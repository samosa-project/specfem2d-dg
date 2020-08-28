! ------------------------------------------------------------ !
! lns_load_background_model                                    !
! ------------------------------------------------------------ !
! Read values of a 2D model and match them on the current mesh.
! Interpolate using a Delaunay triangulation if needed (though quite unstable with ugly meshes).
! Simply match if the model's mesh and the current mesh exactly match.

subroutine lns_load_background_model(nlines_header, nlines_model)

  use constants, only: CUSTOM_REAL, FOUR_THIRDS, NDIM, TINYVAL
  use specfem_par, only: myrank, nglob, USE_DISCONTINUOUS_METHOD, any_elastic, gammaext_dg
  use specfem_par_lns, only: USE_LNS, LNS_rho0, LNS_v0, LNS_p0, &
                             LNS_g, LNS_mu, LNS_eta, LNS_kappa, LNS_E0, LNS_T0, &
                             nabla_v0, sigma_v_0
  
  implicit none
  
  ! Input/output.
  integer, intent(in) :: nlines_header, nlines_model
  
  ! Local variables.
  logical :: meshes_agree
  real(kind=CUSTOM_REAL), parameter :: ZEROcr = 0._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: ONEcr = 1._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: TWOcr = 2._CUSTOM_REAL
  real(kind=CUSTOM_REAL), dimension(NDIM, nlines_model) :: X_m, v_model
  real(kind=CUSTOM_REAL), dimension(nlines_model) :: rho_model, p_model, g_model, gam_model, mu_model, kappa_model
  integer, dimension(nglob) :: idL_m ! Maps (i, j, ispec) to corresponding line of model, if the model agrees with the current mesh.
  
  ! Safeguard.
  if(.not. USE_DISCONTINUOUS_METHOD) then
    write(*,*) "********************************"
    write(*,*) "*            ERROR             *"
    write(*,*) "********************************"
    write(*,*) "* Cannot use the               *"
    write(*,*) "* 'LNS_generalised' model if   *"
    write(*,*) "* not using the DG method.     *"
    write(*,*) "********************************"
    call exit_MPI(myrank, " ")
  endif
  if(.not. USE_LNS) then
    write(*,*) "********************************"
    write(*,*) "*            ERROR             *"
    write(*,*) "********************************"
    write(*,*) "* Cannot use the               *"
    write(*,*) "* 'LNS_generalised' model if   *"
    write(*,*) "* not using the LNS            *"
    write(*,*) "* implementation.              *"
    write(*,*) "********************************"
    call exit_MPI(myrank, " ")
  endif
  if(any_elastic) then
    write(*,*) "********************************"
    write(*,*) "*            ERROR             *"
    write(*,*) "********************************"
    write(*,*) "* Cannot yet use the           *"
    write(*,*) "* 'LNS_generalised' model if   *"
    write(*,*) "* using elastic elements.      *"
    write(*,*) "********************************"
    call exit_MPI(myrank, " ")
  endif
  
  ! Read and store values of model.
  call lns_read_background_model(nlines_header, nlines_model, X_m, &
                                 rho_model, v_model, p_model, g_model, gam_model, mu_model, kappa_model)
  
  ! Initialise initial state registers.
  LNS_rho0   = ZEROcr
  LNS_v0     = ZEROcr
  LNS_E0     = ZEROcr
  ! Initialise initial auxiliary quantities.
  LNS_T0     = ZEROcr
  LNS_p0     = ZEROcr
  nabla_v0   = ZEROcr ! Will be initialised in 'compute_forces_acoustic_LNS_calling_routine.f90'.
  sigma_v_0  = ZEROcr ! Will be initialised in 'compute_forces_acoustic_LNS_calling_routine.f90'.
  ! Physical parameters.
  LNS_g      = ZEROcr
  LNS_mu     = ZEROcr
  LNS_eta    = ZEROcr
  LNS_kappa  = ZEROcr
  
  ! Check whether current mesh agrees with model mesh (in which case, no interpolation is needed).
  meshes_agree = .false.
  call do_meshes_agree(nlines_model, X_m, meshes_agree, idL_m)
  
  if(meshes_agree) then
    if(myrank==0) then
      write(*,*) "> The provided model points exactly agree with current mesh, fast-forwarding by applying directly model to mesh."
    endif
    call apply_model_to_mesh(nlines_model, X_m, &
                             rho_model, v_model, p_model, g_model, gam_model, mu_model, kappa_model, &
                             idL_m)
  else
    ! If mesh does not agree with model, perform linear interpolation.
    if(myrank==0) then
      write(*,*) "> Performing a 2D linear interpolation on the Delaunay triangulation of the provided model points."
    endif
    call delaunay_interp_all_points(nlines_model, X_m, &
                                    rho_model, v_model, p_model, g_model, gam_model, mu_model, kappa_model)
  endif
  
  ! Deduce remaining quantities.    
  LNS_eta = FOUR_THIRDS*LNS_mu
  where(LNS_p0 < TINYVAL) LNS_p0 = ONEcr ! LNS_p0 is uninitialised in solids (remains 0). Crashes compute_E. Hack it.
  where(gammaext_DG < TINYVAL) gammaext_DG = TWOcr ! gammaext_DG is uninitialised in solids (remains 0). Crashes compute_E. Hack it.
  call compute_E(LNS_rho0, LNS_v0, LNS_p0, LNS_E0)
  where(LNS_rho0 < TINYVAL) LNS_rho0 = ONEcr ! LNS_rho0 is uninitialised in solids (remains 0). Crashes compute_T. Hack it.
  call compute_T(LNS_rho0, LNS_v0, LNS_E0, LNS_T0)
  
  ! Safeguards.
  call LNS_prevent_nonsense()
  
  if(any_elastic) then
    ! Update elastic parts.
    call external_DG_update_elastic_from_parfile() ! Update elastic regions by reading parameters directly from parfile (see 'define_external_model.f90').
  endif
  
  ! Eventually save interpolated model.
  !call output_lns_interpolated_model() ! To ASCII (human-readable but very heavy).
  call output_lns_interpolated_model_binary() ! To binary.
  
end subroutine lns_load_background_model


! ------------------------------------------------------------ !
! do_meshes_agree                                              !
! ------------------------------------------------------------ !
! Checks if the given model mesh (xmodel) matches with the current SPECFEM mesh.
! Runs through the current mesh, and for each point of interest loop on all points of the model to find a possible match.
! This is very sub-optimal, but works alright.
! Produces a mapping table along the way, mapping (i, j, ispec) to the corresponding line of model.

subroutine do_meshes_agree(nlines_model, xmodel, meshes_agree, idL_m)
  use constants, only: CUSTOM_REAL, NDIM, NGLLX, NGLLZ, TINYVAL
  use specfem_par, only: nspec, nglob, coord, ispec_is_elastic, ispec_is_acoustic_DG, ibool_before_perio, coord_interface
  use specfem_par_lns, only: norm2r1
  
  implicit none
  
  ! Input/output.
  integer, intent(in) :: nlines_model
  real(kind=CUSTOM_REAL), dimension(NDIM, nlines_model), intent(in) :: xmodel
  logical, intent(out) :: meshes_agree
  integer, dimension(nglob), intent(out) :: idL_m ! Maps (i, j, ispec) to corresponding line of model.
  
  ! Local variables.
  real(kind=CUSTOM_REAL), parameter :: ZEROcr = 0._CUSTOM_REAL
  integer :: ispec, i, j, k, cur_ibool
  logical, dimension(nglob) :: FOUND
  
  FOUND = .false.
  idL_m = 0
  
  do ispec = 1, nspec
    if(ispec_is_elastic(ispec)) then
      do j = 1, NGLLZ; do i = 1, NGLLX
        FOUND(ibool_before_perio(i, j, ispec)) = .true. ! Do not care about elastic points.
      enddo; enddo
    elseif(ispec_is_acoustic_DG(ispec)) then
      ! For DG elements, go through GLL points one by one.
      do j = 1, NGLLZ; do i = 1, NGLLX
        do k = 1, nlines_model
          cur_ibool = ibool_before_perio(i, j, ispec)
          if(norm2r1(xmodel(:, k)-(coord(:, cur_ibool)-(/ZEROcr, coord_interface/)))<=TINYVAL) then
            FOUND(cur_ibool) = .true.
            idL_m(cur_ibool) = k
          endif
        enddo
      enddo; enddo
    endif
  enddo
  ! Meshes agree if all points (in this CPU) were found in the model.
  meshes_agree = all(FOUND)
end subroutine do_meshes_agree


! ------------------------------------------------------------ !
! apply_model_to_mesh                                          !
! ------------------------------------------------------------ !
! Used if the meshes agreed.
! Simply maps the loaded model to the current SPECFEM mesh, using the previously computed mapping.

subroutine apply_model_to_mesh(nlines_model, X_m, &
                               rho_m, v_m, p_m, g_m, gam_m, mu_m, kappa_m, &
                               idL_m)
  use constants, only: CUSTOM_REAL, NDIM, NGLLX, NGLLZ, FOUR_THIRDS, TINYVAL
  use specfem_par_lns, only: LNS_rho0, LNS_v0, LNS_p0, &
                             LNS_g, LNS_mu, LNS_kappa, LNS_c0
  use specfem_par, only: nspec, ispec_is_elastic, rhoext, vpext, myrank, &
                         ispec_is_acoustic_DG, ibool_DG, ibool_before_perio, gammaext_DG, nglob
  
  implicit none
  
  ! Input/output.
  integer, intent(in) :: nlines_model
  real(kind=CUSTOM_REAL), dimension(NDIM, nlines_model), intent(inout) :: X_m
  real(kind=8), dimension(NDIM, nlines_model), intent(in) :: v_m
  real(kind=CUSTOM_REAL), dimension(nlines_model), intent(in) :: rho_m, p_m, g_m, gam_m, mu_m, kappa_m
  integer, dimension(nglob), intent(in) :: idL_m ! Maps (i, j, ispec) to corresponding line of model.
  
  ! Local variables.
  real(kind=CUSTOM_REAL), parameter :: ZEROcr = 0._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: ONEcr = 1._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: TWOcr = 2._CUSTOM_REAL
  integer :: ispec, i, j, v, l_model, iglobDG
  
  do ispec = 1, nspec
    if(ispec_is_elastic(ispec)) then
      cycle ! Elastic regions are treated separately (see call to 'external_DG_update_elastic_from_parfile' in the 'lns_load_background_model' routine).
    elseif(ispec_is_acoustic_DG(ispec)) then
      ! For DG elements, go through GLL points one by one.
      do j = 1, NGLLZ; do i = 1, NGLLX
        l_model = idL_m(ibool_before_perio(i, j, ispec))
        iglobDG = ibool_DG(i, j, ispec)
        LNS_rho0(iglobDG)    = rho_m(l_model)
        do v = 1, NDIM
          LNS_v0(v, iglobDG) = v_m(v, l_model)
        enddo
        LNS_p0(iglobDG)      = p_m(l_model)
        LNS_g(iglobDG)       = g_m(l_model)
        gammaext_DG(iglobDG) = gam_m(l_model)
        LNS_mu(iglobDG)      = mu_m(l_model)
        LNS_kappa(iglobDG)   = kappa_m(l_model)
        ! For the subroutines in 'invert_mass_matrix.f90', one needs to initialise the following:
        rhoext(i, j, ispec) = LNS_rho0(iglobDG)
        LNS_c0(iglobDG) = sqrt(gammaext_DG(iglobDG)*LNS_p0(iglobDG)/LNS_rho0(iglobDG)) ! Take the chance to compute and save c0.
        vpext(i, j, ispec) = LNS_c0(iglobDG)
      enddo; enddo
    else
      ! Neither elastic nor acoustic_dg.
      write(*,*) "********************************"
      write(*,*) "*            ERROR             *"
      write(*,*) "********************************"
      write(*,*) "* See                          *"
      write(*,*) "* lns_load_background_model.f90"
      write(*,*) "********************************"
      call exit_MPI(myrank, " ")
    endif ! Endif on ispec_is_acoustic_DG.
  enddo
end subroutine apply_model_to_mesh


! ------------------------------------------------------------ !
! lns_read_background_model                                    !
! ------------------------------------------------------------ !
! Read and store values of model from the expected file name (BCKGRD_MDL_LNS_FILENAME from specfem_par_lns).
! If the model file is ASCII (BCKGRD_MDL_LNS_is_binary==.false.), it should have the following format:
! --- file begin
!   //header whatever//
!   //header whatever//
!   //header whatever//
!   x[m] z[m] rho0[kg.m^{-3}] v0x[m.s^{-1}] v0z[m.s^{-1}] p0[Pa] g[m.s^{-2}] gamma[1] mu[kg.m^{-1}.s^{-1}] kappa[m.kg.s^{-3}.K{-1}]
! --- file end
! If the model file is binary (BCKGRD_MDL_LNS_is_binary==.true.), the number of lines should have been read before (see nlines_model variable which is an input), and this routine will simply read from the binary stream.
! The Matlab script './utils_new/lns_background_models/write_bg_model.m' allows the printing of models matching the aforementionned requirements.

subroutine lns_read_background_model(nlines_header, nlines_model, X_m, &
                                     rho_model, v_model, p_model, g_model, gam_model, mu_model, kappa_model)
  use constants, only: CUSTOM_REAL, NDIM
  use specfem_par, only: myrank
  use specfem_par_lns, only: BCKGRD_MDL_LNS_is_binary, BCKGRD_MDL_LNS_FILENAME, BCKGRD_MDL_LNS_NCOL
  
  implicit none
  
  ! Input/output.
  integer, intent(in) :: nlines_header, nlines_model
  real(kind=CUSTOM_REAL), dimension(NDIM, nlines_model), intent(out) :: X_m, v_model
  real(kind=CUSTOM_REAL), dimension(nlines_model), intent(out) :: rho_model, p_model, g_model, gam_model, mu_model, kappa_model
  
  ! Local variables.
  integer :: ncolumns_detected
  integer :: i, io
  
  call external_model_DG_testfile(BCKGRD_MDL_LNS_FILENAME) ! Test-reading the file.
  
  if(.not. BCKGRD_MDL_LNS_is_binary) then
    call external_model_DG_only_find_nbcols(BCKGRD_MDL_LNS_FILENAME, nlines_header, ncolumns_detected)
  endif
  
  if(myrank==0) then
    write(*,*) "> Reading general background model file '", trim(BCKGRD_MDL_LNS_FILENAME)
    if(.not. BCKGRD_MDL_LNS_is_binary) then
      write(*,*) "> Detected ",nlines_header," lines for header, and ",nlines_model," lines of model ", &
                 "(total = ",(nlines_header+nlines_model)," lines). Detected ", ncolumns_detected," columns."
    endif
  endif
  
  if(BCKGRD_MDL_LNS_is_binary) then
    ! If binary, loop until whole file is read.
    open(100, file=BCKGRD_MDL_LNS_FILENAME, access='stream', form='unformatted', STATUS="old", action='read', iostat=io)
    if (io/=0) call exit_MPI(myrank, "Error opening background model file.")
    
    ! Read as many n-uplets as needed from the serial binary file.
    do i=1, nlines_model
      read(100, iostat=io) X_m(1, i), X_m(NDIM, i), &
                           rho_model(i), v_model(1, i), v_model(NDIM, i), p_model(i), &
                           g_model(i), gam_model(i), mu_model(i), kappa_model(i)
    enddo

  else
    ! If ASCII, skip header and read detected number of lines.
    OPEN(100, file=BCKGRD_MDL_LNS_FILENAME)
    do i=1, nlines_header
      read(100, *, iostat=io)
      if (io/=0) call exit_MPI(myrank, "Error reading line in background model file.")
    enddo
    do i=1, nlines_model
      if(ncolumns_detected==BCKGRD_MDL_LNS_NCOL) then
        read(100, *, iostat=io) X_m(1, i), X_m(NDIM, i), &
                                rho_model(i), v_model(1, i), v_model(NDIM, i), p_model(i), &
                                g_model(i), gam_model(i), mu_model(i), kappa_model(i)
      else
        write(*,*) "********************************"
        write(*,*) "*            ERROR             *"
        write(*,*) "********************************"
        write(*,*) "* Number of columns in model   *"
        write(*,*) "* file is wrong.               *"
        write(*,*) "********************************"
        call exit_MPI(myrank, " ")
      endif
      if(io/=0) exit
    enddo
    
  endif
  
  close(100)
end subroutine lns_read_background_model


! ------------------------------------------------------------ !
! delaunay_interp_all_points                                   !
! ------------------------------------------------------------ !
! Interpolate the model on the SPECFEM mesh using a Delaunay triangulation.

subroutine delaunay_interp_all_points(nlines_model, X_m, &
                                      rho_m, v_m, p_m, g_m, gam_m, mu_m, kappa_m)
  use constants, only: CUSTOM_REAL, NDIM, NGLLX, NGLLZ, TINYVAL
  !use specfem_par_lns, only: USE_LNS
  use specfem_par, only: myrank, nspec, ispec_is_elastic, ispec_is_acoustic_DG
  
  implicit none
  
  ! Input/output.
  integer, intent(in) :: nlines_model
  real(kind=CUSTOM_REAL), dimension(NDIM, nlines_model), intent(inout) :: X_m
  real(kind=8), dimension(NDIM, nlines_model), intent(in) :: v_m
  real(kind=CUSTOM_REAL), dimension(nlines_model), intent(in) :: rho_m, p_m, g_m, gam_m, mu_m, kappa_m
  
  ! Local variables.
  real(kind=CUSTOM_REAL), parameter :: ZEROcr = 0._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: ONEcr = 1._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: TWOcr = 2._CUSTOM_REAL
  integer :: ispec, i, j
  integer(kind=4) :: tri_num
  integer(kind=4), dimension(3, 2*nlines_model) :: tri_vert, tri_nabe ! Theoretically, tri_vert and tri_nabe are of size (3, tri_num). However, at this point, one does not know tri_num. Nevertheless, an upper bound is 2*nlines_model (cf. dtris2 subroutine in 'table_delaunay').
  
  call delaunay_background_model_2d(nlines_model, X_m, tri_num, tri_vert, tri_nabe)
  
  if(myrank==0) then
    write(*,*) "> > Starting interpolation on the Delaunay triangulation."
  endif
  
  do ispec = 1, nspec
    if(ispec_is_elastic(ispec)) then
      cycle ! Elastic regions are treated separately (see call to 'external_DG_update_elastic_from_parfile' below).
    elseif(ispec_is_acoustic_DG(ispec)) then
      ! For DG elements, go through GLL points one by one.
      do j = 1, NGLLZ; do i = 1, NGLLX
          call delaunay_interpolate_one_point(nlines_model, X_m, tri_num, tri_vert, ispec, i, j, &
                                              rho_m, v_m, p_m, g_m, gam_m, mu_m, kappa_m)
      enddo; enddo
    else
      ! Neither elastic nor acoustic_dg.
      write(*,*) "********************************"
      write(*,*) "*            ERROR             *"
      write(*,*) "********************************"
      write(*,*) "* See                          *"
      write(*,*) "* lns_load_background_model.f90"
      write(*,*) "********************************"
      call exit_MPI(myrank, " ")
    endif ! Endif on ispec_is_acoustic_DG.
  enddo
end subroutine delaunay_interp_all_points


! ------------------------------------------------------------ !
! delaunay_background_model_2d                                 !
! ------------------------------------------------------------ !
! Compute the Delaunay triangulation of a given model.
! This uses the external module 'dtris2' in the 'table_delaunay.f90' file.

subroutine delaunay_background_model_2d(nlines_model, X_m, tri_num, tri_vert, tri_nabe)
  use constants, only: CUSTOM_REAL, NDIM
  use specfem_par, only: myrank
  
  implicit none
  ! Input/output.
  integer, intent(in) :: nlines_model
  real(kind=CUSTOM_REAL), dimension(NDIM, nlines_model), intent(in) :: X_m
  integer(kind=4), intent(out) :: tri_num
  integer(kind=4), dimension(3, 2*nlines_model), intent(out) :: tri_vert, tri_nabe ! Theoretically, tri_vert and tri_nabe are of size (3, tri_num). However, at this point, one does not know tri_num. Nevertheless, an upper bound is 2*nlines_model (cf. dtris2 subroutine in 'table_delaunay').
  ! Local variables.
  integer :: t
  call dtris2(nlines_model, X_m, tri_num, tri_vert, tri_nabe)
  if(myrank==0) then
    write(*,*) "> > Delaunay triangulation of the ",nlines_model," model points: ", tri_num, "triangles found."
    ! Debug: print out triangles.
    if(.false.) then
      do t = 1, tri_num
        write(*,*) "> > Triangle n°", t, " vertices: ", X_m(1:NDIM, tri_vert(1, t))
        write(*,*) "> >                                   ", X_m(1:NDIM, tri_vert(2, t))
        write(*,*) "> >                                   ", X_m(1:NDIM, tri_vert(3, t))
      enddo
    endif
  endif
end subroutine delaunay_background_model_2d


! ------------------------------------------------------------ !
! delaunay_interpolate_one_point                               !
! ------------------------------------------------------------ !
! Given a Delaunay triangulation, interpolate one value.
! Algorithm:
!   for each triangles in the triangulation,
!     if the queried point is in the triangle,
!       perform a barycentric interpolation using the three vertices of the triangle;
!     else, skip to next triangle.
! This type of interpolation is dodgy, because it is only C^1 per triangle, and only C^0 across triangles.
! In other words, derivative discontinuities might (and surely will) occur across the triangle edges.
! It is best to avoid doing this, and generate a model exactly matching the SPECFEM mesh.
! TODO, alternatively: a nicer interpolation algorithm making use of the triangulation and of the neighbouring triangles.

subroutine delaunay_interpolate_one_point(nlines_model, X_m, tri_num, tri_vert, ispec, i, j, &
                                          rho_m, v_m, p_m, g_m, gam_m, mu_m, kappa_m)
  use constants, only: CUSTOM_REAL, NDIM, TINYVAL
  use specfem_par_lns, only: LNS_rho0, LNS_v0, LNS_p0, LNS_g, LNS_mu, LNS_kappa, LNS_c0
  use specfem_par, only: ibool, ibool_DG, coord, coord_interface, gammaext_dg, rhoext, vpext, myrank
  
  implicit none
  
  ! Input/output.
  integer, intent(in) :: nlines_model, ispec, i, j, tri_num
  integer(kind=4), dimension(3, 2*nlines_model), intent(in) :: tri_vert ! Theoretically, tri_vert and tri_nabe are of size (3, tri_num). However, at this point, one does not know tri_num. Nevertheless, an upper bound is 2*nlines_model (cf. dtris2 subroutine in 'table_delaunay').
  real(kind=8), dimension(NDIM, nlines_model), intent(in) :: X_m, v_m
  real(kind=CUSTOM_REAL), dimension(nlines_model), intent(in) :: rho_m, p_m, g_m, gam_m, mu_m, kappa_m
  
  ! Local variables.
  real(kind=CUSTOM_REAL), parameter :: ZEROcr = 0._CUSTOM_REAL
  integer :: t, v, iglobDG
  integer, dimension(3) :: loctri_vertices_ids
  real(kind=8), dimension(3, NDIM) :: local_vertice_list
  real(kind=CUSTOM_REAL), dimension(NDIM) :: point_to_test
  logical found, inTri
  real(kind=CUSTOM_REAL), dimension(3) :: barycor
  
  ! Safety zeroing of local variables.
  t = 0
  v = 0
  iglobDG = 0
  loctri_vertices_ids = 0
  local_vertice_list = 0.
  point_to_test = 0.
  found = .false.
  barycor = 0.
  
  ! Bring back the SPECFEM coordinate to the model convention (where the interface is at altitude zero).
  point_to_test = coord(1:NDIM, ibool(i, j, ispec)) - (/ZEROcr, coord_interface/)
  
  iglobDG = ibool_DG(i, j, ispec)
  
  do t = 1, tri_num
    do v = 1, 3
      local_vertice_list(v, 1:NDIM) = X_m(1:NDIM, tri_vert(v, t))
    enddo
    
    call barycentric_coordinates_2d(local_vertice_list, point_to_test, inTri, barycor)
    
    if(inTri) then
      ! Point is inside triangle.
      loctri_vertices_ids = tri_vert(1:3, t)
      LNS_rho0(iglobDG)    = DOT_PRODUCT(barycor, rho_m(loctri_vertices_ids))
      do v = 1, NDIM
        LNS_v0(v, iglobDG) = DOT_PRODUCT(barycor, v_m(v, loctri_vertices_ids))
      enddo
      LNS_p0(iglobDG)      = DOT_PRODUCT(barycor, p_m(loctri_vertices_ids))
      LNS_g(iglobDG)       = DOT_PRODUCT(barycor, g_m(loctri_vertices_ids))
      gammaext_DG(iglobDG) = DOT_PRODUCT(barycor, gam_m(loctri_vertices_ids))
      LNS_mu(iglobDG)      = DOT_PRODUCT(barycor, mu_m(loctri_vertices_ids))
      LNS_kappa(iglobDG)   = DOT_PRODUCT(barycor, kappa_m(loctri_vertices_ids))
      
      ! Sanity checks.
      if(LNS_rho0(iglobDG) <= TINYVAL) then
        write(*,*) "********************************"
        write(*,*) "*            ERROR             *"
        write(*,*) "********************************"
        write(*,*) "* A negative (<=0) density was *"
        write(*,*) "* found after interpolation.   *"
        write(*,*) "********************************"
        write(*,*) ispec,i,j,&
                   coord(1,iglobDG),coord(2,iglobDG),&
                   LNS_rho0(iglobDG)
        write(*,*) "********************************"
        call exit_MPI(myrank, " ")
      endif
      
      ! For the subroutines in 'invert_mass_matrix.f90', one needs to initialise the following:
      rhoext(i, j, ispec) = LNS_rho0(iglobDG)
      LNS_c0(iglobDG) = sqrt(gammaext_DG(iglobDG)*LNS_p0(iglobDG)/LNS_rho0(iglobDG)) ! Take the chance to compute and save c0.
      vpext(i, j, ispec) = LNS_c0(iglobDG)
      
      return
    else
      ! Point is outside triangle, go to next triangle.
      cycle
    endif
  enddo
  
  if(.not. found) then
    write(*,*) "********************************"
    write(*,*) "*            ERROR             *"
    write(*,*) "********************************"
    write(*,*) "* Physical point               *"
    write(*,*) "* ", point_to_test
    write(*,*) "* (coord_interface removed)    *"
    write(*,*) "* was not found within the     *"
    write(*,*) "* Delaunay triangulation.      *"
    write(*,*) "* 1) Check coord_interface is  *"
    write(*,*) "* set properly wrt the model.  *"
    write(*,*) "* 2) You can look at the       *"
    write(*,*) "* triangles by commenting out  *"
    write(*,*) "* the 'if(.false.)' in the     *"
    write(*,*) "* subroutine                   *"
    write(*,*) "* 'delaunay_background_model_2d'"
    write(*,*) "* (in                          *"
    write(*,*) "* 'lns_load_background_model.f90')."
    write(*,*) "* 3) You can check whether the *"
    write(*,*) "* set of points in your model  *"
    write(*,*) "* is not prone to wrong convex *"
    write(*,*) "* hulls (typically, make sure  *"
    write(*,*) "* it can map a rectangular or  *"
    write(*,*) "* triangular grid).            *"
    write(*,*) "********************************"
    call exit_MPI(myrank, " ")
  endif
  
end subroutine delaunay_interpolate_one_point


! ------------------------------------------------------------ !
! barycentric_coordinates_2d                                   !
! ------------------------------------------------------------ !
! https://en.wikipedia.org/wiki/Barycentric_coordinate_system

subroutine barycentric_coordinates_2d(vlist, P, inTri, barycor)
  use constants, only: CUSTOM_REAL, NDIM, TINYVAL
  use specfem_par, only: myrank
  implicit none
  ! Input/output.
  real(kind=8), dimension(3, NDIM), intent(in) :: vlist
  real(kind=CUSTOM_REAL), dimension(NDIM), intent(in) :: P
  logical, intent(out) :: inTri
  real(kind=CUSTOM_REAL), dimension(3), intent(out) :: barycor
  ! Local variables.
  real(kind=CUSTOM_REAL), parameter :: ZERO = 0._CUSTOM_REAL
  real(kind=CUSTOM_REAL) :: det
  inTri = .false.
  barycor = ZERO
  det = (vlist(1, 1)-vlist(3, 1))*(vlist(2, 2)-vlist(3, 2)) - (vlist(1, 2)-vlist(3, 2))*(vlist(2, 1)-vlist(3, 1))
  if(abs(det)<TINYVAL) then
    write(*,*) "********************************"
    write(*,*) "*            ERROR             *"
    write(*,*) "********************************"
    write(*,*) "* Determinant is very close to *"
    write(*,*) "* zero in this triangle.       *"
    write(*,*) "* Something's wrong I can feel *"
    write(*,*) "* it.                          *"
    write(*,*) "********************************"
    call exit_MPI(myrank, " ")
  endif
  barycor(1) = ((vlist(2, 2)-vlist(3, 2))*(P(1)-vlist(3, 1))+(vlist(3, 1)-vlist(2, 1))*(P(2)-vlist(3, 2))) / det
  if(barycor(1)<-TINYVAL) return ! Point is outside triangle, stop and return.
  barycor(2) = ((vlist(3, 2)-vlist(1, 2))*(P(1)-vlist(3, 1))+(vlist(1, 1)-vlist(3, 1))*(P(2)-vlist(3, 2))) / det
  if(barycor(2)<-TINYVAL) return ! Point is outside triangle, stop and return.
  barycor(3) = 1. - sum(barycor(1:2))
  if(barycor(3)<-TINYVAL) return ! Point is outside triangle, stop and return.
  inTri = .true. ! No return command was hit before, hence point can be considered in triangle.
end subroutine barycentric_coordinates_2d


! ------------------------------------------------------------ !
! output_lns_interpolated_model                                !
! ------------------------------------------------------------ !
! Outputs the loaded model in full to an ASCII file, for checking/plotting purposes.

subroutine output_lns_interpolated_model()
  use constants, only: CUSTOM_REAL, NDIM, NGLLX, NGLLZ
  use specfem_par, only: myrank, nspec, ibool_before_perio, coord, ibool_DG, gammaext_DG
  use specfem_par_lns, only: LNS_rho0, LNS_v0, LNS_p0, LNS_g, LNS_mu, LNS_kappa
  implicit none
  ! Local variables.
  integer :: ispec, i, j, iglobDG
  ! Let all procs write to file (I know this may not be sequential and freak up the file, but I don't care since it's only for debugging purposes).
  open(unit=504,file='OUTPUT_FILES/LNS_GENERAL_MODEL',status='unknown',action='write', position="append")
  do ispec = 1,nspec
    do j = 1,NGLLZ
      do i = 1,NGLLX
        iglobDG = ibool_DG(i, j, ispec)
        write(504,*) coord(1, ibool_before_perio(i, j, ispec)), coord(2, ibool_before_perio(i, j, ispec)),&
                     LNS_rho0(iglobDG), LNS_v0(1, iglobDG), LNS_v0(NDIM, iglobDG), LNS_p0(iglobDG), &
                     LNS_g(iglobDG), gammaext_DG(iglobDG), LNS_mu(iglobDG), LNS_kappa(iglobDG)
      enddo
    enddo
  enddo
  close(504)
  if(myrank==0) then
    write(*,*) "> > Dumped interpolated model (on proc 0) to './OUTPUT_FILES/LNS_GENERAL_MODEL'. ", &
               "Use the following Matlab one-liner to plot:"
    write(*,*) "      OFD=''; [~,lab]=order_bg_model();a=importdata(OFD);x=a(:,1);z=a(:,2);close all;for i=3:10;", &
               "d=[];d.rho=[];d.vel=[];d.pre=a(:,i);[Xi,Yi,Zi]=interpDumps(x,z,d,0,0,1);figure(i);surf(Xi,Yi,Zi.pre);", &
               "view([0,0,1]);shading flat;colorbar;title(lab{i});end;"
  endif
end subroutine output_lns_interpolated_model


! ------------------------------------------------------------ !
! output_lns_interpolated_model_binary                         !
! ------------------------------------------------------------ !
! Same as 'output_lns_interpolated_model', but to a binary file.
! Use the Matlab script './utils_new/lns_background_models/load_bg_model.m' to read the output file.

subroutine output_lns_interpolated_model_binary()
  use constants, only: CUSTOM_REAL, NDIM, NGLLX, NGLLZ, SIZE_DOUBLE
  use specfem_par, only: myrank, nglob, ibool, nspec, ibool_DG, coord, ibool_before_perio, gammaext_DG
  use specfem_par_lns, only: LNS_rho0, LNS_v0, LNS_p0, LNS_g, LNS_mu, LNS_kappa
  implicit none
  ! Local variables.
  character(len=150) :: THEFILE
  integer, parameter :: THEUNIT = 504
  integer :: ispec, i, j, iglobDG, ier, icounter, iglob
  logical, dimension(nglob) :: mask_ibool
  write(THEFILE, "('OUTPUT_FILES/LNS_GENERAL_MODEL_',i3.3,'.bin')") myrank
  open(unit=THEUNIT, file=THEFILE, form='unformatted', access='direct', status='unknown', &
       action='write', recl=(2*NDIM+6)*SIZE_DOUBLE, iostat=ier)
  icounter = 0
  mask_ibool(:) = .false.
  do ispec = 1,nspec
    do j = 1,NGLLZ
      do i = 1,NGLLX
        iglob = ibool(i, j, ispec)
        if (.not. mask_ibool(iglob)) then
          icounter = icounter + 1
          mask_ibool(iglob) = .true.
          iglobDG = ibool_DG(i, j, ispec)
          write(THEUNIT, rec=icounter) dble(coord(:, ibool_before_perio(i, j, ispec))), &
                                       dble(LNS_rho0(iglobDG)), dble(LNS_v0(:, iglobDG)), dble(LNS_p0(iglobDG)), &
                                       dble(gammaext_DG(iglobDG)), dble(LNS_g(iglobDG)), &
                                       dble(LNS_mu(iglobDG)), dble(LNS_kappa(iglobDG))
        endif
      enddo
    enddo
  enddo
  close(THEUNIT)
end subroutine output_lns_interpolated_model_binary



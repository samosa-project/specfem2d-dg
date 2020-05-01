!========================================================================
!
!                   S P E C F E M 2 D  Version 7 . 0
!                   --------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!                and CNRS / University of Marseille, France
!                 (there are currently many more authors!)
! (c) Princeton University and CNRS / University of Marseille, April 2014
!
! This software is a computer program whose purpose is to solve
! the two-dimensional viscoelastic anisotropic or poroelastic wave equation
! using a spectral-element method (SEM).
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
! The full text of the license is available in file "LICENSE".
!
!========================================================================

! ------------------------------------------------------------ !
! lns_load_background_model                                    !
! ------------------------------------------------------------ !
! Read values of model and interpolate on current grid.

subroutine lns_load_background_model(nlines_header, nlines_model)

  use constants, only: CUSTOM_REAL, NDIM
  use specfem_par, only: myrank, USE_DISCONTINUOUS_METHOD, any_elastic
  use specfem_par_lns, only: USE_LNS
  
  implicit none
  
  ! Input/output.
  integer, intent(in) :: nlines_header, nlines_model
  
  ! Local variables.
  real(kind=CUSTOM_REAL), parameter :: ZERO = 0._CUSTOM_REAL
  real(kind=CUSTOM_REAL), dimension(NDIM, nlines_model) :: X_m, v_model
  real(kind=CUSTOM_REAL), dimension(nlines_model) :: rho_model, p_model, g_model, gam_model, mu_model, kappa_model
  
  ! Safeguard.
  if(.not. USE_DISCONTINUOUS_METHOD) then
    write(*,*) "********************************"
    write(*,*) "*            ERROR             *"
    write(*,*) "********************************"
    write(*,*) "* Cannot use the               *"
    write(*,*) "* 'LNS_generalised' model if   *"
    write(*,*) "* not using the DG method.     *"
    write(*,*) "********************************"
    stop
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
    stop
  endif
  if(any_elastic) then
    write(*,*) "********************************"
    write(*,*) "*            ERROR             *"
    write(*,*) "********************************"
    write(*,*) "* Cannot yet use the           *"
    write(*,*) "* 'LNS_generalised' model if   *"
    write(*,*) "* using elastic elements.      *"
    write(*,*) "********************************"
    stop
  endif
  
  ! Read and store values of model.
  call lns_read_background_model(nlines_header, nlines_model, X_m, &
                                 rho_model, v_model, p_model, g_model, gam_model, mu_model, kappa_model &
                                )
  
  ! Check
!  if(.true.) then
!    do i=1, nlines_model
!      write(*,*) X_m(1, i), X_m(2, i), p_model(i)
!    enddo
!    stop
!  endif

  ! Interpolation.
  if(myrank==0) then
    write(*,*) "> Performing a 2D linear interpolation on the Delaunay triangulation of the provided model points."
  endif
  call delaunay_interp_all_points(nlines_model, X_m, &
                                  rho_model, v_model, p_model, g_model, gam_model, mu_model, kappa_model &
                                 )
  
  ! Update elastic parts.
  call external_DG_update_elastic_from_parfile() ! Update elastic regions by reading parameters directly from parfile (see 'define_external_model.f90').
  
  ! Eventually save interpolated model.
  call output_lns_interpolated_model()
  
end subroutine lns_load_background_model


! ------------------------------------------------------------ !
! lns_read_background_model                                    !
! ------------------------------------------------------------ !
! Read and store values of model.
! The model file, if ASCII, should have the following format:
! --- file begin
!   //header whatever//
!   //header whatever//
!   //header whatever//
!   x[m] z[m] rho0[kg.m^{-3}] v0x[m.s^{-1}] v0z[m.s^{-1}] p0[Pa] g[m.s^{-2}] gamma[1] mu[kg.m^{-1}.s^{-1}] kappa[m.kg.s^{-3}.K{-1}]
! --- file end

subroutine lns_read_background_model(nlines_header, nlines_model, X_m, &
                                     rho_model, v_model, p_model, g_model, gam_model, mu_model, kappa_model &
                                    )
  use constants, only: CUSTOM_REAL, NDIM
  !use specfem_par,only: nspec, tau_sigma, tau_epsilon
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
    if (io/=0) stop "Error opening background model file."
    
    ! Read as many n-uplets as needed from the serial binary file.
    do i=1, nlines_model
      read(100, iostat=io) X_m(1, i), X_m(NDIM, i), &
                           rho_model(i), v_model(1, i), v_model(NDIM, i), p_model(i), &
                           g_model(i), gam_model(i), mu_model(i), kappa_model(i)
!      write(*,*) X_m(1, i), X_m(NDIM, i), &
!                 rho_model(i), v_model(1, i), v_model(NDIM, i), p_model(i), &
!                 g_model(i), gam_model(i), mu_model(i), kappa_model(i), io
    enddo

  else
    ! If ASCII, skip header and read detected number of lines.
    OPEN(100, file=BCKGRD_MDL_LNS_FILENAME)
    do i=1, nlines_header
      read(100, *, iostat=io)
      IF (io/=0) stop "Error reading line in background model file."
    enddo
    do i=1, nlines_model
      if(ncolumns_detected==BCKGRD_MDL_LNS_NCOL) then
        read(100, *, iostat=io) X_m(1, i), X_m(NDIM, i), &
                                rho_model(i), v_model(1, i), v_model(NDIM, i), p_model(i), &
                                g_model(i), gam_model(i), mu_model(i), kappa_model(i)
!        if(i<10) write(*,*) i, X_m(1, i), X_m(NDIM, i), p_model(i)
!        if(i>10) stop
      else
        write(*,*) "********************************"
        write(*,*) "*            ERROR             *"
        write(*,*) "********************************"
        write(*,*) "* Number of columns in model   *"
        write(*,*) "* file is wrong.               *"
        write(*,*) "********************************"
        stop
      endif
      if(io/=0) exit
    enddo
    
  endif
  
  close(100)
end subroutine lns_read_background_model

! ------------------------------------------------------------ !
! delaunay_interp_all_points                                   !
! ------------------------------------------------------------ !


subroutine delaunay_interp_all_points(nlines_model, X_m, &
                                      rho_m, v_m, p_m, g_m, gam_m, mu_m, kappa_m &
                                     )
  use constants, only: CUSTOM_REAL, NDIM, NGLLX, NGLLZ, FOUR_THIRDS, TINYVAL
  use specfem_par_lns, only: USE_LNS, LNS_rho0, LNS_v0, LNS_p0, &
                             LNS_g, LNS_mu, LNS_eta, LNS_kappa, LNS_E0, LNS_T0, &
                             nabla_v0, sigma_v_0
  use specfem_par, only: USE_DISCONTINUOUS_METHOD, any_acoustic_DG, myrank, nspec, ispec_is_elastic, &
                         ispec_is_acoustic_DG, coord, ibool_DG, ibool, gammaext_DG
  
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
  
  ! Safeguards.
  if(.not. (USE_DISCONTINUOUS_METHOD .and. USE_LNS)) then
    write(*,*) "********************************"
    write(*,*) "*            ERROR             *"
    write(*,*) "********************************"
    write(*,*) "* Currently cannot use the     *"
    write(*,*) "* 'external_DG' model if not   *"
    write(*,*) "* using the DG method.         *"
    write(*,*) "********************************"
    stop
  endif
  if(.not. any_acoustic_DG) then
    write(*,*) "********************************"
    write(*,*) "*            ERROR             *"
    write(*,*) "********************************"
    write(*,*) "* No acoustic DG elements!     *"
    write(*,*) "* Cannot import                *"
    write(*,*) "* LNS_generalised model.       *"
    write(*,*) "********************************"
    stop
  endif
  
  call delaunay_background_model_2d(nlines_model, X_m, tri_num, tri_vert, tri_nabe)
  
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
  
  if(myrank==0) then
    write(*,*) "> > Starting interpolation on the Delaunay triangulation."
  endif
  
  do ispec = 1, nspec
    if(ispec_is_elastic(ispec)) then
!      write(*,*) 'Element ',ispec,' is elastic.'! DEBUG
      cycle ! Elastic regions are treated separately (see call to 'external_DG_update_elastic_from_parfile' below).
    elseif(ispec_is_acoustic_DG(ispec)) then
      ! For DG elements, go through GLL points one by one.
!      write(*,*) 'Element ',ispec,' is acoustic DG.'! DEBUG
      do j = 1, NGLLZ
        do i = 1, NGLLX
          call delaunay_interpolate_one_point(nlines_model, X_m, tri_num, tri_vert, ispec, i, j, &
                                              rho_m, v_m, p_m, g_m, gam_m, mu_m, kappa_m &
                                             )
        enddo
      enddo
    else
      ! Neither elastic nor acoustic_dg.
      write(*,*) "********************************"
      write(*,*) "*            ERROR             *"
      write(*,*) "********************************"
      write(*,*) "* lns_load_background_model.f90"
      write(*,*) "********************************"
      stop
    endif ! Endif on ispec_is_acoustic_DG.
  enddo
  
  ! Deduce remaining quantities.    
  LNS_eta = FOUR_THIRDS*LNS_mu
  where(LNS_p0 < TINYVAL) LNS_p0 = ONEcr ! LNS_p0 is uninitialised in solids (remains 0). Crashes compute_E. Hack it.
  where(gammaext_DG < TINYVAL) gammaext_DG = TWOcr ! gammaext_DG is uninitialised in solids (remains 0). Crashes compute_E. Hack it.
  call compute_E(LNS_rho0, LNS_v0, LNS_p0, LNS_E0)
  where(LNS_rho0 < TINYVAL) LNS_rho0 = ONEcr ! LNS_rho0 is uninitialised in solids (remains 0). Crashes compute_T. Hack it.
  call compute_T(LNS_rho0, LNS_v0, LNS_E0, LNS_T0)
  
  ! Debug.
  do ispec = 1, nspec; do j = 1, NGLLZ; do i = 1, NGLLX
    if(LNS_E0(ibool_DG(i, j, ispec))<TINYVAL) write(*,*) 'E0<=0 at [',coord(1:NDIM, ibool(i, j, ispec)),']: ', &
                                                         LNS_E0(ibool_DG(i, j, ispec)), gammaext_DG(ibool_DG(i, j, ispec))
  enddo; enddo; enddo
  
  ! Safeguards.
  call LNS_prevent_nonsense()
!  write(*,*) 'min max E0 reldiff', minval(LNS_E0), maxval(LNS_E0), ABS(maxval(LNS_E0)-minval(LNS_E0))/abs(maxval(LNS_E0))
!  write(*,*) 'min max T0 reldiff', minval(LNS_T0), maxval(LNS_T0), ABS(maxval(LNS_T0)-minval(LNS_T0))/abs(maxval(LNS_T0))
!  stop
end subroutine delaunay_interp_all_points


! ------------------------------------------------------------ !
! delaunay_background_model_2d                                 !
! ------------------------------------------------------------ !


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


subroutine delaunay_interpolate_one_point(nlines_model, X_m, tri_num, tri_vert, ispec, i, j, &
                                          rho_m, v_m, p_m, g_m, gam_m, mu_m, kappa_m &
                                         )
  use constants, only: CUSTOM_REAL, NDIM, TINYVAL
  use specfem_par_lns, only: LNS_rho0, LNS_v0, LNS_p0, LNS_g, LNS_mu, LNS_kappa, LNS_c0
  use specfem_par, only: ibool, ibool_DG, coord, coord_interface, gammaext_dg, rhoext, vpext
  
  implicit none
  
  ! Input/output.
  real(kind=CUSTOM_REAL), parameter :: ZEROcr = 0._CUSTOM_REAL
  integer, intent(in) :: nlines_model, ispec, i, j, tri_num
  integer(kind=4), dimension(3, 2*nlines_model), intent(in) :: tri_vert ! Theoretically, tri_vert and tri_nabe are of size (3, tri_num). However, at this point, one does not know tri_num. Nevertheless, an upper bound is 2*nlines_model (cf. dtris2 subroutine in 'table_delaunay').
  real(kind=8), dimension(NDIM, nlines_model), intent(in) :: X_m, v_m
  real(kind=CUSTOM_REAL), dimension(nlines_model), intent(in) :: rho_m, p_m, g_m, gam_m, mu_m, kappa_m
  
  ! Local variables.
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
!    write(*,*) "Triangle ", t, " vertices (IDs [", tri_vert(1:3, t),"]): "
    do v = 1, 3
      local_vertice_list(v, 1:NDIM) = X_m(1:NDIM, tri_vert(v, t))
!      write(*,*) "                           ", local_vertice_list(v, 1:NDIM)
    enddo
    
!    if(point_is_in_triangle(local_vertice_list, point_to_test)) then
!      found = .true.
!      write(*,*) "Point (",point_to_test,") is in triangle n°", t,". Performing barycentric interpolation."
!      container_triangle = t
!      call barycentric_coordinates_2d(local_vertice_list, point_to_test, barycor)
!      ! https://codeplea.com/triangular-interpolation
!      return
!    endif
    
    call barycentric_coordinates_2d(local_vertice_list, point_to_test, inTri, barycor)
    
    if(inTri) then
      ! Point is inside triangle.
      loctri_vertices_ids = tri_vert(1:3, t)
!      write(*,*) "Point (",coord(1:NDIM, ibool(i, j, ispec)),") is     in tri. n°", t," (barcoor. = [", barycor,"])."
!      write(*,*) "Performing barycentric interpolation." ! https://codeplea.com/triangular-interpolation
!      write(*,*) "Vertices of this triangle are model points n°[", loctri_vertices_ids, "]."
!      write(*,*) "Positions of those vertices are: [", local_vertice_list(1, 1:NDIM), "],"
!      write(*,*) "                                 [", local_vertice_list(2, 1:NDIM), "],"
!      write(*,*) "                                 [", local_vertice_list(3, 1:NDIM), "]."
!      write(*,*) "Values of p_model at those points are: [", p_m(loctri_vertices_ids), "]."
!      write(*,*) "Weighted value is (p_model[v1 v2 v3].[barycor(1) barycor(2) barycor(3)]) = ", DOT_PRODUCT(p_m(loctri_vertices_ids), barycor), "."
      LNS_rho0(iglobDG)    = DOT_PRODUCT(barycor, rho_m(loctri_vertices_ids))
      do v = 1, NDIM
        LNS_v0(v, iglobDG) = DOT_PRODUCT(barycor, v_m(v, loctri_vertices_ids))
      enddo
      LNS_p0(iglobDG)      = DOT_PRODUCT(barycor, p_m(loctri_vertices_ids))
      LNS_g(iglobDG)       = DOT_PRODUCT(barycor, g_m(loctri_vertices_ids))
      gammaext_DG(iglobDG) = DOT_PRODUCT(barycor, gam_m(loctri_vertices_ids))
      LNS_mu(iglobDG)      = DOT_PRODUCT(barycor, mu_m(loctri_vertices_ids))
      LNS_kappa(iglobDG)   = DOT_PRODUCT(barycor, kappa_m(loctri_vertices_ids))
      !vpext(i, j, ispec) = vp_model(1)
      !Nsqext(i, j, ispec) = Nsq_model(1)
      !vsext(i, j, ispec) = ZERO
      !Qmu_attenuationext(i, j, ispec) = HUGEVAL
      !QKappa_attenuationext(i, j, ispec) = HUGEVAL
      !Htabext_DG(indglob_DG) = Htab_model(1)
      !tau_sigma(i, j, ispec) = tau_sigma_model(i)
      !tau_epsilon(i, j, ispec) = tau_epsilon_model(i)
      
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
        stop
      endif
      
      ! For the subroutines in 'invert_mass_matrix.f90', one needs to initialise the following:
      rhoext(i, j, ispec) = LNS_rho0(iglobDG)
      LNS_c0(iglobDG) = sqrt(gammaext_DG(iglobDG)*LNS_p0(iglobDG)/LNS_rho0(iglobDG)) ! Take the chance to compute and save c0.
      vpext(i, j, ispec) = LNS_c0(iglobDG)
      
      return
    else
      ! Point is outside triangle, go to next triangle.
!      write(*,*) "Point (",point_to_test,") is not in tri. n°", t," (barcoor. = [", barycor,"])."
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
    stop
  endif
  
end subroutine delaunay_interpolate_one_point


! ------------------------------------------------------------ !
! barycentric_coordinates_2d                                   !
! ------------------------------------------------------------ !
! https://en.wikipedia.org/wiki/Barycentric_coordinate_system

subroutine barycentric_coordinates_2d(vlist, P, inTri, barycor)
  use constants, only: CUSTOM_REAL, NDIM, TINYVAL
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
    stop
  endif
  barycor(1) = ((vlist(2, 2)-vlist(3, 2))*(P(1)-vlist(3, 1))+(vlist(3, 1)-vlist(2, 1))*(P(2)-vlist(3, 2))) / det
  if(barycor(1)<-TINYVAL) return ! Point is outside triangle, stop and return.
  barycor(2) = ((vlist(3, 2)-vlist(1, 2))*(P(1)-vlist(3, 1))+(vlist(1, 1)-vlist(3, 1))*(P(2)-vlist(3, 2))) / det
  if(barycor(2)<-TINYVAL) return ! Point is outside triangle, stop and return.
  barycor(3) = 1. - sum(barycor(1:2))
  if(barycor(3)<-TINYVAL) return ! Point is outside triangle, stop and return.
  inTri = .true. ! No return command was hit before, hence point can be considered in triangle.
  !write(*,*) "Barycentric coordinates = (", l1, l2, l3, ")."
end subroutine barycentric_coordinates_2d


! ------------------------------------------------------------ !
! output_lns_interpolated_model                                !
! ------------------------------------------------------------ !
! Outputs the loaded model in full to a file, for checking/plotting purposes.

subroutine output_lns_interpolated_model()
  use constants, only: CUSTOM_REAL, NDIM, NGLLX, NGLLZ
  use specfem_par, only: myrank, nspec, ibool_before_perio, coord, ibool_DG, gammaext_DG
  use specfem_par_lns, only: LNS_rho0, LNS_v0, LNS_p0, LNS_g, LNS_mu, LNS_kappa
  implicit none
  ! Local variables.
  integer :: ispec, i, j, iglobDG
  ! Let all procs write to file (I know this may not be sequential and freak up the file, but I don't care since it's only for debugging purposes).
  open(unit=504,file='OUTPUT_FILES/LNS_GENERAL_INTERPOLATED_MODEL',status='unknown',action='write', position="append")
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
    write(*,*) "> > Dumped interpolated model (on proc 0) to './OUTPUT_FILES/LNS_GENERAL_INTERPOLATED_MODEL'. ", &
               "Use the following Matlab one-liner to plot:"
    write(*,*) "      OFD=''; [~,lab]=order_bg_model();a=importdata(OFD);x=a(:,1);z=a(:,2);close all;for i=3:10;", &
               "d=[];d.rho=[];d.vel=[];d.pre=a(:,i);[Xi,Yi,Zi]=interpDumps(x,z,d,0,0,1);figure(i);surf(Xi,Yi,Zi.pre);", &
               "view([0,0,1]);shading flat;colorbar;title(lab{i});end;"
  endif
end subroutine output_lns_interpolated_model

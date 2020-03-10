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
  use specfem_par, only: myrank, USE_DISCONTINUOUS_METHOD
  
  implicit none
  
  ! Input/output.
  integer, intent(in) :: nlines_header, nlines_model
  
  ! Local variables.
  real(kind=CUSTOM_REAL), parameter :: ZERO = 0._CUSTOM_REAL
  real(kind=CUSTOM_REAL), dimension(NDIM, nlines_model) :: X_m
  real(kind=CUSTOM_REAL), dimension(nlines_model) :: p_model
  
  ! Safeguard.
  if(.not. USE_DISCONTINUOUS_METHOD) then
    write(*,*) "********************************"
    write(*,*) "*            ERROR             *"
    write(*,*) "********************************"
    write(*,*) "* Currently cannot use the     *"
    write(*,*) "* 'external_DG' model if not   *"
    write(*,*) "* using the DG method.         *"
    write(*,*) "********************************"
    stop
  endif
  
  ! Read and store values of model.
  call lns_read_background_model(nlines_header, nlines_model, X_m, p_model)
  
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
                                  p_model &
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

subroutine lns_read_background_model(nlines_header, nlines_model, X_m, p_model)
  use constants, only: CUSTOM_REAL, NDIM
  !use specfem_par,only: nspec, tau_sigma, tau_epsilon
  use specfem_par, only: myrank
  use specfem_par_lns, only: BCKGRD_MDL_LNS_FILENAME, BCKGRD_MDL_LNS_NCOL
  
  implicit none
  
  ! Input/output.
  integer, intent(in) :: nlines_header, nlines_model
  real(kind=CUSTOM_REAL), dimension(NDIM, nlines_model), intent(out) :: X_m
  real(kind=CUSTOM_REAL), dimension(nlines_model), intent(out) :: p_model
  
  ! Local variables.
  integer :: ncolumns_detected = 3
  integer :: i, io
  
  call external_model_DG_testfile(BCKGRD_MDL_LNS_FILENAME) ! Test-reading the file.
  
  call external_model_DG_only_find_nbcols(BCKGRD_MDL_LNS_FILENAME, nlines_header, ncolumns_detected)
  
  if(myrank==0) then
    write(*,*) "> Reading general background model file '", trim(BCKGRD_MDL_LNS_FILENAME)
    write(*,*) "> Detected ",nlines_header," lines for header, and ",nlines_model," lines of model ", &
               "(total = ",(nlines_header+nlines_model)," lines). Detected ", ncolumns_detected," columns."
  endif
  
  OPEN(100, file=BCKGRD_MDL_LNS_FILENAME)
  
  do i=1, nlines_header
    ! Read and skip in header.
    read(100, *, iostat=io)
    IF (io/=0) stop "Error reading line in background model file."
  enddo
  
  do i=1, nlines_model
    ! Read values.
    if(ncolumns_detected==BCKGRD_MDL_LNS_NCOL) then
      read(100, *, iostat=io) X_m(1, i), X_m(NDIM, i), p_model(i)
      !write(*,*) X_m(1, i), X_m(NDIM, i), p_model(i)
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
  close(100)
end subroutine lns_read_background_model

! ------------------------------------------------------------ !
! delaunay_background_model_2d                                 !
! ------------------------------------------------------------ !


subroutine delaunay_background_model_2d(nlines_model, X_m, tri_num, tri_vert, tri_nabe)
  use constants, only: CUSTOM_REAL, NDIM
  
  implicit none
  ! Input/output.
  integer, intent(in) :: nlines_model
  real(kind=8), dimension(NDIM, nlines_model), intent(in) :: X_m
  integer(kind=4), intent(out) :: tri_num
  integer(kind=4), dimension(3, nlines_model), intent(out) :: tri_vert, tri_nabe
  ! Local variables.
  integer :: t
  call dtris2(nlines_model, X_m, tri_num, tri_vert, tri_nabe)
  
  write(*,*) "> > Delaunay triangulation: ", tri_num, "triangles found."
  
  ! Debug: print out triangles.
  if(.false.) then
    do t = 1, tri_num
      write(*,*) "> > Triangle n°", t, " vertices: ", X_m(1:NDIM, tri_vert(1, t))
      write(*,*) "> >                                   ", X_m(1:NDIM, tri_vert(2, t))
      write(*,*) "> >                                   ", X_m(1:NDIM, tri_vert(3, t))
    enddo
  endif
  
end subroutine delaunay_background_model_2d

! ------------------------------------------------------------ !
! delaunay_interp_all_points                                   !
! ------------------------------------------------------------ !


subroutine delaunay_interp_all_points(nlines_model, X_m, p_m)
  use constants, only: CUSTOM_REAL, NDIM, NGLLX, NGLLZ
  use specfem_par_lns, only: point_is_in_triangle
  use specfem_par,only: nspec, &
        ispec_is_elastic, ispec_is_acoustic_DG
  
  implicit none
  
  ! Input/output.
  integer, intent(in) :: nlines_model
  real(kind=8), dimension(NDIM, nlines_model), intent(in) :: X_m
  real(kind=CUSTOM_REAL), dimension(nlines_model), intent(in) :: p_m
  
  ! Local variables.
  integer :: ispec, i, j
  integer(kind=4) :: tri_num
  integer(kind=4), dimension(3, nlines_model) :: tri_vert, tri_nabe
  
  call delaunay_background_model_2d(nlines_model, X_m, tri_num, tri_vert, tri_nabe)
  
  do ispec = 1, nspec
    if(ispec_is_elastic(ispec)) then
      cycle ! Elastic regions are treated separately (see call to 'external_DG_update_elastic_from_parfile' below).
    else if(ispec_is_acoustic_DG(ispec)) then
        ! For DG elements, go through GLL points one by one.
        do j = 1, NGLLZ
          do i = 1, NGLLX
            call delaunay_interpolate_one_point(nlines_model, X_m, tri_num, tri_vert, ispec, i, j, &
                                                p_m &
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
end subroutine delaunay_interp_all_points


! ------------------------------------------------------------ !
! delaunay_interpolate_one_point                               !
! ------------------------------------------------------------ !


subroutine delaunay_interpolate_one_point(nlines_model, X_m, tri_num, tri_vert, ispec, i, j, &
                                          p_m &
                                         )
  use constants, only: CUSTOM_REAL, NDIM
  use specfem_par_lns, only: point_is_in_triangle
  use specfem_par, only: ibool, coord, pext_DG, coord_interface
  
  implicit none
  
  ! Input/output.
  integer, intent(in) :: nlines_model, ispec, i, j, tri_num
  real(kind=8), dimension(NDIM, nlines_model), intent(in) :: X_m
  integer(kind=4), dimension(3, nlines_model), intent(in) :: tri_vert
  real(kind=CUSTOM_REAL), dimension(nlines_model), intent(in) :: p_m
  
  ! Local variables.
  integer :: t, v
  real(kind=8), dimension(3, NDIM) :: local_vertice_list
  real(kind=CUSTOM_REAL), dimension(NDIM) :: point_to_test
  logical found
  real(kind=CUSTOM_REAL) :: l1, l2, l3
  
  ! Bring back the SPECFEM coordinate to the model convention (where the interface is at altitude zero).
  point_to_test = coord(1:NDIM, ibool(i, j, ispec)) - (/0._CUSTOM_REAL, coord_interface/)
  
  found = .false.
  
  do t = 1, tri_num
    !write(*,*) "Triangle ", t, " vertices: "
    do v = 1, 3
      local_vertice_list(v, 1:NDIM) = X_m(1:NDIM, tri_vert(v, t))
      !write(*,*) "                           ", local_vertice_list(v, 1:NDIM)
    enddo
    
!    if(point_is_in_triangle(local_vertice_list, point_to_test)) then
!      found = .true.
!      write(*,*) "Point (",point_to_test,") is in triangle n°", t,". Performing barycentric interpolation."
!      container_triangle = t
!      call barycentric_coordinates_2d(local_vertice_list, point_to_test, l1, l2, l3)
!      ! https://codeplea.com/triangular-interpolation
!      return
!    endif
    
    call barycentric_coordinates_2d(local_vertice_list, point_to_test, l1, l2, l3)
    
    if(l1>=0. .and. l2>=0. .and. l3>=0.) then
      ! Point is inside triangle.
!      write(*,*) "Point (",coord(1:NDIM, ibool(i, j, ispec)),") is     in tri. n°", t," (barcoor. = [", l1,", ", l2,", ", l3,"])."
!      write(*,*) "Performing barycentric interpolation." ! https://codeplea.com/triangular-interpolation
!      write(*,*) "Vertices of this triangle are model points n°[", tri_vert(1:3, t), "]."
!      write(*,*) "Positions of those vertices are: [", local_vertice_list(1, 1:NDIM), "],"
!      write(*,*) "                                 [", local_vertice_list(2, 1:NDIM), "],"
!      write(*,*) "                                 [", local_vertice_list(3, 1:NDIM), "]."
!      write(*,*) "Values of p_model at those points are: [", p_m(tri_vert(1:3, t)), "]."
!      write(*,*) "Weighted value is (p_model[v1 v2 v3].[l1 l2 l3]) = ", DOT_PRODUCT(p_m(tri_vert(1:3, t)), (/l1, l2, l3/)), "."
      pext_DG(i, j, ispec) = DOT_PRODUCT(p_m(tri_vert(1:3, t)), (/l1, l2, l3/))
      
      return
    else
      ! Point is outside triangle, go to next triangle.
!      write(*,*) "Point (",point_to_test,") is not in tri. n°", t," (barcoor. =", &
!                 " [", l1,", ", l2,", ", l3,"])."
      cycle
    endif
  enddo
  
  if(.not. found) then
    write(*,*) "********************************"
    write(*,*) "*            ERROR             *"
    write(*,*) "********************************"
    write(*,*) "* Point                        *"
    write(*,*) "* ", point_to_test
    write(*,*) "* was not found within the     *"
    write(*,*) "* Delaunay triangulation.      *"
    write(*,*) "* You can look at the          *"
    write(*,*) "* triangles by commenting out  *"
    write(*,*) "* the 'if(.false.)' in the     *"
    write(*,*) "* subroutine                   *"
    write(*,*) "* 'delaunay_background_model_2d'"
    write(*,*) "* (in                          *"
    write(*,*) "* 'lns_load_background_model.f90')."
    write(*,*) "* You can also make sure the   *"
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

subroutine barycentric_coordinates_2d(vlist, P, l1, l2, l3)
  use constants, only: CUSTOM_REAL, NDIM
  implicit none
  ! Input/output.
  real(kind=8), dimension(3, NDIM), intent(in) :: vlist
  real(kind=CUSTOM_REAL), dimension(NDIM), intent(in) :: P
  real(kind=CUSTOM_REAL), intent(out) :: l1, l2, l3
  ! Local variables.
  real(kind=CUSTOM_REAL), parameter :: ZERO = 0._CUSTOM_REAL
  real(kind=CUSTOM_REAL) :: det
  l1 = ZERO
  l2 = ZERO
  l3 = ZERO
  det = (vlist(1, 1)-vlist(3, 1))*(vlist(2, 2)-vlist(3, 2)) - (vlist(1, 2)-vlist(3, 2))*(vlist(2, 1)-vlist(3, 1))
  l1 = ((vlist(2, 2)-vlist(3, 2))*(P(1)-vlist(3, 1))+(vlist(3, 1)-vlist(2, 1))*(P(2)-vlist(3, 2))) / det
  if(l1<0.) return ! Point is outside triangle, stop and return.
  l2 = ((vlist(3, 2)-vlist(1, 2))*(P(1)-vlist(3, 1))+(vlist(1, 1)-vlist(3, 1))*(P(2)-vlist(3, 2))) / det
  if(l2<0.) return ! Point is outside triangle, stop and return.
  l3 = 1. - l1 - l2
  if(l3<0.) return ! Point is outside triangle, stop and return.
  !write(*,*) "Barycentric coordinates = (", l1, l2, l3, ")."
end subroutine barycentric_coordinates_2d


! ------------------------------------------------------------ !
! output_lns_interpolated_model                                !
! ------------------------------------------------------------ !
! Outputs the loaded model in full to a file, for checking/plotting purposes.
! Matlab one-liner plot:
! a=importdata('/home/l.martire/Documents/SPECFEM/specfem-dg-master/EXAMPLES/test_lns_custom_wavefield/OUTPUT_FILES/TESTMODEL'); x=a(:,1); z=a(:,2); d=a(:,3); [Xi, Yi, Zi] = interpDumps(x, z, d, 50, 50); surf(Xi, Yi, Zi); view([0,0,1]); shading flat; colorbar;

subroutine output_lns_interpolated_model()
  use constants, only: CUSTOM_REAL, NDIM, NGLLX, NGLLZ
  use specfem_par, only: myrank, nspec, ibool_before_perio, coord, &
                         pext_DG
  implicit none
  ! Local variables.
  integer :: ispec, i, j
  open(unit=504,file='OUTPUT_FILES/TESTMODEL',status='unknown',action='write', position="append")
  do ispec = 1,nspec
    do j = 1,NGLLZ
      do i = 1,NGLLX
        write(504,*) coord(1, ibool_before_perio(i, j, ispec)), coord(2, ibool_before_perio(i, j, ispec)),&
                     pext_DG(i, j, ispec)
      enddo
    enddo
  enddo
  close(504)
  if(myrank==0) then
    write(*,*) "> > Dumped interpolated model to './OUTPUT_FILES/TESTMODEL'. Use the following Matlab one-liner to plot:"
    write(*,*) "      a=importdata('/home/l.martire/Documents/SPECFEM/specfem-dg-master/EXAMPLES/test_lns_custom_wavefield/", &
               "OUTPUT_FILES/TESTMODEL'); x=a(:,1); z=a(:,2); d=a(:,3); [Xi, Yi, Zi] = interpDumps(x, z, d, 50, 50); ", &
               "surf(Xi, Yi, Zi); view([0,0,1]); shading flat; colorbar;"
  endif
end subroutine output_lns_interpolated_model

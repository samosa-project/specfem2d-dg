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

  use constants,only: CUSTOM_REAL, NDIM, NGLLX,NGLLZ,NDIM,IMAIN,FOUR_THIRDS,TINYVAL,PI
  use specfem_par,only: nspec, ibool, ibool_DG, coord, coord_interface, &
        ispec_is_elastic, ispec_is_acoustic_DG, &
        pext_DG
  
  implicit none
  
  ! Input/output.
  integer, intent(in) :: nlines_header, nlines_model
  
  ! Local variables.
  real(kind=CUSTOM_REAL), parameter :: ZERO = 0._CUSTOM_REAL
  real(kind=CUSTOM_REAL), dimension(NDIM, nlines_model) :: X_m
  real(kind=CUSTOM_REAL), dimension(nlines_model) :: p_model
  
  integer :: i, j, ispec, ii, indglob_DG
  double precision :: z, frac, pii, piim1, piim2!, piip1!,tmp1, gamma_temp,gamma_temp_prev,x,max_z
  
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
  call delaunay_interp_all_points(nlines_model, X_m)
  
  ! Interpolate.
  do ispec = 1, nspec
    if(ispec_is_elastic(ispec)) then
      cycle ! Elastic regions are treated separately (see call to 'external_DG_update_elastic_from_parfile' below).
      
    else if(ispec_is_acoustic_DG(ispec)) then
      ! For DG elements, go through GLL points one by one.
      do j = 1, NGLLZ
        do i = 1, NGLLX
          indglob_DG = ibool_DG(i, j, ispec)
          
          z = coord(2, ibool(i, j, ispec)) ! Get z-coordinate of GLL point.
          z = z - coord_interface ! Get altitude of GLL point.
          
          ! Find first index ii (from bottom) in model such that altitude_model(ii) > z.
          ii = 1
          do while(z >= X_m(NDIM, ii) .and. ii /= nlines_model)
            write(*,*) z, X_m(NDIM, ii), ii
            ii = ii + 1
          enddo
          write(*,*) z, X_m(NDIM, ii), ii
          
          if (ii == 1) then
            ! Altitude of point is <= the altitude of first line of model.
            ! Use first line of model.
            
            pext_DG(i, j, ispec) = p_model(1)
            
          elseif (ii == 2) then
            ! Altitude of point is > 1st line of model, and <= 2nd line of model.
            ! Interpolate using the values at lines 1 (ii-1) and 2 (ii) of model.
            frac = (z-X_m(NDIM, ii-1))/(X_m(NDIM, ii)-X_m(NDIM, ii-1))
            
            pext_DG(i, j, ispec)          = p_model(ii-1)     + frac*(p_model(ii)-p_model(ii-1))

          else
            ! Altitude of point is > (ii-1)-th line of model, and <= (ii)-th line of model.
            ! Interpolate using the values at lines (ii-2), (ii-1), and ii of model.
            pii   = (z-X_m(NDIM, ii-1))*(z-X_m(NDIM, ii-2))/((X_m(NDIM, ii)-X_m(NDIM, ii-1))*(X_m(NDIM, ii)-X_m(NDIM, ii-2)))
            piim1 = (z-X_m(NDIM, ii))  *(z-X_m(NDIM, ii-2))/((X_m(NDIM, ii-1)-X_m(NDIM, ii))*(X_m(NDIM, ii-1)-X_m(NDIM, ii-2)))
            piim2 = (z-X_m(NDIM, ii-1))*(z-X_m(NDIM, ii))  /((X_m(NDIM, ii-2)-X_m(NDIM, ii))*(X_m(NDIM, ii-2)-X_m(NDIM, ii-1)))
            
            pext_DG(i, j, ispec)          = p_model(ii)*pii     + p_model(ii-1)*piim1 + p_model(ii-2)*piim2
          endif ! Endif on ii.
        enddo ! Enddo on i.
      enddo ! Enddo on j.
    else
      write(*,*) "********************************"
      write(*,*) "*            ERROR             *"
      write(*,*) "********************************"
      write(*,*) "* lns_load_background_model.f90"
      write(*,*) "********************************"
      stop
    endif ! Endif on ispec_is_acoustic_DG.
  enddo ! Enddo on ispec.
  
  call external_DG_update_elastic_from_parfile() ! Update elastic regions by reading parameters directly from parfile (see 'define_external_model.f90').
  
end subroutine lns_load_background_model


! ------------------------------------------------------------ !
! lns_read_background_model                                    !
! ------------------------------------------------------------ !
! Read and store values of model.

subroutine lns_read_background_model(nlines_header, nlines_model, X_m, p_model)
  use constants, only: CUSTOM_REAL, NDIM
  !use specfem_par,only: nspec, tau_sigma, tau_epsilon
  
  implicit none
  
  ! Input/output.
  integer, intent(in) :: nlines_header, nlines_model
  real(kind=CUSTOM_REAL), dimension(NDIM, nlines_model), intent(out) :: X_m
  real(kind=CUSTOM_REAL), dimension(nlines_model), intent(out) :: p_model
  
  ! Local variables.
  integer, parameter :: ncolz = 3
  character(len=100), parameter :: BCKGRD_MDL_LNS_FILENAME = './background_model.dat'
  integer :: i, io
  
  !write(*,*) "Opening '", trim(BCKGRD_MDL_LNS_FILENAME), "'."
  OPEN(100, file=BCKGRD_MDL_LNS_FILENAME)
  
  do i=1, nlines_header
    ! Read and skip in header.
    read(100, *, iostat=io)
    IF (io/=0) stop "Error reading line in background model file."
  enddo
  
  do i=1, nlines_model
    ! Read values.
    if(ncolz==3) then
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

subroutine delaunay_background_model_2d(nlines_model, X_m, tri_num, tri_vert, tri_nabe)
  use constants, only: CUSTOM_REAL, NDIM
  
  implicit none
  
  integer, intent(in) :: nlines_model
  real(kind=8), dimension(NDIM, nlines_model), intent(in) :: X_m
  integer(kind=4), intent(out) :: tri_num
  integer(kind=4), dimension(3, nlines_model), intent(out) :: tri_vert, tri_nabe
  call dtris2(nlines_model, X_m, tri_num, tri_vert, tri_nabe)
  
  write(*,*) "Delaunay triangulation: ", tri_num, "triangles found."
end subroutine delaunay_background_model_2d

subroutine delaunay_interp_all_points(nlines_model, X_m)
  use constants, only: CUSTOM_REAL, NDIM
  use specfem_par_lns, only: point_is_in_triangle
  
  implicit none
  
  ! Input/output.
  integer, intent(in) :: nlines_model
  real(kind=8), dimension(NDIM, nlines_model), intent(in) :: X_m
  
  ! l
  integer :: container_triangle
  integer(kind=4) :: tri_num
  integer(kind=4), dimension(3, nlines_model) :: tri_vert, tri_nabe
  real(kind=CUSTOM_REAL), dimension(2), PARAMETER :: point_to_test = (/0.25, 0.25/)
  
  call delaunay_background_model_2d(nlines_model, X_m, tri_num, tri_vert, tri_nabe)
  
!  do ispec = 1, nspec
!    do i = 
!      do j = 
        call delaunay_check_one_point(nlines_model, X_m, tri_num, tri_vert, point_to_test, container_triangle)
!      enddo
!    enddo
!  enddo
end subroutine delaunay_interp_all_points

subroutine delaunay_check_one_point(nlines_model, X_m, tri_num, tri_vert, point_to_test, container_triangle)
  use constants, only: CUSTOM_REAL, NDIM
  use specfem_par_lns, only: point_is_in_triangle
  
  implicit none
  
  ! Input/output.
  integer, intent(in) :: nlines_model
  real(kind=8), dimension(NDIM, nlines_model), intent(in) :: X_m
  integer(kind=4), intent(in) :: tri_num
  integer(kind=4), dimension(3, nlines_model), intent(in) :: tri_vert
  real(kind=CUSTOM_REAL), dimension(2), intent(in) :: point_to_test
  integer, intent(out) :: container_triangle
  
  ! Local variables.
  integer :: t, v
  real(kind=8), dimension(3, NDIM) :: local_vertice_list
  logical found
  
  found = .false.
  do t = 1, tri_num
    !write(*,*) "Triangle ", t, " vertices: "
    do v = 1, 3
      local_vertice_list(v, 1:NDIM) = X_m(1:NDIM, tri_vert(v, t))
      !write(*,*) "                           ", local_vertice_list(v, 1:NDIM)
    enddo
    if(point_is_in_triangle(local_vertice_list, point_to_test)) then
      found = .true.
      write(*,*) "Point (",point_to_test,") is in triangle n°", t,". Performing barycentric interpolation."
      container_triangle = t
      call barycentric_coordinates_2d(local_vertice_list, point_to_test)
      ! https://codeplea.com/triangular-interpolation
      return
    endif
  enddo
  
  if(.not. found) then
    write(*,*) "********************************"
    write(*,*) "*            ERROR             *"
    write(*,*) "********************************"
    write(*,*) "* Point was not found within   *"
    write(*,*) "* the convex hull of points    *"
    write(*,*) "* specified in the background  *"
    write(*,*) "* model file.                  *"
    write(*,*) "********************************"
    stop
  endif
  
end subroutine delaunay_check_one_point

subroutine barycentric_coordinates_2d(vertice_list, P)
  ! https://en.wikipedia.org/wiki/Barycentric_coordinate_system
  use constants, only: CUSTOM_REAL, NDIM
  
  implicit none
  
  !i
  real(kind=8), dimension(3, NDIM), intent(in) :: vertice_list
  real(kind=CUSTOM_REAL), dimension(NDIM), intent(in) :: P  
  
  !l
  real(kind=CUSTOM_REAL) :: x1, x2, x3, y1, y2, y3, x, y, det, l1, l2, l3
  real(kind=CUSTOM_REAL), dimension(2, 2) :: T
  
  x1 = vertice_list(1, 1)
  x2 = vertice_list(2, 1)
  x3 = vertice_list(3, 1)
  y1 = vertice_list(1, 2)
  y2 = vertice_list(2, 2)
  y3 = vertice_list(3, 2)
  x = P(1)
  y = P(2)
  
  T(1, 1) = x1-x3
  T(1, 2) = x2-x3
  T(2, 1) = y1-y3
  T(2, 2) = y2-y3
  det = T(1, 1)*T(2, 2) - T(2, 1)*T(1, 2)
  
  l1 = ((y2-y3)*(x-x3)+(x3-x2)*(y-y3)) / det
  l2 = ((y3-y1)*(x-x3)+(x1-x3)*(y-y3)) / det
  l3 = 1. - l1 - l2
  
  write(*,*) "Barycentric coordinates = (", l1, l2, l3, ")."
  
end subroutine barycentric_coordinates_2d

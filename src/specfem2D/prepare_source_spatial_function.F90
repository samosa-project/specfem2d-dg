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
! prepare_timerun_ssf                                          !
! ------------------------------------------------------------ !
! Wrapper. Allocates and prepares a vector containing the values of the source spatial function.

subroutine prepare_timerun_ssf()
  
  use constants, only: ZERO
  use specfem_par

  implicit none

  if(initialfield) then
    ! One uses an initialfield. Thus, do a dummy allocation.
    allocate(source_spatial_function_DG(1, 1))
  else
    ! One does not use an initial field.
    if(.not. any_acoustic_DG) then
      ! No DG elements exist, ignore call and inform user.
      if(myrank == 0) then
        write(*,*) "********************************"
        write(*,*) "*           WARNING            *"
        write(*,*) "********************************"
        write(*,*) "* No DG elements exist on      *"
        write(*,*) "* CPU ", myrank
        write(*,*) "* The spread source spatial    *"
        write(*,*) "* function(s) cannot be        *"
        write(*,*) "* implemented here. The        *"
        write(*,*) "* classical spatial source     *"
        write(*,*) "* function(s) will be used. If *"
        write(*,*) "* the source is elsewhere (in  *"
        write(*,*) "* another CPU), this poses no  *"
        write(*,*) "* problem.                     *"
        write(*,*) "********************************"
        call flush_IMAIN()
      endif
    else
      ! DG elements exist.
      if(myrank == 0) then
        write(*,*) "********************************"
        write(*,*) "*           WARNING            *"
        write(*,*) "********************************"
        write(*,*) "* This kind of source has to   *"
        write(*,*) "* be tuned not to generate     *"
        write(*,*) "* spurious discontinuities.    *"
        write(*,*) "* Use at your own risk.        *"
        write(*,*) "********************************"
        call flush_IMAIN()
      endif
      ! TODO: Implement a less cumbersome way of storing the values. Indeed, many values far away from the source will be negligible in this array.
      ! For example, instead of storing values over all the points, build the list of indices which associated points have a non-negligible value (value greater than a given threshold), and store only the values of the source spatial function at those points.
      allocate(source_spatial_function_DG(NSOURCES, nglob))
      source_spatial_function_DG(:, :) = ZERO
      call prepare_source_spatial_function_DG() ! Compute the SSF array.
    endif ! Endif on any_acoustic_DG.
  endif ! Endif on initialfield.
  
  call synchronize_all()

end subroutine prepare_timerun_ssf

! ------------------------------------------------------------ !
! prepare_source_spatial_function_DG                           !
! ------------------------------------------------------------ !
! Compute values of the source spatial function at all points.

subroutine prepare_source_spatial_function_DG
  use constants, only: CUSTOM_REAL, NGLLX, NGLLZ, TINYVAL
  use specfem_par, only: coord, ibool_before_perio, IMAIN, ispec_is_acoustic_DG, &
                         myrank, NSOURCES, nspec, SPREAD_SSF_CUSTOM, &!ispec_selected_source, &
                         source_spatial_function_DG, source_type, SPREAD_SSF_SAVE, SPREAD_SSF_SIGMA, x_source, z_source 
  implicit none
  
  ! Local variables.
  real(kind=CUSTOM_REAL), parameter :: ZEROcr = 0._CUSTOM_REAL
  integer :: i_source, iglob_unique, ispec, i, j
  real(kind=CUSTOM_REAL) :: distsqrd
  character(len=33) :: filename ! Used for saving values. Length: 16 for "OUTPUT_FILES/SSF" + 8 for process numbering + 1 + 8 for source numbering.
  
  do i_source = 1, NSOURCES ! Loop on sources.
    if(SPREAD_SSF_SAVE) then
      write(filename, '( "OUTPUT_FILES/SSF", i8.8, "_", i8.8 )' ) i_source, myrank
      open(unit=504,file=filename,status='unknown',action='write', position="append")
    endif
    
    ! 190828: removed this warning. Indeed, only points within DG elements are updated with the SSF. So, this poses no problem when it overlaps with non-DG elements. Worse, it causes some with large sources when procs share DG and non-DG elements.
    !write(*,*) myrank, 'KEK', coord(:, ibool_before_perio(3, 3, ispec_selected_source(i_source))) ! DEBUG
    !if(.not. ispec_is_acoustic_DG(ispec_selected_source(i_source))) then
    !  ! The central point of the source is not in a DG element, ignore call and inform user.
    !  if(myrank == 0) then
    !    write(*,*) "********************************"
    !    write(*,*) "*           WARNING            *"
    !    write(*,*) "********************************"
    !    write(*,*) "* The central point of the     *"
    !    write(*,*) "* source                       *"
    !    write(*,*) "* ", i_source
    !    write(*,*) "* is not in a DG element, a    *"
    !    write(*,*) "* spread source spatial        *"
    !    write(*,*) "* function cannot be           *"
    !    write(*,*) "* implemented. A classical     *"
    !    write(*,*) "* spatial source function will *"
    !    write(*,*) "* be used.                     *"
    !    write(*,*) "********************************"
    !    call flush_IMAIN()
    !  endif
    !else
      ! The central point of the source is indeed in a DG element.
      ! TODO: More cases may exist where the use of a spread source spatial function has to be discarded. Identify those and produce relevant tests and error messages.
      
      !if(source_type(i_source)==2 .and. myrank == 0) then
      !  write(*,*) "********************************"
      !  write(*,*) "*           WARNING            *"
      !  write(*,*) "********************************"
      !  write(*,*) "* You are using source_type=2  *"
      !  write(*,*) "* within a DG simulation and   *"
      !  write(*,*) "* where SPREAD_SSF is          *"
      !  write(*,*) "* activated. This means a      *"
      !  write(*,*) "* hardcoded SSF will be used.  *"
      !  write(*,*) "* Be sure of what you are      *"
      !  write(*,*) "* doing.                       *"
      !  write(*,*) "********************************"
      !  call flush_IMAIN()
      !endif
      
      do ispec = 1, nspec
        if(ispec_is_acoustic_DG(ispec)) then
          do i = 1, NGLLX
            do j = 1, NGLLZ
              iglob_unique = ibool_before_perio(i, j, ispec)
              select case (source_type(i_source))
                case (1)
                  ! If the source is an elastic force or an acoustic pressure.
                  select case(SPREAD_SSF_CUSTOM)
                    case(0)
                      ! Classical Gaussian SSF.
                      distsqrd =   (coord(1, iglob_unique) - x_source(i_source))**2. &
                                 + (coord(2, iglob_unique) - z_source(i_source))**2.
                      source_spatial_function_DG(i_source, iglob_unique) = exp(-distsqrd/(SPREAD_SSF_SIGMA**2.))
                    case(1)
                      ! Plane wave.
                      call CUSTOM_SSF_PLANE_WAVE(coord(:,iglob_unique), source_spatial_function_DG(i_source, iglob_unique))
                    case(99)
                      call CUSTOM_SSF_MOUNTAIN_WAVE(coord(:,iglob_unique), (/x_source(i_source), z_source(i_source)/), &
                                                    source_spatial_function_DG(i_source, iglob_unique))
                  end select ! Endselect on SPREAD_SSF_CUSTOM.
                  
                case default
                  if(myrank == 0) then
                    write(*,*) "********************************"
                    write(*,*) "*            ERROR             *"
                    write(*,*) "********************************"
                    write(*,*) "* This source type as spread   *"
                    write(*,*) "* source spatial function is   *"
                    write(*,*) "* not implemented.             *"
                    write(*,*) "* i_source              = ", i_source
                    write(*,*) "* source_type(i_source) = ", source_type(i_source)
                    write(*,*) "********************************"
                    write(*,*) "* If your source is in an      *"
                    write(*,*) "* elastic media, simply        *"
                    write(*,*) "* disable USE_SPREAD_SSF in    *"
                    write(*,*) "* the parfile.                 *"
                    write(*,*) "* If your source is in an      *"
                    write(*,*) "* acoustic media, it should    *"
                    write(*,*) "* probably not have this type. *"
                    write(*,*) "********************************"
                    stop
                  endif
              end select ! Endselect on source_type.
              
              if(SPREAD_SSF_SAVE) then
                if(source_spatial_function_DG(i_source, iglob_unique)>TINYVAL) then
                  ! Do not save zeros, because that would be silly.
                  ! Do not save negligible values either.
                  write(504,*) coord(1, iglob_unique), coord(2, iglob_unique), source_spatial_function_DG(i_source, iglob_unique)
                endif
              endif
            enddo ! Enddo on j.
          enddo ! Enddo on i.
        endif ! Endif on ispec_is_acoustic_DG(ispec).
      enddo ! Enddo on ispec.
    !endif ! Endif on ispec_is_acoustic_DG(ispec_selected_source(i_source)).
    
    !write(*,*) "MINMAX SSF", minval(source_spatial_function_DG), maxval(source_spatial_function_DG)
    where(source_spatial_function_DG(:, :)<TINYVAL) source_spatial_function_DG=ZEROcr ! Set to zero where value is too small.
    !write(*,*) "MINMAX SSF", minval(source_spatial_function_DG), maxval(source_spatial_function_DG)
    
    if(SPREAD_SSF_CUSTOM/=0) then
      if(myrank == 0) then
        write(*,*) "********************************"
        write(*,*) "*         INFORMATION          *"
        write(*,*) "********************************"
        write(*,*) "* A custom and spread source   *"
        write(*,*) "* spatial function was used.   *"
        write(*,*) "* The actual formulaes are set *"
        write(*,*) "* up directly in               *"
        write(*,*) "* 'prepare_source_spatial_function.f90'. *"
        write(*,*) "********************************"
        call flush_IMAIN()
      endif
    endif
    
    if(SPREAD_SSF_SAVE) then
      close(504)
      if(myrank == 0) then
        write(*,*) "********************************"
        write(*,*) "*         INFORMATION          *"
        write(*,*) "********************************"
        write(*,*) "* The spread source spatial    *"
        write(*,*) "* function's values for source *"
        write(*,*) "* number ", i_source, " at the  *"
        write(*,*) "* mesh's points were saved in  *"
        write(*,*) "* the OUTPUT_FILES folder. Use *"
        write(*,*) "* the Matlab script            *"
        write(*,*) "* '/utils_new/show_SSF.m' to   *"
        write(*,*) "* plot them.                   *"
        write(*,*) "********************************"
        call flush_IMAIN()
      endif
    endif
  enddo ! Enddo on i_source.
end subroutine prepare_source_spatial_function_DG

subroutine CUSTOM_SSF_PLANE_WAVE(X, SSF_out)
  use constants, only: CUSTOM_REAL, NDIM, PI
  use specfem_par, only: CUSTOM_SSF_P00, CUSTOM_SSF_P01
  real(kind=CUSTOM_REAL), dimension(NDIM), intent(in) :: X
  real(kind=CUSTOM_REAL), intent(out) :: SSF_out
  SSF_out = exp( -((X(NDIM) - CUSTOM_SSF_P00 - tan(CUSTOM_SSF_P01*PI/180.)*X(1)) / 2.)**2 )
end subroutine CUSTOM_SSF_PLANE_WAVE

subroutine CUSTOM_SSF_MOUNTAIN_WAVE(X, Xs, SSF_out)
  use constants, only: CUSTOM_REAL, NDIM
  real(kind=CUSTOM_REAL), dimension(NDIM), intent(in) :: X, Xs
  real(kind=CUSTOM_REAL), intent(out) :: SSF_out
  real(kind=CUSTOM_REAL) :: CPAR1, CPAR2, CPAR3
  ! Moutain gravity waves.
  ! Sabatini, R., Snively, J. B., Bailly, C., Hickey, M. P., & Garrison, J. L. (2019). Numerical Modeling of the Propagation of Infrasonic Acoustic Waves Through the Turbulent Field Generated by the Breaking of Mountain Gravity Waves. Geophysical Research Letters, 46(10), 5526–5534. https://doi.org/10.1029/2019GL082456.
  CPAR1 = 15000. ! SIG_X
  CPAR2 =  3000. ! SIG_Z
  CPAR3 = 2. * 3.14159265358979323846 / CPAR1 ! K1
  SSF_out =   exp( -(   ((X(1)-Xs(1))**2/(2.*CPAR1**2)) &
                      + ((X(2)-Xs(2))**2/(2.*CPAR2**2))   ) ) &
            * cos(CPAR3*(X(1)-Xs(1)))
end subroutine CUSTOM_SSF_MOUNTAIN_WAVE








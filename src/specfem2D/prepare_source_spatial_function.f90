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
  
  !integer :: ig, i, j, ispec ! DEBUG

  if(initialfield) then
    ! One uses an initialfield. Thus, do a dummy allocation.
    allocate(source_spatial_function_DG(1, 1))
  else
    ! One does not use an initial field.
    if(.not. any_acoustic_DG) then
      ! No DG elements exist, ignore call and inform user.
      if(myrank == 0) then
        write(IMAIN,*) 'No DG elements exist, spread source spatial function(s) cannot be implemented. ',&
                       'The classical spatial source function(s) will be used.'
        call flush_IMAIN()
      endif
    else
      ! DG elements exist.
      if(myrank == 0) then
        write(IMAIN,*) 'Spread source spatial function(s) is (are) being initialised over the DG elements.'
        write(IMAIN,*) 'WARNING: THIS KIND OF SOURCE HAS TO BE TUNED NOT TO GENERATE SPURIOUS DISCONTINUITIES, ',&
                       'USE AT YOUR OWN RISK.'
        call flush_IMAIN()
      endif
      ! TODO: Implement a less cumbersome way of storing the values. Indeed, many values far away from the source will be negligible in this array.
      ! For example, instead of storing values over all the points, build the list of indices which associated points have a non-negligible value (value greater than a given threshold), and store only the values of the source spatial function at those points.
      allocate(source_spatial_function_DG(NSOURCES, nglob))
      source_spatial_function_DG(:, :) = ZERO
      call prepare_source_spatial_function_DG() ! Compute the SSF array.
      !if(.false.) then ! DEBUG
      !  !write(*, *) ">  proc", myrank, "initialfield", initialfield
      !  !write(*, *) source_spatial_function_DG ! DEBUG
      !  !write(*, *) 'ouloulou nspec', nspec, 'nglob_DG', nglob_DG, "nglob", nglob ! DEBUG
      !  do ispec = 1, nspec
      !    do i = 1, NGLLX
      !      do j = 1, NGLLZ
      !        ig = ibool_before_perio(i, j, ispec)
      !        if(.false. .and. source_spatial_function_DG(1, ig)>9.5d-1) &
      !          write(*, *) ">  proc", myrank, "ig", ig,&
      !                      "xy", coord(1, ig), coord(2, ig), &
      !                      "ssf", source_spatial_function_DG(1, ig) ! DEBUG
      !      enddo
      !    enddo
      !  enddo
      !endif ! Endif on DEBUG.
    endif ! Endif on any_acoustic_DG.
  endif ! Endif on initialfield.
  
  call synchronize_all()

end subroutine prepare_timerun_ssf

! ------------------------------------------------------------ !
! prepare_source_spatial_function_DG                           !
! ------------------------------------------------------------ !
! Compute values of the source spatial function at all points.

subroutine prepare_source_spatial_function_DG
  
  use constants, only: CUSTOM_REAL, NGLLX, NGLLZ
  use specfem_par, only: coord, ibool_before_perio, IMAIN, ispec_is_acoustic_DG, &
                         ispec_selected_source, myrank, NSOURCES, nspec, &
                         source_spatial_function_DG, source_type, SPREAD_SSF_SAVE, SPREAD_SSF_SIGMA, x_source, z_source 
  implicit none
  
  ! Local variables.
  integer :: i_source, iglob_unique, ispec, i, j
  real(kind=CUSTOM_REAL) :: distsqrd
  
  ! Save values.
  character(len=24) :: filename ! 16 for "OUTPUT_FILES/SSF" + N for process numbering.
  
  if(SPREAD_SSF_SAVE) then
    write(filename, '( "OUTPUT_FILES/SSF", i8.8 )' ) myrank
    open(unit=504,file=filename,status='unknown',action='write', position="append")
  endif
  do i_source = 1, NSOURCES ! Loop on sources.
    if(.not. ispec_is_acoustic_DG(ispec_selected_source(i_source))) then
      ! The central point of the source is not in a DG element, ignore call and inform user.
      if(myrank == 0) then
        write(IMAIN,*) 'The central point of the source ', i_source, ' is not in a DG element, a spread source spatial function ',&
                       'cannot be implemented. A classical spatial source function will be used.'
        call flush_IMAIN()
      endif
    else
      ! The central point of the source is indeed in a DG element.
      ! TODO: More cases may exist where the use of a spread source spatial function has to be discarded. Identify those and produce relevant tests and error messages.
      do ispec = 1, nspec
        if(ispec_is_acoustic_DG(ispec)) then
          do i = 1, NGLLX
            do j = 1, NGLLZ
              iglob_unique = ibool_before_perio(i, j, ispec)
              if(source_type(i_source) == 1) then
                ! If the source is an elastic force or an acoustic pressure.
                distsqrd =   (coord(1, iglob_unique) - x_source(i_source))**2. &
                           + (coord(2, iglob_unique) - z_source(i_source))**2.
                source_spatial_function_DG(i_source, iglob_unique) = exp(-distsqrd/(SPREAD_SSF_SIGMA**2.))
                !source_spatial_function_DG(i_source, iglob_unique) = sin(-distsqrd/(SIGMA_SSF**2.)) ! Test purposes.
                !if(distsqrd<SIGMA_SSF**2*log(10.)*8) source_spatial_function_DG(i_source, iglob_unique) = exp(-distsqrd/(SIGMA_SSF**2.)) ! Only set the SSF where SSF(x) > 10^(-8).
              endif ! Endif source_type. ! TODO: Implement the case source_type = 2.
              
              ! Plane waves tests.
              if(.true.) then
                !distsqrd = (coord(2, iglob_unique) - z_source(i_source))**2. ! Horizontal plane wave.
                !source_spatial_function_DG(i_source, iglob_unique) = exp(-distsqrd/0.5)
                distsqrd = (coord(2, iglob_unique) - z_source(i_source) - tan(3.1415*(0.25))*(coord(1, iglob_unique)-0.))**2. ! Oblique plane wave.
                !distsqrd = (coord(2, iglob_unique) - z_source(i_source) - tan(0.)*(coord(1, iglob_unique)-0.))**2. ! Horizontal plane wave.
                if(abs(coord(2, iglob_unique))<15. .and. abs(coord(1, iglob_unique))<35.) then
                  source_spatial_function_DG(i_source, iglob_unique) = exp(-distsqrd/0.1)
                else
                  source_spatial_function_DG(i_source, iglob_unique) = 0.
                endif
              endif
              
              if(SPREAD_SSF_SAVE) then
                write(504,*) coord(1, iglob_unique), coord(2, iglob_unique), source_spatial_function_DG(1, iglob_unique)
              endif
            enddo ! Enddo on j.
          enddo ! Enddo on i.
        endif ! Endif on ispec_is_acoustic_DG(ispec).
      enddo ! Enddo on ispec.
    endif ! Endif on ispec_is_acoustic_DG(ispec_selected_source(i_source)).
  enddo ! Enddo on i_source.
  if(SPREAD_SSF_SAVE) then
    close(504)
    if(myrank == 0) then
      write(IMAIN,*) "The spread source spatial function's values at the mesh's points were saved in the OUTPUT_FILES folder. ",&
                     "Use the Matlab script '/utils_new/show_SSF.m' to plot."
      call flush_IMAIN()
    endif
  endif
end subroutine prepare_source_spatial_function_DG

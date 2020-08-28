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
!=====================================================================

! TODO: Unused, decide what to do with this.

  subroutine compute_add_sources_NS(dt_rhoveloc_acoustic,it,i_stage)

  use constants,only: CUSTOM_REAL,NGLLX,NGLLZ,NDIM

  use specfem_par, only: ispec_is_acoustic,nglob_acoustic,&
                         NSOURCES,source_type,source_time_function,&
                         is_proc_source,ispec_selected_source,&
                         hxis_store,hgammas_store,ibool,myrank
  implicit none

  real(kind=CUSTOM_REAL), dimension(NDIM, nglob_acoustic),intent(inout) :: dt_rhoveloc_acoustic
  integer,intent(in) :: it,i_stage

  !local variables
  integer :: i_source,i,j,iglob
  double precision :: hlagrange

  do i_source= 1,NSOURCES
    ! if this processor core carries the source and the source element is acoustic
    if (is_proc_source(i_source) == 1 .and. ispec_is_acoustic(ispec_selected_source(i_source))) then
      ! collocated force
      ! beware, for acoustic medium, source is: pressure divided by Kappa of the fluid
      ! the sign is negative because pressure p = - Chi_dot_dot therefore we need
      ! to add minus the source to Chi_dot_dot to get plus the source in pressure
      if (source_type(i_source) == 1) then
        ! forward wavefield
        do j = 1,NGLLZ
          do i = 1,NGLLX
            iglob = ibool(i,j,ispec_selected_source(i_source))
            hlagrange = hxis_store(i_source,i) * hgammas_store(i_source,j)
            dt_rhoveloc_acoustic(1, iglob) = dt_rhoveloc_acoustic(1, iglob) + &
                         source_time_function(i_source,it,i_stage)*hlagrange 
            dt_rhoveloc_acoustic(2, iglob) = dt_rhoveloc_acoustic(2, iglob) + &
                         source_time_function(i_source,it,i_stage)*hlagrange 
                                              !  WRITE(*,*) "TOTO", dt_rhoveloc_acoustic(2, iglob)
          enddo
        enddo
      ! moment tensor
      else if (source_type(i_source) == 2) then
         call exit_MPI(myrank,'cannot have moment tensor source in acoustic element')
      endif
    endif ! if this processor core carries the source and the source element is acoustic
  enddo ! do i_source= 1,NSOURCES

  end subroutine compute_add_sources_NS


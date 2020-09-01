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

! for acoustic solver

  subroutine compute_add_sources_acoustic_DG(variable_DG,it,i_stage)

  use constants,only: CUSTOM_REAL,NGLLX,NGLLZ

  use specfem_par, only: ispec_is_acoustic,nglob_DG,&
                         NSOURCES,source_type,source_time_function,&
                         is_proc_source,ispec_selected_source,&
                         hxis_store,hgammas_store,ibool_DG,myrank, &
                         neighbor_DG, neighbor_DG_corner, is_corner, &
                         ispec_is_acoustic, link_DG_CG, coord, myrank
  implicit none

  real(kind=CUSTOM_REAL), dimension(nglob_DG),intent(inout) :: variable_DG
  integer,intent(in) :: it,i_stage

  !local variables
  integer :: i_source,i,j,ispec,iglob,i_neighbor,j_neighbor,ispec_neighbor,k, cpt_loc,i2
  integer, dimension(4) :: iglob_surface
  double precision :: hlagrange

  do i_source= 1,NSOURCES
    ! if this processor core carries the source and the source element is acoustic
    if (is_proc_source(i_source) == 1 .and. ispec_is_acoustic(ispec_selected_source(i_source))) then
      ! collocated force
      ! beware, for acoustic medium, source is: pressure divided by Kappa of the fluid
      ! the sign is negative because pressure p = - Chi_dot_dot therefore we need
      ! to add minus the source to Chi_dot_dot to get plus the source in pressure
      if (ispec_is_acoustic(ispec_selected_source(i_source))) then
      if (source_type(i_source) == 1) then
        ! forward wavefield
        do j = 1,NGLLZ
          do i = 1,NGLLX
            hlagrange = hxis_store(i_source,i) * hgammas_store(i_source,j)
            
            ispec = ispec_selected_source(i_source)
            
            if(.false.) then
            ! Find neighbors if on the element boundary
            !cpt = 1
            cpt_loc = 1
            !if(i == 1 .OR. j == 1 .OR. i == NGLLX .OR. j == NGLLZ) then
            !    cpt = 2
            !    if(is_corner(i,j)) cpt = 4
            !endif
            
            iglob_surface(cpt_loc) = ibool_DG(i,j,ispec)
            
            i_neighbor = neighbor_DG(i,j,ispec,1)
            j_neighbor = neighbor_DG(i,j,ispec,2)
            ispec_neighbor = neighbor_DG(i,j,ispec,3)
            if(i_neighbor > -1) then
            cpt_loc = cpt_loc + 1
            iglob_surface(cpt_loc) = ibool_DG(i_neighbor,j_neighbor,ispec_neighbor)
            endif
            
            if(is_corner(i,j)) then
            
            i_neighbor = neighbor_DG_corner(i,j,ispec,1)
            j_neighbor = neighbor_DG_corner(i,j,ispec,2)
            ispec_neighbor = neighbor_DG_corner(i,j,ispec,3)
            if(i_neighbor > -1) then
            
                cpt_loc = cpt_loc + 1
                iglob_surface(cpt_loc) = ibool_DG(i_neighbor,j_neighbor,ispec_neighbor)
                
                if(neighbor_DG(i_neighbor,j_neighbor,ispec_neighbor,3) /= ispec) then
                    i2 = neighbor_DG(i_neighbor,j_neighbor,ispec_neighbor,1)
                    j_neighbor = neighbor_DG(i_neighbor,j_neighbor,ispec_neighbor,2)
                    ispec_neighbor = neighbor_DG(i_neighbor,j_neighbor,ispec_neighbor,3)
                else
                    !i_neighbor = neighbor_DG_corner(i_neighbor,j_neighbor,ispec_neighbor,1)
                    !j_neighbor = neighbor_DG_corner(i_neighbor,j_neighbor,ispec_neighbor,2)
                    !ispec_neighbor = neighbor_DG_corner(i_neighbor,j_neighbor,ispec_neighbor,3)
                    !WRITE(*,*) "2 ---->", neighbor_DG_corner(i_neighbor,j_neighbor,ispec_neighbor,1)
                    i2 = neighbor_DG_corner(i_neighbor,j_neighbor,ispec_neighbor,1)
                    j_neighbor = neighbor_DG_corner(i_neighbor,j_neighbor,ispec_neighbor,2)
                    ispec_neighbor = neighbor_DG_corner(i_neighbor,j_neighbor,ispec_neighbor,3)
                endif
                    
                if(i2 > -1) then
                        cpt_loc = cpt_loc + 1
                        iglob_surface(cpt_loc) = ibool_DG(i2,j_neighbor,ispec_neighbor)
                endif
            
            endif !if(i_neighbor > -1)
            
            endif !is_corner(i,j)
            
            do k = 1,cpt_loc
            iglob = iglob_surface(k)
            WRITE(*,*) myrank,i,j,ispec,"SOURCE POINTS : ", iglob, coord(:, link_DG_CG(iglob))
            variable_DG(iglob) = variable_DG(iglob) - source_time_function(i_source,it,i_stage)*hlagrange
            enddo
            WRITE(*,*) "----------------"
            endif
            
            variable_DG(ibool_DG(i,j,ispec)) = variable_DG(ibool_DG(i,j,ispec)) &
                - source_time_function(i_source,it,i_stage)*hlagrange
            
          enddo
        enddo
      ! moment tensor
      else if (source_type(i_source) == 2) then
         call exit_MPI(myrank,'cannot have moment tensor source in acoustic element')
      endif
      endif !(ispec_is_acoustic(ispec_selected_source(i_source)))
    endif ! if this processor core carries the source and the source element is acoustic
  enddo ! do i_source= 1,NSOURCES

  end subroutine compute_add_sources_acoustic_DG


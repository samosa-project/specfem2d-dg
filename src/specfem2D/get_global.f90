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

! create a sorted version of the indirect addressing to reduce cache misses

  subroutine get_global()

  use constants,only: NGLLX,NGLLZ

  use specfem_par, only : nspec,ibool,copy_ibool_ori,integer_mask_ibool,SAVE_MODEL,myrank, &
        ! MODIF DG
        ibool_before_perio, integer_mask_ibool_bef, copy_ibool_ori_bef, &
        USE_DISCONTINUOUS_METHOD

  implicit none

  ! local parameters
  integer :: inumber,ispec,i,j, &
        inumber_bef
  integer :: ier
  character(len=150) :: outputname

  ! initializes temporary arrays
  integer_mask_ibool(:) = -1
  copy_ibool_ori(:,:,:) = ibool(:,:,:)
  
  integer_mask_ibool_bef(:) = -1
  ! Save a copy of ibool before renumbering
  ! TODO: Needs to be corrected because useless
  if(USE_DISCONTINUOUS_METHOD) ibool_before_perio(:,:,:) = ibool(:,:,:)
  if(USE_DISCONTINUOUS_METHOD) copy_ibool_ori_bef(:,:,:) = ibool_before_perio(:,:,:)
  
  inumber = 0
  inumber_bef = 0
  
  ! reduce cache misses in all the elements
  ! loop over spectral elements
  do ispec = 1,nspec
    do j = 1,NGLLZ
      do i = 1,NGLLX
      
        if (integer_mask_ibool(copy_ibool_ori(i,j,ispec)) == -1) then
          ! create a new point
          inumber = inumber + 1
          ibool(i,j,ispec) = inumber
          integer_mask_ibool(copy_ibool_ori(i,j,ispec)) = inumber
        else
          ! use an existing point created previously
          ibool(i,j,ispec) = integer_mask_ibool(copy_ibool_ori(i,j,ispec))
        endif
        
        if(USE_DISCONTINUOUS_METHOD) then
          if (integer_mask_ibool_bef(copy_ibool_ori_bef(i,j,ispec)) == -1) then
            ! create a new point
            inumber_bef = inumber_bef + 1
            ibool_before_perio(i,j,ispec) = inumber_bef
            integer_mask_ibool_bef(copy_ibool_ori_bef(i,j,ispec)) = inumber_bef
          else
            ! use an existing point created previously
            ibool_before_perio(i,j,ispec) = integer_mask_ibool_bef(copy_ibool_ori_bef(i,j,ispec))
          endif
        endif
      enddo
    enddo
  enddo
  
  if (trim(SAVE_MODEL) == 'binary') then
    write(outputname,'(a,i6.6,a)') './DATA/proc',myrank,'_NSPEC_ibool.bin'

    open(888,file=trim(outputname),status='unknown',form='unformatted',iostat=ier)
    if (ier /= 0) stop 'Error opening smoothed kernel file'
    write(888) nspec
    write(888) ibool
    close(888)
  endif

  end subroutine get_global

!
!-------------------------------------------------------------------------------------------------
!

! create a sorted version of the indirect addressing to reduce cache misses

  subroutine get_global_indirect_addressing(nspec,nglob,ibool,copy_ibool_ori,integer_mask_ibool)

  use constants,only: NGLLX,NGLLZ

  implicit none

  integer :: nspec,nglob

  integer, dimension(NGLLX,NGLLZ,nspec) :: ibool,copy_ibool_ori
  integer, dimension(nglob) :: integer_mask_ibool

  ! local parameters
  integer :: inumber,ispec,i,j

  ! initializes temporary arrays
  integer_mask_ibool(:) = -1
  copy_ibool_ori(:,:,:) = ibool(:,:,:)

  inumber = 0

  do ispec = 1,nspec
    do j = 1,NGLLZ
      do i = 1,NGLLX
        if (integer_mask_ibool(copy_ibool_ori(i,j,ispec)) == -1) then
          ! create a new point
          inumber = inumber + 1
          ibool(i,j,ispec) = inumber
          integer_mask_ibool(copy_ibool_ori(i,j,ispec)) = inumber
        else
          ! use an existing point created previously
          ibool(i,j,ispec) = integer_mask_ibool(copy_ibool_ori(i,j,ispec))
        endif
      enddo
    enddo
  enddo

  end subroutine get_global_indirect_addressing


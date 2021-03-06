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

#ifdef USE_MPI

  subroutine get_MPI()

! sets up the MPI interface for communication between partitions

  use mpi

  use constants,only: CUSTOM_REAL,IMAIN,NGLLX,NGLLZ

  use specfem_par, only: nspec,ibool,nglob, &
                    ninterface,ibool_interfaces_acoustic,ibool_interfaces_elastic, &
                    ibool_interfaces_poroelastic, &
                    nibool_interfaces_acoustic,nibool_interfaces_elastic, &
                    nibool_interfaces_poroelastic, &
                    ninterface_acoustic,ninterface_elastic,ninterface_poroelastic, &
                    myrank,coord,buffer_send_faces_vector_ac,buffer_recv_faces_vector_ac,&
                    tab_requests_send_recv_acoustic

  implicit none

  !local parameters
  double precision, dimension(:), allocatable :: xp,zp
  double precision, dimension(:), allocatable :: work
  integer, dimension(:), allocatable :: locval
  integer, dimension(:), allocatable :: nibool_interfaces_true
  ! for MPI buffers
  integer, dimension(:), allocatable :: reorder_interface,ind,ninseg,iwork
  integer, dimension(:), allocatable :: ibool_dummy
  
!  integer, dimension(:,:), allocatable :: ibool_interfaces_dummy
  logical, dimension(:), allocatable :: ifseg
  integer :: iinterface,ilocnum
  integer :: num_points1, num_points2
  ! assembly test
  integer :: i,j,ispec,iglob,countval,inum,idomain
  integer :: max_nibool_interfaces,num_nibool,num_interface
  
  real(kind=CUSTOM_REAL), dimension(:),allocatable :: test_flag_cr

  ! gets global indices for points on MPI interfaces
  ! (defined by my_interfaces) between different partitions
  ! and stores them in ibool_interfaces*** & nibool_interfaces*** (number of total points)
  call prepare_assemble_MPI()

  ! sorts ibool comm buffers lexicographically for all MPI interfaces
  num_points1 = 0
  num_points2 = 0
  allocate(nibool_interfaces_true(ninterface))

  do idomain = 1,3

    ! checks number of interface in this domain
    num_interface = 0
    if (idomain == 1) then
      num_interface = ninterface_acoustic
    else if (idomain == 2) then
      num_interface = ninterface_elastic
    else if (idomain == 3) then
      num_interface = ninterface_poroelastic
    endif
    if (num_interface == 0) cycle

    ! loops over interfaces
    do iinterface = 1, ninterface

      ! number of global points in this interface
      num_nibool = 0
      if (idomain == 1) then
        num_nibool = nibool_interfaces_acoustic(iinterface)
      else if (idomain == 2) then
        num_nibool = nibool_interfaces_elastic(iinterface)
      else if (idomain == 3) then
        num_nibool = nibool_interfaces_poroelastic(iinterface)
      endif
      ! checks if anything to sort
      if (num_nibool == 0) cycle

      allocate(xp(num_nibool))
      allocate(zp(num_nibool))
      allocate(locval(num_nibool))
      allocate(ifseg(num_nibool))
      allocate(reorder_interface(num_nibool))
      allocate(ibool_dummy(num_nibool))
      allocate(ind(num_nibool))
      allocate(ninseg(num_nibool))
      allocate(iwork(num_nibool))
      allocate(work(num_nibool))

      ! works with a copy of ibool array
      if (idomain == 1) then
        ibool_dummy(:) = ibool_interfaces_acoustic(1:num_nibool,iinterface)
      else if (idomain == 2) then
        ibool_dummy(:) = ibool_interfaces_elastic(1:num_nibool,iinterface)
      else if (idomain == 3) then
        ibool_dummy(:) = ibool_interfaces_poroelastic(1:num_nibool,iinterface)
      endif

      ! gets x,y,z coordinates of global points on MPI interface
      do ilocnum = 1, num_nibool
        iglob = ibool_dummy(ilocnum)
        xp(ilocnum) = coord(1,iglob)
        zp(ilocnum) = coord(2,iglob)
      enddo

      ! sorts (lexicographically?) ibool_interfaces and updates value
      ! of total number of points nibool_interfaces_true(iinterface)
      call sort_array_coordinates(num_nibool, &
                                xp,zp, &
                                ibool_dummy, &
                                reorder_interface,locval,ifseg, &
                                nibool_interfaces_true(iinterface), &
                                ind,ninseg,iwork,work)

      ! checks that number of MPI points are still the same
      num_points1 = num_points1 + num_nibool
      num_points2 = num_points2 + nibool_interfaces_true(iinterface)
      
      if (num_points1 /= num_points2) then
        write(IMAIN,*) 'Error sorting MPI interface points:',myrank
        write(IMAIN,*) '   domain:',idomain
        write(IMAIN,*) '   interface:',iinterface,num_points1,num_points2
        call exit_MPI(myrank,'error sorting MPI interface')
      endif

      ! stores new order of ibool array
      if (idomain == 1) then
        ibool_interfaces_acoustic(1:num_nibool,iinterface) = ibool_dummy(:)
        !do i = 1, num_nibool 
        !        !WRITE(*,*) myrank, "TEST ACOUSTIC AFTER SORTING", ibool_interfaces_acoustic(i,iinterface), &
        !        !        coord(:,ibool_interfaces_acoustic(i,iinterface))
        !        if(myrank == 0) &
        !        WRITE(*,*) coord(:,ibool_interfaces_acoustic(i,iinterface))
        !enddo
        !stop
      else if (idomain == 2) then
        ibool_interfaces_elastic(1:num_nibool,iinterface) = ibool_dummy(:)
      else if (idomain == 3) then
        ibool_interfaces_poroelastic(1:num_nibool,iinterface) = ibool_dummy(:)
      endif
      
      !stop 'KIKIKIKIKI 10'

      ! cleanup temporary arrays
      deallocate(xp)
      deallocate(zp)
      deallocate(locval)
      deallocate(ifseg)
      deallocate(reorder_interface)
      deallocate(ibool_dummy)
      deallocate(ind)
      deallocate(ninseg)
      deallocate(iwork)
      deallocate(work)
    enddo
  enddo

  ! cleanup
  deallocate(nibool_interfaces_true)

  ! outputs total number of MPI interface points
  call sum_all_i(num_points2,num_points1)
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'total MPI interface points: ',num_points1
    call flush_IMAIN()
  endif

  ! checks interfaces in acoustic domains
  inum = 0
  countval = 0

  if (ninterface_acoustic > 0) then

    ! checks with assembly of test fields
    allocate(test_flag_cr(nglob))
    test_flag_cr(:) = 0._CUSTOM_REAL
    countval = 0
    do ispec = 1, nspec
      ! sets flags on global points
      do j = 1, NGLLZ
        do i = 1, NGLLX
          ! global index
          iglob = ibool(i,j,ispec)

          ! counts number of unique global points to set
          if (nint(test_flag_cr(iglob)) == 0) countval = countval + 1

          ! sets identifier
          test_flag_cr(iglob) = myrank + 1.0
        enddo
      enddo
    enddo

    max_nibool_interfaces = maxval(nibool_interfaces_acoustic(:))

    allocate(tab_requests_send_recv_acoustic(ninterface_acoustic*2))
    allocate(buffer_send_faces_vector_ac(max_nibool_interfaces,ninterface_acoustic))
    allocate(buffer_recv_faces_vector_ac(max_nibool_interfaces,ninterface_acoustic))
    inum = 0
    do iinterface = 1, ninterface
      inum = inum + nibool_interfaces_acoustic(iinterface)
    enddo
  endif
  
  ! note: this mpi reduction awaits information from all processes.
  !          thus, avoid an mpi deadlock in case some of the paritions have no acoustic interface
  call sum_all_i(inum,num_points1)
  if (myrank == 0) then
    write(IMAIN,*) '    DG acoustic interface points: ',num_points1
    call flush_IMAIN()
  endif
  
  ! checks if assembly works
  inum = 0
  if (ninterface_acoustic > 0) then
    ! adds contributions from different partitions to flag arrays
    ! custom_real arrays
    call assemble_MPI_vector_ac(test_flag_cr)

    ! checks number of interface points
    inum = 0
    do iglob= 1,nglob
      ! only counts flags with MPI contributions
      if (nint(test_flag_cr(iglob)) > myrank+1 ) inum = inum + 1
    enddo

    deallocate(tab_requests_send_recv_acoustic)
    deallocate(buffer_send_faces_vector_ac)
    deallocate(buffer_recv_faces_vector_ac)
    deallocate(test_flag_cr)
  endif
  
  ! note: this mpi reduction awaits information from all processes.
  call sum_all_i(inum,num_points2)
  if (myrank == 0) then
    write(IMAIN,*) '       assembly acoustic MPI interface points:',num_points2
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  end subroutine get_MPI

#endif

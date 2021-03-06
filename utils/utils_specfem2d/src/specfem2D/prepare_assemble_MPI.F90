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

!
! This file contains subroutines related to assembling (of the mass matrix, potential_dot_dot and
! accel_elastic, accels_poroelastic, accelw_poroelastic).
! These subroutines are for the most part not used in the sequential version.
!

#ifdef USE_MPI

!-----------------------------------------------
! Determines the points that are on the interfaces with other partitions, to help
! build the communication buffers, and determines which elements are considered 'inner'
! (no points in common with other partitions) and 'outer' (at least one point in common
! with neighbouring partitions).
! We have both acoustic and (poro)elastic buffers, for coupling between acoustic and (poro)elastic elements
! led us to have two sets of communications.
!-----------------------------------------------
  subroutine prepare_assemble_MPI()

  use constants,only: NGLLX,NGLLZ

  use specfem_par, only: nspec,ibool,ibool_DG,nglob_DG,knods, ngnod,nglob, ispec_is_elastic, ispec_is_poroelastic, &
                                ninterface, &
                                my_nelmnts_neighbours, my_interfaces, &
                                ibool_interfaces_acoustic, ibool_interfaces_elastic, &
                                ibool_interfaces_poroelastic, &
                                nibool_interfaces_acoustic, nibool_interfaces_elastic, &
                                nibool_interfaces_poroelastic, &
                                inum_interfaces_acoustic, inum_interfaces_elastic, &
                                inum_interfaces_poroelastic, &
                                ninterface_acoustic, ninterface_elastic, ninterface_poroelastic, &
                                ! MODIF DG
                                ibool_interfaces_acoustic_DG, nibool_interfaces_acoustic_DG, &
                                inum_interfaces_acoustic_DG, ninterface_acoustic_DG, &
                                mask_ispec_inner_outer,nibool_interfaces_ext_mesh, ibool_interfaces_ext_mesh_init, &
                                is_MPI_interface_DG!, coord, myrank, link_DG_CG!, neighbor_DG, neighbor_DG_corner, &
                                !ispec_is_acoustic_coupling_el, ispec_is_acoustic_surface_corner, &
                                !ispec_is_acoustic_surface
  implicit none

  ! local parameters
  integer  :: num_interface
  integer  :: ispec_interface
  logical, dimension(nglob)  :: mask_ibool_acoustic
  logical, dimension(nglob)  :: mask_ibool_elastic
  logical, dimension(nglob)  :: mask_ibool_poroelastic
  logical, dimension(nglob)  :: mask_ibool_ext_mesh
  integer  :: ixmin, ixmax, izmin, izmax, ix, iz
  integer, dimension(ngnod)  :: n
  integer  :: e1, e2, itype, ispec, k, sens, iglob, iglob_DG
  integer  :: nglob_interface_acoustic
  ! MODIF DG
  integer  :: nglob_interface_acoustic_DG
  integer  :: nglob_interface_elastic
  integer  :: nglob_interface_poroelastic
  integer :: npoin_interface_ext_mesh

  ! initializes
  ibool_interfaces_acoustic(:,:) = 0
  nibool_interfaces_acoustic(:) = 0
  ! MODIF DG
  ibool_interfaces_acoustic_DG(:,:) = 0
  nibool_interfaces_acoustic_DG(:) = 0
  ibool_interfaces_elastic(:,:) = 0
  nibool_interfaces_elastic(:) = 0
  ibool_interfaces_poroelastic(:,:) = 0
  nibool_interfaces_poroelastic(:) = 0
  nibool_interfaces_ext_mesh(:) = 0
  ibool_interfaces_ext_mesh_init(:,:) = 0
  
  allocate(is_MPI_interface_DG(nglob_DG))
  is_MPI_interface_DG = .false.

  do num_interface = 1, ninterface
    ! initializes interface point counters
    nglob_interface_acoustic = 0
    ! MODIF DG
    nglob_interface_acoustic_DG = 0
    nglob_interface_elastic = 0
    nglob_interface_poroelastic = 0
    npoin_interface_ext_mesh = 0
    mask_ibool_acoustic(:) = .false.
    mask_ibool_elastic(:) = .false.
    mask_ibool_poroelastic(:) = .false.
    mask_ibool_ext_mesh(:) = .false.

    do ispec_interface = 1, my_nelmnts_neighbours(num_interface)
      ! element id
      ispec = my_interfaces(1,ispec_interface,num_interface)
      ! type of interface: 1 = common point, 2 = common edge
      itype = my_interfaces(2,ispec_interface,num_interface)
      ! element control node ids
      do k = 1, ngnod
        n(k) = knods(k,ispec)
      enddo
      ! common node ids
      e1 = my_interfaces(3,ispec_interface,num_interface)
      e2 = my_interfaces(4,ispec_interface,num_interface)

      call get_edge(ngnod, n, itype, e1, e2, ixmin, ixmax, izmin, izmax, sens)

      do iz = izmin, izmax, sens
        do ix = ixmin, ixmax, sens
          ! global index
          iglob = ibool(ix,iz,ispec)
          
          !WRITE(100,*) myrank, ix, iz,ispec, coord(:,iglob)
          !WRITE(*,*) myrank,"COORD",coord(:,iglob)
          
          ! MODIF DG
          iglob_DG = ibool_DG(ix,iz,ispec)

          if (.not. mask_ibool_ext_mesh(iglob)) then
            ! masks point as being accounted for
            mask_ibool_ext_mesh(iglob) = .true.
            ! adds point to interface
            npoin_interface_ext_mesh = npoin_interface_ext_mesh + 1
            ibool_interfaces_ext_mesh_init(npoin_interface_ext_mesh,num_interface) = iglob
          endif

          ! checks to which material this common interface belongs
          if (ispec_is_elastic(ispec)) then
            ! elastic element
            if (.not. mask_ibool_elastic(iglob)) then
              mask_ibool_elastic(iglob) = .true.
              nglob_interface_elastic = nglob_interface_elastic + 1
              ibool_interfaces_elastic(nglob_interface_elastic,num_interface) = iglob
            endif
          else if (ispec_is_poroelastic(ispec)) then
            ! poroelastic element
            if (.not. mask_ibool_poroelastic(iglob)) then
              mask_ibool_poroelastic(iglob) = .true.
              nglob_interface_poroelastic = nglob_interface_poroelastic + 1
              ibool_interfaces_poroelastic(nglob_interface_poroelastic,num_interface) = iglob
            endif
          else
            ! acoustic element
            nglob_interface_acoustic_DG = nglob_interface_acoustic_DG + 1
            ! MODIF DG
            ibool_interfaces_acoustic_DG(nglob_interface_acoustic_DG,num_interface) = iglob_DG
            is_MPI_interface_DG(iglob_DG) = .true.
            
            ! If at convex corner "x" one has two MPI neighbors (N1 and N2)
            ! => Thus one needs to be able to recover this info 
            !    by increasing the number of acoustic points on the interface
            ! (!) But this is valid only when it is not an outside or coupled edge
            !
            !          N1
            !   ---------
            !          x| N2
            !           |
            !           |
            !
            !if(neighbor_DG(ix,iz,ispec,1) == -1 .AND. &
            !   neighbor_DG_corner(ix,iz,ispec,1) == -1 .AND. &
            !   ispec_is_acoustic_coupling_el(ix,iz,ispec,3) < 0 .AND. &
            !   .not. ispec_is_acoustic_surface_corner(ix,iz,ispec) .AND. &
            !   .not. ispec_is_acoustic_surface(ix,iz,ispec)) then
            !    nglob_interface_acoustic_DG = nglob_interface_acoustic_DG + 1
            !    ibool_interfaces_acoustic_DG(nglob_interface_acoustic_DG,num_interface) = iglob_DG
            !    if(myrank == 0) &
            !    WRITE(*,*) myrank,"WHATISCORNER", iglob_DG, coord(:, iglob), ix,iz,ispec
            !endif               
           
            ! acoustic element
            if (.not. mask_ibool_acoustic(iglob)) then
              mask_ibool_acoustic(iglob) = .true.
              nglob_interface_acoustic = nglob_interface_acoustic + 1
              ibool_interfaces_acoustic(nglob_interface_acoustic,num_interface) = iglob
            endif
          endif
        enddo
      enddo

    enddo
    
    ! stores counters for interface points
    nibool_interfaces_acoustic(num_interface) = nglob_interface_acoustic
    ! MODIF DG
    nibool_interfaces_acoustic_DG(num_interface) = nglob_interface_acoustic_DG
    nibool_interfaces_elastic(num_interface) = nglob_interface_elastic
    nibool_interfaces_poroelastic(num_interface) = nglob_interface_poroelastic
    nibool_interfaces_ext_mesh(num_interface) = npoin_interface_ext_mesh
    ! sets inner/outer element flags
    do ispec = 1, nspec
      do iz = 1, NGLLZ
        do ix = 1, NGLLX

           if (mask_ibool_acoustic(ibool(ix,iz,ispec)) &
            .or. mask_ibool_elastic(ibool(ix,iz,ispec)) &
            .or. mask_ibool_poroelastic(ibool(ix,iz,ispec))) then
               mask_ispec_inner_outer(ispec) = .true.
          endif

        enddo
      enddo
    enddo

  enddo
  
  ! sets number of interfaces for each material domain
  ninterface_acoustic = 0
  ! MODIF DG
  ninterface_acoustic_DG = 0
  ninterface_elastic =  0
  ninterface_poroelastic =  0

  ! loops over all MPI interfaces
  do num_interface = 1, ninterface
    ! sets acoustic MPI interface (local) indices in range [1,ninterface_acoustic]
    if (nibool_interfaces_acoustic(num_interface) > 0) then
      ninterface_acoustic = ninterface_acoustic + 1
      inum_interfaces_acoustic(ninterface_acoustic) = num_interface
    endif
    ! MODIF DG
    if (nibool_interfaces_acoustic_DG(num_interface) > 0) then
      ninterface_acoustic_DG = ninterface_acoustic_DG + 1
      inum_interfaces_acoustic_DG(ninterface_acoustic_DG) = num_interface
    endif
    ! elastic
    if (nibool_interfaces_elastic(num_interface) > 0) then
      ninterface_elastic = ninterface_elastic + 1
      inum_interfaces_elastic(ninterface_elastic) = num_interface
    endif
    ! poroelastic
    if (nibool_interfaces_poroelastic(num_interface) > 0) then
      ninterface_poroelastic = ninterface_poroelastic + 1
      inum_interfaces_poroelastic(ninterface_poroelastic) = num_interface
    endif
  enddo

  end subroutine prepare_assemble_MPI


!-----------------------------------------------
! Get the points (ixmin, ixmax, izmin and izmax) on an node/edge for one element.
! 'sens' is used to have DO loops with increment equal to 'sens' (-/+1).
!-----------------------------------------------
  subroutine get_edge ( ngnod, n, itype, e1, e2, ixmin, ixmax, izmin, izmax, sens )

  use constants,only: NGLLX,NGLLZ

  implicit none

  integer, intent(in)  :: ngnod
  integer, dimension(ngnod), intent(in)  :: n
  integer, intent(in)  :: itype, e1, e2
  integer, intent(out)  :: ixmin, ixmax, izmin, izmax
  integer, intent(out)  :: sens

  if (itype == 1) then

    ! common single point

    ! checks which corner point is given
    if (e1 == n(1)) then
        ixmin = 1
        ixmax = 1
        izmin = 1
        izmax = 1
    endif
    if (e1 == n(2)) then
        ixmin = NGLLX
        ixmax = NGLLX
        izmin = 1
        izmax = 1
    endif
    if (e1 == n(3)) then
        ixmin = NGLLX
        ixmax = NGLLX
        izmin = NGLLZ
        izmax = NGLLZ
    endif
    if (e1 == n(4)) then
        ixmin = 1
        ixmax = 1
        izmin = NGLLZ
        izmax = NGLLZ
    endif
    sens = 1

  else if (itype == 2) then

    ! common edge

    ! checks which edge and corner points are given
    if (e1 ==  n(1)) then
        ixmin = 1
        izmin = 1
        if (e2 == n(2)) then
           ixmax = NGLLX
           izmax = 1
           sens = 1
        endif
        if (e2 == n(4)) then
           ixmax = 1
           izmax = NGLLZ
           sens = 1
        endif
     endif
     if (e1 == n(2)) then
        ixmin = NGLLX
        izmin = 1
        if (e2 == n(3)) then
           ixmax = NGLLX
           izmax = NGLLZ
           sens = 1
        endif
        if (e2 == n(1)) then
           ixmax = 1
           izmax = 1
           sens = -1
        endif
     endif
     if (e1 == n(3)) then
        ixmin = NGLLX
        izmin = NGLLZ
        if (e2 == n(4)) then
           ixmax = 1
           izmax = NGLLZ
           sens = -1
        endif
        if (e2 == n(2)) then
           ixmax = NGLLX
           izmax = 1
           sens = -1
        endif
     endif
     if (e1 == n(4)) then
        ixmin = 1
        izmin = NGLLZ
        if (e2 == n(1)) then
           ixmax = 1
           izmax = 1
           sens = -1
        endif
        if (e2 == n(3)) then
           ixmax = NGLLX
           izmax = NGLLZ
           sens = 1
        endif
     endif

  else

    stop 'ERROR get_edge unknown type'

  endif

  end subroutine get_edge

#endif

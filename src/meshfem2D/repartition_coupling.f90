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


  !--------------------------------------------------
  ! repartitioning: coupled acoustic/elastic elements are transferred to the same partition
  !--------------------------------------------------

  subroutine acoustic_elastic_repartitioning(elmnts_l, nbmodels, phi_material, num_material, nproc)

  use part_unstruct_par,only: nelmnts,edges_coupled,nedges_coupled,part

  implicit none
  include "constants.h"

  integer, dimension(0:NCORNERS*nelmnts-1), intent(in)  :: elmnts_l
  integer, intent(in)  :: nproc, nbmodels
  double precision, dimension(nbmodels), intent(in)  :: phi_material
  integer, dimension(1:nelmnts), intent(in)  :: num_material

  ! local parameters
  integer, dimension(:), allocatable  :: xadj_l
  integer, dimension(:), allocatable  :: adjncy_l
  logical, dimension(nbmodels)  :: is_acoustic, is_elastic
  integer  :: i, ier, num_edge
  integer  :: el, el_adj
  logical  :: is_repartitioned

  allocate(xadj_l(0:nelmnts))
  allocate(adjncy_l(0:MAX_NEIGHBORS*nelmnts-1))

  is_acoustic(:) = .false.
  is_elastic(:) = .false.

  do i = 1, nbmodels
     if (phi_material(i) >= 1.d0) then
        is_acoustic(i) = .true.
     endif
     if (phi_material(i) < TINYVAL) then
        is_elastic(i) = .true.
     endif
  enddo

  ! determines maximum neighbors based on 2 common nodes (common edge)
  call mesh2dual_ncommonnodes(elmnts_l, 2, xadj_l, adjncy_l)

  nedges_coupled = 0
  do el = 0, nelmnts-1
     if (is_acoustic(num_material(el+1))) then
        do el_adj = xadj_l(el), xadj_l(el+1) - 1
           if (is_elastic(num_material(adjncy_l(el_adj)+1))) then
              nedges_coupled = nedges_coupled + 1
           endif
        enddo
     endif
  enddo

  print *, 'nedges_coupled (acoustic/elastic)', nedges_coupled

  allocate(edges_coupled(2,nedges_coupled),stat=ier)
  if (ier /= 0) stop 'Error allocating array edges_coupled'

  nedges_coupled = 0
  do el = 0, nelmnts-1
     if (is_acoustic(num_material(el+1))) then
        do el_adj = xadj_l(el), xadj_l(el+1) - 1
           if (is_elastic(num_material(adjncy_l(el_adj)+1))) then
              nedges_coupled = nedges_coupled + 1
              edges_coupled(1,nedges_coupled) = el
              edges_coupled(2,nedges_coupled) = adjncy_l(el_adj)
           endif

        enddo
     endif
  enddo

  do i = 1, nedges_coupled*nproc
     is_repartitioned = .false.
     do num_edge = 1, nedges_coupled
        if (part(edges_coupled(1,num_edge)) /= part(edges_coupled(2,num_edge))) then
           if (part(edges_coupled(1,num_edge)) < part(edges_coupled(2,num_edge))) then
              part(edges_coupled(2,num_edge)) = part(edges_coupled(1,num_edge))
           else
              part(edges_coupled(1,num_edge)) = part(edges_coupled(2,num_edge))
           endif

           is_repartitioned = .true.
        endif

     enddo
     if (.not. is_repartitioned) then
        exit
     endif
  enddo

  deallocate(xadj_l,adjncy_l)

  end subroutine acoustic_elastic_repartitioning

!
!---------------------------------------------------------------------------------------
!

  !--------------------------------------------------
  ! repartitioning: coupled acoustic/poroelastic elements are transferred to the same partition
  !--------------------------------------------------

  subroutine acoustic_poro_repartitioning(elmnts_l, nbmodels, phi_material, num_material, nproc)

  use part_unstruct_par,only: nelmnts,edges_acporo_coupled,nedges_acporo_coupled,part

  implicit none
  include "constants.h"

  integer, dimension(0:NCORNERS*nelmnts-1), intent(in)  :: elmnts_l
  integer, intent(in)  :: nproc, nbmodels
  double precision, dimension(nbmodels), intent(in)  :: phi_material
  integer, dimension(1:nelmnts), intent(in)  :: num_material

  ! local parameters
  integer, dimension(:), allocatable  :: xadj_l
  integer, dimension(:), allocatable  :: adjncy_l
  logical, dimension(nbmodels)  :: is_acoustic,is_poroelastic
  integer  :: i, ier, num_edge
  integer  :: el, el_adj
  logical  :: is_repartitioned

  allocate(xadj_l(0:nelmnts))
  allocate(adjncy_l(0:MAX_NEIGHBORS*nelmnts-1))

  is_acoustic(:) = .false.
  is_poroelastic(:) = .false.

  do i = 1, nbmodels
     if (phi_material(i) >=1.d0) then
        is_acoustic(i) = .true.
     endif
     if (phi_material(i) <1.d0 .and. phi_material(i) > TINYVAL) then
        is_poroelastic(i) = .true.
     endif
  enddo

  ! determines maximum neighbors based on 2 common nodes (common edge)
  call mesh2dual_ncommonnodes(elmnts_l, 2, xadj_l, adjncy_l)

  nedges_acporo_coupled = 0
  do el = 0, nelmnts-1
     if (is_acoustic(num_material(el+1))) then
        do el_adj = xadj_l(el), xadj_l(el+1) - 1
           if (is_poroelastic(num_material(adjncy_l(el_adj)+1))) then
              nedges_acporo_coupled = nedges_acporo_coupled + 1
           endif

        enddo
     endif
  enddo

  print *, 'nedges_coupled (acoustic/poroelastic)', nedges_acporo_coupled

  allocate(edges_acporo_coupled(2,nedges_acporo_coupled),stat=ier)
  if (ier /= 0) stop 'Error allocating array edges_acporo_coupled'

  nedges_acporo_coupled = 0
  do el = 0, nelmnts-1
     if (is_acoustic(num_material(el+1))) then
        do el_adj = xadj_l(el), xadj_l(el+1) - 1
           if (is_poroelastic(num_material(adjncy_l(el_adj)+1))) then
              nedges_acporo_coupled = nedges_acporo_coupled + 1
              edges_acporo_coupled(1,nedges_acporo_coupled) = el
              edges_acporo_coupled(2,nedges_acporo_coupled) = adjncy_l(el_adj)
           endif

        enddo
     endif
  enddo

  do i = 1, nedges_acporo_coupled*nproc
     is_repartitioned = .false.
     do num_edge = 1, nedges_acporo_coupled
        if (part(edges_acporo_coupled(1,num_edge)) /= part(edges_acporo_coupled(2,num_edge))) then
           if (part(edges_acporo_coupled(1,num_edge)) < part(edges_acporo_coupled(2,num_edge))) then
              part(edges_acporo_coupled(2,num_edge)) = part(edges_acporo_coupled(1,num_edge))
           else
              part(edges_acporo_coupled(1,num_edge)) = part(edges_acporo_coupled(2,num_edge))
           endif
           is_repartitioned = .true.
        endif

     enddo
     if (.not. is_repartitioned) then
        exit
     endif
  enddo

  deallocate(xadj_l,adjncy_l)

  end subroutine acoustic_poro_repartitioning

!
!---------------------------------------------------------------------------------------
!

  !--------------------------------------------------
  ! repartitioning: coupled poroelastic/elastic elements are transferred to the same partition
  !--------------------------------------------------

  subroutine poro_elastic_repartitioning(elmnts_l, nbmodels, phi_material, num_material, nproc)

  use part_unstruct_par,only: nelmnts,nedges_elporo_coupled,edges_elporo_coupled,part

  implicit none
  include "constants.h"

  integer, dimension(0:NCORNERS*nelmnts-1), intent(in)  :: elmnts_l
  integer, intent(in)  :: nproc, nbmodels
  double precision, dimension(nbmodels), intent(in)  :: phi_material
  integer, dimension(1:nelmnts), intent(in)  :: num_material

  ! local parameters
  integer, dimension(:), allocatable  :: xadj_l
  integer, dimension(:), allocatable  :: adjncy_l
  logical, dimension(nbmodels)  :: is_elastic,is_poroelastic
  integer  :: i, ier, num_edge
  integer  :: el, el_adj
  logical  :: is_repartitioned

  allocate(xadj_l(0:nelmnts))
  allocate(adjncy_l(0:MAX_NEIGHBORS*nelmnts-1))

  is_elastic(:) = .false.
  is_poroelastic(:) = .false.

  do i = 1, nbmodels
     if (phi_material(i) < TINYVAL) then
        is_elastic(i) = .true.
     endif
     if (phi_material(i) <1.d0 .and. phi_material(i) > TINYVAL) then
        is_poroelastic(i) = .true.
     endif
  enddo

  ! determines maximum neighbors based on 2 common nodes (common edge)
  call mesh2dual_ncommonnodes(elmnts_l, 2, xadj_l, adjncy_l)

  nedges_elporo_coupled = 0
  do el = 0, nelmnts-1
     if (is_poroelastic(num_material(el+1))) then
        do el_adj = xadj_l(el), xadj_l(el+1) - 1
           if (is_elastic(num_material(adjncy_l(el_adj)+1))) then
              nedges_elporo_coupled = nedges_elporo_coupled + 1
           endif

        enddo
     endif
  enddo

  print *, 'nedges_coupled (poroelastic/elastic)', nedges_elporo_coupled

  allocate(edges_elporo_coupled(2,nedges_elporo_coupled),stat=ier)
  if (ier /= 0) stop 'Error allocating array edges_elporo_coupled'

  nedges_elporo_coupled = 0
  do el = 0, nelmnts-1
     if (is_poroelastic(num_material(el+1))) then
        do el_adj = xadj_l(el), xadj_l(el+1) - 1
           if (is_elastic(num_material(adjncy_l(el_adj)+1))) then
              nedges_elporo_coupled = nedges_elporo_coupled + 1
              edges_elporo_coupled(1,nedges_elporo_coupled) = el
              edges_elporo_coupled(2,nedges_elporo_coupled) = adjncy_l(el_adj)
           endif

        enddo
     endif
  enddo

  do i = 1, nedges_elporo_coupled*nproc
     is_repartitioned = .false.
     do num_edge = 1, nedges_elporo_coupled
        if (part(edges_elporo_coupled(1,num_edge)) /= part(edges_elporo_coupled(2,num_edge))) then
           if (part(edges_elporo_coupled(1,num_edge)) < part(edges_elporo_coupled(2,num_edge))) then
              part(edges_elporo_coupled(2,num_edge)) = part(edges_elporo_coupled(1,num_edge))
           else
              part(edges_elporo_coupled(1,num_edge)) = part(edges_elporo_coupled(2,num_edge))
           endif
           is_repartitioned = .true.
        endif

     enddo
     if (.not. is_repartitioned) then
        exit
     endif
  enddo

  deallocate(xadj_l,adjncy_l)

  end subroutine poro_elastic_repartitioning

!
!---------------------------------------------------------------------------------------
!

  !--------------------------------------------------
  ! repartitioning: coupled periodic elements are transferred to the same partition
  !--------------------------------------------------

  subroutine periodic_edges_repartitioning(elmnts_l,nnodes,nodes_coords,PERIODIC_HORIZ_DIST)

  use part_unstruct_par,only: nelmnts,part, nzread
  
  use parameter_file_par,only:periodic_BC_search_bandwidth,read_external_mesh

  implicit none
  include "constants.h"

  integer, dimension(0:NCORNERS*nelmnts-1), intent(in) :: elmnts_l

  integer :: nnodes
  double precision, dimension(2,nnodes) :: nodes_coords
  double precision :: PERIODIC_HORIZ_DIST

  ! local parameters
  logical, dimension(0:nelmnts-1) :: is_periodic

  integer :: el,el2,icorner,icorner2,num_node,num_node2,ifirst_partition_found,startPointOtherLoop

  double precision :: xtol,xtypdist
  double precision :: x,y,x2,y2
  
  double precision :: xmin, xmax, band ! Acceleration of periodic elements finding.
  logical accelerate_finding, debug_finding

! set up a local geometric tolerance by computing the typical horizontal size of an element.
! the sqrt() assumes that the geometrical model is 'not too elongated' and thus 'not too far from a square'
! and thus contains more or less the same number of points along X and Y. If this is not the case i.e.
! if the model is very elongated then this trick will work anyway because we just want to have a rough idea
! of a typical length in the mesh, even if it is not very accurate it will work anyway.
  xtypdist = (maxval(nodes_coords(1,:)) - minval(nodes_coords(1,:))) / sqrt(dble(nnodes))

! define a tolerance, small with respect to the minimum size
  xtol = 1.d-4 * xtypdist

! detect the points that are on the same horizontal line (i.e. at the same height Z)
! and that have a value of the horizontal coordinate X that differs by exactly the periodicity length;
! if so, make them all have the same global number, which will then implement periodic boundary conditions automatically.
! We select the smallest value of iglob and assign it to all the points that are the same due to periodicity,
! this way the maximum value of the ibool() array will remain as small as possible.
!
! *** IMPORTANT: this simple algorithm will be slow for large meshes because it has a cost of NGLOB^2 / 2
! (where NGLOB is the number of points per MPI slice, not of the whole mesh though). This could be
! reduced to O(NGLOB log(NGLOB)) by using a quicksort algorithm on the coordinates of the points to detect the multiples
! (as implemented in routine createnum_fast() elsewhere in the code). This could be done one day if needed instead
! of the very simple double loop below.

  print *,'start detecting points for periodic boundary conditions (the current algorithm can be slow and could be improved)...'

  is_periodic(:) = .false.
  
  accelerate_finding = .true. ! Set this to .false. to go back to the original O(N^2) search algorithm.
  if(read_external_mesh) then
    ! For now, determination of band is not working well for external meshes. Thus, use old algorithm.
    accelerate_finding = .false.
    write(*,*) "> Using external mesh, thus using raw algorithm."
  endif
  
  debug_finding = .false.
  
  if(accelerate_finding) then
    xmin = minval(nodes_coords(1,:))
    xmax = maxval(nodes_coords(1,:))
    ! Band must be greater than the maximum horizontal dx. The larger the band the better (less chances of not finding all points), but also slower. (xmax-xmin)/(number of points along x) the minimum. Other choices can be hard-coded here.
    if(abs(periodic_BC_search_bandwidth)<TINYVAL) then
      ! If left to a small value, must mean an external mesh is being used.
      ! Thus, choose a band only based of (xmin, xmax)
      band = 0.04*(xmax-xmin)
    else
      ! Else, it has been set before, so leave the band as such.
      band = periodic_BC_search_bandwidth
    endif
  endif
  
! loop on all the elements
  do el = 0, nelmnts-2 ! we stop one element before the end in order for the second loop to be OK in all cases
    
    if(debug_finding .and. mod(el, floor(nelmnts/100.))==0) then
      ! DEBUG, only here to make sure periodic points detection does not crash on huge meshes.
      if(el>0) write(*,*) "> DEBUG OF PERIODIC ELEMENTS' DETECTION: ELEMENT ", el, ", ", (100.*el/nelmnts), '% DONE.'
    endif
    
    ! Skip elements too far from left/right boundaries.
    if(accelerate_finding) then
      if(el==0) then
        write(*,*) "> Using acceleration method (skip points too far from left/right boundaries)."
        write(*,*) "> Band: ", band, " ( ",(100.*band/(xmax-xmin)),"% of horizontal span)."
      endif
      if(      nodes_coords(1,elmnts_l(NCORNERS*el) + 1) > xmin+band & ! Corner n°0 is laterally too far from left boundary.
         .and. nodes_coords(1,elmnts_l(NCORNERS*el) + 1) < xmax-band & ! Corner n°0 is laterally too far from right boundary.
        ) then
        cycle
      endif
    endif
    
    if(read_external_mesh) then
      if(el==0) then
        write(*,*) "> For some reason, when using external meshes, periodic point finding with raw algorithm crashes."
        write(*,*) "> Thus cycle through ALL points twice. Might be very long, prepare to sit there for a long time."
        write(*,*) "> Maybe go jogging or something meanwhile."
      endif
      startPointOtherLoop = 0
    else
      startPointOtherLoop = el+1
    endif
    
    do el2 = startPointOtherLoop, nelmnts-1
    !do el2 = el+1, nelmnts-1
      if (is_periodic(el2)) cycle
      ! it is sufficient to loop on the four corners to determine if this element has at least one periodic point
      do icorner = 0,NCORNERS-1
        num_node = elmnts_l(icorner + NCORNERS*el) + 1 ! the plus one is because elmnts_l() starts at zero
        x = nodes_coords(1,num_node)
        y = nodes_coords(2,num_node)
        do icorner2 = 0,NCORNERS-1
          num_node2 = elmnts_l(icorner2 + NCORNERS*el2) + 1 ! the plus one is because elmnts_l() starts at zero
          x2 = nodes_coords(1,num_node2)
          y2 = nodes_coords(2,num_node2)
          ! if the two points are at the same height Y
          !if(y>0. .and. y2>0. .and. &
          !  ((abs(x-75.)<xtol .and. abs(x2+75.)<xtol) .or. (abs(x+75.)<xtol .and. abs(x2-75.)<xtol))) &
          !    write(*,*) 'points (',x,y,') and (',x2,y2,') are considered' ! debug
          if (abs(y2 - y) < xtol) then
            !if(y>0. .and. y2>0. .and. &
            !  ((abs(x-75.)<xtol .and. abs(x2+75.)<xtol) .or. (abs(x+75.)<xtol .and. abs(x2-75.)<xtol))) &
            !    write(*,*) 'close enough on yZ' ! debug
            ! if in addition their X coordinates differ by exactly the periodicity distance
            if (abs(abs(x2 - x) - PERIODIC_HORIZ_DIST) < xtol) then
              ! then these two elements are in contact by a periodic edge
              is_periodic(el) = .true.
              is_periodic(el2) = .true.
              !write(*,*) '>> considered the same (elements ',el,' and ',el2,' are periodic)' ! debug
              goto 100
            else! debug
            !if(y>0. .and. y2>0. .and. &
            !  ((abs(x-75.)<xtol .and. abs(x2+75.)<xtol) .or. (abs(x+75.)<xtol .and. abs(x2-75.)<xtol))) &
            !      write(*,*) '>> considered NOT the same' ! debug
            endif
          endif
        enddo
      enddo
 100  continue
    enddo
  enddo
  
  print *,'done detecting points for periodic boundary conditions.'
  print *,'number of periodic elements found and grouped in the same partition: ',count(is_periodic)
  !stop 'kek' ! DEBUGff
  
  if(mod(count(is_periodic),2)/=0) then
    write(*,*) "********************************"
    write(*,*) "*            ERROR             *"
    write(*,*) "********************************"
    write(*,*) "* Odd number of periodic       *"
    write(*,*) "* elements found, this is      *"
    write(*,*) "* impossible.                  *"
    write(*,*) "********************************"
    if(accelerate_finding) then
      write(*,*) "* Acceleration method was      *"
      write(*,*) "* used, this may be the cause. *"
      write(*,*) "* Try to go to code            *"
      write(*,*) "* (meshfem2D/repartition_coupling.f90)"
      write(*,*) "* in order to increase band    *"
      write(*,*) "* size, or to deactivate       *"
      write(*,*) "* acceleration method.         *"
      write(*,*) "********************************"
    endif
    stop
  endif
  if(count(is_periodic)/=2*nzread) then
    write(*,*) "********************************"
    write(*,*) "*           WARNING            *"
    write(*,*) "********************************"
    write(*,*) "* Number of periodic elements  *"
    write(*,*) "* should (except maybe a few   *"
    write(*,*) "* specific cases) equal 2      *"
    write(*,*) "* times the number of spectral *"
    write(*,*) "* elements along Z. This is    *"
    write(*,*) "* not the case, be careful.    *"
    write(*,*) "* n_periodic = ", count(is_periodic)
    write(*,*) "* 2*nz_layer = ", (2*nzread)
    write(*,*) "********************************"
  endif

! loop on all the elements to find the first partition that contains a periodic element
  ifirst_partition_found = -1
  do el = 0, nelmnts-1
    ! TODO: maybe the acceleration trick above can be used here too, since we only look for periodic elements.
    if (is_periodic(el)) then
      ifirst_partition_found = part(el)
      exit
    endif
  enddo
  if (ifirst_partition_found < 0) stop 'error: no periodic element found, even though ADD_PERIODIC_CONDITIONS is set'

! loop on all the elements to move all periodic elements to the first partition found
  do el = 0, nelmnts-1
    ! TODO: maybe the acceleration trick above can be used here too, since we only look for periodic elements.
    if (is_periodic(el)) part(el) = ifirst_partition_found
  enddo

  end subroutine periodic_edges_repartitioning


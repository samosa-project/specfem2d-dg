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


  subroutine setup_mesh()

! creates mesh related properties, local to global mesh numbering and node locations

  use specfem_par

  implicit none
  
  ! generate the global numbering
  call setup_mesh_numbering()

  ! sets point coordinates
  call setup_mesh_coordinates()

  ! sets material properties on node points
  call setup_mesh_properties()

  ! for periodic edges
  call setup_mesh_periodic_edges()
  
  ! for DG runs
  call setup_mesh_DG()
  
  ! for acoustic forcing
  call setup_mesh_acoustic_forcing_edges()

  ! reads in external models and re-assigns material properties
  call setup_mesh_external_models()

  ! synchronizes all processes
  call synchronize_all()

  ! performs basic checks on parameters read
  all_anisotropic = .false.
  if (count(ispec_is_anisotropic(:) .eqv. .true.) == nspec) all_anisotropic = .true.

! absorbing boundaries work, but not perfect for anisotropic
!  if (all_anisotropic .and. anyabs) &
!    call exit_MPI(myrank,'Cannot put absorbing boundaries if anisotropic materials along edges')

  if (ATTENUATION_VISCOELASTIC_SOLID .and. all_anisotropic) then
    call exit_MPI(myrank,'Cannot turn attenuation on in anisotropic materials')
  endif

  ! synchronizes all processes
  call synchronize_all()

  ! global domain flags
  ! (sets global flag for all slices)
  call any_all_l(any_elastic, ELASTIC_SIMULATION)
  call any_all_l(any_poroelastic, POROELASTIC_SIMULATION)
  call any_all_l(any_acoustic, ACOUSTIC_SIMULATION)
  call any_all_l(any_gravitoacoustic, GRAVITOACOUSTIC_SIMULATION)

  ! check for acoustic
  if (ATTENUATION_VISCOELASTIC_SOLID .and. .not. ELASTIC_SIMULATION) &
    call exit_MPI(myrank,'currently cannot have attenuation if acoustic/poroelastic simulation only')

  ! sets up domain coupling, i.e. edge detection for domain coupling
  call get_coupling_edges()
  
  call setup_mesh_surface_DG_coupling()

  ! synchronizes all processes
  call synchronize_all()

  end subroutine setup_mesh

!
!-----------------------------------------------------------------------------------
!

  subroutine setup_mesh_numbering()

  use specfem_par

  implicit none

  ! local parameters
  integer :: ier
  ! to count the number of degrees of freedom
  integer :: count_nspec_acoustic_total,nspec_total,nglob_total
  integer :: nb_acoustic_DOFs,nb_elastic_DOFs
  double precision :: ratio_1DOF,ratio_2DOFs

  ! "slow and clean" or "quick and dirty" version
  if (FAST_NUMBERING) then
    call createnum_fast()
  else
    call createnum_slow()
  endif

  ! gets total numbers for all slices
  call sum_all_i(count_nspec_acoustic,count_nspec_acoustic_total)
  call sum_all_i(nspec,nspec_total)
  call sum_all_i(nglob,nglob_total)

  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'Total number of elements: ',nspec_total
    write(IMAIN,*) 'decomposed as follows:'
    write(IMAIN,*)
    write(IMAIN,*) 'Total number of elastic/visco/poro elements: ',nspec_total - count_nspec_acoustic_total
    write(IMAIN,*) 'Total number of acoustic elements: ',count_nspec_acoustic_total
    write(IMAIN,*)
#ifdef USE_MPI
    write(IMAIN,*) 'Approximate total number of grid points in the mesh'
    write(IMAIN,*) '(with a few duplicates coming from MPI buffers): ',nglob_total
#else
    write(IMAIN,*) 'Exact total number of grid points in the mesh: ',nglob_total
#endif

    ! percentage of elements with 2 degrees of freedom per point
    ratio_2DOFs = (nspec_total - count_nspec_acoustic_total) / dble(nspec_total)
    ratio_1DOF  = count_nspec_acoustic_total / dble(nspec_total)

    nb_acoustic_DOFs = nint(nglob_total*ratio_1DOF)

    ! elastic elements have two degrees of freedom per point
    nb_elastic_DOFs  = nint(nglob_total*ratio_2DOFs*2)

    if (P_SV) then
      write(IMAIN,*)
      write(IMAIN,*) 'Approximate number of acoustic degrees of freedom in the mesh: ',nb_acoustic_DOFs
      write(IMAIN,*) 'Approximate number of elastic degrees of freedom in the mesh: ',nb_elastic_DOFs
      write(IMAIN,*) '  (there are 2 degrees of freedom per point for elastic elements)'
      write(IMAIN,*)
      write(IMAIN,*) 'Approximate total number of degrees of freedom in the mesh'
      write(IMAIN,*) '(sum of the two values above): ',nb_acoustic_DOFs + nb_elastic_DOFs
      write(IMAIN,*)
      write(IMAIN,*) ' (for simplicity viscoelastic or poroelastic elements, if any,'
      write(IMAIN,*) '  are counted as elastic in the above three estimates;'
      write(IMAIN,*) '  in reality they have more degrees of freedom)'
      write(IMAIN,*)
    endif
    call flush_IMAIN()
  endif

  ! allocate temporary arrays
  allocate(integer_mask_ibool(nglob),stat=ier)
  if (ier /= 0 ) stop 'error allocating integer_mask_ibool'
  allocate(copy_ibool_ori(NGLLX,NGLLZ,nspec),stat=ier)
  if (ier /= 0 ) stop 'error allocating copy_ibool_ori'
  
  ! allocate temporary arrays
  allocate(integer_mask_ibool_bef(nglob),stat=ier)
  allocate(copy_ibool_ori_bef(NGLLX,NGLLZ,nspec),stat=ier)

  ! reduce cache misses by sorting the global numbering in the order in which it is accessed in the time loop.
  ! this speeds up the calculations significantly on modern processors
  call get_global()

  ! synchronizes all processes
  call synchronize_all()
  
  end subroutine setup_mesh_numbering

!
!-----------------------------------------------------------------------------------
!


  subroutine setup_mesh_coordinates()

  use specfem_par

  implicit none

  ! local parameters
  ! Jacobian matrix and determinant
  double precision :: xixl,xizl,gammaxl,gammazl,jacobianl
  double precision :: xi,gamma,x,z

  integer :: i,j,ispec,iglob,ier

  ! to help locate elements with a negative Jacobian using OpenDX
  logical :: found_a_negative_jacobian
  
  integer :: i2,j2,ispec2,iglob2,it_corner
  logical :: found
  
  integer :: nb

  ! allocate other global arrays
  allocate(coord(NDIM,nglob),stat=ier)
  if (ier /= 0) stop 'Error allocating coord array'

  is_corner = .false.

  ! Small trick just to find forcing element for solid part
  nb = 0

  ! sets the coordinates of the points of the global grid
  found_a_negative_jacobian = .false.
  do ispec = 1,nspec
    do j = 1,NGLLZ
      do i = 1,NGLLX
        if (AXISYM) then
          if (is_on_the_axis(ispec)) then
            xi = xiglj(i)
          else
            xi = xigll(i)
          endif
        else
          xi = xigll(i)
        endif
        gamma = zigll(j)

        call recompute_jacobian(xi,gamma,x,z,xixl,xizl,gammaxl,gammazl, &
                                jacobianl,coorg,knods,ispec,ngnod,nspec,npgeo, &
                                .false.)

        if (jacobianl <= ZERO) found_a_negative_jacobian = .true.

        ! coordinates of global nodes
        iglob = ibool(i,j,ispec)
        coord(1,iglob) = x
        coord(2,iglob) = z
        
        if(z == 0) nb = nb + 1
        
        xix(i,j,ispec) = xixl
        xiz(i,j,ispec) = xizl
        gammax(i,j,ispec) = gammaxl
        gammaz(i,j,ispec) = gammazl
        jacobian(i,j,ispec) = jacobianl
        
        ! Quick use of the loops
        if( ispec == 1 .AND. (j == 1 .AND. i == 1) &
                .OR. (j == NGLLZ .AND. i == 1) &
                .OR. (j == 1 .AND. i == NGLLX) &
                .OR. (j == NGLLZ .AND. i == NGLLX) ) is_corner(i,j) = .true.

      enddo
    enddo
  enddo
  
  ! Small trick just to find forcing element for solid part
  allocate(forcing_solid(nb, 3))
  nb_forcing_solid = nb
  nb = 0
  do ispec = 1,nspec
    do j = 1,NGLLZ
      do i = 1,NGLLX
      
        if(coord(2,ibool(i,j,ispec)) == 0) then
                nb = nb + 1
                forcing_solid(nb,1) = i
                forcing_solid(nb,2) = j
                forcing_solid(nb,3) = ispec
        endif
      
      enddo
    enddo
  enddo
  
  if(.false.) then
  ! Create neighboring elements for DG flux comm. 
   do ispec = 1,nspec
    do j = 1,NGLLZ
      do i = 1,NGLLX
      
      ! Initialization
      neighbor_DG(i,j,ispec,:)        = -1
      neighbor_DG_corner(i,j,ispec,:) = -1
      
      ! Ignore interior element nodes
      if(j > 1 .AND. j < NGLLZ .AND. i > 1 .AND. i < NGLLX) cycle
      
      iglob = ibool(i,j,ispec)
      
      found = .false.
      it_corner = 0
      ! Find neighbor by going through every mesh node
      do ispec2 = 1,nspec
      
            ! Ignore current element we are working on
            if(ispec2 == ispec) cycle
            ! If we already found the neighbor we exit the loop
            if(found) exit
            
            do j2 = 1,NGLLZ
              ! If we already found the neighbor we exit the loop
              if(found) exit
              ! If we already found the first corner neighbor we exit the loop
              if(it_corner == 1) exit
              
              do i2 = 1,NGLLX
              
              ! Ignore interior element nodes
              if(j2 > 1 .AND. j2 < NGLLZ .AND. i2 > 1 .AND. i2 < NGLLX) cycle
              
              iglob2 = ibool(i2,j2,ispec2)
              
              ! If same coordinates => Neighbor found
              !if(( coord(1,iglob) == coord(1,iglob2) .AND. coord(2,iglob) == coord(2,iglob2) .AND. &
              if( iglob == iglob2 .AND. &
              ! If must not be the opposite corner
                 (i+i2+j+j2 /= (NGLLX + NGLLZ + 2) .OR. .not. is_corner(i,j)) ) then
              
                if(it_corner == 0) then
                neighbor_DG(i,j,ispec,1) = i2
                neighbor_DG(i,j,ispec,2) = j2
                neighbor_DG(i,j,ispec,3) = ispec2
                
                !WRITE(*,*) "CLASSICAL FINDING 1", i, j, ispec, coord(:,ibool(i,j,ispec))
                !WRITE(*,*) "CLASSICAL FINDING 2", i2, j2, ispec2, coord(:,ibool(i2,j2,ispec2))
                !WRITE(*,*) "*************"
                !if(coord(2,ibool(i,j,ispec)) /= coord(2,ibool(i2,j2,ispec2)) .OR. &
                !     coord(1,ibool(i,j,ispec)) /= coord(1,ibool(i2,j2,ispec2))   ) then
                !     WRITE(*,*) "ERROR FINDING 1", i, j, ispec, coord(:,ibool(i,j,ispec))
                !WRITE(*,*) "ERROR FINDING 2", i2, j2, ispec2, coord(:,ibool(i2,j2,ispec2))
                !WRITE(*,*) "*************"
                !endif
                endif
                
                ! If at corner => needs two neighbors
                if( is_corner(i,j) .AND. it_corner == 0) then
                
                        it_corner = 1
                        
                ! Second corner neighbor
                elseif(it_corner == 2) then
                
                        found = .true.
                        neighbor_DG_corner(i,j,ispec,1) = i2
                        neighbor_DG_corner(i,j,ispec,2) = j2
                        neighbor_DG_corner(i,j,ispec,3) = ispec2
                        if(neighbor_DG(i,j,ispec,3) == neighbor_DG_corner(i,j,ispec,3)) then
                        WRITE(*,*) "CORNER FINDING 1", i, j, ispec, coord(:,ibool(i,j,ispec))
                        WRITE(*,*) "CORNER FINDING 2", i2, j2, ispec2, coord(:,ibool(i2,j2,ispec2))
                        WRITE(*,*) "CORNER FINDING 3", &
                        neighbor_DG(i,j,ispec,1), neighbor_DG(i,j,ispec,2), neighbor_DG(i,j,ispec,3)
                        WRITE(*,*) "*************"
                        stop
                        endif
                        it_corner = 0
                endif
                
                ! Normal interior point
                if(.not. is_corner(i,j)) found = .true.
                
                ! If found a neighbor leave the current loop
                exit
              endif
              
              enddo
            enddo
            
            ! If we went through all the nodes of the previous element where we found a neighbor
            ! We notify that the research should continue normally
            if(it_corner == 1) it_corner = 2
            
      enddo
      
      enddo
    enddo
    
  enddo
  
  !stop
  endif
! create an OpenDX file containing all the negative elements displayed in red, if any
! this allows users to locate problems in a mesh based on the OpenDX file created at the second iteration
! do not create OpenDX files if no negative Jacobian has been found, or if we are running in parallel
! (because writing OpenDX routines is much easier in serial)
  if (found_a_negative_jacobian .and. NPROC == 1) then
    call save_openDX_jacobian(nspec,npgeo,ngnod,knods,coorg,xigll,zigll,AXISYM,is_on_the_axis,xiglj)
  endif

  ! stop the code at the first negative element found, because such a mesh cannot be computed
  if (found_a_negative_jacobian) then
    do ispec = 1,nspec
      do j = 1,NGLLZ
        do i = 1,NGLLX
          if (AXISYM) then
            if (is_on_the_axis(ispec)) then
              xi = xiglj(i)
            else
              xi = xigll(i)
            endif
          else
            xi = xigll(i)
          endif
          gamma = zigll(j)

          call recompute_jacobian(xi,gamma,x,z,xixl,xizl,gammaxl,gammazl, &
                                  jacobianl,coorg,knods,ispec,ngnod,nspec,npgeo, &
                                  .true.)
        enddo
      enddo
    enddo
  endif

  ! synchronizes all processes
  call synchronize_all()

  end subroutine setup_mesh_coordinates

!
!-----------------------------------------------------------------------------------
!


  subroutine setup_mesh_DG()

  use specfem_par

  implicit none

  ! local parameters

  integer :: i,j,ispec,iglob
  
  integer :: i2,j2,ispec2,iglob2,it_corner
  logical :: found

  is_corner = .false.
  
  ! Create neighboring elements for DG flux comm. 
   do ispec = 1,nspec
   
    do j = 1,NGLLZ
      do i = 1,NGLLX
      
      ! Setup is_corner table link
        if( ispec == 1 .AND. (j == 1 .AND. i == 1) &
                .OR. (j == NGLLZ .AND. i == 1) &
                .OR. (j == 1 .AND. i == NGLLX) &
                .OR. (j == NGLLZ .AND. i == NGLLX) ) is_corner(i,j) = .true.
      
      ! Initialization
      neighbor_DG(i,j,ispec,:)        = -1
      neighbor_DG_corner(i,j,ispec,:) = -1
      
      ! Skip non acoustic (thus non fluid) elements
      if (.not. ispec_is_acoustic(ispec)) cycle
      
      ! Ignore interior element nodes
      if(j > 1 .AND. j < NGLLZ .AND. i > 1 .AND. i < NGLLX) cycle
      
      ! If already found neighbor => cycle
      if(.not. is_corner(i,j) .AND. neighbor_DG(i,j,ispec,1) > -1) cycle
      
      iglob = ibool(i,j,ispec)
      
      found = .false.
      it_corner = 0
      ! Find neighbor by going through every mesh node
      do ispec2 = 1,nspec
      
            ! Skip non acoustic (thus non fluid) elements
            if (.not. ispec_is_acoustic(ispec2)) cycle
      
            ! Ignore current element we are working on
            if(ispec2 == ispec) cycle
            ! If we already found the neighbor we exit the loop
            if(found) exit
            
            do j2 = 1,NGLLZ
            
              ! If we already found the neighbor we exit the loop
              if(found) exit
              
              ! If we already found the first corner neighbor we exit the loop
              if(it_corner == 1) exit
              
              do i2 = 1,NGLLX
              
              ! Ignore interior element nodes
              if(j2 > 1 .AND. j2 < NGLLZ .AND. i2 > 1 .AND. i2 < NGLLX) cycle
              
              iglob2 = ibool(i2,j2,ispec2)
              
              ! If same coordinates => Neighbor found
              !if(( coord(1,iglob) == coord(1,iglob2) .AND. coord(2,iglob) == coord(2,iglob2) .AND. &
              if( iglob == iglob2 .AND. &
              ! If must not be the opposite corner
                 (i+i2+j+j2 /= (NGLLX + NGLLZ + 2) .OR. .not. is_corner(i,j)) ) then
              
                if(it_corner == 0) then
                neighbor_DG(i,j,ispec,1) = i2
                neighbor_DG(i,j,ispec,2) = j2
                neighbor_DG(i,j,ispec,3) = ispec2
                
                !neighbor_DG(i2,j2,ispec2,1) = i
                !neighbor_DG(i2,j2,ispec2,2) = j
                !neighbor_DG(i2,j2,ispec2,3) = ispec
                !if(coord(1,ibool_before_perio(i,j,ispec)) == 10.) then
                !WRITE(*,*) "CLASSICAL FINDING 1", i, j, ispec, coord(:,ibool_before_perio(i,j,ispec))
                !WRITE(*,*) "CLASSICAL FINDING 2", i2, j2, ispec2, coord(:,ibool_before_perio(i2,j2,ispec2))
                !WRITE(*,*) "*************"
                !endif
                endif
                
                ! If at corner => needs two neighbors
                if( is_corner(i,j) .AND. it_corner == 0) then
                
                        it_corner = 1
                        
                ! Second corner neighbor
                elseif(it_corner == 2) then
                
                        found = .true.
                        neighbor_DG_corner(i,j,ispec,1) = i2
                        neighbor_DG_corner(i,j,ispec,2) = j2
                        neighbor_DG_corner(i,j,ispec,3) = ispec2
                        !if(coord(1,ibool_before_perio(i,j,ispec)) == 14000.) then
                        !WRITE(*,*) "CORNER FINDING 1", i, j, ispec, coord(:,ibool_before_perio(i,j,ispec))
                        !WRITE(*,*) "CORNER FINDING 2", i2, j2, ispec2, coord(:,ibool_before_perio(i2,j2,ispec2))
                        !WRITE(*,*) myrank, "neighbor_DG_corner(i,j,ispec,1)", neighbor_DG_corner(i,j,ispec,:), &
                        !        neighbor_DG(i,j,ispec,:)
                        !WRITE(*,*) "*************"
                        !endif
                        it_corner = 0
                endif
                
                ! Normal interior point
                if(.not. is_corner(i,j)) found = .true.
                
                ! If found a neighbor leave the current loop
                exit
              endif
              
              enddo
            enddo
            
            ! If we went through all the nodes of the previous element where we found a neighbor
            ! We notify that the research should continue normally
            if(it_corner == 1) it_corner = 2
            
      enddo
      
      enddo
    enddo
    
  enddo
  
  !stop
  
  call setup_mesh_surface_DG()
  
  end subroutine setup_mesh_DG
  
!
!-----------------------------------------------------------------------------------
!


  subroutine setup_mesh_surface_DG()

  use constants,only: CUSTOM_REAL,NGLLX,NGLLZ
  use specfem_par,only: neighbor_DG, neighbor_DG_corner, is_corner, &
        ACOUSTIC_FORCING, numacforcing, ispec_is_acoustic_forcing, nelem_acforcing, &
        acoustic_surface, ispec_is_acoustic_surface, ispec_is_acoustic_surface_corner, &
        nelem_acoustic_surface!, &
        !coord, ibool_before_perio,myrank, nspec

  implicit none

  ! local parameters

  integer :: i, j, ispec, numelem
  
  !!!!!!!!!!!!!!!!!!!!!!!
   ! Find forcing elements
   ispec_is_acoustic_forcing = .false.
   if(ACOUSTIC_FORCING) then
   do ispec = 1,nelem_acforcing
   
        numelem = numacforcing(ispec)
        do i = 1,NGLLX
        do j = 1,NGLLZ
                if(neighbor_DG(i,j,numelem,3) == -1) &
                        ispec_is_acoustic_forcing(i,j,numelem) = .true.
        enddo
        enddo
   
   enddo
   endif
   
   !!!!!!!!!!!!!!!!!!!!!!!
   ! Find surface elements
   ispec_is_acoustic_surface = .false.
   ispec_is_acoustic_surface_corner = .false.
   do ispec = 1, nelem_acoustic_surface
        
        numelem = acoustic_surface(1,ispec)
        !WRITE(*,*) "TEST numelem", numelem
        do i = 1,NGLLX
        do j = 1,NGLLZ
        
                if(neighbor_DG(i,j,numelem,3) == -1 .AND. &
                (i == 1 .OR. i == NGLLX .OR. j == 1 .OR. j == NGLLZ)) &!then
                        ispec_is_acoustic_surface(i,j,numelem) = .true.
                        !WRITE(*,*) myrank,"BOUNDARY", i, j, numelem, coord(:,ibool_before_perio(i,j,ispec))
                !endif
                if(is_corner(i,j) .AND. neighbor_DG_corner(i,j,numelem,3) == -1) &
                        ispec_is_acoustic_surface_corner(i,j,numelem) = .true.
        enddo
        enddo
        
   enddo
  
  !stop
  !endif

  end subroutine setup_mesh_surface_DG

!
!-----------------------------------------------------------------------------------
!

subroutine setup_mesh_surface_DG_coupling()

  use constants,only: CUSTOM_REAL,NGLLX,NGLLZ
  use specfem_par,only: num_fluid_solid_edges, &
        fluid_solid_acoustic_ispec, fluid_solid_elastic_ispec, ivalue, jvalue, &
        ispec_is_acoustic_coupling_el, fluid_solid_elastic_iedge, &
        fluid_solid_acoustic_iedge, ivalue_inverse, jvalue_inverse, &
        ispec_is_acoustic_surface, ispec_is_acoustic_surface_corner, is_corner, &
        coord_interface, ibool_before_perio, coord
        !, coord, ibool_before_perio

  implicit none

  ! local parameters

  integer :: i, j, i_el, j_el, ispec, ipoin1D, ispec_elastic, iedge_elastic, &
        ispec_acoustic, iedge_acoustic
   
   !!!!!!!!!!!!!!!!!!!!!!!
   ! Find coupling elements
   coord_interface = 10**8
   if(num_fluid_solid_edges < 1) coord_interface = 0
   ispec_is_acoustic_coupling_el = -1
   do ispec = 1,num_fluid_solid_edges

        !do i = 1,NGLLX
        !do j = 1,NGLLZ
        !        if(neighbor_DG(i,j,numelem,3) == -1) &
        !                ispec_is_acoustic_coupling(i,j,numelem) = .true.
        !enddo
        !enddo
        
        ! get the edge of the acoustic element
        ispec_acoustic = fluid_solid_acoustic_ispec(ispec)
        iedge_acoustic = fluid_solid_acoustic_iedge(ispec)
        !WRITE(*,*) "ispec_acoustic", ispec_acoustic
        ! get the corresponding edge of the elastic element
        ispec_elastic = fluid_solid_elastic_ispec(ispec)
        iedge_elastic = fluid_solid_elastic_iedge(ispec)

        ! implement 1D coupling along the edge
        do ipoin1D = 1,NGLLX

        ! get point values for the acoustic side
        i = ivalue(ipoin1D,iedge_acoustic)
        j = jvalue(ipoin1D,iedge_acoustic)
        
        ispec_is_acoustic_surface(i,j,ispec_acoustic) = .true.
        if(is_corner(i,j)) &
                ispec_is_acoustic_surface_corner(i,j,ispec_acoustic) = .true.

        ! get point values for the elastic side
        i_el = ivalue_inverse(ipoin1D,iedge_elastic)
        j_el = jvalue_inverse(ipoin1D,iedge_elastic)
        
        ispec_is_acoustic_coupling_el(i,j,ispec_acoustic,1) = i_el
        ispec_is_acoustic_coupling_el(i,j,ispec_acoustic,2) = j_el
        ispec_is_acoustic_coupling_el(i,j,ispec_acoustic,3) = ispec_elastic
        
        if(coord(2, ibool_before_perio(i,j,ispec_acoustic)) < coord_interface) &
        coord_interface = coord(2, ibool_before_perio(i,j,ispec_acoustic))
        
        !WRITE(*,*) ">>>> ", i, j,ispec_acoustic, i_el, j_el, ispec_elastic, &
        !        ispec_is_acoustic_coupling_el(i,j,ispec_acoustic,:)
        
        enddo

   enddo
   
   !stop

  end subroutine setup_mesh_surface_DG_coupling

!
!-----------------------------------------------------------------------------------
!


  subroutine setup_mesh_properties()

#ifdef USE_MPI
  use mpi
#endif

  use specfem_par
  use specfem_par_movie

  implicit none

  ! local parameters
  double precision :: xmin,xmax,zmin,zmax
  double precision :: xmin_local,xmax_local,zmin_local,zmax_local
  integer :: i,n

  ! determines mesh dimensions
  xmin_local = minval(coord(1,:))
  xmax_local = maxval(coord(1,:))
  zmin_local = minval(coord(2,:))
  zmax_local = maxval(coord(2,:))

  ! collect min/max
  call min_all_all_dp(xmin_local, xmin)
  call max_all_all_dp(xmax_local, xmax)
  call min_all_all_dp(zmin_local, zmin)
  call max_all_all_dp(zmax_local, zmax)

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'Xmin,Xmax of the whole mesh = ',xmin,xmax
    write(IMAIN,*) 'Zmin,Zmax of the whole mesh = ',zmin,zmax
    write(IMAIN,*)
  endif

  ! checks that no source is located outside the mesh
  if (myrank == 0) then
    do i = 1,NSOURCES
      if (x_source(i) < xmin) stop 'error: at least one source has x < xmin of the mesh'
      if (x_source(i) > xmax) stop 'error: at least one source has x > xmax of the mesh'

      if (z_source(i) < zmin) stop 'error: at least one source has z < zmin of the mesh'
      if (z_source(i) > zmax) stop 'error: at least one source has z > zmax of the mesh'
    enddo
  endif

  ! saves the grid of points in a file
  if (output_grid_ASCII .and. myrank == 0) then
     write(IMAIN,*)
     write(IMAIN,*) 'Saving the grid in an ASCII text file...'
     write(IMAIN,*)
     open(unit=55,file='OUTPUT_FILES/ASCII_dump_of_grid_points.txt',status='unknown')
     write(55,*) nglob
     do n = 1,nglob
        write(55,*) (coord(i,n), i = 1,NDIM)
     enddo
     close(55)
  endif

  ! plots the GLL mesh in a Gnuplot file
  if (output_grid_Gnuplot .and. myrank == 0) then
    call plot_gll()
  endif

  ! synchronizes all processes
  call synchronize_all()

  end subroutine setup_mesh_properties

!
!-----------------------------------------------------------------------------------
!

  subroutine setup_mesh_periodic_edges()

  use constants,only: IMAIN,NGLLX,NGLLZ,HUGEVAL
  use specfem_par

  implicit none

  ! local parameters
  integer :: ispec,i,j,iglob,iglob2,ier
  double precision :: xmaxval,xminval,ymaxval,yminval,xtol,xtypdist
  integer :: counter, nb, iglob_nb, iglob2_nb
  logical :: find
  integer, dimension(nspec*NGLLX*NGLLZ,3) :: iglob_perio

  ! Save ibool for DG normal computations
  ibool_before_perio = ibool

! allocate an array to make sure that an acoustic free surface is not enforced on periodic edges
  allocate(this_ibool_is_a_periodic_edge(NGLOB),stat=ier)
  if (ier /= 0) stop 'Error allocating periodic edge array'

  this_ibool_is_a_periodic_edge(:) = .false.

! periodic conditions: detect common points between left and right edges and replace one of them with the other
  if (ADD_PERIODIC_CONDITIONS) then
    ! user output
    if (myrank == 0) then
      write(IMAIN,*)
      write(IMAIN,*) 'implementing periodic boundary conditions'
      write(IMAIN,*) 'in the horizontal direction with a periodicity distance of ',PERIODIC_HORIZ_DIST,' m'
      if (PERIODIC_HORIZ_DIST <= 0.d0) stop 'PERIODIC_HORIZ_DIST should be greater than zero when using ADD_PERIODIC_CONDITIONS'
      write(IMAIN,*)
      write(IMAIN,*) '*****************************************************************'
      write(IMAIN,*) '*****************************************************************'
      write(IMAIN,*) '**** BEWARE: because of periodic conditions, values computed ****'
      write(IMAIN,*) '****         by check_grid() below will not be reliable       ****'
      write(IMAIN,*) '*****************************************************************'
      write(IMAIN,*) '*****************************************************************'
      write(IMAIN,*)
    endif

    ! set up a local geometric tolerance
    xtypdist = +HUGEVAL

    do ispec = 1,nspec

      xminval = +HUGEVAL
      yminval = +HUGEVAL
      xmaxval = -HUGEVAL
      ymaxval = -HUGEVAL

      nb = 0
      ! only loop on the four corners of each element to get a typical size
      do j = 1,NGLLZ,NGLLZ-1
        do i = 1,NGLLX,NGLLX-1
          iglob = ibool(i,j,ispec)
          xmaxval = max(coord(1,iglob),xmaxval)
          xminval = min(coord(1,iglob),xminval)
          ymaxval = max(coord(2,iglob),ymaxval)
          yminval = min(coord(2,iglob),yminval)
        enddo
      enddo

      ! compute the minimum typical "size" of an element in the mesh
      xtypdist = min(xtypdist,xmaxval-xminval)
      xtypdist = min(xtypdist,ymaxval-yminval)

    enddo
    
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
    if (myrank == 0) then
      write(IMAIN,*) 'start detecting points for periodic boundary conditions '// &
                     '(the current algorithm can be slow and could be improved)...'
    endif

    if(.true.) then
    counter = 0
    do iglob = 1,NGLOB-1
      !find = .false.
      do iglob2 = iglob + 1,NGLOB
        ! check if the two points have the exact same Z coordinate
        if (abs(coord(2,iglob2) - coord(2,iglob)) < xtol) then
         ! find = .true.
          ! if so, check if their X coordinate differs by exactly the periodicity distance
          if (abs(abs(coord(1,iglob2) - coord(1,iglob)) - PERIODIC_HORIZ_DIST) < xtol) then
            ! if so, they are the same point, thus replace the highest value of ibool with the lowest
            ! to make them the same global point and thus implement periodicity automatically
            counter = counter + 1
            this_ibool_is_a_periodic_edge(iglob) = .true.
            this_ibool_is_a_periodic_edge(iglob2) = .true.
            do ispec = 1,nspec
              do j = 1,NGLLZ
                do i = 1,NGLLX
                  if (ibool(i,j,ispec) == iglob2) ibool(i,j,ispec) = iglob
                enddo
              enddo
            enddo
          endif
        endif
      !  if(find) exit
      enddo
     ! if(find) exit
    enddo
    
    else
    
     do ispec = 1,nspec

      nb = 0
      ! only loop on the four corners of each element to get a typical size
      do j = 1,NGLLZ
        do i = 1,NGLLX
          iglob = ibool(i,j,ispec)
          
          if(abs(coord(1,iglob)) < xtol .OR. &
                abs(coord(1,iglob) - PERIODIC_HORIZ_DIST) < xtol) then
                nb = nb + 1
                iglob_perio(nb,1) = i
                iglob_perio(nb,2) = j
                iglob_perio(nb,3) = ispec
                !iglob_perio(nb,4) = ispec
          endif
          
        enddo
      enddo

    enddo
    
    counter = 0
    do iglob_nb = 1,nb-1
        find = .false.
        iglob = ibool(iglob_perio(iglob_nb,1), iglob_perio(iglob_nb,2), iglob_perio(iglob_nb,3))
        do iglob2_nb = iglob_nb + 1,nb
        iglob2 = ibool(iglob_perio(iglob2_nb,1), iglob_perio(iglob2_nb,2), iglob_perio(iglob2_nb,3))
        ! check if the two points have the exact same Z coordinate
        if (abs(coord(2,iglob2) - coord(2,iglob)) < xtol) then
          find = .true.
          ! if so, check if their X coordinate differs by exactly the periodicity distance
          if (abs(abs(coord(1,iglob2) - coord(1,iglob)) - PERIODIC_HORIZ_DIST) < xtol) then
            ! if so, they are the same point, thus replace the highest value of ibool with the lowest
            ! to make them the same global point and thus implement periodicity automatically
            counter = counter + 1
            this_ibool_is_a_periodic_edge(iglob) = .true.
            this_ibool_is_a_periodic_edge(iglob2) = .true.
            ibool(iglob_perio(iglob2_nb,1), iglob_perio(iglob2_nb,2), iglob_perio(iglob2_nb,3)) = iglob
          endif
        endif
        if(find) exit
      enddo
      if(find) exit
    enddo
    
    endif

    if (myrank == 0) write(IMAIN,*) 'done detecting points for periodic boundary conditions.'

    if (counter > 0) write(IMAIN,*) 'implemented periodic conditions on ',counter,' grid points on proc ',myrank

  endif ! of if (ADD_PERIODIC_CONDITIONS)

  end subroutine setup_mesh_periodic_edges

!
!-----------------------------------------------------------------------------------
!

  subroutine setup_mesh_acoustic_forcing_edges()

! acoustic forcing edge detection

  use specfem_par

  implicit none

  ! local parameters
  integer :: ipoin1D

  ! acoustic forcing edge detection
  ! the elements forming an edge are already known (computed in meshfem2D),
  ! the common nodes forming the edge are computed here
  if (ACOUSTIC_FORCING) then

    ! user output
    if (myrank == 0) then
      write(IMAIN,*)
      write(IMAIN,*) 'Acoustic forcing simulation'
      write(IMAIN,*)
      write(IMAIN,*) 'Beginning of acoustic forcing edge detection'
      call flush_IMAIN()
    endif

    ! define i and j points for each edge
    do ipoin1D = 1,NGLLX

      ivalue(ipoin1D,IBOTTOM) = NGLLX - ipoin1D + 1
      ivalue_inverse(ipoin1D,IBOTTOM) = ipoin1D
      jvalue(ipoin1D,IBOTTOM) = NGLLZ
      jvalue_inverse(ipoin1D,IBOTTOM) = NGLLZ

      ivalue(ipoin1D,IRIGHT) = 1
      ivalue_inverse(ipoin1D,IRIGHT) = 1
      jvalue(ipoin1D,IRIGHT) = NGLLZ - ipoin1D + 1
      jvalue_inverse(ipoin1D,IRIGHT) = ipoin1D

      ivalue(ipoin1D,ITOP) = ipoin1D
      ivalue_inverse(ipoin1D,ITOP) = NGLLX - ipoin1D + 1
      jvalue(ipoin1D,ITOP) = 1
      jvalue_inverse(ipoin1D,ITOP) = 1

      ivalue(ipoin1D,ILEFT) = NGLLX
      ivalue_inverse(ipoin1D,ILEFT) = NGLLX
      jvalue(ipoin1D,ILEFT) = ipoin1D
      jvalue_inverse(ipoin1D,ILEFT) = NGLLZ - ipoin1D + 1

    enddo

  endif ! if (ACOUSTIC_FORCING)

  ! synchronizes all processes
  call synchronize_all()

  end subroutine setup_mesh_acoustic_forcing_edges


!
!-----------------------------------------------------------------------------------
!

  subroutine setup_mesh_external_models()

! external models

  use specfem_par

  implicit none

  ! local parameters
  integer :: nspec_ext,ier

  ! allocates material arrays
  if (assign_external_model) then
    nspec_ext = nspec
  else
    ! dummy allocations
    nspec_ext = 1
  endif

  allocate(vpext(NGLLX,NGLLZ,nspec_ext), &
           vsext(NGLLX,NGLLZ,nspec_ext), &
           rhoext(NGLLX,NGLLZ,nspec_ext), &
           gravityext(NGLLX,NGLLZ,nspec_ext), &
           Nsqext(NGLLX,NGLLZ,nspec_ext), &
           QKappa_attenuationext(NGLLX,NGLLZ,nspec_ext), &
           Qmu_attenuationext(NGLLX,NGLLZ,nspec_ext), &
           c11ext(NGLLX,NGLLZ,nspec_ext), &
           c13ext(NGLLX,NGLLZ,nspec_ext), &
           c15ext(NGLLX,NGLLZ,nspec_ext), &
           c33ext(NGLLX,NGLLZ,nspec_ext), &
           c35ext(NGLLX,NGLLZ,nspec_ext), &
           c55ext(NGLLX,NGLLZ,nspec_ext), &
           c12ext(NGLLX,NGLLZ,nspec_ext), &
           c23ext(NGLLX,NGLLZ,nspec_ext), &
           c25ext(NGLLX,NGLLZ,nspec_ext), &
           c22ext(NGLLX,NGLLZ,nspec_ext),stat=ier)
  ! MODIF DG
  allocate(gammaext(NGLLX,NGLLZ,nspec_ext), &
           windxext(NGLLX,NGLLZ,nspec_ext), &
           windzext(NGLLX,NGLLZ,nspec_ext), &
           muext(NGLLX,NGLLZ,nspec_ext), &
           etaext(NGLLX,NGLLZ,nspec_ext), &
           pext_DG(NGLLX,NGLLZ,nspec_ext), &
           gammaext_DG(nglob_DG), stat=ier) 
           
  if (ier /= 0) stop 'Error allocating external model arrays'

  ! reads in external models
  if (assign_external_model) then
    if (myrank == 0) then
      write(IMAIN,*) 'Assigning an external velocity and density model'
      call flush_IMAIN()
    endif
    call read_external_model()
  endif

  ! synchronizes all processes
  call synchronize_all()

  end subroutine setup_mesh_external_models


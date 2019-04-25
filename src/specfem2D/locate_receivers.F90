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

!----
!---- locate_receivers finds the correct position of the receivers
!----

  subroutine locate_receivers(ibool,coord,nspec,nglob,xigll,zigll, &
                              nrec,nrecloc,recloc,which_proc_receiver,NPROC,myrank, &
                              st_xval,st_zval,ispec_selected_rec, &
                              xi_receiver,gamma_receiver,station_name,network_name, &
                              x_source,z_source, &
                              coorg,knods,ngnod,npgeo, &
                              x_final_receiver, z_final_receiver)

  use constants,only: NDIM,NGLLX,NGLLZ,MAX_LENGTH_STATION_NAME,MAX_LENGTH_NETWORK_NAME, &
    IMAIN,HUGEVAL,TINYVAL,NUM_ITER

  use specfem_par, only : AXISYM,is_on_the_axis,xiglj,ispec_is_acoustic,USE_TRICK_FOR_BETTER_PRESSURE&
                          , ispec_is_elastic, ispec_is_acoustic_dg

#ifdef USE_MPI
  use mpi
#endif

  implicit none

  integer :: nrec,nspec,nglob,ngnod,npgeo
  integer, intent(in)  :: NPROC, myrank

  integer :: knods(ngnod,nspec)
  double precision :: coorg(NDIM,npgeo)

  integer, dimension(NGLLX,NGLLZ,nspec) :: ibool

! array containing coordinates of the points
  double precision :: coord(NDIM,nglob)

  integer :: irec,i,j,ispec,iglob,iter_loop,ix_initial_guess,iz_initial_guess

  double precision :: x_source,z_source,dist_squared,stele,stbur
  double precision, dimension(nrec)  :: distance_receiver
  double precision :: xi,gamma,dx,dz,dxi,dgamma

! Gauss-Lobatto-Legendre points of integration
  double precision :: xigll(NGLLX)
  double precision :: zigll(NGLLZ)

  double precision :: x,z,xix,xiz,gammax,gammaz,jacobian

! use dynamic allocation
  double precision :: distmin_squared
  double precision, dimension(:), allocatable :: final_distance

! receiver information
  integer  :: nrecloc
  integer, dimension(nrec) :: ispec_selected_rec, recloc
  double precision, dimension(nrec) :: xi_receiver,gamma_receiver

! station information for writing the seismograms
  character(len=MAX_LENGTH_STATION_NAME), dimension(nrec) :: station_name
  character(len=MAX_LENGTH_NETWORK_NAME), dimension(nrec) :: network_name

  double precision, dimension(nrec) :: st_xval,st_zval

! tangential detection
  double precision, dimension(nrec)  :: x_final_receiver, z_final_receiver

  double precision, dimension(nrec,NPROC)  :: gather_final_distance
  double precision, dimension(nrec,NPROC)  :: gather_xi_receiver, gather_gamma_receiver

  integer, dimension(nrec), intent(inout)  :: which_proc_receiver
  integer, dimension(:,:), allocatable  :: gather_ispec_selected_rec
  integer  :: ier

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) '********************'
    write(IMAIN,*) ' locating receivers'
    write(IMAIN,*) '********************'
    write(IMAIN,*)
    write(IMAIN,*) 'reading receiver information from the DATA/STATIONS file'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  open(unit= 1,file='DATA/STATIONS',status='old',action='read')

! allocate memory for arrays using number of stations
  allocate(final_distance(nrec))

! loop on all the stations
  do irec= 1,nrec

    ! set distance to huge initial value
    distmin_squared = HUGEVAL

    read(1,*) station_name(irec),network_name(irec),st_xval(irec),st_zval(irec),stele,stbur

    ! check that station is not buried, burial is not implemented in current code
    if (abs(stbur) > TINYVAL) call exit_MPI(myrank,'stations with non-zero burial not implemented yet')

    ! compute distance between source and receiver
    distance_receiver(irec) = sqrt((st_zval(irec)-z_source)**2 + (st_xval(irec)-x_source)**2)

    do ispec= 1,nspec
      
      ! The test below fails if stations are demanded too close to the fluid-solid interface in DG simulations. See ranting about that issue below, near user output.
      ! A classical test would be to check if the station is within some element using the fact that they are convex quadrilaterals. Principle, for each element:
      ! 1) Compute area of the quadrilateral (maybe it is already stored somewhere though).
      ! 2) Compute sum of the areas of the four triangles formed by each side of the quadrilateral and the station.
      ! 3a) If the area of the quadrilateral and the sum of the four triangles' areas are equal (within some tolerance, typically TINYVAL), then the station is in this element. Being on the edge of the element will also verify this condition, though very close to the tolerance.
      ! 3b) If the area of the quadrilateral is SMALLER than the sum of the four triangles' areas, then the station is outside the element.
      ! 3c) [area of the quadrilateral > sum of the four triangles' areas] cannot happen.
      ! See e.g. https://i.stack.imgur.com/zFEPn.jpg for an illustration.
      ! 
      ! In fine, we leave the original code because we do not have time to thouroughly test our new method.
      
      ! loop only on points inside the element
      ! exclude edges to ensure this point is not shared with other elements
      do j = 2,NGLLZ-1
        do i = 2,NGLLX-1

          iglob = ibool(i,j,ispec)

          !  we compare squared distances instead of distances themselves to significantly speed up calculations
          dist_squared = (st_xval(irec)-dble(coord(1,iglob)))**2 + (st_zval(irec)-dble(coord(2,iglob)))**2

          ! keep this point if it is closer to the receiver
          if (dist_squared < distmin_squared) then
            distmin_squared = dist_squared
            ispec_selected_rec(irec) = ispec
            ix_initial_guess = i
            iz_initial_guess = j
          endif

        enddo
      enddo

    ! end of loop on all the spectral elements
    enddo


    ! ****************************************
    ! find the best (xi,gamma) for each receiver
    ! ****************************************
    ! use initial guess in xi and gamma

    if (AXISYM) then ! TODO
      if (is_on_the_axis(ispec_selected_rec(irec))) then
        xi = xiglj(ix_initial_guess)
      else
        xi = xigll(ix_initial_guess)
      endif
    else
      xi = xigll(ix_initial_guess)
    endif
    gamma = zigll(iz_initial_guess)

    ! iterate to solve the non linear system
    do iter_loop = 1,NUM_ITER
      ! compute coordinates of the new point and derivatives dxi/dx, dxi/dz
      call recompute_jacobian(xi,gamma,x,z,xix,xiz,gammax,gammaz,jacobian, &
                  coorg,knods,ispec_selected_rec(irec),ngnod,nspec,npgeo, &
                  .true.)

      ! compute distance to target location
      dx = - (x - st_xval(irec))
      dz = - (z - st_zval(irec))

      ! compute increments
      dxi  = xix*dx + xiz*dz
      dgamma = gammax*dx + gammaz*dz

      ! update values
      xi = xi + dxi
      gamma = gamma + dgamma

      ! impose that we stay in that element
      ! (useful if user gives a receiver outside the mesh for instance)
      ! we can go slightly outside the [1,1] segment since with finite elements
      ! the polynomial solution is defined everywhere
      ! this can be useful for convergence of itertive scheme with distorted elements
      if (xi > 1.10d0) xi = 1.10d0
      if (xi < -1.10d0) xi = -1.10d0
      if (gamma > 1.10d0) gamma = 1.10d0
      if (gamma < -1.10d0) gamma = -1.10d0

    ! end of non linear iterations
    enddo

    ! compute final coordinates of point found
    call recompute_jacobian(xi,gamma,x,z,xix,xiz,gammax,gammaz,jacobian, &
                coorg,knods,ispec_selected_rec(irec),ngnod,nspec,npgeo, &
                .true.)

    ! store xi,gamma of point found
    xi_receiver(irec) = xi
    gamma_receiver(irec) = gamma

    ! compute final distance between asked and found
    final_distance(irec) = sqrt((st_xval(irec)-x)**2 + (st_zval(irec)-z)**2)

    x_final_receiver(irec) = x
    z_final_receiver(irec) = z

  enddo

  ! close receiver file
  close(1)

  ! select one mesh slice for each receiver
  allocate(gather_ispec_selected_rec(nrec,NPROC),stat=ier)
  if (ier /= 0) call exit_MPI(myrank,'Error allocating gather array')

#ifdef USE_MPI
  ! gathers infos onto master process
  call MPI_GATHER(final_distance(1),nrec,MPI_DOUBLE_PRECISION,&
        gather_final_distance(1,1),nrec,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_GATHER(xi_receiver(1),nrec,MPI_DOUBLE_PRECISION,&
        gather_xi_receiver(1,1),nrec,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_GATHER(gamma_receiver(1),nrec,MPI_DOUBLE_PRECISION,&
        gather_gamma_receiver(1,1),nrec,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_GATHER(ispec_selected_rec(1),nrec,MPI_INTEGER,&
        gather_ispec_selected_rec(1,1),nrec,MPI_INTEGER,0,MPI_COMM_WORLD,ier)

  if (myrank == 0) then
    ! selects best slice which minimum distance to receiver location
    do irec = 1, nrec
      which_proc_receiver(irec:irec) = minloc(gather_final_distance(irec,:)) - 1
    enddo
  endif
  call MPI_BCAST(which_proc_receiver(1),nrec,MPI_INTEGER,0,MPI_COMM_WORLD,ier)

#else
  ! serial
  gather_final_distance(:,1) = final_distance(:)

  gather_xi_receiver(:,1) = xi_receiver(:)
  gather_gamma_receiver(:,1) = gamma_receiver(:)
  gather_ispec_selected_rec(:,1) = ispec_selected_rec(:)
  which_proc_receiver(:) = 0
#endif

  if (USE_TRICK_FOR_BETTER_PRESSURE) then
    do irec= 1,nrec
      if (which_proc_receiver(irec) == myrank) then
        if (.not. ispec_is_acoustic(ispec_selected_rec(irec))) then
          call exit_MPI(myrank,'USE_TRICK_FOR_BETTER_PRESSURE : receivers must be in acoustic elements')
        endif
      endif
    enddo
  endif
  
  ! counts local receivers in this slice
  nrecloc = 0
  do irec = 1, nrec
    if (which_proc_receiver(irec) == myrank) then
      nrecloc = nrecloc + 1
      recloc(nrecloc) = irec
    endif
  enddo
  
  ! user output
  if (myrank == 0) then
    do irec = 1, nrec
      write(IMAIN,*)
      write(IMAIN,*) 'Station # ',irec,'    ',network_name(irec),station_name(irec)

      if (gather_final_distance(irec,which_proc_receiver(irec)+1) == HUGEVAL) &
        call exit_MPI(myrank,'Error locating receiver')

      write(IMAIN,*) '            original x: ',sngl(st_xval(irec))
      write(IMAIN,*) '            original z: ',sngl(st_zval(irec))
      write(IMAIN,*) '  distance from source: ',sngl(distance_receiver(irec))
      write(IMAIN,*) 'closest estimate found: ',sngl(gather_final_distance(irec,which_proc_receiver(irec)+1)), &
                    ' m away'
      write(IMAIN,*) ' in element ',gather_ispec_selected_rec(irec,which_proc_receiver(irec)+1)
      
      
      ! ISSUE: if stations are demanded too close to the fluid-solid interface in DG simulations, this routine  will sometimes WRONGLY detect in which material the station actually is.
      ! Why you ask? It is because it tries to find the closest INTERIOR point to the station. If it happens that the first closest interior point if the fluid one instead of the solid one, the software believes the station is in a fluid element, when it is not.
      ! In the case of 5 GLL points, the key value is the altitude midpoint between 0.1727*DZ_{fluid} and 0.1727*DZ_{solid}, that is z_k = [z_interface + 0.5*0.1727*(DZ_{fluid}-DZ_{solid})]. If z_{station} > z_k, station will be believed to be in material above. If z_{station} < z_k, station will be believed to be in material below. This is the case wherever the interface is.
      ! This poses very few problems in classic SPECFEM since quantites recorded are the same in all materials. In SPECFEM-DG, this can lead to nonsense.
      ! Two solutions. (A) Always make sure that:
      !                    1) Stations in fluid are such that z > 0.5*0.1727*(DZ_{fluid}-DZ_{solid}), and
      !                    2) Stations in solid are such that z < 0.5*0.1727*(DZ_{fluid}-DZ_{solid}),
      !                    where DZ_{fluid} is the vertical dimension of the first layer of fluid elements, and DZ_{solid} is the vertical dimension of the first layer of solid elements.
      !                (B) Implement an actual test checking if the station is in either element or not, instead of the trick actually implemented. This is easy since we know the elements are convex quadrilaterals. See above for a somewhat complete pseudo-code.
      
      ! DEBUGGING REMAINS, maybe useful in the future.
      !write(IMAIN,*) '   > which is [elastic?', ispec_is_elastic(ispec_selected_rec(irec)),&! DEBUG
      !               '], [acoustic?', ispec_is_acoustic(ispec_selected_rec(irec)),&! DEBUG
      !               ' | acoustic DG?', ispec_is_acoustic_DG(ispec_selected_rec(irec)),']'! DEBUG
      !write(IMAIN,*) '   > which is [elastic?', ispec_is_elastic(gather_ispec_selected_rec(irec,which_proc_receiver(irec)+1)),&! DEBUG
      !               '], [acoustic?', ispec_is_acoustic(gather_ispec_selected_rec(irec,which_proc_receiver(irec)+1)),&! DEBUG
      !               ' | acoustic DG?', ispec_is_acoustic_DG(gather_ispec_selected_rec(irec,which_proc_receiver(irec)+1)),']'! DEBUG
      ! gather_ispec_selected_rec(irec,which_proc_receiver(irec)+1) and ispec_selected_rec(irec) seem to contain the same information.! DEBUG
      !write(IMAIN,*) ' xminmaxzminmax of that element:', &! DEBUG
      !               minval(coord(1, pack(ibool(:,:,ispec_selected_rec(irec)), .true.))), &! DEBUG
      !               maxval(coord(1, pack(ibool(:,:,ispec_selected_rec(irec)), .true.))), &! DEBUG
      !               minval(coord(2, pack(ibool(:,:,ispec_selected_rec(irec)), .true.))), &! DEBUG
      !               maxval(coord(2, pack(ibool(:,:,ispec_selected_rec(irec)), .true.)))! DEBUG
      
      
      write(IMAIN,*) ' at process ', which_proc_receiver(irec)
      write(IMAIN,*) ' at xi,gamma coordinates = ',gather_xi_receiver(irec,which_proc_receiver(irec)+1),&
                                  gather_gamma_receiver(irec,which_proc_receiver(irec)+1)
      write(IMAIN,*)
    enddo

    write(IMAIN,*)
    write(IMAIN,*) 'end of receiver detection'
    write(IMAIN,*)
    call flush_IMAIN()

  endif
  
  ! RE-OUTPUT TO USER TO MAKE SURE THE MATERIAL OF EACH STATION IS RIGHT (SAFEGUARD)
  !if(myrank == 0) then
  call synchronize_all() ! This is a trick to force those prints to be put in the same place, after previous prints.
  do irec = 1, nrec
    if(myrank == which_proc_receiver(irec)) then
      ! The CPU in which the station is (or, to be more exact, in which the station is believed to be in) has the last word.
      ! Thus, believe it.
      write(IMAIN,*) 'Station #',irec,": closest element match is [on CPU",myrank,"] & ", & ! DEBUG
                     '[el.?', ispec_is_elastic(ispec_selected_rec(irec)),'] ', & ! DEBUG
                     '[ac.?', ispec_is_acoustic(ispec_selected_rec(irec)),'] ', & ! DEBUG
                     '[ac. DG?', ispec_is_acoustic_DG(ispec_selected_rec(irec)),']' ! DEBUG
    endif
  enddo
  call flush_IMAIN()

  ! deallocate arrays
  deallocate(final_distance)
  deallocate(gather_ispec_selected_rec)

  end subroutine locate_receivers


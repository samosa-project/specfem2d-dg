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

  subroutine compute_forces_acoustic_DG_main()

  use specfem_par

  implicit none

  ! local parameters
  ! for rk44
  !double precision :: weight_rk
  
  real(kind=CUSTOM_REAL) :: timelocal
  
  ! Parameters
  real(kind=CUSTOM_REAL), parameter :: ZEROl = 0._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: ONEl  = 1._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: TWOl  = 2._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: SIXl  = 6._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: HALFl = 0.5_CUSTOM_REAL
  
  real(kind=CUSTOM_REAL), dimension(5) :: rk4a, rk4b, rk4c
  double precision, dimension(5) :: rk4a_d, rk4b_d, rk4c_d
  
  !integer :: i, j, ispec, numelem
  
  real(kind=CUSTOM_REAL), parameter :: threshold = 0.0000001_CUSTOM_REAL
  !integer :: i ,j, ispec
  ! checks if anything to do in this slice
  if (.not. any_acoustic) return
  
  ! Intialization
  if(it == 1 .AND. i_stage == 1) then
  
    ! Sloper limiter initialization
    call setUpVandermonde()
  
    resu_rhovx = ZEROl
    resu_rhovz = ZEROl
    resu_rho   = ZEROl
    resu_E     = ZEROl
  
    if(.not. trim(MODEL) == 'external') then
        deallocate(gravityext)
        allocate(gravityext(NGLLX,NGLLZ,nspec)) 
    endif
  
    !gravityext = 9.831
   
    call initial_condition_DG()
    
    p_DG_init = (gammaext_DG - ONE)*( E_DG &
        - (HALF)*rho_DG*( veloc_x_DG**2 + veloc_z_DG**2 ) )
    
    call prepare_MPI_DG()
    
    !WRITE(*,*) it,"MAXVAL START ", maxval(rhovx_DG), minval(rhovx_DG)
    !WRITE(*,*) it,"MAXVAL START", maxval(rhovz_DG), minval(rhovz_DG)
    
  endif
  
  rk4a_d(1) = 0d0
  rk4a_d(2) = -567301805773.0/1357537059087.0
  rk4a_d(3) = -2404267990393.0/2016746695238.0
  rk4a_d(4) = -3550918686646.0/2091501179385.0
  rk4a_d(5) = -1275806237668.0/842570457699.0
    
  rk4b_d(1) = 1432997174477.0/9575080441755.0 
  rk4b_d(2) = 5161836677717.0/13612068292357.0 
  rk4b_d(3) = 1720146321549.0/2090206949498.0 
  rk4b_d(4) = 3134564353537.0/4481467310338.0 
  rk4b_d(5) = 2277821191437.0/14882151754819.0
    
  rk4c_d(1) = 0d0
  rk4c_d(2) = 1432997174477.0/9575080441755.0 
  rk4c_d(3) = 2526269341429.0/6820363962896.0 
  rk4c_d(4) = 2006345519317.0/3224310063776.0 
  rk4c_d(5) = 2802321613138.0/2924317926251.0
  
  ! RK3
  !rk4a_d(1) = 0d0
  !rk4a_d(2) = -567301805773.0/1357537059087.0
  !rk4a_d(3) = -2404267990393.0/2016746695238.0
  !rk4a_d(4) = -3550918686646.0/2091501179385.0
  !rk4a_d(5) = -1275806237668.0/842570457699.0
    
  !rk4b_d(1) = 1432997174477.0/9575080441755.0 
  !rk4b_d(2) = 5161836677717.0/13612068292357.0 
  !rk4b_d(3) = 1720146321549.0/2090206949498.0 
  !rk4b_d(4) = 3134564353537.0/4481467310338.0 
  !rk4b_d(5) = 2277821191437.0/14882151754819.0
    
  !rk4c_d(1) = 0d0
  !rk4c_d(2) = 1d0
  !rk4c_d(3) = 1d0/2d0
  !rk4c_d(4) = 0d0
  !rk4c_d(5) = 0d0
  
  rk4a = real(rk4a_d, kind=CUSTOM_REAL)
  rk4b = real(rk4b_d, kind=CUSTOM_REAL)
  rk4c = real(rk4c_d, kind=CUSTOM_REAL)
  
  timelocal = (it-1)*deltat + rk4c(i_stage)*deltat
  
  if(myrank == 0) &
  WRITE(*,*) "iter", it, i_stage, timelocal
 ! WRITE(*,*) ">> rho_DG", minval(rho_DG), maxval(rho_DG)
 ! WRITE(*,*) ">> rhovx ", maxval(rhovx_DG), minval(rhovx_DG)
  !if(i_stage == 1) &
    ! assembling potential_dot_dot or b_potential_dot_dot for acoustic elements
#ifdef USE_MPI
  if (NPROC > 1 .and. ninterface_acoustic > 0) then
    !WRITE(*,*) myrank, "rho_DG ***********************"
    call assemble_MPI_vector_DG(rho_DG, buffer_DG_rho_P)
    !WRITE(*,*) myrank, "rhovx_DG ***********************"
    call assemble_MPI_vector_DG(rhovx_DG, buffer_DG_rhovx_P)
    !WRITE(*,*) myrank, "rhovz_DG ***********************"
    call assemble_MPI_vector_DG(rhovz_DG, buffer_DG_rhovz_P)
    !WRITE(*,*) myrank, "E_DG ***********************"
    call assemble_MPI_vector_DG(E_DG, buffer_DG_E_P)
  endif
#endif
  
  call compute_forces_acoustic_DG(rho_DG, rhovx_DG, rhovz_DG, E_DG, &
        dot_rho, dot_rhovx, dot_rhovz, dot_E, &
        timelocal)
  
  if (time_stepping_scheme == 3) then
    
    ! Inverse mass matrix
    dot_rho(:)   = dot_rho(:) * rmass_inverse_acoustic_DG(:)
    dot_rhovx(:) = dot_rhovx(:) * rmass_inverse_acoustic_DG(:)
    dot_rhovz(:) = dot_rhovz(:) * rmass_inverse_acoustic_DG(:)
    dot_E(:)     = dot_E(:) * rmass_inverse_acoustic_DG(:)
  
    !WRITE(*,*) ">> 2.5rho_DG", minval(dot_rho), maxval(dot_rho)
    !WRITE(*,*) ">> 2.5rhovx ", maxval(dot_rhovx), minval(dot_rhovx)
    !WRITE(*,*) ">> 2.5rhovz ", maxval(dot_rhovz), minval(dot_rhovz)
    !WRITE(*,*) ">> 2.5E ", maxval(dot_E), minval(dot_E)
  
    ! RK5-low dissipation Update
    resu_rho = rk4a(i_stage)*resu_rho + deltat*dot_rho
    resu_rhovx = rk4a(i_stage)*resu_rhovx + deltat*dot_rhovx
    resu_rhovz = rk4a(i_stage)*resu_rhovz + deltat*dot_rhovz
    resu_E = rk4a(i_stage)*resu_E + deltat*dot_E
    
    !WRITE(*,*) ">> 2.5rho_DG", minval(resu_rho), maxval(resu_rho)
    !WRITE(*,*) ">> 2.5rhovx ", maxval(resu_rhovx), minval(resu_rhovx)
    !WRITE(*,*) ">> 2.5rhovz ", maxval(resu_rhovz), minval(resu_rhovz)
    !WRITE(*,*) ">> 2.5E ", maxval(resu_E), minval(resu_E)
    
    ! SSP-RK
    !resu_rho   = rho_DG + deltat*dot_rho
    !resu_rhovx = rhovx_DG + deltat*dot_rhovx
    !resu_rhovz = rhovz_DG + deltat*dot_rhovz
    !resu_E     = E_DG + deltat*dot_E
    
   ! do ispec = 1,nspec
!
!      ! first double loop over GLL points to compute and store gradients
!      do j = 1,NGLLZ
!        do i = 1,NGLLX
!                if(abs(resu_rho(ibool_DG(i,j,ispec))) < threshold) resu_rho(ibool_DG(i,j,ispec)) = ZEROl
!                if(abs(resu_rhovx(ibool_DG(i,j,ispec))) < threshold) resu_rhovx(ibool_DG(i,j,ispec)) = ZEROl
!                if(abs(resu_rhovz(ibool_DG(i,j,ispec))) < threshold) resu_rhovz(ibool_DG(i,j,ispec)) = ZEROl
!                if(abs(resu_E(ibool_DG(i,j,ispec))) < threshold) resu_E(ibool_DG(i,j,ispec)) = ZEROl
!        enddo
!      enddo
!      
!    enddo
    
    
    rho_DG   = rho_DG + rk4b(i_stage)*resu_rho
    rhovx_DG   = rhovx_DG + rk4b(i_stage)*resu_rhovx
    rhovz_DG   = rhovz_DG + rk4b(i_stage)*resu_rhovz
    E_DG   = E_DG + rk4b(i_stage)*resu_E  
    
    !if(i_stage == 2) then
    !rho_DG   = 0.5*(rho_DG + resu_rho + deltat*dot_rho)
    !rhovx_DG   = 0.5*(rhovx_DG + resu_rhovx + deltat*dot_rhovx)
    !rhovz_DG   = 0.5*(rhovz_DG + resu_rhovz + deltat*dot_rhovz)
    !E_DG   =  0.5*(E_DG + resu_E + deltat*dot_E)
    !endif
    
  endif
  
  !stop
  ! FLux limiter for shocks on constitutive variables
  !veloc_x_DG = rhovx_DG/rho_DG
  !veloc_z_DG = rhovz_DG/rho_DG
  !call SlopeLimit1(veloc_x_DG)
  !call SlopeLimit1(veloc_z_DG)
  
  call SlopeLimit1(rhovx_DG)
  call SlopeLimit1(rhovz_DG)
  call SlopeLimit1(rho_DG)
  
  !rhovx_DG = veloc_x_DG*rho_DG
  !rhovz_DG = veloc_z_DG*rho_DG
  
  call SlopeLimit1(E_DG)

  end subroutine compute_forces_acoustic_DG_main

  subroutine prepare_MPI_DG()
    
    use constants,only: CUSTOM_REAL,NGLLX,NGLLZ
    
    use specfem_par

    use mpi

    implicit none 
    
    include "precision.h"
    
    integer :: iinterface, ipoin, num_interface, iglob, cpt, &
        i, j, ispec, iglob_DG, nb_values, ier
    integer, dimension(nglob_DG, 3) :: link_ij_iglob
    integer, dimension(nglob_DG) :: MPI_iglob
    
    integer, dimension(max_ibool_interfaces_size_ac,ninterface_acoustic) :: &
        buffer_recv_faces_vector_DG_i, &
        buffer_send_faces_vector_DG_i, &
        buffer_recv_faces_vector_DG_j, &
        buffer_send_faces_vector_DG_j
    
    integer :: i_2, j_2
    !character(len=100) file_name
    
    !write(file_name,"('MPI_',i3.3,'.txt')") myrank
    ! Open output forcing file
    !open(100,file=file_name,form='formatted')
  
    allocate(MPI_transfer(nglob_DG, 2, 4))
    allocate(diag_MPI(max_ibool_interfaces_size_ac,ninterface_acoustic))
    
    do ispec = 1, nspec
    
    do i = 1, NGLLX
    do j = 1, NGLLZ
    
        iglob_DG = ibool_DG(i,j,ispec)
    
        link_ij_iglob(iglob_DG,1) = i
        link_ij_iglob(iglob_DG,2) = j
        link_ij_iglob(iglob_DG,3) = ispec
    
    enddo
    enddo
    
    enddo
    
    buffer_recv_faces_vector_DG_i = -1
    buffer_send_faces_vector_DG_i = -1  
    buffer_recv_faces_vector_DG_j = -1
    buffer_send_faces_vector_DG_j = -1   
    
    ! MPI SEND INFO ABOUT DIAG ELEMENT OR NOT
    do iinterface = 1, ninterface_acoustic_DG

    ! gets interface index in the range of all interfaces [1,ninterface]
    num_interface = inum_interfaces_acoustic_DG(iinterface)

    ! loops over all interface points
    do ipoin = 1, nibool_interfaces_acoustic_DG(num_interface)
    
        iglob = ibool_interfaces_acoustic_DG(ipoin,num_interface)
        
        i = link_ij_iglob(iglob,1)
        j = link_ij_iglob(iglob,2)
        
        buffer_send_faces_vector_DG_i(ipoin,num_interface) = i
        buffer_send_faces_vector_DG_j(ipoin,num_interface) = j
        
    enddo
    
    enddo
    
    do iinterface = 1, ninterface_acoustic_DG

    ! gets interface index in the range of all interfaces [1,ninterface]
    num_interface = inum_interfaces_acoustic_DG(iinterface)
        
    nb_values = nibool_interfaces_acoustic_DG(num_interface)
    
    call MPI_ISEND( buffer_send_faces_vector_DG_i(1,iinterface), &
             nb_values, MPI_INTEGER, &
             my_neighbours(num_interface), 14, MPI_COMM_WORLD, &
             tab_requests_send_recv_DG(iinterface), ier)
        
    if (ier /= MPI_SUCCESS) then
      call exit_MPI(myrank,'MPI_ISEND unsuccessful in assemble_MPI_vector_start')
    endif

    ! starts a non-blocking receive
    call MPI_IRECV ( buffer_recv_faces_vector_DG_i(1,iinterface), &
             nb_values, MPI_INTEGER, &
             my_neighbours(num_interface), 14, MPI_COMM_WORLD, &
             tab_requests_send_recv_DG(ninterface_acoustic_DG+iinterface), ier)
     
    if (ier /= MPI_SUCCESS) then
      call exit_MPI(myrank,'MPI_IRECV unsuccessful in assemble_MPI_vector')
    endif
    
    call MPI_ISEND( buffer_send_faces_vector_DG_j(1,iinterface), &
             nb_values, MPI_INTEGER, &
             my_neighbours(num_interface), 14, MPI_COMM_WORLD, &
             tab_requests_send_recv_DG(ninterface_acoustic_DG*2+iinterface), ier)
        
    if (ier /= MPI_SUCCESS) then
      call exit_MPI(myrank,'MPI_ISEND unsuccessful in assemble_MPI_vector_start')
    endif

    ! starts a non-blocking receive
    call MPI_IRECV ( buffer_recv_faces_vector_DG_j(1,iinterface), &
             nb_values, MPI_INTEGER, &
             my_neighbours(num_interface), 14, MPI_COMM_WORLD, &
             tab_requests_send_recv_DG(ninterface_acoustic_DG*3+iinterface), ier)
     
    if (ier /= MPI_SUCCESS) then
      call exit_MPI(myrank,'MPI_IRECV unsuccessful in assemble_MPI_vector')
    endif
    
    enddo
    
    ! waits for MPI requests to complete (recv)
    ! each wait returns once the specified MPI request completed
    do iinterface = 1, ninterface_acoustic_DG*4
      call MPI_Wait(tab_requests_send_recv_DG(iinterface), &
                  MPI_STATUS_IGNORE, ier)
    enddo
    
    MPI_transfer = -1
    MPI_iglob    = 0
    
    do iinterface = 1, ninterface_acoustic_DG

    ! gets interface index in the range of all interfaces [1,ninterface]
    num_interface = inum_interfaces_acoustic_DG(iinterface)

    ! loops over all interface points
    do ipoin = 1, nibool_interfaces_acoustic_DG(num_interface)
    
        iglob = ibool_interfaces_acoustic_DG(ipoin,num_interface)
        
        i = link_ij_iglob(iglob,1)
        j = link_ij_iglob(iglob,2)
        ispec = link_ij_iglob(iglob,3)
        
        i_2     = buffer_recv_faces_vector_DG_i(ipoin, iinterface)
        j_2     = buffer_recv_faces_vector_DG_j(ipoin, iinterface)
        
        !WRITE(100,*) myrank,i,j,ispec,coord(:,ibool(i,j,ispec))
        
        if( (neighbor_DG(i,j,ispec,3) == -1 .OR. &
                neighbor_DG_corner(i,j,ispec,3) == -1) .AND. (i == i_2 .OR. j == j_2) ) then
        
        !WRITE(*,*) myrank,"------------>", i,j,ispec,i_2,j_2,coord(:,ibool(i,j,ispec))
        !WRITE(*,*) "-----", neighbor_DG(i,j,ispec,3), neighbor_DG_corner(i,j,ispec,3)
        !WRITE(*,*) "***********"
        
        MPI_iglob(iglob) = MPI_iglob(iglob) + 1
        cpt = MPI_iglob(iglob)
        MPI_transfer(iglob,cpt,1) = ipoin 
        MPI_transfer(iglob,cpt,2) = num_interface 
        
        MPI_transfer(iglob,cpt,3) = i_2
        MPI_transfer(iglob,cpt,4) = j_2
       
       endif
        
    enddo
    
    enddo
    
    !close(100)
    
    !stop
    
  end subroutine prepare_MPI_DG

  ! Setup Vandermonde matrices and derivation matrix for slope limiter
  ! purposes
  subroutine setUpVandermonde()
    
    use constants,only: CUSTOM_REAL,NGLLX,NGLLZ
    
    use specfem_par,only: Vandermonde, invVandermonde, Drx, Drz

    implicit none 
    
    call compute_Vander_matrices(NGLLX*NGLLZ,0._CUSTOM_REAL,0._CUSTOM_REAL,&
       Vandermonde, invVandermonde, Drx, Drz )
       
       !WRITE(*,*) "maxval Vandermonde", maxval(Vandermonde), minval(Vandermonde)
       
  end subroutine setUpVandermonde
  
  subroutine compute_Vander_matrices(np,alpha,beta,&
       V1D, V1D_inv, Drx, Drz )

    use constants,only: CUSTOM_REAL,NGLLX,NGLLZ,GAUSSALPHA,GAUSSBETA
    use specfem_par,only: xigll, zigll!, hprime_xx, hprime_zz

    implicit none 

  ! Parameters
  real(kind=CUSTOM_REAL), parameter :: ZEROl = 0._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: ONEl  = 1._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: TWOl  = 2._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: SIXl  = 6._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: HALFl = 0.5_CUSTOM_REAL

  !double precision, parameter :: zero=0.d0,one=1.d0,two=2.d0

  integer np
  real(kind=CUSTOM_REAL) alpha,beta
  real(kind=CUSTOM_REAL), parameter :: threshold = 0.00001_CUSTOM_REAL
  real(kind=CUSTOM_REAL) V1D(np,np), V1D_inv(np, np),&
         Drx(np, np), Drz(np, np)!,  diagMass(np,np)!, V1D_der_x(np, np),&
          !V1D_der_z(np, np), V1D_t(np, np)!, diagMass(np,np)!, R(np, np)
  double precision, dimension(np,np) :: V1D_d, V1D_inv_d
  !double precision, dimension(NGLLX) :: hx_d, hprimex_d
  !double precision, dimension(NGLLZ) :: hz_d, hprimez_d   
  !real(kind=CUSTOM_REAL), dimension(NGLLX) :: hx, hprimex
  !real(kind=CUSTOM_REAL), dimension(NGLLZ) :: hz, hprimez  

  double precision :: p,pd,pm1,pdm1,pm2,pdm2,p2,pd2

  integer i,j,k,errorflag,l,m,n
  
  !double precision, external :: hgll

  !!!!!!!!!!!!!!!!!!!!!!!
  !!!! => NEEDS TO BE FINISHED
  !!!!!!!!!!!!!!!!!!!!!!!
  ! Init
  V1D     = ZEROl
  !V1D_t   = ZEROl
  V1D_inv = ZEROl
  Drx     = ZEROl
  Drz     = ZEROl
  V1D_d   = 0d0
  
  ! NEW MATRIcES
 ! do k=1,np
 k = 0
  do m = 1,NGLLX
    do n = 1,NGLLZ
    k = k+1
  l = 0
  do i = 1,NGLLX
    do j = 1,NGLLZ
      l = l + 1
      !WRITE(*,*) "1*********"
      !call lagrange_any(xigll(j),NGLLX,xigll,hx_d,hprimex_d)
      !hx_d(i) = hgll(i-1,xigll(j-1),xigll,NGLLX) !lagrange_deriv_GLL(i1-1,i2-1,xigll,NGLLX)
      call jacobf(p,pd,pm1,pdm1,pm2,pdm2,i-1,GAUSSALPHA,GAUSSBETA,xigll(m))
      call jacobf(p2,pd2,pm1,pdm1,pm2,pdm2,j-1,GAUSSALPHA,GAUSSBETA,zigll(n))
      !call jacobf(p,pd,pm1,pdm1,pm2,pdm2,m-1,GAUSSALPHA,GAUSSBETA,xigll(i))
      !call jacobf(p2,pd2,pm1,pdm1,pm2,pdm2,n-1,GAUSSALPHA,GAUSSBETA,zigll(j))
  
      !WRITE(*,*) "TEST LAGRANGE h", h(i2)
      !WRITE(*,*) "TEST LAGRANGE hprime", hprimex_d(i), xigll(j),NGLLX
      !WRITE(*,*) "TEST LAGRANGE hprime_xx(i2,i1)", hprime_xx(i2,i1)
      !WRITE(*,*) "*********"
      !hprimewgll_xx(i2,i1) = wxgll(i2) * hprime_xx(i2,i1)
      !Drx(k,j) = REAL(hprimex_d(i), kind=CUSTOM_REAL)
      !V1D(k,j) = REAL(hx_d(i), kind=CUSTOM_REAL)
      !V1D_der_x(k,l) = REAL(pd, kind=CUSTOM_REAL)
      !V1D_der_z(k,l) = REAL(pd2, kind=CUSTOM_REAL)
      V1D(k,l) = REAL(p*p2, kind=CUSTOM_REAL)
      if(abs(V1D(k,l)) < threshold) V1D(k,l) = ZEROl
      !WRITE(*,*) "V1D(",k,",",l,") = ", V1D(k,l)
      V1D_d(k,l) = p*p2
      !V1D_t(l,k) = V1D(k,l)
    enddo
  enddo
  !enddo
    enddo
  enddo
  
  call FINDInv(V1D_d, V1D_inv_d, np, errorflag)
  do j=1,np
    do i=1,np
        V1D_inv(i,j) = REAL(V1D_inv_d(i,j), kind=CUSTOM_REAL)
        if(abs(V1D_inv(i,j)) < threshold) V1D_inv(i,j) = ZEROl
    enddo
  enddo
  
  !do j=1,np
  !  do i=1,np
  !      do k=1,np
  !      Drx(i,j) = Drx(i,j) + V1D_der_x(i,k)*V1D_inv(k,j)
  !      !WRITE(*,*) "A(",k,",",l,") = ", V1D_der_x(i,k), REAL(V1D_inv(k,j), kind=CUSTOM_REAL)
  !      Drz(i,j) = Drz(i,j) + V1D_der_z(i,k)*V1D_inv(k,j)
  !      enddo
  !     ! WRITE(*,*) "Dr", Dr(i,j) 
  !  enddo
  !enddo
  
  !TEST INVERSION
  !R = ZEROl
  !do k=1,np
  !      do l=1,np
  !              do m=1,np
  !              R(k,l)  = R(k,l) + V1D(k,m)*V1D_inv(m,l)
  !              enddo
  !              WRITE(*,*) "R(k,l)",k,l,R(k,l)
  !    enddo
  !    
  !enddo

  end subroutine compute_Vander_matrices
  
  !Subroutine to find the inverse of a square matrix
!Author : Louisda16th a.k.a Ashwith J. Rego
!Reference : Algorithm has been well explained in:
!http://math.uww.edu/~mcfarlat/inverse.htm 
!http://www.tutor.ms.unimelb.edu.au/matrix/matrix_inverse.html

SUBROUTINE FINDInv(matrix, inverse, n, errorflag)

use constants,only: CUSTOM_REAL

IMPLICIT NONE
!Declarations
INTEGER, INTENT(IN) :: n
INTEGER, INTENT(OUT) :: errorflag !Return error status. -1 for error, 0 for normal
double precision, INTENT(IN), DIMENSION(n,n) :: matrix !Input matrix
double precision, INTENT(OUT), DIMENSION(n,n) :: inverse !Inverted matrix

LOGICAL :: FLAG = .TRUE.
INTEGER :: i, j, k
double precision :: m
double precision, DIMENSION(n,2*n) :: augmatrix !augmented matrix

!Augment input matrix with an identity matrix
DO i = 1, n
DO j = 1, 2*n
IF (j <= n ) THEN
augmatrix(i,j) = matrix(i,j)
ELSE IF ((i+n) == j) THEN
augmatrix(i,j) = 1
Else
augmatrix(i,j) = 0
ENDIF
END DO
END DO

!Reduce augmented matrix to upper traingular form
DO k =1, n-1
IF (augmatrix(k,k) == 0) THEN
        FLAG = .FALSE.
        DO i = k+1, n
                        IF (augmatrix(i,k) /= 0) THEN
                        DO j = 1,2*n
                                augmatrix(k,j) = augmatrix(k,j)+augmatrix(i,j)
                        END DO
                        FLAG = .TRUE.
                        EXIT
                ENDIF
                IF (FLAG .EQV. .FALSE.) THEN
                        PRINT*, "Matrix is non - invertible"
                        inverse = 0
                        errorflag = -1
                        return
                ENDIF
        END DO
ENDIF

DO j = k+1, n 
        m = augmatrix(j,k)/augmatrix(k,k)
        DO i = k, 2*n
                augmatrix(j,i) = augmatrix(j,i) - m*augmatrix(k,i)
        END DO
END DO

END DO

!Test for invertibility
DO i = 1, n
IF (augmatrix(i,i) == 0) THEN
PRINT*, "Matrix is non - invertible"
inverse = 0
errorflag = -1
return
ENDIF
END DO

!Make diagonal elements as 1
DO i = 1 , n
m = augmatrix(i,i)
DO j = i , (2 * n) 
augmatrix(i,j) = (augmatrix(i,j) / m)
END DO
END DO

!Reduced right side half of augmented matrix to identity matrix
DO k = n-1, 1, -1
DO i =1, k
m = augmatrix(i,k+1)
DO j = k, (2*n)
augmatrix(i,j) = augmatrix(i,j) -augmatrix(k+1,j) * m
END DO
END DO
END DO 

!store answer
DO i =1, n
DO j = 1, n
inverse(i,j) = augmatrix(i,j+n)
END DO
END DO
errorflag = 0
RETURN
END

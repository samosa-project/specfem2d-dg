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
  
  ! SLOPE LIMITER
  
  subroutine SlopeLimit1_old(Q, timelocal, type_var)
  
    use constants,only: CUSTOM_REAL,NGLLX,NGLLZ

    use specfem_par, only: nspec, nglob_DG, ibool_DG, Vandermonde, invVandermonde, &
        neighbor_DG, neighbor_DG_corner,is_corner, NGLLX, NGLLZ, &
        max_interface_size, ninterface, NPROC, MPI_transfer, ninterface_acoustic!, myrank!, ibool, ibool_before_perio!, coord, i_stage&
        !ispec_is_acoustic_surface, ispec_is_acoustic_surface_corner!,ibool_before_perio,

    implicit none 
  
    integer type_var
    real(kind=CUSTOM_REAL) timelocal
    real(kind=CUSTOM_REAL) :: rho_DG_P, rhovx_DG_P, rhovz_DG_P, E_DG_P, &
                veloc_x_DG_P, veloc_z_DG_P, p_DG_P, e1_DG_P
    
    integer :: iglob, iglob2, ispec, j, k, l, h, prod_ind
    integer, parameter :: NGLL = NGLLX*NGLLZ
    real(kind=CUSTOM_REAL) :: vkm1_x,vkp1_x,vkm1_z,vkp1_z!, vkm1_x_save
    real(kind=CUSTOM_REAL), dimension(nglob_DG) :: Q, ulimit
    real(kind=CUSTOM_REAL), dimension(NGLL) :: uh, uavg, uhl, ulimit_temp
    real(kind=CUSTOM_REAL), dimension(1:nspec,1:NGLL) :: ul
    real(kind=CUSTOM_REAL), dimension(1:nspec) :: v!, vk
    real(kind=CUSTOM_REAL), dimension(nglob_DG) :: v_MPI
    real(kind=CUSTOM_REAL), dimension(NGLLX*max_interface_size,ninterface) :: buffer_v
    
    integer :: ipoin,num_interface
    !logical, dimension(NGLLX, NGLLZ) :: ispec_bound
    logical ispec_bound_MPI, has_been_activated
    !real(kind=CUSTOM_REAL), dimension(1:NGLLX,1:NGLLZ) :: val_in_elmnt
    
    logical :: activate
    
    ! Parameters
    real(kind=CUSTOM_REAL), parameter :: ZEROl = 0._CUSTOM_REAL
    
    has_been_activated = .false.
    
    ulimit = ZEROl
    buffer_v = ZEROl

    do ispec=1,nspec
    
    uh = ZEROl
    ! Compute cell averages
    ! Compute => uh = invV*u
    do k=1,NGLL
        prod_ind = 0
        do j=1,NGLLX
          do h=1,NGLLZ
                iglob2 = ibool_DG(j,h,ispec)
                prod_ind = prod_ind + 1
                uh(k)  = uh(k) + &
                        invVandermonde(k,prod_ind)*Q(iglob2)
               ! val_in_elmnt(j,h) = Q(iglob2)
          enddo
      enddo
    enddo
    
    ! Extract linear polynomial
    uhl = uh
    uhl(3:NGLLX) = ZEROl
    uhl(NGLLX+2:NGLL) = ZEROl
    
    ul(ispec,:) = ZEROl
    ! => ul = V*uhl;
    ! Compute => uavg = V*uh
    do k=1,NGLL
        do j=1,NGLL
                ul(ispec,k)  = ul(ispec,k) + Vandermonde(k,j)*uhl(j)
        enddo
    enddo
    
    ! Remove high order coef (> 0)
    uh(2:NGLL) = ZEROl
    
    ! Init
    uavg = ZEROl
    ! Compute => uavg = V*uh
    do k=1,NGLL
        do j=1,NGLL
                uavg(k)  = uavg(k) + Vandermonde(k,j)*uh(j)
        enddo
    enddo
    ! Store cell average
   v(ispec) = uavg(1)
   
   ! Change Quentin 14 mars 17:16
   do k=1,NGLLX
          do l=1,NGLLZ
          v_MPI(ibool_DG(k,l,ispec)) = v(ispec)
          enddo
   enddo
   
   enddo ! do ispec=1,NSPEC
   
   ! find cell averages
   !vk = v
   
   if (NPROC > 1 .and. ninterface_acoustic > 0) then
     ! Change Quentin 14 mars 17:16
     call assemble_MPI_vector_DG(v_MPI, buffer_v)
   endif
   
   !nb_mod = 0
    do ispec=1,nspec
    
     ispec_bound_MPI = .false.
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !
     !  LEFT NEIGHBOR
     !
     ! MPI neighbor
     ipoin         = -1
     num_interface = -1
     if(NPROC > 1) then
        iglob         = ibool_DG(1,NGLLZ/2,ispec)
        ipoin         = MPI_transfer(iglob, 1,1)
        num_interface = MPI_transfer(iglob, 1,2)
        
        if(ipoin > -1) vkm1_x = buffer_v(ipoin,num_interface)
     endif
     
     if(ipoin == -1 .AND. neighbor_DG(1,NGLLZ/2,ispec,3) == - 1) then
     ! Test with boundary conditions
     call boundary_condition_DG(1, NGLLZ/2, ispec, timelocal, rho_DG_P, rhovx_DG_P, rhovz_DG_P, E_DG_P, &
                veloc_x_DG_P, veloc_z_DG_P, p_DG_P, e1_DG_P)
    
     if(type_var == 1) vkm1_x = rho_DG_P
     if(type_var == 2) vkm1_x = rhovx_DG_P
     if(type_var == 3) vkm1_x = rhovz_DG_P
     if(type_var == 4) vkm1_x = E_DG_P
     
     ispec_bound_MPI = .true.
     endif
     
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !
     !  RIGHT NEIGHBOR
     !
     ! MPI neighbor
     ipoin         = -1
     num_interface = -1
     if(NPROC > 1) then
        iglob         = ibool_DG(NGLLX,NGLLZ/2,ispec)
        ipoin         = MPI_transfer(iglob, 1,1)
        num_interface = MPI_transfer(iglob, 1,2)
        
        if(ipoin > -1) vkp1_x = buffer_v(ipoin,num_interface)
     endif
     
     if(ipoin == -1 .AND. neighbor_DG(NGLLX,NGLLZ/2,ispec,3) == - 1) then
     call boundary_condition_DG(NGLLX, NGLLZ/2, ispec, timelocal, rho_DG_P, rhovx_DG_P, rhovz_DG_P, E_DG_P, &
                veloc_x_DG_P, veloc_z_DG_P, p_DG_P, e1_DG_P)
    
     if(type_var == 1) vkp1_x = rho_DG_P
     if(type_var == 2) vkp1_x = rhovx_DG_P
     if(type_var == 3) vkp1_x = rhovz_DG_P
     if(type_var == 4) vkp1_x = E_DG_P
     
     ispec_bound_MPI = .true.
     endif
     
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !
     !  BOTTOM NEIGHBOR
     !
     ipoin         = -1
     num_interface = -1
     if(NPROC > 1) then
        iglob         = ibool_DG(NGLLX/2, 1, ispec)
        ipoin         = MPI_transfer(iglob, 1,1)
        num_interface = MPI_transfer(iglob, 1,2)
        
        if(ipoin > -1) vkm1_z = buffer_v(ipoin,num_interface)
     endif
     
     if(ipoin == -1 .AND. neighbor_DG(NGLLX/2,1,ispec,3) == - 1) then
     call boundary_condition_DG(NGLLX/2, 1, ispec, timelocal, rho_DG_P, rhovx_DG_P, rhovz_DG_P, E_DG_P, &
                veloc_x_DG_P, veloc_z_DG_P, p_DG_P, e1_DG_P)
    
     if(type_var == 1) vkm1_z = rho_DG_P
     if(type_var == 2) vkm1_z = rhovx_DG_P
     if(type_var == 3) vkm1_z = rhovz_DG_P
     if(type_var == 4) vkm1_z = E_DG_P
     
     ispec_bound_MPI = .true.
     endif
     
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !
     !  TOP NEIGHBOR
     !
     ipoin         = -1
     num_interface = -1
     if(NPROC > 1) then
        iglob         = ibool_DG(NGLLX/2, NGLLZ, ispec)
        ipoin         = MPI_transfer(iglob, 1,1)
        num_interface = MPI_transfer(iglob, 1,2)
        
        if(ipoin > -1) then 
                vkp1_z = buffer_v(ipoin,num_interface)
        endif
     endif
     
     if(ipoin == -1 .AND. neighbor_DG(NGLLX/2,NGLLZ,ispec,3) == - 1) then
     call boundary_condition_DG(NGLLX/2, NGLLZ, ispec, timelocal, rho_DG_P, rhovx_DG_P, rhovz_DG_P, E_DG_P, &
                veloc_x_DG_P, veloc_z_DG_P, p_DG_P, e1_DG_P)
    
     if(type_var == 1) vkp1_z = rho_DG_P
     if(type_var == 2) vkp1_z = rhovx_DG_P
     if(type_var == 3) vkp1_z = rhovz_DG_P
     if(type_var == 4) vkp1_z = E_DG_P
     
     ispec_bound_MPI = .true.
     endif
     
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !
     ! NEIGHBOR IN THE MPI DOMAIN
     !
     
     if(neighbor_DG(1,NGLLZ/2,ispec,3) > - 1) vkm1_x = v(neighbor_DG(1,NGLLZ/2,ispec,3))
     
     if(neighbor_DG(NGLLX,NGLLZ/2,ispec,3) > - 1) vkp1_x = v(neighbor_DG(NGLLX,NGLLZ/2,ispec,3))
     
     if(neighbor_DG(NGLLX/2,1,ispec,3) > - 1) vkm1_z = v(neighbor_DG(NGLLX/2,1,ispec,3))
     
     if(neighbor_DG(NGLLX/2,NGLLZ,ispec,3) > - 1) vkp1_z = v(neighbor_DG(NGLLX/2,NGLLZ,ispec,3)) !v(neighbor_DG(NGLLX/2,NGLLZ,ispec,3))
     
     if(.not. ispec_bound_MPI) then
     ! Compute Limited flux
     call SlopeLimitLin(ulimit_temp, ispec, ul(ispec,:),v(ispec),&
        vkm1_x,vkp1_x,vkm1_z,vkp1_z,activate,type_var)
     endif
     
     prod_ind = 0
     do k=1,NGLLX
       do l=1,NGLLZ
        iglob = ibool_DG(k,l,ispec)
        prod_ind = prod_ind + 1
        
        if(.not. activate .OR. ispec_bound_MPI) then
        !if(.not. activate) then
                ulimit(iglob) = Q(iglob)
        else
                ulimit(iglob) = ulimit_temp(prod_ind)
                has_been_activated = .true.
        endif
        
       enddo
     enddo
     
   enddo
   
   if(has_been_activated) stop 'TOTT'!WRITE(*,*) "SLOPE LIMITER HAS BEEN ACTIVATED", timelocal, type_var
   
   ! Update solution
   Q(1:nglob_DG) = ulimit(1:nglob_DG)
   
  end subroutine SlopeLimit1_old
  
  subroutine SlopeLimit1(Q, timelocal, type_var)
  
    use constants,only: CUSTOM_REAL,NGLLX,NGLLZ

    use specfem_par, only: nspec, nglob_DG, ibool_DG, Vandermonde, invVandermonde, &
        neighbor_DG, neighbor_DG_corner,i_stage,is_corner, ibool, coord,NGLLX, NGLLZ, &
        max_interface_size, ninterface, NPROC, MPI_transfer, myrank, ibool_before_perio!, &
        !ispec_is_acoustic_surface, ispec_is_acoustic_surface_corner!,ibool_before_perio,

    implicit none 
  
    integer type_var
    real(kind=CUSTOM_REAL) timelocal
    real(kind=CUSTOM_REAL) :: rho_DG_P, rhovx_DG_P, rhovz_DG_P, E_DG_P, &
                veloc_x_DG_P, veloc_z_DG_P, p_DG_P, e1_DG_P
    
    integer :: iglob, iglob2, ispec, j, k, l, h, prod_ind
    integer, parameter :: NGLL = NGLLX*NGLLZ
    real(kind=CUSTOM_REAL) :: vkm1_x,vkp1_x,vkm1_z,vkp1_z!, vkm1_x_save
    real(kind=CUSTOM_REAL), dimension(nglob_DG) :: Q, ulimit
    real(kind=CUSTOM_REAL), dimension(NGLL) :: uh, uavg, uhl, ulimit_temp
    real(kind=CUSTOM_REAL), dimension(1:nspec,1:NGLL) :: ul
    real(kind=CUSTOM_REAL), dimension(1:nspec) :: v!, vk
    real(kind=CUSTOM_REAL), dimension(nglob_DG) :: v_MPI
    real(kind=CUSTOM_REAL), dimension(NGLLX*max_interface_size,ninterface) :: buffer_v
    
    integer, dimension(nglob_DG) :: MPI_iglob
    integer :: ipoin,num_interface
    !logical, dimension(NGLLX, NGLLZ) :: ispec_bound
    logical is_ispec_onbound
    !real(kind=CUSTOM_REAL), dimension(1:NGLLX,1:NGLLZ) :: val_in_elmnt
    
    logical :: activate, activate_total
    
    ! Parameters
    real(kind=CUSTOM_REAL), parameter :: ZEROl = 0._CUSTOM_REAL
    
    !character(len=100) file_name
    
    !write(file_name,"('isMPI_slope_',i3.3,'.txt')") myrank
  ! Open output forcing file
    !open(10,file=file_name,form='formatted')

    MPI_iglob = 1

    ulimit = ZEROl
    buffer_v = ZEROl

    do ispec=1,nspec
    
    uh = ZEROl
    ! Compute cell averages
    ! Compute => uh = invV*u
    do k=1,NGLL
        prod_ind = 0
        do j=1,NGLLX
          do h=1,NGLLZ
                iglob2 = ibool_DG(j,h,ispec)
                prod_ind = prod_ind + 1
                uh(k)  = uh(k) + &
                        invVandermonde(k,prod_ind)*Q(iglob2)
               ! val_in_elmnt(j,h) = Q(iglob2)
          enddo
      enddo
    enddo
    
    ! Extract linear polynomial
    uhl = uh
    uhl(3:NGLLX) = ZEROl
    uhl(NGLLX+2:NGLL) = ZEROl
    
    ul(ispec,:) = ZEROl
    ! => ul = V*uhl;
    ! Compute => uavg = V*uh
    do k=1,NGLL
        do j=1,NGLL
                ul(ispec,k)  = ul(ispec,k) + Vandermonde(k,j)*uhl(j)
        enddo
    enddo
    
    ! Remove high order coef (> 0)
    uh(2:NGLL) = ZEROl
    
    ! Init
    uavg = ZEROl
    ! Compute => uavg = V*uh
    do k=1,NGLL
        do j=1,NGLL
                uavg(k)  = uavg(k) + Vandermonde(k,j)*uh(j)
        enddo
    enddo
    ! Store cell average
   v(ispec) = uavg(1)
   !WRITE(*,*) "maxval/minval", maxval(Vandermonde), minval(Vandermonde), &
   !     maxval(invVandermonde), minval(invVandermonde)
   !WRITE(*,*) ispec, "cell average v(ispec)", v(ispec)!, maxval(val_in_elmnt),minval(val_in_elmnt)!, minval(Q(ibool_DG(:,1,ispec))), &
    !    WRITE(*,*) ispec, "cell average v(ispec)", v(ispec), maxval(Q(ibool_DG(:,1,ispec))), minval(Q(ibool_DG(:,1,ispec))), &
    !    maxval(Q(ibool_DG(:,5,ispec))), minval(Q(ibool_DG(:,1,ispec)))
   do k=1,NGLLX
          do l=1,NGLLZ
          v_MPI(ibool_DG(k,l,ispec)) = v(ispec)
          enddo
   enddo
   
   enddo ! do ispec=1,NSPEC
   
   ! find cell averages
   !vk = v
   
   call assemble_MPI_vector_DG(v_MPI, buffer_v)
   
   activate_total = .false.
   !nb_mod = 0
    do ispec=1,nspec
    
     ! Test with boundary conditions
     call boundary_condition_DG(1, NGLLZ/2, ispec, timelocal, rho_DG_P, rhovx_DG_P, rhovz_DG_P, E_DG_P, &
                veloc_x_DG_P, veloc_z_DG_P, p_DG_P, e1_DG_P)
    
     if(type_var == 1) vkm1_x = rho_DG_P
     if(type_var == 2) vkm1_x = rhovx_DG_P
     if(type_var == 3) vkm1_x = rhovz_DG_P
     if(type_var == 4) vkm1_x = E_DG_P
     
     call boundary_condition_DG(NGLLX, NGLLZ/2, ispec, timelocal, rho_DG_P, rhovx_DG_P, rhovz_DG_P, E_DG_P, &
                veloc_x_DG_P, veloc_z_DG_P, p_DG_P, e1_DG_P)
    
     if(type_var == 1) vkp1_x = rho_DG_P
     if(type_var == 2) vkp1_x = rhovx_DG_P
     if(type_var == 3) vkp1_x = rhovz_DG_P
     if(type_var == 4) vkp1_x = E_DG_P
     
     call boundary_condition_DG(NGLLX/2, 1, ispec, timelocal, rho_DG_P, rhovx_DG_P, rhovz_DG_P, E_DG_P, &
                veloc_x_DG_P, veloc_z_DG_P, p_DG_P, e1_DG_P)
    
     if(type_var == 1) vkm1_z = rho_DG_P
     if(type_var == 2) vkm1_z = rhovx_DG_P
     if(type_var == 3) vkm1_z = rhovz_DG_P
     if(type_var == 4) vkm1_z = E_DG_P
     
     call boundary_condition_DG(NGLLX/2, NGLLZ, ispec, timelocal, rho_DG_P, rhovx_DG_P, rhovz_DG_P, E_DG_P, &
                veloc_x_DG_P, veloc_z_DG_P, p_DG_P, e1_DG_P)
    
     if(type_var == 1) vkp1_z = rho_DG_P
     if(type_var == 2) vkp1_z = rhovx_DG_P
     if(type_var == 3) vkp1_z = rhovz_DG_P
     if(type_var == 4) vkp1_z = E_DG_P
    
     is_ispec_onbound = .true.
    
   ! apply slope limiter to selected elements
     ! Find neighbors
     !vkm1_x = 0.
     !do k=1,NGLLZ
     !!vkm1_x = vkm1_x + (1./NGLLZ)*Q(ibool_DG(1,k,ispec))!v(ispec)
     !enddo
     !vkm1_x = Q(ibool_DG(1,NGLLZ/2,ispec))
     !vkm1_x = v(ispec)
     if(neighbor_DG(1,NGLLZ/2,ispec,3) > - 1) then
        vkm1_x = v(neighbor_DG(1,NGLLZ/2,ispec,3))
        is_ispec_onbound = .false.
     endif
     
     !vkp1_x = 0.
     !do k=1,NGLLZ
     !vkp1_x = vkp1_x + (1./NGLLZ)*Q(ibool_DG(NGLLX,k,ispec))!v(ispec)
     !enddo
     !vkp1_x = Q(ibool_DG(NGLLX,NGLLZ/2,ispec))!v(ispec)
     !vkp1_x = v(ispec)
     if(neighbor_DG(NGLLX,NGLLZ/2,ispec,3) > - 1) then
        vkp1_x = v(neighbor_DG(NGLLX,NGLLZ/2,ispec,3))
        is_ispec_onbound = .false.
     endif
     
     !vkm1_z = 0.
     !do k=1,NGLLX
     !vkm1_z = vkm1_z + (1./NGLLX)*Q(ibool_DG(k,1,ispec))!v(ispec)
     !enddo
     !vkm1_z = Q(ibool_DG(NGLLX/2,1,ispec))!v(ispec)
     !vkm1_z = v(ispec)
     if(neighbor_DG(NGLLX/2,1,ispec,3) > - 1) then
        vkm1_z = v(neighbor_DG(NGLLX/2,1,ispec,3))
        is_ispec_onbound = .false.
     endif
     
     !vkp1_z = 0.
     !do k=1,NGLLX
     !vkp1_z = vkp1_z + (1./NGLLX)*Q(ibool_DG(k,NGLLZ,ispec))!v(ispec)
     !enddo
     !vkp1_z = Q(ibool_DG(NGLLX/2,NGLLZ,ispec))!v(ispec)
     !vkp1_z = v(ispec)
     if(neighbor_DG(NGLLX/2,NGLLZ,ispec,3) > - 1) then
        vkp1_z = v(neighbor_DG(NGLLX/2,NGLLZ,ispec,3))
        is_ispec_onbound = .false.
     endif
     
     ipoin         = -1
     num_interface = -1
     if(NPROC > 1) then
        iglob         = ibool_DG(1,NGLLZ/2,ispec)
        ipoin         = MPI_transfer(iglob, 1,1)
        num_interface = MPI_transfer(iglob, 1,2)
        
        if(ipoin > -1) then
                vkm1_x = buffer_v(ipoin,num_interface)
                is_ispec_onbound = .false.
        endif
     endif
     
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !
     !  RIGHT NEIGHBOR
     !
     ! MPI neighbor
     ipoin         = -1
     num_interface = -1
     if(NPROC > 1) then
        iglob         = ibool_DG(NGLLX,NGLLZ/2,ispec)
        ipoin         = MPI_transfer(iglob, 1,1)
        num_interface = MPI_transfer(iglob, 1,2)
        
        if(ipoin > -1) then
                vkp1_x = buffer_v(ipoin,num_interface)
                is_ispec_onbound = .false.
        endif
     endif
     
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !
     !  BOTTOM NEIGHBOR
     !
     ipoin         = -1
     num_interface = -1
     if(NPROC > 1) then
        iglob         = ibool_DG(NGLLX/2, 1, ispec)
        ipoin         = MPI_transfer(iglob, 1,1)
        num_interface = MPI_transfer(iglob, 1,2)
        
        if(ipoin > -1) then
                vkm1_z = buffer_v(ipoin,num_interface)
                is_ispec_onbound = .false.
        endif
     endif
     
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !
     !  TOP NEIGHBOR
     !
     ipoin         = -1
     num_interface = -1
     if(NPROC > 1) then
        iglob         = ibool_DG(NGLLX/2, NGLLZ, ispec)
        ipoin         = MPI_transfer(iglob, 1,1)
        num_interface = MPI_transfer(iglob, 1,2)
        
        if(ipoin > -1) then 
                vkp1_z = buffer_v(ipoin,num_interface)
                is_ispec_onbound = .false.
        endif
     endif
     
     ! Compute Limited flux
     activate = .false.
     if(.not. is_ispec_onbound) &
     call SlopeLimitLin(ulimit_temp, ispec, ul(ispec,:),v(ispec),&
        vkm1_x,vkp1_x,vkm1_z,vkp1_z,activate,type_var)
     
     if(activate) activate_total = activate
     
     !if(ispec_bound) &
     !WRITE(*,*) ">>>>", coord(:,ibool(1,1,ispec)), ispec_bound
     
     prod_ind = 0
     do k=1,NGLLX
       do l=1,NGLLZ
        iglob = ibool_DG(k,l,ispec)
        prod_ind = prod_ind + 1
        
        if(.not. is_ispec_onbound) &
        ulimit(iglob) = ulimit_temp(prod_ind)
        
        !if(ulimit_temp(prod_ind) /= Q(iglob)) &
        !WRITE(*,*) "DIFF", ulimit_temp(prod_ind), Q(iglob)
        
        !if(abs(ulimit_temp(prod_ind) - Q(iglob)) < 1.*10**(-3)) ulimit(iglob) = Q(iglob)
        
        ! No slope limiter on the outer boundary
        !if(ispec_bound(k,l))  ulimit(iglob) = Q(iglob)
        if(.not. activate) ulimit(iglob) = Q(iglob)
        !if(coord(2,ibool(k,l,ispec)) == 0 .OR. &
        !     coord(2,ibool(k,l,ispec)) == 1   )  ulimit(iglob) = Q(iglob)
        
        !if(ispec == 1) WRITE(*,*) "Qispec1", k,l,ispec, Q(iglob),ulimit_temp(prod_ind)
        !if(abs(Q(iglob)) > 1d-8) then
        !if(abs(ulimit(iglob)/Q(iglob)) < 0.9999 .OR. abs(ulimit(iglob)/Q(iglob)) > 1.0001) &
        !        WRITE(*,*) coord(:,ibool_before_perio(k,l,ispec)),ispec,neighbor_DG(k,l,ispec,:),"MEGA DIF :", k,l,ispec, &
        !        ulimit(iglob), Q(iglob), abs(ulimit(iglob)/Q(iglob)), vkm1_x,vkm1_z, vkp1_x, vkp1_z, v(ispec)
        !endif
       enddo
     enddo
   enddo
   
   if(activate_total) WRITE(*,*) "SLOPE LIMITER HAS BEEN ACTIVATED", timelocal, type_var
   
   !close(10)
   !stop
   !WRITE(*,*) "MODIFIED : ", nb_mod
   ! Update solution
   Q(1:nglob_DG) = ulimit(1:nglob_DG)
   
  end subroutine SlopeLimit1
  
  subroutine SlopeLimitLin(ulimit, ispec, ul, v0, vkm1_x,vkp1_x,vkm1_z,vkp1_z,activate,type_unknown)

    use constants,only: CUSTOM_REAL,NGLLX,NGLLZ

    use specfem_par, only: ibool_before_perio, coord, &
        xix, xiz, gammax, gammaz, hprime_xx, hprime_zz!, Drx, Drz

    implicit none 

        !double precision, dimension(0:NGLOB_DG+1) :: ulimit
        integer, parameter :: NGLL = NGLLX*NGLLZ
        
        real(kind=CUSTOM_REAL), dimension(1:NGLL) :: x0, ul, xl, ux, uz, ulimit
        real(kind=CUSTOM_REAL), dimension(1:NGLL) :: z0, zl
        real(kind=CUSTOM_REAL), dimension(1:3,1)  :: minmod
        real(kind=CUSTOM_REAL) :: hx, hz, v0, vkm1_x,vkp1_x,vkm1_z,vkp1_z, &
                xixl, xizl, gammaxl, gammazl, du_dxi, du_dgamma
        real(kind=CUSTOM_REAL), dimension(1)  :: ulimit_temp
        integer :: ispec, iglob, j, k, m, n, i, type_unknown
        real(kind=CUSTOM_REAL), dimension(1:NGLLX,1:NGLLZ) :: ul_loc
        
        ! Parameters
        real(kind=CUSTOM_REAL), parameter :: ZEROl = 0._CUSTOM_REAL
        real(kind=CUSTOM_REAL), parameter :: ONEl  = 1._CUSTOM_REAL
        real(kind=CUSTOM_REAL), parameter :: TWOl  = 2._CUSTOM_REAL
        real(kind=CUSTOM_REAL), parameter :: HALFl = ONEl/TWOl
        
        real(kind=CUSTOM_REAL), parameter :: gradient_factor = ONEl/TWOl

        logical :: activate_temp, activate

        ! function ulimit = SlopeLimitLin(ul,xl,vm1,v0,vp1);
        ! Purpose: Apply slopelimited on linear function ul(Np,1) on x(Np,1)
        !          (vm1,v0,vp1) are cell averages left, center, and right
        k = 0
        do m=1,NGLLX
          do n=1,NGLLZ
                k = k + 1
                iglob = ibool_before_perio(m,n,ispec)
                xl(k) = real(coord(1,iglob), kind=CUSTOM_REAL)
                zl(k) = real(coord(2,iglob), kind=CUSTOM_REAL)
          enddo
        enddo
        
        ! Assume that x, z > 0
        hx  = maxval(xl) - minval(xl)
        !WRITE(*,*) "hx",hx
        hz  = maxval(zl) - minval(zl)
        
        ! Coordinates of the element centeroid
        x0 = minval(xl) + hx*HALFl
        z0 = minval(zl) + hz*HALFl
        
        ! Local version of the limited flux
        m = 0
        do i=1,NGLLX
          do j=1,NGLLZ
                m = m + 1
                ul_loc(i,j) = ul(m)
          enddo
        enddo
              
        ! Compute 1st order derivatives at center        
        ux = ZEROl
        uz = ZEROl
        m = 0
        do i=1,NGLLX
          do j=1,NGLLZ
                m = m + 1
                
                du_dxi    = ZEROl
                du_dgamma = ZEROl
                ! first double loop over GLL points to compute and store gradients
                ! we can merge the two loops because NGLLX == NGLLZ
                do k = 1,NGLLX
                  du_dxi    = du_dxi + ul_loc(k,j) * hprime_xx(i,k)
                  du_dgamma = du_dgamma + ul_loc(i,k) * hprime_zz(j,k)
                enddo

                xixl = xix(i,j,ispec)
                xizl = xiz(i,j,ispec)
                gammaxl = gammax(i,j,ispec)
                gammazl = gammaz(i,j,ispec)

                ! derivatives of potential
                ux(m) = du_dxi * xixl + du_dgamma * gammaxl
                uz(m) = du_dxi * xizl + du_dgamma * gammazl
                
          enddo
        enddo
        
        ulimit = v0
        if(.false.) then
        minmod(1,1) = 0
        do k=1,NGLLX*NGLLZ
                minmod(1,1) = minmod(1,1) + (1/NGLLX*NGLLZ)*ux(k)
        enddo
        endif
        minmod(1,1) = ux(NGLLX*NGLLZ/2)!+ux(NGLLZ))!0.5*(ux(1)+ux(NGLLX*(NGLLZ-1)+1))!
        minmod(2,1) = (vkp1_x-v0)/(hx*gradient_factor)
        minmod(3,1) = (v0-vkm1_x)/(hx*gradient_factor)
        
        activate = .false.
        
        call minmod_computeB(ulimit_temp, minmod, 3, 1,hx, activate_temp)
        
        if(activate_temp) activate = .true.
        
        ulimit = ulimit + (xl-x0)*ulimit_temp(1) 
        
        if(.false.) then 
        minmod(1,1) = 0
        do k=1,NGLLX*NGLLZ
                minmod(1,1) = minmod(1,1) + (1/NGLLX*NGLLZ)*uz(k)
        enddo
        endif
        minmod(1,1) = uz(NGLLX*NGLLZ/2)!*(uz(1)+uz(NGLLX*(NGLLZ-1)+1)) :0.5*(uz(1)+uz(NGLLX*(NGLLZ-1)+1))!
        minmod(2,1) = (vkp1_z-v0)/(hz*gradient_factor)
        minmod(3,1) = (v0-vkm1_z)/(hz*gradient_factor)
        
        call minmod_computeB(ulimit_temp, minmod, 3, 1,hz, activate_temp)
        
        if(activate_temp) activate = .true.
        
        ulimit = ulimit + (zl-z0)*ulimit_temp(1)
        
        ! Positivity preserving ?
        if(.false.) then
        do k=1,NGLLX*NGLLZ
        if((type_unknown == 1 .OR. type_unknown == 4) &
                 .AND. ulimit(k) < 1d-10) ulimit(k) = 1d-10
        enddo
        endif
        
  end subroutine SlopeLimitLin
  
  subroutine minmod_computeB(R, v, n, m, h, activate)

     use constants,only: CUSTOM_REAL
     use specfem_par,only: MINMOD_FACTOR

     implicit none 

        ! Usually 
        ! m = NSPEC
        ! n = 3
        integer :: n, m
        real(kind=CUSTOM_REAL), dimension(1:n,1:m) :: v
        real(kind=CUSTOM_REAL), dimension(1:m) :: R
        
        real(kind=CUSTOM_REAL) :: h, M_param
        
        real(kind=CUSTOM_REAL), parameter :: ZEROl = 0._CUSTOM_REAL
        real(kind=CUSTOM_REAL), parameter :: ONEl  = 1._CUSTOM_REAL
        
        logical activate
        
        ! function mfunc = minmod(v)
        ! Purpose: Implement the midmod function v is a vector

        M_param = MINMOD_FACTOR!0.004!1.
        
        ! Find regions where limiter is needed
        if(abs(v(1,1)) > M_param*h**2) then
                call minmod_compute(R, v, n, m)
                activate = .true.
                !WRITE(*,*) ">>>", M_param*h**2, abs(v(1,1))
                !WRITE(*,*) "ispec modified", ispec, R, v(1,1)
        else
                R = v(1,1)
                activate = .false.
        endif
        
   end subroutine minmod_computeB
  
  subroutine minmod_compute(R, v, n, m)

     use constants,only: CUSTOM_REAL

     implicit none 

        ! Usually 
        ! m = NSPEC
        ! n = 3
        integer :: n, m, k, j
        real(kind=CUSTOM_REAL), dimension(1:n,1:m) :: v
        real(kind=CUSTOM_REAL), dimension(1:m) :: R, s
        
        real(kind=CUSTOM_REAL), parameter :: ZEROl = 0._CUSTOM_REAL
        real(kind=CUSTOM_REAL), parameter :: ONEl  = 1._CUSTOM_REAL
        
        ! function mfunc = minmod(v)
        ! Purpose: Implement the midmod function v is a vector

        R = ZEROl
        !R = v(1,:)
        
        ! Find regions where limiter is needed
        !if(abs(v(1,1)) > M_param*h_param**2) then
        
        ! Sum of lines
        s = ZEROl
        ! Skip through columns
        do k=1,m
             ! Skip through lines
             do j=1,n
             s(k) = s(k) + sign(ONEl,v(j,k))/REAL(n, kind=CUSTOM_REAL)
             enddo
        enddo
        
        do k=1,m
          if(abs(s(k)) == ONEl) then
                R(k) = s(k)*minval(abs(v(:,k)))
          endif
        enddo
        
        !endif
        
   end subroutine minmod_compute
  

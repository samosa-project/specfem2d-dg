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
  
! ------------------------------------------------------------ !
! SlopeLimit1                                                  !
! ------------------------------------------------------------ !
! MUSCL slope limiter.
! See details in (Brissaud et al., 2017, Section 3.4).
! Brissaud, Q., Martin, R., Garcia, R. F., & Komatitsch, D. (2017). Hybrid Galerkin numerical modelling of elastodynamics and compressible Navier-Stokes couplings: Applications to seismo-gravito acoustic waves. Geophysical Journal International, 210(2), 1047–1069. https://doi.org/10.1093/gji/ggx185

  subroutine SlopeLimit1(Q, timelocal, type_var)
  
  use constants,only: CUSTOM_REAL,NGLLX,NGLLZ

  use specfem_par, only: nspec, nglob_DG, ibool_DG, Vandermonde, invVandermonde, &
      neighbor_DG_iface, NGLLX, NGLLZ, &
      max_interface_size, ninterface, NPROC, MPI_transfer_iface, &
      ispec_is_acoustic_DG

  implicit none 

  integer type_var
  real(kind=CUSTOM_REAL) timelocal
  real(kind=CUSTOM_REAL) :: rho_DG_P, rhovx_DG_P, rhovz_DG_P, E_DG_P, &
              veloc_x_DG_P, veloc_z_DG_P, p_DG_P, e1_DG_P
  
  integer :: iglob, iglob2, ispec, j, k, l, h, prod_ind
  integer, parameter :: NGLL = NGLLX*NGLLZ
  real(kind=CUSTOM_REAL) :: vkm1_x,vkp1_x,vkm1_z,vkp1_z
  real(kind=CUSTOM_REAL), dimension(nglob_DG) :: Q, ulimit
  real(kind=CUSTOM_REAL), dimension(NGLL) :: uh, uavg, uhl, ulimit_temp
  real(kind=CUSTOM_REAL), dimension(1:nspec,1:NGLL) :: ul
  real(kind=CUSTOM_REAL), dimension(1:nspec) :: v
  real(kind=CUSTOM_REAL), dimension(nglob_DG) :: v_MPI
  real(kind=CUSTOM_REAL), dimension(NGLLX*max_interface_size,ninterface) :: buffer_v
  
  integer :: ipoin,num_interface, count_bound
  logical is_ispec_onbound
  
  logical :: activate, activate_total
  
  ! Parameters
  real(kind=CUSTOM_REAL), parameter :: ZEROl = 0._CUSTOM_REAL
  
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
            uh(k) = uh(k) + &
                    invVandermonde(k,prod_ind)*Q(iglob2)
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
    do k=1,NGLLX
      do l=1,NGLLZ
        v_MPI(ibool_DG(k,l,ispec)) = v(ispec)
      enddo
    enddo
  enddo ! do ispec=1,NSPEC
 
  ! find cell averages
  if(NPROC > 1) &
  call assemble_MPI_vector_DG(v_MPI, buffer_v)

  ulimit = 1.

  activate_total = .false.
  !nb_mod = 0
  do ispec=1,nspec
  
    if(.not. ispec_is_acoustic_DG(ispec)) cycle

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

    count_bound = 0

    ! apply slope limiter to selected elements
    ! Find neighbors
    !vkm1_x = 0.
    !do k=1,NGLLZ
    !!vkm1_x = vkm1_x + (1./NGLLZ)*Q(ibool_DG(1,k,ispec))!v(ispec)
    !enddo
    !vkm1_x = Q(ibool_DG(1,NGLLZ/2,ispec))
    !vkm1_x = v(ispec)
    !if(neighbor_DG(1,NGLLZ/2,ispec,3) > - 1) then
    if(neighbor_DG_iface(NGLLX/2,1,ispec,3) > - 1) then
      vkm1_x = v(neighbor_DG_iface(NGLLX/2,1,ispec,3))
      is_ispec_onbound = .false.
      count_bound = count_bound + 1
    endif

    !vkp1_x = 0.
    !do k=1,NGLLZ
    !vkp1_x = vkp1_x + (1./NGLLZ)*Q(ibool_DG(NGLLX,k,ispec))!v(ispec)
    !enddo
    !vkp1_x = Q(ibool_DG(NGLLX,NGLLZ/2,ispec))!v(ispec)
    !vkp1_x = v(ispec)
    !if(neighbor_DG(NGLLX,NGLLZ/2,ispec,3) > - 1) then
    if(neighbor_DG_iface(NGLLX/2,2,ispec,3) > - 1) then
      vkp1_x = v(neighbor_DG_iface(NGLLX/2,2,ispec,3))
      is_ispec_onbound = .false.
      count_bound = count_bound + 1
    endif

    !vkm1_z = 0.
    !do k=1,NGLLX
    !vkm1_z = vkm1_z + (1./NGLLX)*Q(ibool_DG(k,1,ispec))!v(ispec)
    !enddo
    !vkm1_z = Q(ibool_DG(NGLLX/2,1,ispec))!v(ispec)
    !vkm1_z = v(ispec)
    !if(neighbor_DG(NGLLX/2,1,ispec,3) > - 1) then
    if(neighbor_DG_iface(NGLLX/2,3,ispec,3) > - 1) then
      vkm1_z = v(neighbor_DG_iface(NGLLX/2,3,ispec,3))
      is_ispec_onbound = .false.
      count_bound = count_bound + 1
    endif

    !vkp1_z = 0.
    !do k=1,NGLLX
    !vkp1_z = vkp1_z + (1./NGLLX)*Q(ibool_DG(k,NGLLZ,ispec))!v(ispec)
    !enddo
    !vkp1_z = Q(ibool_DG(NGLLX/2,NGLLZ,ispec))!v(ispec)
    !vkp1_z = v(ispec)
    !if(neighbor_DG(NGLLX/2,NGLLZ,ispec,3) > - 1) then
    if(neighbor_DG_iface(NGLLX/2,4,ispec,3) > - 1) then
      vkp1_z = v(neighbor_DG_iface(NGLLX/2,4,ispec,3))
      is_ispec_onbound = .false.
      count_bound = count_bound + 1
    endif
    
    is_ispec_onbound = .true.
    if(count_bound == 4) is_ispec_onbound = .false.
    
    count_bound = 0
    ! Left neighbour.
    ipoin         = -1
    num_interface = -1
    if(NPROC > 1) then
      iglob         = ibool_DG(1,NGLLZ/2,ispec)
      ipoin         = MPI_transfer_iface(NGLLX/2, 1, ispec, 1)
      num_interface = MPI_transfer_iface(NGLLX/2, 1, ispec, 2)
      if(ipoin > -1) then
        vkm1_x = buffer_v(ipoin,num_interface)
        is_ispec_onbound = .false.
        count_bound = count_bound + 1
      endif
    endif
    
    ! Right neighbour.
    ipoin         = -1
    num_interface = -1
    if(NPROC > 1) then
      iglob         = ibool_DG(NGLLX,NGLLZ/2,ispec)
      ipoin         = MPI_transfer_iface(NGLLX/2, 2, ispec, 1)
      num_interface = MPI_transfer_iface(NGLLX/2, 2, ispec, 2)
      if(ipoin > -1) then
        vkp1_x = buffer_v(ipoin,num_interface)
        is_ispec_onbound = .false.
        count_bound = count_bound + 1
      endif
    endif
    
    ! Bottom neighbour.
    ipoin         = -1
    num_interface = -1
    if(NPROC > 1) then
    iglob         = ibool_DG(NGLLX/2, 1, ispec)
    ipoin         = MPI_transfer_iface(NGLLX/2, 3, ispec, 1)
    num_interface = MPI_transfer_iface(NGLLX/2, 3, ispec, 2)
    if(ipoin > -1) then
        vkm1_z = buffer_v(ipoin,num_interface)
        is_ispec_onbound = .false.
        count_bound = count_bound + 1
    endif
    endif
    
    ! Top neighbour.
    ipoin         = -1
    num_interface = -1
    if(NPROC > 1) then
    iglob         = ibool_DG(NGLLX/2, NGLLZ, ispec)
    ipoin         = MPI_transfer_iface(NGLLX/2, 4, ispec, 1)
    num_interface = MPI_transfer_iface(NGLLX/2, 4, ispec, 2)
    if(ipoin > -1) then 
        vkp1_z = buffer_v(ipoin,num_interface)
        is_ispec_onbound = .false.
        count_bound = count_bound + 1
    endif
    endif
    
    !is_ispec_onbound = .true.
    if(count_bound == 4) is_ispec_onbound = .false.

    ! Compute Limited flux
    activate = .false.
    if(.not. is_ispec_onbound) then
      call SlopeLimitLin(ulimit_temp, ispec, ul(ispec,:),v(ispec),&
                         vkm1_x,vkp1_x,vkm1_z,vkp1_z,activate,type_var)
    
    endif

    if(activate) activate_total = activate

    prod_ind = 0
    do k=1,NGLLX
      do l=1,NGLLZ
        iglob = ibool_DG(k,l,ispec)
        prod_ind = prod_ind + 1
        if(.not. is_ispec_onbound .AND. activate) then
          ulimit(iglob) = ulimit_temp(prod_ind)
        endif

        ! No slope limiter on the outer boundary
        !if(ispec_bound(k,l))  ulimit(iglob) = Q(iglob)
        if(.not. activate) then
          ulimit(iglob) = Q(iglob)
        endif

      enddo
    enddo
  enddo

  if(activate_total) WRITE(*,*) "SLOPE LIMITER HAS BEEN ACTIVATED", timelocal, type_var

  Q(1:nglob_DG) = ulimit(1:nglob_DG)
   
end subroutine SlopeLimit1

! ------------------------------------------------------------ !
! SlopeLimitLin                                                !
! ------------------------------------------------------------ !
! Linear slope limiter.
! Called by SlopeLimit1.

subroutine SlopeLimitLin(ulimit, ispec, ul, v0, vkm1_x,vkp1_x,vkm1_z,vkp1_z,activate,type_unknown)

  use constants,only: CUSTOM_REAL,NGLLX,NGLLZ
  
  use specfem_par, only: ibool_before_perio, coord, xix, xiz, gammax, gammaz, hprime_xx, hprime_zz
  
  implicit none

  real(kind=CUSTOM_REAL), parameter :: ZEROl = 0._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: ONEl  = 1._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: TWOl  = 2._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: HALFl = ONEl/TWOl
  integer, parameter :: NGLL = NGLLX*NGLLZ
  real(kind=CUSTOM_REAL), dimension(1:NGLL) :: x0, ul, xl, ux, uz, ulimit
  real(kind=CUSTOM_REAL), dimension(1:NGLL) :: z0, zl
  real(kind=CUSTOM_REAL), dimension(1:3,1)  :: minmod
  real(kind=CUSTOM_REAL) :: hx, hz, v0, vkm1_x,vkp1_x,vkm1_z,vkp1_z, &
                            xixl, xizl, gammaxl, gammazl, du_dxi, du_dgamma
  real(kind=CUSTOM_REAL), dimension(1)  :: ulimit_temp
  integer :: ispec, iglob, j, k, m, n, i, type_unknown
  real(kind=CUSTOM_REAL), dimension(1:NGLLX,1:NGLLZ) :: ul_loc
  real(kind=CUSTOM_REAL), parameter :: gradient_factor = ONEl
  logical :: activate_temp, activate
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
  minmod(1,1) = uz(NGLLX*NGLLZ/2)
  minmod(2,1) = (vkp1_z-v0)/(hz*gradient_factor)
  minmod(3,1) = (v0-vkm1_z)/(hz*gradient_factor)
  
  call minmod_computeB(ulimit_temp, minmod, 3, 1, hz, activate_temp)
  
  if(activate_temp) activate = .true.
  
  ulimit = ulimit + (zl-z0)*ulimit_temp(1)
  
  ! Positivity preserving?
  if(.false.) then
    do k=1,NGLLX*NGLLZ
      if((type_unknown == 1 .OR. type_unknown == 4) &
         .AND. ulimit(k) < 1d-10) ulimit(k) = 1d-10
    enddo
  endif
        
end subroutine SlopeLimitLin


! ------------------------------------------------------------ !
! minmod_compute                                               !
! ------------------------------------------------------------ !
! Another minmod function.

subroutine minmod_compute(R, v, n, m)
  use constants, only: CUSTOM_REAL
  implicit none
  real(kind=CUSTOM_REAL), parameter :: ZEROl = 0._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: ONEl  = 1._CUSTOM_REAL
  integer :: n, m, k, j
  real(kind=CUSTOM_REAL), dimension(1:n,1:m) :: v
  real(kind=CUSTOM_REAL), dimension(1:m) :: R, s
  R = ZEROl
  s = ZEROl
  do k=1,m
    do j=1,n
      s(k) = s(k) + sign(ONEl,v(j,k))/REAL(n, kind=CUSTOM_REAL)
    enddo
  enddo
  do k=1,m
    if(abs(s(k)) == ONEl) then
      R(k) = s(k)*minval(abs(v(:,k)))
    endif
  enddo
end subroutine minmod_compute


! ------------------------------------------------------------ !
! minmod_computeB                                              !
! ------------------------------------------------------------ !
! Another minmod function.
subroutine minmod_computeB(R, v, n, m, h, activate)
  use constants, only: CUSTOM_REAL
  use specfem_par, only: MINMOD_FACTOR
  implicit none
  integer :: n, m
  real(kind=CUSTOM_REAL), dimension(1:n,1:m) :: v
  real(kind=CUSTOM_REAL), dimension(1:m) :: R
  real(kind=CUSTOM_REAL) :: h, M_param
  real(kind=CUSTOM_REAL), parameter :: ZEROl = 0._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: ONEl  = 1._CUSTOM_REAL
  logical :: activate
  M_param = MINMOD_FACTOR
  ! Find regions where limiter is needed
  if(abs(v(1,1)) > M_param*h**2) then
    call minmod_compute(R, v, n, m)
    activate = .true.
  else
    R = v(1,1)
    activate = .false.
  endif
end subroutine minmod_computeB
  

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


  subroutine compute_forces_acoustic_DG(rho_DG_main, rhovx_DG_main, rhovz_DG_main, E_DG_main, &
        dot_rho, dot_rhovx, dot_rhovz, dot_E, timelocal)


! compute forces in the acoustic elements in forward simulation and in adjoint simulation in adjoint inversion

  use constants,only: CUSTOM_REAL,NGLLX,NGLLZ,gamma_euler

  use specfem_par, only: nglob_DG,nspec,ispec_is_acoustic, &
                         xix,xiz,gammax,gammaz,jacobian, &
                         hprimewgll_xx, &
                         hprimewgll_zz,wxgll,wzgll, &
                         normal_DG, normal_DG_corner, ibool_DG, weight_DG, weight_DG_corner, &
                         neighbor_DG, ispec_is_acoustic_forcing, ispec_is_acoustic_surface, &
                         ACOUSTIC_FORCING, neighbor_DG_corner, is_corner, &
                         ispec_is_acoustic_surface_corner, it, &
                         potential_dphi_dx_DG, potential_dphi_dz_DG, ibool, &
                         hprime_xx, hprime_zz, elastic_tensor, ispec_is_acoustic_coupling_el, veloc_elastic,&
                         coord!, &
                         !rho0_DG, E0_DG!, &!this_ibool_is_a_periodic_edge, &
                         !ibool, coord, ibool_before_perio!, ADD_PERIODIC_CONDITIONS, rmass_inverse_acoustic_DG! &
                         !,coord, ibool!, link_DG_CG

  implicit none

  real(kind=CUSTOM_REAL), dimension(nglob_DG) :: rho_DG_main, rhovx_DG_main, rhovz_DG_main, E_DG_main
  real(kind=CUSTOM_REAL), dimension(nglob_DG) :: rho_DG, rhovx_DG, rhovz_DG, E_DG
  real(kind=CUSTOM_REAL), dimension(nglob_DG) :: dot_rho, dot_rhovx, dot_rhovz, dot_E
  
  ! local parameters
  integer :: ispec,i,j,k,iglob, it_corner
  integer :: ifirstelem,ilastelem

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: temp_rho_1, temp_rho_2, &
        temp_rhovx_1, temp_rhovx_2, temp_rhovz_1, temp_rhovz_2, &
        temp_E_1, temp_E_2, &
        temp_rho_gravi, temp_rhovx_gravi, temp_rhovz_gravi, temp_E_gravi

  ! Jacobian matrix and determinant
  real(kind=CUSTOM_REAL) :: xixl,xizl,gammaxl,gammazl,jacobianl

  real(kind=CUSTOM_REAL) :: lambda, nx, nz, weight, &
        temp_unknown_M, temp_unknown_P, temp_unknown2_M, temp_unknown2_P, &
        temp_unknown, temp_unknown2, &
        flux_x, flux_z, flux_n, jump, &
        rho_DG_P, veloc_x_DG_P, veloc_z_DG_P, &
        E_DG_P, p_DG_P, rhovx_DG_P, rhovz_DG_P, timelocal!, &
        !templ_rhovx_gravi, templ_rhovz_gravi, templ_rho_gravi, templ_E_gravi
        
  integer :: ind, iglobM, iglobP
  
  integer, dimension(3) :: neighbor
  
  ! Local
  real(kind=CUSTOM_REAL), dimension(nglob_DG) :: veloc_x_DG, veloc_z_DG, p_DG
  integer :: i_el, j_el, ispec_el
  real(kind=CUSTOM_REAL) :: veloc_x, veloc_z
  
  ! Viscosity
  real(kind=CUSTOM_REAL) :: dux_dxi, dux_dgamma, duz_dxi, duz_dgamma
  real(kind=CUSTOM_REAL) :: dux_dx, dux_dz, duz_dx, duz_dz
  real(kind=CUSTOM_REAL) :: wxl, wzl
  real(kind=CUSTOM_REAL), parameter :: mu_visco = 0.00
  
  ! Parameters
  real(kind=CUSTOM_REAL), parameter :: ZERO = 0._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: ONE  = 1._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: TWO  = 2._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: HALF = 0.5_CUSTOM_REAL
  
  real(kind=CUSTOM_REAL), parameter :: threshold = 0.000000001_CUSTOM_REAL
  
  ! Temporary
  real(kind=CUSTOM_REAL), dimension(nglob_DG) :: rho0_DG, E0_DG, p0_DG
  logical, parameter :: ONLY_PERTURBATION = .false.
  
  ifirstelem = 1
  ilastelem = nspec
  
  if(ONLY_PERTURBATION) then
        
        rho0_DG  = 1.
        p0_DG   = 1.
        E0_DG    = p0_DG/(gamma_euler-1.)
  
        rho_DG   = rho0_DG + rho_DG_main
        rhovx_DG = rhovx_DG_main
        rhovz_DG = rhovz_DG_main
        E_DG     = E0_DG + E_DG_main
        
  else
  
        rho_DG   = rho_DG_main
        rhovx_DG = rhovx_DG_main
        rhovz_DG = rhovz_DG_main
        E_DG     = E_DG_main
        
  endif
  
  ! Init auxiliary unknwons
  veloc_x_DG = rhovx_DG/rho_DG
  veloc_z_DG = rhovz_DG/rho_DG
  p_DG       = (gamma_euler - ONE)*( E_DG &
        - (HALF)*rho_DG*( veloc_x_DG**2 + veloc_z_DG**2 ) )
  ! Test with potential temperature instead E_DG = Theta_DG
  !p0_DG      = 1.
  !p_DG       = p0_DG * (2.*E_DG/p0_DG)**(gamma_euler)
  
  if(ONLY_PERTURBATION) then
        !p_DG = (gamma_euler - ONE)*( E_DG &
        !        - (HALF*rho_DG)*( (rhovx_DG/rho_DG)**2 &
        !        + (rhovz_DG/rho_DG)**2 ) ) - p0_DG 
        p_DG = (gamma_euler - ONE)*( E_DG_main &
                - (HALF*rho_DG_main)*( veloc_x_DG**2 &
                + veloc_z_DG**2 ) )  
  endif
        
  ! Initialization
  dot_rho   = ZERO
  dot_rhovx = ZERO
  dot_rhovz = ZERO
  dot_E     = ZERO
        
  WRITE(*,*) it,"MAXVAL ", maxval(rho_DG), minval(rho_DG)
  WRITE(*,*) it,"MAXVAL ", maxval(rhovx_DG), minval(rhovx_DG)
  WRITE(*,*) it,"MAXVAL ", maxval(rhovz_DG), minval(rhovz_DG)
  WRITE(*,*) it,"MAXVAL ", maxval(E_DG), minval(E_DG)
  !WRITE(*,*) "hprime_xx", hprime_xx
  !WRITE(*,*) "hprime_zz", hprime_zz

  ! Test Source
  !WRITE(*,*) "TESST", size(ibool_DG(1,:,1)), ibool_DG(2,2,45)
  !p_DG(ibool_DG(3,3,45)) = p_DG(ibool_DG(2,2,45)) + 1000.*exp(-((timelocal-0.05)/0.1)**2)

! loop over spectral elements
  do ispec = ifirstelem,ilastelem

    ! acoustic spectral element
    if (ispec_is_acoustic(ispec)) then

      ! first double loop over GLL points to compute and store gradients
      do j = 1,NGLLZ
        do i = 1,NGLLX
          
          iglob = ibool_DG(i,j,ispec)
          
          !p_DG(iglob) = p_DG(iglob) &
        !- (gamma_euler - ONE)*rho_DG(iglob)*potential_dphi_dz_DG(ibool(i,j,ispec))*coord(2, ibool(i,j,ispec))
          !WRITE(*,*) "potential_dphi_dz_DG(ibool(i,j,ispec))",potential_dphi_dz_DG(ibool(i,j,ispec))

          jacobianl = real(jacobian(i,j,ispec), kind=CUSTOM_REAL)
        
          xixl = real(xix(i,j,ispec), kind=CUSTOM_REAL)
          xizl = real(xiz(i,j,ispec), kind=CUSTOM_REAL)
          gammaxl = real(gammax(i,j,ispec), kind=CUSTOM_REAL)
          gammazl = real(gammaz(i,j,ispec), kind=CUSTOM_REAL)
          
          wzl = real(wzgll(j), kind=CUSTOM_REAL)
          wxl = real(wxgll(i), kind=CUSTOM_REAL)
        
          !!!!!!!!!!!!!!!!!!!!!!!!
          ! Inviscid stress tensor
          temp_unknown = rhovx_DG(iglob)
          temp_unknown2 = rhovz_DG(iglob)
          temp_rho_1(i,j) = wzl * jacobianl * (xixl * temp_unknown + xizl * temp_unknown2) 
          temp_rho_2(i,j) = wxl * jacobianl * (gammaxl * temp_unknown + gammazl * temp_unknown2) 
          
          temp_unknown = rho_DG(iglob)*veloc_x_DG(iglob)**2 + p_DG(iglob)
          ! MODIF GRAVI
          !temp_unknown = temp_unknown &
          !      - rho_DG(iglob)*potential_dphi_dx_DG(ibool(i,j,ispec))*coord(2,ibool(i,j,ispec))
          temp_unknown2 = rho_DG(iglob)*veloc_x_DG(iglob)*veloc_z_DG(iglob)
          temp_rhovx_1(i,j) = wzl * jacobianl * (xixl * temp_unknown + xizl * temp_unknown2) 
          temp_rhovx_2(i,j) = wxl * jacobianl * (gammaxl * temp_unknown + gammazl * temp_unknown2) 
          !WRITE(*,*) "p_DG(iglob)", p_DG(iglob)
          
          temp_unknown = rho_DG(iglob)*veloc_x_DG(iglob)*veloc_z_DG(iglob)
          temp_unknown2 = rho_DG(iglob)*veloc_z_DG(iglob)**2 + p_DG(iglob)
          !temp_unknown2 = temp_unknown2 &
          !      - rho_DG(iglob)*potential_dphi_dx_DG(ibool(i,j,ispec))*coord(2,ibool(i,j,ispec))
          temp_rhovz_1(i,j) = wzl * jacobianl * (xixl * temp_unknown + xizl * temp_unknown2) 
          temp_rhovz_2(i,j) = wxl * jacobianl * (gammaxl * temp_unknown + gammazl * temp_unknown2) 
          
          temp_unknown = veloc_x_DG(iglob)*(E_DG(iglob) + p_DG(iglob))
          !if(ONLY_PERTURBATION) temp_unknown = temp_unknown + veloc_x_DG(iglob)*p0_DG(iglob)
          temp_unknown2 = veloc_z_DG(iglob)*(E_DG(iglob) + p_DG(iglob))
          !if(ONLY_PERTURBATION) temp_unknown2 = temp_unknown2 + veloc_z_DG(iglob)*p0_DG(iglob)
          temp_E_1(i,j) = wzl * jacobianl * (xixl * temp_unknown + xizl * temp_unknown2) 
          temp_E_2(i,j) = wxl * jacobianl * (gammaxl * temp_unknown + gammazl * temp_unknown2) 
          
          ! Test with potential temperature instead E_DG = rho_DG*Theta_DG
          !temp_unknown  = veloc_x_DG(iglob)*E_DG(iglob)
          !temp_unknown2 = veloc_z_DG(iglob)*E_DG(iglob)
          !temp_E_1(i,j) = wzl * jacobianl * (xixl * temp_unknown + xizl * temp_unknown2) 
          !temp_E_2(i,j) = wxl * jacobianl * (gammaxl * temp_unknown + gammazl * temp_unknown2)
          
          !!!!!!!!!!!!!!!!!!!!!!!
          ! Viscous stress tensor
          dux_dxi    = ZERO
          dux_dgamma = ZERO
          duz_dxi    = ZERO
          duz_dgamma = ZERO
          ! first double loop over GLL points to compute and store gradients
          ! we can merge the two loops because NGLLX == NGLLZ
          do k = 1,NGLLX
                dux_dxi    = dux_dxi    + veloc_x_DG(ibool_DG(k,j,ispec)) * real(hprime_xx(i,k), kind=CUSTOM_REAL)
                dux_dgamma = dux_dgamma + veloc_x_DG(ibool_DG(i,k,ispec)) * real(hprime_zz(j,k), kind=CUSTOM_REAL)
                duz_dxi    = duz_dxi    + veloc_z_DG(ibool_DG(k,j,ispec)) * real(hprime_xx(i,k), kind=CUSTOM_REAL)
                duz_dgamma = duz_dgamma + veloc_z_DG(ibool_DG(i,k,ispec)) * real(hprime_zz(j,k), kind=CUSTOM_REAL)
          enddo

          ! derivatives of potential
          dux_dx = ( dux_dxi * xixl + dux_dgamma * gammaxl )
          dux_dz = ( dux_dxi * xizl + dux_dgamma * gammazl )
          duz_dx = ( duz_dxi * xixl + duz_dgamma * gammaxl )
          duz_dz = ( duz_dxi * xizl + duz_dgamma * gammazl )
          
          temp_unknown = mu_visco*( TWO*dux_dx - (TWO/3.)*(dux_dx + duz_dz) )
          temp_unknown2 = mu_visco*( dux_dz + duz_dx )
          temp_rhovx_1(i,j) = temp_rhovx_1(i,j) - wzl * jacobianl * (xixl * temp_unknown + xizl * temp_unknown2) 
          temp_rhovx_2(i,j) = temp_rhovx_2(i,j) - wxl * jacobianl * (gammaxl * temp_unknown + gammazl * temp_unknown2) 
          
          temp_E_1(i,j) = temp_E_1(i,j) - wzl * jacobianl * (xixl * veloc_x_DG(iglob) * temp_unknown &
                + xizl * veloc_z_DG(iglob) * temp_unknown2) 
          
          temp_unknown = mu_visco*( dux_dz + duz_dx )
          temp_unknown2 = mu_visco*( TWO*duz_dz - (TWO/3.)*(dux_dx + duz_dz) )
          temp_rhovz_1(i,j) = temp_rhovz_1(i,j) - wzl * jacobianl * (xixl * temp_unknown + xizl * temp_unknown2) 
          temp_rhovz_2(i,j) = temp_rhovz_2(i,j) - wxl * jacobianl * (gammaxl * temp_unknown + gammazl * temp_unknown2) 
          
          temp_E_2(i,j) = temp_E_2(i,j) + wxl * jacobianl * (gammaxl * veloc_x_DG(iglob) * temp_unknown &
                + gammazl * veloc_z_DG(iglob) * temp_unknown2) 
          
          ! Gravity potentials
          temp_rho_gravi(i,j)   = ZERO
          temp_rhovx_gravi(i,j) = -rho_DG(iglob)*potential_dphi_dx_DG(ibool(i,j,ispec))* jacobianl
          temp_rhovz_gravi(i,j) = -rho_DG(iglob)*potential_dphi_dz_DG(ibool(i,j,ispec))* jacobianl
          temp_E_gravi(i,j)     = -rho_DG(iglob)*(veloc_x_DG(iglob)*potential_dphi_dx_DG(ibool(i,j,ispec)) + &
                veloc_z_DG(iglob)*potential_dphi_dz_DG(ibool(i,j,ispec)))* jacobianl
                !WRITE(*,*) "potential_dphi_dx_DG(ibool(i,j,ispec))", potential_dphi_dx_DG(ibool(i,j,ispec))
          
        enddo
      enddo
      !WRITE(*,*) "dot_rhovx1", maxval(dot_rhovx), minval(dot_rhovx)
!    
! second double-loop over GLL to compute all the terms
!
      do j = 1,NGLLZ
        do i = 1,NGLLX
          iglob = ibool_DG(i,j,ispec)
          ! along x direction and z direction
          ! and assemble the contributions
            do k = 1,NGLLX
            
              dot_rho(iglob) = dot_rho(iglob) + &
                       (temp_rho_1(k,j) * real(hprimewgll_xx(k,i), kind=CUSTOM_REAL) + &
                       temp_rho_2(i,k) * real(hprimewgll_zz(k,j), kind=CUSTOM_REAL))
                       
              dot_rhovx(iglob) = dot_rhovx(iglob) + &
                       (temp_rhovx_1(k,j) * real(hprimewgll_xx(k,i), kind=CUSTOM_REAL) + &
                       temp_rhovx_2(i,k) * real(hprimewgll_zz(k,j), kind=CUSTOM_REAL))
              
              !WRITE(*,*) i,j,ispec,"hprimewgll_xx(k,i)", hprimewgll_xx(k,i)*wzgll(j)
              !WRITE(*,*) "hprimewgll_zz(k,i)", hprimewgll_zz(k,j)*wzgll(i)
                       
              dot_rhovz(iglob) = dot_rhovz(iglob) + &
                       (temp_rhovz_1(k,j) * real(hprimewgll_xx(k,i), kind=CUSTOM_REAL) + &
                       temp_rhovz_2(i,k) * real(hprimewgll_zz(k,j), kind=CUSTOM_REAL))
              
              dot_E(iglob) = dot_E(iglob) + &
                       (temp_E_1(k,j) * real(hprimewgll_xx(k,i), kind=CUSTOM_REAL) + &
                       temp_E_2(i,k) * real(hprimewgll_zz(k,j), kind=CUSTOM_REAL))
                       
            enddo
            
            wzl = real(wzgll(j), kind=CUSTOM_REAL)
            wxl = real(wxgll(i), kind=CUSTOM_REAL)
            
            dot_rho(iglob)   = dot_rho(iglob)   + temp_rho_gravi(i,j) * wxl * wzl
            dot_rhovx(iglob) = dot_rhovx(iglob) + temp_rhovx_gravi(i,j) * wxl * wzl
            dot_rhovz(iglob) = dot_rhovz(iglob) + temp_rhovz_gravi(i,j) * wxl * wzl
            dot_E(iglob)     = dot_E(iglob)     + temp_E_gravi(i,j) * wxl * wzl

        enddo ! second loop over the GLL points
      enddo
      
      !WRITE(*,*) "1 - dot_rhovx2", maxval(dot_rhovx), minval(dot_rhovx)
      
    it_corner = 0

    do  i = 1,NGLLX
      
      j = 1
      do while (j <= NGLLZ)
      
        ! We skip interior points
        if(i > 1 .AND. i < NGLLX .AND. j > 1 .AND. j < NGLLZ) then
                j = j + 1
                cycle
        endif
        
        ! Recover neighbor location
        neighbor = neighbor_DG(i,j,ispec,:)
      
        !WRITE(*,*) "CURRENT WORK : ", i, j, ispec, iglobM, iglobP, &
        !neighbor(3),neighbor(2),neighbor(1)
      
        ! Recover normal and right coeficients
        if(j == 1) ind = 1
        if(i == 1) ind = 2
        if(j == NGLLZ) ind = 3
        if(i == NGLLX) ind = 4
        
        if(is_corner(i,j) .AND. it_corner == 0) then
                if(j == 1 .AND. i == 1) ind = 1
                if(j == NGLLZ .AND. i == 1) ind = 2
                if(j == 1 .AND. i == NGLLX) ind = 1
                if(j == NGLLZ .AND. i == NGLLX) ind = 3
        endif
        
        ! Needs x2 points at corners to correctly map edges
        ! => 2 step on the same point
        if(it_corner == 1) then
                neighbor = neighbor_DG_corner(i,j,ispec,:)
                if(j == 1) ind = 1
                if(i == 1) ind = 2
                if(j == NGLLZ) ind = 3
                if(i == NGLLX) ind = 4
                it_corner = 2
        endif
        
        !if(is_corner(i,j) .AND. neighbor_DG_corner(i,j,ispec,3) == neighbor_DG(i,j,ispec,3)) then
        !        WRITE(*,*) "i, j, ispec", i, j, ispec
        !        WRITE(*,*) "neighbor corner",neighbor_DG_corner(i,j,ispec,:)
        !        WRITE(*,*) "neighbor", neighbor_DG(i,j,ispec,:)
        !        WRITE(*,*) "nx/nz", nx, nz
        !        WRITE(*,*) "***********", it_corner
        !endif
        
        !nx = normal_DG(ispec, ind, 1)
        !nz = normal_DG(ispec, ind, 2)
        if((is_corner(i,j) .AND. it_corner == 2) .OR. .not. is_corner(i,j)) then
                nx     = normal_DG(i,j,ispec, 1)
                nz     = normal_DG(i,j,ispec, 2)
                weight = weight_DG(i, j, ispec)
                !if(is_corner(i,j) .AND. it_corner == 2) then
                !WRITE(*,*) "i, j, ispec", i, j, ispec
                !WRITE(*,*) "neighbor", neighbor
                !WRITE(*,*) "nx/nz", nx, nz
                !WRITE(*,*) "***********", weight_DG(i, j, ispec)
                !endif
        else
                nx     = normal_DG_corner(i,j,ispec, 1)
                nz     = normal_DG_corner(i,j,ispec, 2)
                weight = weight_DG_corner(i, j, ispec)
                
                !WRITE(*,*) "i, j, ispec", i, j, ispec
                !WRITE(*,*) "neighbor", neighbor
                !WRITE(*,*) "nx/nz", nx, nz
                !WRITE(*,*) "***********", weight
        endif
        
        
        
        ! Interior point
        iglobM = ibool_DG(i,j,ispec)
        
        ! If at corner notify that we will need to go again 
        if( is_corner(i,j) .AND. it_corner == 0) it_corner = 1
        
        !if(ispec == 6400 .OR. ispec == 1) then
        !!WRITE(*,*) "TEST : ", iglobP, iglobM, i, j, ispec, ispec_is_acoustic_surface(i,j,ispec), &
        !!nx, nz, weight
        !WRITE(*,*) "coord : ", coord(:,ibool(i,j,ispec)), nx, nz, weight, i, j, ispec
        !WRITE(*,*) "iglobM : ", iglobM 
        !WRITE(*,*) "ispec : ", i,j,ispec , nx, nz, weight
        !WRITE(*,*) "******"
        !endif
        rho_DG_P     = ZERO
        rhovx_DG_P   = ZERO
        rhovz_DG_P   = ZERO
        E_DG_P       = ZERO
        veloc_x_DG_P = ZERO
        veloc_z_DG_P = ZERO
        p_DG_P       = ZERO
        ! Define neighbor unknown in order to easily treat B.C.
        ! If surface element => we impose P value
        if ( (ispec_is_acoustic_surface(i,j,ispec) .AND. it_corner < 2) .OR. &
                (ispec_is_acoustic_surface_corner(i,j,ispec) .AND. it_corner == 2) ) then
        ! Forcing boundary conditions
        if(ACOUSTIC_FORCING .AND. ispec_is_acoustic_forcing(i,j,ispec)) then
        
                rho_DG_P     = ONE
                veloc_x_DG_P = ZERO
                veloc_z_DG_P = ZERO
                p_DG_P       = ONE
                E_DG_P       = p_DG_P/(gamma_euler - ONE) &
                        + HALF*rho_DG_P*( veloc_x_DG_P**2 + veloc_z_DG_P**2 )
                rhovx_DG_P   = rho_DG_P*veloc_x_DG_P
                rhovz_DG_P   = rho_DG_P*veloc_z_DG_P
                
        ! Elastic coupling        
        elseif(ispec_is_acoustic_coupling_el(i,j,ispec,3) >= 0 .AND. it_corner < 2) then
        
               !WRITE(*,*) "ispec_is_acoustic_coupling_el", i,j,ispec, ispec_is_acoustic_coupling_el(i,j,ispec,:), &
               ! nx, nz, weight  
                
               ! Coordinates of elastic element
               i_el     = ispec_is_acoustic_coupling_el(i,j,ispec,1)
               j_el     = ispec_is_acoustic_coupling_el(i,j,ispec,2)
               ispec_el = ispec_is_acoustic_coupling_el(i,j,ispec,3)
               
               !WRITE(*,*) "ispec_is_acoustic_coupling_el", i,j,ispec,coord(:,ibool(i,j,ispec)), nx, nz, &
               ! elastic_tensor(i_el,j_el,ispec_el,4)*nz
               
               iglob = ibool(i_el,j_el,ispec_el)
               ! Elastic velocities
               veloc_x = veloc_elastic(1,iglob)
               veloc_z = veloc_elastic(2,iglob)
        
               rho_DG_P = 1.!rho_DG(iglobM)
               
               ! Continuity of the velocity
               veloc_x_DG_P = 1.*veloc_x
               veloc_z_DG_P = 1.*veloc_z!-0.2*((timelocal-0.1)/0.1)*(1/0.1)*exp(-((timelocal-0.1)/0.1)**2)!
               
               p_DG_P = 1.!-elastic_tensor(i_el,j_el,ispec_el,4)
               
               !WRITE(*,*) "VELOC ", veloc_x, veloc_z
               
               ! Continuity of the normal velocity
               !veloc_x_DG_P = 0!veloc_x
               !veloc_z_DG_P = veloc_z_DG(iglobM) + (1/nz)*(veloc_x*nx + veloc_z*nz)
               
               ! Continuity of the normal stress tensor (w/ non linear terms)
               !p_DG_P       = (-1/(nx + nz))*( elastic_tensor(i_el,j_el,ispec_el,1)*nx &
               ! + elastic_tensor(i_el,j_el,ispec_el,2)*nz &
               ! + elastic_tensor(i_el,j_el,ispec_el,3)*nx &
               ! + elastic_tensor(i_el,j_el,ispec_el,4)*nz &
               ! - nx*veloc_x_DG_P**2 - nz*veloc_z_DG_P**2 &! &
               ! - (rho_DG_P*veloc_x_DG_P*veloc_z_DG_P)*(nx + nz) )
               
               ! (w/out non linear terms)
               !p_DG_P = p_DG(iglobM) -(1/nz)*( elastic_tensor(i_el,j_el,ispec_el,3)*nx &
               ! + elastic_tensor(i_el,j_el,ispec_el,4)*nz ) + p0_DG(iglobM)
               
               E_DG_P       = p_DG_P/(gamma_euler - 1.) + HALF*rho_DG_P*( veloc_x_DG_P**2 + veloc_z_DG_P**2 )
               rhovx_DG_P   = rho_DG_P*veloc_x_DG_P
               rhovz_DG_P   = rho_DG_P*veloc_z_DG_P
               
               if(ONLY_PERTURBATION) then
                        rho_DG_P   = 0.
                        rhovx_DG_P = rho0_DG(iglobM)*veloc_x_DG_P
                        rhovz_DG_P = rho0_DG(iglobM)*veloc_z_DG_P
                        E_DG_P     = (p_DG_P)/(gamma_euler - 1.) + HALF*rho_DG_P*( veloc_x_DG_P**2 + veloc_z_DG_P**2 )
               endif
               
               ! Notify that we already went through this node
               !not_treated_coupling = .false.
               
               if(.false.) WRITE(*,*) elastic_tensor, coord
               !if(elastic_tensor(i_el,j_el,ispec_el,1) > 0.000001) &
               !WRITE(*,*) "ACOUSTIC >>",i_el,j_el,ispec_el,elastic_tensor(i_el,j_el,ispec_el,:)
               !if(elastic_tensor(i_el,j_el,ispec_el,4) > 0) &
               !if(veloc_z_DG_P > 0.00001) WRITE(*,*) i,j,ispec,E_DG_P, rhovz_DG_P, veloc_z_DG_P, &
                !elastic_tensor(i_el,j_el,ispec_el,:), nx, nz
               !WRITE(*,*) "elastic_tensor", nx, nz, p_DG_P, veloc_x, veloc_z, elastic_tensor(i_el,j_el,ispec_el,4)
        ! Classical boundary conditions
        else
         ! WRITE(*,*) "ispec_is_acoustic_coupling_el2", i,j,ispec, neighbor(:), &
         !       nx, nz, weight  
                
                call boundary_condition_DG(i, j, ispec, timelocal, rho_DG_P, rhovx_DG_P, rhovz_DG_P, E_DG_P, &
                veloc_x_DG_P, veloc_z_DG_P, p_DG_P)
                
                if(ONLY_PERTURBATION) then
                        rhovx_DG_P = rhovx_DG_P+rho0_DG(iglobP)*veloc_x
                        rhovz_DG_P = rhovz_DG_P+rho0_DG(iglobP)*veloc_z
                        
               !         WRITE(*,*) "FRONTEER", i,j,ispec,  rho_DG_P, rhovx_DG_P, rhovz_DG_P, E_DG_P, &
               ! veloc_x_DG_P, veloc_z_DG_P, p_DG_P, coord(:,ibool(i,j,ispec))
                endif
                
                rho_DG_P     = rho_DG(iglobM)
                veloc_x_DG_P = veloc_x_DG(iglobM)
                veloc_z_DG_P = veloc_z_DG(iglobM)
                E_DG_P       = E_DG(iglobM)
                p_DG_P       = p_DG(iglobM)
                rhovx_DG_P   = rho_DG_P*veloc_x_DG_P
                rhovz_DG_P   = rho_DG_P*veloc_z_DG_P
                
                !WRITE(*,*) "TEST", p_DG_P
                !if(veloc_z_DG_P > 0) WRITE(*,*) veloc_z_DG_P 
                
                !if(coord(1,ibool(i,j,ispec)) == 4800. &
                !        .AND. coord(2,ibool(i,j,ispec)) == 4800.) then
                !WRITE(*,*) "VERIF BC", ispec, i, j, coord(:, ibool(i,j,ispec))
                !WRITE(*,*) "coord : ", nx, nz, weight
                !WRITE(*,*) "iglobM : ", iglobM 
                !WRITE(*,*) "ispec : ", i,j,ispec 
                !WRITE(*,*) "*********"
                !endif
        endif
        ! Not an outside edge
        else
                iglobP       = ibool_DG(neighbor(1),neighbor(2),neighbor(3))
                !WRITE(*,*) "coord : ", coord(:,ibool_before_perio(i,j,ispec)), nx, nz, weight, i, j, ispec
                !WRITE(*,*) "iglobP : ", iglobP, coord(:,ibool_before_perio(neighbor(1),neighbor(2),neighbor(3))), &
                !        neighbor(1),neighbor(2),neighbor(3), &
                !        ibool(neighbor(1),neighbor(2),neighbor(3)), ibool(i,j,ispec)
                !WRITE(*,*) "*********"
                rho_DG_P     = rho_DG(iglobP)
                veloc_x_DG_P = veloc_x_DG(iglobP)
                veloc_z_DG_P = veloc_z_DG(iglobP)
                E_DG_P       = E_DG(iglobP)
                p_DG_P       = p_DG(iglobP)
                rhovx_DG_P   = rhovx_DG(iglobP)
                rhovz_DG_P   = rhovz_DG(iglobP)
                
                !if(ONLY_PERTURBATION) then
                !        rho_DG_P     = rho_DG(iglobP) + rho0_DG(iglobM)
                !        veloc_x_DG_P = veloc_x_DG(iglobP)
                !        veloc_z_DG_P = veloc_z_DG(iglobP)
                !        E_DG_P       = E_DG(iglobP) + E0_DG(iglobM)
                !        p_DG_P       = p_DG(iglobP) + p0_DG(iglobM)
                !        rhovx_DG_P   = rhovx_DG(iglobP)
                !        rhovz_DG_P   = rhovz_DG(iglobP)
                !endif
        endif
        
        if(ONLY_PERTURBATION) then
                        rho_DG_P   = rho0_DG(iglobM) + rho_DG_P
                        E_DG_P     = E0_DG(iglobM) + E_DG_P
                        ! - p_DG_0(ibool_DG(i,j,ispec))
                        !WRITE(*,*) "rho_DG_P", rho_DG_P, E_DG_P
        endif
        
        !WRITE(*,*) "*********"
        
        !if(ispec_is_acoustic_coupling_el(i,j,ispec,3) >= 0 .AND. p_DG_P > 0) WRITE(*,*) "TOTO >> ", p_DG_P, p_DG(iglobM)
        
        ! Approximate local maximum linearized acoustic wave speed
        ! Lax-Friedrich
        lambda = max( sqrt(veloc_x_DG(iglobM)**2+veloc_z_DG(iglobM)**2) + &
                sqrt(abs(gamma_euler*p_DG(iglobM)/rho_DG(iglobM))), &
                sqrt(veloc_x_DG_P**2+veloc_z_DG_P**2) + &
                sqrt(abs(gamma_euler*p_DG_P/rho_DG_P)) )
        ! Roe
        !lambda = max( abs(veloc_x_DG(iglobM)*nx+veloc_z_DG(iglobM)*nz) + &
        !        sqrt(abs(gamma_euler*p_DG(iglobM)/rho_DG(iglobM))), &
        !        abs(veloc_x_DG_P*nx + veloc_z_DG_P*nz) + &
        !        sqrt(abs(gamma_euler*p_DG_P/rho_DG_P)) )
        
        !if ((ispec_is_acoustic_surface(i,j,ispec) .AND. it_corner < 2) .OR. &
        !        (ispec_is_acoustic_surface_corner(i,j,ispec) .AND. it_corner == 2)) &
        !                lambda = sqrt(veloc_x_DG(iglobM)**2+veloc_z_DG(iglobM)**2) + &
        !        sqrt(abs(gamma_euler*p_DG(iglobM)/rho_DG(iglobM)))
        !lambda = max(lambda, 1._CUSTOM_REAL)
      
        !WRITE(*,*) "CHECK NORMALS", i, j, ispec, nx, nz
      
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Mass conservation equation
        temp_unknown_M  = rhovx_DG(iglobM)
        temp_unknown_P  = rhovx_DG_P
        
        temp_unknown2_M = rhovz_DG(iglobM)
        temp_unknown2_P = rhovz_DG_P
        
        ! compute dot product
        flux_x = temp_unknown_M + temp_unknown_P
        flux_z = temp_unknown2_M + temp_unknown2_P
        flux_n = flux_x*nx + flux_z*nz
        jump   = rho_DG(iglobM) - rho_DG_P
        !if(abs(jump) > 0) WRITE(*,*) it,timelocal,"JUMP > 0", jump, i, j, ispec, rho_DG(iglobM), rho_DG_P
        !if(abs(flux_x) > 0) WRITE(*,*) "flux_x > 0", flux_x
        !if(abs(jump) < threshold) jump = ZERO
        !if(ispec_is_acoustic_coupling_el(i,j,ispec,3) < 0) then
        !if ( (ispec_is_acoustic_surface(i,j,ispec) .AND. it_corner < 2) .OR. &
        !        (ispec_is_acoustic_surface_corner(i,j,ispec) .AND. it_corner == 2) ) then
        !        dot_rho(iglobM) = dot_rho(iglobM) - weight*(temp_unknown_P*nx + temp_unknown2_P*nz)
        !else
                dot_rho(iglobM) = dot_rho(iglobM) - weight*(flux_n + lambda*jump)*HALF
        !endif
        !else        
        !        dot_rho(iglobM) = dot_rho(iglobM) - 0.
        !endif
        
        
        !!!!!!!!!!!!!!!!!!!!!
        ! x-Momentum equation
        temp_unknown_M = rho_DG(iglobM)*veloc_x_DG(iglobM)**2 + p_DG(iglobM)
        !temp_unknown_M = temp_unknown_M - rho_DG(iglobM)*potential_dphi_dx_DG(ibool(i,j,ispec))*p_DG_P
        temp_unknown_P = rho_DG_P*veloc_x_DG_P**2 + p_DG_P
        !temp_unknown_P = temp_unknown_P - rho_DG_P*potential_dphi_dx_DG(ibool(i,j,ispec))*p_DG_P
        
        temp_unknown2_M = rho_DG(iglobM)*veloc_x_DG(iglobM)*veloc_z_DG(iglobM)
        temp_unknown2_P = rho_DG_P*veloc_x_DG_P*veloc_z_DG_P
        
        ! compute dot product
        flux_x = temp_unknown_M + temp_unknown_P
        flux_z = temp_unknown2_M + temp_unknown2_P
        flux_n = flux_x*nx + flux_z*nz
        jump   = rhovx_DG(iglobM) - rhovx_DG_P
        !if(abs(jump) < threshold) jump = ZERO
        !if(ispec_is_acoustic_coupling_el(i,j,ispec,3) < 0) then
        !        if ( (ispec_is_acoustic_surface(i,j,ispec) .AND. it_corner < 2) .OR. &
        !        (ispec_is_acoustic_surface_corner(i,j,ispec) .AND. it_corner == 2) ) then
        !        dot_rhovx(iglobM) = dot_rhovx(iglobM) - weight*(temp_unknown_P*nx + temp_unknown2_P*nz)
        !else
                dot_rhovx(iglobM) = dot_rhovx(iglobM) - weight*(flux_n + lambda*jump)*HALF
        !endif
        !else    
        !        dot_rhovx(iglobM) = dot_rhovx(iglobM) &
        !                - weight*( rho_DG_P*veloc_x_DG_P**2*nx + rho_DG_P*veloc_x_DG_P*veloc_z_DG_P*nz )
        !endif
        !WRITE(*,*) timelocal,"jump dot_rhovx", jump, flux_n, weight, temp_unknown_M, temp_unknown_P, nx, nz, i, j, ispec
        !!!!!!!!!!!!!!!!!!!!!
        ! z-Momentum equation
        temp_unknown_M  = rho_DG(iglobM)*veloc_x_DG(iglobM)*veloc_z_DG(iglobM)
        temp_unknown_P  = rho_DG_P*veloc_x_DG_P*veloc_z_DG_P
        
        temp_unknown2_M = rho_DG(iglobM)*veloc_z_DG(iglobM)**2 + p_DG(iglobM)
        !temp_unknown2_M = temp_unknown2_M - rho_DG(iglobM)*potential_dphi_dx_DG(ibool(i,j,ispec))*p_DG(iglobM)
        temp_unknown2_P = rho_DG_P*veloc_z_DG_P**2 + p_DG_P
        !temp_unknown2_P = temp_unknown2_P - rho_DG_P*potential_dphi_dx_DG(ibool(i,j,ispec))*p_DG_P
          
        ! compute dot product
        flux_x = temp_unknown_M + temp_unknown_P
        flux_z = temp_unknown2_M + temp_unknown2_P
        flux_n = flux_x*nx + flux_z*nz
        jump   = rhovz_DG(iglobM) - rhovz_DG_P
        !if(abs(jump) < threshold) jump = ZERO
        !if(abs(flux_z) > 0) WRITE(*,*) "flux_z > 0", i, j, ispec ,rmass_inverse_acoustic_DG(iglobM),&
        !weight*(flux_n + lambda*jump)*HALF, dot_rhovz(iglobM), &
        !dot_rhovz(iglobM) - weight*(flux_n + lambda*jump)*HALF, &
        !(dot_rhovz(iglobM) - weight*(flux_n + lambda*jump)*HALF)*rmass_inverse_acoustic_DG(iglobM)
        !WRITE(*,*) "jump dot_rhovz", jump, flux_n, E_DG(iglobM), E_DG_P, E_DG(iglobM) - E_DG_P
        !if(ispec_is_acoustic_coupling_el(i,j,ispec,3) < 0) then
        !        if ( (ispec_is_acoustic_surface(i,j,ispec) .AND. it_corner < 2) .OR. &
        !        (ispec_is_acoustic_surface_corner(i,j,ispec) .AND. it_corner == 2) ) then
        !        dot_rhovz(iglobM) = dot_rhovz(iglobM) - weight*(temp_unknown_P*nx + temp_unknown2_P*nz)
        ! else 
                dot_rhovz(iglobM) = dot_rhovz(iglobM) - weight*(flux_n + lambda*jump)*HALF
       !endif
        !else
        !        dot_rhovz(iglobM) = dot_rhovz(iglobM) &
        !                - weight*( rho_DG_P*veloc_z_DG_P**2*nz + rho_DG_P*veloc_x_DG_P*veloc_z_DG_P*nx )
        !endif
         
        !!!!!!!!!!!!!!!!!
        ! Energy equation
        temp_unknown_M = veloc_x_DG(iglobM)*(E_DG(iglobM) + p_DG(iglobM))
        temp_unknown_P = veloc_x_DG_P*(E_DG_P + p_DG_P)
        
        temp_unknown2_M = veloc_z_DG(iglobM)*(E_DG(iglobM) + p_DG(iglobM))
        temp_unknown2_P = veloc_z_DG_P*(E_DG_P + p_DG_P)
        
        if(ONLY_PERTURBATION) then
                temp_unknown_M = temp_unknown_M + veloc_x_DG(iglobM)*p0_DG(iglobM)
                temp_unknown_P = temp_unknown_P + veloc_x_DG_P*p0_DG(iglobM)
                
                temp_unknown2_M = temp_unknown2_M + veloc_z_DG(iglobM)*p0_DG(iglobM)
                temp_unknown2_P = temp_unknown2_P + veloc_z_DG_P*p0_DG(iglobM)
        endif
        ! Test with potential temperature instead E_DG = rho_DG*Theta_DG
        !temp_unknown_M = veloc_x_DG(iglobM)*E_DG(iglobM)
        !temp_unknown_P = veloc_x_DG_P*E_DG_P
        !
        !temp_unknown2_M = veloc_z_DG(iglobM)*E_DG(iglobM)
        !temp_unknown2_P = veloc_z_DG_P*E_DG_P
        
        ! compute dot product
        flux_x = temp_unknown_M + temp_unknown_P
        flux_z = temp_unknown2_M + temp_unknown2_P
        flux_n = flux_x*nx + flux_z*nz
        jump   = E_DG(iglobM) - E_DG_P
        !if(abs(jump) < threshold) jump = ZERO
        !if(ispec_is_acoustic_coupling_el(i,j,ispec,3) >= 0 .AND. abs(flux_n) > 1d-4) &
        !WRITE(*,*) "jump E", jump, flux_n, E_DG(iglobM), E_DG_P, E_DG(iglobM) - E_DG_P
        !if(ispec_is_acoustic_coupling_el(i,j,ispec,3) < 0) then
        !        if ( (ispec_is_acoustic_surface(i,j,ispec) .AND. it_corner < 2) .OR. &
        !        (ispec_is_acoustic_surface_corner(i,j,ispec) .AND. it_corner == 2) ) then
        !        dot_E(iglobM) = dot_E(iglobM) - weight*(temp_unknown_P*nx + temp_unknown2_P*nz)
        !        else
                dot_E(iglobM) = dot_E(iglobM) - weight*(flux_n + lambda*jump)*HALF
        !        endif
       ! else 
        !        dot_E(iglobM) = dot_E(iglobM) - 0.
        !endif
        ! Increment NGLLZ counter
        j = j + 1
        
        ! If at corner and first step => go again 
        if(it_corner == 1) j = j - 1
        ! Reset corner notification
        if(it_corner == 2) it_corner = 0

      enddo
    enddo
    
    !WRITE(*,*) "FIN ELMT", ispec
    
    endif ! end of test if acoustic element
    
  enddo
  
  end subroutine compute_forces_acoustic_DG
  
  ! Compute initial condition
  subroutine initial_condition_DG()

! compute forces in the acoustic elements in forward simulation and in adjoint simulation in adjoint inversion

  use constants,only: CUSTOM_REAL,NGLLX,NGLLZ,gamma_euler,PI

  use specfem_par, only: nglob_DG, ibool_DG, &
        rho_DG, rhovx_DG, rhovz_DG, E_DG, nspec, &
        potential_dphi_dx_DG, potential_dphi_dz_DG

  implicit none
  
  real(kind=CUSTOM_REAL), dimension(nglob_DG) :: veloc_x_DG, veloc_z_DG, p_DG
  
  ! Parameters
  real(kind=CUSTOM_REAL), parameter :: ZEROl = 0._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: ONEl  = 1._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: TWOl  = 2._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: SIXl  = 6._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: HALFl = 0.5_CUSTOM_REAL
  
  integer :: ispec, iglob, i, j
  
  real(kind=CUSTOM_REAL) :: rho_DG_P, rhovx_DG_P, rhovz_DG_P, E_DG_P, &
      veloc_x_DG_P, veloc_z_DG_P, p_DG_P
  
  do ispec = 1,nspec

        ! first double loop over GLL points to compute and store gradients
        do j = 1,NGLLZ
          do i = 1,NGLLX
  
           iglob = ibool_DG(i,j,ispec)
  
           call boundary_condition_DG(i, j, ispec, ZEROl, rho_DG_P, rhovx_DG_P, rhovz_DG_P, E_DG_P, &
                veloc_x_DG_P, veloc_z_DG_P, p_DG_P)
  
           rho_DG(iglob)     = rho_DG_P
           veloc_x_DG(iglob) = veloc_x_DG_P
           veloc_z_DG(iglob) = veloc_z_DG_P
           p_DG(iglob)       = p_DG_P
           rhovx_DG(iglob)   = rhovx_DG_P
           rhovz_DG(iglob)   = rhovz_DG_P
           E_DG(iglob)       = E_DG_P
           
          enddo
       enddo
       
       call setup_gravity_potential(potential_dphi_dx_DG, potential_dphi_dz_DG, ispec)
       
  enddo  
  
  end subroutine initial_condition_DG
  
  ! Compute gravity potential derivatives
  subroutine setup_gravity_potential(potential_dphi_dx_DG,potential_dphi_dz_DG,ispec)

! compute forces in the acoustic elements in forward simulation and in adjoint simulation in adjoint inversion

  use constants,only: CUSTOM_REAL, NGLLX, NGLLZ

  use specfem_par, only: xix, xiz, gammax, gammaz, coord, ibool_before_perio, &
        jacobian, hprime_xx, hprime_zz, ibool, nglob

  implicit none
  
  ! Parameters
  real(kind=CUSTOM_REAL), parameter :: ZEROl = 0._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: ONEl  = 1._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: TWOl  = 2._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: SIXl  = 6._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: HALFl = 0.5_CUSTOM_REAL
  
  integer :: ispec, i, j, k
  
  real(kind=CUSTOM_REAL) :: dphi_dxi, dphi_dgamma, dphi_dxl, dphi_dzl, &
        xixl, xizl, gammaxl, gammazl, &
        jacobianl
        
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: potential_phi_DG
        
  real(kind=CUSTOM_REAL), dimension(nglob) :: potential_dphi_dx_DG, potential_dphi_dz_DG
  
    do i=1,NGLLX
        do j=1,NGLLZ
        potential_phi_DG(i,j) = 1.*( real(coord(2,ibool_before_perio(i,j,ispec)),kind=CUSTOM_REAL) &
              +  real(coord(1,ibool_before_perio(i,j,ispec)),kind=CUSTOM_REAL) )
        enddo
    enddo
  
    do i=1,NGLLX
        do j=1,NGLLZ
  
          ! Gravity potential derivatives
          jacobianl = jacobian(i,j,ispec)
        
          xixl = xix(i,j,ispec)
          xizl = xiz(i,j,ispec)
          gammaxl = gammax(i,j,ispec)
          gammazl = gammaz(i,j,ispec)
           
          dphi_dxi    = ZEROl
          dphi_dgamma = ZEROl
          
          ! Type of potential
          !potential_phi_DG = coord(2,ibool_before_perio(i,j,ispec))
          !potential_phi_DG = coord(1,ibool_before_perio(i,j,ispec)) &
          !      + coord(2,ibool_before_perio(i,j,ispec))
          
          ! first double loop over GLL points to compute and store gradients
          ! we can merge the two loops because NGLLX == NGLLZ
          do k = 1,NGLLX
            dphi_dxi    = dphi_dxi    + potential_phi_DG(k,j) * hprime_xx(i,k)
            dphi_dgamma = dphi_dgamma + potential_phi_DG(i,k) * hprime_zz(j,k)
            ! Test w/ Well-balanced condition (Based on hydrisattic solution)
            !dphi_dxi    = dphi_dxi    + exp(-potential_phi_DG(k,j)/(1.)) * hprime_xx(i,k)
            !dphi_dgamma = dphi_dgamma + exp(-potential_phi_DG(i,k)/(1.)) * hprime_zz(j,k)
          enddo

          ! derivatives of potential
          dphi_dxl = dphi_dxi * xixl + dphi_dgamma * gammaxl
          dphi_dzl = dphi_dxi * xizl + dphi_dgamma * gammazl 
          potential_dphi_dx_DG(ibool(i,j,ispec)) = dphi_dxl
          !potential_dphi_dx_DG(ibool(i,j,ispec)) = -exp(potential_phi_DG(i,j))*dphi_dxl
          !if(abs( potential_dphi_dx_DG(ibool(i,j,ispec))) < 0.0000000001) potential_dphi_dx_DG(ibool(i,j,ispec)) = 0.
          potential_dphi_dz_DG(ibool(i,j,ispec)) = dphi_dzl
          !potential_dphi_dz_DG(ibool(i,j,ispec)) = -exp(potential_phi_DG(i,j))*dphi_dzl
          !WRITE(*,*) "dphi_dzl", dphi_dzl, i, j, ispec,  coord(2,ibool_before_perio(i,j,ispec))
          
     enddo
  enddo   
          
  end subroutine setup_gravity_potential
  
  ! Compute initial condition
  subroutine boundary_condition_DG(i, j, ispec, timelocal, rho_DG_P, rhovx_DG_P, rhovz_DG_P, E_DG_P, &
      veloc_x_DG_P, veloc_z_DG_P, p_DG_P)

! compute forces in the acoustic elements in forward simulation and in adjoint simulation in adjoint inversion

  use constants,only: CUSTOM_REAL,gamma_euler,PI

  use specfem_par, only: ibool_before_perio, coord

  implicit none
  
  integer :: i, j, ispec
  
  real(kind=CUSTOM_REAL) :: rho_DG_P, rhovx_DG_P, rhovz_DG_P, E_DG_P, &
      veloc_x_DG_P, veloc_z_DG_P, p_DG_P, timelocal
  
  ! Parameters
  real(kind=CUSTOM_REAL), parameter :: ZEROl = 0._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: ONEl  = 1._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: TWOl  = 2._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: SIXl  = 6._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: HALFl = 0.5_CUSTOM_REAL
  
  real(kind=CUSTOM_REAL) :: x, x0, z, z0, r, beta
  real(kind=CUSTOM_REAL) :: to, perio, lambdo, xo
  ! Hydrostatic solution
  real(kind=CUSTOM_REAL) :: RR, p0, rho0
  ! Density current
  !real(kind=CUSTOM_REAL) :: cp, cnu, exner, RR, p0, rho0, rs, theta, theta0
  ! Linear mountains
  !real(kind=CUSTOM_REAL) :: cp, cnu, exner, RR, p0, rho0, rs, theta, theta0, Nsq
  !real(kind=CUSTOM_REAL) :: Tl, p0,  RR
  !real(kind=CUSTOM_REAL) :: Tl, Tu, p0, rho0, RR
  !real(kind=CUSTOM_REAL) :: Tl, Tu, p0, rho0, RR
  !real(kind=CUSTOM_REAL) :: Tl, Tu, p0, rho0, RR, rs!, theta, thetaprime, pibar
  !real(kind=CUSTOM_REAL) :: Tl, Tu, p0, RR, rs, theta0, Nsq, rho0
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! ISENTROPIC VORTEX
  x0 = 5.
  z0 = 5.
  beta = 5.
        
        x = real(coord(1,ibool_before_perio(i,j,ispec)), kind=CUSTOM_REAL)
        z = real(coord(2,ibool_before_perio(i,j,ispec)), kind=CUSTOM_REAL)
        r = sqrt( (x - timelocal - x0)**2 + (z - z0)**2 )
  
        rho_DG_P     = &
                ( ONEl - ((gamma_euler-1)/(16.*gamma_euler*PI**2))*(beta**2)*&
                exp( TWOl*( ONEl - r**2 ) ) )**(1/(gamma_euler - 1))
        
        veloc_x_DG_P = ONEL - beta*exp(ONEl - r**2)*( z - z0 )/(2*PI)
        veloc_z_DG_P = beta*exp(ONEl - r**2)*( x - x0 )/(2*PI)
        
        p_DG_P       = rho_DG_P**gamma_euler
   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! BASIC FORCING
   rho_DG_P     = ONEl*exp(-z/10000)
   veloc_x_DG_P = ZEROl
   veloc_z_DG_P = ZEROl
   if(z == ZEROl) &
   veloc_z_DG_P = -0.2*((timelocal-0.1)/0.1)*(1/0.1)*exp(-((timelocal-0.1)/0.1)**2)*&
        sin(x*(2*PI)/(1000))
        
   ! Acoustic plane wave
   to = 10.
   perio = 10.     
   if(z == ZEROl) &     
   veloc_z_DG_P = (&
                  - (2d0/(perio/4d0))*((timelocal-(to-perio/4d0))/(perio/4d0))* &
                           (exp(-((timelocal-(to-perio/4d0))/(perio/4d0))**2)) &
                  + (2d0/(perio/4d0))*((timelocal-(to+perio/4d0))/(perio/4d0))* &
                           (exp(-((timelocal-(to+perio/4d0))/(perio/4d0))**2)) ) 
   lambdo = 500
   xo = 1000
   ! Simple gravity wave
   if(z == ZEROl) &     
   veloc_z_DG_P = (&
                  - (2d0/(perio/4d0))*((timelocal-(to-perio/4d0))/(perio/4d0))* &
                           (exp(-((timelocal-(to-perio/4d0))/(perio/4d0))**2)) &
                  + (2d0/(perio/4d0))*((timelocal-(to+perio/4d0))/(perio/4d0))* &
                           (exp(-((timelocal-(to+perio/4d0))/(perio/4d0))**2)) ) &
                   * ( exp(-((x-(xo-lambdo/4))/(lambdo/4))**2) - &
                          exp(-((x-(xo+lambdo/4))/(lambdo/4))**2) ) 
   
   !WRITE(*,*) x,z,"veloc_z_DG_P",veloc_z_DG_P
                           
   p_DG_P = -rho_DG_P*9.81!(340**2)*rho_DG_P/gamma_euler 
   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Rayleigh-Taylor instability (case 5.3 paper)
   !Tl = TWOl
   !Tu = ONEl
   !p0   = ONEl
   !rho0 = ONEl
   !RR   = 1.!8.3145_CUSTOM_REAL
   !!if(z <= 1._CUSTOM_REAL) &
   !p_DG_P = p0*exp(-(z-ONEl/(TWOl))/(RR*Tl))
   !if(z > ONEl/TWOl) &
  ! p_DG_P = p0*exp(-(z-ONEl/(TWOl))/(RR*Tu))
  ! !if(z <= 1._CUSTOM_REAL) &
  ! rho_DG_P = p_DG_p/(RR*Tl)
  ! if(z > ONEl/TWOl) &
  ! rho_DG_P = p_DG_p/(RR*Tu)
  ! veloc_x_DG_P = ZEROl
  ! veloc_z_DG_P = ZEROl
   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Steady shear flow
   !rho_DG_P = ONEl
   !RR   = 1.!8.3145_CUSTOM_REAL
  ! 
  ! veloc_x_DG_P = ZEROl
  ! veloc_z_DG_P = (z/10)**2
  ! 
  ! E_DG_P = (10.+2.*0.01*x)/(gamma_euler-1.) + ((z/10)**4)/2.
   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Pressure pulse
   !rho_DG_P = ONEl
  ! 
  ! veloc_x_DG_P = ZEROl
  ! veloc_z_DG_P = ZEROl
  ! 
  ! E_DG_P = (12.)/(gamma_euler-1.) + &
  !      (rho_DG_P/2.)*exp( -4.*( cos(PI*x)**2 + cos(PI*z)**2 ) )
   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! RTI 2nd test
   !p0   = ONEl
   !rho0 = ONEl
   !RR = 6.
   !! epsilonx
   !Tl = 0.1_CUSTOM_REAL*sqrt( 1.4_CUSTOM_REAL/ TWOl )
   !! epsilony
   !Tu = -Tl*RR/16.
  !
  ! p_DG_P = TWOl - z
  ! if(z >= ONEl/TWOl) p_DG_P = TWOl - TWOL*z
  !  
  ! rho_DG_P  = ONEl
  ! if(z >= ONEl/TWOl) rho_DG_P = TWOl
  ! 
  ! veloc_x_DG_P = Tu*sin(8.*PI*x)*cos(PI*z)*(sin(PI*z))**(RR-1)
  ! 
  ! veloc_z_DG_P = -Tl*cos(8.*PI*x)*(sin(PI*z))**(RR)
   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! HYDROSTATIC SOLUTION 
   RR = 0.0
   rho0 = 1.
   p0   = 1.
   rho_DG_P = rho0*exp(-(x + z)*rho0/p0)
   p_DG_P   = p0*exp(-(x + z)*rho0/p0) &
        + RR*exp(-100._CUSTOM_REAL*(rho0/p0)*( (x - 0.3_CUSTOM_REAL)**2 &
        + (z - 0.3_CUSTOM_REAL)**2 ))
   veloc_x_DG_P = ZEROl
   veloc_z_DG_P = ZEROl
   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Thermal bubble
   ! Location
   !x0 = 5.
   !z0 = 2.5
   !r = sqrt( (x - x0)**2 + (z - z0)**2 )
   !rs = 2.5
   !! initial cond.
   !Tl = 2.
   !Tu = 1.
   !p0   = ONEl/1.
   !rho0 = ONEl
   !RR   = 2.!8.3145_CUSTOM_REAL
   !
   !theta = 300
   !thetaprime = 0
   !if(r <= rs) thetaprime = 150.25*(1 + cos(PI*r/rs))
   !p_DG_P = p0
   !pibar = (-9.81/(7*theta))*z
   !E_DG_P = rho0*5.*(theta + thetaprime)*pibar
   !p_DG_P = p0*exp(-(r-rs)/(RR*Tl))
   !if(r <= rs) &
   !p_DG_P = p0*exp(-(r-rs)/(RR*Tu))
   !if(abs(p_DG_P) > 10000000) WRITE(*,*) "p_DG_P", p_DG_P
   !:WRITE(*,*) "TOITOTTO",p_DG_P,x,x0,z,z0,r,p0,RR,Tl,Tu,exp(-(r-250)/(RR*Tu))
   !if(z <= 1._CUSTOM_REAL) &
   !rho_DG_P = p_DG_p/(RR*Tl)
   !if(r > rs) &
   !rho_DG_P = p_DG_p/(RR*Tu)
   !rho_DG_P = rho0
   !veloc_x_DG_P = ZEROl
   !veloc_z_DG_P = ZEROl
   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! TEST BOUNDARY     
   !Tl = 1.
   !p0   = ONEl/1.
   !RR   = 2.!8.3145_CUSTOM_REAL 
   !p_DG_P   = p0*exp(-z/(RR*Tl))    
   !rho_DG_P = p_DG_P/(RR*Tl)
   !veloc_x_DG_P = ZEROl
   !veloc_z_DG_P = ZEROl
   !WRITE(*,*) "p_DG_P",p_DG_P,rho_DG_P
   
   !p_DG_P   = 1.
   !rho_DG_P = 1.
   !veloc_x_DG_P = ZEROl
   !veloc_z_DG_P = ZEROl
   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! PERTURBATION INITIALIZATION
   !if(ONLY_PERTURBATION) then
   !     rho_DG_P     = ZEROl
   !     veloc_x_DG_P = ZEROl
   !     veloc_z_DG_P = ZEROl
   !     p_DG_P       = ZEROl
   !     if(z == ZEROl) &
   !        p_DG_P = -0.2*((timelocal-0.2)/0.1)*(1/0.1)*exp(-((timelocal-0.2)/0.1)**2)!*&
                !sin(x*(2*PI)/(0.25))
   !endif
        
        E_DG_P       = p_DG_P/(gamma_euler - ONEl) &
                        + rho_DG_P*HALFl*( veloc_x_DG_P**2 + veloc_z_DG_P**2 ) 
                        
        ! Test with potential temperature instead E_DG = Theta_DG
        ! => We know E_DG = Theta_DG and search for P_DG
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! RISING THERMAL BUBBLE WITH POTENTIAL TEMPERATURE
        !x0 = 5.
        !z0 = 2.5
        !r = sqrt( (x - x0)**2 + (z - z0)**2 )
        !rs = 2.5
        !! initial cond.
        !p0   = ONEl
        !RR   = 2.!8.3145_CUSTOM_REAL
        !!
        !E_DG_P  = 300.
        !if(r <= rs) E_DG_P  = E_DG_P + (273.65/2.)*(1 + cos(PI*r/rs))
        !!if(r <= rs) E_DG_P  = E_DG_P + 100.5*exp(-(r/1)**2)
        !rho_DG_P = (p0/(RR*E_DG_P))*(9.81*z/(7*300.))**(5/RR)+1.
        !
        !E_DG_P = E_DG_P*rho_DG_P
        !
        !p_DG_P       = p0 * (RR*E_DG_P/p0)**gamma_euler
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! FORCING WITH POTENTIAL TEMPERATURE
        !E_DG_P  = 0.
       ! 
       ! rho_DG_P = 0.
       ! 
       ! p_DG_P    = 0.
       ! 
       ! veloc_x_DG_P = ZEROl
       ! veloc_z_DG_P = ZEROl
       ! if(z == ZEROl) &
       ! veloc_z_DG_P = -0.2*((timelocal-0.02)/0.1)*(1/0.1)*exp(-((timelocal-0.02)/0.1)**2)
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! INERTIA GRAVITY WAVES
        !x0 = 33333
        !z0 = 2.5
        !r = sqrt( (x - x0)**2 + (z - z0)**2 )
        !rs = 2.5
        !! initial cond.
        !p0   = 100000
        !rho0 = ONEl
        !
        !cp  = 7.
        !cnu = 5.
        !RR   = cp - cnu
        !
        !Nsq  = 0.01**2
       ! 
       ! theta0 = 300.
       ! theta = theta0*exp(Nsq*z/9.81)
       ! 
       ! veloc_x_DG_P = 20.
       ! veloc_z_DG_P = ZEROl
       ! theta  = theta + (273.15+0.01)*sin(PI*z/10000)/(1 + ((x-x0)/1600)**2)
       ! exner = 1 + ((9.81**2)/(cp*theta0*Nsq))*(exp( -Nsq*z/(9.81) ) - 0.2)
       ! rho_DG_P = (p0/(RR*theta))*( exner )**(5/RR)
       ! 
       ! 
       ! E_DG_P  = rho_DG_P*cnu*theta + (rho_DG_P/2.)*abs(veloc_x_DG_P)
        
        !WRITE(*,*) x,z,"exner", theta, exner, E_DG_P, ((9.81**2)/(cp*theta0*Nsq))*(exp( -Nsq*z/(9.81) ) - 1.)
        !E_DG_P  = theta0*exp( -Nsq*z/(9.81) ) &
        !        + (273.15+0.01)*( sin(PI*z/10000)/(1 + ((x-x0)/5000)**2) )
        !if(r <= rs) E_DG_P  = E_DG_P + 100.5*exp(-(r/1)**2)
        !WRITE(*,*) "xz", x, z
        !WRITE(*,*) "E_DG", x, z, E_DG_P, (p0/(RR*E_DG_P)),  0 + &
        !        ((9.81**2)/(7.*theta0*Nsq)),(exp( -Nsq*z/(9.81) ) - 1.) , &
        !        1 + ((9.81**2)/(7.*theta0*Nsq))*(exp( -Nsq*z/(9.81) ) - 1.)
        !(p0/(RR*E_DG_P)), ( 1 + &
        !        ((9.81**2)/(7.*theta0*Nsq))*(E_DG_P/theta0 - 1.) )**(5/RR)
        
        
        !E_DG_P = E_DG_P*rho_DG_P
        
   !     p_DG_P       = p0 * (RR*E_DG_P/p0)**gamma_euler
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! DENSITY CURRENT
        !x0 = 0.
        !z0 = 3000.
        !r = sqrt( ((x - x0)/4000.)**2 + ((z - z0)/2000.)**2 )
        !rs = 1
        !! initial cond.
        !p0   = 10000
        !rho0 = ONEl
        !theta0 = 273.15-15.
        !cp  = 7.
        !cnu = 5.
        !RR   = cp - cnu!8.3145_CUSTOM_REAL
        !
        !theta  = 300.
        !if(r <= rs) theta  = theta + (theta0/2.)*(1+cos(PI*r))
        !
        !veloc_x_DG_P = ZEROl
        !veloc_z_DG_P = ZEROl
        !
        !exner = -9.81*z/(7*theta) + 1
        !!WRITE(*,*) "THETA", theta, exner
        !rho_DG_P = (p0/(RR*theta))*exner**(cnu/RR) 
        !E_DG_P   = rho_DG_P*cnu*PI*theta
       ! 
       ! p_DG_P     = (gamma_euler - ONEl)*( E_DG_P &
       ! - (HALFl)*rho_DG_P*( veloc_x_DG_P**2 + veloc_z_DG_P**2 ) )
                  
        rhovx_DG_P   = rho_DG_P*veloc_x_DG_P
        rhovz_DG_P   = rho_DG_P*veloc_z_DG_P
  
  end subroutine boundary_condition_DG
  
  ! SLOPE LIMITER
  
  subroutine SlopeLimit1(Q)
  
    use constants,only: CUSTOM_REAL,NGLLX,NGLLZ

    use specfem_par, only: nspec, nglob_DG, ibool_DG, Vandermonde, invVandermonde, &
        neighbor_DG

    implicit none 
  
    integer :: iglob, iglob2, ispec, j, k, l, h, prod_ind
    integer, parameter :: NGLL = NGLLX*NGLLZ
    real(kind=CUSTOM_REAL) :: vkm1_x,vkp1_x,vkm1_z,vkp1_z
    real(kind=CUSTOM_REAL), dimension(nglob_DG) :: Q, ulimit
    real(kind=CUSTOM_REAL), dimension(NGLL) :: uh, uavg, uhl, ulimit_temp
    real(kind=CUSTOM_REAL), dimension(1:nspec,1:NGLL) :: ul
    real(kind=CUSTOM_REAL), dimension(1:nspec) :: v, vk
    !real(kind=CUSTOM_REAL), dimension(1:NGLLX,1:NGLLZ) :: val_in_elmnt
    
    ! Parameters
    real(kind=CUSTOM_REAL), parameter :: ZEROl = 0._CUSTOM_REAL

    ulimit = ZEROl

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
    uhl(3:NGLLZ) = ZEROl
    uhl(NGLLZ+2:NGLL) = ZEROl
    
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
   !WRITE(*,*) ispec, "cell average v(ispec)", v(ispec), maxval(val_in_elmnt),minval(val_in_elmnt)!, minval(Q(ibool_DG(:,1,ispec))), &
      !  maxval(Q(ibool_DG(:,1,ispec))), minval(Q(ibool_DG(:,5,ispec))), &
      !  maxval(Q(ibool_DG(:,5,ispec)))
   
   enddo ! do ispec=1,NSPEC
   
   ! find cell averages
   vk = v
   
    do ispec=1,nspec
    
   ! apply slope limiter to selected elements
     ! Find neighbors
     vkm1_x = v(ispec)
     if(neighbor_DG(1,NGLLZ/2,ispec,3) > - 1) vkm1_x = v(neighbor_DG(1,NGLLZ/2,ispec,3))
     vkp1_x = v(ispec)
     if(neighbor_DG(NGLLX,NGLLZ/2,ispec,3) > - 1) vkp1_x = v(neighbor_DG(NGLLX,NGLLZ/2,ispec,3))
     vkm1_z = v(ispec)
     if(neighbor_DG(NGLLX/2,1,ispec,3) > - 1) vkm1_z = v(neighbor_DG(NGLLX/2,1,ispec,3))
     vkp1_z = v(ispec)
     if(neighbor_DG(NGLLX/2,NGLLZ,ispec,3) > - 1) vkp1_z = v(neighbor_DG(NGLLX/2,NGLLZ,ispec,3))
     
     !WRITE(*,*) "TEST >> ", neighbor_DG(1,NGLLZ/2,ispec,3), neighbor_DG(NGLLX,NGLLZ/2,ispec,3), &
     !   neighbor_DG(NGLLX/2,1,ispec,3), neighbor_DG(NGLLX/2,NGLLZ,ispec,3)
    ! WRITE(*,*) "2TEST >> ",   vkm1_x, vkp1_x, vkm1_z, vkp1_z
     ! Compute Limited flux
     call SlopeLimitLin(ulimit_temp, ispec, ul(ispec,:),vk(ispec),vkm1_x,vkp1_x,vkm1_z,vkp1_z)
     
     prod_ind = 0
     do k=1,NGLLX
       do l=1,NGLLZ
        iglob = ibool_DG(k,l,ispec)
        prod_ind = prod_ind + 1
        ulimit(iglob) = ulimit_temp(prod_ind)
       enddo
     enddo
   enddo
   ! Update solution
   Q(1:nglob_DG) = ulimit(1:nglob_DG)
   
  end subroutine SlopeLimit1
  
  subroutine SlopeLimitLin(ulimit, ispec, ul, v0, vkm1_x,vkp1_x,vkm1_z,vkp1_z)
  
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
        integer :: ispec, iglob, j, k, m, n, i
        real(kind=CUSTOM_REAL), dimension(1:NGLLX,1:NGLLZ) :: ul_loc
        
        ! Parameters
        real(kind=CUSTOM_REAL), parameter :: ZEROl = 0._CUSTOM_REAL
        real(kind=CUSTOM_REAL), parameter :: ONEl  = 1._CUSTOM_REAL
        real(kind=CUSTOM_REAL), parameter :: TWOl  = 2._CUSTOM_REAL
        real(kind=CUSTOM_REAL), parameter :: HALFl = ONEl/TWOl
        
        real(kind=CUSTOM_REAL), parameter :: gradient_factor = ONEl/TWOl

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
        x0 = minval(xl) + hx/TWOl
        z0 = minval(zl) + hz/TWOl
        
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
                ux(m) = ( du_dxi * xixl + du_dgamma * gammaxl )
                uz(m) = ( du_dxi * xizl + du_dgamma * gammazl )
                
          enddo
        enddo
        
        ! Limit function
        !do k=1,NGLL
        ! m = 0
        ! do n=1,NGLLX
        !  do j=1,NGLLZ
        !        m = m + 1
        !        xixl = xix(n,j,ispec)
        !        xizl = xiz(n,j,ispec)
        !        gammaxl = gammax(n,j,ispec)
        !        gammazl = gammaz(n,j,ispec)
        !  
        !        ux(k)  = ux(k) &
        !                + (TWOl/hNx(k))*(xixl*Drx(k,m) &
        !                + gammaxl*Drz(k,m))*ul(m) 
        !                
        !        uz(k)  = uz(k) &
        !                + (TWOl/hNz(k))*(xizl*Drx(k,m) &
        !                + gammazl*Drz(k,m))*ul(m)
        !  enddo
        ! enddo
        !enddo

        ulimit = v0

        minmod(1,1) = HALFl*(ux(1)+ux(NGLLZ))
        minmod(2,1) = (vkp1_x-v0)/(hx*gradient_factor)
        minmod(3,1) = (v0-vkm1_x)/(hx*gradient_factor)
        
        !WRITE(*,*) ispec,"TEST MINMOD", v0, &
        !        ux(1), (vkp1_x-v0)/hx, (v0-vkm1_x)/hx, uz(1), (vkp1_z-v0)/hz, (v0-vkm1_z)/hz
        
       ! stop
        
        call minmod_computeB(ulimit_temp, minmod, 3, 1,ispec,hx)
        !WRITE(*,*) "ux", ux
        !WRITE(*,*) "ulimit_tempx", ulimit_temp
        
        ulimit = ulimit + (xl-x0)*ulimit_temp(1) 
        
        minmod(1,1) = HALFl*(uz(1)+uz(NGLLX*(NGLLZ-1)+1))
        minmod(2,1) = (vkp1_z-v0)/(hz*gradient_factor)
        minmod(3,1) = (v0-vkm1_z)/(hz*gradient_factor)
        
        call minmod_computeB(ulimit_temp, minmod, 3, 1,ispec,hz)
        !WRITE(*,*) "uz", uz
        !WRITE(*,*) "ulimit_tempz", ulimit_temp
        !WRITE(*,*) "********"
        
        ulimit = ulimit + (zl-z0)*ulimit_temp(1)
        
  end subroutine SlopeLimitLin
  
  subroutine minmod_computeB(R, v, n, m, ispec, h)

     use constants,only: CUSTOM_REAL

     implicit none 

        ! Usually 
        ! m = NSPEC
        ! n = 3
        integer :: n, m, ispec
        real(kind=CUSTOM_REAL), dimension(1:n,1:m) :: v
        real(kind=CUSTOM_REAL), dimension(1:m) :: R
        
        real(kind=CUSTOM_REAL) :: h, M_param
        
        real(kind=CUSTOM_REAL), parameter :: ZEROl = 0._CUSTOM_REAL
        real(kind=CUSTOM_REAL), parameter :: ONEl  = 1._CUSTOM_REAL
        
        ! function mfunc = minmod(v)
        ! Purpose: Implement the midmod function v is a vector

        M_param = 1.

        R = ZEROl
        
        ! Find regions where limiter is needed
        if(abs(v(1,1)) > M_param*h**2) then
                call minmod_compute(R, v, n, m,ispec)
        else
                R = v(1,1)
        endif
        
   end subroutine minmod_computeB
  
  subroutine minmod_compute(R, v, n, m,ispec)

     use constants,only: CUSTOM_REAL

     implicit none 

        ! Usually 
        ! m = NSPEC
        ! n = 3
        integer :: n, m, ispec, k, j
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
             !WRITE(*,*) "=====>", k, n, j, m, v(k,j), s(j)
             enddo
        enddo
        
        !if(~isempty(ids))
        do k=1,m
          if(abs(s(k)) == ONEl) then
                R(k) = s(k)*minval(abs(v(:,k)))
                !WRITE(*,*) "CHOSEN =>", ispec, v
          endif
        enddo
        
        !endif
        
   end subroutine minmod_compute

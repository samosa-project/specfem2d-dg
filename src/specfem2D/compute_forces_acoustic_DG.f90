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
! compute_forces_acoustic_DG                                   !
! ------------------------------------------------------------ !
! TODO: Description.
! compute forces in the acoustic elements in forward simulation and in adjoint simulation in adjoint inversion

  subroutine compute_forces_acoustic_DG(rho_DG_main, rhovx_DG_main, rhovz_DG_main, E_DG_main, &
        T_DG_main, V_DG_main, e1_DG, &
        dot_rho, dot_rhovx, dot_rhovz, dot_E, dot_e1, timelocal)

  use constants,only: CUSTOM_REAL,NGLLX,NGLLZ,gamma_euler

  use specfem_par, only: nglob_DG,nspec, ispec_is_acoustic_DG,&
                         xix,xiz,gammax,gammaz,jacobian, &
                         hprimewgll_xx, &
                         hprimewgll_zz,wxgll,wzgll, &
                         ibool_DG, &
                         it,potential_dphi_dx_DG, potential_dphi_dz_DG, ibool, &
                         DIR_RIGHT, DIR_LEFT, DIR_UP, DIR_DOWN, &
                         myrank, &
                         i_stage, p_DG_init, gammaext_DG, muext, etaext, kappa_DG,tau_epsilon, tau_sigma, &
                         rhovx_init, rhovz_init, E_init, rho_init, &
                         CONSTRAIN_HYDROSTATIC, TYPE_SOURCE_DG, &
                         link_iface_ijispec, nx_iface, nz_iface, weight_iface, neighbor_DG_iface,&
                         ABC_STRETCH &! Stretching-based absorbing conditions.
                         , coord, ibool_before_perio!,cnu
                         
  implicit none

  real(kind=CUSTOM_REAL), dimension(nglob_DG) :: rho_DG_main, rhovx_DG_main, rhovz_DG_main, E_DG_main, e1_DG
  real(kind=CUSTOM_REAL), dimension(2, nglob_DG), intent(in) :: T_DG_main
  real(kind=CUSTOM_REAL), dimension(2, 2, nglob_DG), intent(in) :: V_DG_main
  real(kind=CUSTOM_REAL), dimension(nglob_DG) :: rho_DG, rhovx_DG, rhovz_DG, E_DG
  real(kind=CUSTOM_REAL), dimension(2, nglob_DG) :: T_DG
  real(kind=CUSTOM_REAL), dimension(2, 2, nglob_DG) :: V_DG
  real(kind=CUSTOM_REAL), dimension(nglob_DG), intent(out) :: dot_rho, dot_rhovx, dot_rhovz, dot_E, dot_e1
  
  ! Local variables.
  integer :: ispec,i, j,k, iglob
  integer :: ifirstelem,ilastelem

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: temp_rho_1, temp_rho_2, &
        temp_rhovx_1, temp_rhovx_2, temp_rhovz_1, temp_rhovz_2, &
        temp_E_1, temp_E_2, &
        !temp_rho_gravi, &
        temp_rhovx_gravi, temp_rhovz_gravi, temp_E_gravi!, temp_e1

  ! Jacobian matrix and determinant
  real(kind=CUSTOM_REAL) :: xixl,xizl,gammaxl,gammazl,jacobianl

  real(kind=CUSTOM_REAL) :: lambda, nx, nz, weight, &
        temp_unknown_M, temp_unknown_P, temp_unknown2_M, temp_unknown2_P, &
        temp_unknown, temp_unknown2, &
        flux_n, jump, &
        rho_DG_P, veloc_x_DG_P, veloc_z_DG_P, &
        E_DG_P, p_DG_P, rhovx_DG_P, rhovz_DG_P, timelocal, &
        Tx_DG_P, Tz_DG_P, Vxx_DG_P, Vzz_DG_P, Vxz_DG_P, Vzx_DG_P, T_P, &
        ! TEST
        gamma_P
    
  real(kind=CUSTOM_REAL) :: dT_dx, dT_dz
  !real(kind=CUSTOM_REAL) :: veloc_n_M, veloc_n_P
        
  integer :: iglobM, iglobP
  
  integer, dimension(3) :: neighbor
  
  ! Local
  real(kind=CUSTOM_REAL), dimension(nglob_DG) :: veloc_x_DG, veloc_z_DG, p_DG
  
  ! Viscosity
  real(kind=CUSTOM_REAL) :: dux_dx, dux_dz, duz_dx, duz_dz
  real(kind=CUSTOM_REAL) :: wxl, wzl
  
  ! Parameters
  real(kind=CUSTOM_REAL), parameter :: ZERO = 0._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: ONE  = 1._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: TWO  = 2._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: HALF = 0.5_CUSTOM_REAL
  
  real(kind=CUSTOM_REAL), parameter :: threshold = 0.000000001_CUSTOM_REAL
  
  ! Temporary
  logical, parameter :: ONLY_PERTURBATION = .false.
  
  logical :: exact_interface_flux
  
  ! For better CONSTRAIN_HYDROSTATIC switches.
  integer :: cnst_hdrsttc
  
  ! TEST ABSORB
  real(kind=CUSTOM_REAL) :: maxval_rho,maxval_rhovx,maxval_rhovz,maxval_E
  logical :: ABSORB_BC

  integer :: iface1, iface, iface1_neighbor, iface_neighbor, ispec_neighbor
  
  ! TEST STRETCH
  real(kind=CUSTOM_REAL) :: coef_stretch_x_ij, coef_stretch_x_ij_prime, &
                            coef_stretch_z_ij, coef_stretch_z_ij_prime
  real(kind=CUSTOM_REAL) :: viscous_tens_11, viscous_tens_12, viscous_tens_22
  real(kind=CUSTOM_REAL) :: nx_unit, nz_unit
  
  ! For more convinient CONSTRAIN_HYDROSTATIC switches.
  ! TODO: Replace the CONSTRAIN_HYDROSTATIC switches using this variable.
  if(CONSTRAIN_HYDROSTATIC) then
    cnst_hdrsttc=ONE
  else
    cnst_hdrsttc=ZERO
  endif
  
  ifirstelem = 1
  ilastelem = nspec
  
  rho_DG   = rho_DG_main
  rhovx_DG = rhovx_DG_main
  rhovz_DG = rhovz_DG_main
  E_DG     = E_DG_main
  
  ABSORB_BC = .false.
  if(ABSORB_BC) then
    maxval_rho   = maxval(rho_DG-rho_init)
    maxval_rhovx = maxval(rhovx_DG-rhovx_init)
    maxval_rhovz = maxval(rhovz_DG-rhovz_init)
    maxval_E     = maxval(E_DG-E_init)
  endif
  
  T_DG = T_DG_main
  V_DG = V_DG_main
  
  ! Initialise auxiliary unknowns.
  veloc_x_DG = rhovx_DG/rho_DG
  veloc_z_DG = rhovz_DG/rho_DG
  p_DG       = (gammaext_DG - ONE)*( E_DG &
               - (HALF)*rho_DG*( veloc_x_DG**2 + veloc_z_DG**2 ) )
        
  ! Initialisation.
  dot_rho   = ZERO
  dot_rhovx = ZERO
  dot_rhovz = ZERO
  dot_E     = ZERO
  dot_e1    = ZERO
  
  ! Add force source.
  if(TYPE_SOURCE_DG == 1) then
    call compute_add_sources_acoustic_DG_spread(dot_rho, it, i_stage)   
  elseif(TYPE_SOURCE_DG == 2) then
    call compute_add_sources_acoustic_DG_spread(dot_rhovx, it, i_stage)
    call compute_add_sources_acoustic_DG_spread(dot_rhovz, it, i_stage)
  elseif(TYPE_SOURCE_DG == 3) then
    call compute_add_sources_acoustic_DG_spread(dot_E, it, i_stage)
  endif
  
  ! TODO: introduce a verbosity parameter in order to prevent unwanted flooding of the terminal.
  if(myrank == 0 .AND. mod(it, 50) == 0) then
    write(*,"(a)")                 "               | max                     | min"
    WRITE(*,"(a,e24.16,a,e24.16)") " rho           |", maxval(rho_DG), " |", minval(rho_DG)
    WRITE(*,"(a,e24.16,a,e24.16)") " rhovx         |", maxval(rhovx_DG), " |", minval(rhovx_DG)
    WRITE(*,"(a,e24.16,a,e24.16)") " rhovz         |", maxval(rhovz_DG), " |", minval(rhovz_DG)
    WRITE(*,"(a,e24.16,a,e24.16)") " E             |", maxval(E_DG), " |", minval(E_DG)
    WRITE(*,"(a,e23.16,a)")        "Ratio |p-p_{init}|/p_{init}:", maxval(abs((p_DG-p_DG_init)/p_DG_init)), "."
  endif
  
  do ispec = ifirstelem, ilastelem ! Loop over elements.
    ! acoustic spectral element
    if (ispec_is_acoustic_DG(ispec)) then !if (ispec_is_acoustic(ispec)) then
      ! --------------------------- !
      ! First set of loops: compute !
      ! volumic contributions.      !
      ! --------------------------- !
      do j = 1, NGLLZ
        do i = 1, NGLLX
          iglob = ibool_DG(i, j, ispec)
          jacobianl = jacobian(i, j, ispec)
          xixl = xix(i, j, ispec)
          xizl = xiz(i, j, ispec)
          gammaxl = gammax(i, j, ispec)
          gammazl = gammaz(i, j, ispec)
          
          if(ABC_STRETCH) then
            ! Here is updated the operator \nabla and the jacobian.
            ! \partial_x becomes \ya_x\partial_x, and since \partial_x=(\partial_x\xi)\partial_\xi+(\partial_x\eta)\partial_\eta, only updating \partial_x\xi and \partial_x\eta is enough. Idem for \partial_z.
            ! Hence, only updating xix to \ya_x * xix, xiz to \ya_z * xiz, etc. is enough to update the operator.
            call virtual_stretch(i, j, ispec, coef_stretch_x_ij, coef_stretch_z_ij)
            xixl = coef_stretch_x_ij * xixl
            xizl = coef_stretch_z_ij * xizl
            gammaxl = coef_stretch_x_ij * gammaxl
            gammazl = coef_stretch_z_ij * gammazl
            ! Add jacobian of stretching into the integrand (artifically).
            ! TODO: Do that more clearly.
            jacobianl = coef_stretch_x_ij*coef_stretch_z_ij*jacobianl
          endif
          
          wzl = wzgll(j)
          wxl = wxgll(i)
          
          ! Inviscid stress tensor's contributions.
          temp_unknown = rhovx_DG(iglob)
          temp_unknown2 = rhovz_DG(iglob)
          
          temp_rho_1(i, j) = wzl * jacobianl * (xixl * temp_unknown + xizl * temp_unknown2)
          temp_rho_2(i, j) = wxl * jacobianl * (gammaxl * temp_unknown + gammazl * temp_unknown2)
          
          if(.not. CONSTRAIN_HYDROSTATIC) then
            temp_unknown = rho_DG(iglob)*veloc_x_DG(iglob)**2 + p_DG(iglob)
          else
            temp_unknown = rho_DG(iglob)*veloc_x_DG(iglob)**2 + (p_DG(iglob) - p_DG_init(iglob))
          endif
          temp_unknown2 = rho_DG(iglob)*veloc_x_DG(iglob)*veloc_z_DG(iglob)
          
          temp_rhovx_1(i, j) = wzl * jacobianl * (xixl * temp_unknown + xizl * temp_unknown2) 
          temp_rhovx_2(i, j) = wxl * jacobianl * (gammaxl * temp_unknown + gammazl * temp_unknown2) 
          
          temp_unknown = rho_DG(iglob)*veloc_x_DG(iglob)*veloc_z_DG(iglob)
          if(.not. CONSTRAIN_HYDROSTATIC) then
            temp_unknown2 = rho_DG(iglob)*veloc_z_DG(iglob)**2 + p_DG(iglob)
          else
            temp_unknown2 = rho_DG(iglob)*veloc_z_DG(iglob)**2 + (p_DG(iglob) - p_DG_init(iglob))
          endif
          
          temp_rhovz_1(i, j) = wzl * jacobianl * (xixl * temp_unknown + xizl * temp_unknown2) 
          temp_rhovz_2(i, j) = wxl * jacobianl * (gammaxl * temp_unknown + gammazl * temp_unknown2) 
          
          if(.not. CONSTRAIN_HYDROSTATIC) then
            temp_unknown = veloc_x_DG(iglob)*(E_DG(iglob) + p_DG(iglob))
            temp_unknown2 = veloc_z_DG(iglob)*(E_DG(iglob) + p_DG(iglob))
          else          
            temp_unknown = veloc_x_DG(iglob)*(E_DG(iglob) + (p_DG(iglob) - p_DG_init(iglob)))
            temp_unknown2 = veloc_z_DG(iglob)*(E_DG(iglob) + (p_DG(iglob) - p_DG_init(iglob)))
          endif
          
          temp_E_1(i, j) = wzl * jacobianl * (xixl * temp_unknown + xizl * temp_unknown2) 
          temp_E_2(i, j) = wxl * jacobianl * (gammaxl * temp_unknown + gammazl * temp_unknown2) 
          
          ! Add the viscous stress tensor's contributions.
          dux_dx = V_DG(1, 1, iglob)
          dux_dz = V_DG(1, 2, iglob)
          duz_dx = V_DG(2, 1, iglob)
          duz_dz = V_DG(2, 2, iglob)
          if(muext(i, j, ispec) > 0 .OR. &
             etaext(i, j, ispec) > 0 .OR. &
             kappa_DG(i, j, ispec) > 0) then
          
            dT_dx  = T_DG(1, iglob)
            dT_dz  = T_DG(2, iglob)
            
            ! This vector, [temp_unknown, temp_unknown2], is the first line of the viscous Navier-Stokes tensor (\Sigma_v).
            temp_unknown = muext(i, j, ispec)*TWO*dux_dx + (etaext(i, j, ispec) - (TWO/3.)*muext(i, j, ispec))*(dux_dx + duz_dz) 
            temp_unknown2 = muext(i, j, ispec)*( dux_dz + duz_dx )
            temp_rhovx_1(i, j) = temp_rhovx_1(i, j) - wzl * jacobianl * (xixl * temp_unknown + xizl * temp_unknown2) 
            temp_rhovx_2(i, j) = temp_rhovx_2(i, j) - wxl * jacobianl * (gammaxl * temp_unknown + gammazl * temp_unknown2) 
            
            ! Use the values stored in temp_unknown and temp_unknown2 to compute the x component of the first part of the viscous energy vector (\Sigma_v\cdot\vect{v}).
            temp_unknown  = veloc_x_DG(iglob)*temp_unknown + veloc_z_DG(iglob)*temp_unknown2
            temp_E_1(i, j) = temp_E_1(i, j) - wzl * jacobianl * (xixl * temp_unknown) 
            temp_E_2(i, j) = temp_E_2(i, j) - wxl * jacobianl * (gammaxl * temp_unknown) 
            
            ! This vector, [temp_unknown, temp_unknown2], is the second line of the viscous Navier-Stokes tensor (\Sigma_v).
            temp_unknown = muext(i, j, ispec)*( dux_dz + duz_dx )
            temp_unknown2 = muext(i, j, ispec)*TWO*duz_dz + (etaext(i, j, ispec) - (TWO/3.)*muext(i, j, ispec))*(dux_dx + duz_dz) 
            temp_rhovz_1(i, j) = temp_rhovz_1(i, j) - wzl * jacobianl * (xixl * temp_unknown + xizl * temp_unknown2) 
            temp_rhovz_2(i, j) = temp_rhovz_2(i, j) - wxl * jacobianl * (gammaxl * temp_unknown + gammazl * temp_unknown2) 
            
            ! Use the values stored in temp_unknown and temp_unknown2 to compute the z component of the first part of the viscous energy vector (\Sigma_v\cdot\vect{v}).
            temp_unknown2  = veloc_x_DG(iglob)*temp_unknown + veloc_z_DG(iglob)*temp_unknown2
            temp_E_1(i, j) = temp_E_1(i, j) - wzl * jacobianl * (xizl * temp_unknown2) 
            temp_E_2(i, j) = temp_E_2(i, j) - wxl * jacobianl * (gammazl * temp_unknown2) 
            
            ! Add the heat contributions (second part of the viscous energy tensor).
            temp_unknown  = kappa_DG(i, j, ispec)*dT_dx
            temp_unknown2 = kappa_DG(i, j, ispec)*dT_dz
            temp_E_1(i, j) = temp_E_1(i, j) - wzl * jacobianl * (xixl * temp_unknown + xizl * temp_unknown2) 
            temp_E_2(i, j) = temp_E_2(i, j) - wxl * jacobianl * (gammaxl * temp_unknown + gammazl * temp_unknown2) 
          endif
          
          
          ! Gravity contributions (separated from the rest).
          !temp_rho_gravi(i, j) = ZERO ! Save one operation.
          temp_rhovx_gravi(i, j) = -rho_DG(iglob)*potential_dphi_dx_DG(ibool(i, j, ispec))* jacobianl
          if(.not. CONSTRAIN_HYDROSTATIC) then
            temp_rhovz_gravi(i, j) = -rho_DG(iglob)*potential_dphi_dz_DG(ibool(i, j, ispec))* jacobianl
          else
            temp_rhovz_gravi(i, j) = -(rho_DG(iglob) - rho_init(iglob)) * potential_dphi_dz_DG(ibool(i, j, ispec)) * jacobianl 
          endif
          if(.not. CONSTRAIN_HYDROSTATIC) then
            temp_E_gravi(i, j) = -rho_DG(iglob)*(veloc_x_DG(iglob)*potential_dphi_dx_DG(ibool(i, j, ispec)) + &
                                veloc_z_DG(iglob)*potential_dphi_dz_DG(ibool(i, j, ispec)))* jacobianl
          else
            ! By constraining hydrostaticity, remark that an additional term appears (p_0 \nabla\cdot v'). We add it here despite the fact that it has nothing to do with gravity.
            temp_E_gravi(i, j) = &
                                -(rho_DG(iglob) - rho_init(iglob))*(veloc_x_DG(iglob)*potential_dphi_dx_DG(ibool(i, j, ispec)) + &
                                veloc_z_DG(iglob)*potential_dphi_dz_DG(ibool(i, j, ispec)))* jacobianl         
            !temp_E_gravi(i, j) = temp_E_gravi(i, j) - p_DG_init(iglob)*(dux_dx + duz_dz)* jacobianl
          endif
          temp_E_gravi(i, j) = temp_E_gravi(i, j) - jacobianl * (p_DG_init(iglob)*gammaext_DG(iglob)) &
                              * ( (tau_epsilon(i, j, ispec)/tau_sigma(i, j, ispec)) - ONE ) &
                              * ( dux_dx + duz_dz - e1_DG(iglob))/(gammaext_DG(iglob) - ONE)
          
          ! By constraining hydrostaticity, remark that an additional term appears ($p_0 \nabla\cdot v'$). Add it here despite the fact that it has nothing to do with gravity.
          if(CONSTRAIN_HYDROSTATIC) then
            temp_E_gravi(i, j) = temp_E_gravi(i, j) - p_DG_init(iglob)*(dux_dx + duz_dz)*jacobianl
          endif
          
          ! When using a virtual stretching method, more terms appear ($\Sigma\cdot(\nabla\cdot\Ya)$). Add them here despite the fact that they have nothing to do with gravity.
          if(ABC_STRETCH .and. .false.) then
            viscous_tens_11 =   TWO*muext(i, j, ispec)*dux_dx & ! Viscous tensor, line 1, column 1, term 1
                              + (etaext(i, j, ispec)-(TWO/3.)*muext(i, j, ispec)) & ! Viscous tensor, line 1, column 1, term 2 (1 of 2)
                                * (dux_dx+duz_dz) ! Viscous tensor, line 1, column 1, term 2 (2 of 2)
            viscous_tens_12 = muext(i, j, ispec)*(duz_dx+dux_dz) ! Visous tensor, line 1, column 2
            viscous_tens_22 =   TWO*muext(i, j, ispec)*duz_dz & ! Viscous tensor, line 2, column 2, term 1
                              + (etaext(i, j, ispec)-(TWO/3.)*muext(i, j, ispec)) & ! Viscous tensor, line 2, column 2, term 2 (1 of 2)
                                * (dux_dx+duz_dz) ! Viscous tensor, line 2, column 2, term 2 (2 of 2)
            
            call virtual_stretch_prime(i, j, ispec, coef_stretch_x_ij_prime, coef_stretch_z_ij_prime)
            
            if(.false. .and. i==3 .and. j==3 &
               .and. abs(coord(1, ibool_before_perio(i, j, ispec)))<5.&
               .and. abs(coord(2, ibool_before_perio(i, j, ispec))-21.)<2.) then
              write(*,*) 'omegalul', coord(1, ibool_before_perio(i, j, ispec)), coord(2, ibool_before_perio(i, j, ispec)) &
                         , coef_stretch_x_ij_prime, coef_stretch_z_ij_prime &
                         , ((rhovx_DG(iglob)*coef_stretch_x_ij_prime&
                             +rhovz_DG(iglob)*coef_stretch_z_ij_prime&
                            )* jacobianl*wxl*wzl)&
                         , rhovx_DG(iglob), rhovz_DG(iglob)
            endif
            if(.false. .and. i==3 .and. j==3 &
               .and. abs(coord(1, ibool_before_perio(i, j, ispec)))<5.&
               .and. abs(coord(2, ibool_before_perio(i, j, ispec))-21.)<2.) then
              write(*,*) 'omegalul', coord(1, ibool_before_perio(i, j, ispec)), coord(2, ibool_before_perio(i, j, ispec)) &
                         , coef_stretch_x_ij_prime, coef_stretch_z_ij_prime &
                         , (  (   (  E_DG(iglob) + p_DG(iglob)-cnst_hdrsttc*p_DG_init(iglob)&
                                                      + viscous_tens_11 &
                                                     ) * veloc_x_DG(iglob) &
                                                   + ( viscous_tens_12 ) * veloc_z_DG(iglob) &
                                                 ) * coef_stretch_x_ij_prime & ! End of z contribution.
                                               + (   (  E_DG(iglob) + p_DG(iglob)-cnst_hdrsttc*p_DG_init(iglob)&
                                                      + viscous_tens_22 &
                                                     ) * veloc_z_DG(iglob) &
                                                   + ( viscous_tens_12 ) * veloc_x_DG(iglob) &
                                                 ) * coef_stretch_z_ij_prime  & ! End of x contribution.
                                              ) * jacobianl*wxl*wzl
            endif
            
            !temp_rho_gravi(i, j) = temp_rho_gravi(i, j) + (  rhovx_DG(iglob)*coef_stretch_x_ij_prime &
            !                                               + rhovz_DG(iglob)*coef_stretch_z_ij_prime &
            !                                              ) * jacobianl
            temp_rhovx_gravi(i, j) = temp_rhovx_gravi(i, j) + (  (  rho_DG(iglob)&
                                                                    *veloc_x_DG(iglob)**2 & ! Inviscid tensor, line 1, column 1
                                                                    + p_DG(iglob)-cnst_hdrsttc*p_DG_init(iglob) & ! Inviscid tensor, line 1, column 1
                                                                  - viscous_tens_11 & ! Viscous tensor, line 1, column 1
                                                                 ) * coef_stretch_x_ij_prime  & ! End of x contribution.
                                                               + (  rho_DG(iglob)&
                                                                    *veloc_x_DG(iglob)*veloc_z_DG(iglob) & ! Inviscid tensor, line 1, column 2
                                                                  - viscous_tens_12 & ! Visous tensor, line 1, column 2
                                                                 ) * coef_stretch_z_ij_prime & ! End of z contribution.
                                                              ) * jacobianl
            temp_rhovz_gravi(i, j) = temp_rhovz_gravi(i, j) + (  (  rho_DG(iglob)&
                                                                    *veloc_x_DG(iglob)*veloc_z_DG(iglob) & ! Inviscid tensor, line 2, column 1
                                                                  - viscous_tens_12 & ! Visous tensor, line 2, column 1
                                                                 ) * coef_stretch_x_ij_prime & ! End of z contribution.
                                                               !+ (  (rho_DG(iglob)-cnst_hdrsttc*rho_init(iglob))&
                                                               + (  rho_DG(iglob)&
                                                                     *veloc_z_DG(iglob)**2  & ! Inviscid tensor
                                                                    + p_DG(iglob)-cnst_hdrsttc*p_DG_init(iglob) & ! Inviscid tensor, line 2, column 2
                                                                  - viscous_tens_22 & ! Viscous tensor, line 2, column 2
                                                                 ) * coef_stretch_z_ij_prime  & ! End of x contribution.
                                                              ) * jacobianl
            temp_E_gravi(i, j) = temp_E_gravi(i, j) + (  (   (  E_DG(iglob) + p_DG(iglob)-cnst_hdrsttc*p_DG_init(iglob)&
                                                              + viscous_tens_11 &
                                                             ) * veloc_x_DG(iglob) &
                                                           + ( viscous_tens_12 ) * veloc_z_DG(iglob) &
                                                         ) * coef_stretch_x_ij_prime & ! End of z contribution.
                                                       + (   ( viscous_tens_12 ) * veloc_x_DG(iglob) &
                                                           + (  E_DG(iglob) + p_DG(iglob)-cnst_hdrsttc*p_DG_init(iglob)&
                                                              + viscous_tens_22 &
                                                             ) * veloc_z_DG(iglob) &
                                                         ) * coef_stretch_z_ij_prime  & ! End of x contribution.
                                                      ) * jacobianl
          endif
          
          ! Memory variable evolution. TODO: Describe more precisely.
          dot_e1(iglob) = dot_e1(iglob) - (ONE/tau_sigma(i, j, ispec)) &
                          *( (ONE - (tau_sigma(i, j, ispec)/tau_epsilon(i, j, ispec))) * (dux_dx + duz_dz) + e1_DG(iglob) )
        enddo
      enddo
      
      ! Assemble the contributions previously computed, and add gravity's contribution.
      ! The integration by quadrature on the GLL points leads to three sums. See in particular Komatitsch (Méthodes spectrales et éléments spectraux pour l'équation de l'élastodynamique 2D et 3D en milieu hétérogène), Annexe 3.A.
      do j = 1, NGLLZ
        do i = 1, NGLLX
          iglob = ibool_DG(i, j, ispec)
          do k = 1, NGLLX
            dot_rho(iglob) = dot_rho(iglob) + &
                             (temp_rho_1(k, j) * real(hprimewgll_xx(k, i), kind=CUSTOM_REAL) + &
                              temp_rho_2(i, k) * real(hprimewgll_zz(k, j), kind=CUSTOM_REAL)   )
            dot_rhovx(iglob) = dot_rhovx(iglob) + &
                               (temp_rhovx_1(k, j) * real(hprimewgll_xx(k, i), kind=CUSTOM_REAL) + &
                                temp_rhovx_2(i, k) * real(hprimewgll_zz(k, j), kind=CUSTOM_REAL)   )
            dot_rhovz(iglob) = dot_rhovz(iglob) + &
                               (temp_rhovz_1(k, j) * real(hprimewgll_xx(k, i), kind=CUSTOM_REAL) + &
                                temp_rhovz_2(i, k) * real(hprimewgll_zz(k, j), kind=CUSTOM_REAL)   )
            dot_E(iglob) = dot_E(iglob) + &
                           (temp_E_1(k, j) * real(hprimewgll_xx(k, i), kind=CUSTOM_REAL) + &
                            temp_E_2(i, k) * real(hprimewgll_zz(k, j), kind=CUSTOM_REAL)   )
          enddo ! Enddo on k.
          
          wzl = real(wzgll(j), kind=CUSTOM_REAL)
          wxl = real(wxgll(i), kind=CUSTOM_REAL)
          
          ! Gravity terms.
          ! TODO: consider renaming the "temp_*_gravi" variables since some contain terms which are not related to gravity.
          !dot_rho(iglob)   = dot_rho(iglob)   + temp_rho_gravi(i, j) * wxl * wzl ! Save one operation.
          dot_rhovx(iglob) = dot_rhovx(iglob) + temp_rhovx_gravi(i, j) * wxl * wzl
          dot_rhovz(iglob) = dot_rhovz(iglob) + temp_rhovz_gravi(i, j) * wxl * wzl
          dot_E(iglob)     = dot_E(iglob)     + temp_E_gravi(i, j) * wxl * wzl
        enddo
      enddo
      
      ! --------------------------- !
      ! Second set of loops: add    !
      ! fluxes between elements,    !
      ! but only on exterior        !
      ! points.                     !
      ! --------------------------- !
      do  iface = 1, 4 
       do  iface1 = 1, NGLLX
          i = link_iface_ijispec(iface1,iface,ispec,1)
          j = link_iface_ijispec(iface1,iface,ispec,2)
          
          ! Step 1: prepare the normals' parameters (nx, nz, weight, etc.).
          ! Interior point
          iglobM = ibool_DG(i, j, ispec)
          
          ! TEST WITH IFACE FORMULATION
          nx     = nx_iface(iface, ispec)
          nz     = nz_iface(iface, ispec)
          weight = weight_iface(iface1,iface, ispec)
          neighbor = -1
          if(neighbor_DG_iface(iface1, iface, ispec, 3) > -1) then
            iface1_neighbor = neighbor_DG_iface(iface1, iface, ispec, 1)
            iface_neighbor  = neighbor_DG_iface(iface1, iface, ispec, 2)
            ispec_neighbor = neighbor_DG_iface(iface1, iface, ispec, 3)
            neighbor(1) = link_iface_ijispec(iface1_neighbor, iface_neighbor, ispec_neighbor,1)
            neighbor(2) = link_iface_ijispec(iface1_neighbor, iface_neighbor, ispec_neighbor,2)
            neighbor(3) = ispec_neighbor
          endif
          
          if(ABC_STRETCH) then
            ! Add jacobian of stretching into the integrand (artifically).
            ! TODO: Do that more clearly.
            call virtual_stretch(i, j, ispec, coef_stretch_x_ij, coef_stretch_z_ij)
            weight=coef_stretch_x_ij*coef_stretch_z_ij*weight
          endif
          
          ! Step 2: knowing the normals' parameters, compute now the fluxes.
          !TEST STRETCH
          nx_unit=nx
          nz_unit=nz
          if(ABC_STRETCH .and. .false.) then
            call virtual_stretch(i, j, ispec, coef_stretch_x_ij, coef_stretch_z_ij)
            nx = coef_stretch_x_ij * nx
            nz = coef_stretch_z_ij * nz
          endif
          
          rho_DG_P     = ZERO
          rhovx_DG_P   = ZERO
          rhovz_DG_P   = ZERO
          E_DG_P       = ZERO
          veloc_x_DG_P = ZERO
          veloc_z_DG_P = ZERO
          p_DG_P       = ZERO
          
          iglobP = 1
          if(neighbor(1) > -1) then
            iglobP = ibool_DG(neighbor(1), neighbor(2), neighbor(3))
          endif
          
          exact_interface_flux = .false. ! Reset this variable to .false.: by default, the fluxes have to be computed (jump!=0). In some specific cases (assigned during the call to compute_interface_unknowns), the flux can be exact (jump==0).
          call compute_interface_unknowns(i, j, ispec, rho_DG_P, rhovx_DG_P, &
                  rhovz_DG_P, E_DG_P, veloc_x_DG_P, veloc_z_DG_P, p_DG_P, T_P, &
                  Tx_DG_P, Tz_DG_P, Vxx_DG_P, Vzz_DG_P, Vzx_DG_P, Vxz_DG_P, gamma_P,&
                  neighbor, &
                  exact_interface_flux, &
                  rho_DG(iglobM), E_DG(iglobM), rhovx_DG(iglobM), rhovz_DG(iglobM), &
                  V_DG(:,:,iglobM), T_DG(:,iglobM), &
                  rho_DG(iglobP), E_DG(iglobP), rhovx_DG(iglobP), rhovz_DG(iglobP), &
                  V_DG(:,:,iglobP), T_DG(:,iglobP), &
                  nx, nz, weight, timelocal, iface1, iface)
                  !TEST STRETCH
                  !nx_unit, nz_unit, weight, timelocal, iface1, iface)
          
          ! Recover an approximate local maximum linearized acoustic wave speed. See for example Hesthaven (doi.org/10.1007/9780387720678), page 208.
          lambda = 0.
          jump   = 0.
          !gamma_P = gammaext_DG(iglobM) ! DEBUG
          !veloc_n_M = sqrt(veloc_x_DG(iglobM)**2 + veloc_z_DG(iglobM)**2) !DEBUG
          !veloc_n_P = sqrt(veloc_x_DG_P**2 + veloc_z_DG_P**2) !DEBUG
          
          !GOODVERSION
          !veloc_n_M = abs(veloc_x_DG(iglobM)*nx + veloc_z_DG(iglobM)*nz)
          !veloc_n_P = abs(veloc_x_DG_P*nx + veloc_z_DG_P*nz)
          !TEST STRETCH: use unit normal because we might be interested by the conventionnal fluid speed rather than the artificial one
          !veloc_n_M = abs(veloc_x_DG(iglobM)*nx_unit + veloc_z_DG(iglobM)*nz_unit)
          !veloc_n_P = abs(veloc_x_DG_P*nx_unit + veloc_z_DG_P*nz_unit)
          
          ! Save some operations.
          !lambda = max( veloc_n_M + sqrt(abs(gammaext_DG(iglobM)*p_DG(iglobM)/rho_DG(iglobM))), &
          !              veloc_n_P + sqrt(abs(gamma_P*p_DG_P/rho_DG_P)) )
          lambda = max(abs(veloc_x_DG(iglobM)*nx + veloc_z_DG(iglobM)*nz) &
                       + sqrt(abs(gammaext_DG(iglobM)*p_DG(iglobM)/rho_DG(iglobM))), &
                       abs(veloc_x_DG_P*nx + veloc_z_DG_P*nz) &
                       + sqrt(abs(gamma_P*p_DG_P/rho_DG_P)))
          
          ! Viscous stress tensor's contributions (already under the form of the average mean flux).
          dux_dx = ZERO
          dux_dz = ZERO
          duz_dx = ZERO
          duz_dz = ZERO
          dT_dx = ZERO
          dT_dz = ZERO
          if(muext(i, j, ispec) > 0 .OR. &
             etaext(i, j, ispec) > 0 .OR. &
             kappa_DG(i, j, ispec)  > 0) then
            dux_dx = HALF*(V_DG(1, 1, iglobM) + Vxx_DG_P)
            dux_dz = HALF*(V_DG(1, 2, iglobM) + Vxz_DG_P)
            duz_dx = HALF*(V_DG(2, 1, iglobM) + Vzx_DG_P)
            duz_dz = HALF*(V_DG(2, 2, iglobM) + Vzz_DG_P)
            dT_dx = HALF*(T_DG(1, iglobM) + Tx_DG_P)
            dT_dz = HALF*(T_DG(2, iglobM) + Tz_DG_P)
          endif
          
          ! Mass conservation equation's contributions (which are fully inviscid).
          !temp_unknown_M  = rhovx_DG(iglobM)
          !temp_unknown_P  = rhovx_DG_P
          !temp_unknown2_M = rhovz_DG(iglobM)
          !temp_unknown2_P = rhovz_DG_P
          ! Dot product.
          flux_n = (rhovx_DG(iglobM)+rhovx_DG_P)*nx + (rhovz_DG(iglobM)+rhovz_DG_P)*nz
          !flux_n = (temp_unknown_M+temp_unknown_P)*nx + (temp_unknown2_M+temp_unknown2_P)*nz ! [5 operations + 1 affectation], instead of [5 operations + 3 affectations].
          jump   = rho_DG(iglobM) - rho_DG_P
          ! Add flux' contribution.
          if(exact_interface_flux) then
            jump = ZERO
          endif
          dot_rho(iglobM) = dot_rho(iglobM) - weight*(flux_n + lambda*jump)*HALF
          
          ! x-Momentum equation's inviscid contributions.
          if(.not. CONSTRAIN_HYDROSTATIC) then
            temp_unknown_M = rho_DG(iglobM)*veloc_x_DG(iglobM)**2 + p_DG(iglobM)
            temp_unknown_P = rho_DG_P*veloc_x_DG_P**2 + p_DG_P
            temp_unknown2_M = rho_DG(iglobM)*veloc_x_DG(iglobM)*veloc_z_DG(iglobM)
            temp_unknown2_P = rho_DG_P*veloc_x_DG_P*veloc_z_DG_P
          else        
            temp_unknown_M = rho_DG(iglobM)*veloc_x_DG(iglobM)**2 + (p_DG(iglobM) - p_DG_init(iglobM))
            temp_unknown_P = rho_DG_P*veloc_x_DG_P**2 + (p_DG_P - p_DG_init(iglobM))
            temp_unknown2_M = rho_DG(iglobM)*veloc_x_DG(iglobM)*veloc_z_DG(iglobM)
            temp_unknown2_P = rho_DG_P*veloc_x_DG_P*veloc_z_DG_P
          endif
          ! Dot product.
          flux_n = (temp_unknown_M+temp_unknown_P)*nx + (temp_unknown2_M+temp_unknown2_P)*nz ! [5 operations + 1 affectation], instead of [5 operations + 3 affectations].
          jump   = rhovx_DG(iglobM) - rhovx_DG_P
          ! Add flux' contribution.
          if(exact_interface_flux) then
            jump = ZERO
          endif
          dot_rhovx(iglobM) = dot_rhovx(iglobM) - weight*(flux_n + lambda*jump)*HALF

          ! x-Momentum equation's viscous contributions.
          ! The vector [temp_unknown, temp_unknown2] represents the mean average flux at the boundary of the x-momentum.
          ! Recall: dux_dx, duz_dx, dux_dz, and duz_dz already contain the 0.5 factor to put the flux under mean average form.
          temp_unknown = muext(i, j, ispec)*TWO*dux_dx + (etaext(i, j, ispec) - (TWO/3.)*muext(i, j, ispec))*(dux_dx + duz_dz) 
          temp_unknown2 = muext(i, j, ispec)*( dux_dz + duz_dx )
          ! Dot product.
          flux_n = temp_unknown*nx + temp_unknown2*nz ! [3 operations + 1 affectation], instead of [3 operations + 3 affectations].
          dot_rhovx(iglobM) = dot_rhovx(iglobM) + weight*flux_n
          ! The computed values contained in the variables temp_unknown and temp_unknown2 can be used to compute the energy's x component of the mean average flux at the boundary. Thus, we add this contribution here.
          dot_E(iglobM)     = dot_E(iglobM) &
                              + weight * HALF * (  (veloc_x_DG(iglobM) + veloc_x_DG_P) * temp_unknown &
                                                  +(veloc_z_DG(iglobM) + veloc_z_DG_P) * temp_unknown2 )*nx
          
          ! z-Momentum equation's inviscid contributions.
          temp_unknown_M  = rho_DG(iglobM)*veloc_x_DG(iglobM)*veloc_z_DG(iglobM)
          temp_unknown_P  = rho_DG_P*veloc_x_DG_P*veloc_z_DG_P
          if(.not. CONSTRAIN_HYDROSTATIC) then        
            temp_unknown2_M = rho_DG(iglobM)*veloc_z_DG(iglobM)**2 + p_DG(iglobM)
            temp_unknown2_P = rho_DG_P*veloc_z_DG_P**2 + p_DG_P
          else
            temp_unknown2_M = rho_DG(iglobM)*veloc_z_DG(iglobM)**2 + (p_DG(iglobM) - p_DG_init(iglobM))
            temp_unknown2_P = rho_DG_P*veloc_z_DG_P**2 + (p_DG_P - p_DG_init(iglobM))
          endif
          ! Dot product.
          flux_n = (temp_unknown_M+temp_unknown_P)*nx + (temp_unknown2_M+temp_unknown2_P)*nz ! [5 operations + 1 affectation], instead of [5 operations + 3 affectations].
          jump   = rhovz_DG(iglobM) - rhovz_DG_P
          ! Add flux' contribution.
          if(exact_interface_flux) then
            jump = ZERO
          endif
          dot_rhovz(iglobM) = dot_rhovz(iglobM) - weight*(flux_n + lambda*jump)*HALF
          
          ! z-Momentum equation's viscous contributions.
          ! The vector [temp_unknown, temp_unknown2] represents the mean average flux at the boundary of the z-momentum.
          ! Recall: dux_dx, duz_dx, dux_dz, and duz_dz already contain the 0.5 factor to put the flux under mean average form.
          temp_unknown = muext(i, j, ispec)*( dux_dz + duz_dx )
          temp_unknown2 = muext(i, j, ispec)*TWO*duz_dz + (etaext(i, j, ispec) - (TWO/3.)*muext(i, j, ispec))*(dux_dx + duz_dz)
          ! Dot product.
          flux_n = temp_unknown*nx + temp_unknown2*nz ! [3 operations + 1 affectation], instead of [3 operations + 3 affectations].
          dot_rhovz(iglobM) = dot_rhovz(iglobM) + weight*flux_n
          ! The computed values contained in the variables temp_unknown and temp_unknown2 can be used to compute the energy's z component of the mean average flux at the boundary. Thus, we add this contribution here.
          dot_E(iglobM)     = dot_E(iglobM) &
                              + weight * HALF * (  (veloc_x_DG(iglobM) + veloc_x_DG_P) * temp_unknown &
                                                  +(veloc_z_DG(iglobM) + veloc_z_DG_P) * temp_unknown2 )*nz

          ! Energy equation's fully inviscid contributions.
          if(.not. CONSTRAIN_HYDROSTATIC) then
            temp_unknown_M = veloc_x_DG(iglobM)*(E_DG(iglobM) + p_DG(iglobM))
            temp_unknown_P = veloc_x_DG_P*(E_DG_P + p_DG_P)
            temp_unknown2_M = veloc_z_DG(iglobM)*(E_DG(iglobM) + p_DG(iglobM))
            temp_unknown2_P = veloc_z_DG_P*(E_DG_P + p_DG_P)
          else        
            temp_unknown_M = veloc_x_DG(iglobM)*(E_DG(iglobM) + (p_DG(iglobM)- p_DG_init(iglobM)))
            temp_unknown_P = veloc_x_DG_P*(E_DG_P + (p_DG_P - p_DG_init(iglobM)))
            temp_unknown2_M = veloc_z_DG(iglobM)*(E_DG(iglobM) + (p_DG(iglobM)- p_DG_init(iglobM)))
            temp_unknown2_P = veloc_z_DG_P*(E_DG_P + (p_DG_P - p_DG_init(iglobM)))
          endif        
          ! Dot product.
          flux_n = (temp_unknown_M+temp_unknown_P)*nx + (temp_unknown2_M+temp_unknown2_P)*nz ! [5 operations + 1 affectation], instead of [5 operations + 3 affectations].
          jump   = E_DG(iglobM) - E_DG_P
          ! Add flux' contribution.
          if(exact_interface_flux) then
            jump = ZERO
          endif
          dot_E(iglobM) = dot_E(iglobM) - weight*(flux_n + lambda*jump)*HALF
          
          ! Energy equation's heat flux' contribution (last remaining term, viscous).
          ! Recall: dT_dx, and dT_dx already contain the 0.5 factor to put the flux under mean average form.
          dot_E(iglobM) = dot_E(iglobM) &
                          + weight*( kappa_DG(i, j, ispec)*( dT_dx*nx + dT_dx*nz ) )
          
        enddo
      enddo
    endif ! End of test if acoustic element
  enddo ! End of loop on elements.
  
  end subroutine compute_forces_acoustic_DG
  
! ------------------------------------------------------------ !
! compute_viscous_tensors                                      !
! ------------------------------------------------------------ !
! Computes the values of the auxiliary viscous tensors \mathcal{T} and \mathcal{V} at every GLL point of the mesh.
! See doi:10.1016/j.jcp.2007.12.009, section 4.3.2.
   
  subroutine compute_viscous_tensors(T_DG, V_DG, rho_DG, rhovx_DG, rhovz_DG, E_DG, timelocal)

! ?? compute forces in the acoustic elements in forward simulation and in adjoint simulation in adjoint inversion
  
  use constants,only: CUSTOM_REAL,NGLLX,NGLLZ,gamma_euler

  use specfem_par, only: nglob_DG,nspec, ispec_is_acoustic_DG,&
                         xix,xiz,gammax,gammaz,jacobian, &
                         wxgll,wzgll, ibool_DG, &
                         hprimewgll_zz, hprimewgll_xx, &
                         hprime_xx, hprime_zz, rmass_inverse_acoustic_DG, &
                         cnu, &
                         rhovx_init, rhovz_init, rho_init, T_init, CONSTRAIN_HYDROSTATIC, &
                         link_iface_ijispec, nx_iface, nz_iface, weight_iface, neighbor_DG_iface,&
                         ABC_STRETCH ! Stretching-based absorbing conditions.
                         
  implicit none
  
  ! local parameters
  integer :: ispec,i, j,k, iglob, iglobM, iglobP
  real(kind=CUSTOM_REAL) :: rho_DG_P, rhovx_DG_P, rhovz_DG_P, &
        E_DG_P, veloc_x_DG_P, veloc_z_DG_P, p_DG_P, T_P, &
        Tx_DG_P, Tz_DG_P, Vxx_DG_P, Vzz_DG_P, Vzx_DG_P, Vxz_DG_P, &
        flux_n, flux_x, flux_z, nx, nz, timelocal, weight, gamma_P
  logical :: exact_interface_flux
  integer, dimension(nglob_DG) :: MPI_iglob
  integer, dimension(3) :: neighbor
  integer :: iface1, iface, iface1_neighbor, iface_neighbor, ispec_neighbor

  ! Jacobian matrix and determinant
  real(kind=CUSTOM_REAL) :: xixl,xizl,gammaxl,gammazl,jacobianl

  real(kind=CUSTOM_REAL) :: temp_Tx, temp_Tz, temp_Vxx, temp_Vzx, temp_Vxz, temp_Vzz
  
  ! Local variables.
  real(kind=CUSTOM_REAL), dimension(2, 2, nglob_DG), intent(out) :: V_DG
  real(kind=CUSTOM_REAL), dimension(2, nglob_DG), intent(out) :: T_DG
  real(kind=CUSTOM_REAL), dimension(nglob_DG) :: &
         rho_DG, rhovx_DG, rhovz_DG, E_DG, veloc_x_DG, veloc_z_DG, T, &
         grad_Tx, grad_Tz, grad_Vxx, grad_Vzz, grad_Vxz, grad_Vzx
  
  ! Viscosity
  real(kind=CUSTOM_REAL) :: dux_dxi, dux_dgamma, duz_dxi, duz_dgamma, dT_dxi, dT_dgamma
  real(kind=CUSTOM_REAL) :: dux_dx, dux_dz, duz_dx, duz_dz, dT_dx, dT_dz
  real(kind=CUSTOM_REAL) :: wxl, wzl
  
  ! Parameters
  real(kind=CUSTOM_REAL), parameter :: ZERO = 0._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: ONE  = 1._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: TWO  = 2._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: HALF = 0.5_CUSTOM_REAL
  
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: temp_Tx_1, temp_Tx_2, &
        temp_Tz_1, temp_Tz_2, temp_Vxx_1, temp_Vxx_2, &
        temp_Vxz_1, temp_Vxz_2, temp_Vzx_1, temp_Vzx_2, temp_Vzz_1, temp_Vzz_2
  
  logical :: ADD_SURFACE_TERMS
  real(kind=CUSTOM_REAL) :: vx_init, vz_init
  
!  integer :: coef_surface
  
  ADD_SURFACE_TERMS = .true. ! TODO: Decide to set a value (.true. or .false.), or to introduce a parameter in the parfile.
  
  !T_init = (E_DG/rho_DG - 0.5*((rhovx_DG/rho_DG)**2 + (rhovz_DG/rho_DG)**2))/(cnu)
  
  MPI_iglob = 1
  
  veloc_x_DG = rhovx_DG/rho_DG
  veloc_z_DG = rhovz_DG/rho_DG
  T = (E_DG/rho_DG - 0.5*(veloc_x_DG**2 + veloc_z_DG**2))/(cnu)
  
  grad_Tx = ZERO
  grad_Tz = ZERO
  grad_Vxx = ZERO
  grad_Vzz = ZERO
  grad_Vxz = ZERO
  grad_Vzx = ZERO

  do ispec = 1, nspec ! Loop over elements.
    ! acoustic spectral element
    !if (ispec_is_acoustic(ispec)) then
    if (ispec_is_acoustic_DG(ispec)) then
      
      ! --------------------------- !
      ! Volume terms.               !
      ! --------------------------- !
      do j = 1, NGLLZ
        do i = 1, NGLLX
          iglob = ibool_DG(i, j, ispec)
          jacobianl = jacobian(i, j, ispec)
          xixl = xix(i, j, ispec)
          xizl = xiz(i, j, ispec)
          gammaxl = gammax(i, j, ispec)
          gammazl = gammaz(i, j, ispec)
          
          wzl = wzgll(j)
          wxl = wxgll(i)
          
          if(ADD_SURFACE_TERMS) then
            ! In that case, we want to compute \int_{\Omega^k} \mathcal{T}\Phi d\Omega^k = \int_{\Omega^k} (\nabla T)\Phi d\Omega^k as:
            ! - \int_{\Omega^k} T\nabla\Phi d\Omega^k + \int_{\partial\Omega^k} T\Phi n d\Gamma^k.
            ! The idea is to store in:
            !   grad_Tx,
            !   grad_Tz,
            !   grad_Vxx,
            !   grad_Vxz,
            !   grad_Vzx,
            !   grad_Vzz
            ! the values of the approximated integrals:
            !  $\int \partial_xT\Phi_x d\Omega^k$,
            !  $\int \partial_zT\Phi_z d\Omega^k$,
            !  $\int \partial_xVx\Phi_x d\Omega^k$,
            !  etc.
            ! and "desintegrate" those values by multiplying the obtained vector by the inverse mass matrix. As explained before, the integrals are computed by using the divergence theorem and taking into account the surface terms.
            
            ! Viscous stress tensors
            if(.not. CONSTRAIN_HYDROSTATIC) then
              temp_Tx_1(i, j)  = wzl * jacobianl * (xixl * T(iglob)) 
              temp_Tz_1(i, j)  = wzl * jacobianl * (xizl * T(iglob)) 
              temp_Vxx_1(i, j) = wzl * jacobianl * (xixl * veloc_x_DG(iglob))
              temp_Vxz_1(i, j) = wzl * jacobianl * (xizl * veloc_x_DG(iglob))
              temp_Vzx_1(i, j) = wzl * jacobianl * (xixl * veloc_z_DG(iglob)) 
              temp_Vzz_1(i, j) = wzl * jacobianl * (xizl * veloc_z_DG(iglob)) 
              
              temp_Tx_2(i, j)  = wxl * jacobianl * (gammaxl * T(iglob)) 
              temp_Tz_2(i, j)  = wxl * jacobianl * (gammazl * T(iglob)) 
              temp_Vxx_2(i, j) = wxl * jacobianl * (gammaxl * veloc_x_DG(iglob)) 
              temp_Vxz_2(i, j) = wxl * jacobianl * (gammazl * veloc_x_DG(iglob)) 
              temp_Vzx_2(i, j) = wxl * jacobianl * (gammaxl * veloc_z_DG(iglob)) 
              temp_Vzz_2(i, j) = wxl * jacobianl * (gammazl * veloc_z_DG(iglob)) 
            else
              vx_init = rhovx_init(iglob)/rho_init(iglob)
              vz_init = rhovz_init(iglob)/rho_init(iglob)
              
              temp_Tx_1(i, j)  = wzl * jacobianl * (xixl * (T(iglob) - T_init(iglob)))
              temp_Tz_1(i, j)  = wzl * jacobianl * (xizl * (T(iglob) - T_init(iglob)))
              temp_Vxx_1(i, j) = wzl * jacobianl * (xixl * (veloc_x_DG(iglob) - vx_init))
              temp_Vxz_1(i, j) = wzl * jacobianl * (xizl * (veloc_x_DG(iglob) - vx_init))
              temp_Vzx_1(i, j) = wzl * jacobianl * (xixl * (veloc_z_DG(iglob) - vz_init))
              temp_Vzz_1(i, j) = wzl * jacobianl * (xizl * (veloc_z_DG(iglob) - vz_init))
              
              temp_Tx_2(i, j)  = wxl * jacobianl * (gammaxl * (T(iglob) - T_init(iglob)))
              temp_Tz_2(i, j)  = wxl * jacobianl * (gammazl * (T(iglob) - T_init(iglob)))
              temp_Vxx_2(i, j) = wxl * jacobianl * (gammaxl * (veloc_x_DG(iglob) - vx_init))
              temp_Vxz_2(i, j) = wxl * jacobianl * (gammazl * (veloc_x_DG(iglob) - vx_init))
              temp_Vzx_2(i, j) = wxl * jacobianl * (gammaxl * (veloc_z_DG(iglob) - vz_init))
              temp_Vzz_2(i, j) = wxl * jacobianl * (gammazl * (veloc_z_DG(iglob) - vz_init))
            endif
          else
            ! In that case, we want to compute \int_{\Omega^k} \mathcal{T}\Phi d\Omega^k = \int_{\Omega^k} (\nabla T)\Phi d\Omega^k directly as it.
            ! The idea is to store in:
            !   grad_Tx,
            !   grad_Tz,
            !   grad_Vxx,
            !   grad_Vxz,
            !   grad_Vzx,
            !   grad_Vzz
            ! the actual values of the quantities:
            !   \partial_xT,
            !   \partial_zT,
            !   \partial_xVx,
            !   etc.
            ! which is immediate through the SEM formulation.
            
            dux_dxi    = ZERO
            dux_dgamma = ZERO
            duz_dxi    = ZERO
            duz_dgamma = ZERO
            dT_dxi     = ZERO
            dT_dgamma  = ZERO
            
            ! first double loop over GLL points to compute and store gradients
            ! we can merge the two loops because NGLLX == NGLLZ
            do k = 1, NGLLX
              if(.not. CONSTRAIN_HYDROSTATIC) then
                dux_dxi    = dux_dxi    + veloc_x_DG(ibool_DG(k, j,ispec)) * real(hprime_xx(i, k), kind=CUSTOM_REAL)
                dux_dgamma = dux_dgamma + veloc_x_DG(ibool_DG(i, k, ispec)) * real(hprime_zz(j,k), kind=CUSTOM_REAL)
                duz_dxi    = duz_dxi    + veloc_z_DG(ibool_DG(k, j,ispec)) * real(hprime_xx(i, k), kind=CUSTOM_REAL)
                duz_dgamma = duz_dgamma + veloc_z_DG(ibool_DG(i, k, ispec)) * real(hprime_zz(j,k), kind=CUSTOM_REAL)
                dT_dxi     = dT_dxi    + T(ibool_DG(k, j,ispec)) * real(hprime_xx(i, k), kind=CUSTOM_REAL)
                dT_dgamma  = dT_dgamma + T(ibool_DG(i, k, ispec)) * real(hprime_zz(j,k), kind=CUSTOM_REAL)
              else
                vx_init = rhovx_init(ibool_DG(k, j,ispec))/rho_init(ibool_DG(k, j,ispec))
                dux_dxi = dux_dxi + (veloc_x_DG(ibool_DG(k, j,ispec)) - vx_init) &
                                    * real(hprime_xx(i, k), kind=CUSTOM_REAL)
                vx_init = rhovx_init(ibool_DG(i, k, ispec))/rho_init(ibool_DG(i, k, ispec))       
                dux_dgamma = dux_dgamma + (veloc_x_DG(ibool_DG(i, k, ispec)) - vx_init) &
                                          * real(hprime_zz(j,k), kind=CUSTOM_REAL)
                vz_init = rhovz_init(ibool_DG(k, j,ispec))/rho_init(ibool_DG(k, j,ispec))
                duz_dxi = duz_dxi + (veloc_z_DG(ibool_DG(k, j,ispec)) - vz_init) &
                                    * real(hprime_xx(i, k), kind=CUSTOM_REAL)
                vz_init = rhovz_init(ibool_DG(i, k, ispec))/rho_init(ibool_DG(i, k, ispec))
                duz_dgamma = duz_dgamma + (veloc_z_DG(ibool_DG(i, k, ispec)) - vz_init) &
                                          * real(hprime_zz(j,k), kind=CUSTOM_REAL)
                dT_dxi = dT_dxi + (T(ibool_DG(k, j,ispec)) - T_init(ibool_DG(k, j,ispec))) &
                                  * real(hprime_xx(i, k), kind=CUSTOM_REAL)
                dT_dgamma = dT_dgamma + (T(ibool_DG(i, k, ispec)) - T_init(ibool_DG(i, k, ispec))) &
                                         * real(hprime_zz(j,k), kind=CUSTOM_REAL)
              endif
            enddo

            ! Derivatives of velocities.
            dux_dx = ( dux_dxi * xixl + dux_dgamma * gammaxl )
            dux_dz = ( dux_dxi * xizl + dux_dgamma * gammazl )
            duz_dx = ( duz_dxi * xixl + duz_dgamma * gammaxl )
            duz_dz = ( duz_dxi * xizl + duz_dgamma * gammazl )
            ! Derivatives of temperature.
            dT_dx = ( dT_dxi * xixl + dT_dgamma * gammaxl )
            dT_dz = ( dT_dxi * xizl + dT_dgamma * gammazl )
            
            temp_Tx  = dT_dx * jacobianl
            temp_Tz  = dT_dz * jacobianl
            temp_Vxx = dux_dx * jacobianl
            temp_Vxz = dux_dz * jacobianl
            temp_Vzx = duz_dx * jacobianl
            temp_Vzz = duz_dz * jacobianl
            
            grad_Tx(iglob)  = dT_dx
            grad_Tz(iglob)  = dT_dz
            grad_Vxx(iglob) = dux_dx
            grad_Vxz(iglob) = dux_dz
            grad_Vzx(iglob) = duz_dx
            grad_Vzz(iglob) = duz_dz
          endif
        enddo
      enddo
      
      if(ADD_SURFACE_TERMS) then
        do j = 1, NGLLZ
          do i = 1, NGLLX
            iglob = ibool_DG(i, j, ispec)
            ! along x direction and z direction
            ! and assemble the contributions
            do k = 1, NGLLX
              grad_Tx(iglob) = grad_Tx(iglob) - &
                               (temp_Tx_1(k, j) * real(hprimewgll_xx(k, i), kind=CUSTOM_REAL) + &
                               temp_Tx_2(i, k) * real(hprimewgll_zz(k, j), kind=CUSTOM_REAL))
              grad_Tz(iglob) = grad_Tz(iglob) - &
                               (temp_Tz_1(k, j) * real(hprimewgll_xx(k, i), kind=CUSTOM_REAL) + &
                               temp_Tz_2(i, k) * real(hprimewgll_zz(k, j), kind=CUSTOM_REAL))
              grad_Vxx(iglob) = grad_Vxx(iglob) - &
                                (temp_Vxx_1(k, j) * real(hprimewgll_xx(k, i), kind=CUSTOM_REAL) + &
                                temp_Vxx_2(i, k) * real(hprimewgll_zz(k, j), kind=CUSTOM_REAL))
              grad_Vxz(iglob) = grad_Vxz(iglob) - &
                                (temp_Vxz_1(k, j) * real(hprimewgll_xx(k, i), kind=CUSTOM_REAL) + &
                                temp_Vxz_2(i, k) * real(hprimewgll_zz(k, j), kind=CUSTOM_REAL))
              grad_Vzx(iglob) = grad_Vzx(iglob) - &
                                (temp_Vzx_1(k, j) * real(hprimewgll_xx(k, i), kind=CUSTOM_REAL) + &
                                temp_Vzx_2(i, k) * real(hprimewgll_zz(k, j), kind=CUSTOM_REAL))
              grad_Vzz(iglob) = grad_Vzz(iglob) - &
                                (temp_Vzz_1(k, j) * real(hprimewgll_xx(k, i), kind=CUSTOM_REAL) + &
                                temp_Vzz_2(i, k) * real(hprimewgll_zz(k, j), kind=CUSTOM_REAL))
            enddo
          enddo
        enddo
          
      ! --------------------------- !
      ! Interface terms.            !
      ! --------------------------- !
        do  iface = 1, 4 
         do  iface1 = 1, NGLLX
         
            i = link_iface_ijispec(iface1,iface,ispec,1)
            j = link_iface_ijispec(iface1,iface,ispec,2)

            ! Interior point
            iglobM = ibool_DG(i, j, ispec)
            
            rho_DG_P     = ZERO
            rhovx_DG_P   = ZERO
            rhovz_DG_P   = ZERO
            E_DG_P       = ZERO
            veloc_x_DG_P = ZERO
            veloc_z_DG_P = ZERO
            p_DG_P       = ZERO
            T_P          = ZERO
            Vxx_DG_P     = ZERO
            Vzz_DG_P     = ZERO
            Vzx_DG_P     = ZERO
            Vxz_DG_P     = ZERO
            
            ! TEST WITH IFACE FORMULATION
            nx     = nx_iface(iface, ispec)
            nz     = nz_iface(iface, ispec)
            
            weight = weight_iface(iface1,iface, ispec)
            neighbor = -1
            if(neighbor_DG_iface(iface1, iface, ispec, 3) > -1) then
              iface1_neighbor = neighbor_DG_iface(iface1, iface, ispec, 1)
              iface_neighbor  = neighbor_DG_iface(iface1, iface, ispec, 2)
              ispec_neighbor = neighbor_DG_iface(iface1, iface, ispec, 3)
              neighbor(1) = link_iface_ijispec(iface1_neighbor, iface_neighbor, ispec_neighbor,1)
              neighbor(2) = link_iface_ijispec(iface1_neighbor, iface_neighbor, ispec_neighbor,2)
              neighbor(3) = ispec_neighbor
            endif
            
            iglobP = 1
            if(neighbor(1) > -1) then
              iglobP = ibool_DG(neighbor(1), neighbor(2), neighbor(3))
            endif
            
            exact_interface_flux = .false. ! Reset this variable to .false.: by default, the fluxes have to be computed (jump!=0). In some specific cases (assigned during the call to compute_interface_unknowns), the flux can be exact (jump==0).
            call compute_interface_unknowns(i, j, ispec, rho_DG_P, rhovx_DG_P, &
                    rhovz_DG_P, E_DG_P, veloc_x_DG_P, veloc_z_DG_P, p_DG_P, T_P, &
                    Tx_DG_P, Tz_DG_P, Vxx_DG_P, Vzz_DG_P, Vzx_DG_P, Vxz_DG_P, gamma_P, &
                    neighbor,&
                    exact_interface_flux, &
                    rho_DG(iglobM), E_DG(iglobM), rhovx_DG(iglobM), rhovz_DG(iglobM), &
                    V_DG(:,:,iglobM), T_DG(:,iglobM), &
                    rho_DG(iglobP), E_DG(iglobP), rhovx_DG(iglobP), rhovz_DG(iglobP), &
                    V_DG(:,:,iglobP), T_DG(:,iglobP), &
                    nx, nz, weight, timelocal,iface1, iface)

            vx_init = rhovx_init(iglobM)/rho_init(iglobM)
            vz_init = rhovz_init(iglobM)/rho_init(iglobM)

            ! Dot product.
            flux_x = T(iglobM) + T_P
            if(CONSTRAIN_HYDROSTATIC) flux_x = flux_x - 2*T_init(iglobM)
            flux_n = flux_x*nx
            grad_Tx(iglobM) = grad_Tx(iglobM) + weight*flux_n*HALF
            
            ! Dot product.
            flux_z = T(iglobM) + T_P
            if(CONSTRAIN_HYDROSTATIC) flux_z = flux_z - 2*T_init(iglobM)
            flux_n = flux_z*nz
            grad_Tz(iglobM) = grad_Tz(iglobM) + weight*flux_n*HALF
            
            ! Dot product.
            flux_x = veloc_x_DG(iglobM) + veloc_x_DG_P
            if(CONSTRAIN_HYDROSTATIC) flux_x = flux_x - 2*vx_init
            flux_n = flux_x*nx
            grad_Vxx(iglobM) = grad_Vxx(iglobM) + weight*flux_n*HALF
            
            ! Dot product.
            flux_z = veloc_x_DG(iglobM) + veloc_x_DG_P
            if(CONSTRAIN_HYDROSTATIC) flux_z = flux_z - 2*vx_init
            flux_n = flux_z*nz
            grad_Vxz(iglobM) = grad_Vxz(iglobM) + weight*flux_n*HALF
            
            ! Dot product.
            flux_x = veloc_z_DG(iglobM) + veloc_z_DG_P
            if(CONSTRAIN_HYDROSTATIC) flux_x = flux_x - 2*vz_init
            flux_n = flux_x*nx
            grad_Vzx(iglobM) = grad_Vzx(iglobM) + weight*flux_n*HALF
            
            ! Dot product.
            flux_z = veloc_z_DG(iglobM) + veloc_z_DG_P
            if(CONSTRAIN_HYDROSTATIC) flux_z = flux_z - 2*vz_init
            flux_n = flux_z*nz
            grad_Vzz(iglobM) = grad_Vzz(iglobM) + weight*flux_n*HALF
            
          enddo
        enddo
      endif ! Endif on ADD_SURFACE_TERMS.
    endif ! End of test if acoustic element.
  enddo ! End of loop over elements.
  
  if(ADD_SURFACE_TERMS) then
    ! "Desintegrate".
    grad_Tx(:)  = grad_Tx(:) * rmass_inverse_acoustic_DG(:)
    grad_Tz(:)  = grad_Tz(:) * rmass_inverse_acoustic_DG(:)
    grad_Vxx(:) = grad_Vxx(:) * rmass_inverse_acoustic_DG(:)
    grad_Vxz(:) = grad_Vxz(:) * rmass_inverse_acoustic_DG(:)
    grad_Vzx(:) = grad_Vzx(:) * rmass_inverse_acoustic_DG(:)
    grad_Vzz(:) = grad_Vzz(:) * rmass_inverse_acoustic_DG(:)
  endif
  
  ! Store in variables intended for output.
  T_DG(1, :) = grad_Tx
  T_DG(2, :) = grad_Tz
  V_DG(1, 1, :) = grad_Vxx
  V_DG(1, 2, :) = grad_Vxz
  V_DG(2, 1, :) = grad_Vzx
  V_DG(2, 2, :) = grad_Vzz
  
  ! --------------------------- !
  ! Virtual mesh stretching.    !
  ! --------------------------- !
  if(ABC_STRETCH) then
    ! TODO: Maybe, this has to be put in the developements above rather than "hardcoded" here.
    call virtual_stretch_auxi_visc_tensors_DG(T_DG, V_DG)
  endif
  
  end subroutine compute_viscous_tensors

! ------------------------------------------------------------ !
! virtual_stretch_auxi_visc_tensors_DG                         !
! ------------------------------------------------------------ !
! TODO: This is a test routine.
! Multiplies the terms in the auxiliary viscous tensors \mathcal{T} and \mathcal{V} by the corresponding stretching functions. The stretching consists in replacing the derivatives \partial_x and \partial_z by R_x\partial_x and R_y\partial_y.
! TODO: Instead of recomputing the coefficients at each iteration, a 2 * nglob_DG matrix containing them for all mesh points can be computed and stored before execution, and called to here.

  subroutine virtual_stretch_auxi_visc_tensors_DG(T_DG, V_DG)

  use specfem_par, only: ispec_is_acoustic, nspec, nglob_DG, ibool_DG
  use constants,only: CUSTOM_REAL,NGLLX,NGLLZ

  implicit none
  
  ! Input/Output
  real(kind=CUSTOM_REAL), dimension(2, 2, nglob_DG), intent(inout) :: V_DG
  real(kind=CUSTOM_REAL), dimension(2, nglob_DG), intent(inout) :: T_DG
  
  ! Local
  integer :: i, j, ispec, iglob
  real(kind=CUSTOM_REAL) :: coef_stretch_x, coef_stretch_z
  
  do ispec = 1, nspec
    if (ispec_is_acoustic(ispec)) then
      do j = 1, NGLLZ
        do i = 1, NGLLX
          iglob = ibool_DG(i, j, ispec)
          call virtual_stretch(i, j, ispec, coef_stretch_x, coef_stretch_z) ! Get coef_stretch_x and coef_stretch_z based on the element's coordinates.
          T_DG(1, iglob) = coef_stretch_x * T_DG(1, iglob)
          T_DG(2, iglob) = coef_stretch_z * T_DG(2, iglob)
          V_DG(1, :, iglob) = coef_stretch_x * V_DG(1, :, iglob)
          V_DG(2, :, iglob) = coef_stretch_z * V_DG(2, :, iglob)
        enddo
      enddo
     endif
   enddo
  end subroutine virtual_stretch_auxi_visc_tensors_DG

! ------------------------------------------------------------ !
! virtual_stretch                                              !
! ------------------------------------------------------------ !
! Computes the value of the stretching coefficients.
! TODO: Instead of a function to be called, build a vector that will only need to be computed once and to which simple and less expensive memory calls can be made.

  subroutine virtual_stretch(i, j, ispec, coef_stretch_x, coef_stretch_z)
  
  use specfem_par, only: ibool_DG, ibool_before_perio, coord, ABC_STRETCH_LBUF
  use constants,only: CUSTOM_REAL,NGLLX,NGLLZ

  implicit none
  
  ! Input/Output
  integer, intent(in) :: i, j, ispec
  real(kind=CUSTOM_REAL), intent(out) :: coef_stretch_x, coef_stretch_z
  
  ! Local
  integer iglob
  real(kind=CUSTOM_REAL), parameter :: ONE = 1._CUSTOM_REAL
  real(kind=CUSTOM_REAL) :: eps_l, p, q ! Arina's stretching.
  real(kind=CUSTOM_REAL) :: C_1, C_2 ! Arina's damping.
  real(kind=CUSTOM_REAL) :: beta, sigma_max ! Richards' damping.
  real(kind=CUSTOM_REAL) :: zmax, z, z_l ! Variables used for testing.
  
  ! Coefficients for the stretching function.
  ! Arina's stretching
  eps_l = 1.0d-4 ! 1.d-4 in Arina's paper.
  p = 3.25d0
  q = 1.75d0
  ! Arina's damping coefficients.
  C_1 = 0.0d0 ! 0 in Arina's paper, 0<=C_1<=0.1 in Wasistho's.
  C_2 = 13.0d0 ! 13 in Arina's paper, 10<=C_2<=20 in Wasistho's.
  ! Richards' damping coefficients.
  beta = 2.0d0 ! 4 in Richards' paper.
  sigma_max = -1.0d0
  
  zmax = 20. ! Domain top coordinate.
  
  iglob = ibool_DG(i, j, ispec)
  coef_stretch_x = ONE ! By default, mesh is not stretched.
  coef_stretch_z = ONE ! By default, mesh is not stretched.
  z = coord(2, ibool_before_perio(i, j, ispec)) ! Absolute coordinate.
  z_l = (1. - ((zmax - z)/ABC_STRETCH_LBUF)) ! Relative buffer coordinate.
  if(z_l > 0. .AND. z_l <= 1.) then
    ! Note: new expressions can be implemented here. Some compatibility conditions should be respected. The function has be 1 when z_l==0. The derivative of the function should be 0 when z_l==0 and when z_l==1.
    !write(*, *) "z", z, "z_l", z_l ! DEBUG
    !coef_stretch_z = 1. + z_l**2.! Stretching function.
    !coef_stretch_z = 1. + 5.*z_l**2.*(z_l-1.)**2.! Stretching function.
    !coef_stretch_z = 1. + 1. - (1. - eps_l) * (1. - z_l**p)**q! Stretching function.
    !coef_stretch_z = 1. - 1.*z_l ! Stretching function.
    ! Arina's stretching.
    coef_stretch_z = 1. - (1. - eps_l) * (1. - (1. - z_l)**p)**q ! Stretching function.
    ! Arina's damping.
    !coef_stretch_z = (1.0d0 - C_1*z_l**2.0d0)*(1.0d0 - ( 1.0d0 - exp(C_2*(z_l)**2.0d0) )/( 1.0d0 - exp(C_2) ))
    ! Richards' damping.
    !coef_stretch_z = 1.0d0 + sigma_max * z_l ** beta
    !write(*, *) "z", z, "coef_stretch_z", coef_stretch_z ! DEBUG
  endif
  
  end subroutine virtual_stretch

! ------------------------------------------------------------ !
! virtual_stretch_prime                                        !
! ------------------------------------------------------------ !
! Computes the value of the derivative of the stretching coefficients.
! TODO: Instead of a function to be called, build a vector that will only need to be computed once and to which simple and less expensive memory calls can be made.

  subroutine virtual_stretch_prime(i, j, ispec, coef_stretch_x_prime, coef_stretch_z_prime)
  
  use specfem_par, only: ibool_DG, ibool_before_perio, coord, ABC_STRETCH_LBUF
  use constants,only: CUSTOM_REAL,NGLLX,NGLLZ

  implicit none
  
  ! Input/Output
  integer, intent(in) :: i, j, ispec
  real(kind=CUSTOM_REAL), intent(out) :: coef_stretch_x_prime, coef_stretch_z_prime
  
  ! Local
  integer iglob
  real(kind=CUSTOM_REAL), parameter :: ZERO = 0._CUSTOM_REAL
  real(kind=CUSTOM_REAL) :: eps_l, p, q
  real(kind=CUSTOM_REAL) :: zmax, z, z_l ! Variables used for testing.
  
  stop 'should not be called'
  
  ! Coefficients for the stretching function.
  eps_l = 1.0d-4
  p = 3.25d0
  q = 1.75d0
  
  zmax = 50. ! Domain top coordinate.
  
  iglob = ibool_DG(i, j, ispec)
  coef_stretch_x_prime = ZERO ! By default, mesh is not stretched.
  coef_stretch_z_prime = ZERO ! By default, mesh is not stretched.
  z = coord(2, ibool_before_perio(i, j, ispec)) ! Absolute coordinate.
  z_l = (1. - ((zmax - z)/ABC_STRETCH_LBUF)) ! Relative buffer coordinate.
  if(z_l > 0. .AND. z_l <= 1.) then
    !write(*, *) "z", z, "z_l", z_l ! DEBUG
    coef_stretch_z_prime = - (1. - eps_l) * p * q * (1.-z_l)**(p-1.) * (1.-(1.-z_l)**p)**(q-1.) ! Stretching function derivative.
    !coef_stretch_z_prime = -1. ! Stretching function derivative.
    !write(*, *) "z", z, "coef_stretch_z", coef_stretch_z ! DEBUG
  endif
  
  end subroutine virtual_stretch_prime

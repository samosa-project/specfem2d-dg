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
                         ABC_STRETCH, stretching_ya, &!stretching_buffer,&! Stretching-based absorbing conditions.
                         ABC_STRETCH_LEFT, ABC_STRETCH_RIGHT, ABC_STRETCH_TOP, ABC_STRETCH_BOTTOM,&
                         ABC_STRETCH_LEFT_LBUF, ABC_STRETCH_RIGHT_LBUF, ABC_STRETCH_TOP_LBUF, ABC_STRETCH_BOTTOM_LBUF,&
                         mesh_xmin, mesh_xmax, mesh_zmin, mesh_zmax,&
                         coord, ibool_before_perio,stretching_buffer!,cnu
                         
  implicit none

  ! Input/output.
  real(kind=CUSTOM_REAL), dimension(nglob_DG) :: rho_DG_main, rhovx_DG_main, rhovz_DG_main, E_DG_main, e1_DG
  real(kind=CUSTOM_REAL), dimension(2, nglob_DG), intent(in) :: T_DG_main
  real(kind=CUSTOM_REAL), dimension(2, 2, nglob_DG), intent(in) :: V_DG_main
  real(kind=CUSTOM_REAL), dimension(nglob_DG), intent(out) :: dot_rho, dot_rhovx, dot_rhovz, dot_E, dot_e1
  
  ! Local variables.
  real(kind=CUSTOM_REAL), parameter :: ZERO = 0._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: ONE  = 1._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: TWO  = 2._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: HALF = 0.5_CUSTOM_REAL
  !real(kind=CUSTOM_REAL), parameter :: threshold = 0.000000001_CUSTOM_REAL ! Unused.
  real(kind=CUSTOM_REAL), dimension(nglob_DG) :: rho_DG, rhovx_DG, rhovz_DG, E_DG
  real(kind=CUSTOM_REAL), dimension(2, nglob_DG) :: T_DG
  real(kind=CUSTOM_REAL), dimension(2, 2, nglob_DG) :: V_DG
  integer :: ispec, i, j, k, iglob, iglob_unique
  real(kind=CUSTOM_REAL), dimension(NGLLX, NGLLZ) :: temp_rho_1, temp_rho_2, &
                                                     temp_rhovx_1, temp_rhovx_2, temp_rhovz_1, temp_rhovz_2, &
                                                     temp_E_1, temp_E_2, &
                                                     temp_nondiv_rho, &
                                                     temp_nondiv_rhovx, temp_nondiv_rhovz, temp_nondiv_E!, temp_e1
  real(kind=CUSTOM_REAL) :: xixl,xizl,gammaxl,gammazl,jacobianl ! Jacobian matrix and determinant.
  real(kind=CUSTOM_REAL) :: lambda, nx, nz, weight, &
                            temp_unknown_M, temp_unknown_P, temp_unknown2_M, temp_unknown2_P, &
                            temp_unknown, temp_unknown2, &
                            flux_n, jump, &
                            rho_DG_P, veloc_x_DG_P, veloc_z_DG_P, &
                            E_DG_P, p_DG_P, rhovx_DG_P, rhovz_DG_P, timelocal, &
                            Tx_DG_P, Tz_DG_P, Vxx_DG_P, Vzz_DG_P, Vxz_DG_P, Vzx_DG_P, T_P, &
                            ! TEST
                            gamma_P,&
                            e1_DG_P
  real(kind=CUSTOM_REAL) :: dT_dx, dT_dz
  !real(kind=CUSTOM_REAL) :: veloc_n_M, veloc_n_P
  integer :: iglobM, iglobP
  integer, dimension(3) :: neighbor
  real(kind=CUSTOM_REAL), dimension(nglob_DG) :: veloc_x_DG, veloc_z_DG, p_DG
  real(kind=CUSTOM_REAL) :: dux_dx, dux_dz, duz_dx, duz_dz ! Viscosity
  real(kind=CUSTOM_REAL) :: wxl, wzl ! Integration weigths.
  logical :: exact_interface_flux
  integer :: cnst_hdrsttc ! For computationnaly better CONSTRAIN_HYDROSTATIC switches.
  integer :: iface1, iface, iface1_neighbor, iface_neighbor, ispec_neighbor
  real(kind=CUSTOM_REAL) :: ya_x_l, ya_z_l ! Stretching absorbing boundary conditions.
  
  ! TESTS
  real(kind=CUSTOM_REAL) :: x,z ! Artifical advection and a priori damping.
  real(kind=CUSTOM_REAL) :: maxval_rho,maxval_rhovx,maxval_rhovz,maxval_E ! ABSORB
  logical :: ABSORB_BC ! ABSORB
  
  ! For more convinient CONSTRAIN_HYDROSTATIC switches.
  ! TODO: Replace the CONSTRAIN_HYDROSTATIC switches using this variable.
  if(CONSTRAIN_HYDROSTATIC) then
    cnst_hdrsttc=ONE
  else
    cnst_hdrsttc=ZERO
  endif
  
  rho_DG   = rho_DG_main
  rhovx_DG = rhovx_DG_main
  rhovz_DG = rhovz_DG_main
  E_DG     = E_DG_main
  
  ! TODO: remove?
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
  
  do ispec = 1, nspec ! Loop over elements.
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
          
          if(ABC_STRETCH .and. stretching_buffer(ibool_before_perio(i, j, ispec))>0) then
          !if(ABC_STRETCH) then
            ! Here are updated the operator \nabla and the jacobian, but only if stretching is activated and we are in a buffer.
            ! \partial_x becomes \ya_x\partial_x, and since \partial_x=(\partial_x\xi)\partial_\xi+(\partial_x\eta)\partial_\eta, only updating \partial_x\xi and \partial_x\eta is enough. Idem for \partial_z. Hence, only updating xix to \ya_x * xix, xiz to \ya_z * xiz, etc. is enough to update the operator.
            ! The jacobian of the stretching transformation is updated following the same rationale.
            iglob_unique = ibool_before_perio(i, j, ispec);
            ya_x_l=stretching_ya(1, iglob_unique)
            ya_z_l=stretching_ya(2, iglob_unique)
            xixl = ya_x_l * xixl
            xizl = ya_z_l * xizl
            gammaxl = ya_x_l * gammaxl
            gammazl = ya_z_l * gammazl
            ! TODO: This is a bit rough, do it more clearly.
            jacobianl = ya_x_l*ya_z_l*jacobianl
          endif
          
          wzl = wzgll(j)
          wxl = wxgll(i)
          
          ! Inviscid stress tensor's contributions.
          ! Mass conservation.
          temp_unknown = rhovx_DG(iglob)
          temp_unknown2 = rhovz_DG(iglob)
          temp_rho_1(i, j) = wzl * jacobianl * (xixl * temp_unknown + xizl * temp_unknown2)
          temp_rho_2(i, j) = wxl * jacobianl * (gammaxl * temp_unknown + gammazl * temp_unknown2)
          
          ! x-Momentum.
          if(.not. CONSTRAIN_HYDROSTATIC) then
            temp_unknown = rho_DG(iglob)*veloc_x_DG(iglob)**2 + p_DG(iglob)
          else
            temp_unknown = rho_DG(iglob)*veloc_x_DG(iglob)**2 + (p_DG(iglob) - p_DG_init(iglob))
          endif
          temp_unknown2 = rho_DG(iglob)*veloc_x_DG(iglob)*veloc_z_DG(iglob)
          temp_rhovx_1(i, j) = wzl * jacobianl * (xixl * temp_unknown + xizl * temp_unknown2) 
          temp_rhovx_2(i, j) = wxl * jacobianl * (gammaxl * temp_unknown + gammazl * temp_unknown2) 
          
          ! z-Momentum.
          temp_unknown = rho_DG(iglob)*veloc_x_DG(iglob)*veloc_z_DG(iglob)
          if(.not. CONSTRAIN_HYDROSTATIC) then
            temp_unknown2 = rho_DG(iglob)*veloc_z_DG(iglob)**2 + p_DG(iglob)
          else
            temp_unknown2 = rho_DG(iglob)*veloc_z_DG(iglob)**2 + (p_DG(iglob) - p_DG_init(iglob))
          endif
          temp_rhovz_1(i, j) = wzl * jacobianl * (xixl * temp_unknown + xizl * temp_unknown2) 
          temp_rhovz_2(i, j) = wxl * jacobianl * (gammaxl * temp_unknown + gammazl * temp_unknown2) 
          
          ! Energy.
          if(.not. CONSTRAIN_HYDROSTATIC) then
            temp_unknown = veloc_x_DG(iglob)*(E_DG(iglob) + p_DG(iglob))
            temp_unknown2 = veloc_z_DG(iglob)*(E_DG(iglob) + p_DG(iglob))
          else          
            temp_unknown = veloc_x_DG(iglob)*(E_DG(iglob) + (p_DG(iglob) - p_DG_init(iglob)))
            temp_unknown2 = veloc_z_DG(iglob)*(E_DG(iglob) + (p_DG(iglob) - p_DG_init(iglob)))
          endif
          temp_E_1(i, j) = wzl * jacobianl * (xixl * temp_unknown + xizl * temp_unknown2) 
          temp_E_2(i, j) = wxl * jacobianl * (gammaxl * temp_unknown + gammazl * temp_unknown2) 
          
          ! Viscous stress tensor's contributions.
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
            ! When CONSTRAIN_HYDROSTATIC==.true., it is the first line of the perturbed tensor (\Sigma'_v).
            temp_unknown = muext(i, j, ispec)*TWO*dux_dx + (etaext(i, j, ispec) - (TWO/3.)*muext(i, j, ispec))*(dux_dx + duz_dz) 
            temp_unknown2 = muext(i, j, ispec)*( dux_dz + duz_dx )
            temp_rhovx_1(i, j) = temp_rhovx_1(i, j) - wzl * jacobianl * (xixl * temp_unknown + xizl * temp_unknown2) 
            temp_rhovx_2(i, j) = temp_rhovx_2(i, j) - wxl * jacobianl * (gammaxl * temp_unknown + gammazl * temp_unknown2) 
            
            ! Use the values stored in temp_unknown and temp_unknown2 to compute the x component of the first part of the viscous energy vector (\Sigma_v\cdot\vect{v}, or \Sigma'_v\cdot\vect{v} if CONSTRAIN_HYDROSTATIC==.true.).
            temp_unknown  = veloc_x_DG(iglob)*temp_unknown + veloc_z_DG(iglob)*temp_unknown2
            temp_E_1(i, j) = temp_E_1(i, j) - wzl * jacobianl * (xixl * temp_unknown) 
            temp_E_2(i, j) = temp_E_2(i, j) - wxl * jacobianl * (gammaxl * temp_unknown) 
            
            ! This vector, [temp_unknown, temp_unknown2], is the second line of the viscous Navier-Stokes tensor (\Sigma_v).
            ! When CONSTRAIN_HYDROSTATIC==.true., it is the first line of the perturbed tensor (\Sigma'_v).
            temp_unknown = muext(i, j, ispec)*( dux_dz + duz_dx )
            temp_unknown2 = muext(i, j, ispec)*TWO*duz_dz + (etaext(i, j, ispec) - (TWO/3.)*muext(i, j, ispec))*(dux_dx + duz_dz)
            temp_rhovz_1(i, j) = temp_rhovz_1(i, j) - wzl * jacobianl * (xixl * temp_unknown + xizl * temp_unknown2)
            temp_rhovz_2(i, j) = temp_rhovz_2(i, j) - wxl * jacobianl * (gammaxl * temp_unknown + gammazl * temp_unknown2)
            
            ! Use the values stored in temp_unknown and temp_unknown2 to compute the z component of the first part of the viscous energy vector (\Sigma_v\cdot\vect{v}, or \Sigma'_v\cdot\vect{v} if CONSTRAIN_HYDROSTATIC==.true.).
            temp_unknown2  = veloc_x_DG(iglob)*temp_unknown + veloc_z_DG(iglob)*temp_unknown2
            temp_E_1(i, j) = temp_E_1(i, j) - wzl * jacobianl * (xizl * temp_unknown2)
            temp_E_2(i, j) = temp_E_2(i, j) - wxl * jacobianl * (gammazl * temp_unknown2)
            
            ! TODO: When CONSTRAIN_HYDROSTATIC==.true., doesn't the energy contribution lack the term \Sigma_{v,0}{\vect{v}'} here? As follows:
            !if(.false. .and. CONSTRAIN_HYDROSTATIC) then
            !  ! DV0X_DZ should be set to $\partial_z{v_{0,x}}$ here.
            !  temp_unknown = muext(i, j, ispec)*DV0X_DZ*(veloc_z_DG(iglob)-rhovz_init(iglob)/rho_init(iglob))
            !  temp_unknown2 = muext(i, j, ispec)*DV0X_DZ*(veloc_x_DG(iglob)-rhovx_init(iglob)/rho_init(iglob))
            !  temp_E_1(i, j) = temp_E_1(i, j) + wzl*jacobianl*(xixl*temp_unknown + xizl*temp_unknown2) 
            !  temp_E_2(i, j) = temp_E_2(i, j) + wxl*jacobianl*(gammaxl*temp_unknown + gammazl*temp_unknown2) 
            !endif
            
            ! Add the heat contributions (second part of the viscous energy tensor, \kappa\nabla T, or \kappa\nabla{T'} if CONSTRAIN_HYDROSTATIC==.true.).
            temp_unknown  = kappa_DG(i, j, ispec)*dT_dx
            temp_unknown2 = kappa_DG(i, j, ispec)*dT_dz
            temp_E_1(i, j) = temp_E_1(i, j) - wzl * jacobianl * (xixl * temp_unknown + xizl * temp_unknown2) 
            temp_E_2(i, j) = temp_E_2(i, j) - wxl * jacobianl * (gammaxl * temp_unknown + gammazl * temp_unknown2) 
          endif
          
          ! Gravity (and in fact, everything not inside the divergence operator in the strong form, except sources) contributions (separated from the rest).
          temp_nondiv_rho(i, j) = ZERO
          temp_nondiv_rhovx(i, j) = -rho_DG(iglob)*potential_dphi_dx_DG(ibool(i, j, ispec))* jacobianl
          if(.not. CONSTRAIN_HYDROSTATIC) then
            temp_nondiv_rhovz(i, j) = -rho_DG(iglob)*potential_dphi_dz_DG(ibool(i, j, ispec))* jacobianl
          else
            temp_nondiv_rhovz(i, j) = -(rho_DG(iglob) - rho_init(iglob)) * potential_dphi_dz_DG(ibool(i, j, ispec)) * jacobianl 
          endif
          if(.not. CONSTRAIN_HYDROSTATIC) then
            temp_nondiv_E(i, j) = -rho_DG(iglob)*(veloc_x_DG(iglob)*potential_dphi_dx_DG(ibool(i, j, ispec)) + &
                                veloc_z_DG(iglob)*potential_dphi_dz_DG(ibool(i, j, ispec)))* jacobianl
          else
            ! By constraining hydrostaticity, remark that an additional term appears outside the divergence ($p_0\nabla\cdot{\vect{v}'}$).
            temp_nondiv_E(i, j) = &
                                -(rho_DG(iglob) - rho_init(iglob))*(veloc_x_DG(iglob)*potential_dphi_dx_DG(ibool(i, j, ispec)) + &
                                veloc_z_DG(iglob)*potential_dphi_dz_DG(ibool(i, j, ispec)))* jacobianl         
            temp_nondiv_E(i, j) = temp_nondiv_E(i, j) - p_DG_init(iglob)*(dux_dx + duz_dz)* jacobianl
          endif
          ! Memory variable evolution.
          temp_nondiv_E(i, j) = temp_nondiv_E(i, j) - jacobianl * (p_DG_init(iglob)*gammaext_DG(iglob)) &
                              * ( (tau_epsilon(i, j, ispec)/tau_sigma(i, j, ispec)) - ONE ) &
                              * ( dux_dx + duz_dz - e1_DG(iglob))/(gammaext_DG(iglob) - ONE)
          
          ! TEST ARTIFICIAL ADVECTION on top buffer
          if(.false.) then
            x=coord(1, ibool_before_perio(i, j, ispec))
            z=coord(2, ibool_before_perio(i, j, ispec))
            if(     (ABC_STRETCH_LEFT   .and. x < mesh_xmin + ABC_STRETCH_LEFT_LBUF) & ! left stretching and in left buffer zone
               .or. (ABC_STRETCH_RIGHT  .and. x > mesh_xmax - ABC_STRETCH_RIGHT_LBUF) & ! right stretching and in right buffer zone
               .or. (ABC_STRETCH_BOTTOM .and. z < mesh_zmin + ABC_STRETCH_BOTTOM_LBUF) & ! bottom stretching and in bottom buffer zone
               .or. (ABC_STRETCH_TOP    .and. z > mesh_zmax - ABC_STRETCH_TOP_LBUF)) then ! top stretching and in top buffer zone
              !temp_nondiv_rho(i, j) = temp_nondiv_rho(i, j) + ((z-mesh_zmax)/ABC_STRETCH_LBUF+ONE)*10.*rho_DG(iglob)
              !temp_nondiv_rhovx(i, j) = temp_nondiv_rhovx(i, j) + ((z-mesh_zmax)/ABC_STRETCH_LBUF+ONE)*10.*rhovx_DG(iglob)
              !temp_nondiv_rhovz(i, j) = temp_nondiv_rhovz(i, j) + ((z-mesh_zmax)/ABC_STRETCH_LBUF+ONE)*10.*rhovz_DG(iglob)
              !temp_nondiv_E(i, j) = temp_nondiv_E(i, j) + ((z-mesh_zmax)/ABC_STRETCH_LBUF+ONE)*10.*E_DG(iglob)
              temp_nondiv_rho(i, j) = temp_nondiv_rho(i, j) + 10.*rho_DG(iglob)
              temp_nondiv_rhovx(i, j) = temp_nondiv_rhovx(i, j) + 10.*rhovx_DG(iglob)
              temp_nondiv_rhovz(i, j) = temp_nondiv_rhovz(i, j) + 10.*rhovz_DG(iglob)
              temp_nondiv_E(i, j) = temp_nondiv_E(i, j) + 10.*E_DG(iglob)
            endif
          endif
          ! TEST PRIORI DAMPING on top buffer
          if(.false.) then
            x=coord(1, ibool_before_perio(i, j, ispec))
            z=coord(2, ibool_before_perio(i, j, ispec))
            if(     (ABC_STRETCH_LEFT   .and. x < mesh_xmin + ABC_STRETCH_LEFT_LBUF) & ! left stretching and in left buffer zone
               .or. (ABC_STRETCH_RIGHT  .and. x > mesh_xmax - ABC_STRETCH_RIGHT_LBUF) & ! right stretching and in right buffer zone
               .or. (ABC_STRETCH_BOTTOM .and. z < mesh_zmin + ABC_STRETCH_BOTTOM_LBUF) & ! bottom stretching and in bottom buffer zone
               .or. (ABC_STRETCH_TOP    .and. z > mesh_zmax - ABC_STRETCH_TOP_LBUF)) then ! top stretching and in top buffer zone
              call boundary_condition_DG(i, j, ispec, timelocal, rho_DG_P, rhovx_DG_P, rhovz_DG_P, E_DG_P, &
                      veloc_x_DG_P, veloc_z_DG_P, p_DG_P, e1_DG_P)
              
              if(ABC_STRETCH_TOP) then
                if((z-mesh_zmax)/ABC_STRETCH_TOP_LBUF+ONE>ZERO .and. (z-mesh_zmax)/ABC_STRETCH_TOP_LBUF+ONE<=ONE) then
                  temp_nondiv_rho(i, j) = temp_nondiv_rho(i, j)&
                                         + ((z-mesh_zmax)/ABC_STRETCH_TOP_LBUF+ONE)**2.*10.*( rho_DG(iglob) - rho_DG_P)
                  temp_nondiv_rhovx(i, j) = temp_nondiv_rhovx(i, j)&
                                           + ((z-mesh_zmax)/ABC_STRETCH_TOP_LBUF+ONE)**2.*10.*( rhovx_DG(iglob) - rhovx_DG_P)
                  temp_nondiv_rhovz(i, j) = temp_nondiv_rhovz(i, j)&
                                           + ((z-mesh_zmax)/ABC_STRETCH_TOP_LBUF+ONE)**2.*10.*( rhovz_DG(iglob) - rhovz_DG_P)
                  temp_nondiv_E(i, j) = temp_nondiv_E(i, j)&
                                       + ((z-mesh_zmax)/ABC_STRETCH_TOP_LBUF+ONE)**2.*10.*( E_DG(iglob) - E_DG_P)
                endif
              endif
            endif
          endif
          ! END OF TESTS
          
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
          
          ! Terms outside divergence operator.
          dot_rho(iglob)   = dot_rho(iglob)   + temp_nondiv_rho(i, j) * wxl * wzl
          dot_rhovx(iglob) = dot_rhovx(iglob) + temp_nondiv_rhovx(i, j) * wxl * wzl
          dot_rhovz(iglob) = dot_rhovz(iglob) + temp_nondiv_rhovz(i, j) * wxl * wzl
          dot_E(iglob)     = dot_E(iglob)     + temp_nondiv_E(i, j) * wxl * wzl
          
          !dot_E(iglob) = 0. ! DEBUG: DEACTIVATE ENERGY INTERNAL (ALL EXCEPT FLUXES) EVOLUTION.
        enddo
      enddo
      
      !DEBUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUG
      !DEBUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUG
      !DEBUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUG
      if(.false. .and. timelocal>=2e-5 .and. timelocal<8e-5) then
        do j = 1, NGLLZ
          do i = 1, NGLLX
            if(coord(1, ibool(i, j, ispec)) > 59. &
               .and. coord(1, ibool(i, j, ispec)) < 60. &
               .and. coord(2, ibool(i, j, ispec)) > 59.5) then ! DEBUG
              iglob = ibool_DG(i, j, ispec)
              write(*,*) timelocal, coord(1, ibool(i, j, ispec)), coord(2, ibool(i, j, ispec)), veloc_z_DG(iglob)
            endif
          enddo
        enddo
      endif
      
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
          
          ! Step 2: knowing the normals' parameters, compute now the fluxes.
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
          lambda = ZERO
          jump   = ZERO
          !veloc_n_M = abs(veloc_x_DG(iglobM)*nx + veloc_z_DG(iglobM)*nz)
          !veloc_n_P = abs(veloc_x_DG_P*nx + veloc_z_DG_P*nz)
          ! Save some operations.
          !lambda = max( veloc_n_M + sqrt(abs(gammaext_DG(iglobM)*p_DG(iglobM)/rho_DG(iglobM))), &
          !              veloc_n_P + sqrt(abs(gamma_P*p_DG_P/rho_DG_P)) )
          lambda = max(abs(veloc_x_DG(iglobM)*nx + veloc_z_DG(iglobM)*nz) &
                       + sqrt(abs(gammaext_DG(iglobM)*p_DG(iglobM)/rho_DG(iglobM))), &
                       abs(veloc_x_DG_P*nx + veloc_z_DG_P*nz) &
                       + sqrt(abs(gamma_P*p_DG_P/rho_DG_P)))
          
          if(ABC_STRETCH .and. stretching_buffer(ibool_before_perio(i, j, ispec))>0) then
            ! Update flux with stretching components. It is quite ugly to implement stretching like this (since stretching has nothing to do with the normals), but at least it is quick and does the job. I am sorry.
            ! TODO: Do it more clearly.
            iglob_unique=ibool_before_perio(i, j, ispec)
            !weight=weight*(stretching_ya(1,iglob_unique)*stretching_ya(2,iglob_unique))
            nx=stretching_ya(1,iglob_unique)*stretching_ya(2,iglob_unique)*nx
            nz=stretching_ya(1,iglob_unique)*stretching_ya(2,iglob_unique)*nz
            !lambda = max(abs(veloc_x_DG(iglobM)*nx + veloc_z_DG(iglobM)*nz) &
            !             *(stretching_ya(1,iglob_unique)*stretching_ya(2,iglob_unique)) &
            !             + sqrt(abs(gammaext_DG(iglobM)*p_DG(iglobM)/rho_DG(iglobM))), &
            !             abs(veloc_x_DG_P*nx + veloc_z_DG_P*nz) &
            !             *(stretching_ya(1,iglob_unique)*stretching_ya(2,iglob_unique)) &
            !             + sqrt(abs(gamma_P*p_DG_P/rho_DG_P)))!&
            !            *(stretching_ya(1,iglob_unique)*stretching_ya(2,iglob_unique))
          endif
          
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
          ! 1st line of \Sigma_v, or \Sigma'_v if CONSTRAIN_HYDROSTATIC==.true..
          ! The vector [temp_unknown, temp_unknown2] represents the mean average flux at the boundary of the x-momentum.
          ! Recall: dux_dx, duz_dx, dux_dz, and duz_dz already contain the 0.5 factor to put the flux under mean average form.
          temp_unknown = muext(i, j, ispec)*TWO*dux_dx + (etaext(i, j, ispec) - (TWO/3.)*muext(i, j, ispec))*(dux_dx + duz_dz) 
          temp_unknown2 = muext(i, j, ispec)*( dux_dz + duz_dx )
          ! Dot product.
          flux_n = temp_unknown*nx + temp_unknown2*nz ! [3 operations + 1 affectation], instead of [3 operations + 3 affectations].
          dot_rhovx(iglobM) = dot_rhovx(iglobM) + weight*flux_n
          ! The computed values contained in the variables temp_unknown and temp_unknown2 can be used to compute the energy's x component of the mean average flux at the boundary. Thus, we add this contribution here. (\Sigma_v\cdot\vect{v}, or \Sigma'_v\cdot\vect{v} if CONSTRAIN_HYDROSTATIC==.true..)
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
          ! 2nd line of \Sigma_v, or \Sigma'_v if CONSTRAIN_HYDROSTATIC==.true..
          ! The vector [temp_unknown, temp_unknown2] represents the mean average flux at the boundary of the z-momentum.
          ! Recall: dux_dx, duz_dx, dux_dz, and duz_dz already contain the 0.5 factor to put the flux under mean average form.
          temp_unknown = muext(i, j, ispec)*( dux_dz + duz_dx )
          temp_unknown2 = muext(i, j, ispec)*TWO*duz_dz + (etaext(i, j, ispec) - (TWO/3.)*muext(i, j, ispec))*(dux_dx + duz_dz)
          ! Dot product.
          flux_n = temp_unknown*nx + temp_unknown2*nz ! [3 operations + 1 affectation], instead of [3 operations + 3 affectations].
          dot_rhovz(iglobM) = dot_rhovz(iglobM) + weight*flux_n
          ! The computed values contained in the variables temp_unknown and temp_unknown2 can be used to compute the energy's z component of the mean average flux at the boundary. Thus, we add this contribution here. (\Sigma_v\cdot\vect{v}, or \Sigma'_v\cdot\vect{v} if CONSTRAIN_HYDROSTATIC==.true..)
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
          
          ! TODO: When CONSTRAIN_HYDROSTATIC==.true., doesn't the energy contribution lack the term \Sigma_{v,0}{\vect{v}'} here? As follows:
          !if(.false. .and. CONSTRAIN_HYDROSTATIC) then
          !  ! DV0X_DZ should be set to $\partial_z{v_{0,x}}$ here.
          !  temp_unknown = muext(i, j, ispec)*DV0X_DZ*(veloc_z_DG(iglob)-rhovz_init(iglob)/rho_init(iglob))
          !  temp_unknown2 = muext(i, j, ispec)*DV0X_DZ*(veloc_x_DG(iglob)-rhovx_init(iglob)/rho_init(iglob))
          !  dot_E(iglobM) = dot_E(iglobM)+weight*(temp_unknown*nx+temp_unknown2*nz)
          !endif
          
          ! Energy equation's heat flux' contribution (last remaining term, viscous).
          ! Recall: dT_dx already contains the 0.5 factor to put the flux under mean average form.
          dot_E(iglobM) = dot_E(iglobM) &
                          + weight*( kappa_DG(i, j, ispec)*( dT_dx*nx + dT_dz*nz ) )
          
          !dot_E(iglobM) = 0. ! DEBUG: DEACTIVATE ENERGY FLUXES.
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
                         ibool_before_perio,&
                         rhovx_init, rhovz_init, rho_init, T_init, CONSTRAIN_HYDROSTATIC, &
                         link_iface_ijispec, nx_iface, nz_iface, weight_iface, neighbor_DG_iface,&
                         ABC_STRETCH,stretching_ya,stretching_buffer ! Stretching-based absorbing conditions.
                         
  implicit none
  
  ! local parameters
  integer :: ispec,i, j,k, iglob, iglobM, iglobP, iglob_unique
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
  real(kind=CUSTOM_REAL) :: ya_x_l, ya_z_l
  
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
          
          if(ABC_STRETCH .and. stretching_buffer(ibool_before_perio(i, j, ispec))>0) then
            ! See beginning of subroutine compute_forces_acoustic_DG for detailed explanations.
            iglob_unique=ibool_before_perio(i, j, ispec)
            ya_x_l=stretching_ya(1, iglob_unique)
            ya_z_l=stretching_ya(2, iglob_unique)
            xixl = ya_x_l * xixl
            xizl = ya_z_l * xizl
            gammaxl = ya_x_l * gammaxl
            gammazl = ya_z_l * gammazl
            ! TODO: Do that more clearly.
            jacobianl = ya_x_l*ya_z_l*jacobianl
          endif
          
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
          
            if(ABC_STRETCH .and. stretching_buffer(ibool_before_perio(i, j, ispec))>0) then
              ! Update flux with stretching components. See explanation in the surface terms part in the subroutine above.
              ! TODO: Do that more clearly.
              iglob_unique=ibool_before_perio(i, j, ispec)
              weight=stretching_ya(1,iglob_unique)*stretching_ya(2,iglob_unique)*weight
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
end subroutine compute_viscous_tensors

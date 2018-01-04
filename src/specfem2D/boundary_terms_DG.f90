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
! initial_condition_DG                                         !
! ------------------------------------------------------------ !
! Computes initial conditions.

  subroutine initial_condition_DG()

  ! compute forces in the acoustic elements in forward simulation and in adjoint simulation in adjoint inversion

  use constants,only: CUSTOM_REAL,NGLLX,NGLLZ,PI

  use specfem_par, only: ibool_DG, &
        rho_DG, rhovx_DG, rhovz_DG, E_DG, e1_DG, nspec

  implicit none
  
  ! Parameters.
  real(kind=CUSTOM_REAL), parameter :: ZEROl = 0._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: ONEl  = 1._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: TWOl  = 2._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: SIXl  = 6._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: HALFl = 0.5_CUSTOM_REAL
  
  integer :: ispec, iglob, i, j
  
  real(kind=CUSTOM_REAL) :: rho_DG_P, rhovx_DG_P, rhovz_DG_P, E_DG_P, &
      veloc_x_DG_P, veloc_z_DG_P, p_DG_P, e1_DG_P
  
    
  do ispec = 1, nspec ! Loop over elements.
    ! Double loop over GLL points. ??[to compute and store gradients]
    do j = 1, NGLLZ
      do i = 1, NGLLX
       iglob = ibool_DG(i, j, ispec)
       call boundary_condition_DG(i, j, ispec, ZEROl, rho_DG_P, rhovx_DG_P, rhovz_DG_P, E_DG_P, &
            veloc_x_DG_P, veloc_z_DG_P, p_DG_P, e1_DG_P)

       rho_DG(iglob)   = rho_DG_P
       rhovx_DG(iglob) = rhovx_DG_P
       rhovz_DG(iglob) = rhovz_DG_P
       E_DG(iglob)     = E_DG_P
       e1_DG(iglob)    = e1_DG_P
      enddo
    enddo
  enddo  
  
  if(.false.) then
    ! TODO: Why 'if(.false.)'?
    call recompute_density()
  endif
  
  end subroutine initial_condition_DG
  
! ------------------------------------------------------------ !
! recompute_density                                            !
! ------------------------------------------------------------ !
  
  subroutine recompute_density()

  use constants,only: CUSTOM_REAL, NGLLX, NGLLZ

  use specfem_par, only: xix, xiz, &
        gammaext_DG,gammax, gammaz,ibool_DG, gravityext, &
        E_DG, rhovx_DG, rhovz_DG, rho_DG, pext_DG, windxext, nspec, rhoext, &
        nglob_DG, hprime_xx, hprime_zz , ispec_is_acoustic 
        
  implicit none
  
  ! local parameters
  integer :: ispec,i,j,k,iglob

  ! Jacobian matrix and determinant
  real(kind=CUSTOM_REAL) :: xixl,xizl,gammaxl,gammazl
  
  ! Viscosity
  real(kind=CUSTOM_REAL) :: duz_dxi, duz_dgamma, duz_dz
  
  ! Parameters
  real(kind=CUSTOM_REAL), parameter :: ZERO = 0._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: ONE  = 1._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: TWO  = 2._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: HALF = 0.5_CUSTOM_REAL
  
  real(kind=CUSTOM_REAL), dimension(nglob_DG) :: p_DG

  p_DG = (gammaext_DG - ONE)*( E_DG &
        - (HALF/rho_DG)*( rhovx_DG**2 + rhovz_DG**2 ) )
  
  do ispec = 1, nspec ! Loop over elements.
    ! Acoustic spectral element.
    if (ispec_is_acoustic(ispec)) then
    !if (ispec_is_acoustic_DG(ispec)) then
      do j = 1, NGLLZ
        do i = 1, NGLLX
          iglob = ibool_DG(i, j, ispec)
          duz_dxi    = ZERO
          duz_dgamma = ZERO
          xixl = xix(i, j, ispec)
          xizl = xiz(i, j, ispec)
          gammaxl = gammax(i, j, ispec)
          gammazl = gammaz(i, j, ispec)
          
          ! first double loop over GLL points to compute and store gradients
          ! we can merge the two loops because NGLLX == NGLLZ
          do k = 1, NGLLX
            duz_dxi    = duz_dxi    + log(p_DG(ibool_DG(k,j,ispec))) * real(hprime_xx(i,k), kind=CUSTOM_REAL)
            duz_dgamma = duz_dgamma + log(p_DG(ibool_DG(i,k,ispec))) * real(hprime_zz(j,k), kind=CUSTOM_REAL)
          enddo
          ! derivatives of velocities
          duz_dz = ( duz_dxi * xizl + duz_dgamma * gammazl )
          rhoext(i, j, ispec) = -duz_dz*p_DG(iglob)/gravityext(i, j, ispec)
          rho_DG(iglob)     = rhoext(i, j, ispec)
          E_DG(iglob)         = pext_DG(i, j, ispec)/(gammaext_DG(iglob) - 1.) &
                                + 0.5*rho_DG(iglob)*( windxext(i, j, ispec)**2 )
          rhovx_DG(iglob)    = rho_DG(iglob)*windxext(i, j, ispec)
         enddo
       enddo
     endif
   enddo
  
  end subroutine recompute_density
  
! ------------------------------------------------------------ !
! boundary_condition_DG                                        !
! ------------------------------------------------------------ !
! Sets all variables of the model at initial time. Those are
! also used as the far-field model. This subroutine also takes
! care of the various forcings.

  subroutine boundary_condition_DG(i, j, ispec, timelocal, rho_DG_P, rhovx_DG_P, rhovz_DG_P, E_DG_P, &
      veloc_x_DG_P, veloc_z_DG_P, p_DG_P, e1_DG_P)

  use constants,only: CUSTOM_REAL, gamma_euler, PI

  use specfem_par, only: ibool_before_perio, ibool_DG, coord, MODEL, &
        rhoext, windxext, pext_DG, gravityext, gammaext_DG, &
        etaext, muext, coord_interface, kappa_DG, cp, cnu, &
        tau_epsilon, tau_sigma, &
        gravity_cte_DG, dynamic_viscosity_cte_DG, thermal_conductivity_cte_DG, tau_eps_cte_DG, tau_sig_cte_DG, SCALE_HEIGHT, &
        USE_ISOTHERMAL_MODEL, potential_dphi_dx_DG, potential_dphi_dz_DG, ibool, &
        surface_density, sound_velocity, wind, TYPE_FORCING, &
        forcing_initial_time, main_time_period, forcing_initial_loc, main_spatial_period

  implicit none
  
  integer :: i, j, ispec
  
  real(kind=CUSTOM_REAL) :: rho_DG_P, rhovx_DG_P, rhovz_DG_P, E_DG_P, &
      veloc_x_DG_P, veloc_z_DG_P, p_DG_P, e1_DG_P, timelocal
  
  ! Parameters
  real(kind=CUSTOM_REAL), parameter :: ZEROl = 0._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: ONEl  = 1._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: TWOl  = 2._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: SIXl  = 6._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: HALFl = 0.5_CUSTOM_REAL
  
  real(kind=CUSTOM_REAL) :: x, z, H, G!, A
  !real(kind=CUSTOM_REAL) :: x0, z0, r, beta
  ! Forcing
  real(kind=CUSTOM_REAL) :: to, perio, lambdo, xo
  
  real(kind=CUSTOM_REAL) :: mu_visco, eta_visco
  ! Hydrostatic solution
  !real(kind=CUSTOM_REAL) :: RR, p0, rho0
  ! Density current
  !real(kind=CUSTOM_REAL) :: cp, cnu, exner, RR, p0, rho0, rs, theta, theta0
  ! Linear mountains
  !real(kind=CUSTOM_REAL) :: cp, cnu, exner, RR, p0, rho0, rs, theta, theta0, Nsq
  real(kind=CUSTOM_REAL) :: VELOC_TSUNAMI
  ! Thermal bubble
  real(kind=CUSTOM_REAL) :: Tl, Tu, rho0, p0, RR
  !real(kind=CUSTOM_REAL) :: Tl, Tu, p0, rho0, RR, rs!, theta, thetaprime, pibar
  !real(kind=CUSTOM_REAL) :: Tl, Tu, p0, RR, rs, theta0, Nsq, rho0
  !real(kind=CUSTOM_REAL) :: cp, cnu, R, P0, N, theta0, exner, theta
  !real(kind=CUSTOM_REAL) :: x0, z0, r, rs, p0, RR, cnu, cp, theta0, exner
  !real(kind=CUSTOM_REAL) :: x0, z0, r, theta, RR, p0!, exner, theta0
  !real(kind=CUSTOM_REAL) :: &
  !      lz, vs, vs_x, vs_z, Ms, p_1, rho_1, veloc_x_1, veloc_z_1, p_2, rho_2, veloc_x_2, veloc_z_2 
  
  ! Coordinate of the solid-fluid interface
  !coord_interface = 304.!250.
  
  x = real(coord(1, ibool_before_perio(i, j, ispec)), kind=CUSTOM_REAL)
  z = real(coord(2, ibool_before_perio(i, j, ispec)), kind=CUSTOM_REAL)
  z = z - coord_interface
  
  if(trim(MODEL) == 'external') then
    ! If an external model data file is given for initial conditions, use it.
    rho_DG_P     = rhoext(i, j, ispec)
    p_DG_P       = pext_DG(i, j, ispec)
    !p_DG_P       = rho_DG_P * gravityext(i, j, ispec) * Htabext_DG(ibool_DG(i, j, ispec)) ! Hydrostaticity: impose p = \rho * g * H, where the RHS is gotten from the external model data file.
    veloc_x_DG_P = windxext(i, j, ispec) ! Read horizontal wind.
    veloc_z_DG_P = ZEROl ! Impose vertical wind to zero.

    if(timelocal == 0) then
      !gravityext(i, j, ispec) = G!9.831
      !gammaext_DG(ibool_DG(i, j, ispec))   = 1.4!(rho_DG_P/p_DG_P)*(vpext(i, j, ispec)**2)!gamma_euler!
      !muext(i, j, ispec)  = 1.0d-05
      !etaext(i, j, ispec)  = (4/3)*muext(i, j, ispec)
      !kappa_DG(i, j, ispec) = 0.025!kappa_DG(i, j, ispec) / 10.!0.2
      !cp = 7/2
      !cnu = 5/2

      !tau_epsilon(i, j, ispec) = 1.
      !tau_sigma(i, j, ispec)   = 1.!0.4013

      potential_dphi_dx_DG(ibool(i, j, ispec)) = ZEROl
      potential_dphi_dz_DG(ibool(i, j, ispec)) = gravityext(i, j, ispec)
    endif
  
  else
    ! If no external model data file is given (no initial conditions were specified), build model.
    ! Set gravity, viscosity coefficients, relaxation times, and gravity potential (Phi) derivatives.
    if(timelocal == 0) then
      ! > Isobaric case. Since we need to stay hydrostatic, the gravity field needs to stay 0.
      gravityext(i, j, ispec) = ZEROl !real(coord(1,ibool_before_perio(i, j, ispec)), kind=CUSTOM_REAL)
      ! > Isothermal case.
      if(USE_ISOTHERMAL_MODEL) then
        gravityext(i, j, ispec) = real(gravity_cte_DG, kind=CUSTOM_REAL)
      endif
      !WRITE(*,*) ">>>", gravity_cte_DG, gravityext(i, j, ispec), abs(real(gravity_cte_DG, kind=CUSTOM_REAL)- gravityext(i, j, ispec)) ! DEBUG
      !cp = 29 !7!1010!7/2 ! DEBUG
      !cnu = 20.8 !5!718!5/2 ! DEBUG
      gammaext_DG(ibool_DG(i, j, ispec)) = cp/cnu!1.33076167! 1.4!gamma_euler
      mu_visco  = dynamic_viscosity_cte_DG!1.092656e-05
      eta_visco = (4.0d0/3.0d0)*dynamic_viscosity_cte_DG!2.67e-05
      muext(i, j, ispec) = mu_visco
      etaext(i, j, ispec) = eta_visco
      !kappa_DG(i, j, ispec) = 0.0
      kappa_DG(i, j, ispec) = thermal_conductivity_cte_DG!0.!0.05!4.79046750E-04
      tau_epsilon(i, j, ispec) = tau_eps_cte_DG!1 !1.5!2.!1.5
      tau_sigma(i, j, ispec)   = tau_sig_cte_DG!1 !1.!1./(8*PI**2)!0.025!1
      
      potential_dphi_dx_DG(ibool(i, j, ispec)) = ZEROl
      potential_dphi_dz_DG(ibool(i, j, ispec)) = gravityext(i, j, ispec)
    endif
    
    ! Set auxiliary variable.
    e1_DG_P = ZEROl
    
    ! Set density.
    rho_DG_P = surface_density
    if(USE_ISOTHERMAL_MODEL) then
      H = SCALE_HEIGHT!9350!10000
      G = gravityext(i, j, ispec)
      rho_DG_P = surface_density*exp(-z/H)
    endif
    
    ! Set pressure.
    ! > Isobaric case.
    p_DG_P = (sound_velocity**2)*rho_DG_P/gammaext_DG(ibool_DG(i, j, ispec)) ! Acoustic only (under ideal gas hypothesis): p = c^2 * \rho / \gamma.
    ! > Isothermal case.
    if(USE_ISOTHERMAL_MODEL) then
      ! If we use an isothermal model, replace the value just computed by the one under hydrostatic hypothesis: p = \rho * g * H.
      p_DG_P = rho_DG_P*G*H
    endif
    
    ! Set wind.
    veloc_x_DG_P = wind ! Read horizontal wind, from specfem2_par. TODO: why not windxext(i, j, ispec) as above?
    veloc_z_DG_P = ZEROl ! Impose vertical wind to zero.
  endif
  
  ! --------------------------- !
  ! Rayleigh-Taylor instability !
  ! (case 5.3 paper).           !
  ! --------------------------- !
  ! TODO: do something instead of this 'if(.false.)'.
  if(.false.) then
    Tl = 2.!TWOl
    Tu = 1.
    p0   = ONEl
    rho0 = ONEl
    RR   = 1.!8.3145_CUSTOM_REAL
    gammaext_DG(ibool_DG(i, j, ispec)) = 1.4
    gravityext(i, j, ispec) = 1
    !if(z <= 1._CUSTOM_REAL) &
    p_DG_P = p0*exp(-(z-ONEl/(TWOl))/(RR*Tl))
    if(z > ONEl/TWOl) &
    p_DG_P = p0*exp(-(z-ONEl/(TWOl))/(RR*Tu))
    !!if(z <= 1._CUSTOM_REAL) &
    rho_DG_P = p_DG_p/(RR*Tl)
    if(z > ONEl/TWOl) &
    rho_DG_P = p_DG_p/(RR*Tu)
    veloc_x_DG_P = ZEROl
    veloc_z_DG_P = ZEROl
  endif
       
  !call boundary_forcing_DG(timelocal, rho_DG_P, veloc_x_DG_P, veloc_z_DG_P, E_DG_P)

  ! Set (overwrite) velocities according to forcing.
 
  ! Forcing parameters from par_file
  lambdo = main_spatial_period
  xo     = forcing_initial_loc
  perio  = main_time_period
  to     = forcing_initial_time

  ! --------------------------- !
  ! Acoustic plane wave         !
  ! forcing.                    !
  ! --------------------------- !
  ! TODO: read the parameters from parfile.
  if(z == ZEROl .AND. TYPE_FORCING == 1) then
    veloc_z_DG_P = 0.01*(&
                   - (2d0/(perio/4d0))*((timelocal-(to-perio/4d0))/(perio/4d0))* &
                         (exp(-((timelocal-(to-perio/4d0))/(perio/4d0))**2)) &
                   + (2d0/(perio/4d0))*((timelocal-(to+perio/4d0))/(perio/4d0))* &
                         (exp(-((timelocal-(to+perio/4d0))/(perio/4d0))**2)) )
  endif

  ! --------------------------- !
  ! Gravity wave forcing.       !
  ! --------------------------- !
  ! TODO: read the parameters from parfile.
  if(z == ZEROl .AND. TYPE_FORCING == 2) &     
  veloc_z_DG_P = 0.001*(&
                - (2d0/(perio/4d0))*((timelocal-(to-perio/4d0))/(perio/4d0))* &
                         (exp(-((timelocal-(to-perio/4d0))/(perio/4d0))**2)) &
                + (2d0/(perio/4d0))*((timelocal-(to+perio/4d0))/(perio/4d0))* &
                         (exp(-((timelocal-(to+perio/4d0))/(perio/4d0))**2)) ) &
                 * ( exp(-((x-(xo-lambdo/4))/(lambdo/4))**2) - &
                        exp(-((x-(xo+lambdo/4))/(lambdo/4))**2) ) 
       
  ! --------------------------- !
  ! Tsunami forcing.            !
  ! --------------------------- !
  ! TODO: do something instead of this 'if(.false.)'.
  if(.false.) then
    ! TODO: read the parameters from parfile.
    VELOC_TSUNAMI = 200
    perio = 1000
    lambdo = 50000
    if(timelocal < perio) then
        veloc_z_DG_P = (1d0/perio)*exp( -((x-lambdo)/(sqrt(2d0)*perio))**2 )
    else
        veloc_z_DG_P = &
                2d0*((VELOC_TSUNAMI)/(sqrt(2d0)*perio))*(((x-lambdo)-VELOC_TSUNAMI*(timelocal - perio))/(sqrt(2d0)*perio))&
                *exp( -(((x-lambdo)-VELOC_TSUNAMI*(timelocal - perio))/(sqrt(2d0)*perio))**2 )
    endif
  endif

  ! Set energy.
  E_DG_P = p_DG_P/(gammaext_DG(ibool_DG(i, j, ispec)) - ONEl) &
                      + rho_DG_P*HALFl*( veloc_x_DG_P**2 + veloc_z_DG_P**2 ) 

  ! Set momenta.
  rhovx_DG_P = rho_DG_P*veloc_x_DG_P
  rhovz_DG_P = rho_DG_P*veloc_z_DG_P
   
  end subroutine boundary_condition_DG

! ------------------------------------------------------------ !
! compute_interface_unknowns                                   !
! ------------------------------------------------------------ !
! From the coordinates of a GLL point in an element (local coordinates (i, j) and global element number ispec), and its neighbour's identifier (neighbor), compute the values of the constitutive variables at the neighbour.
! Variables ending in "_P" (for "plus") are output-intended and are the sought values, or exterior values.
! Variables ending in "_iM" (for "minus") are input-intended and should correspond to the interior values.
! Variables ending in "_iP" are input-intended, and correspond to the exterior values, if known. Remark that those variables are only used if neighbor(3) != -1.
 
  subroutine compute_interface_unknowns(i, j, ispec, rho_DG_P, rhovx_DG_P, &
                rhovz_DG_P, E_DG_P, veloc_x_DG_P, veloc_z_DG_P, p_DG_P, T_P, &
                Tx_DG_P, Tz_DG_P, Vxx_DG_P, Vzz_DG_P, Vzx_DG_P, Vxz_DG_P, gamma_P, &
                neighbor, &
                exact_interface_flux, &
                rho_DG_iM, E_DG_iM, rhovx_DG_iM, rhovz_DG_iM, &
                V_DG_iM, T_DG_iM, &
                rho_DG_iP, E_DG_iP, rhovx_DG_iP, rhovz_DG_iP, &
                V_DG_iP, T_DG_iP, &
                nx, nz, weight, timelocal, &
                iface1, iface)
  
! compute forces in the acoustic elements in forward simulation and in adjoint simulation in adjoint inversion
  use constants,only: CUSTOM_REAL,NGLLX,NGLLZ,gamma_euler

  use specfem_par,only:  ibool_DG, &
                         ispec_is_acoustic_forcing, &
                         ACOUSTIC_FORCING, &
                          ispec_is_acoustic_coupling_el, ispec_is_acoustic_coupling_ac, potential_dot_dot_acoustic, veloc_elastic,&
                         DIR_RIGHT, DIR_LEFT, DIR_UP, DIR_DOWN, &
                         buffer_DG_rho_P, buffer_DG_rhovx_P, buffer_DG_rhovz_P, buffer_DG_E_P, NPROC, &
                         buffer_DG_Vxx_P, buffer_DG_Vzz_P, buffer_DG_Vxz_P, buffer_DG_Vzx_P, buffer_DG_Tz_P, buffer_DG_Tx_P, &
                         p_DG_init, gammaext_DG, muext, etaext, kappa_DG, ibool, cnu, &
                         buffer_DG_gamma_P, coord,  &
                         rho_init, rhovx_init, rhovz_init, E_init, &
                         potential_dot_dot_acoustic, &
                         veloc_vector_acoustic_DG_coupling, MPI_transfer_iface
                         
  implicit none
  
  ! TODO: Why not loading it from specfem_par? All calls seem to use the one loaded from specfem_par.
  
  integer, intent(in) :: i, j, ispec, iface1, iface
  
  integer, dimension(3), intent(in) :: neighbor
  
  real(kind=CUSTOM_REAL), intent(out) :: rho_DG_P, rhovx_DG_P, rhovz_DG_P, &
        E_DG_P, veloc_x_DG_P, veloc_z_DG_P, p_DG_P, T_P, &
        Tx_DG_P, Tz_DG_P, Vxx_DG_P, Vzz_DG_P, Vzx_DG_P, Vxz_DG_P
        
  real(kind=CUSTOM_REAL), intent(in) :: timelocal
        
  logical, intent(inout) :: exact_interface_flux
  
  real(kind=CUSTOM_REAL) :: tx, tz, normal_v, tangential_v, &
        nx, nz, veloc_x, veloc_z, weight, &
        veloc_x_DG_iM, veloc_z_DG_iM, p_DG_iM, gamma_P, e1_DG_P
        
  real(kind=CUSTOM_REAL), intent(in) :: rho_DG_iM, E_DG_iM, rhovx_DG_iM, rhovz_DG_iM, &
                rho_DG_iP, E_DG_iP, rhovx_DG_iP, rhovz_DG_iP
  real(kind=CUSTOM_REAL), dimension(2), intent(in) :: T_DG_iP, T_DG_iM
  real(kind=CUSTOM_REAL), dimension(2, 2), intent(in) :: V_DG_iP, V_DG_iM
        
  ! Local variables     
  integer :: iglobM, i_el, j_el, ispec_el, i_ac, j_ac, ispec_ac, iglob, iglobP, ipoin, num_interface
  real(kind=CUSTOM_REAL), dimension(2, 2) :: trans_boundary
  
  ! Parameters
  real(kind=CUSTOM_REAL), parameter :: ZERO = 0._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: ONE  = 1._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: TWO  = 2._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: HALF = 0.5_CUSTOM_REAL
  
  ! Characteristic based BC
  real(kind=CUSTOM_REAL) :: p_b, rho_b, s_b, rho_inf, v_b_x, v_b_z, c_b, un_b, &
  rlambda_max, rlambda_min, c_in, c_inf, un_in, un_inf, &
  deltaZ1, deltaZ2star, p_n, a_n, alpha0
  
  iglobM = ibool_DG(i, j, ispec)
  
  ! Deduce velocities from momenta and density (the two latter being constitutive variables, and thus passed to the routine as input parameters).
  veloc_x_DG_iM = rhovx_DG_iM/rho_DG_iM
  veloc_z_DG_iM = rhovz_DG_iM/rho_DG_iM
  ! Deduce pressure from energy (consitutive), density (consitutive), and velocities (deduced above).
  p_DG_iM       = (gammaext_DG(iglobM) - ONE)*( E_DG_iM &
        - (HALF)*rho_DG_iM*( veloc_x_DG_iM**2 + veloc_z_DG_iM**2 ) )
  
  Vxx_DG_P = -V_DG_iM(1, 1)
  Vzz_DG_P = -V_DG_iM(2, 2)
  Vxz_DG_P = -V_DG_iM(1, 2)
  Vzx_DG_P = -V_DG_iM(2, 1)
  Tx_DG_P  = -T_DG_iM(1)
  Tz_DG_P  = -T_DG_iM(2)
  
  gamma_P = gammaext_DG(iglobM)
  
  e1_DG_P = ZERO
  
  exact_interface_flux = .false. ! If set to .true., an exact flux formula will be used. If set to .false., the Lax-Friedrich approximation for the flux will be used.
  MPI_change_cpt = .false.
  
  if(neighbor(3) == -1) then
    ! --------------------------- !
    ! neighbor(3) == -1.          !
    ! --------------------------- !
    
    ipoin         = -1
    num_interface = -1
    if(NPROC > 1) then
      ipoin         = MPI_transfer_iface(iface1, iface, ispec, 1)
      num_interface = MPI_transfer_iface(iface1, iface, ispec, 2)
    endif                  
    
    if(ipoin > -1) then
      ! --------------------------- !
      ! neighbor(3) == -1,          !
      !   not diagonal corner       !
      !   element and not outside   !
      !   element.                  !
      ! --------------------------- !
      
      if(ispec_is_acoustic_coupling_ac(ibool_DG(i, j, ispec)) >= 0 .AND. .false.) then
        ! --------------------------- !
        ! COUPLING ACOUSTIC           !
        ! POTENTIAL - FLUID           !
        ! --------------------------- !
        ! TODO: Decide what to do with this case (we are in a 'if(.false.)').
        
        ! We already know the "real" flux at boundary.
        exact_interface_flux = .true. ! If set to .true., an exact flux formula will be used. If set to .false., the Lax-Friedrich approximation for the flux will be used.

        ! Coordinates of elastic element
        iglob = ibool(i, j, ispec)

        ! Only for density
        call boundary_condition_DG(i, j, ispec, timelocal, rho_DG_P, rhovx_DG_P, rhovz_DG_P, E_DG_P, &
        veloc_x_DG_P, veloc_z_DG_P, p_DG_P, e1_DG_P)

        ! derivatives of potential
        veloc_x = veloc_vector_acoustic_DG_coupling(iglob, 1)
        veloc_z = veloc_vector_acoustic_DG_coupling(iglob, 2)

        ! Tangential vector
        ! Since only bottom topography nz > 0
        tx = -nz
        tz = nx

        ! Normal velocity of the solid perturbation and tangential velocity of the fluid flow
        normal_v     = veloc_x*nx + veloc_x*nz
        tangential_v = veloc_x*tx + veloc_x*tz

        ! Set the matrix for the transformation from normal/tangential coordinates to mesh coordinates.
        trans_boundary(1, 1) =  tz
        trans_boundary(1, 2) = -nz
        trans_boundary(2, 1) = -tx
        trans_boundary(2, 2) =  nx
        trans_boundary = trans_boundary/(nx * tz - tx * nz)

        ! From free slip and normal velocity continuity
        veloc_x_DG_P = trans_boundary(1, 1)*normal_v + trans_boundary(1, 2)*tangential_v!veloc_elastic(1,iglob)
        veloc_z_DG_P = trans_boundary(2, 1)*normal_v + trans_boundary(2, 2)*tangential_v

        ! From traction continuity
        p_DG_P = p_DG_init(iglobM) - potential_dot_dot_acoustic(iglob)

        E_DG_P       = p_DG_P/(gammaext_DG(iglobM) - 1.) + HALF*rho_DG_P*( veloc_x_DG_P**2 + veloc_z_DG_P**2 )

        rhovx_DG_P   = rho_DG_P*veloc_x_DG_P
        rhovz_DG_P   = rho_DG_P*veloc_z_DG_P

        T_P = (E_DG_iM/rho_DG_iM - 0.5*(veloc_x_DG_iM**2 + veloc_z_DG_iM**2))/(cnu)
        
      else ! if(ispec_is_acoustic_coupling_ac(i, j, ispec))
        ! --------------------------- !
        ! CLASSICAL MPI NEIGHBOR (NOT !
        ! ACOUSTIC)                   !
        ! --------------------------- !
        
        ! Get all values from the buffer.
        rho_DG_P     = buffer_DG_rho_P(ipoin, num_interface)
        E_DG_P       = buffer_DG_E_P(ipoin, num_interface)
        rhovx_DG_P   = buffer_DG_rhovx_P(ipoin, num_interface)
        rhovz_DG_P   = buffer_DG_rhovz_P(ipoin, num_interface)
        gamma_P = buffer_DG_gamma_P(ipoin,num_interface)
        if(muext(i, j, ispec) > 0 .OR. &
           etaext(i, j, ispec) > 0 .OR. &
           kappa_DG(i, j, ispec) > 0) then
          ! If there is viscosity, get the values from the buffer.
          ! If not, values as initialised above are kept.
          Vxx_DG_P = buffer_DG_Vxx_P(ipoin, num_interface)
          Vzz_DG_P = buffer_DG_Vzz_P(ipoin, num_interface)
          Vxz_DG_P = buffer_DG_Vxz_P(ipoin, num_interface)
          Vzx_DG_P = buffer_DG_Vzx_P(ipoin, num_interface)
          Tx_DG_P = buffer_DG_Tx_P(ipoin, num_interface)
          Tz_DG_P = buffer_DG_Tz_P(ipoin, num_interface)
        endif
        
        ! Deduce velocities, pressure, and temperature.
        veloc_x_DG_P = rhovx_DG_P/rho_DG_P
        veloc_z_DG_P = rhovz_DG_P/rho_DG_P
        p_DG_P = (gamma_P - ONE)*( E_DG_P &
                 - (HALF)*rho_DG_P*( veloc_x_DG_P**2 + veloc_z_DG_P**2 ) )
        T_P = (E_DG_P/rho_DG_P - 0.5*((rhovx_DG_P/rho_DG_P)**2 + (rhovz_DG_P/rho_DG_P)**2))/(cnu)

      endif
                    
    elseif(ACOUSTIC_FORCING .AND. ispec_is_acoustic_forcing(i, j, ispec)) then
      ! --------------------------- !
      ! neighbor(3) == -1,          !
      !   acoustic forcing.         !
      ! --------------------------- !
      stop 'ACOUSTIC_FORCING obsolete for DG simulations'
         
    elseif(ispec_is_acoustic_coupling_el(i, j, ispec, 3) >= 0) then
      ! --------------------------- !
      ! neighbor(3) == -1,          !
      !   elastic coupling.         !
      ! --------------------------- !
      ! TODO: Clean up all this part.
      
      exact_interface_flux = .false. ! If set to .true., an exact flux formula will be used. If set to .false., the Lax-Friedrich approximation for the flux will be used.

      ! Coordinates of elastic element
      i_el     = ispec_is_acoustic_coupling_el(i, j, ispec, 1)
      j_el     = ispec_is_acoustic_coupling_el(i, j, ispec, 2)
      ispec_el = ispec_is_acoustic_coupling_el(i, j, ispec, 3)

      iglob = ibool(i_el, j_el, ispec_el)
      
      call boundary_condition_DG(i, j, ispec, timelocal, rho_DG_P, rhovx_DG_P, rhovz_DG_P, E_DG_P, &
                                 veloc_x_DG_P, veloc_z_DG_P, p_DG_P, e1_DG_P)

      rho_DG_P = rho_DG_iM 

      ! Elastic velocities
      veloc_x = veloc_elastic(1, iglob)
      veloc_z = veloc_elastic(2, iglob)

      ! Tangential vector
      ! Since only bottom topography nz > 0
      tx = -nz
      tz =  nx

      ! Normal velocity of the solid perturbation and tangential velocity of the fluid flow
      normal_v     = (veloc_x*nx + veloc_z*nz)
      tangential_v = veloc_x_DG_P*tx + veloc_z_DG_P*tz
      !normal_v     = -(veloc_x_DG_iM*nx + veloc_z_DG_iM*nz) + 2*(veloc_x*nx + veloc_z*nz) 
      !tangential_v = veloc_x_DG_iM*tx + veloc_z_DG_iM*tz
      !tangential_v    = (veloc_x*tx + veloc_z*tz) 

      ! Set the matrix for the transformation from normal/tangential coordinates to mesh coordinates.
      trans_boundary(1, 1) =  tz
      trans_boundary(1, 2) = -nz
      trans_boundary(2, 1) = -tx
      trans_boundary(2, 2) =  nx
      trans_boundary = trans_boundary/(nx * tz - tx * nz)

      ! Convert (back) the velocity components from normal/tangential coordinates to mesh coordinates.
      veloc_x_DG_P = trans_boundary(1, 1)*normal_v + trans_boundary(1, 2)*tangential_v!veloc_elastic(1,iglob)
      veloc_z_DG_P = trans_boundary(2, 1)*normal_v + trans_boundary(2, 2)*tangential_v

      ! No stress continuity.
      p_DG_P = p_DG_iM

      ! Deduce energy.
      E_DG_P = p_DG_P/(gammaext_DG(iglobM) - 1.) + HALF*rho_DG_P*( veloc_x_DG_P**2 + veloc_z_DG_P**2 )

      ! Recompute momenta in order not to keep those coming from the call to boundary_condition_DG above.
      rhovx_DG_P = rho_DG_P*veloc_x_DG_P
      rhovz_DG_P = rho_DG_P*veloc_z_DG_P
      
      ! Deduce temperature.
      T_P = (E_DG_iM/rho_DG_iM - 0.5*(veloc_x_DG_iM**2 + veloc_z_DG_iM**2))/(cnu)

    else
      ! --------------------------- !
      ! neighbor(3) == -1,          !
      !   classical boundary        !
      !   conditions.               !
      ! --------------------------- !
    
      exact_interface_flux = .false. ! If set to .true., an exact flux formula will be used. If set to .false., the Lax-Friedrich approximation for the flux will be used.

      ! Get "background" (or "far-field", or "unperturbed") values.
      call boundary_condition_DG(i, j, ispec, timelocal, rho_DG_P, rhovx_DG_P, rhovz_DG_P, E_DG_P, &
                                 veloc_x_DG_P, veloc_z_DG_P, p_DG_P, e1_DG_P)
      
      ! Set the velocity. --------- !
      ! Convert the velocity components from mesh coordinates to normal/tangential coordinates.
      ! This is more practical to set the boundary conditions.
      ! The setting of the boundary conditions is thus made here.
      tx = -nz
      tz =  nx
      normal_v = (veloc_x_DG_P*nx + veloc_z_DG_P*nz)
      tangential_v = (veloc_x_DG_P*tx + veloc_z_DG_P*tz)
      
      ! Notes:
      !   At that point, normal_v and tangential_v are set to their "background" (or "far-field", or "unperturbed") values.
      !   Treatment of boundary conditions based on normal/tangential velocities should be done here. The "free slip" condition and the "normal velocity continuity" conditions can be set here.
      
      ! Set the matrix for the transformation from normal/tangential coordinates to mesh coordinates.
      trans_boundary(1, 1) =  tz
      trans_boundary(1, 2) = -nz
      trans_boundary(2, 1) = -tx
      trans_boundary(2, 2) =  nx
      trans_boundary = trans_boundary/(nx * tz - tx * nz)
     
      ! Convert (back) the velocity components from normal/tangential coordinates to mesh coordinates.
      veloc_x_DG_P = trans_boundary(1, 1)*normal_v + trans_boundary(1, 2)*tangential_v
      veloc_z_DG_P = trans_boundary(2, 1)*normal_v + trans_boundary(2, 2)*tangential_v
      
      ! Set the energy. ----------- !
      E_DG_P = p_DG_P/(gammaext_DG(iglobM) - ONE) &
               + rho_DG_P*HALF*( veloc_x_DG_P**2 + veloc_z_DG_P**2 )
      
      ! Set the momenta. ---------- !
      rhovx_DG_P   = rho_DG_P*veloc_x_DG_P
      rhovz_DG_P   = rho_DG_P*veloc_z_DG_P
      
      
      ! TODO: Decide what to do with those cases below.
      ! Those are test cases, intended to test characteristic-based boundary conditions.
      if(coord(2, ibool(i, j, ispec)) == 20000. .AND. .false.) then
        ! --------------------------- !
        ! neighbor(3) == -1,          !
        !   classical boundary        !
        !   conditions,               !
        !     characteristic based BC,!
        !     top boundary.           !
        ! --------------------------- !
        
        un_in       = veloc_x_DG_iM*nx + veloc_z_DG_iM*nz
        un_inf      = veloc_x_DG_P*nx  + veloc_z_DG_P*nz
        c_in        = sqrt(gammaext_DG(iglobM)*p_DG_iM/rho_DG_iM)
        c_inf       = sqrt(gammaext_DG(iglobM)*p_DG_P/rho_DG_P)
        
        rlambda_max = un_in  - 2.*c_in/(gammaext_DG(iglobM) - 1.)
        rlambda_min = un_inf + 2.*c_inf/(gammaext_DG(iglobM) - 1.)
        
        un_b = (rlambda_min + rlambda_max)/2.
        c_b  = (gammaext_DG(iglobM) - 1.)*(rlambda_min - rlambda_max)/4.
        
        ! Boudnary values for velocity
        v_b_x  = veloc_x_DG_P + ( un_b - un_inf )*nx
        v_b_z  = veloc_z_DG_P + ( un_b - un_inf )*nz
        
        rho_inf = rho_DG_P
        s_b = (c_in**2)/( gammaext_DG(iglobM)*(rho_DG_iM**(gammaext_DG(iglobM) - 1.)) )
        
        ! Boundary values for density and pressure
        rho_b = ((c_b**2)/(gammaext_DG(iglobM)*s_b))**(1/(gammaext_DG(iglobM)-1.))
        p_b    = rho_b*(c_b**2)/gammaext_DG(iglobM)
        
        rho_DG_P      = -rho_DG_iM + rho_b
        veloc_x_DG_P  = v_b_x!veloc_x_DG_P -2.*( normal_v + veloc_x_DG_P*nx + veloc_z_DG_P*nz )*nx
        veloc_z_DG_P  = v_b_z!veloc_z_DG_P -2.*( normal_v + veloc_x_DG_P*nx + veloc_z_DG_P*nz )*nz
        p_DG_P = p_b
        E_DG_P = p_DG_P/(gammaext_DG(iglobM) - ONE) &
                + rho_DG_P*HALF*( veloc_x_DG_P**2 + veloc_z_DG_P**2 )
        
        rhovx_DG_P   = -rhovx_DG_iM + 2.*rho_DG_P*veloc_x_DG_P
        rhovz_DG_P   = -rhovz_DG_iM + 2.*rho_DG_P*veloc_z_DG_P
        rho_DG_P     = -rho_DG_iM + 2.*rho_DG_P
        !veloc_x_DG_P = -veloc_x_DG_iM + 2.*veloc_x_DG_P
        !veloc_z_DG_P = -veloc_z_DG_iM + 2.*veloc_z_DG_P     
        E_DG_P       = -E_DG_iM + 2.*E_DG_P
      endif
      if(coord(2, ibool(i, j, ispec)) == 0. .AND. .false.) then
        ! --------------------------- !
        ! neighbor(3) == -1,          !
        !   classical boundary        !
        !   conditions,               !
        !     characteristic based BC,!
        !     bottom boundary.        !
        ! --------------------------- !
        
        p_n = (gammaext_DG(iglobM) - ONE)*( E_init(iglobM) &
                - (HALF/rho_init(iglobM))*( rhovx_init(iglobM)**2 + rhovz_init(iglobM)**2 ) )
        a_n = sqrt(gammaext_DG(iglobM)*p_n/rho_init(iglobM))
        deltaZ1 = 2*p_DG_P - p_DG_iM - p_n &
                - rho_init(iglobM)*a_n&
                        *((veloc_x_DG_iM - rhovx_init(iglobM)/rho_init(iglobM))*nx &
                                +(veloc_z_DG_iM-rhovz_init(iglobM)/rho_init(iglobM))*nz)

        deltaZ2star = (p_DG_iM - p_n) &
           + rho_init(iglobM)*a_n*((veloc_x_DG_iM - rhovx_init(iglobM)/rho_init(iglobM))*nx &
                                 + (veloc_z_DG_iM-rhovz_init(iglobM)/rho_init(iglobM))*nz)

       ! rho_DG_P = rho_DG_iM
        alpha0 = 0.0001          
        
        p_DG_P = p_n - (alpha0*deltaZ1 + deltaZ2star)/2
        
        normal_v     = ((rhovx_init(iglobM)/rho_init(iglobM))*nx+(rhovz_init(iglobM)/rho_init(iglobM))*nz) &
                - (alpha0*deltaZ1 - deltaZ2star)/(2.*a_n*rho_init(iglobM))
        tangential_v = 0.

        ! Set the matrix for the transformation from normal/tangential coordinates to mesh coordinates.
        trans_boundary(1, 1) =  tz
        trans_boundary(1, 2) = -nz
        trans_boundary(2, 1) = -tx
        trans_boundary(2, 2) =  nx
        trans_boundary = trans_boundary/(nx * tz - tx * nz)

        ! From free slip and normal velocity continuity
        veloc_x_DG_P = trans_boundary(1, 1)*normal_v + trans_boundary(1, 2)*tangential_v!veloc_elastic(1,iglob)
        veloc_z_DG_P = trans_boundary(2, 1)*normal_v + trans_boundary(2, 2)*tangential_v  
      endif
      
      ! Set the temperature. ------ !
      T_P = (E_DG_P/rho_DG_P - 0.5*(veloc_x_DG_P**2 + veloc_z_DG_P**2))/(cnu)
      !T_P = (E_DG_iM/rho_DG_iM - 0.5*(veloc_x_DG_iM**2 + veloc_z_DG_iM**2))/(cnu) ! DEBUG
      !Tx_DG_P = -T_DG_iM(1) ! DEBUG
      !Tz_DG_P = -T_DG_iM(2) ! DEBUG
    endif
  elseif(ispec_is_acoustic_coupling_ac(ibool_DG(i, j, ispec)) >= 0 .AND. .false.) then
    ! --------------------------- !
    ! neighbor(3) != -1,          !
    !   acoustic coupling.        !
    ! --------------------------- !
    ! TODO: Decide what to do with this case.

    !!elseif(.false.) then

    !!WRITE(*,*) ">>>>>>>>>>>> TOTO" ! DEBUG

    exact_interface_flux = .true. ! If set to .true., an exact flux formula will be used. If set to .false., the Lax-Friedrich approximation for the flux will be used.

    ! Coordinates of elastic element
    i_ac     = neighbor(1)
    j_ac     = neighbor(2)
    ispec_ac = neighbor(3)

    iglob = ibool(i_ac,j_ac,ispec_ac)!ibool(i, j, ispec)
    if(.false.) then
      WRITE(*,*) ibool(i_ac,j_ac,ispec_ac), xixl, xizl, gammaxl, gammazl, duz_dxi, duz_dgamma, k ! DEBUG
    endif
    ! Only for density
    call boundary_condition_DG(i_ac, j_ac, ispec_ac, timelocal, rho_DG_P, rhovx_DG_P, rhovz_DG_P, E_DG_P, &
    veloc_x_DG_P, veloc_z_DG_P, p_DG_P, e1_DG_P)

    !!rho_DG_P = rho_DG_iM

    !!duz_dxi    = ZERO
    !!duz_dgamma = ZERO

    !!xixl = xix(i_ac,j_ac,ispec_ac)
    !!xizl = xiz(i_ac,j_ac,ispec_ac)
    !!gammaxl = gammax(i_ac,j_ac,ispec_ac)
    !!gammazl = gammaz(i_ac,j_ac,ispec_ac)

    ! first double loop over GLL points to compute and store gradients
    ! we can merge the two loops because NGLLX == NGLLZ
    !do k = 1, NGLLX
    !         duz_dxi    = duz_dxi    &
    !                 + potential_dot_acoustic(ibool(k,j_ac,ispec_ac)) * real(hprime_xx(i_ac,k), kind=CUSTOM_REAL)
    !         duz_dgamma = duz_dgamma &
    !                 + potential_dot_acoustic(ibool(i_ac,k,ispec_ac)) * real(hprime_zz(j_ac,k), kind=CUSTOM_REAL)
    !enddo

    !!! Derivatives of potential.
    !!veloc_x = ( duz_dxi * xixl + duz_dgamma * gammaxl )
    !!veloc_z = ( duz_dxi * xizl + duz_dgamma * gammazl )

    veloc_x = veloc_vector_acoustic_DG_coupling(iglob, 1)
    veloc_z = veloc_vector_acoustic_DG_coupling(iglob, 2)

    ! Tangential vector, since only bottom topography nz > 0.
    tx = -nz
    tz =  nx

    ! Normal velocity of the solid perturbation and tangential velocity of the fluid flow
    normal_v     = veloc_x*nx + veloc_x*nz
    tangential_v = veloc_x*tx + veloc_x*tz

    ! Set the matrix for the transformation from normal/tangential coordinates to mesh coordinates.
    trans_boundary(1, 1) =  tz
    trans_boundary(1, 2) = -nz
    trans_boundary(2, 1) = -tx
    trans_boundary(2, 2) =  nx
    trans_boundary = trans_boundary/(nx * tz - tx * nz)

    ! From free slip and normal velocity continuity
    veloc_x_DG_P = trans_boundary(1, 1)*normal_v + trans_boundary(1, 2)*tangential_v!veloc_elastic(1,iglob)
    veloc_z_DG_P = trans_boundary(2, 1)*normal_v + trans_boundary(2, 2)*tangential_v

    ! From traction continuity
    p_DG_P = p_DG_init(iglobM) - potential_dot_dot_acoustic(iglob)

    E_DG_P       = p_DG_P/(gammaext_DG(iglobM) - 1.) + HALF*rho_DG_P*( veloc_x_DG_P**2 + veloc_z_DG_P**2 )

    rhovx_DG_P   = rho_DG_P*veloc_x_DG_P
    rhovz_DG_P   = rho_DG_P*veloc_z_DG_P

    T_P = (E_DG_iM/rho_DG_iM - 0.5*(veloc_x_DG_iM**2 + veloc_z_DG_iM**2))/(cnu)
  
  else
    ! --------------------------- !
    ! neighbor(3) != -1,          !
    !   not an outside edge,      !
    !   thus values are known.    !
    ! --------------------------- !
    iglobP = ibool_DG(neighbor(1), neighbor(2), neighbor(3))
    
    gamma_P = gammaext_DG(iglobP)
    
    rho_DG_P     = rho_DG_iP
    E_DG_P       = E_DG_iP
    rhovx_DG_P   = rhovx_DG_iP
    rhovz_DG_P   = rhovz_DG_iP
    
    veloc_x_DG_P = rhovx_DG_P/rho_DG_P
    veloc_z_DG_P = rhovz_DG_P/rho_DG_P
    p_DG_P       = (gamma_P - ONE)*( E_DG_P &
                   - (HALF)*rho_DG_P*( veloc_x_DG_P**2 + veloc_z_DG_P**2 ) )
    
    if(muext(i, j, ispec) > 0 .OR. &
       etaext(i, j, ispec) > 0 .OR. &
       kappa_DG(i, j, ispec) > 0) then
      ! Viscosity  
      Vxx_DG_P = V_DG_iP(1, 1)!,iglobP)
      Vzz_DG_P = V_DG_iP(2, 2)!,iglobP)
      Vxz_DG_P = V_DG_iP(1, 2)!,iglobP)
      Vzx_DG_P = V_DG_iP(2, 1)!,iglobP)
      Tx_DG_P = T_DG_iP(1)!,iglobP)
      Tz_DG_P = T_DG_iP(2)!,iglobP)
    endif
    T_P = (E_DG_P/rho_DG_P - 0.5*((rhovx_DG_P/rho_DG_P)**2 + (rhovz_DG_P/rho_DG_P)**2))/(cnu)
  endif
  
  !close(10)
  
  ! Temperature computation
  !T_P = (E_DG_P/rho_DG_P - 0.5*((rhovx_DG_P/rho_DG_P)**2 + (rhovz_DG_P/rho_DG_P)**2))/(cnu)
  ! (E_DG_iM/rho_DG_iM - 0.5*(veloc_x_DG_iM**2 + veloc_z_DG_iM**2))/(cnu)
  
  end subroutine compute_interface_unknowns

! ------------------------------------------------------------ !
! damp_solution_DG                                             !
! ------------------------------------------------------------ !
! TODO: This is a test routine.

  subroutine damp_solution_DG(rho_DG, rhovx_DG, rhovz_DG, E_DG, timelocal)

  use specfem_par, only: ispec_is_acoustic, nspec, coord, ibool_DG, nglob_DG, &
        ibool_before_perio
  use constants,only: CUSTOM_REAL,NGLLX,NGLLZ

  implicit none
  
  real(kind=CUSTOM_REAL), dimension(nglob_DG) :: rho_DG, rhovx_DG, rhovz_DG, E_DG
  
  ! Local
  real(kind=CUSTOM_REAL) :: L_buffer_absorb, z, sigma, zmax, z_l
  ! Arina's damping coefficients.
  real(kind=CUSTOM_REAL) :: C_1, C_2
  ! Richards' damping coefficients.
  real(kind=CUSTOM_REAL) :: beta, sigma_max
  
  real(kind=CUSTOM_REAL) :: timelocal, rho_DG_P, rhovx_DG_P, rhovz_DG_P, E_DG_P, &
                            veloc_x_DG_P, veloc_z_DG_P, p_DG_P, e1_DG_P
  integer :: i, j, ispec, ifirstelem, ilastelem, ibool
  
  ! Domain-related quantities.
  L_buffer_absorb = 20.0d0 ! Length of the buffer in which the solution is damped.
  zmax            = 50.0d0
  
  ! Arina's damping coefficients.
  C_1 = 0.0d0 ! 0 in Arina's paper, 0<=C_1<=0.1 in Wasistho's.
  C_2 = 13.0d0 ! 13 in Arina's paper, 10<=C_2<=20 in Wasistho's.
  ! Richards' damping coefficients.
  beta = 4.0d0 ! 4 in Richards' paper.
  sigma_max = -1.0d0
  
  ifirstelem = 1
  ilastelem = nspec
  
  do ispec = ifirstelem, ilastelem
    if (ispec_is_acoustic(ispec)) then
      do j = 1, NGLLZ
        do i = 1, NGLLX
          ibool = ibool_DG(i, j, ispec)
          
          ! Absolute vertical position of the GLL point.
          z = coord(2, ibool_before_perio(i, j, ispec))
          ! Relative vertical position of the GLL point with respect to:
          ! 1) the upper limit of the simulated domain, and
          ! 2) the wanted buffer length.
          ! z_l = 1 on the outer limit of the buffer, z_l = 0 at the inner limit of the buffer, z_l > 1 outside domain, and z_l < 0 inside the domain (before the buffer).
          z_l = (1.0d0 - ((zmax - z)/L_buffer_absorb))
          
          !sigma = sigma_max*(abs( 1. - (z - L_buffer_absorb)/L_buffer_absorb ))**beta
          
          if(z_l > 0.0d0 .AND. z_l <= 1.0d0) then
            ! Arina's damping.
            sigma = (1.0d0 - C_1*z_l**2.0d0)*(1.0d0 - ( 1.0d0 - exp(C_2*(z_l)**2.0d0) )/( 1.0d0 - exp(C_2) ))
            ! Richards' damping.
            !sigma = 1.0d0 + sigma_max * z_l ** beta
            
            ! Note: sigma is a function that goes gradually from 1 to 0 in the buffer.
            
            call boundary_condition_DG(i, j, ispec, timelocal, rho_DG_P, rhovx_DG_P, rhovz_DG_P, E_DG_P, &
                    veloc_x_DG_P, veloc_z_DG_P, p_DG_P, e1_DG_P)
            
            !WRITE(*, *) timelocal, z,z_l,abs(z - zmax),abs(z - zmax)/L_buffer_absorb,&
            ! sigma!, rho_DG(ibool) , rho_DG_P, sigma*( rho_DG(ibool) -  rho_DG_P), &
                   ! rho_DG(ibool) - sigma*( rho_DG(ibool) -  rho_DG_P) ! DEBUG
             
            !if(sigma <= 1.0d-4) WRITE(*, *) sigma, z_l, (1-exp(C_2*(z_l)**2.0d0))/(1-exp(C_2)) ! DEBUG
            
            ! 1st approach (set only the perturbation to zero).
            ! For a variable v: v_{tot} = v_{init} + v'. Then do:
            ! v_{tot}^{n+1} = v_{init} + sigma * (v_{tot, computed} - v_{init})
            rho_DG(ibool)   = rho_DG_P + sigma*( rho_DG(ibool) - rho_DG_P)
            rhovx_DG(ibool) = rhovx_DG_P + sigma*( rhovx_DG(ibool) - rhovx_DG_P)
            rhovz_DG(ibool) = rhovz_DG_P + sigma*( rhovz_DG(ibool) - rhovz_DG_P)
            E_DG(ibool)     = E_DG_P + sigma*( E_DG(ibool) - E_DG_P)
            
            ! 2nd approach (do not use since in particular it sets rho to zero, which is not physically acceptable).
            !rho_DG(ibool)   = rho_DG(ibool) * sigma
            !rhovx_DG(ibool) = rhovx_DG(ibool) * sigma
            !rhovz_DG(ibool) = rhovz_DG(ibool) * sigma
            !E_DG(ibool)     = E_DG(ibool) * sigma
          endif
        enddo
      enddo
     endif
   enddo
  end subroutine damp_solution_DG

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
        rho_DG, rhovx_DG, rhovz_DG, E_DG, e1_DG, nspec!, ibool_before_perio, coord

  implicit none
  
  ! Parameters.
  real(kind=CUSTOM_REAL), parameter :: ZEROl = 0._CUSTOM_REAL
  !real(kind=CUSTOM_REAL), parameter :: ONEl  = 1._CUSTOM_REAL
  !real(kind=CUSTOM_REAL), parameter :: TWOl  = 2._CUSTOM_REAL
  !real(kind=CUSTOM_REAL), parameter :: SIXl  = 6._CUSTOM_REAL
  !real(kind=CUSTOM_REAL), parameter :: HALFl = 0.5_CUSTOM_REAL
  
  integer :: ispec, iglob, i, j
  
  real(kind=CUSTOM_REAL) :: rho_DG_P, rhovx_DG_P, rhovz_DG_P, E_DG_P, &
      veloc_x_DG_P, veloc_z_DG_P, p_DG_P, e1_DG_P
  
    
  do ispec = 1, nspec ! Loop over elements.
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
       !if(abs(coord(2,ibool_before_perio(i,j,ispec)))<1.) then ! DEBUG
       !  write(*,*) coord(:,ibool_before_perio(i,j,ispec)), & ! DEBUG
       !             rho_DG_P, rhovx_DG_P, rhovz_DG_P, E_DG_P ! DEBUG
       !endif ! DEBUG
      enddo
    enddo
  enddo  
  
  !if(.false.) then
  !  ! TODO: Why 'if(.false.)'?
  !  call recompute_density()
  !endif
  
  end subroutine initial_condition_DG
  
! ------------------------------------------------------------ !
! recompute_density                                            !
! ------------------------------------------------------------ !
! TODO: Description.
  
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

  use constants,only: CUSTOM_REAL, gamma_euler, PI, HUGEVAL, TINYVAL, FOUR_THIRDS

  use specfem_par, only: ibool_before_perio, ibool_DG, coord, &
        rhoext, windxext, pext_DG, gravityext, gammaext_DG, &
        etaext, muext, coord_interface, kappa_DG, cp, c_V, &
        tau_epsilon, tau_sigma, &
        gravity_cte_DG, dynamic_viscosity_cte_DG, thermal_conductivity_cte_DG, tau_eps_cte_DG, tau_sig_cte_DG, SCALE_HEIGHT, &
        USE_ISOTHERMAL_MODEL, potential_dphi_dx_DG, potential_dphi_dz_DG, ibool, &
        surface_density, sound_velocity, wind, TYPE_FORCING, &
        !forcing_initial_time, main_time_period, forcing_initial_loc, main_spatial_period,&
        assign_external_model!,myrank!, &
        !DT,&!, XPHASE_RANDOMWALK, TPHASE_RANDOMWALK, PHASE_RANDOMWALK_LASTTIME,& ! Microbarom forcing.
        !EXTERNAL_FORCING_MAXTIME,EXTERNAL_FORCING, EXTFORC_MAP_ibbp_TO_LOCAL,& ! External forcing.
        !EXTFORC_MINX, EXTFORC_MAXX,EXTFORC_FILEDT! External forcing.

  implicit none
  
  integer :: i, j, ispec
  
  real(kind=CUSTOM_REAL) :: timelocal
  real(kind=CUSTOM_REAL), intent(out) :: rho_DG_P, rhovx_DG_P, rhovz_DG_P, E_DG_P, &
      veloc_x_DG_P, veloc_z_DG_P, p_DG_P, e1_DG_P
  
  ! Parameters
  real(kind=CUSTOM_REAL), parameter :: ZEROl = 0._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: ONEl  = 1._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: TWOl  = 2._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: SIXl  = 6._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: HALFl = 0.5_CUSTOM_REAL
  
  real(kind=CUSTOM_REAL) :: x, z, H, G!, A
  !real(kind=CUSTOM_REAL) :: x0, z0, r, beta
  ! Forcing
  !real(kind=CUSTOM_REAL) :: to, perio, lambdo, xo
  !real(kind=CUSTOM_REAL) :: mu_visco, eta_visco
  ! Hydrostatic solution
  !real(kind=CUSTOM_REAL) :: RR, p0, rho0
  ! Density current
  !real(kind=CUSTOM_REAL) :: cp, c_V, exner, RR, p0, rho0, rs, theta, theta0
  ! Linear mountains
  !real(kind=CUSTOM_REAL) :: cp, c_V, exner, RR, p0, rho0, rs, theta, theta0, Nsq
  ! Tsunami
  !real(kind=CUSTOM_REAL) :: VELOC_TSUNAMI
  ! Microbaroms.
  !real(kind=CUSTOM_REAL) :: MICROBAROM_AMPLITUDE, MICROBAROM_MAXTIME, MICROBAROM_RANGE
  !real(kind=CUSTOM_REAL) :: UNIFORM1, UNIFORM2, NORMAL1, NORMAL2
  ! External focing.
  !integer :: externalforcingid
  
  ! Thermal bubble
  real(kind=CUSTOM_REAL) :: Tl, Tu, rho0, p0, RR
  !real(kind=CUSTOM_REAL) :: Tl, Tu, p0, rho0, RR, rs!, theta, thetaprime, pibar
  !real(kind=CUSTOM_REAL) :: Tl, Tu, p0, RR, rs, theta0, Nsq, rho0
  !real(kind=CUSTOM_REAL) :: cp, c_V, R, P0, N, theta0, exner, theta
  !real(kind=CUSTOM_REAL) :: x0, z0, r, rs, p0, RR, c_V, cp, theta0, exner
  !real(kind=CUSTOM_REAL) :: x0, z0, r, theta, RR, p0!, exner, theta0
  !real(kind=CUSTOM_REAL) :: &
  !      lz, vs, vs_x, vs_z, Ms, p_1, rho_1, veloc_x_1, veloc_z_1, p_2, rho_2, veloc_x_2, veloc_z_2 
  
  ! Coordinate of the solid-fluid interface
  !coord_interface = 304.!250.
  
  x = real(coord(1, ibool_before_perio(i, j, ispec)), kind=CUSTOM_REAL)
  z = real(coord(2, ibool_before_perio(i, j, ispec)), kind=CUSTOM_REAL)
  z = z - coord_interface
  
  if(assign_external_model) then
    ! If an external model data file is given for initial conditions, use it.
    rho_DG_P     = rhoext(i, j, ispec)
    p_DG_P       = pext_DG(i, j, ispec)
    !p_DG_P       = rho_DG_P * gravityext(i, j, ispec) * Htabext_DG(ibool_DG(i, j, ispec)) ! Hydrostaticity: impose p = \rho * g * H, where the RHS is gotten from the external model data file.
    veloc_x_DG_P = windxext(i, j, ispec) ! Read horizontal wind.
    veloc_z_DG_P = ZEROl ! Impose vertical wind to zero.
    if(timelocal == 0) then
      !gravityext(i, j, ispec) = G!9.831
      !gammaext_DG(ibool_DG(i, j, ispec)) = 1.4!(rho_DG_P/p_DG_P)*(vpext(i, j, ispec)**2)!gamma_euler!
      !muext(i, j, ispec) = 1.0d-05
      !etaext(i, j, ispec) = (4/3)*muext(i, j, ispec)
      !kappa_DG(i, j, ispec) = 0.025!kappa_DG(i, j, ispec) / 10.!0.2
      !cp = 7/2
      !c_V = 5/2
      !tau_epsilon(i, j, ispec) = 1.
      !tau_sigma(i, j, ispec) = 1.!0.4013
      potential_dphi_dx_DG(ibool(i, j, ispec)) = ZEROl
      potential_dphi_dz_DG(ibool(i, j, ispec)) = gravityext(i, j, ispec)
    endif
    
    ! Set auxiliary variable.
    e1_DG_P = ZEROl
  else
    ! If no external model data file is given (no initial conditions were specified), build model.
    ! Set gravity, viscosity coefficients, relaxation times, and gravity potential (Phi) derivatives.
    if(abs(timelocal)<TINYVAL) then
      ! > Isobaric case. Since we need to stay hydrostatic, the gravity field needs to stay 0.
      gravityext(i, j, ispec) = ZEROl !real(coord(1,ibool_before_perio(i, j, ispec)), kind=CUSTOM_REAL)
      ! > Isothermal case.
      if(USE_ISOTHERMAL_MODEL) then
        gravityext(i, j, ispec) = real(gravity_cte_DG, kind=CUSTOM_REAL)
      endif
      gammaext_DG(ibool_DG(i, j, ispec)) = cp/c_V
      muext(i, j, ispec) = dynamic_viscosity_cte_DG
      etaext(i, j, ispec) = FOUR_THIRDS*dynamic_viscosity_cte_DG
      kappa_DG(i, j, ispec) = thermal_conductivity_cte_DG
      tau_epsilon(i, j, ispec) = tau_eps_cte_DG
      tau_sigma(i, j, ispec)   = tau_sig_cte_DG
      
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
    veloc_x_DG_P = wind ! Read horizontal wind (from specfem2_par), which is the scalar value read from parfile.
    veloc_z_DG_P = ZEROl ! Impose vertical wind to zero.
  endif ! Endif on assign_external_model.
  
  ! TEST VISCOSITY in top buffer
  if(.false.) then
    call change_visco(i, j, ispec,&
                      real(coord(1, ibool_before_perio(i, j, ispec)), kind=CUSTOM_REAL),&
                      real(coord(2, ibool_before_perio(i, j, ispec)), kind=CUSTOM_REAL)) ! See 'prepare_stretching.f90'.
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
  
  ! If bottom forcing is activated, set (overwrite) velocities.
  if(TYPE_FORCING>0 .and. abs(z)<TINYVAL) then
    call forcing_DG(i, j, ispec, timelocal, veloc_z_DG_P) ! Used the subroutine 'forcing_DG' defined below, instead of crowding this one.
    ! Read forcing parameters from parameter file.
    !lambdo = main_spatial_period
    !xo     = forcing_initial_loc
    !perio  = main_time_period
    !to     = forcing_initial_time
    !! --------------------------- !
    !! Time Gaussian derivative    !
    !! (acoustic plane wave        !
    !! forcing).                   !
    !! --------------------------- !
    !!if(TYPE_FORCING == 1) then
    !!  call forcing_DG(i, j, ispec, timelocal, veloc_z_DG_P)
    !!endif
    !if(TYPE_FORCING == 1) then
    !  veloc_z_DG_P = 0.01*(&
    !                 - (2d0/(perio/4d0))*((timelocal-(to-perio/4d0))/(perio/4d0))* &
    !                       (exp(-((timelocal-(to-perio/4d0))/(perio/4d0))**2)) &
    !                 + (2d0/(perio/4d0))*((timelocal-(to+perio/4d0))/(perio/4d0))* &
    !                       (exp(-((timelocal-(to+perio/4d0))/(perio/4d0))**2)) )
    !endif
    !! --------------------------- !
    !! Time and space Gaussian     !
    !! derivative (gravity wave    !
    !! forcing.                    !
    !! --------------------------- !
    !if(TYPE_FORCING == 2) then
    !  veloc_z_DG_P = 0.001*(&
    !                - (2d0/(perio/4d0))*((timelocal-(to-perio/4d0))/(perio/4d0))* &
    !                         (exp(-((timelocal-(to-perio/4d0))/(perio/4d0))**2)) &
    !                + (2d0/(perio/4d0))*((timelocal-(to+perio/4d0))/(perio/4d0))* &
    !                         (exp(-((timelocal-(to+perio/4d0))/(perio/4d0))**2)) ) &
    !                 * ( exp(-((x-(xo-lambdo/4))/(lambdo/4))**2) - &
    !                        exp(-((x-(xo+lambdo/4))/(lambdo/4))**2) )
    !endif
    !! --------------------------- !
    !! Hardcoded custom forcing    !
    !! (custom forcing case for    !
    !! user).                      !
    !! --------------------------- !
    !! The Matlab script 'utils_new/forcings_test.m' contains some tests for some of those forcings, in order to make hardcoding a little easier.
    !if(TYPE_FORCING == 9) then
    !  ! Tsunami forcing.
    !  if(.false.) then        
    !    VELOC_TSUNAMI = 200.
    !    perio = 1000.
    !    lambdo = 50000.
    !    if(timelocal < perio) then
    !        veloc_z_DG_P = (1d0/perio)*exp( -((x-lambdo)/(sqrt(2d0)*perio))**2 )
    !    else
    !        veloc_z_DG_P = &
    !                2d0*((VELOC_TSUNAMI)/(sqrt(2d0)*perio))*(((x-lambdo)-VELOC_TSUNAMI*(timelocal - perio))/(sqrt(2d0)*perio))&
    !                *exp( -(((x-lambdo)-VELOC_TSUNAMI*(timelocal - perio))/(sqrt(2d0)*perio))**2 )
    !    endif
    !  endif
    !  ! Analytic microbarom forcing.
    !  if(.true.) then
    !    perio  = main_time_period
    !    lambdo = main_spatial_period
    !    MICROBAROM_AMPLITUDE = 1. ! Microbarom amplitude. Unit: m.
    !    MICROBAROM_RANGE = 80.*lambdo ! Range around x=0 to which impose microbaroms. Unit: m. Be careful with apodisation below.
    !    MICROBAROM_MAXTIME = 10.5 * perio ! Microbarom active from t=0 to t=MICROBAROM_MAXTIME. Unit: s. Be careful with apodisation below.
    !    if(timelocal==0) then
    !      ! Start random phase walk.
    !      XPHASE_RANDOMWALK = 0.
    !      TPHASE_RANDOMWALK = 0.
    !      PHASE_RANDOMWALK_LASTTIME=0.
    !    endif
    !    if(timelocal<MICROBAROM_MAXTIME) then
    !      if(timelocal>=PHASE_RANDOMWALK_LASTTIME+DT) then
    !        ! Update the random walk only once par time step.
    !        call random_number(UNIFORM1)
    !        call random_number(UNIFORM2)
    !        NORMAL1 = (PI*DT/(0.2*2.*perio)) * sqrt(-2.*log(UNIFORM1))*cos(2.*PI*UNIFORM2) ! Box-Muller method to generate a 1nd N(0, \sigma^2) random variable.
    !        NORMAL2 = (PI*DT/(0.2*2.*perio)) * sqrt(-2.*log(UNIFORM1))*sin(2.*PI*UNIFORM2) ! Box-Muller method to generate a 2nd N(0, \sigma^2) random variable.
    !        XPHASE_RANDOMWALK = XPHASE_RANDOMWALK + NORMAL1
    !        TPHASE_RANDOMWALK = TPHASE_RANDOMWALK + NORMAL2
    !        PHASE_RANDOMWALK_LASTTIME = timelocal
    !      endif
    !      if(abs(x)<MICROBAROM_RANGE) then
    !        veloc_z_DG_P = MICROBAROM_AMPLITUDE & ! Amplitude.
    !        * sin(2.*PI*x/lambdo+XPHASE_RANDOMWALK) &
    !        !* 0.25*(1.-erf((x-MICROBAROM_RANGE+10.*lambdo)/(5.*lambdo)))*& ! Spatial apodisation. TODO: find why intrinsic functions won't compile.
    !        !       (1.+erf((x+MICROBAROM_RANGE-10.*lambdo)/(5.*lambdo)))& ! Spatial apodisation, continued.
    !        * sin(2.*PI*timelocal/perio+TPHASE_RANDOMWALK) ! &
    !        !* 0.5*(1.-erf((timelocal-MICROBAROM_MAXTIME+perio)/(0.5*perio))) ! Temporal apodisation. TODO: find why intrinsic functions won't compile.
    !        ! Spatial apodisation over 10 periods.
    !        if(abs(x)>MICROBAROM_RANGE-10.*lambdo) then
    !          veloc_z_DG_P = veloc_z_DG_P * (1.-((-abs(x)+MICROBAROM_RANGE)/(10.*lambdo)-1.)**2.)
    !        endif
    !        ! Temporal beginning apodisation over 1 period (split into 0.5 and 0.5 for smoothness).
    !        if(timelocal<0.5*perio) then
    !          veloc_z_DG_P = veloc_z_DG_P * (2.*(timelocal/perio)**2.)
    !        endif
    !        if(timelocal>=0.5*perio .and. timelocal<=perio) then
    !          veloc_z_DG_P = veloc_z_DG_P * ((4.*timelocal)/perio-2.*(timelocal/perio)**2.-1.)
    !        endif
    !        ! Temporal end apodisation over 1.5 periods.
    !        if(timelocal>MICROBAROM_MAXTIME-1.5*perio) then
    !          veloc_z_DG_P = veloc_z_DG_P * (1.-((timelocal-MICROBAROM_MAXTIME)/(1.5*perio)+1.)**2.)
    !        endif
    !      endif
    !    endif
    !  endif
    !  
    !  ! Plane wave forcing.
    !  if(.false.) then
    !    lambdo = main_spatial_period
    !    perio  = main_time_period
    !    if(timelocal<2. .and. abs(x)<20.d3) then
    !      veloc_z_DG_P = sin(2.*PI*timelocal/perio)
    !    endif
    !  endif
    !  
    !  ! Output a warning.
    !  if(timelocal == 0 .and. myrank==0) then
    !    write(*,*) "********************************"
    !    write(*,*) "*           WARNING            *"
    !    write(*,*) "********************************"
    !    write(*,*) "* A hardcoded bottom forcing   *"
    !    write(*,*) "* is being used. Use at your   *"
    !    write(*,*) "* own risk. See                *"
    !    write(*,*) "* 'boundary_terms_DG.f90'.     *"
    !    write(*,*) "********************************"
    !    call flush_IMAIN()
    !  endif
    !endif ! Endif on TYPEFORCING==9
    !! --------------------------- !
    !! Forcing read from file.     !
    !! --------------------------- !
    !! The Matlab script 'utils_new/forcings.m' can be used to generate the external bottom forcing file.
    !if(TYPE_FORCING == 10) then
    !  !write(*,*) "KEK", floor(timelocal/DT+1)
    !  !write(*,*) ibool_before_perio(i,j,ispec)
    !  !write(*,*) EXTFORC_MAP_ibbp_TO_LOCAL(1)
    !  if(timelocal<EXTERNAL_FORCING_MAXTIME) then
    !    externalforcingid=EXTFORC_MAP_ibbp_TO_LOCAL(ibool_before_perio(i,j,ispec))
    !    if(externalforcingid/=HUGE(0)) then
    !      ! Getting time step ID:
    !      veloc_z_DG_P=EXTERNAL_FORCING(floor(timelocal/EXTFORC_FILEDT+1),externalforcingid)
    !      ! Note: with tests done during loading, EXTFORC_FILEDT is either DT or DT/stage_time_scheme, which is consistent with the line above.
    !      !if(abs(timelocal-3.5)<0.1) then
    !      !  write(*,*) veloc_z_DG_P, rho_DG_P
    !      !endif
    !      !if(abs(timelocal-0.16)<0.2 .and. abs(abs(x)-50.)<0.4) then
    !      !  write(*,*) timelocal, x, veloc_z_DG_P
    !      !endif
    !    else
    !      if(x>=EXTFORC_MINX .and. x<=EXTFORC_MAXX) then
    !        ! This block is entered if forcing should happen (t < max_t and min_x < x < max_x), but can't be read (externalforcingid==HUGE(0)).
    !        ! Prompt an error.
    !        write(*,*) "********************************"
    !        write(*,*) "*            ERROR             *"
    !        write(*,*) "********************************"
    !        write(*,*) "* External bottom forcing      *"
    !        write(*,*) "* should happen, but does not  *"
    !        write(*,*) "* happen. Maybe some points    *"
    !        write(*,*) "* were not paired well, or     *"
    !        write(*,*) "* time to index conversion     *"
    !        write(*,*) "* fails.                       *"
    !        write(*,*) "* x_min = ", EXTFORC_MINX
    !        write(*,*) "* x     = ", x
    !        write(*,*) "* x_max = ", EXTFORC_MAXX
    !        write(*,*) "* t     = ", timelocal
    !        write(*,*) "* t_max = ", EXTERNAL_FORCING_MAXTIME
    !        write(*,*) "********************************"
    !        stop
    !      endif
    !    endif
    !  endif
    !endif ! Endif on TYPEFORCING==10
  endif

  ! Set energy.
  E_DG_P = p_DG_P/(gammaext_DG(ibool_DG(i, j, ispec)) - ONEl) &
                      + rho_DG_P*HALFl*( veloc_x_DG_P**2 + veloc_z_DG_P**2 ) 

  ! Set momenta.
  rhovx_DG_P = rho_DG_P*veloc_x_DG_P
  rhovz_DG_P = rho_DG_P*veloc_z_DG_P
   
end subroutine boundary_condition_DG

! ------------------------------------------------------------ !
! forcing_DG                                                   !
! ------------------------------------------------------------ !
! Computes a forcing.
! IN:
!   i, j, ispec, indices of the considered point
!   current_time, current time
! OUT:
!   forced_SF, forced scalar field
   
subroutine forcing_DG(i, j, ispec, current_time, forced_SF)
  use constants ! TODO: select variables to use.
  use specfem_par ! TODO: select variables to use.
  use specfem_par_LNS ! TODO: select variables to use.
  
  implicit none
  
  ! Input/Output.
  integer, intent(in) :: i, j, ispec
  real(kind=CUSTOM_REAL), intent(in) :: current_time
  real(kind=CUSTOM_REAL), intent(out) :: forced_SF
  
  ! Local.
  real(kind=CUSTOM_REAL), dimension(NDIM) :: x
  real(kind=CUSTOM_REAL) :: L_0,P_0,t_0,x_0
  integer externalforcingid ! TYPE_FORCING==10.
  
  if(TYPE_FORCING==0) return ! Safeguard in case routine is accidentally called, to prevent useless computations.
  
  ! Get point position.
  x = real(coord(:, ibool_before_perio(i, j, ispec)), kind=CUSTOM_REAL)
  x(NDIM) = x(NDIM) - coord_interface ! Remove eventual non-zero altitude of interface.
  
  if(abs(x(NDIM))<TINYVAL) then
    ! For now, forcings only occur on bottom boundary. This test thus enables faster treatment.
    
    L_0 = main_spatial_period ! TODO: replace "_0" variable by its longer name. Reduces readability, but will increase computation speed.
    P_0 = main_time_period ! TODO: replace "_0" variable by its longer name. Reduces readability, but will increase computation speed.
    t_0 = forcing_initial_time ! TODO: replace "_0" variable by its longer name. Reduces readability, but will increase computation speed.
    x_0 = forcing_initial_loc ! TODO: replace "_0" variable by its longer name. Reduces readability, but will increase computation speed.
    
    select case (TYPE_FORCING)
      case (99)
        ! Testing purposes case.
        if(abs(x(NDIM))<TINYVAL) then
          ! If z==0.
          !forced_SF = sin(2.*PI*x(1)/10.)
          forced_SF = cos(2.*PI*current_time/0.08)-1.
        else
          ! Leave everything unchanged and leave early.
          return
        endif
      
      case (1)
        ! Time Gaussian derivative (acoustic plane wave forcing).
        ! Old code, which is actually a Gaussian second derivative.
        !forced_SF = 0.01*(&
        !               - (2d0/(P_0/4d0))*((current_time-(t_0-P_0/4d0))/(P_0/4d0))* &
        !                     (exp(-((current_time-(t_0-P_0/4d0))/(P_0/4d0))**2)) &
        !               + (2d0/(P_0/4d0))*((current_time-(t_0+P_0/4d0))/(P_0/4d0))* &
        !                     (exp(-((current_time-(t_0+P_0/4d0))/(P_0/4d0))**2)) )
        ! Actual Gaussian derivative.
        forced_SF = -(2.*(PI/P_0)**2) * (current_time-t0-t_0) * exp(-((current_time-t0-t_0)*PI/P_0)**2)
      
      case (2)
        ! Time and space Gaussian derivative (gravity wave forcing).
        ! Old code. TODO: recenter peak by taking starting delay ("negative times") into account.
        forced_SF = 0.001*(&
                      - (2d0/(P_0/4d0))*((current_time-(t_0-P_0/4d0))/(P_0/4d0))* &
                               (exp(-((current_time-(t_0-P_0/4d0))/(P_0/4d0))**2)) &
                      + (2d0/(P_0/4d0))*((current_time-(t_0+P_0/4d0))/(P_0/4d0))* &
                               (exp(-((current_time-(t_0+P_0/4d0))/(P_0/4d0))**2)) ) &
                       * ( exp(-((x(1)-(x_0-L_0/4))/(L_0/4))**2) - &
                              exp(-((x(1)-(x_0+L_0/4))/(L_0/4))**2) )
      
      case (3)
        ! Gaussian second derivative.
        forced_SF = (1.-2.*(PI*(current_time-t0-t_0)/P_0)**2) * exp(-((current_time-t0-t_0)*PI/P_0)**2)
      
      case (9)
        ! Hardcoded forcing.
        
        !if(.false.) then ! Tsunami forcing.
        !  VELOC_TSUNAMI = 200.
        !  P_0 = 1000.
        !  L_0 = 50000.
        !  if(timelocal < P_0) then
        !      forced_SF = (1d0/P_0)*exp( -((x-L_0)/(sqrt(2d0)*P_0))**2 )
        !  else
        !      forced_SF = &
        !              2d0*((VELOC_TSUNAMI)/(sqrt(2d0)*P_0))*(((x-L_0)-VELOC_TSUNAMI*(current_time - P_0))/(sqrt(2d0)*P_0))&
        !              *exp( -(((x-L_0)-VELOC_TSUNAMI*(current_time - P_0))/(sqrt(2d0)*P_0))**2 )
        !  endif
        !endif
        
        !if(.true.) then ! Analytic microbarom forcing.
        !  MICROBAROM_AMPLITUDE = 1. ! Microbarom amplitude. Unit: m.
        !  MICROBAROM_RANGE = 80.*L_0 ! Range around x=0 to which impose microbaroms. Unit: m. Be careful with apodisation below.
        !  MICROBAROM_MAXTIME = 10.5 * P_0 ! Microbarom active from t=0 to t=MICROBAROM_MAXTIME. Unit: s. Be careful with apodisation below.
        !  if(timelocal==0) then
        !    ! Start random phase walk.
        !    XPHASE_RANDOMWALK = 0.
        !    TPHASE_RANDOMWALK = 0.
        !    PHASE_RANDOMWALK_LASTTIME=0.
        !  endif
        !  if(timelocal<MICROBAROM_MAXTIME) then
        !    if(timelocal>=PHASE_RANDOMWALK_LASTTIME+DT) then
        !      ! Update the random walk only once par time step.
        !      call random_number(UNIFORM1)
        !      call random_number(UNIFORM2)
        !      NORMAL1 = (PI*DT/(0.2*2.*P_0)) * sqrt(-2.*log(UNIFORM1))*cos(2.*PI*UNIFORM2) ! Box-Muller method to generate a 1nd N(0, \sigma^2) random variable.
        !      NORMAL2 = (PI*DT/(0.2*2.*P_0)) * sqrt(-2.*log(UNIFORM1))*sin(2.*PI*UNIFORM2) ! Box-Muller method to generate a 2nd N(0, \sigma^2) random variable.
        !      XPHASE_RANDOMWALK = XPHASE_RANDOMWALK + NORMAL1
        !      TPHASE_RANDOMWALK = TPHASE_RANDOMWALK + NORMAL2
        !      PHASE_RANDOMWALK_LASTTIME = timelocal
        !    endif
        !    if(abs(x)<MICROBAROM_RANGE) then
        !      forced_SF = MICROBAROM_AMPLITUDE & ! Amplitude.
        !      * sin(2.*PI*x/L_0+XPHASE_RANDOMWALK) &
        !      !* 0.25*(1.-erf((x-MICROBAROM_RANGE+10.*L_0)/(5.*L_0)))*& ! Spatial apodisation. TODO: find why intrinsic functions won't compile.
        !      !       (1.+erf((x+MICROBAROM_RANGE-10.*L_0)/(5.*L_0)))& ! Spatial apodisation, continued.
        !      * sin(2.*PI*timelocal/P_0+TPHASE_RANDOMWALK) ! &
        !      !* 0.5*(1.-erf((timelocal-MICROBAROM_MAXTIME+P_0)/(0.5*P_0))) ! Temporal apodisation. TODO: find why intrinsic functions won't compile.
        !      ! Spatial apodisation over 10 periods.
        !      if(abs(x)>MICROBAROM_RANGE-10.*L_0) then
        !        forced_SF = forced_SF * (1.-((-abs(x)+MICROBAROM_RANGE)/(10.*L_0)-1.)**2.)
        !      endif
        !      ! Temporal beginning apodisation over 1 period (split into 0.5 and 0.5 for smoothness).
        !      if(timelocal<0.5*P_0) then
        !        forced_SF = forced_SF * (2.*(timelocal/P_0)**2.)
        !      endif
        !      if(timelocal>=0.5*P_0 .and. timelocal<=P_0) then
        !        forced_SF = forced_SF * ((4.*timelocal)/P_0-2.*(timelocal/P_0)**2.-1.)
        !      endif
        !      ! Temporal end apodisation over 1.5 periods.
        !      if(timelocal>MICROBAROM_MAXTIME-1.5*P_0) then
        !        forced_SF = forced_SF * (1.-((timelocal-MICROBAROM_MAXTIME)/(1.5*P_0)+1.)**2.)
        !      endif
        !    endif
        !  endif
        !endif
        
        ! Output a warning.
        if(current_time == 0 .and. myrank==0) then
          write(*,*) "********************************"
          write(*,*) "*           WARNING            *"
          write(*,*) "********************************"
          write(*,*) "* A hardcoded bottom forcing   *"
          write(*,*) "* is being used. Use at your   *"
          write(*,*) "* own risk. See                *"
          write(*,*) "* 'boundary_terms_DG.f90'.     *"
          write(*,*) "********************************"
          call flush_IMAIN()
        endif
      
      case (10)
        ! --------------------------- !
        ! Forcing read from file.     !
        ! --------------------------- !
        ! The Matlab script 'utils_new/forcings.m' can be used to generate the external bottom forcing file.
        !write(*,*) "KEK", floor(timelocal/DT+1)
        !write(*,*) ibool_before_perio(i,j,ispec)
        !write(*,*) EXTFORC_MAP_ibbp_TO_LOCAL(ibool_before_perio(i,j,ispec)), i, j, ispec ! DEBUG
        if(current_time<EXTERNAL_FORCING_MAXTIME) then
          externalforcingid=EXTFORC_MAP_ibbp_TO_LOCAL(ibool_before_perio(i,j,ispec))
          if(externalforcingid/=HUGE(0)) then
            ! If there exists an id, retrieve forcing:
            forced_SF=EXTERNAL_FORCING(floor(current_time/EXTFORC_FILEDT+1),externalforcingid)
          else ! Else on externalforcingid.
            if(x(1)>=EXTFORC_MINX .and. x(1)<=EXTFORC_MAXX) then
              ! This block is entered if forcing should happen (t < max_t and min_x < x < max_x), but can't be read (externalforcingid==HUGE(0)).
              ! Prompt an error.
              write(*,*) "********************************"
              write(*,*) "*            ERROR             *"
              write(*,*) "********************************"
              write(*,*) "* External bottom forcing      *"
              write(*,*) "* should happen, but does not  *"
              write(*,*) "* happen. Maybe some points    *"
              write(*,*) "* were not paired well, or     *"
              write(*,*) "* time to index conversion     *"
              write(*,*) "* fails.                       *"
              write(*,*) "* x_min = ", EXTFORC_MINX, "(included)"
              write(*,*) "* x     = ", x
              write(*,*) "* x_max = ", EXTFORC_MAXX, "(included)"
              write(*,*) "* t     = ", current_time
              write(*,*) "* t_max = ", EXTERNAL_FORCING_MAXTIME, "(excluded)"
              write(*,*) "********************************"
              stop
            endif
          endif ! Endif on externalforcingid.
        else ! Else on current_time.
          ! Leave everything unchanged and leave early.
          return
        endif ! Endif on current_time.
      
      case (0)
        ! Safeguard in case routine is accidentally called, to prevent useless computations.
        return
      
      case default
        ! TYPE_FORCING not implemented.
        write(*,*) "********************************"
        write(*,*) "*            ERROR             *"
        write(*,*) "********************************"
        write(*,*) "* Forcing unknown or not yet   *"
        write(*,*) "* implemented:                 *"
        write(*,*) "* TYPE_FORCING = ", TYPE_FORCING
        write(*,*) "********************************"
        stop
    end select ! Endselect on TYPEFORCING.
  else ! Else on abs(z).
    ! Leave everything unchanged and leave early.
    return
  endif ! Endif on abs(z).
end subroutine forcing_DG

! ------------------------------------------------------------ !
! prepare_external_forcing                                     !
! ------------------------------------------------------------ !
! Reads EXTERNAL_FORCING_FILENAME, save values of bottom forcing in an array, and map mesh coordinates to that array.

subroutine prepare_external_forcing()

  use constants,only: CUSTOM_REAL, HUGEVAL, TINYVAL, NDIM

  use specfem_par, only: EXTERNAL_FORCING_FILENAME, EXTERNAL_FORCING_MAXTIME,&
                         EXTFORC_MINX, EXTFORC_MAXX, &
                         EXTERNAL_FORCING, EXTFORC_MAP_ibbp_TO_LOCAL, DT,&
                         EXTFORC_FILEDT,&
                         nglob_DG,coord,nspec,NGLLX,NGLLZ,myrank,&!NPROC,&
                         only_DG_acoustic,ibool_before_perio,stage_time_scheme,&
                         ADD_PERIODIC_CONDITIONS

  implicit none

  ! Local variables.
  logical :: fileexists, counting_nx
  real(kind=CUSTOM_REAL) :: t, x, val, tmp_t_old
  real(kind=CUSTOM_REAL), dimension(NDIM) :: Xspcfm
  integer :: io, istat,NSTEPFORCING,nx,it,ix,ibbp,ispec,i,j,nx_paired,nx_paired_tot
  
  if(.not. only_DG_acoustic) then
    write(*,*) "********************************"
    write(*,*) "*            ERROR             *"
    write(*,*) "********************************"
    write(*,*) "* External bottom forcing can  *"
    write(*,*) "* only be used on a full DG    *"
    write(*,*) "* simulation, which is now not *"
    write(*,*) "* the case.                    *"
    write(*,*) "********************************"
    stop
  endif
  ! Check existence of model file.
  fileexists=.false.
  INQUIRE(File=EXTERNAL_FORCING_FILENAME, Exist=fileexists)
  if(.not. fileexists) then
    write(*,*) "********************************"
    write(*,*) "*            ERROR             *"
    write(*,*) "********************************"
    write(*,*) "* external_bottom_forcing file *"
    write(*,*) "* does not exist in folder.    *"
    write(*,*) "********************************"
    stop
  endif
  
  ! Read once first to find out forcing stopping time, number of points at z=0, and associate each index to a iglob_before_perio.
  ! TODO: There must be a way to do this more efficiently.
  allocate(EXTFORC_MAP_ibbp_TO_LOCAL(nglob_DG),stat=istat)
  if(istat/=0) then
    write(*,*) "********************************"
    write(*,*) "*            ERROR             *"
    write(*,*) "********************************"
    write(*,*) "* Allocation of                *"
    write(*,*) "* EXTFORC_MAP_ibbp_TO_LOCAL    *"
    write(*,*) "* has failed.                  *"
    write(*,*) "********************************"
    stop
  endif
  if(ADD_PERIODIC_CONDITIONS) then
    write(*,*) "********************************"
    write(*,*) "*           WARNING            *"
    write(*,*) "********************************"
    write(*,*) "* Periodic boundary conditions *"
    write(*,*) "* are activated. External      *"
    write(*,*) "* bottom forcing can be used   *"
    write(*,*) "* with those, but it must be   *"
    write(*,*) "* prepared beforehand.         *"
    write(*,*) "********************************"
  endif
  EXTFORC_MAP_ibbp_TO_LOCAL=HUGE(0) ! Safeguard.
  
  if(myrank==0) then
    write(*,*) "  First scan: counting number of points, figuring file time step, and pairing file points to SPECFEM points."
  endif
  tmp_t_old=0.
  counting_nx=.true.
  nx=0 ! Counter for number of spatial points found in file.
  nx_paired=0 ! Counter for points actually paired with their SPECFEM-DG counterparts. Safeguard.
  EXTFORC_MINX=HUGEVAL ! Start at +inf, will be updated
  EXTFORC_MAXX=-HUGEVAL ! Start at -inf, will be updated
  OPEN(100, file=EXTERNAL_FORCING_FILENAME)
  DO
    READ(100,*,iostat=io) t, x, val ! In fact, we do not care about val here.
    !write(*,*) t, tmp_t_old
    if(counting_nx .and. t>tmp_t_old) then
      !nx=nx+1 ! Add last one.
      counting_nx=.false. ! Deactivate counting.
      EXTFORC_FILEDT=t-tmp_t_old ! Save file dt.
    endif
    if(counting_nx) then ! This if is entered only for first time step (see else below).
      nx=nx+1 ! Count.
      ! Find the ibool_before_perio corresponding to the current x.
      do ispec = 1, nspec
        ! For DG elements, go through GLL points one by one.
        do j = 1, NGLLZ
          do i = 1, NGLLX
            ibbp=ibool_before_perio(i, j, ispec)
            Xspcfm = coord(:, ibbp) ! Get coordinates of GLL point.
            if(abs(Xspcfm(NDIM))>TINYVAL) then ! If z is too far away from 0, cycle point.
              cycle
            else
              !if(.true.) then! .and. abs(Xspcfm(1))<650.) then ! DEBUG
              !  if(ispec==nspec) stop
              !  write(*,*) Xspcfm(1)
              !endif
              !write(*,*) "abs(x-coord(1, ibbp))", abs(x-coord(1, ibbp)) ! DEBUG
              !if(.true. .and. abs(x-0.)<650. .and. abs(Xspcfm(1)-0.)<650.) then ! DEBUG
              !  write(*,*) t, x, Xspcfm(1), abs(x-Xspcfm(1)),&
              !             nx, nx_paired
              !endif
              if(     abs(x-Xspcfm(1))<1e-6 &
                 .or. abs(x-Xspcfm(1))/(abs(x)+abs(Xspcfm(1)))<1e-5) then ! 1e-6 is a relative empiric threshold based on precision we can obtain by generating forcings with our Matlab script. Basically, we chose '%.8e' as printing format for x (in the script), which means we have 8 significant digits, thus theoretically we have a relative precision of 1e-8. We chose the threshold to be 1e-6 in order to a have a safety margin. It is acceptable because we expect the GLL points never to be (relatively) 1e-7 apart.
                EXTFORC_MAP_ibbp_TO_LOCAL(ibbp)=nx
                !write(*,*) EXTFORC_MAP_ibbp_TO_LOCAL(ibbp) ! DEBUG
                nx_paired=nx_paired+1
              endif
            endif
          enddo ! Enddo on i.
        enddo ! Enddo on j.
      enddo ! Enddo on ispec.
    endif ! Endif on counting_nx.
    if(x<=EXTFORC_MINX) then
      EXTFORC_MINX=x
    endif
    if(x>=EXTFORC_MAXX) then
      EXTFORC_MAXX=x
    endif
    tmp_t_old=t
    IF (io/=0) EXIT
  ENDDO
  close(100)
  EXTERNAL_FORCING_MAXTIME=t ! Save maximum time.
  
  !write(*,*) "nx_paired, nx_paired_tot", nx_paired, nx_paired_tot
    
  if(myrank==0) then
    write(*,*) "  File:", nx, " mesh points were found."
  endif
  write(*,*) "  CPU ", myrank, ": ", nx_paired, " SPECFEM mesh points were paired."
  
  call sum_all_all_i(nx_paired, nx_paired_tot)
  if(myrank==0) then
    write(*,*) "  Across CPUs:       ", nx_paired_tot, " SPECFEM mesh points were paired."
    write(*,*) "  Because of the DG implementation, duplicates can occur without problem."
  endif

  if(myrank==0 .and. nx_paired_tot<nx) then
    ! TODO: work on that condition. For instance, find exactly the number of expected duplicates.
    write(*,*) "********************************"
    write(*,*) "*            ERROR             *"
    write(*,*) "********************************"
    write(*,*) "* Not all spatial points read  *"
    write(*,*) "* from external file were      *"
    write(*,*) "* paired with SPECFEM          *"
    write(*,*) "* counterparts:                *"
    write(*,*) "* proc      = ", myrank
    write(*,*) "* nx        = ", nx
    write(*,*) "* nx_paired = ", nx_paired_tot
    write(*,*) "* Consider recompiling the     *"
    write(*,*) "* external file with matching  *"
    write(*,*) "* spatial mesh (see            *"
    write(*,*) "* 'utils_new/forcings.m').     *"
    write(*,*) "* Using mesh points directly   *"
    write(*,*) "* read from databases could    *"
    write(*,*) "* also be a safer method. See  *"
    write(*,*) "* also the import routine,     *"
    write(*,*) "* 'prepare_external_forcing'   *"
    write(*,*) "* in 'boundary_terms_DG.f90'.  *"
    write(*,*) "********************************"
    stop
  endif ! Endif on nx_paired and nx.
  
  !write(*,*) "observed dt from file: ", file_dt, DT/stage_time_scheme ! DEBUG
  if(EXTFORC_FILEDT==DT) then ! Normal sampling (ok files).
    NSTEPFORCING = int(floor(EXTERNAL_FORCING_MAXTIME/DT+2))
  else
    write(*,*) "********************************"
    write(*,*) "*            ERROR             *"
    write(*,*) "********************************"
    write(*,*) "* File time step does not      *"
    write(*,*) "* match simulation time step   *"
    write(*,*) "* or simulation time substeps  *"
    write(*,*) "* (RK substeps).               *"
    write(*,*) "* file_dt              = ", EXTFORC_FILEDT
    write(*,*) "* DT                   = ", DT
    write(*,*) "* DT/stage_time_scheme = ", DT/stage_time_scheme
    write(*,*) "********************************"
    stop
  endif
  
  !write(*,*) t, DT, NSTEPFORCING, EXTFORC_MINX, EXTFORC_MAXX, nx,stage_time_scheme ! DEBUG
  
  allocate(EXTERNAL_FORCING(NSTEPFORCING, nx),stat=istat)
  if(istat/=0) then
    write(*,*) "********************************"
    write(*,*) "*            ERROR             *"
    write(*,*) "********************************"
    write(*,*) "* Allocation of                *"
    write(*,*) "* EXTERNAL_FORCING has failed. *"
    write(*,*) "********************************"
    stop
  endif
  
  ! Re-read and save values.
  if(myrank==0) then
    write(*,*) "  Second scan: saving values of forcing for each point and each time step."
  endif
  it=1
  ix=1
  OPEN(100, file=EXTERNAL_FORCING_FILENAME)
  DO
    READ(100,*,iostat=io) t, x, val
    if(.false. .and. it>1245) then
      write(*,*) it, t, ix, x, val ! DEBUG
    endif
    EXTERNAL_FORCING(it, ix)=val
    ix=ix+1
    if(mod(ix,nx+1)==0) then
      ix=1
      it=it+1
    endif
    IF (io/=0) EXIT
  ENDDO
  close(100)
  
  ! DEBUG
  !write(*,*) EXTERNAL_FORCING(floor((DT+2*DT/stage_time_scheme)/DT+1),:)
  !if(.false.) then
  !  do ispec = 1, nspec
  !    do j = 1, NGLLZ
  !      do i = 1, NGLLX
  !        ibbp=ibool_before_perio(i, j, ispec)
  !        z = coord(2, ibbp)
  !        if(z/=0.) then
  !          cycle
  !        else
  !          write(*,*) coord(1,ibbp), z, ibbp, EXTFORC_MAP_ibbp_TO_LOCAL(ibbp)
  !        endif
  !      enddo ! Enddo on i.
  !    enddo ! Enddo on j.
  !  enddo ! Enddo on ispec.
  !endif
  !if(.false.) then
  !  write(*,*) "this must be 0"
  !  write(*,*) EXTERNAL_FORCING(1,:)
  !  write(*,*) "this must be non 0"
  !  write(*,*) EXTERNAL_FORCING(2,:)
  !  write(*,*) "this must be non 0"
  !  write(*,*) EXTERNAL_FORCING(3,:)
  !endif
  !stop
end subroutine prepare_external_forcing

! ------------------------------------------------------------ !
! compute_interface_unknowns                                   !
! ------------------------------------------------------------ !
! From the coordinates of a GLL point in an element (local coordinates (i, j) and global element number ispec), and its neighbour's identifier (neighbor), compute the values of the constitutive variables at the neighbour.
! Variables ending in "_P" (for "plus") are output-intended and are the sought values, or exterior values.
! Variables ending in "_iM" (for "minus") are input-intended and should correspond to the interior values.
! Variables ending in "_iP" are input-intended, and correspond to the exterior values, if known. Remark that those variables are only used if neighbor(3) != -1.
! Variable "exact_interface_flux" is output-intended. If set to .true., an exact flux formula will be used. If set to .false., the Lax-Friedrich approximation for the flux will be used.
 
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
                         p_DG_init, gammaext_DG, muext, etaext, kappa_DG, ibool, c_V, &
                         buffer_DG_gamma_P, coord,  &
                         rho_init, rhovx_init, rhovz_init, E_init, &
                         potential_dot_dot_acoustic, &
                         veloc_vector_acoustic_DG_coupling, MPI_transfer_iface,&
                         ibool_before_perio, ABC_STRETCH, ABC_STRETCH_LEFT,&
                         ABC_STRETCH_RIGHT, ABC_STRETCH_LEFT_LBUF, ABC_STRETCH_RIGHT_LBUF,&
                         mesh_xmin, mesh_xmax
                         
  implicit none
  
  integer, intent(in) :: i, j, ispec, iface1, iface
  
  integer, dimension(3), intent(in) :: neighbor
  
  real(kind=CUSTOM_REAL), intent(out) :: rho_DG_P, rhovx_DG_P, rhovz_DG_P, &
                                         E_DG_P, veloc_x_DG_P, veloc_z_DG_P, p_DG_P, T_P, &
                                         Tx_DG_P, Tz_DG_P, Vxx_DG_P, Vzz_DG_P, Vzx_DG_P, Vxz_DG_P
        
  real(kind=CUSTOM_REAL), intent(in) :: timelocal
        
  logical, intent(out) :: exact_interface_flux
  
  real(kind=CUSTOM_REAL), intent(in) :: nx, nz
  
  real(kind=CUSTOM_REAL) :: tx, tz, normal_v, tangential_v, &
                            veloc_x, veloc_z, weight, &
                            veloc_x_DG_iM, veloc_z_DG_iM, p_DG_iM, gamma_P, e1_DG_P
        
  real(kind=CUSTOM_REAL), intent(in) :: rho_DG_iM, E_DG_iM, rhovx_DG_iM, rhovz_DG_iM, &
                                        rho_DG_iP, E_DG_iP, rhovx_DG_iP, rhovz_DG_iP
  real(kind=CUSTOM_REAL), dimension(2), intent(in) :: T_DG_iP, T_DG_iM
  real(kind=CUSTOM_REAL), dimension(2, 2), intent(in) :: V_DG_iP, V_DG_iM
  
  ! Local variables.
  integer :: iglobM, i_el, j_el, ispec_el, i_ac, j_ac, ispec_ac, iglob, iglobP, ipoin, num_interface
  real(kind=CUSTOM_REAL), dimension(2, 2) :: trans_boundary
  
  ! Parameters
  real(kind=CUSTOM_REAL), parameter :: ZERO = 0._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: ONE  = 1._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: TWO  = 2._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: HALF = 0.5_CUSTOM_REAL
  
  real(kind=CUSTOM_REAL) :: x ! For coupling deactivation in buffers.
  
  ! Characteristic based BC
  real(kind=CUSTOM_REAL) :: p_b, rho_b, s_b, rho_inf, v_b_x, v_b_z, c_b, un_b, &
                            rlambda_max, rlambda_min, c_in, c_inf, un_in, un_inf, &
                            deltaZ1, deltaZ2star, p_n, a_n, alpha0
  
  iglobM = ibool_DG(i, j, ispec)
  
  if(.false.) weight=1. ! Horrible hack to get rid of the warning "Warning: Unused dummy argument 'weight'". I'm sorry. TODO: Remove weight from arguments in all calls to compute_interface_unknowns.
  
  ! Deduce velocities from momenta and density (the two latter being constitutive variables, and thus passed to the routine as input parameters).
  veloc_x_DG_iM = rhovx_DG_iM/rho_DG_iM
  veloc_z_DG_iM = rhovz_DG_iM/rho_DG_iM
  ! Deduce pressure from energy (consitutive), density (consitutive), and velocities (deduced above).
  p_DG_iM       = (gammaext_DG(iglobM) - ONE)*( E_DG_iM &
        - (HALF)*rho_DG_iM*( veloc_x_DG_iM**2 + veloc_z_DG_iM**2 ) )
  
  ! Extract derivatives.
  Tx_DG_P  = -T_DG_iM(1)
  Tz_DG_P  = -T_DG_iM(2)
  Vxx_DG_P = -V_DG_iM(1, 1)
  Vzz_DG_P = -V_DG_iM(2, 2)
  Vxz_DG_P = -V_DG_iM(1, 2)
  Vzx_DG_P = -V_DG_iM(2, 1)
  
  gamma_P = gammaext_DG(iglobM)
  
  e1_DG_P = ZERO
  
  exact_interface_flux = .false. ! Unless otherwise specified, use the Lax-Friedrich approximation.
  
  if(neighbor(3) == -1 ) then
    ! --------------------------- !
    ! neighbor(3) == -1, edge is  !
    ! an 'outside' edge.          !
    ! --------------------------- !
    ! That is, either its corresponding edge is in another partition, or the edge is a coupling (with another material) edge (also in another partition, by construction), or the edge is on the outer boundary of the computational domain.
    
    ipoin         = -1 ! By default, specify that for the triplet (iface1,iface,ispec), the values should not be sought in another partition. That is, either the edge is a coupling (with another material) edge, or an edge on the outer boundary of the computational domain.
    num_interface = -1 ! Initialised to this value, but should not be used anywhere as is.
    if(NPROC > 1) then
      ipoin         = MPI_transfer_iface(iface1, iface, ispec, 1)
      num_interface = MPI_transfer_iface(iface1, iface, ispec, 2)
    endif                  
    
    if(ipoin > -1) then
      ! --------------------------- !
      ! ipoin > -1, values should   !
      ! be sought in another        !
      ! partition, thus on an       !
      ! interface, thus DG buffers  !
      ! should be used.             !
      ! --------------------------- !
      
      if(ispec_is_acoustic_coupling_ac(ibool_DG(i, j, ispec)) >= 0 .AND. .false.) then
        ! --------------------------- !
        ! MPI acoustic potential      !
        ! neighbour.                  !
        ! --------------------------- !
        ! TODO: Decide what to do with this case (we are in an 'if(.false.)').
        
        ! We already know the "real" flux at boundary.
        exact_interface_flux = .true.

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

        T_P = (E_DG_iM/rho_DG_iM - 0.5*(veloc_x_DG_iM**2 + veloc_z_DG_iM**2))/c_V
        
      else ! Happens everytime (since 'if' above is an 'if(.false.)').
        ! --------------------------- !
        ! MPI acoustic DG neighbour   !
        ! (not acoustic potential     !
        ! neighbour).                 !
        ! --------------------------- !
        
        ! Get all values from the MPI buffers.
        rho_DG_P     = buffer_DG_rho_P(ipoin, num_interface)
        E_DG_P       = buffer_DG_E_P(ipoin, num_interface)
        rhovx_DG_P   = buffer_DG_rhovx_P(ipoin, num_interface)
        rhovz_DG_P   = buffer_DG_rhovz_P(ipoin, num_interface)
        gamma_P      = buffer_DG_gamma_P(ipoin,num_interface)
        if(muext(i, j, ispec) > 0 .OR. &
           etaext(i, j, ispec) > 0 .OR. &
           kappa_DG(i, j, ispec) > 0) then
          ! If there is viscosity, get the values from the MPI buffers.
          ! If not, values as initialised above are kept.
          Vxx_DG_P = buffer_DG_Vxx_P(ipoin, num_interface)
          Vzz_DG_P = buffer_DG_Vzz_P(ipoin, num_interface)
          Vxz_DG_P = buffer_DG_Vxz_P(ipoin, num_interface)
          Vzx_DG_P = buffer_DG_Vzx_P(ipoin, num_interface)
          Tx_DG_P  = buffer_DG_Tx_P(ipoin, num_interface)
          Tz_DG_P  = buffer_DG_Tz_P(ipoin, num_interface)
        endif
        
        ! Deduce velocities, pressure, and temperature.
        veloc_x_DG_P = rhovx_DG_P/rho_DG_P
        veloc_z_DG_P = rhovz_DG_P/rho_DG_P
        p_DG_P = (gamma_P - ONE)*( E_DG_P &
                 - (HALF)*rho_DG_P*( veloc_x_DG_P**2 + veloc_z_DG_P**2 ) )
        T_P = (E_DG_P/rho_DG_P - 0.5*((rhovx_DG_P/rho_DG_P)**2 + (rhovz_DG_P/rho_DG_P)**2))/c_V

      endif
                    
    elseif(ACOUSTIC_FORCING .AND. ispec_is_acoustic_forcing(i, j, ispec)) then
      ! --------------------------- !
      ! ipoin == -1                 !
      !   (+) acoustic forcing      !
      ! --------------------------- !
      write(*,*) "********************************"
      write(*,*) "*            ERROR             *"
      write(*,*) "********************************"
      write(*,*) "* ACOUSTIC_FORCING obsolete    *"
      write(*,*) "* for DG simulations.          *"
      write(*,*) "********************************"
      stop
         
    elseif(ispec_is_acoustic_coupling_el(i, j, ispec, 3) >= 0) then
      ! --------------------------- !
      ! ipoin == -1                 !
      !   (+) elastic coupling      !
      ! --------------------------- !
      ! TODO: Clean up all this part.
      
      exact_interface_flux = .false.

      ! Coordinates of elastic element
      i_el     = ispec_is_acoustic_coupling_el(i, j, ispec, 1)
      j_el     = ispec_is_acoustic_coupling_el(i, j, ispec, 2)
      ispec_el = ispec_is_acoustic_coupling_el(i, j, ispec, 3)

      iglob = ibool(i_el, j_el, ispec_el)
      
      call boundary_condition_DG(i, j, ispec, timelocal, rho_DG_P, rhovx_DG_P, rhovz_DG_P, E_DG_P, &
                                 veloc_x_DG_P, veloc_z_DG_P, p_DG_P, e1_DG_P)

      rho_DG_P = rho_DG_iM
      
      ! Get elastic medium velocities.
      veloc_x = veloc_elastic(1, iglob)
      veloc_z = veloc_elastic(2, iglob)
      ! Tangential vector. It is assumed that we only have the elastic media under the DG medium, hence we always have nz>=0.
      tx = -nz
      tz =  nx
      ! Get the normal velocity of the solid perturbation (normal continuity) and the tangential velocity of the fluid flow (free slip).
      normal_v     = (veloc_x*nx + veloc_z*nz)
      tangential_v = veloc_x_DG_P*tx + veloc_z_DG_P*tz
      !normal_v     = -(veloc_x_DG_iM*nx + veloc_z_DG_iM*nz) + 2*(veloc_x*nx + veloc_z*nz) ! DEBUG
      !tangential_v = veloc_x_DG_iM*tx + veloc_z_DG_iM*tz ! DEBUG
      !tangential_v    = (veloc_x*tx + veloc_z*tz) ! DEBUG
      ! Set the matrix for the transformation from normal/tangential coordinates to mesh coordinates.
      trans_boundary(1, 1) =  tz
      trans_boundary(1, 2) = -nz
      trans_boundary(2, 1) = -tx
      trans_boundary(2, 2) =  nx
      trans_boundary = trans_boundary/(nx * tz - tx * nz)
      
      ! QUICK HACK: DEACTIVATE COUPLING IN BUFFER ZONES.
      x = coord(1, ibool_before_perio(i, j, ispec))
      if(ABC_STRETCH .and. &
         (     (ABC_STRETCH_LEFT   .and. x < mesh_xmin + ABC_STRETCH_LEFT_LBUF) & ! left stretching and in left buffer zone ! TODO: use stretching_buffer variable
          .or. (ABC_STRETCH_RIGHT  .and. x > mesh_xmax - ABC_STRETCH_RIGHT_LBUF) & ! right stretching and in right buffer zone ! TODO: use stretching_buffer variable
         ) &
        ) then
        
        if(ABC_STRETCH_LEFT) x = (x - mesh_xmin)/ABC_STRETCH_LEFT_LBUF ! x is a local buffer coordinate now (1 at beginning, 0 at end).
        if(ABC_STRETCH_RIGHT) x = (mesh_xmax - x)/ABC_STRETCH_RIGHT_LBUF ! x is a local buffer coordinate now (1 at beginning, 0 at end).
        
        if(x<0.1 .and. .false.) then !DEBUG
          write(*,*) x, 'b4:', veloc_x_DG_P
        endif
        
        ! Convert (back) the velocity components from normal/tangential coordinates to mesh coordinates.
        veloc_x_DG_P = x*(trans_boundary(1, 1)*normal_v + trans_boundary(1, 2)*tangential_v)&
                       +(1.-x)*veloc_x_DG_P
        veloc_z_DG_P = x*(trans_boundary(2, 1)*normal_v + trans_boundary(2, 2)*tangential_v)&
                       +(1.-x)*veloc_z_DG_P
        ! This gradually deactivates coupling in the horizontal buffers.
        ! TODO: This is a hack. A better method is to be preferred.
        
        if(x<0.1 .and. .false.) then !DEBUG
          write(*,*) x, 'aftR:', veloc_x_DG_P, x, (1.-x)
        endif
      else
        ! Convert (back) the velocity components from normal/tangential coordinates to mesh coordinates.
        veloc_x_DG_P = trans_boundary(1, 1)*normal_v + trans_boundary(1, 2)*tangential_v
        veloc_z_DG_P = trans_boundary(2, 1)*normal_v + trans_boundary(2, 2)*tangential_v
      endif
      
      !if(abs(coord(1, ibool(i, j, ispec))) < 1. &
      !   .and. abs(coord(2, ibool(i, j, ispec))) < 1e-6 &
      !   .and. timelocal>0.78) then ! DEBUG
      !  write(*,*)  coord(2, ibool(i, j, ispec)), veloc_z, veloc_z_DG_P ! DEBUG
      !endif
      
      ! No stress continuity.
      p_DG_P = p_DG_iM

      ! Deduce energy.
      E_DG_P = p_DG_P/(gammaext_DG(iglobM) - 1.) + HALF*rho_DG_P*( veloc_x_DG_P**2 + veloc_z_DG_P**2 )

      ! Recompute momenta in order not to keep those coming from the call to boundary_condition_DG above.
      rhovx_DG_P = rho_DG_P*veloc_x_DG_P
      rhovz_DG_P = rho_DG_P*veloc_z_DG_P
      
      ! Deduce temperature.
      T_P = (E_DG_iM/rho_DG_iM - 0.5*(veloc_x_DG_iM**2 + veloc_z_DG_iM**2))/c_V
      
      !if(      coord(2,ibool_before_perio(i,j,ispec))<1. & ! DEBUG
      !   .and. coord(2,ibool_before_perio(i,j,ispec))>=0. & ! DEBUG
      !   .and. abs(coord(1,ibool_before_perio(i,j,ispec)))<2.) then ! DEBUG
      !  write(*,*) timelocal, coord(:,ibool_before_perio(i,j,ispec)), & ! DEBUG
      !             !p_DG_iM, p_DG_P ! DEBUG
      !             !rhovx_DG_P, rhovz_DG_P ! DEBUG
      !             !veloc_x_DG_P, veloc_z_DG_P ! DEBUG
      !             !trans_boundary ! DEBUG
      !endif ! DEBUG

    else
      ! --------------------------- !
      ! ipoin == -1                 !
      !   classical outer boundary  !
      !   conditions                !
      ! --------------------------- !
    
      exact_interface_flux = .false.

      ! Get "background" (or "far-field", or "unperturbed") values.
      call boundary_condition_DG(i, j, ispec, timelocal, rho_DG_P, rhovx_DG_P, rhovz_DG_P, E_DG_P, &
                                 veloc_x_DG_P, veloc_z_DG_P, p_DG_P, e1_DG_P)
      
      ! Set the velocity.
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
      
      ! Deduce energy.
      E_DG_P = p_DG_P/(gammaext_DG(iglobM) - ONE) &
               + rho_DG_P*HALF*( veloc_x_DG_P**2 + veloc_z_DG_P**2 )
      
      ! Deduce momenta.
      rhovx_DG_P   = rho_DG_P*veloc_x_DG_P
      rhovz_DG_P   = rho_DG_P*veloc_z_DG_P
      
      ! The cases below are test cases, intended to test characteristic-based boundary conditions.
      ! TODO: Decide what to do with those.
      if(coord(2, ibool(i, j, ispec)) == 20000. .AND. .false.) then
        ! --------------------------- !
        ! ipoin == -1                 !
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
        ! ipoin == -1           !
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
      
      ! Deduce temperature.
      T_P = (E_DG_P/rho_DG_P - 0.5*(veloc_x_DG_P**2 + veloc_z_DG_P**2))/c_V
      !T_P = (E_DG_iM/rho_DG_iM - 0.5*(veloc_x_DG_iM**2 + veloc_z_DG_iM**2))/c_V ! DEBUG
      !Tx_DG_P = -T_DG_iM(1) ! DEBUG
      !Tz_DG_P = -T_DG_iM(2) ! DEBUG
    endif
  elseif(ispec_is_acoustic_coupling_ac(ibool_DG(i, j, ispec)) >= 0 .AND. .false.) then
    ! --------------------------- !
    ! neighbor(3) /= -1           !
    !   (+) acoustic coupling.    !
    ! --------------------------- !
    ! TODO: Decide what to do with this case.
    !!elseif(.false.) then
    exact_interface_flux = .true.
    ! Coordinates of elastic element
    i_ac     = neighbor(1)
    j_ac     = neighbor(2)
    ispec_ac = neighbor(3)
    iglob = ibool(i_ac,j_ac,ispec_ac)!ibool(i, j, ispec)
    !if(.false.) then ! DEBUG
    !  WRITE(*,*) ibool(i_ac,j_ac,ispec_ac), xixl, xizl, gammaxl, gammazl, duz_dxi, duz_dgamma, k
    !endif
    ! Only for density
    call boundary_condition_DG(i_ac, j_ac, ispec_ac, timelocal, rho_DG_P, &
                               rhovx_DG_P, rhovz_DG_P, E_DG_P, &
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
    T_P = (E_DG_iM/rho_DG_iM - 0.5*(veloc_x_DG_iM**2 + veloc_z_DG_iM**2))/c_V
  else ! Happens everytime (since 'else if' above is an 'if(.false.)').
    ! --------------------------- !
    ! neighbor(3) /= -1           !
    !   (+) not an 'outside' edge !
    !       (thus values are      !
    !       already known)        !
    ! --------------------------- !
    exact_interface_flux = .false.
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
      Vxx_DG_P = V_DG_iP(1, 1)
      Vzz_DG_P = V_DG_iP(2, 2)
      Vxz_DG_P = V_DG_iP(1, 2)
      Vzx_DG_P = V_DG_iP(2, 1)
      Tx_DG_P = T_DG_iP(1)
      Tz_DG_P = T_DG_iP(2)
    endif
    T_P = (E_DG_P/rho_DG_P - 0.5*((rhovx_DG_P/rho_DG_P)**2 + (rhovz_DG_P/rho_DG_P)**2))/c_V
  endif
  
  !close(10)
  
  ! Temperature computation
  !T_P = (E_DG_P/rho_DG_P - 0.5*((rhovx_DG_P/rho_DG_P)**2 + (rhovz_DG_P/rho_DG_P)**2))/c_V
  ! (E_DG_iM/rho_DG_iM - 0.5*(veloc_x_DG_iM**2 + veloc_z_DG_iM**2))/c_V
  
  end subroutine compute_interface_unknowns

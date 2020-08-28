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
  use constants, only: CUSTOM_REAL,NGLLX,NGLLZ,PI
  use specfem_par, only: ibool_DG, rho_DG, rhovx_DG, rhovz_DG, E_DG, e1_DG, nspec
  implicit none
  ! Input/output.
  ! N. A.
  ! Local Variables.
  real(kind=CUSTOM_REAL), parameter :: ZEROl = 0._CUSTOM_REAL
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
      enddo
    enddo
  enddo
end subroutine initial_condition_DG


!********************************!
! TODO: Marked for deletion.     !
!********************************!
! ------------------------------------------------------------ !
! recompute_density                                            !
! ------------------------------------------------------------ !
! TODO: Description.
  
  subroutine recompute_density()

  use constants, only: CUSTOM_REAL, NGLLX, NGLLZ

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
!********************************!
! TODO: Marked for deletion.     !
!********************************!

! ------------------------------------------------------------ !
! boundary_condition_DG                                        !
! ------------------------------------------------------------ !
! Sets all variables of the model at initial time.
! Those are also used as the far-field model.
! This subroutine also takes care of the various forcings.

subroutine boundary_condition_DG(i, j, ispec, timelocal, rho_DG_P, rhovx_DG_P, rhovz_DG_P, E_DG_P, &
                                 veloc_x_DG_P, veloc_z_DG_P, p_DG_P, e1_DG_P)
  use constants, only: CUSTOM_REAL, gamma_euler, PI, HUGEVAL, TINYVAL, FOUR_THIRDS
  use specfem_par, only: ibool_before_perio, ibool_DG, coord, &
        rhoext, windxext, pext_DG, gravityext, gammaext_DG, Bxext, Bzext, N0ext, &
        etaext, muext, coord_interface, kappa_DG, cp, c_V, &
        tau_epsilon, tau_sigma, &
        gravity_cte_DG, dynamic_viscosity_cte_DG, thermal_conductivity_cte_DG, tau_eps_cte_DG, tau_sig_cte_DG, SCALE_HEIGHT, &
        USE_ISOTHERMAL_MODEL, potential_dphi_dx_DG, potential_dphi_dz_DG, ibool, &
        surface_density, sound_velocity, wind, TYPE_FORCING, &
        assign_external_model, IONOSPHERIC_COUPLING

  implicit none
  
  ! Input/output.
  real(kind=CUSTOM_REAL), intent(out) :: rho_DG_P, rhovx_DG_P, rhovz_DG_P, E_DG_P, &
                                         veloc_x_DG_P, veloc_z_DG_P, p_DG_P, e1_DG_P
  
  ! Local Variables.
  integer :: i, j, ispec
  real(kind=CUSTOM_REAL) :: timelocal
  real(kind=CUSTOM_REAL), parameter :: ZEROl = 0._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: ONEl  = 1._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: TWOl  = 2._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: SIXl  = 6._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: HALFl = 0.5_CUSTOM_REAL
  real(kind=CUSTOM_REAL) :: x, z, H, G
  real(kind=CUSTOM_REAL) :: Tl, Tu, rho0, p0, RR ! Variables for the thermal bubble illustrative case.
  
  x = real(coord(1, ibool_before_perio(i, j, ispec)), kind=CUSTOM_REAL)
  z = real(coord(2, ibool_before_perio(i, j, ispec)), kind=CUSTOM_REAL)
  z = z - coord_interface
  
  if(assign_external_model) then
    ! If an external model data file is given for initial conditions, use it.
    rho_DG_P     = rhoext(i, j, ispec)
    p_DG_P       = pext_DG(i, j, ispec)
    veloc_x_DG_P = windxext(i, j, ispec) ! Read horizontal wind.
    veloc_z_DG_P = ZEROl ! Impose vertical wind to zero.
    if(abs(timelocal)<TINYVAL) then
      ! Initialise gravity only once.
      potential_dphi_dx_DG(ibool(i, j, ispec)) = ZEROl
      potential_dphi_dz_DG(ibool(i, j, ispec)) = gravityext(i, j, ispec)
    endif
    
  else
    ! If no external model data file is given (no initial conditions were specified), build model.
    ! Set gravity, viscosity coefficients, relaxation times, and gravity potential (Phi) derivatives.
    
    if(abs(timelocal)<TINYVAL) then
      gravityext(i, j, ispec) = ZEROl ! Isobaric case. Since we need to stay hydrostatic, the gravity field needs to stay 0.
      if(USE_ISOTHERMAL_MODEL) then   
        ! Isothermal case.
        gravityext(i, j, ispec) = real(gravity_cte_DG, kind=CUSTOM_REAL)
      endif
      
      gammaext_DG(ibool_DG(i, j, ispec)) = cp/c_V
      muext(i, j, ispec) = dynamic_viscosity_cte_DG
      etaext(i, j, ispec) = FOUR_THIRDS*dynamic_viscosity_cte_DG
      kappa_DG(i, j, ispec) = thermal_conductivity_cte_DG
      tau_epsilon(i, j, ispec) = tau_eps_cte_DG
      tau_sigma(i, j, ispec)   = tau_sig_cte_DG
      
      if(IONOSPHERIC_COUPLING) then
        Bxext(ibool_DG(i, j, ispec)) = ZEROl
        Bzext(ibool_DG(i, j, ispec)) = ZEROl
        N0ext(i, j, ispec) = ZEROl
      endif
      
      potential_dphi_dx_DG(ibool(i, j, ispec)) = ZEROl
      potential_dphi_dz_DG(ibool(i, j, ispec)) = gravityext(i, j, ispec)
    endif
    
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
  
  e1_DG_P = ZEROl ! Auxiliary variable for vibrationnal attenuation.
  
  ! Rayleigh-Taylor instability illustrative case.
  if(.false.) then
    Tl = 2.
    Tu = 1.
    p0   = ONEl
    rho0 = ONEl
    RR   = 1.!8.3145_CUSTOM_REAL
    gammaext_DG(ibool_DG(i, j, ispec)) = 1.4
    gravityext(i, j, ispec) = 1.
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
  
  ! If bottom forcing is activated, set (overwrite) velocities.
  if(TYPE_FORCING>0 .and. abs(z)<TINYVAL) then
    ! Used the routine 'forcing_DG' defined below, instead of crowding this one.
    call forcing_DG(i, j, ispec, timelocal, veloc_z_DG_P)
  endif

  ! Set energy.
  E_DG_P =   p_DG_P/(gammaext_DG(ibool_DG(i, j, ispec)) - ONEl) &
           + rho_DG_P*HALFl*( veloc_x_DG_P**2 + veloc_z_DG_P**2 ) 

  ! Set momenta.
  rhovx_DG_P = rho_DG_P*veloc_x_DG_P
  rhovz_DG_P = rho_DG_P*veloc_z_DG_P
   
end subroutine boundary_condition_DG


! ------------------------------------------------------------ !
! forcing_DG                                                   !
! ------------------------------------------------------------ !
! Set the bottom forcing for some variable (usually vertical velocity).
   
subroutine forcing_DG(i, j, ispec, current_time, forced_SF)
  use constants, only: NDIM, CUSTOM_REAL, PI, TINYVAL
  use specfem_par, only: myrank, t0, coord, coord_interface, &
                         external_forcing, external_forcing_maxtime, extforc_filedt, extforc_maxx, extforc_minx, &
                         extforc_map_ibbp_to_local, forcing_dg_factor, &
                         TYPE_FORCING, forcing_initial_loc, forcing_initial_time, main_spatial_period, main_time_period, &
                         ibool_before_perio
  
  implicit none
  
  ! Input/Output.
  integer, intent(in) :: i, j, ispec
  real(kind=CUSTOM_REAL), intent(in) :: current_time
  real(kind=CUSTOM_REAL), intent(out) :: forced_SF
  
  ! Local Variables.
  real(kind=CUSTOM_REAL), dimension(NDIM) :: x
  integer externalforcingid ! For TYPE_FORCING==10.
  
  if(TYPE_FORCING==0) return ! Safeguard in case routine is accidentally called, to prevent useless computations.
  
  ! Get point position.
  x = real(coord(:, ibool_before_perio(i, j, ispec)), kind=CUSTOM_REAL)
  x(NDIM) = x(NDIM) - coord_interface ! Remove eventual non-zero altitude of interface.
  
  if(abs(x(NDIM))<TINYVAL) then
    ! For now, forcings only occur on bottom boundary. This test thus enables faster treatment.
    
    select case (TYPE_FORCING)
      case (99)
        ! Testing purposes case. This is the place to shove your hard-coded expressions.
        if(abs(x(NDIM))<TINYVAL) then
          ! If z==0.
          !forced_SF = sin(2.*PI*x(1)/10.)
          forced_SF = cos(2.*PI*current_time/0.08)-1.
        else
          ! Leave everything unchanged and leave early.
          return
        endif
      
      case (1)
        ! Time Gaussian derivative.
        forced_SF = -   (2.*(PI/main_time_period)**2) &
                      * (current_time-t0-forcing_initial_time) &
                      * exp(-((current_time-t0-forcing_initial_time)*PI/main_time_period)**2)
        ! Matlab one-liner: TT=???; T0=???; GGG = -(2.*(pi/TT)^2) * (t-T0) .* exp(-((t-T0)*pi/TT).^2)
      
      case (2)
        ! Time and space Gaussian derivative (gravity wave forcing).
        ! Old code. TODO: recenter peak by taking starting delay ("negative times") into account.
        forced_SF = 0.001 &
                    * ( &
                        - (2./(main_time_period/4.)) &
                          *((current_time-(forcing_initial_time-main_time_period/4.))/(main_time_period/4.)) &
                          *(exp(-((current_time-(forcing_initial_time-main_time_period/4.))/(main_time_period/4.))**2)) &
                        + (2./(main_time_period/4.)) &
                          *((current_time-(forcing_initial_time+main_time_period/4.))/(main_time_period/4.)) &
                          *(exp(-((current_time-(forcing_initial_time+main_time_period/4.))/(main_time_period/4.))**2)) &
                      ) &
                    * (   exp(-((x(1)-(forcing_initial_loc-main_spatial_period/4))/(main_spatial_period/4))**2) &
                        - exp(-((x(1)-(forcing_initial_loc+main_spatial_period/4))/(main_spatial_period/4))**2) )
      
      case (3)
        ! Gaussian second derivative.
        forced_SF =   (1.-2.*(PI*(current_time-t0-forcing_initial_time)/main_time_period)**2) &
                    * exp(-((current_time-t0-forcing_initial_time)*PI/main_time_period)**2)
      
      case (4)
        ! Gaussian pulse.
        forced_SF =   exp(-((current_time-t0-forcing_initial_time)*PI/main_time_period)**2) &
                    * exp(-((x(1)-forcing_initial_loc)*pi/main_spatial_period)**2)
      
      case (9)
        ! Hardcoded forcing.
        
        forced_SF = 0.
        
        !if(.false.) then ! Tsunami forcing.
        !  VELOC_TSUNAMI = 200.
        !  main_time_period = 1000.
        !  main_spatial_period = 50000.
        !  if(timelocal < main_time_period) then
        !      forced_SF = (1d0/main_time_period)*exp( -((x-main_spatial_period)/(sqrt(2.)*main_time_period))**2 )
        !  else
        !      forced_SF = &
        !              2.*((VELOC_TSUNAMI)/(sqrt(2.)*main_time_period))*(((x-main_spatial_period)-VELOC_TSUNAMI*(current_time - main_time_period))/(sqrt(2.)*main_time_period))&
        !              *exp( -(((x-main_spatial_period)-VELOC_TSUNAMI*(current_time - main_time_period))/(sqrt(2.)*main_time_period))**2 )
        !  endif
        !endif
        
        !if(.true.) then ! Analytic microbarom forcing.
        !  MICROBAROM_AMPLITUDE = 1. ! Microbarom amplitude. Unit: m.
        !  MICROBAROM_RANGE = 80.*main_spatial_period ! Range around x=0 to which impose microbaroms. Unit: m. Be careful with apodisation below.
        !  MICROBAROM_MAXTIME = 10.5 * main_time_period ! Microbarom active from t=0 to t=MICROBAROM_MAXTIME. Unit: s. Be careful with apodisation below.
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
        !      NORMAL1 = (PI*DT/(0.2*2.*main_time_period)) * sqrt(-2.*log(UNIFORM1))*cos(2.*PI*UNIFORM2) ! Box-Muller method to generate a 1nd N(0, \sigma^2) random variable.
        !      NORMAL2 = (PI*DT/(0.2*2.*main_time_period)) * sqrt(-2.*log(UNIFORM1))*sin(2.*PI*UNIFORM2) ! Box-Muller method to generate a 2nd N(0, \sigma^2) random variable.
        !      XPHASE_RANDOMWALK = XPHASE_RANDOMWALK + NORMAL1
        !      TPHASE_RANDOMWALK = TPHASE_RANDOMWALK + NORMAL2
        !      PHASE_RANDOMWALK_LASTTIME = timelocal
        !    endif
        !    if(abs(x)<MICROBAROM_RANGE) then
        !      forced_SF = MICROBAROM_AMPLITUDE & ! Amplitude.
        !      * sin(2.*PI*x/main_spatial_period+XPHASE_RANDOMWALK) &
        !      !* 0.25*(1.-erf((x-MICROBAROM_RANGE+10.*main_spatial_period)/(5.*main_spatial_period)))*& ! Spatial apodisation.
        !      !       (1.+erf((x+MICROBAROM_RANGE-10.*main_spatial_period)/(5.*main_spatial_period)))& ! Spatial apodisation, continued.
        !      * sin(2.*PI*timelocal/main_time_period+TPHASE_RANDOMWALK) ! &
        !      !* 0.5*(1.-erf((timelocal-MICROBAROM_MAXTIME+main_time_period)/(0.5*main_time_period))) ! Temporal apodisation.
        !      ! Spatial apodisation over 10 periods.
        !      if(abs(x)>MICROBAROM_RANGE-10.*main_spatial_period) then
        !        forced_SF = forced_SF * (1.-((-abs(x)+MICROBAROM_RANGE)/(10.*main_spatial_period)-1.)**2.)
        !      endif
        !      ! Temporal beginning apodisation over 1 period (split into 0.5 and 0.5 for smoothness).
        !      if(timelocal<0.5*main_time_period) then
        !        forced_SF = forced_SF * (2.*(timelocal/main_time_period)**2.)
        !      endif
        !      if(timelocal>=0.5*main_time_period .and. timelocal<=main_time_period) then
        !        forced_SF = forced_SF * ((4.*timelocal)/main_time_period-2.*(timelocal/main_time_period)**2.-1.)
        !      endif
        !      ! Temporal end apodisation over 1.5 periods.
        !      if(timelocal>MICROBAROM_MAXTIME-1.5*main_time_period) then
        !        forced_SF = forced_SF * (1.-((timelocal-MICROBAROM_MAXTIME)/(1.5*main_time_period)+1.)**2.)
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
        ! Forcing read from a file.
        ! The Matlab script 'utils_new/forcings.m' can be used to generate the external bottom forcing file.
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
    ! We're not at the bottom of the mesh, leave everything unchanged and leave early.
    return
  endif ! Endif on abs(z).
  
  ! Scale using factor provided in parfile.
  forced_SF = FORCING_DG_FACTOR * forced_SF
  
end subroutine forcing_DG


! ------------------------------------------------------------ !
! prepare_external_forcing                                     !
! ------------------------------------------------------------ !
! Reads EXTERNAL_FORCING_FILENAME, save values of bottom forcing in an array, and map mesh coordinates to that array.
! This routine is only called by the 'forcing_DG' routine.

subroutine prepare_external_forcing()

  use constants, only: CUSTOM_REAL, HUGEVAL, TINYVAL, NDIM

  use specfem_par, only: EXTERNAL_FORCING_FILENAME, EXTERNAL_FORCING_MAXTIME,&
                         EXTFORC_MINX, EXTFORC_MAXX, &
                         EXTERNAL_FORCING, EXTFORC_MAP_ibbp_TO_LOCAL, DT,&
                         EXTFORC_FILEDT,&
                         nglob_DG,coord,nspec,NGLLX,NGLLZ,myrank,&!NPROC,&
                         only_DG_acoustic,ibool_before_perio,stage_time_scheme,&
                         ADD_PERIODIC_CONDITIONS

  implicit none

  ! Input/output.
  ! N. A.
  
  ! Local Variables.
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
  EXTFORC_MAP_ibbp_TO_LOCAL = HUGE(0) ! Safeguard.
  
  if(myrank==0) then
    write(*,*) "  First scan: counting number of points, figuring file time step, and pairing file points to SPECFEM points."
  endif
  
  tmp_t_old    = 0.
  counting_nx  = .true.
  nx           = 0 ! Counter for number of spatial points found in file.
  nx_paired    = 0 ! Counter for points actually paired with their SPECFEM-DG counterparts. Safeguard.
  EXTFORC_MINX = HUGEVAL ! Start at +inf, will be updated
  EXTFORC_MAXX = -HUGEVAL ! Start at -inf, will be updated
  
  OPEN(100, file=EXTERNAL_FORCING_FILENAME)
  DO
    READ(100,*,iostat=io) t, x, val ! In fact, we do not care about val here.
    if(counting_nx .and. t>tmp_t_old) then
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
              if(     abs(x-Xspcfm(1))<1e-6 &
                 .or. abs(x-Xspcfm(1))/(abs(x)+abs(Xspcfm(1)))<1e-5) then ! 1e-6 is a relative empiric threshold based on precision we can obtain by generating forcings with our Matlab script. Basically, we chose '%.8e' as printing format for x (in the script), which means we have 8 significant digits, thus theoretically we have a relative precision of 1e-8. We chose the threshold to be 1e-6 in order to a have a safety margin. It is acceptable because we expect the GLL points never to be (relatively) 1e-7 apart.
                EXTFORC_MAP_ibbp_TO_LOCAL(ibbp)=nx
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
    EXTERNAL_FORCING(it, ix)=val
    ix=ix+1
    if(mod(ix,nx+1)==0) then
      ix=1
      it=it+1
    endif
    IF (io/=0) EXIT
  ENDDO
  close(100)
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
                                      Vix_P, Viz_P, Ni_P, &
                                      neighbor, &
                                      exact_interface_flux, &
                                      rho_DG_iM, E_DG_iM, rhovx_DG_iM, rhovz_DG_iM, &
                                      V_DG_iM, T_DG_iM, &
                                      Vix_iM, Viz_iM, Ni_iM, &
                                      rho_DG_iP, E_DG_iP, rhovx_DG_iP, rhovz_DG_iP, &
                                      V_DG_iP, T_DG_iP, &
                                      Vix_iP, Viz_iP, Ni_iP, &
                                      nx, nz, weight, timelocal, &
                                      iface1, iface)
  
  use constants, only: CUSTOM_REAL,NGLLX,NGLLZ,gamma_euler

  use specfem_par, only: ibool_DG, &
                         ispec_is_acoustic_forcing, &
                         ACOUSTIC_FORCING, &
                          ispec_is_acoustic_coupling_el, ispec_is_acoustic_coupling_ac, potential_dot_dot_acoustic, veloc_elastic,&
                         DIR_RIGHT, DIR_LEFT, DIR_UP, DIR_DOWN, &
                         buffer_DG_rho_P, buffer_DG_rhovx_P, buffer_DG_rhovz_P, buffer_DG_E_P, buffer_DG_Ni_P, NPROC, &
                         buffer_DG_Vxx_P, buffer_DG_Vzz_P, buffer_DG_Vxz_P, buffer_DG_Vzx_P, buffer_DG_Tz_P, buffer_DG_Tx_P, &
                         p_DG_init, gammaext_DG, muext, etaext, kappa_DG, ibool, c_V, &
                         buffer_DG_gamma_P, coord,  &
                         rho_init, rhovx_init, rhovz_init, E_init, &
                         potential_dot_dot_acoustic, &
                         veloc_vector_acoustic_DG_coupling, MPI_transfer_iface,&
                         ibool_before_perio, ABC_STRETCH, ABC_STRETCH_LEFT,&
                         ABC_STRETCH_RIGHT, ABC_STRETCH_LEFT_LBUF, ABC_STRETCH_RIGHT_LBUF,&
                         mesh_xmin, mesh_xmax, sigma_elastic, vpext, vsext, rhoext, &
                         assign_external_model, kmato, poroelastcoef, density, surface_density, sound_velocity, IONOSPHERIC_COUPLING
                         
  implicit none
  
  ! Input/output.
  integer, intent(in) :: i, j, ispec, iface1, iface
  integer, dimension(3), intent(in) :: neighbor
  real(kind=CUSTOM_REAL), intent(out) :: rho_DG_P, rhovx_DG_P, rhovz_DG_P, &
                                         E_DG_P, veloc_x_DG_P, veloc_z_DG_P, p_DG_P, T_P, &
                                         Tx_DG_P, Tz_DG_P, Vxx_DG_P, Vzz_DG_P, Vzx_DG_P, Vxz_DG_P, &
                                         Vix_P, Viz_P, Ni_P
  real(kind=CUSTOM_REAL), intent(in) :: timelocal      
  logical, intent(out) :: exact_interface_flux
  real(kind=CUSTOM_REAL), intent(in) :: nx, nz
  real(kind=CUSTOM_REAL), intent(in) :: rho_DG_iM, E_DG_iM, rhovx_DG_iM, rhovz_DG_iM, &
                                        rho_DG_iP, E_DG_iP, rhovx_DG_iP, rhovz_DG_iP
  real(kind=CUSTOM_REAL), intent(in) :: Vix_iP, Viz_iP, Vix_iM, Viz_iM, Ni_iM, Ni_iP ! Variables for electron density.
  real(kind=CUSTOM_REAL), dimension(2), intent(in) :: T_DG_iP, T_DG_iM ! Variables for electron density.
  real(kind=CUSTOM_REAL), dimension(2, 2), intent(in) :: V_DG_iP, V_DG_iM ! Variables for electron density.
  
  ! Local Variables.
  real(kind=CUSTOM_REAL), parameter :: ZERO = 0._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: ONE  = 1._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: TWO  = 2._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: HALF = 0.5_CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: FOUR_THIRDS = 4._CUSTOM_REAL/3._CUSTOM_REAL
  real(kind=CUSTOM_REAL) :: tx, tz, normal_v, tangential_v, &
                            veloc_x, veloc_z, weight, &
                            veloc_x_DG_iM, veloc_z_DG_iM, p_DG_iM, gamma_P, e1_DG_P
  integer :: iglobM, i_el, j_el, ispec_el, i_ac, j_ac, ispec_ac, iglob, iglobP, ipoin, num_interface
  real(kind=CUSTOM_REAL), dimension(2, 2) :: trans_boundary
  real(kind=CUSTOM_REAL), dimension(4) :: sumtau, invsumtau, tauplus, tauac
  real(kind=CUSTOM_REAL), dimension(2) :: v_interface
  real(kind=CUSTOM_REAL) :: x ! Variable for coupling deactivation in buffers.
  real(kind=CUSTOM_REAL) :: p_b, rho_b, s_b, rho_inf, v_b_x, v_b_z, c_b, un_b, &
                            rlambda_max, rlambda_min, c_in, c_inf, un_in, un_inf, &
                            deltaZ1, deltaZ2star, p_n, a_n, alpha0 ! Variables for characteristic-based BCs.
  real(kind=CUSTOM_REAL) :: rhol, cpl, csl, lambdal_unrelaxed_elastic, mul_unrelaxed_elastic, cpal, rhoal

  iglobM = ibool_DG(i, j, ispec)
  
  if(.false.) weight=1. ! Horrible hack to get rid of the warning "Warning: Unused dummy argument 'weight'". I'm sorry. TODO: Remove weight from arguments in all calls to compute_interface_unknowns.
  
  ! Deduce velocities from momenta and density (the two latter being constitutive variables, and thus passed to the routine as input parameters).
  veloc_x_DG_iM = rhovx_DG_iM/rho_DG_iM
  veloc_z_DG_iM = rhovz_DG_iM/rho_DG_iM
  ! Deduce pressure from energy (consitutive), density (consitutive), and velocities (deduced above).
  p_DG_iM       =   (gammaext_DG(iglobM) - ONE) &
                  * ( E_DG_iM - (HALF)*rho_DG_iM*( veloc_x_DG_iM**2 + veloc_z_DG_iM**2 ) )
  
  ! Extract derivatives.
  Tx_DG_P  = -T_DG_iM(1)
  Tz_DG_P  = -T_DG_iM(2)
  Vxx_DG_P = -V_DG_iM(1, 1)
  Vzz_DG_P = -V_DG_iM(2, 2)
  Vxz_DG_P = -V_DG_iM(1, 2)
  Vzx_DG_P = -V_DG_iM(2, 1)

  ! Electron density.
  Vix_P = Vix_iM
  Viz_P = Viz_iM
  Ni_P  = Ni_iM
  
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
        if(IONOSPHERIC_COUPLING) Ni_P     = buffer_DG_Ni_P(ipoin, num_interface)

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
         
    elseif(ispec_is_acoustic_coupling_el(i, j, ispec, 3) >= 0 .AND. abs(nx) < 1.) then ! nx < 1 to remove outer edges from coupling
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
      
      ! Elastic parameters
      lambdal_unrelaxed_elastic = poroelastcoef(1, 1, kmato(ispec_el))
      mul_unrelaxed_elastic = poroelastcoef(2, 1, kmato(ispec_el))
      rhol  = density(1, kmato(ispec_el))
      cpl = sqrt((lambdal_unrelaxed_elastic + TWO*mul_unrelaxed_elastic)/rhol)
      csl = sqrt(mul_unrelaxed_elastic/rhol)
      rhoal = surface_density ! This is an approximation, valid unless we have very large topographies (>few km).
      cpal  = sound_velocity ! This is an approximation, valid unless we have very large topographies (>few km).
      if(assign_external_model) then
        rhol = rhoext(i_el,j_el,ispec_el)
        cpl  = vpext(i_el,j_el,ispec_el)
        csl  = vsext(i_el,j_el,ispec_el)
        rhoal = rhoext(i,j,ispec)
        cpal  = vpext(i,j,ispec)
        ! Yes, both fluid and solid external values are stored in the same array, this is not an error.
      endif
      
      ! From Terrana (2017) eq. 56
      ! Terrana, S., Vilotte, J. P., & Guillot, L. (2017). A spectral hybridizable discontinuous Galerkin method for elastic–acoustic wave propagation. Geophysical Journal International, 213(1), 574-602.
      tauac(1) = nx**2
      tauac(2) = nz*nx
      tauac(3) = nz*nx
      tauac(4) = nz**2
      tauac = rhoal*cpal*tauac
      tauplus(1) = cpl*nx**2 + csl*nz**2
      tauplus(2) = (cpl - csl)*nx*nz
      tauplus(3) = (cpl - csl)*nx*nz
      tauplus(4) = csl*nx**2 + cpl*nz**2
      
      tauplus = rhol*tauplus
      sumtau    = tauplus + tauac
      invsumtau(1) = sumtau(4) 
      invsumtau(2) = -sumtau(2) 
      invsumtau(3) = -sumtau(3) 
      invsumtau(4) = sumtau(1) 
      invsumtau    = invsumtau/(sumtau(1)*sumtau(4) - sumtau(2)*sumtau(3))
      
      v_interface(1) = tauac(1)*veloc_x_DG_iM + tauac(2)*veloc_z_DG_iM &
        + tauplus(1)*veloc_elastic(1, iglob) + tauplus(2)*veloc_elastic(2, iglob) &
        - (p_DG_iM- p_DG_init(ibool_DG(i,j,ispec)))*nx &
        + sigma_elastic(1,1,iglob)*nx + sigma_elastic(1,2,iglob)*nz
        
      v_interface(2) = tauac(3)*veloc_x_DG_iM + tauac(4)*veloc_z_DG_iM &
        + tauplus(3)*veloc_elastic(1, iglob) + tauplus(4)*veloc_elastic(2, iglob) &
        - (p_DG_iM- p_DG_init(ibool_DG(i,j,ispec)))*nz &
        + sigma_elastic(2,1,iglob)*nx + sigma_elastic(2,2,iglob)*nz
        
      v_interface(1) = v_interface(1)*invsumtau(1) + v_interface(2)*invsumtau(2)
      v_interface(2) = v_interface(1)*invsumtau(3) + v_interface(2)*invsumtau(4)
      
      ! Get elastic medium velocities.
      veloc_x = 0*veloc_elastic(1, iglob) + 1*v_interface(1)
      veloc_z = 0*veloc_elastic(2, iglob) + 1*v_interface(2)
      ! Tangential vector. It is assumed that we only have the elastic media under the DG medium, hence we always have nz>=0.
      tx = -nz
      tz =  nx
      ! Get the normal velocity of the solid perturbation (normal continuity) and the tangential velocity of the fluid flow (free slip).
      normal_v     = (veloc_x*nx + veloc_z*nz)
      tangential_v = veloc_x_DG_iM*tx + veloc_z_DG_iM*tz
      ! Set the matrix for the transformation from normal/tangential coordinates to mesh coordinates.
      trans_boundary(1, 1) =  tz
      trans_boundary(1, 2) = -nz
      trans_boundary(2, 1) = -tx
      trans_boundary(2, 2) =  nx
      trans_boundary = trans_boundary/(nx * tz - tx * nz)
      
      ! Deactivate air/ground coupling in the buffer regions.
      ! This method gradually deactivates the coupling in the horizontal fluid ABC buffers.
      ! This is somewhat of a hack, but which works to an appreciable level.
      ! A mathematically sounder method would be preferred. For instance, write out the system of equations at the boundary, and apply the real stretching method (see Martire GJI 2020) to the whole set.
      x = coord(1, ibool_before_perio(i, j, ispec))
      if(ABC_STRETCH .and. &
         (     (ABC_STRETCH_LEFT   .and. x < mesh_xmin + ABC_STRETCH_LEFT_LBUF) & ! left stretching and in left buffer zone
          .or. (ABC_STRETCH_RIGHT  .and. x > mesh_xmax - ABC_STRETCH_RIGHT_LBUF) & ! right stretching and in right buffer zone
         ) &
        ) then
        
        if(ABC_STRETCH_LEFT) x = (x - mesh_xmin)/ABC_STRETCH_LEFT_LBUF ! x is a local buffer coordinate now (1 at beginning, 0 at end).
        if(ABC_STRETCH_RIGHT) x = (mesh_xmax - x)/ABC_STRETCH_RIGHT_LBUF ! x is a local buffer coordinate now (1 at beginning, 0 at end).
        
        ! Convert (back) the velocity components from normal/tangential coordinates to mesh coordinates.
        veloc_x_DG_P = x*(trans_boundary(1, 1)*normal_v + trans_boundary(1, 2)*tangential_v)&
                       +(1.-x)*veloc_x_DG_P
        veloc_z_DG_P = x*(trans_boundary(2, 1)*normal_v + trans_boundary(2, 2)*tangential_v)&
                       +(1.-x)*veloc_z_DG_P
      else
        ! Convert (back) the velocity components from normal/tangential coordinates to mesh coordinates.
        veloc_x_DG_P = trans_boundary(1, 1)*normal_v + trans_boundary(1, 2)*tangential_v
        veloc_z_DG_P = trans_boundary(2, 1)*normal_v + trans_boundary(2, 2)*tangential_v
      endif
      
      ! Stress continuity.
      p_DG_P = p_DG_iM + 1*( (tauac(1)*(veloc_x_DG_iM - veloc_x) + tauac(2)*(veloc_z_DG_iM - veloc_z)) * nx  &
        + (tauac(3)*(veloc_x_DG_iM - veloc_x) + tauac(4)*(veloc_z_DG_iM - veloc_z)) * nz )

      ! Deduce energy.
      E_DG_P = p_DG_P/(gammaext_DG(iglobM) - 1.) + HALF*rho_DG_P*( veloc_x_DG_P**2 + veloc_z_DG_P**2 )

      ! Recompute momenta in order not to keep those coming from the call to boundary_condition_DG above.
      rhovx_DG_P = rho_DG_P*veloc_x_DG_P
      rhovz_DG_P = rho_DG_P*veloc_z_DG_P
      
      ! Deduce temperature.
      T_P = (E_DG_iM/rho_DG_iM - 0.5*(veloc_x_DG_iM**2 + veloc_z_DG_iM**2))/c_V

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
    endif
    
  elseif(ispec_is_acoustic_coupling_ac(ibool_DG(i, j, ispec)) >= 0 .AND. .false.) then
    ! --------------------------- !
    ! neighbor(3) /= -1           !
    !   (+) acoustic coupling.    !
    ! --------------------------- !
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
    p_DG_P       =   (gamma_P - ONE) &
                   * ( E_DG_P - (HALF)*rho_DG_P*( veloc_x_DG_P**2 + veloc_z_DG_P**2 ) )
    
    ! Electron density.
    if(IONOSPHERIC_COUPLING) then
      Vix_P = Vix_iP
      Viz_P = Viz_iP
      Ni_P = Ni_iP
    endif
    
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
  
end subroutine compute_interface_unknowns


! ------------------------------------------------------------ !
! compute_forces_acoustic_LNS_main                              !
! ------------------------------------------------------------ !
! TODO: Description.

subroutine compute_forces_acoustic_LNS_main()
  
  use constants ! TODO: select variables to use.
  use specfem_par ! TODO: select variables to use.
  use specfem_par_LNS ! TODO: select variables to use.
  !use constants, only: &
  !  CUSTOM_REAL, &
  !  rk4a_d, rk4b_d, rk4c_d, &
  !  ls33rk_a, ls33rk_b, ls33rk_c, &
  !  HALF, ONE, NGLLX, NGLLZ
  !use specfem_par, only:&! kmato,&
  !  deltat, nglob_DG,stage_time_scheme,&
  !  !gravityext, muext, etaext, kappa_DG,& ! TODO: make sure those aren't allocated
  !  !tau_epsilon, tau_sigma, & ! TODO: make sure those aren't allocated
  !  ispec_is_acoustic_coupling_ac,&
  !  rmass_inverse_acoustic_DG,&
  !  assign_external_model, any_acoustic_DG, only_DG_acoustic, &
  !  !c_V,&!,gammaext_DG,&
  !  it, i_stage,&
  !  myrank, nspec, nproc, ibool_DG, ibool_before_perio, &
  !  ninterface_acoustic, &
  !  time_stepping_scheme, &
  !  !CONSTRAIN_HYDROSTATIC, &
  !  coord!, &
  !  !T_DG, V_DG ! Use those to store respectively \nabla T' and \nabla v'.
  !use specfem_par_LNS, only: &
  !  LNS_g, LNS_mu, LNS_eta, LNS_kappa,&
  !  LNS_rho0, LNS_v0, LNS_E0,&
  !  LNS_drho, LNS_rho0dv, LNS_dE,&
  !  aux_drho, aux_rho0dv, aux_dE,&
  !  RHS_drho, RHS_rho0dv, RHS_dE,&
  !  nabla_v0, LNS_T0, sigma_v_0, LNS_p0,&
  !  LNS_dm, LNS_dp, LNS_dT,&
  !  buffer_LNS_drho_P, buffer_LNS_rho0dv_P, buffer_LNS_dE_P,&
  !  LNS_VERBOSE,&
  !  SPACEDIM,&
  !  LNS_dummy_1d, LNS_dummy_2d!, LNS_dummy_3d

  implicit none

  ! Local.
  real(kind=CUSTOM_REAL), parameter :: ZEROcr = 0._CUSTOM_REAL
  real(kind=CUSTOM_REAL) :: timelocal
  real(kind=CUSTOM_REAL), dimension(stage_time_scheme) :: scheme_A, scheme_B, scheme_C
  real(kind=CUSTOM_REAL), dimension(SPACEDIM,nglob_DG) :: LNS_dv
  real(kind=CUSTOM_REAL), dimension(SPACEDIM,SPACEDIM,nglob_DG) :: nabla_dv
  real(kind=CUSTOM_REAL), dimension(3,nglob_DG) :: sigma_dv
  real(kind=CUSTOM_REAL), dimension(SPACEDIM,nglob_DG) :: nabla_dT
  integer :: i,j,ispec,ier
  logical CHECK_NONPOSITIVITY, CHECK_NONPOSITIVITY_ON_ALL_PROCS, CHECK_NONPOSITIVITY_FIND_POINT
  logical switch_gradient
  
  ! Checks if anything has to be done.
  if (.not. any_acoustic_DG) then
    return
  endif
  
  switch_gradient = .false. ! Switch to activate the use of the "desintegration method" for gradient computation methods, and to desactivate to use the SEM definition of the gradient.
  
  ! Those are debug switches. They are computationaly very heavy and should not be used on every simulation.
  ! CHECK_NONPOSITIVITY_ON_ALL_PROCS=.true. and CHECK_NONPOSITIVITY_FIND_POINT=.true. are particularly heavy.
  CHECK_NONPOSITIVITY=.false.              ! Set to .true. to enable nonpositivity checking.
  CHECK_NONPOSITIVITY_ON_ALL_PROCS=.false. ! Only used if CHECK_NONPOSITIVITY==.true.. Set to .false. for checking only on proc 0. Set to .true. for checking on all procs.
  CHECK_NONPOSITIVITY_FIND_POINT=.false.   ! Only used if CHECK_NONPOSITIVITY==.true.. Set to .true. to find where nonpositivity was encountered.
  
  ! Stop if time stepping scheme is not implemented.
  if(time_stepping_scheme/=3 .and. time_stepping_scheme/=4 .and. myrank==0) then
    write(*,*) "********************************"
    write(*,*) "*            ERROR             *"
    write(*,*) "********************************"
    write(*,*) "* This time_stepping_scheme is *"
    write(*,*) "* not implemented for DG yet.  *"
    write(*,*) "* time_stepping_scheme ", time_stepping_scheme
    write(*,*) "********************************"
    stop
  endif
  ! Load LSRK coefficients.
  if (time_stepping_scheme == 3) then
    ! 5 stages.
    scheme_A = real(rk4a_d, kind=CUSTOM_REAL)
    scheme_B = real(rk4b_d, kind=CUSTOM_REAL)
    scheme_C = real(rk4c_d, kind=CUSTOM_REAL)
  else if (time_stepping_scheme == 4) then
    ! 3 stages.
    scheme_A = real(ls33rk_a, kind=CUSTOM_REAL)
    scheme_B = real(ls33rk_b, kind=CUSTOM_REAL)
    scheme_C = real(ls33rk_c, kind=CUSTOM_REAL)
  endif
  ! Compute current time.
  timelocal = (it-1)*deltat + scheme_C(i_stage)*deltat
  
  ! Intialisation.
  if(it == 1 .and. i_stage == 1) then
    !if(USE_SLOPE_LIMITER) then
    !  ! The Vandermonde matrices are only used when the slope limiter is.
    !  call setUpVandermonde()
    !endif
    
    ! Initialise auxiliary registers.
    aux_drho    = ZEROcr
    aux_rho0dv = ZEROcr
    aux_dE      = ZEROcr
    RHS_drho    = ZEROcr
    RHS_rho0dv = ZEROcr
    RHS_dE      = ZEROcr
    ! Initialise state registers. Note: since constitutive variables are perturbations, they are necessarily zero at start.
    LNS_drho   = ZEROcr
    LNS_rho0dv = ZEROcr
    LNS_dE     = ZEROcr
    
    ! Prepare MPI buffers.
#ifdef USE_MPI
    if(NPROC>1) then
      call prepare_MPI_DG()
      ! TODO: buffers will need to be different than DG-FNS's ones, thus one would have to call a dedicated function.
    endif
#endif
    
    !write(*,*) kmato ! DEBUG
    
    ! Allocate acoustic coupling array.
    allocate(ispec_is_acoustic_coupling_ac(nglob_DG), stat=ier)
    if (ier /= 0) then
      stop "Error allocating 'ispec_is_acoustic_coupling_ac' arrays (see 'compute_forces_acoustic_LNS_calling_routine.F90')."
    endif
    ispec_is_acoustic_coupling_ac = -1
    
    if(.not. only_DG_acoustic) then
      call find_DG_acoustic_coupling()
    endif
    
    !write(*,*) kmato ! DEBUG
    
    ! Physical parameters.
    if(.not. assign_external_model) then
      !deallocate(gravityext, muext, etaext, kappa_DG)
      !allocate(gravityext(NGLLX, NGLLZ, nspec), &
      !         etaext(NGLLX, NGLLZ, nspec), &
      !         muext(NGLLX, NGLLZ, nspec), &
      !         kappa_DG(NGLLX, NGLLZ, nspec))
      !deallocate(tau_epsilon, tau_sigma) ! Temporary hack until boundary_condition_DG is modified.
      !allocate(tau_epsilon(NGLLX, NGLLZ, nspec), tau_sigma(NGLLX, NGLLZ, nspec)) ! Temporary hack until boundary_condition_DG is modified.
    endif
    ! Send the (:,:,:) classic registers depending on (i,j,ispec) into one (:) register depending on (iglob), in order to be able to use vector calculus.
    ! This is only done here because we use the subroutine "boundary_condition_DG" in "initial_state_LNS", which needs the classic registers.
    allocate(LNS_g(nglob_DG), &
             LNS_mu(nglob_DG), &
             LNS_eta(nglob_DG), &
             LNS_kappa(nglob_DG))
    LNS_g=9.81
    LNS_mu=0.
    LNS_eta=0.!(4./3.)*LNS_mu
    LNS_kappa=0.
    !deallocate(gravityext, muext, etaext, kappa_DG) ! Not really optimal, but less invasive.
    
    ! Initialise initial state.
    !write(*,*) "call initial_state_LNS" ! DEBUG
    call initial_state_LNS()
    ! Initialise T_0 and p_0.
    !call compute_p(LNS_rho0, LNS_v0, LNS_E0, LNS_p0) ! In fact, LNS_p0 was computed in initial_state_LNS, but we recompute p_0 based on energy to be coherent with rest of program.
    call compute_T(LNS_rho0, LNS_v0, LNS_E0, LNS_T0)
    !LNS_T0=LNS_p0/(LNS_rho0*287.05684504212125397) ! Approximation assuming air molar mass is constant at 0.028964.
    
    where(LNS_p0 <= 0._CUSTOM_REAL) LNS_p0 = 1._CUSTOM_REAL ! When doing elastic-DG simulations, LNS_p0 = 0 in elastic elements and 1/LNS_p0 will not be properly defined (but should not happen). Use a hack if problems occur.
    ! Initialise \nabla v_0.
    call compute_gradient_TFSF(LNS_v0, LNS_dummy_1d, .true., .false., switch_gradient, nabla_v0, LNS_dummy_2d) ! Dummy variables are not optimal, but prevent from duplicating subroutines.
    ! Initialise \Sigma_v_0.
    call LNS_compute_viscous_stress_tensor(nabla_v0, sigma_v_0)
    
#ifdef USE_MPI
    if (NPROC > 1 .and. ninterface_acoustic > 0) then
      !call assemble_MPI_vector_DG(gammaext_DG, buffer_LNS_gamma_P)
      ! TODO: call a dedicated routine.
    endif
#endif
  endif ! Endif on (it == 1) and (i_stage == 1).
  
  if(LNS_VERBOSE>=1 .and. myrank == 0 .AND. mod(it, LNS_MODPRINT)==0) then
    WRITE(*,*) "****************************************************************"
    WRITE(*,"(a,i9,a,i1,a,e23.16,a)") "Iteration", it, ", stage ", i_stage, ", local time", timelocal, "."
  endif
  
#ifdef USE_MPI
  if (NPROC > 1 .and. ninterface_acoustic > 0) then
    call assemble_MPI_vector_DG(LNS_drho, buffer_LNS_drho_P)
    do i=1,SPACEDIM
      call assemble_MPI_vector_DG(LNS_rho0dv(i,:), buffer_LNS_rho0dv_P(i,:,:))
    enddo
    call assemble_MPI_vector_DG(LNS_dE, buffer_LNS_dE_P)
  endif
#endif

  ! Local Discontinuous Galerkin for viscous fluxes.
  !if((maxval(LNS_mu) > 0 .OR. maxval(LNS_eta) > 0 .OR. maxval(LNS_kappa) > 0) .OR. CONSTRAIN_HYDROSTATIC) then
  if(     maxval(LNS_mu) > 0. &
     .OR. maxval(LNS_eta) > 0. &
     .OR. maxval(LNS_kappa) > 0.) then
#ifdef USE_MPI
    if (NPROC > 1 .and. ninterface_acoustic > 0) then
      !call assemble_MPI_vector_DG(T_DG(1, :), buffer_LNS_Tx_P)
      !call assemble_MPI_vector_DG(T_DG(2, :), buffer_LNS_Tz_P)
      !call assemble_MPI_vector_DG(V_DG(1, 1, :), buffer_LNS_Vxx_P)
      !call assemble_MPI_vector_DG(V_DG(2, 2, :), buffer_LNS_Vzz_P)
      !call assemble_MPI_vector_DG(V_DG(1, 2, :), buffer_LNS_Vxz_P)
      !call assemble_MPI_vector_DG(V_DG(2, 1, :), buffer_LNS_Vzx_P)
    endif
#endif
    !call compute_viscous_tensors(T_DG, V_DG, LNS_drho, LNS_rho0dv(1,:), LNS_rho0dv(SPACEDIM,:), LNS_dE, timelocal)
  endif
  
  ! Precompute momentum perturbation and velocity perturbation.
  LNS_dv=ZEROcr
  do i=1,SPACEDIM
    LNS_dm(i,:) = LNS_rho0dv(i,:)+LNS_drho*LNS_v0(i,:)
    where(LNS_rho0/=ZEROcr) LNS_dv(i,:)=LNS_rho0dv(i,:)/LNS_rho0 ! 'where(...)' as safeguard, as rho0=0 should not happen.
  enddo
  
  ! Recompute temperature and pressure perturbation.
  call compute_dT(LNS_rho0+LNS_drho, LNS_v0+LNS_dv, LNS_E0+LNS_dE, LNS_dT)
  call compute_dp(LNS_rho0+LNS_drho, LNS_v0+LNS_dv, LNS_E0+LNS_dE, LNS_dp)
  
  !write(*,*)it, i_stage, "rho0", LNS_rho0(1), "E0", LNS_E0(1), "T0", LNS_T0(1), "p0", LNS_p0(1), &
  !             "gamma", gammaext_DG(1), "c_V", c_V, &
  !             "dp", LNS_dp(1)! DEBUG
  !stop
  
  ! Precompute gradients.
  call compute_gradient_TFSF(LNS_dv, LNS_dT, .true., .true., switch_gradient, nabla_dv, nabla_dT) ! Dummy variables are not optimal, but prevent from duplicating
  ! Precompute \Sigma_v'.
  call LNS_compute_viscous_stress_tensor(nabla_dv, sigma_dv)
  
  ! Compute RHS.
  call compute_forces_acoustic_LNS(LNS_drho, LNS_rho0dv, LNS_dE, & ! Constitutive variables.
                                   LNS_dm, LNS_dp, LNS_dT, nabla_dT, sigma_dv, & ! Precomputed quantities.
                                   RHS_drho, RHS_rho0dv, RHS_dE, & ! Output.
                                   timelocal) ! Time.
  
  if (time_stepping_scheme == 3 .or. time_stepping_scheme == 4) then
    ! Inverse mass matrix multiplication, in order to obtain actual RHS.
    RHS_drho(:)            = RHS_drho(:)*rmass_inverse_acoustic_DG(:) ! RHS = A^{-1}*b
    RHS_dE(:)              = RHS_dE(:)  *rmass_inverse_acoustic_DG(:)
    ! Update the auxiliary register.
    ! Note: no need to zero it beforehand if scheme_A(1) is 0.
    aux_drho               = scheme_A(i_stage)*aux_drho + deltat*RHS_drho ! U_{i} = a_{i}*U_{i-1} + dt*RHS
    aux_dE                 = scheme_A(i_stage)*aux_dE   + deltat*RHS_dE
    ! Update the state register.
    LNS_drho               = LNS_drho + scheme_B(i_stage)*aux_drho ! Y^{n+1} = Y^{n+1} + b_i*U_{i}
    LNS_dE                 = LNS_dE   + scheme_B(i_stage)*aux_dE
    ! Group momentum treatment in loop.
    do i=1,SPACEDIM
      RHS_rho0dv(i,:) = RHS_rho0dv(i,:)*rmass_inverse_acoustic_DG(:)
      aux_rho0dv(i,:) = scheme_A(i_stage)*aux_rho0dv(i,:) + deltat*RHS_rho0dv(i,:)
      LNS_rho0dv(i,:) = LNS_rho0dv(i,:) + scheme_B(i_stage)*aux_rho0dv(i,:)
    enddo
    
    ! Eventually check non-positivity.
    if(CHECK_NONPOSITIVITY) then
      if(     CHECK_NONPOSITIVITY_ON_ALL_PROCS &
         .or. ((.not. CHECK_NONPOSITIVITY_ON_ALL_PROCS) .and. myrank==0) &
        ) then
        ! If:    we check on all procs,
        !     or we check only on proc 0 and we are on proc 0.
        if(minval(LNS_drho) < 10d-14) then
          WRITE(*,*) "***************************************************************"
          WRITE(*,*) "* CAREFUL, VERY SMALL DENSITY: ", minval(LNS_drho), "    *"
          if(CHECK_NONPOSITIVITY_FIND_POINT) then
            ! Find where density is low.
            do ispec = 1,nspec
              do j = 1,NGLLZ
                do i = 1,NGLLX
                  !write(*,*) LNS_drho(ibool_DG(ispec, i, j))
                  !write(*,*) ispec, i, j, ibool_DG(ispec, i, j)
                  if(LNS_drho(ibool_DG(i, j, ispec))==minval(LNS_drho)) then
                    WRITE(*, *) "* Element", ispec, ", GLL", i, j, ".         *"
                    write(*, *) "* Coords", coord(1, ibool_before_perio(i, j, ispec)), &
                                coord(SPACEDIM, ibool_before_perio(i, j, ispec)), ".*"
                    !call virtual_stretch_prime(i, j, ispec, coef_stretch_x_ij_prime, coef_stretch_z_ij_prime)
                    !write(*, *) coef_stretch_x_ij_prime, coef_stretch_z_ij_prime
                  endif
                enddo
              enddo
            enddo
          endif ! Endif on CHECK_NONPOSITIVITY_FIND_POINT.
          WRITE(*,*) "***************************************************************"
        endif ! Endif on minval(LNS_drho).
      endif ! Endif on CHECK_NONPOSITIVITY_ON_ALL_PROCS.
    endif ! Endif on CHECK_NONPOSITIVITY.
  endif
  
  ! --------------------------- !
  ! Remove high-order           !
  ! coefficients of the         !
  ! solution.                   !
  ! --------------------------- !
  !if(USE_SLOPE_LIMITER) then
  !  ! rho.
  !  if(CONSTRAIN_HYDROSTATIC) then
  !    veloc_x = LNS_drho - LNS_rho0
  !    call SlopeLimit1(veloc_x, timelocal, 1)
  !    LNS_drho = veloc_x + LNS_rho0
  !  else
  !    call SlopeLimit1(LNS_drho, timelocal, 1)
  !  endif
  !  ! rhovx.
  !  if(CONSTRAIN_HYDROSTATIC) then
  !    veloc_x = LNS_rho0dv(1,:) - LNS_v0(1,:)
  !    call SlopeLimit1(veloc_x, timelocal, 1)
  !    LNS_rho0dv(1,:) = veloc_x + LNS_v0(1,:)
  !  else
  !    call SlopeLimit1(LNS_rho0dv(1,:), timelocal, 2)
  !  endif
  !  ! rhovz.
  !  if(CONSTRAIN_HYDROSTATIC) then
  !    veloc_x = LNS_rho0dv(SPACEDIM,:) - LNS_v0(SPACEDIM,:)
  !    call SlopeLimit1(veloc_x, timelocal, 1)
  !    LNS_rho0dv(SPACEDIM,:) = veloc_x + LNS_v0(SPACEDIM,:)
  !  else
  !    call SlopeLimit1(LNS_rho0dv(SPACEDIM,:), timelocal, 3)
  !  endif
  !  ! E.
  !  if(CONSTRAIN_HYDROSTATIC) then
  !    veloc_x = LNS_dE - LNS_E0
  !    call SlopeLimit1(veloc_x, timelocal, 1)
  !    LNS_dE = veloc_x + LNS_E0
  !  else
  !    call SlopeLimit1(LNS_dE, timelocal, 4)
  !  endif
  !endif
end subroutine compute_forces_acoustic_LNS_main





















! ------------------------------------------------------------ !
! initial_state_LNS                                            !
! ------------------------------------------------------------ !
! Computes initial state.

subroutine initial_state_LNS()
  use constants ! TODO: select variables to use.
  use specfem_par ! TODO: select variables to use.
  use specfem_par_LNS ! TODO: select variables to use.

  implicit none
  
  ! Input/Output.
  ! N/A.
  
  ! Local.
  real(kind=CUSTOM_REAL), parameter :: ZEROcr = 0._CUSTOM_REAL
  integer :: ispec, iglob, i, j
  !real(kind=CUSTOM_REAL) :: rho_P, veloc_x_P, veloc_z_P, E_P
  
  do ispec = 1, nspec
    do j = 1, NGLLZ
      do i = 1, NGLLX
        iglob = ibool_DG(i, j, ispec)
        !write(*,*) "call background_physical_parameters" ! DEBUG
        call background_physical_parameters(i, j, ispec, ZEROcr, LNS_rho0(iglob), LNS_v0(:,iglob), LNS_E0(iglob), LNS_p0(iglob))
      enddo
    enddo
  enddo
end subroutine initial_state_LNS

! ------------------------------------------------------------ !
! initial_state_LNS                                            !
! ------------------------------------------------------------ !
! Affects values of background state. May thus be used as initialiser (if time is 0), for far-field boundary conditions, or for bottom forcings.

subroutine background_physical_parameters(i, j, ispec, timelocal, out_rho, out_v, out_E, out_p)
  use constants, only: CUSTOM_REAL
  use specfem_par ! TODO: select variables to use. , only: ibool_before_perio, ibool_DG, coord, &
        !rhoext, windxext, pext_DG, gravityext, gammaext_DG, &
        !etaext, muext, coord_interface, kappa_DG, cp, c_V, &
        !tau_epsilon, tau_sigma, &
        !gravity_cte_DG, dynamic_viscosity_cte_DG, thermal_conductivity_cte_DG, tau_eps_cte_DG, tau_sig_cte_DG, SCALE_HEIGHT, &
        !USE_ISOTHERMAL_MODEL, potential_dphi_dx_DG, potential_dphi_dz_DG, ibool, &
        !surface_density, sound_velocity, wind, TYPE_FORCING, &
        !forcing_initial_time, main_time_period, forcing_initial_loc, main_spatial_period,&
        !assign_external_model,myrank, &
        !DT, XPHASE_RANDOMWALK, TPHASE_RANDOMWALK, PHASE_RANDOMWALK_LASTTIME,& ! Microbarom forcing.
        !EXTERNAL_FORCING_MAXTIME,EXTERNAL_FORCING, EXTFORC_MAP_ibbp_TO_LOCAL,& ! External forcing.
        !EXTFORC_MINX, EXTFORC_MAXX,EXTFORC_FILEDT! External forcing.
  use specfem_par_LNS ! TODO: select variables to use.

  implicit none
  
  ! Input/Output.
  integer, intent(in) :: i, j, ispec
  real(kind=CUSTOM_REAL), intent(out) :: out_rho, out_E, out_p
  real(kind=CUSTOM_REAL), dimension(2), intent(out) :: out_v
  real(kind=CUSTOM_REAL), intent(in) :: timelocal
  
  ! Local.
  integer :: iglob
  real(kind=CUSTOM_REAL), parameter :: ZEROcr = 0._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: HALFcr = 0.5_CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: ONEcr = 1._CUSTOM_REAL
  real(kind=CUSTOM_REAL) :: z, H!, G!, A
  
  iglob = ibool_DG(i, j, ispec)
  if(assign_external_model) then
    ! If an external model data file is given for initial conditions, read from it.
    out_rho = rhoext(i, j, ispec)
    out_p = pext_DG(i, j, ispec)
    out_v(1) = windxext(i, j, ispec)
  else
    ! If no external model data file is given (no initial conditions were specified), build model.
    
    ! > If initialisation (condition on timelocal), set gravity, viscosity coefficients (mu and eta), and gamma.
    !write(*,*) "kek" ! DEBUG
    if(abs(timelocal)<TINYVAL) then
      ! We are at t=0.
      if(USE_ISOTHERMAL_MODEL) then
        ! > Isothermal case.
        LNS_g(iglob) = real(gravity_cte_DG, kind=CUSTOM_REAL)
      else
        ! > Isobaric case. Since we need to stay hydrostatic, the gravity field needs to stay 0.
        LNS_g(iglob) = ZEROcr
      endif
      gammaext_DG(iglob) = cp/c_V
      LNS_mu(iglob) = dynamic_viscosity_cte_DG
      LNS_eta(iglob) = (4./3.)*dynamic_viscosity_cte_DG
      LNS_kappa(iglob) = thermal_conductivity_cte_DG
    endif
    
    if(USE_ISOTHERMAL_MODEL) then
      ! > Set density.
      H = SCALE_HEIGHT ! Also for pressure, below.
      out_rho = surface_density*exp(-z/H)
      ! > Set pressure.
      out_p = out_rho*LNS_g(iglob)*H
    else
      ! > Set density.
      out_rho = surface_density
      ! > Set pressure.
      !write(*,*) sound_velocity ! DEBUG
      !write(*,*) out_rho ! DEBUG
      !write(*,*) gammaext_DG(ibool_DG(i, j, ispec)) ! DEBUG
      out_p = (sound_velocity**2)*out_rho/gammaext_DG(iglob) ! Acoustic only (under ideal gas hypothesis): p = c^2 * \rho / \gamma.
    endif
    
    ! > Set wind.
    out_v(1) = wind ! Read horizontal wind from the scalar value read from parfile.
  endif ! Endif on assign_external_model.
  
  ! Impose vertical wind to zero.
  out_v(SPACEDIM) = ZEROcr 
  
  !! > Set gravity potentials.
  !if(timelocal == ZEROcr) then
  !  potential_dphi_dx_DG(ibool(i, j, ispec)) = ZEROcr
  !  potential_dphi_dz_DG(ibool(i, j, ispec)) = gravityext(i, j, ispec)
  !endif
  
  ! Set energy based on pressure.
  out_E =   out_p/(gammaext_DG(iglob) - ONEcr) &
                  + out_rho*HALFcr*( out_v(1)**2 + out_v(SPACEDIM)**2 )
  
end subroutine background_physical_parameters





















! ------------------------------------------------------------ !
! LNS_compute_viscous_stress_tensor                            !
! ------------------------------------------------------------ !
! Computes the viscous Navier-Stokes stress tensor \Sigma_v for a given velocity v.
! IN:
!   \nabla v, gradient of v, (1,1) being dxvx, (1,2) dzvx, (2,1) dxvz, and (2,2) dzvz.
! OUT:
!   \Sigma_v(v), as a size 3 vector since the stress tensor is symmetric: 1<->(1,1), 2<->(1,2)&(2,1), and 3<->(2,2).

subroutine LNS_compute_viscous_stress_tensor(nabla_v,&
                                             sigma_v)
  use constants ! TODO: select variables to use.
  use specfem_par ! TODO: select variables to use.
  use specfem_par_LNS ! TODO: select variables to use.
  
  implicit none
  
  ! Input/Output.
  real(kind=CUSTOM_REAL), dimension(SPACEDIM, SPACEDIM, nglob_DG), intent(in) :: nabla_v
  real(kind=CUSTOM_REAL), dimension(3, nglob_DG), intent(out) :: sigma_v
  
  ! Local.
  real(kind=CUSTOM_REAL), parameter :: TWOcr = 2._CUSTOM_REAL
  real(kind=CUSTOM_REAL), dimension(nglob_DG) :: LNS_lambda
  
  LNS_lambda = LNS_eta-(2._CUSTOM_REAL/3._CUSTOM_REAL)*LNS_mu
  
  ! Explicit.
  !sigma_v(1,:) = TWOcr*LNS_mu*nabla_v(1,1,:) + LNS_lambda*(nabla_v(1,1,:)+nabla_v(2,2,:))
  !sigma_v(2,:) = LNS_mu*(nabla_v(1,2,:)+nabla_v(2,1,:))
  !sigma_v(3,:) = TWOcr*LNS_mu*nabla_v(2,2,:) + LNS_lambda*(nabla_v(1,1,:)+nabla_v(2,2,:))
  
  ! Compact.
  sigma_v(1,:) = (TWOcr*LNS_mu+LNS_lambda)*nabla_v(1,1,:) + LNS_lambda*nabla_v(SPACEDIM,2,:)
  sigma_v(2,:) = LNS_mu*(nabla_v(1,2,:)+nabla_v(SPACEDIM,1,:))
  sigma_v(3,:) = (TWOcr*LNS_mu+LNS_lambda)*nabla_v(SPACEDIM,2,:) + LNS_lambda*nabla_v(1,1,:)
  
end subroutine LNS_compute_viscous_stress_tensor





















! ------------------------------------------------------------ !
! compute_gradient_TFSF                                        !
! ------------------------------------------------------------ !
! Computes the gradient of a tensor field and/or of a scalar field SF.
! IN:
!   TF, a tensor field
!   SF, a scalar field
!   swTF, a switch to activate computation for TF
!   swSF, a switch to activate computation for SF
!   swMETHOD, a switch to activate the use of the "desintegration method" for gradient computation methods, and to desactivate to use the SEM definition of the gradient
! OUT:
!   nabla_TF, gradient of TF, (1,1) being dxTFx, (1,2) dzTFx, (2,1) dxTFz, and (2,2) dzTFz. Only makes sense if swTF is active.
!   nabla_SF, gradient of SF, (1) being dxSF, and (2) dzSF. Only makes sense if swSF is active.
! NOTES:
!   See doi:10.1016/j.jcp.2007.12.009, section 4.3.2 for the "desintegration method".
   
subroutine compute_gradient_TFSF(TF, SF, swTF, swSF, swMETHOD, nabla_TF, nabla_SF)
  use constants ! TODO: select variables to use.
  use specfem_par ! TODO: select variables to use.
  use specfem_par_LNS ! TODO: select variables to use. SPACEDIM
  
  implicit none
  
  ! Input/Output.
  real(kind=CUSTOM_REAL), dimension(nglob_DG) :: SF
  real(kind=CUSTOM_REAL), dimension(SPACEDIM, nglob_DG) :: TF
  logical, intent(in) :: swTF, swSF, swMETHOD
  real(kind=CUSTOM_REAL), dimension(SPACEDIM, nglob_DG), intent(out) :: nabla_SF
  real(kind=CUSTOM_REAL), dimension(SPACEDIM, SPACEDIM, nglob_DG), intent(out) :: nabla_TF
  
  ! Local.
  real(kind=CUSTOM_REAL), parameter :: ZEROcr = 0._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: HALFcr = 0.5_CUSTOM_REAL
  real(kind=CUSTOM_REAL) :: SF_P ! When swMETHOD==.true., this variable is used to store the value of the scalar field SF across the element's boundary, in order to compute the flux.
  real(kind=CUSTOM_REAL), dimension(SPACEDIM) :: TF_P, n_out ! When swMETHOD==.true., those variables are used to store the value of the tensor field TF across the element's boundary, in order to compute the flux.
  integer :: ispec,i,j,k,iglob, iglobM, iglobP!, iglob_unique
  real(kind=CUSTOM_REAL) :: flux_n, flux_x, flux_z, nx, nz, weight!rho_DG_P, rhovx_DG_P, rhovz_DG_P, &
        !E_DG_P&!, p_DG_P, &
        !Tx_DG_P, Tz_DG_P, Vxx_DG_P, Vzz_DG_P, Vzx_DG_P, Vxz_DG_P, &
        !timelocal,gamma_P
  logical :: exact_interface_flux
  !integer, dimension(nglob_DG) :: MPI_iglob
  integer, dimension(3) :: neighbor
  integer :: iface1, iface, iface1_neighbor, iface_neighbor, ispec_neighbor
  real(kind=CUSTOM_REAL) :: xixl,xizl,gammaxl,gammazl,jacLoc ! Jacobian matrix and determinant
  !real(kind=CUSTOM_REAL) :: temp_SFx, temp_SFz, temp_TFxx, temp_TFzx, temp_TFxz, temp_TFzz
  !real(kind=CUSTOM_REAL), dimension(nglob_DG) :: &
  !       rho_DG, rhovx_DG, rhovz_DG, E_DG, veloc_x_DG, veloc_z_DG, T, &
  !       nabla_SF(1,:), nabla_SF(SPACEDIM,:), nabla_TF(1,1,:), nabla_TF(SPACEDIM,SPACEDIM,:), nabla_TF(1,SPACEDIM,:), nabla_TF(SPACEDIM,1,:)
  
  !real(kind=CUSTOM_REAL), dimension(SPACEDIM,nglob_DG) :: locgrad_SF
  !real(kind=CUSTOM_REAL), dimension(SPACEDIM,SPACEDIM,nglob_DG) :: locgrad_TF
  
  !real(kind=CUSTOM_REAL) :: dxiTF(1), dgamTF(1), dxiTF(SPACEDIM), dgamTF(SPACEDIM), dSF_dxi, dSF_dgamma ! Derivatives in \Lambda.
  real(kind=CUSTOM_REAL) :: dSF_dxi, dSF_dgamma ! Derivatives in \Lambda.
  real(kind=CUSTOM_REAL), dimension(SPACEDIM) :: dxiTF, dgamTF ! Derivatives in \Lambda.
  !real(kind=CUSTOM_REAL) :: dTFx_dx, dTFx_dz, dTFz_dx, dTFz_dz, dSF_dx, dSF_dz
  real(kind=CUSTOM_REAL) :: wxl, wzl ! Quadrature weights.
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: temp_SFx_1, temp_SFx_2, &
        temp_SFz_1, temp_SFz_2, temp_TFxx_1, temp_TFxx_2, &
        temp_TFxz_1, temp_TFxz_2, temp_TFzx_1, temp_TFzx_2, temp_TFzz_1, temp_TFzz_2
  !real(kind=CUSTOM_REAL) :: vx_init, vz_init
  !real(kind=CUSTOM_REAL) :: ya_x_l, ya_z_l
  
  if(swMETHOD) then
    ! "Desintegrate."
    veloc_x_DG = rhovx_DG/rho_DG
    !veloc_z_DG = rhovz_DG/rho_DG
    !T = (E_DG/rho_DG - HALFcr*(veloc_x_DG**2 + veloc_z_DG**2))/c_V
  endif
  
  if(swSF) then
    nabla_SF = ZEROcr
  endif
  if(swTF) then
    nabla_TF = ZEROcr
  endif

  do ispec = 1, nspec ! Loop over elements.
    ! acoustic spectral element
    !if (ispec_is_acoustic(ispec)) then
    if (ispec_is_acoustic_DG(ispec)) then
      
      ! --------------------------- !
      ! Volume terms.               !
      ! --------------------------- !
      do j = 1, NGLLZ
        do i = 1, NGLLX
          iglob = ibool_DG(i,j,ispec)
          jacLoc = jacobian(i,j,ispec)
          xixl = xix(i,j,ispec)
          xizl = xiz(i,j,ispec)
          gammaxl = gammax(i,j,ispec)
          gammazl = gammaz(i,j,ispec)
          
          !if(ABC_STRETCH .and. stretching_buffer(ibool_before_perio(i,j,ispec))>0) then
          !  ! See beginning of subroutine compute_forces_acoustic_DG for detailed explanations.
          !  iglob_unique=ibool_before_perio(i,j,ispec)
          !  ya_x_l=stretching_ya(1, iglob_unique)
          !  ya_z_l=stretching_ya(2, iglob_unique)
          !  xixl = ya_x_l * xixl
          !  xizl = ya_z_l * xizl
          !  gammaxl = ya_x_l * gammaxl
          !  gammazl = ya_z_l * gammazl
          !  ! TODO: Do that more clearly.
          !  jacLoc = ya_x_l*ya_z_l*jacLoc
          !endif
          
          wzl = wzgll(j)
          wxl = wxgll(i)
          
          if(swMETHOD) then
            ! In that case, we want to compute \int_{\Omega} (\nabla f)\Phi as:
            ! - \int_{\Omega} T\nabla\Phi + \int_{\partial\Omega} f\Phi.
            ! The idea is to store in:
            !   nabla_SF(1,:),
            !   nabla_SF(SPACEDIM,:),
            !   nabla_TF(1,1,:),
            !   nabla_TF(1,SPACEDIM,:),
            !   nabla_TF(SPACEDIM,1,:),
            !   nabla_TF(SPACEDIM,SPACEDIM,:)
            ! the values of the approximated integrals:
            !  $\int \partial_xT\Phi_x$,
            !  $\int \partial_zT\Phi_z$,
            !  $\int \partial_xVx\Phi_x$,
            !  etc.
            ! and "desintegrate" those values by multiplying the obtained vector by the inverse mass matrix. As explained before, the integrals are computed by using the divergence theorem and taking into account the surface terms.
            
            ! Compute inner contributions.
            !if(.not. CONSTRAIN_HYDROSTATIC) then
            if(swSF) then
              temp_SFx_1(i,j)  = wzl * jacLoc * (xixl * SF(iglob))
              temp_SFz_1(i,j)  = wzl * jacLoc * (xizl * SF(iglob))
              temp_SFx_2(i,j)  = wxl * jacLoc * (gammaxl * SF(iglob))
              temp_SFz_2(i,j)  = wxl * jacLoc * (gammazl * SF(iglob))
            endif
            if(swTF) then
              temp_TFxx_1(i,j) = wzl * jacLoc * (xixl * TF(1,iglob))
              temp_TFxz_1(i,j) = wzl * jacLoc * (xizl * TF(1,iglob))
              temp_TFzx_1(i,j) = wzl * jacLoc * (xixl * TF(SPACEDIM,iglob))
              temp_TFzz_1(i,j) = wzl * jacLoc * (xizl * TF(SPACEDIM,iglob))
              temp_TFxx_2(i,j) = wxl * jacLoc * (gammaxl * TF(1,iglob))
              temp_TFxz_2(i,j) = wxl * jacLoc * (gammazl * TF(1,iglob))
              temp_TFzx_2(i,j) = wxl * jacLoc * (gammaxl * TF(SPACEDIM,iglob))
              temp_TFzz_2(i,j) = wxl * jacLoc * (gammazl * TF(SPACEDIM,iglob))
            endif
            !else
            !  vx_init = rhovx_init(iglob)/rho_init(iglob)
            !  vz_init = rhovz_init(iglob)/rho_init(iglob)
            !  
            !  temp_SFx_1(i,j)  = wzl * jacLoc * (xixl * (SF(iglob) - T_init(iglob)))
            !  temp_SFz_1(i,j)  = wzl * jacLoc * (xizl * (SF(iglob) - T_init(iglob)))
            !  temp_TFxx_1(i,j) = wzl * jacLoc * (xixl * (TF(1,iglob) - vx_init))
            !  temp_TFxz_1(i,j) = wzl * jacLoc * (xizl * (TF(1,iglob))) ! Some hypothesis on initial velocity is used here.
            !  temp_TFzx_1(i,j) = wzl * jacLoc * (xixl * (TF(SPACEDIM,iglob) - vz_init))
            !  temp_TFzz_1(i,j) = wzl * jacLoc * (xizl * (TF(SPACEDIM,iglob) - vz_init))
            !  
            !  temp_SFx_2(i,j)  = wxl * jacLoc * (gammaxl * (SF(iglob) - T_init(iglob)))
            !  temp_SFz_2(i,j)  = wxl * jacLoc * (gammazl * (SF(iglob) - T_init(iglob)))
            !  temp_TFxx_2(i,j) = wxl * jacLoc * (gammaxl * (TF(1,iglob) - vx_init))
            !  temp_TFxz_2(i,j) = wxl * jacLoc * (gammazl * (TF(1,iglob))) ! Some hypothesis on initial velocity is used here.
            !  temp_TFzx_2(i,j) = wxl * jacLoc * (gammaxl * (TF(SPACEDIM,iglob) - vz_init))
            !  temp_TFzz_2(i,j) = wxl * jacLoc * (gammazl * (TF(SPACEDIM,iglob) - vz_init))
            !endif
          else
            ! In that case, we want to compute \int_{\Omega} (\nabla f)\Phi directly as it.
            ! The idea is to store in:
            !   nabla_SF(1,:),
            !   nabla_SF(SPACEDIM,:),
            !   nabla_TF(1,1,:),
            !   nabla_TF(1,SPACEDIM,:),
            !   nabla_TF(SPACEDIM,1,:),
            !   nabla_TF(SPACEDIM,SPACEDIM,:)
            ! the actual values of the quantities:
            !   \partial_xT,
            !   \partial_zT,
            !   \partial_xVx,
            !   etc.
            ! which is immediate through the SEM formulation.
            
            if(swSF) then
              dSF_dxi     = ZEROcr
              dSF_dgamma  = ZEROcr
            endif
            if(swTF) then
              dxiTF  = ZEROcr
              dgamTF = ZEROcr
            endif
            
            ! Compute derivatives in unit element \Lambda.
            ! Note: we can merge the two loops because NGLLX=NGLLZ.
            do k = 1, NGLLX
              !if(.not. CONSTRAIN_HYDROSTATIC) then
              if(swSF) then
                dSF_dxi     = dSF_dxi     + SF(ibool_DG(k,j,ispec)) * real(hprime_xx(i,k), kind=CUSTOM_REAL)
                dSF_dgamma  = dSF_dgamma  + SF(ibool_DG(i,k,ispec)) * real(hprime_zz(j,k), kind=CUSTOM_REAL)
              endif
              if(swTF) then
                dxiTF(1)         = dxiTF(1)         + TF(1,ibool_DG(k,j,ispec)) * real(hprime_xx(i,k), kind=CUSTOM_REAL)
                dgamTF(1)        = dgamTF(1)        + TF(1,ibool_DG(i,k,ispec)) * real(hprime_zz(j,k), kind=CUSTOM_REAL)
                dxiTF(SPACEDIM)  = dxiTF(SPACEDIM)  + TF(SPACEDIM,ibool_DG(k,j,ispec)) * real(hprime_xx(i,k), kind=CUSTOM_REAL)
                dgamTF(SPACEDIM) = dgamTF(SPACEDIM) + TF(SPACEDIM,ibool_DG(i,k,ispec)) * real(hprime_zz(j,k), kind=CUSTOM_REAL)
              endif
              !else
              !  vx_init = rhovx_init(ibool_DG(k,j,ispec))/rho_init(ibool_DG(k,j,ispec))
              !  dxiTF(1) = dxiTF(1) + (TF(1,ibool_DG(k,j,ispec)) - vx_init) &
              !                      * real(hprime_xx(i,k), kind=CUSTOM_REAL)
              !  vx_init = rhovx_init(ibool_DG(i,k,ispec))/rho_init(ibool_DG(i,k,ispec))       
              !  dgamTF(1) = dgamTF(1) + (TF(1,ibool_DG(i,k,ispec)) - vx_init) &
              !                            * real(hprime_zz(j,k), kind=CUSTOM_REAL)
              !  vz_init = rhovz_init(ibool_DG(k,j,ispec))/rho_init(ibool_DG(k,j,ispec))
              !  dxiTF(SPACEDIM) = dxiTF(SPACEDIM) + (TF(SPACEDIM,ibool_DG(k,j,ispec)) - vz_init) &
              !                      * real(hprime_xx(i,k), kind=CUSTOM_REAL)
              !  vz_init = rhovz_init(ibool_DG(i,k,ispec))/rho_init(ibool_DG(i,k,ispec))
              !  dgamTF(SPACEDIM) = dgamTF(SPACEDIM) + (TF(SPACEDIM,ibool_DG(i,k,ispec)) - vz_init) &
              !                            * real(hprime_zz(j,k), kind=CUSTOM_REAL)
              !  dSF_dxi = dSF_dxi + (SF(ibool_DG(k,j,ispec)) - T_init(ibool_DG(k,j,ispec))) &
              !                    * real(hprime_xx(i,k), kind=CUSTOM_REAL)
              !  dSF_dgamma = dSF_dgamma + (SF(ibool_DG(i,k,ispec)) - T_init(ibool_DG(i,k,ispec))) &
              !                           * real(hprime_zz(j,k), kind=CUSTOM_REAL)
              !endif
            enddo ! Enddo on k.

            ! Compute derivatives in element \Omega using the chain rule.
            if(swSF) then
              !dSF_dx   = dSF_dxi * xixl + dSF_dgamma * gammaxl
              !dSF_dz   = dSF_dxi * xizl + dSF_dgamma * gammazl
              !temp_SFx = dSF_dx * jacLoc
              !temp_SFz = dSF_dz * jacLoc
              !nabla_SF(1,iglob) = dSF_dx
              !nabla_SF(SPACEDIM,iglob) = dSF_dz
              nabla_SF(1,iglob) = dSF_dxi * xixl + dSF_dgamma * gammaxl
              nabla_SF(SPACEDIM,iglob) = dSF_dxi * xizl + dSF_dgamma * gammazl
            endif
            if(swTF) then
              !dTFx_dx   = dxiTF(1) * xixl + dgamTF(1) * gammaxl
              !dTFx_dz   = dxiTF(1) * xizl + dgamTF(1) * gammazl
              !dTFz_dx   = dxiTF(SPACEDIM) * xixl + dgamTF(SPACEDIM) * gammaxl
              !dTFz_dz   = dxiTF(SPACEDIM) * xizl + dgamTF(SPACEDIM) * gammazl
              !temp_TFxx = dTFx_dx * jacLoc
              !temp_TFxz = dTFx_dz * jacLoc
              !temp_TFzx = dTFz_dx * jacLoc
              !temp_TFzz = dTFz_dz * jacLoc
              !nabla_TF(1,1,iglob) = dTFx_dx
              !nabla_TF(1,SPACEDIM,iglob) = dTFx_dz
              !nabla_TF(SPACEDIM,1,iglob) = dTFz_dx
              !nabla_TF(SPACEDIM,SPACEDIM,iglob) = dTFz_dz
              do k=1,SPACEDIM ! Re-using index k for spatial dimension.
                nabla_TF(k,1,iglob) = dxiTF(k) * xixl + dgamTF(k) * gammaxl
                nabla_TF(k,SPACEDIM,iglob) = dxiTF(k) * xizl + dgamTF(k) * gammazl
                !nabla_TF(1,1,iglob) = dxiTF(1) * xixl + dgamTF(1) * gammaxl
                !nabla_TF(1,SPACEDIM,iglob) = dxiTF(1) * xizl + dgamTF(1) * gammazl
                !nabla_TF(SPACEDIM,1,iglob) = dxiTF(SPACEDIM) * xixl + dgamTF(SPACEDIM) * gammaxl
                !nabla_TF(SPACEDIM,SPACEDIM,iglob) = dxiTF(SPACEDIM) * xizl + dgamTF(SPACEDIM) * gammazl
              enddo
            endif
          endif
        enddo
      enddo
      
      if(swMETHOD) then
        ! "Desintegrate".
        do j = 1, NGLLZ
          do i = 1, NGLLX
            iglob = ibool_DG(i,j,ispec)
            ! Assemble the contributions using the inner contributions.
            do k = 1, NGLLX
              if(swSF) then
              nabla_SF(1,iglob) = nabla_SF(1,iglob) - &
                               (temp_SFx_1(k,j) * real(hprimewgll_xx(k,i), kind=CUSTOM_REAL) + &
                                temp_SFx_2(i,k) * real(hprimewgll_zz(k,j), kind=CUSTOM_REAL))
              nabla_SF(SPACEDIM,iglob) = nabla_SF(SPACEDIM,iglob) - &
                               (temp_SFz_1(k,j) * real(hprimewgll_xx(k,i), kind=CUSTOM_REAL) + &
                                temp_SFz_2(i,k) * real(hprimewgll_zz(k,j), kind=CUSTOM_REAL))
              endif
              if(swTF) then
              nabla_TF(1,1,iglob) = nabla_TF(1,1,iglob) - &
                                (temp_TFxx_1(k,j) * real(hprimewgll_xx(k,i), kind=CUSTOM_REAL) + &
                                 temp_TFxx_2(i,k) * real(hprimewgll_zz(k,j), kind=CUSTOM_REAL))
              nabla_TF(1,SPACEDIM,iglob) = nabla_TF(1,SPACEDIM,iglob) - &
                                (temp_TFxz_1(k,j) * real(hprimewgll_xx(k,i), kind=CUSTOM_REAL) + &
                                 temp_TFxz_2(i,k) * real(hprimewgll_zz(k,j), kind=CUSTOM_REAL))
              nabla_TF(SPACEDIM,1,iglob) = nabla_TF(SPACEDIM,1,iglob) - &
                                (temp_TFzx_1(k,j) * real(hprimewgll_xx(k,i), kind=CUSTOM_REAL) + &
                                 temp_TFzx_2(i,k) * real(hprimewgll_zz(k,j), kind=CUSTOM_REAL))
              nabla_TF(SPACEDIM,SPACEDIM,iglob) = nabla_TF(SPACEDIM,SPACEDIM,iglob) - &
                                (temp_TFzz_1(k,j) * real(hprimewgll_xx(k,i), kind=CUSTOM_REAL) + &
                                 temp_TFzz_2(i,k) * real(hprimewgll_zz(k,j), kind=CUSTOM_REAL))
              endif
            enddo
          enddo
        enddo
        
        ! Interface terms.
        do  iface = 1, 4 
          do  iface1 = 1, NGLLX
            i = link_iface_ijispec(iface1,iface,ispec,1)
            j = link_iface_ijispec(iface1,iface,ispec,2)
            ! Interior point
            iglobM = ibool_DG(i,j,ispec)
            !rho_DG_P     = ZEROcr
            !rhovx_DG_P   = ZEROcr
            !rhovz_DG_P   = ZEROcr
            !E_DG_P       = ZEROcr
            SF_P         = ZEROcr
            TF_P         = ZEROcr
            !p_DG_P       = ZEROcr
            !Vxx_DG_P     = ZEROcr
            !Vzz_DG_P     = ZEROcr
            !Vzx_DG_P     = ZEROcr
            !Vxz_DG_P     = ZEROcr
            nx = nx_iface(iface, ispec)
            nz = nz_iface(iface, ispec)
            
            weight = weight_iface(iface1,iface,ispec)
            neighbor = -1
            if(neighbor_DG_iface(iface1, iface, ispec, 3) > -1) then
              iface1_neighbor = neighbor_DG_iface(iface1, iface, ispec, 1)
              iface_neighbor  = neighbor_DG_iface(iface1, iface, ispec, 2)
              ispec_neighbor = neighbor_DG_iface(iface1, iface, ispec, 3)
              neighbor(1) = link_iface_ijispec(iface1_neighbor, iface_neighbor, ispec_neighbor,1)
              neighbor(2) = link_iface_ijispec(iface1_neighbor, iface_neighbor, ispec_neighbor,2)
              neighbor(3) = ispec_neighbor
            endif
          
            !if(ABC_STRETCH .and. stretching_buffer(ibool_before_perio(i,j,ispec))>0) then
            !  ! Update flux with stretching components. See explanation in the surface terms part in the subroutine above.
            !  ! TO DO: Do that more clearly.
            !  iglob_unique=ibool_before_perio(i,j,ispec)
            !  weight=stretching_ya(1,iglob_unique)*stretching_ya(2,iglob_unique)*weight
            !endif
            
            iglobP = 1
            if(neighbor(1) > -1) then
              iglobP = ibool_DG(neighbor(1), neighbor(2), neighbor(3))
            endif
            
            exact_interface_flux = .false. ! Reset this variable to .false.: by default, the fluxes have to be computed (jump!=0). In some specific cases (assigned during the call to LNS_get_interfaces_unknowns), the flux can be exact (jump==0).
          ! TODO: dedicated routine.
          !call compute_interface_unknowns(i,j,ispec, drho_P, rho0dv_P(1), &
          !        rho0dv_P(2), dE_P, TF_P(1), TF_P(SPACEDIM), in_dp_P, T_P, &
          !        Tx_DG_P, Tz_DG_P, Vxx_DG_P, Vzz_DG_P, Vzx_DG_P, Vxz_DG_P, gamma_P,&
          !        neighbor, &
          !        exact_interface_flux, &
          !        cv_drho(iglobM), cv_dE(iglobM), cv_rho0dv(1,iglobM), cv_rho0dv(SPACEDIM,iglobM), &
          !        V_DG(:,:,iglobM), T_DG(:,iglobM), &
          !        cv_drho(iglobP), cv_dE(iglobP), cv_rho0dv(1,iglobP), cv_rho0dv(SPACEDIM,iglobP), &
          !        V_DG(:,:,iglobP), T_DG(:,iglobP), &
          !        nx, nz, weight, currentTime, iface1, iface)
          !        !TEST STRETCH
          !        !nx_unit, nz_unit, weight, currentTime, iface1, iface)
          !call LNS_get_interfaces_unknowns(i, j, ispec, &
          !        drho_P, rho0dv_P(1), &
          !        rho0dv_P(2), dE_P, TF_P(1), TF_P(SPACEDIM), in_dp_P, T_P, &
          !        Tx_DG_P, Tz_DG_P, Vxx_DG_P, Vzz_DG_P, Vzx_DG_P, Vxz_DG_P, gamma_P,&
          !        neighbor, &
          !        exact_interface_flux, &
          !        cv_drho(iglobM), cv_dE(iglobM), cv_rho0dv(1,iglobM), cv_rho0dv(SPACEDIM,iglobM), &
          !        V_DG(:,:,iglobM), T_DG(:,iglobM), &
          !        cv_drho(iglobP), cv_dE(iglobP), cv_rho0dv(1,iglobP), cv_rho0dv(SPACEDIM,iglobP), &
          !        V_DG(:,:,iglobP), T_DG(:,iglobP), &
          !        nx, nz, weight, currentTime, iface1, iface)
            
            call LNS_get_interfaces_unknowns(i, j, ispec, iface1, iface, neighbor, LNS_dummy_1d(1), & ! Point identifier (input).
                  LNS_dummy_1d(1), LNS_dummy_2d(:,1), LNS_dummy_1d(1), & ! Input constitutive variables, "M" side.
                  LNS_dummy_1d(1), LNS_dummy_2d(:,1), LNS_dummy_1d(1), & ! Input constitutive variables, "P" side.
                  LNS_dummy_1d(1), & ! Input other variable, "M" side.
                  !V_DG(:,:,iglobM), T_DG(:,iglobM), & ! Input derivatives, "M" side. MIGHT NEED.
                  !V_DG(:,:,iglobP), T_DG(:,iglobP), & ! Input derivatives, "M" side. MIGHT NEED.
                  n_out, & ! Normal vector (input).
                  exact_interface_flux, & ! Switch to disable jump in some cases (output).
                  LNS_dummy_1d(1), LNS_dummy_2d(:,1), LNS_dummy_1d(1), & ! Output constitutive variables.
                  !Tx_DG_P, Tz_DG_P, Vxx_DG_P, Vzz_DG_P, Vzx_DG_P, Vxz_DG_P,& ! Output derivatives. MIGHT NEED.
                  LNS_dummy_2d(:,1), LNS_dummy_1d(1), LNS_dummy_2d(:,1), LNS_dummy_2d(:,1), LNS_dummy_1d(1:3), & ! Output other variables.
                  SF_P, TF_P, .true.) ! We know that the scalar field will always be temperature, and that the tensor field will always be the velocity, so we use the already-built LNS_get_interfaces_unknowns routine to get them. The switch is here to prevent unecessary quantites to be computed during that specific call. TODO: if other fields need to be computed, one will have to use a dedicated routine.

            !vx_init = rhovx_init(iglobM)/rho_init(iglobM)
            !vz_init = rhovz_init(iglobM)/rho_init(iglobM)

            ! Dot products.
            if(swSF) then
              flux_x = SF(iglobM) + SF_P ! Once multiplied with HALFcr below, will represent flux along x of the scalar field SF.
              !if(CONSTRAIN_HYDROSTATIC) flux_x = flux_x - 2*T_init(iglobM)
              flux_n = flux_x*n_out(1) ! Recall: n_out(1)=n_x.
              nabla_SF(1,iglobM) = nabla_SF(1,iglobM) + weight*flux_n*HALFcr
              flux_z = SF(iglobM) + SF_P ! Once multiplied with HALFcr below, will represent flux along z of the scalar field SF.
              !if(CONSTRAIN_HYDROSTATIC) flux_z = flux_z - 2*T_init(iglobM)
              flux_n = flux_z*n_out(SPACEDIM) ! Recall: n_out(SPACEDIM)=n_z.
              nabla_SF(SPACEDIM,iglobM) = nabla_SF(SPACEDIM,iglobM) + weight*flux_n*HALFcr
            endif
            if(swTF) then
              flux_x = TF(1,iglobM) + TF_P(1) ! Once multiplied with HALFcr below, will represent flux along x of the x-component of the tensor field TF.
              !if(CONSTRAIN_HYDROSTATIC) flux_x = flux_x - 2*vx_init
              flux_n = flux_x*n_out(1) ! Recall: n_out(1)=n_x.
              nabla_TF(1,1,iglobM) = nabla_TF(1,1,iglobM) + weight*flux_n*HALFcr
              flux_z = TF(1,iglobM) + TF_P(1) ! Once multiplied with HALFcr below, will represent flux along z of the x-component of the tensor field TF.
              !if(CONSTRAIN_HYDROSTATIC) flux_z = flux_z - 2*vx_init
              flux_n = flux_z*n_out(SPACEDIM) ! Recall: n_out(SPACEDIM)=n_z.
              nabla_TF(1,SPACEDIM,iglobM) = nabla_TF(1,SPACEDIM,iglobM) + weight*flux_n*HALFcr
              flux_x = TF(SPACEDIM,iglobM) + TF_P(SPACEDIM) ! Once multiplied with HALFcr below, will represent flux along x of the z-component of the tensor field TF.
              !if(CONSTRAIN_HYDROSTATIC) flux_x = flux_x - 2*vz_init
              flux_n = flux_x*n_out(1) ! Recall: n_out(1)=n_x.
              nabla_TF(SPACEDIM,1,iglobM) = nabla_TF(SPACEDIM,1,iglobM) + weight*flux_n*HALFcr
              flux_z = TF(SPACEDIM,iglobM) + TF_P(SPACEDIM) ! Once multiplied with HALFcr below, will represent flux along z of the z-component of the tensor field TF.
              !if(CONSTRAIN_HYDROSTATIC) flux_z = flux_z - 2*vz_init
              flux_n = flux_z*n_out(SPACEDIM) ! Recall: n_out(SPACEDIM)=n_z.
              nabla_TF(SPACEDIM,SPACEDIM,iglobM) = nabla_TF(SPACEDIM,SPACEDIM,iglobM) + weight*flux_n*HALFcr
            endif
          enddo ! Enddo on iface1.
        enddo ! Enddo on iface.
      endif ! Endif on swMETHOD.
    endif ! End of test if acoustic element.
  enddo ! Enddo on ispec.
  
  if(swMETHOD) then
    ! "Desintegrate".
    if(swSF) then
      nabla_SF(1,:)  = nabla_SF(1,:) * rmass_inverse_acoustic_DG(:)
      nabla_SF(SPACEDIM,:)  = nabla_SF(SPACEDIM,:) * rmass_inverse_acoustic_DG(:)
    endif
    if(swTF) then
      do i=1,SPACEDIM
        nabla_TF(i,1,:) = nabla_TF(i,1,:) * rmass_inverse_acoustic_DG(:)
        nabla_TF(i,SPACEDIM,:) = nabla_TF(i,SPACEDIM,:) * rmass_inverse_acoustic_DG(:)
        !nabla_TF(1,1,:) = nabla_TF(1,1,:) * rmass_inverse_acoustic_DG(:)
        !nabla_TF(1,SPACEDIM,:) = nabla_TF(1,SPACEDIM,:) * rmass_inverse_acoustic_DG(:)
        !nabla_TF(SPACEDIM,1,:) = nabla_TF(SPACEDIM,1,:) * rmass_inverse_acoustic_DG(:)
        !nabla_TF(SPACEDIM,SPACEDIM,:) = nabla_TF(SPACEDIM,SPACEDIM,:) * rmass_inverse_acoustic_DG(:)
      enddo
    endif
  endif
  
  !! Store in variables intended for output.
  !if(swSF) then
  !  nabla_SF(1, :) = nabla_SF(1,:)
  !  nabla_SF(2, :) = nabla_SF(SPACEDIM,:)
  !endif
  !if(swTF) then
  !  nabla_TF(1, 1, :) = nabla_TF(1,1,:) ! dxTFx
  !  nabla_TF(1, SPACEDIM, :) = nabla_TF(1,SPACEDIM,:) ! dzTFx
  !  nabla_TF(SPACEDIM, 1, :) = nabla_TF(SPACEDIM,1,:) ! dxTFz
  !  nabla_TF(SPACEDIM, SPACEDIM, :) = nabla_TF(SPACEDIM,SPACEDIM,:) ! dzTFz
  !endif
end subroutine compute_gradient_TFSF




















! ------------------------------------------------------------ !
! compute_p                                                    !
! ------------------------------------------------------------ !
! Computes pressure from constitutive variables.
subroutine compute_p(in_rho, in_v, in_E, out_p)
  use constants, only: CUSTOM_REAL
  use specfem_par, only: gammaext_DG, nglob_DG
  use specfem_par_LNS, only: SPACEDIM
  implicit none
  ! Input/Output.
  real(kind=CUSTOM_REAL), dimension(nglob_DG), intent(in) :: in_rho, in_E
  real(kind=CUSTOM_REAL), dimension(SPACEDIM, nglob_DG), intent(in) :: in_v
  real(kind=CUSTOM_REAL), dimension(nglob_DG), intent(out) :: out_p
  ! Local.
  real(kind=CUSTOM_REAL), parameter :: ONEcr = 1._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: HALFcr = 0.5_CUSTOM_REAL
  out_p =   (gammaext_DG - ONEcr) &
          * (in_E - HALFcr * in_rho * (  in_v(1,:)**2 &
                                       + in_v(SPACEDIM,:)**2))
end subroutine compute_p
! ------------------------------------------------------------ !
! compute_dp                                                   !
! ------------------------------------------------------------ !
! Computes pressure perturbation from constitutive variables.
! Note: this could have been done nearly inline by using the subroutine compute_p, but defining this function enables one to use less RAM.
subroutine compute_dp(in_rho, in_v, in_E, out_dp)
  use constants, only: CUSTOM_REAL
  use specfem_par, only: gammaext_DG, nglob_DG
  use specfem_par_LNS, only: LNS_p0, SPACEDIM
  implicit none
  ! Input/Output.
  real(kind=CUSTOM_REAL), dimension(nglob_DG), intent(in) :: in_rho, in_E
  real(kind=CUSTOM_REAL), dimension(SPACEDIM, nglob_DG), intent(in) :: in_v
  real(kind=CUSTOM_REAL), dimension(nglob_DG), intent(out) :: out_dp
  ! Local.
  real(kind=CUSTOM_REAL), parameter :: ONEcr = 1._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: HALFcr = 0.5_CUSTOM_REAL
  out_dp =   (gammaext_DG - ONEcr) &
           * (in_E - HALFcr * in_rho * (  in_v(1,:)**2 &
                                        + in_v(SPACEDIM,:)**2)) &
           - LNS_p0
end subroutine compute_dp
! ------------------------------------------------------------ !
! compute_dp_i                                                 !
! ------------------------------------------------------------ !
! Same as compute_dp, but point by point (unvectorised).
subroutine compute_dp_i(in_rho, in_v, in_E, out_p, iglob)
  use constants, only: CUSTOM_REAL
  use specfem_par, only: gammaext_DG
  use specfem_par_LNS, only: LNS_p0, SPACEDIM
  implicit none
  ! Input/Output.
  real(kind=CUSTOM_REAL), intent(in) :: in_rho, in_E
  real(kind=CUSTOM_REAL), dimension(SPACEDIM), intent(in) :: in_v
  integer, intent(in) :: iglob
  real(kind=CUSTOM_REAL), intent(out) :: out_p
  ! Local.
  real(kind=CUSTOM_REAL), parameter :: ONEcr = 1._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: HALFcr = 0.5_CUSTOM_REAL
  out_p =   (gammaext_DG(iglob) - ONEcr) &
          * (in_E - HALFcr * in_rho * (  in_v(1)**2 &
                                       + in_v(SPACEDIM)**2)) &
          - LNS_p0(iglob)
end subroutine compute_dp_i
! ------------------------------------------------------------ !
! compute_T                                                    !
! ------------------------------------------------------------ !
! Computes temperature from constitutive variables.
subroutine compute_T(in_rho, in_v, in_E, out_T)
  use constants, only: CUSTOM_REAL
  use specfem_par, only: c_V, nglob_DG
  use specfem_par_LNS, only: SPACEDIM
  implicit none
  ! Input/Output.
  real(kind=CUSTOM_REAL), dimension(nglob_DG), intent(in) :: in_rho, in_E
  real(kind=CUSTOM_REAL), dimension(SPACEDIM, nglob_DG), intent(in) :: in_v
  real(kind=CUSTOM_REAL), dimension(nglob_DG), intent(out) :: out_T
  ! Local.
  real(kind=CUSTOM_REAL), parameter :: HALFcr = 0.5_CUSTOM_REAL
  out_T=(in_E/in_rho - HALFcr*(in_v(1,:)**2+in_v(SPACEDIM,:)**2))/c_V
end subroutine compute_T
! ------------------------------------------------------------ !
! compute_T_i                                                  !
! ------------------------------------------------------------ !
! Same as compute_T, but point by point (unvectorised).
subroutine compute_T_i(in_rho, in_v, in_E, out_T)
  use constants, only: CUSTOM_REAL
  use specfem_par, only: c_V
  use specfem_par_LNS, only: SPACEDIM
  implicit none
  ! Input/Output.
  real(kind=CUSTOM_REAL), intent(in) :: in_rho, in_E
  real(kind=CUSTOM_REAL), dimension(SPACEDIM), intent(in) :: in_v
  real(kind=CUSTOM_REAL), intent(out) :: out_T
  ! Local.
  real(kind=CUSTOM_REAL), parameter :: HALFcr = 0.5_CUSTOM_REAL
  out_T=(in_E/in_rho - HALFcr*(in_v(1)**2+in_v(SPACEDIM)**2))/c_V
end subroutine compute_T_i
! ------------------------------------------------------------ !
! compute_dT                                                   !
! ------------------------------------------------------------ !
! Computes temperature perturbation from constitutive variables.
! Note: this could have been done nearly inline by using the subroutine compute_T, but defining this function enables one to use less RAM.
subroutine compute_dT(in_rho, in_v, in_E, out_dT)
  use constants, only: CUSTOM_REAL
  use specfem_par, only: c_V, nglob_DG
  use specfem_par_LNS, only: LNS_T0, SPACEDIM
  implicit none
  ! Input/Output.
  real(kind=CUSTOM_REAL), dimension(nglob_DG), intent(in) :: in_rho, in_E
  real(kind=CUSTOM_REAL), dimension(SPACEDIM, nglob_DG), intent(in) :: in_v
  real(kind=CUSTOM_REAL), dimension(nglob_DG), intent(out) :: out_dT
  ! Local.
  real(kind=CUSTOM_REAL), parameter :: HALFcr = 0.5_CUSTOM_REAL
  out_dT=  (in_E/in_rho - HALFcr*(in_v(1,:)**2+in_v(SPACEDIM,:)**2))/c_V &
         - LNS_T0
end subroutine compute_dT

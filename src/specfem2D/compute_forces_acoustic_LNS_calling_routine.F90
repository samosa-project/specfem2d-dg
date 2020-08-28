! ------------------------------------------------------------ !
! compute_forces_acoustic_LNS_main                             !
! ------------------------------------------------------------ !
! Main routine taking care of the DG elements to be solved with the LNS (Linear Navier-Stokes) module.
! The initialisation of the background parameters (see routine 'initial_state_LNS()' below) is done before the call to this routine, in
! 'prepare_timerun.f90'.
! This routine is called by iterate_time().

subroutine compute_forces_acoustic_LNS_main()
  use constants
  use specfem_par
  use specfem_par_LNS

  implicit none

  ! Local.
  real(kind=CUSTOM_REAL), parameter :: ZEROcr = 0._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: ONEcr = 1._CUSTOM_REAL
  real(kind=CUSTOM_REAL) :: timelocal
  real(kind=CUSTOM_REAL), dimension(NDIM, NDIM, nglob_DG) :: nabla_dv
  integer :: ier, i_aux
  logical check_linear_hypothesis, check_initial_conditions
  
  ! Checks if anything has to be done.
  if(.not. any_acoustic_DG) return
  
  ! Switch dictating whether the linear hypothesis should be checked at every iteration (set to .true.) or not (set to .false).
  check_linear_hypothesis = .true.
  ! Switch dictating whether the initial conditions should be checked at the first iteration (set to .true.) or not (set to .false).
  check_initial_conditions = .false.
  
  ! Compute current time.
  timelocal = (it-1)*deltat + LNS_scheme_C(i_stage)*deltat
  
  ! Intialisation.
  if(it == 1 .and. i_stage == 1) then
    ! Initialise state registers.
    ! Note: since constitutive variables are perturbations, they are necessarily zero at start.
    ! Allocations occur in 'prepare_timerun_wavefields.f90'.
    LNS_drho   = ZEROcr
    LNS_rho0dv = ZEROcr
    LNS_dE     = ZEROcr
    
    ! Initialise state auxiliary quantities.
    LNS_dT     = ZEROcr
    LNS_dp     = ZEROcr
    LNS_dv     = ZEROcr
    LNS_dm     = ZEROcr
    if(LNS_viscous) then
      sigma_dv   = ZEROcr
      nabla_dT   = ZEROcr
    endif
    
    ! Initialise auxiliary registers (for RK).
    aux_drho   = ZEROcr
    aux_rho0dv = ZEROcr
    aux_dE     = ZEROcr
    RHS_drho   = ZEROcr
    RHS_rho0dv = ZEROcr
    RHS_dE     = ZEROcr
    if(LNS_avib) then
      ! Those variables are deallocated before even starting time iterations in initialise_LNS below.
      LNS_e1 = ZEROcr
      aux_e1 = ZEROcr
      RHS_e1 = ZEROcr
    endif
    
    ! Prepare MPI buffers.
#ifdef USE_MPI
    if(NPROC>1) then
      call prepare_MPI_DG()
    endif
#endif
    
    ! Allocate acoustic coupling array.
    allocate(ispec_is_acoustic_coupling_ac(nglob_DG), stat=ier)
    if (ier /= 0) then
      call exit_MPI(myrank, "Error allocating 'ispec_is_acoustic_coupling_ac' arrays.")
    endif
    ispec_is_acoustic_coupling_ac = -1
    if(.not. only_DG_acoustic) then
      call find_DG_acoustic_coupling()
    endif
  
    ! Check initial conditions.
    if(check_initial_conditions) call LNS_prevent_nonsense()

#ifdef USE_MPI
    ! Fill the MPI buffers related to background model, but only once at the first iteration.
    call LNS_fill_MPI_buffers_var1_bckg0(.false.)
#endif
    
  endif ! Endif on (it == 1) and (i_stage == 1).
  
  if(LNS_VERBOSE>=1 .and. myrank == 0 .and. mod(it, LNS_MODPRINT)==0) then
    WRITE(*, *) "****************************************************************"
    WRITE(*, "(a,i9,a,i1,a,e23.16,a)") "Iteration", it, ", stage ", i_stage, ", local time", timelocal, "."
  endif
  
  ! Precompute momentum perturbation, velocity perturbation,
  ! temperature perturbation, and pressure perturbation.
  LNS_dv = ZEROcr
  do i_aux = 1, NDIM
    LNS_dm(i_aux, :) = LNS_rho0dv(i_aux, :) + LNS_drho*LNS_v0(i_aux, :)
    where(LNS_rho0/=ZEROcr) LNS_dv(i_aux, :) = LNS_rho0dv(i_aux, :)/LNS_rho0 ! 'where(...)' as safeguard, as rho0=0 should never happen.
  enddo
  call compute_dp(LNS_rho0+LNS_drho, LNS_v0+LNS_dv, LNS_E0+LNS_dE, LNS_dp)
  call compute_dT(LNS_rho0+LNS_drho, LNS_v0+LNS_dv, LNS_E0+LNS_dE, LNS_dT)
  
  ! Precompute gradients, only if viscosity (classical or vibrational) exists whatsoever.
  if(LNS_viscous .or. LNS_avib) then
    ! Note: this will compute nabla_dv only if viscosity is activated OR vibrational attenuation is needed.
    call compute_gradient_TFSF(LNS_dv, LNS_dT, (LNS_viscous .or. LNS_avib), LNS_viscous, &
                               LNS_switch_gradient, nabla_dv, nabla_dT, timelocal)
    ! Precompute \Sigma_v' from \nabla v', only if viscosity exists whatsoever. This is not needed when only vibrational attenuation is needed.
    if(LNS_viscous) call LNS_compute_viscous_stress_tensor(nabla_dv, sigma_dv)
  endif ! Endif on (LNS_viscous .or. LNS_avib).

#ifdef USE_MPI
  ! Fill MPI buffers related to the variables.
  call LNS_fill_MPI_buffers_var1_bckg0(.true.)
#endif
  
  ! Compute right-hand side (RHS) of the differential system.
  call compute_forces_acoustic_LNS(LNS_drho, LNS_rho0dv, LNS_dE, LNS_e1, & ! Constitutive variables.
                                   LNS_dm, LNS_dp, nabla_dT, nabla_dv, sigma_dv, & ! Precomputed quantities. sigma_dv is sent even if viscosity is deactivated, but in that case it should be zero and unused in the subroutine.
                                   RHS_drho, RHS_rho0dv, RHS_dE, RHS_e1, & ! Output.
                                   timelocal) ! Time.
  
  ! Inverse mass matrix multiplication, in order to obtain actual RHS.
  RHS_drho               = RHS_drho * rmass_inverse_acoustic_DG ! RHS = A^{-1}*b
  RHS_dE                 = RHS_dE   * rmass_inverse_acoustic_DG
  ! Update the auxiliary register. Note: no need to zero it beforehand if LNS_scheme_A(1) is 0.
  aux_drho               = LNS_scheme_A(i_stage)*aux_drho + deltat*RHS_drho ! U_{i} = a_{i}*U_{i-1} + dt*RHS
  aux_dE                 = LNS_scheme_A(i_stage)*aux_dE   + deltat*RHS_dE
  ! Update the state register.
  LNS_drho               = LNS_drho + LNS_scheme_B(i_stage)*aux_drho ! Y^{n+1} = Y^{n+1} + b_i*U_{i}
  LNS_dE                 = LNS_dE   + LNS_scheme_B(i_stage)*aux_dE
  
  ! Group treatment for momentum in a loop on spatial dimension.
  do i_aux=1,NDIM
    RHS_rho0dv(i_aux,:) = RHS_rho0dv(i_aux,:) * rmass_inverse_acoustic_DG(:)
    aux_rho0dv(i_aux,:) = LNS_scheme_A(i_stage)*aux_rho0dv(i_aux,:) + deltat*RHS_rho0dv(i_aux,:)
    LNS_rho0dv(i_aux,:) = LNS_rho0dv(i_aux,:) + LNS_scheme_B(i_stage)*aux_rho0dv(i_aux,:)
  enddo
  
  ! Vibrational attenuation.
  if(LNS_avib) then
    RHS_e1 = RHS_e1 * rmass_inverse_acoustic_DG
    aux_e1 = LNS_scheme_A(i_stage)*aux_e1 + deltat*RHS_e1
    LNS_e1 = LNS_e1 + LNS_scheme_B(i_stage)*aux_e1
  endif
  
  ! Check linear hypothesis.
  if(check_linear_hypothesis) call LNS_warn_nonsense()
  
#if 0
  ! Test: posteriori damping in absorbing boundary condition buffers.
  if(.false.) then
    if(ABC_STRETCH) then
      if(it==1 .and. i_stage==1 .and. myrank==0) then
        write(*,*) "********************************"
        write(*,*) "*           WARNING            *"
        write(*,*) "********************************"
        write(*,*) "* A posteriori damping of the  *"
        write(*,*) "* solution activated. Solution *"
        write(*,*) "* is being damped in the       *"
        write(*,*) "* absorbing buffers.           *"
        write(*,*) "********************************"
      endif
      call damp_solution_LNS(LNS_drho, LNS_rho0dv, LNS_dE) ! See 'prepare_stretching.f90'.
    endif
  endif
#endif

end subroutine compute_forces_acoustic_LNS_main


! ------------------------------------------------------------ !
! LNS_fill_MPI_buffers_var1_bckg0                              !
! ------------------------------------------------------------ !
! Fills MPI buffers with either the background state or the constitutive variables.
! Called exclusively by the 'compute_forces_acoustic_LNS_main' routine above.

subroutine LNS_fill_MPI_buffers_var1_bckg0(var1_bckg0)
  use constants, only: NDIM
  use specfem_par, only: NPROC, ninterface_acoustic, gammaext_DG, buffer_DG_gamma_P
  use specfem_par_lns, only: NVALSIGMA, LNS_viscous, &
                             LNS_drho, LNS_rho0dv, LNS_dE, nabla_dT, sigma_dv, &
                             buffer_LNS_drho_P, buffer_LNS_rho0dv_P, buffer_LNS_dE_P, &
                             buffer_LNS_nabla_dT, buffer_LNS_sigma_dv, &
                             LNS_rho0, LNS_v0, LNS_E0, LNS_p0, LNS_kappa, sigma_v_0, &
                             buffer_LNS_rho0, buffer_LNS_E0, buffer_LNS_p0, buffer_LNS_kappa, buffer_LNS_v0, buffer_sigma_v0_P
  ! Input/Output.
  logical, intent(in) :: var1_bckg0
  ! Local.
  integer :: i_aux
  if (NPROC > 1 .and. ninterface_acoustic > 0) then
    if(var1_bckg0) then
      ! Fill MPI buffers related to the variables.
      ! Assemble state buffers.
      call assemble_MPI_vector_DG(LNS_drho, buffer_LNS_drho_P)
      do i_aux=1,NDIM
        call assemble_MPI_vector_DG(LNS_rho0dv(i_aux, :), buffer_LNS_rho0dv_P(i_aux, :, :))
      enddo
      call assemble_MPI_vector_DG(LNS_dE, buffer_LNS_dE_P)
      ! Assemble viscous buffers.
      if(LNS_viscous) then ! Check if viscosity exists whatsoever.
        do i_aux=1,NDIM
          call assemble_MPI_vector_DG(nabla_dT(i_aux, :), buffer_LNS_nabla_dT(i_aux, :, :))
        enddo
        do i_aux=1,NVALSIGMA
          call assemble_MPI_vector_DG(sigma_dv(i_aux, :), buffer_LNS_sigma_dv(i_aux, :, :))
        enddo
      endif
    else
      ! Fill the MPI buffers related to background model.
      call assemble_MPI_vector_DG(LNS_rho0, buffer_LNS_rho0)
      do i_aux=1,NDIM
        call assemble_MPI_vector_DG(LNS_v0(i_aux,:), buffer_LNS_v0(i_aux, :, :))
      enddo
      call assemble_MPI_vector_DG(LNS_E0, buffer_LNS_E0)
      call assemble_MPI_vector_DG(LNS_p0, buffer_LNS_p0)
      call assemble_MPI_vector_DG(gammaext_DG, buffer_DG_gamma_P)
      if(LNS_viscous) then
        call assemble_MPI_vector_DG(LNS_kappa, buffer_LNS_kappa)
        call assemble_MPI_vector_DG(sigma_v_0, buffer_sigma_v0_P)
      else
        ! Take this opportunity to deallocate unneeded buffers.
        deallocate(buffer_LNS_kappa)
        deallocate(buffer_sigma_v0_P)
      endif
    endif
  endif
end subroutine LNS_fill_MPI_buffers_var1_bckg0


! ------------------------------------------------------------ !
! damp_solution_LNS                                            !
! ------------------------------------------------------------ !
! Function used to damp the constitutive variables in pre-defined buffers.
! Possibly called in the 'compute_forces_acoustic_LNS_main' routine (see above), but commented out for now.

subroutine damp_solution_LNS(drho, rho0dv, dE)

  use specfem_par, only: nspec, coord, ibool_DG, nglob_DG, &
        ibool_before_perio, &
        ABC_STRETCH_TOP_LBUF, ABC_STRETCH_LEFT_LBUF, ABC_STRETCH_BOTTOM_LBUF, ABC_STRETCH_RIGHT_LBUF, &
        ispec_is_acoustic_DG, stretching_buffer, &
        mesh_xmin, mesh_xmax, mesh_zmin, mesh_zmax
  
  use constants, only: CUSTOM_REAL,NGLLX,NGLLZ,NDIM

  implicit none
  
  ! Input/output.
  real(kind=CUSTOM_REAL), dimension(nglob_DG), intent(inout) :: drho, dE
  real(kind=CUSTOM_REAL), dimension(NDIM, nglob_DG), intent(inout) :: rho0dv
  
  ! Local
  integer :: iglob
  real(kind=CUSTOM_REAL) :: sigma
  integer iglob_unique,ispec,i,j
  real(kind=CUSTOM_REAL), parameter :: ZERO = 0._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: ONE  = 1._CUSTOM_REAL
  real(kind=CUSTOM_REAL) :: x, z
  real(kind=CUSTOM_REAL) :: r_l
  
  do ispec = 1, nspec
    if(ispec_is_acoustic_DG(ispec)) then
      if(any(stretching_buffer(pack(ibool_before_perio(:,:,ispec),.true.))>0)) then
        do j = 1, NGLLZ
          do i = 1, NGLLX
            iglob_unique = ibool_before_perio(i, j, ispec)
            iglob = ibool_DG(i, j, ispec)
            x = coord(1, iglob_unique)
            z = coord(2, iglob_unique)
            sigma = ONE ! In case something bad happens.
            ! Load damping value.
            if(ibits(stretching_buffer(ibool_before_perio(i,j,ispec)),0,1)==1) then
              r_l = (z - mesh_zmax)/ABC_STRETCH_TOP_LBUF + ONE
              if(r_l>ZERO .and. r_l<=ONE) call damp_function(r_l, sigma)
            endif
            if(ibits(stretching_buffer(ibool_before_perio(i,j,ispec)),1,1)==1) then
              r_l = ONE - (x - mesh_xmin)/ABC_STRETCH_LEFT_LBUF
              if(r_l>ZERO .and. r_l<=ONE) call damp_function(r_l, sigma)
            endif
            if(ibits(stretching_buffer(ibool_before_perio(i,j,ispec)),2,1)==1) then
              r_l = ONE - (z - mesh_zmin)/ABC_STRETCH_BOTTOM_LBUF
              if(r_l>ZERO .and. r_l<=ONE) call damp_function(r_l, sigma)
            endif
            if(ibits(stretching_buffer(ibool_before_perio(i,j,ispec)),3,1)==1) then
              r_l = (x - mesh_xmax)/ABC_STRETCH_RIGHT_LBUF + ONE
              if(r_l>ZERO .and. r_l<=ONE) call damp_function(r_l, sigma)
            endif
            ! Quiet state is zero.
            ! Damp perturbation.
            drho(iglob)     = sigma*drho(iglob)
            rho0dv(:,iglob) = sigma*rho0dv(:,iglob)
            dE(iglob)       = sigma*dE(iglob)
          enddo
        enddo
      endif
    endif
  enddo
end subroutine damp_solution_LNS


! ------------------------------------------------------------ !
! initial_state_LNS                                            !
! ------------------------------------------------------------ !
! Computes the initial state for the LNS module.
! This routine is called by the 'prepare_timerun' routine ('prepare_timerun.f90'), well before the 'iterate_time' routine
! ('iterate_time.f90'), and therefore before the 'compute_forces_acoustic_LNS_main' routine (above).

subroutine initial_state_LNS()
  use constants, only: CUSTOM_REAL, NGLLX, NGLLZ, TINYVAL, HUGEVAL
  use specfem_par, only: MODEL, ibool_DG, nspec, coord, myrank, ibool_before_perio, deltat, &
                         gravityext, muext, kappa_DG, tau_epsilon, tau_sigma, &
                         NPROC, gammaext_DG, ispec_is_acoustic_DG
  use specfem_par_LNS, only: LNS_E0, LNS_p0, LNS_rho0, LNS_v0, LNS_T0, LNS_mu, &
                             buffer_LNS_sigma_dv, buffer_LNS_nabla_dT, sigma_dv, &
                             LNS_eta, LNS_kappa, LNS_g,LNS_c0, LNS_dummy_1d, LNS_dummy_2d, &
                             VALIDATION_MMS, &
                             LNS_e1, RHS_e1, aux_e1, LNS_avib, LNS_avib_taueps, LNS_avib_tausig, &
                             LNS_viscous, sigma_v_0, nabla_v0, LNS_switch_gradient

  implicit none
  
  ! Input/Output.
  ! N/A.
  
  ! Local.
  real(kind=CUSTOM_REAL), parameter :: ZEROcr = 0._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: ONEcr = 1._CUSTOM_REAL
  integer :: ispec, iglob, i, j
  
  ! Initialisation of the background state.
  if(trim(MODEL)=='LNS_generalised') then
    ! If a generalised model was loaded, do not re-initialise the 'LNS_*0' variables.
  else
    ! If no external model was loaded, or if the classical stratified DG model was loaded, specify the 'LNS_*0' variables.
    ! Initialise initial state registers.
    LNS_rho0   = ZEROcr
    LNS_v0     = ZEROcr
    LNS_E0     = ZEROcr
    ! Initialise initial auxiliary quantities.
    LNS_T0     = ZEROcr
    LNS_p0     = ZEROcr
    nabla_v0   = ZEROcr
    sigma_v_0  = ZEROcr
    ! Physical parameters.
    LNS_g      = ZEROcr
    LNS_mu     = ZEROcr
    LNS_eta    = ZEROcr
    LNS_kappa  = ZEROcr
    do ispec = 1, nspec
      if(ispec_is_acoustic_DG(ispec)) then
        do j = 1, NGLLZ
          do i = 1, NGLLX
            iglob = ibool_DG(i, j, ispec)
            call set_fluid_properties(i, j, ispec)
            call background_physical_parameters(i, j, ispec, ZEROcr, LNS_rho0(iglob), &
                                                .true., LNS_v0(:,iglob), &
                                                .true., LNS_E0(iglob), &
                                                .true., LNS_p0(iglob))
          enddo
        enddo
      endif
    enddo

    ! Recompute and save globally c0.
    LNS_c0 = ZEROcr
    where(LNS_rho0 > TINYVAL) LNS_c0 = sqrt(gammaext_DG*LNS_p0/LNS_rho0)

    ! Initialise T_0.
    call compute_T(LNS_rho0, LNS_v0, LNS_E0, LNS_T0)
  endif
  
  ! Deallocate old DG variables in all cases.
  ! Not really optimal, but less invasive.
  if(allocated(gravityext)) deallocate(gravityext)
  if(allocated(muext)) deallocate(muext)
  if(allocated(kappa_DG)) deallocate(kappa_DG)
  if(allocated(tau_epsilon)) deallocate(tau_epsilon)
  if(allocated(tau_sigma)) deallocate(tau_sigma)
  if(LNS_avib) then
    ! Make sure taueps tausig nonzero and equal in elastic elements.
    where(LNS_avib_taueps==ZEROcr) LNS_avib_taueps = HUGEVAL
    where(LNS_avib_tausig==ZEROcr) LNS_avib_tausig = HUGEVAL
    if(min(minval(LNS_avib_taueps), minval(LNS_avib_tausig)) < 0.2_CUSTOM_REAL*deltat) then
      write(*,*) "********************************"
      write(*,*) "*            ERROR             *"
      write(*,*) "********************************"
      write(*,*) "* Some relaxation time for     *"
      write(*,*) "* vibrational attenuation is   *"
      write(*,*) "* lower than 0.2*dt (empirical *"
      write(*,*) "* threshold). You cannot hope  *"
      write(*,*) "* to model a process with such *"
      write(*,*) "* a characteristic time with   *"
      write(*,*) "* such a dt. Set dt lower than *"
      write(*,*) "* ", min(minval(LNS_avib_taueps), minval(LNS_avib_tausig))/0.2_CUSTOM_REAL
      write(*,*) "* or consider relaxation times *"
      write(*,*) "* higher than                  *"
      write(*,*) "* ", 0.2*deltat
      write(*,*) "* If you're feeling extra      *"
      write(*,*) "* zealous, you may deactivate  *"
      write(*,*) "* this stop in the source      *"
      write(*,*) "* code, at                     *"
      write(*,*) "* 'compute_forces_acoustic_LNS_calling_routine.f90'."
      write(*,*) "********************************"
      call exit_MPI(myrank, '')
    else
      write(*,*) "********************************"
      write(*,*) "*           WARNING            *"
      write(*,*) "********************************"
      write(*,*) "* Activated vibrational        *"
      write(*,*) "* attenuation modelling.       *"
      write(*,*) "* EXPERIMENTAL, as it's NOT    *"
      write(*,*) "* VALIDATED. Be sure to know   *"
      write(*,*) "* what you're doing.           *"
      write(*,*) "********************************"
    endif
  else
    deallocate(LNS_avib_taueps, LNS_avib_tausig)
    deallocate(LNS_e1, RHS_e1, aux_e1)
  endif
  
  ! Detect if viscosity exists somewhere on this CPU, and set the LNS_viscous switch accordingly.
  if(     maxval(LNS_mu) > TINYVAL &
     .or. maxval(LNS_eta) > TINYVAL &
     .or. maxval(LNS_kappa) > TINYVAL) then
    LNS_viscous = .true.
  else
    LNS_viscous = .false.
    deallocate(LNS_mu, LNS_eta, LNS_kappa) ! Ambitious deallocate to free memory in the inviscid case.
#ifdef USE_MPI
    if(NPROC > 1) then
      deallocate(buffer_LNS_nabla_dT, buffer_LNS_sigma_dv) ! Even more ambitious deallocate.
    endif
#endif
    deallocate(sigma_dv) ! Even more ambitious deallocate. nabla_dv may be added here if it was declared in specfem_par_lns.
  endif
  
  ! Initialise the gradient of background velocity, \nabla v_0.
  call compute_gradient_TFSF(LNS_v0, LNS_dummy_1d, .true., .false., LNS_switch_gradient, nabla_v0, LNS_dummy_2d, ZEROcr) ! Dummy variables are not optimal, but prevent from duplicating subroutines.
  if(LNS_viscous) then ! Check if viscosity exists whatsoever.
    ! Initialise \Sigma_v_0.
    call LNS_compute_viscous_stress_tensor(nabla_v0, sigma_v_0)
  endif

  ! DEBUG. Print the model to some file.
  if(.false. .and. myrank==0) then
    open(unit=504,file='OUTPUT_FILES/TESTMODEL',status='unknown',action='write', position="append")
    do ispec = 1,nspec
      do j = 1,NGLLZ
        do i = 1,NGLLX
          iglob = ibool_DG(i, j, ispec)
          write(504,*) coord(1, ibool_before_perio(i, j, ispec)), coord(2, ibool_before_perio(i, j, ispec)),&
                       !LNS_rho0(iglob)
                       !LNS_p0(iglob)
                       gammaext_DG(iglob)
                       
        enddo
      enddo
    enddo
    close(504)
    call exit_MPI(myrank, 'STOPPING')
    ! Matlab one-liner plot:
    !a=importdata("../OUTPUT_FILES/TESTMODEL");X=a(:,1);Y=a(:,2);V=a(:,3);scatter(X,Y,20,V,'filled'); colorbar
  endif
    
  ! Eventually initialise MMS validation coefficients.
  IF(VALIDATION_MMS) call initialise_VALIDATION_MMS()
  
end subroutine initial_state_LNS


! ------------------------------------------------------------ !
! background_physical_parameters                               !
! ------------------------------------------------------------ !
! Initialises the values for the background state. May thus be used as initialiser (if time is 0), or for far-field boundary
! conditions or for bottom forcings. Both reasons rely on the specified background parameters, or initial state.
! Notes:
!   *) This model-building routine builds essentially the same model as the 'boundary_condition_DG' routine (in
!      'boundary_terms_DG.f90') does.
!   *) This routine is called for two reasons at two different places:
!      *) In 'initial_state_LNS' to initialise the simulation.
!      *) In 'LNS_get_interfaces_unknowns' (in 'compute_forces_acoustic_LNS.f90') to set the far-field outer boundary
!         conditions.

subroutine background_physical_parameters(i, j, ispec, timelocal, out_rho, swComputeV, out_v, swComputeE, out_E, swComputeP, out_p)
  use constants, only: CUSTOM_REAL, TINYVAL, NDIM
  use specfem_par, only: MODEL, assign_external_model, coord, coord_interface, gammaext_DG, ibool_DG, ibool_before_perio, &
                         pext_dg, rhoext, SCALE_HEIGHT, sound_velocity, surface_density, TYPE_FORCING, USE_ISOTHERMAL_MODEL, &
                         wind, windxext, myrank
  use specfem_par_LNS, only: LNS_g, LNS_rho0, LNS_v0, LNS_p0

  implicit none
  
  ! Input/Output.
  integer, intent(in) :: i, j, ispec
  real(kind=CUSTOM_REAL), intent(out) :: out_rho, out_E, out_p
  real(kind=CUSTOM_REAL), dimension(NDIM), intent(out) :: out_v
  real(kind=CUSTOM_REAL), intent(in) :: timelocal
  logical, intent(in) :: swComputeV, swComputeE, swComputeP
  
  ! Local.
  integer :: iglob
  real(kind=CUSTOM_REAL), parameter :: ZEROcr = 0._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: HALFcr = 0.5_CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: ONEcr = 1._CUSTOM_REAL
  real(kind=CUSTOM_REAL) :: z, H
  
  if(swComputeE .and. (.not. swComputeV)) then
    write(*,*) "********************************"
    write(*,*) "*            ERROR             *"
    write(*,*) "********************************"
    write(*,*) "* Cannot ask to compute energy *"
    write(*,*) "* without computing velocity.  *"
    write(*,*) "* See the                      *"
    write(*,*) "* background_physical_parameters"
    write(*,*) "* routine.                     *"
    write(*,*) "********************************"
    call exit_MPI(myrank, " ")
  endif
  if(swComputeE .and. (.not. swComputeP)) then
    write(*,*) "********************************"
    write(*,*) "*            ERROR             *"
    write(*,*) "********************************"
    write(*,*) "* Cannot ask to compute energy *"
    write(*,*) "* without computing pressure.  *"
    write(*,*) "* See the                      *"
    write(*,*) "* background_physical_parameters"
    write(*,*) "* routine.                     *"
    write(*,*) "********************************"
    call exit_MPI(myrank, " ")
  endif
  
  iglob = ibool_DG(i, j, ispec)
  
  if(assign_external_model) then
    ! If an external model data file is given for initial conditions, read from it.
    if(trim(MODEL)=='LNS_generalised') then
      ! If a generalised model was loaded, return the 'LNS_*0' variables as they were initialised in 'lns_load_background_model.f90'.
      ! > Set density.
      out_rho = LNS_rho0(iglob)
      ! > Set pressure.
      if(swComputeP) then
        out_p = LNS_p0(iglob)
      endif
      ! > Set wind.
      if(swComputeV) then
        out_v = LNS_v0(1:NDIM, iglob)
      endif
    else
      ! If the classical stratified DG model was loaded, use it.
    
      ! > Set density.
      out_rho = rhoext(i, j, ispec)
      
      ! > Set pressure.
      if(swComputeP) then
        out_p = pext_DG(i, j, ispec)
      endif
      
      ! > Set wind.
      if(swComputeV) then
        out_v(1) = windxext(i, j, ispec)
        out_v(NDIM) = ZEROcr ! One could set vertical wind here, too.
      endif
    endif
  else
    ! If no external model data file is given (no initial conditions were specified), build model.
    
    if(USE_ISOTHERMAL_MODEL) then
      ! > Set density.
      H = SCALE_HEIGHT ! Also for pressure, below.
      z = real(coord(NDIM, ibool_before_perio(i, j, ispec)), kind=CUSTOM_REAL) - coord_interface
      out_rho = surface_density*exp(-z/H)
      ! > Set pressure.
      if(swComputeP) then
        out_p = out_rho*LNS_g(iglob)*H
      endif
    else
      ! > Set density.
      out_rho = surface_density
      ! > Set pressure.
      if(swComputeP) then
        out_p = (sound_velocity**2)*out_rho/gammaext_DG(iglob) ! Acoustic only (under ideal gas hypothesis): p = c^2 * \rho / \gamma.
      endif
    endif
    
    ! > Set wind.
    if(swComputeV) then
      out_v(1) = wind ! Re-read horizontal wind from the scalar value read from parfile.
      out_v(NDIM) = ZEROcr ! One could set vertical wind here, too.
    endif
  endif ! Endif on assign_external_model.
  
  ! Eventually set vertical wind.
  if(swComputeV) then
    if(TYPE_FORCING/=0 .and. timelocal>ZEROcr) then
      ! If forcing exists, apply it, but only after very first iteration.
      ! This should pose no problem, since anyhow the forcing should not be significantly non-zero at the very first iteration.
      call forcing_DG(i, j, ispec, timelocal, out_v(NDIM))
    else
      ! Else, force it to be zero.
      out_v(NDIM) = ZEROcr ! One could set vertical wind here, too.
    endif
  endif
  
  ! Set energy based on pressure.
  if(swComputeE) then
    call compute_E_i(out_rho, out_v, out_p, out_E, iglob)
  endif
end subroutine background_physical_parameters


! ------------------------------------------------------------ !
! set_fluid_properties                                         !
! ------------------------------------------------------------ !
! Set the fluid's intrinsic properties.
! Notes:
!   *) This routine sets essentially the same values as the 'boundary_condition_DG' (in 'boundary_terms_DG.f90') routine does.
!   *) Called exclusively by the 'initial_state_LNS' routine above.

subroutine set_fluid_properties(i, j, ispec)
  use constants, only: CUSTOM_REAL, TINYVAL, NDIM, FOUR_THIRDS
  use specfem_par, only: assign_external_model, cp, c_v, dynamic_viscosity_cte_DG, etaext, &
                         gammaext_DG, gravityext, gravity_cte_DG, ibool_DG, kappa_DG, muext, &
                         thermal_conductivity_cte_DG, USE_ISOTHERMAL_MODEL, &
                         tau_eps_cte_DG, tau_sig_cte_DG
  use specfem_par_LNS, only: isClose, LNS_eta, LNS_kappa, LNS_g, LNS_mu, LNS_avib, LNS_avib_taueps, LNS_avib_tausig

  implicit none
  
  ! Input/Output.
  integer, intent(in) :: i, j, ispec
  
  ! Local.
  integer :: iglob
  real(kind=CUSTOM_REAL), parameter :: ZEROcr = 0._CUSTOM_REAL
  
  iglob = ibool_DG(i, j, ispec)
  
  ! All these variables are allocated in setup_mesh.f90.
  
  if(assign_external_model) then
    ! If an external model data file is given for initial conditions, read from it.
    LNS_g(iglob) = gravityext(i, j, ispec)
    ! gammaext_DG unchanged, since already defined (read from atmospheric_model.dat file).
    LNS_mu(iglob) = muext(i, j, ispec)
    LNS_eta(iglob) = etaext(i, j, ispec)
    LNS_kappa(iglob) = kappa_DG(i, j, ispec)
  else ! Else on assign_external_model.
    ! If no external model data file is given (no initial conditions were specified), build model.
    if(USE_ISOTHERMAL_MODEL) then
      ! > Isothermal case.
      LNS_g(iglob) = real(gravity_cte_DG, kind=CUSTOM_REAL)
    else
      ! > Isobaric case. Since we need to stay hydrostatic, the gravity field needs to stay 0.
      LNS_g(iglob) = ZEROcr
    endif
    gammaext_DG(iglob) = cp/c_V
    LNS_mu(iglob) = dynamic_viscosity_cte_DG
    LNS_eta(iglob) = FOUR_THIRDS*dynamic_viscosity_cte_DG
    LNS_kappa(iglob) = thermal_conductivity_cte_DG
    
    ! Vibrational attenuation.
    if(.not. isClose(tau_eps_cte_DG, tau_sig_cte_DG, 1.0e-3_CUSTOM_REAL)) then
      LNS_avib = .true.
    else
      LNS_avib = .false.
    endif
    if(LNS_avib) then
      LNS_avib_taueps(iglob) = tau_eps_cte_DG
      LNS_avib_tausig(iglob) = tau_sig_cte_DG
    endif
  endif ! Endif on assign_external_model.
end subroutine set_fluid_properties


! ------------------------------------------------------------ !
! LNS_compute_viscous_stress_tensor                            !
! ------------------------------------------------------------ !
! Computes the viscous Navier-Stokes stress tensor \Sigma_v for a given velocity v.
! IN:
!   \nabla v, gradient of v, (1, 1) being dxvx, (1, 2) dzvx, (2, 1) dxvz, and (2, 2) dzvz.
! OUT:
!   \Sigma_v(v), as a size 3 vector since the stress tensor is symmetric: 1<->(1, 1), 2<->(1, 2)&(2, 1), and 3<->(2, 2).

subroutine LNS_compute_viscous_stress_tensor(nabla_v, sigma_v)
  use constants, only: CUSTOM_REAL, NDIM
  use specfem_par, only: nglob_DG
  use specfem_par_LNS, only: LNS_eta, LNS_mu
  
  implicit none
  
  ! Input/Output.
  real(kind=CUSTOM_REAL), dimension(NDIM, NDIM, nglob_DG), intent(in) :: nabla_v
  real(kind=CUSTOM_REAL), dimension(3, nglob_DG), intent(out) :: sigma_v
  
  ! Local.
  real(kind=CUSTOM_REAL), parameter :: TWOcr = 2._CUSTOM_REAL
  real(kind=CUSTOM_REAL), dimension(nglob_DG) :: LNS_lambda
  
  LNS_lambda = LNS_eta - (2._CUSTOM_REAL/3._CUSTOM_REAL)*LNS_mu
  
  ! Explicit.
  !sigma_v(1,:) = TWOcr*LNS_mu*nabla_v(1,1,:) + LNS_lambda*(nabla_v(1,1,:)+nabla_v(2,2,:))
  !sigma_v(2,:) = LNS_mu*(nabla_v(1,2,:)+nabla_v(2,1,:))
  !sigma_v(3,:) = TWOcr*LNS_mu*nabla_v(2,2,:) + LNS_lambda*(nabla_v(1,1,:)+nabla_v(2,2,:))
  
  ! Compact.
  sigma_v(1, :) = (TWOcr*LNS_mu+LNS_lambda)*nabla_v(1, 1, :) + LNS_lambda*nabla_v(NDIM, 2, :)
  sigma_v(2, :) = LNS_mu*(nabla_v(1, 2, :)+nabla_v(NDIM, 1, :))
  sigma_v(3, :) = (TWOcr*LNS_mu+LNS_lambda)*nabla_v(NDIM, 2, :) + LNS_lambda*nabla_v(1, 1, :)
  
end subroutine LNS_compute_viscous_stress_tensor


! ------------------------------------------------------------ !
! compute_gradient_TFSF                                        !
! ------------------------------------------------------------ !
! Computes the gradient of a tensor field TF and/or of a scalar field SF.
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
   
subroutine compute_gradient_TFSF(TF, SF, swTF, swSF, swMETHOD, nabla_TF, nabla_SF, timelocal)
  use constants, only: CUSTOM_REAL, NGLLX, NGLLZ, TINYVAL, NDIM
  use specfem_par, only: gammax, gammaz, hprime_xx, hprime_zz, hprimewgll_xx, hprimewgll_zz, &
                         ibool_DG, ispec_is_acoustic_DG, jacobian, link_iface_ijispec, neighbor_dg_iface, nglob_DG, nspec, &
                         nx_iface, nz_iface, rmass_inverse_acoustic_DG, weight_iface, wxgll, wzgll, xix, xiz, &
                         ABC_STRETCH, stretching_buffer, stretching_ya, ibool_before_perio
  use specfem_par_LNS, only: LNS_dE, LNS_dp, LNS_drho, LNS_dummy_1d, LNS_dummy_2d, LNS_rho0dv, LNS_v0
  
  implicit none
  
  ! Input/Output.
  real(kind=CUSTOM_REAL), dimension(NDIM, nglob_DG) :: TF
  real(kind=CUSTOM_REAL), dimension(nglob_DG) :: SF
  logical, intent(in) :: swTF, swSF, swMETHOD
  real(kind=CUSTOM_REAL), dimension(NDIM, NDIM, nglob_DG), intent(out) :: nabla_TF
  real(kind=CUSTOM_REAL), dimension(NDIM, nglob_DG), intent(out) :: nabla_SF
  real(kind=CUSTOM_REAL), intent(in) :: timelocal
  
  ! Local.
  real(kind=CUSTOM_REAL), parameter :: ZEROcr = 0._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: HALFcr = 0.5_CUSTOM_REAL
  real(kind=CUSTOM_REAL) :: SF_P ! When swMETHOD==.true., this variable is used to store the value of the scalar field SF across the element's boundary, in order to compute the flux.
  real(kind=CUSTOM_REAL), dimension(NDIM) :: TF_P, n_out ! When swMETHOD==.true., those variables are used to store the value of the tensor field TF across the element's boundary, in order to compute the flux.
  integer :: ispec,i,j,k,iglob, iglobP
  real(kind=CUSTOM_REAL) :: halfWeight
  logical :: exact_interface_flux
  integer, dimension(3) :: neighbor
  integer :: neighbour_type, iface1, iface, iface1_neighbor, iface_neighbor, ispec_neighbor ,SPCDM
  real(kind=CUSTOM_REAL), dimension(NDIM,NDIM) :: DXiEta_L ! Derivatives in \Lambda.
  real(kind=CUSTOM_REAL), dimension(NDIM) :: Jac_WzWx_L ! Jacobian multiplied by the quadrature weights.
  real(kind=CUSTOM_REAL), dimension(NDIM) :: dSF_dXiEta ! Derivatives in \Lambda, only used if swMETHOD==.false..
  real(kind=CUSTOM_REAL), dimension(NDIM) :: dTF_dXi, dTF_dEta ! Derivatives in \Lambda.
  real(kind=CUSTOM_REAL), dimension(NGLLX, NGLLZ, NDIM) :: dSF_dXZ_ctrbXi, dSF_dXZ_ctrbEta
  real(kind=CUSTOM_REAL), dimension(NDIM, NGLLX, NGLLZ, NDIM) :: dTF_dXZ_ctrbXi, dTF_dXZ_ctrbEta
  
  if(.not. (swSF .or. swTF)) then
    write(*,*) "********************************"
    write(*,*) "*           WARNING            *"
    write(*,*) "********************************"
    write(*,*) "* compute_gradient_TFSF is     *"
    write(*,*) "* being called with both       *"
    write(*,*) "* switches to .false.. We      *"
    write(*,*) "* cannot see why such a thing  *"
    write(*,*) "* would be needed, be careful. *"
    write(*,*) "********************************"
    return
  endif
  
  if(swMETHOD) then
    ! "Desintegrate."
    n_out = ZEROcr
  endif
  
  if(swSF) then
    nabla_SF = ZEROcr
  endif
  if(swTF) then
    nabla_TF = ZEROcr
  endif

  do ispec = 1, nspec ! Loop over elements.
    if (ispec_is_acoustic_DG(ispec)) then
      ! --------------------------- !
      ! Volume terms.               !
      ! --------------------------- !
      do j = 1, NGLLZ
        do i = 1, NGLLX
          iglob     = ibool_DG(i,j,ispec)
          DXiEta_L(1,    1)    = xix(i,j,ispec) ! = \partial_x\xi from report
          DXiEta_L(1,    NDIM) = xiz(i,j,ispec) ! = \partial_z\xi from report
          DXiEta_L(NDIM, 1)    = gammax(i,j,ispec) ! = \partial_x\eta from report
          DXiEta_L(NDIM, NDIM) = gammaz(i,j,ispec) ! = \partial_z\eta from report
          
          if(ABC_STRETCH .and. stretching_buffer(ibool_before_perio(i,j,ispec))>0) then
            ! See beginning of subroutine compute_forces_acoustic_LNS for detailed explanations.
            ! If you change anything, remember to do it also in the 'compute_forces_acoustic_LNS' subroutine.
            do k = 1, NDIM
              DXiEta_L(:, k) = stretching_ya(k, ibool_before_perio(i, j, ispec))*DXiEta_L(:, k)
            enddo
          endif
          
          Jac_WzWx_L(1)    = wzgll(j)*jacobian(i,j,ispec)
          Jac_WzWx_L(NDIM) = wxgll(i)*jacobian(i,j,ispec)
          ! In Jac_WzWx_L, notice how 1 is z and 2 is x. This comes from the use of the chain rule on the bilinear form, and reorganisation. See Martire's PhD thesis manuscript.
          
          if(swMETHOD) then
            ! In that case, we want to compute \int_{\Omega} (\nabla f)\Phi as:
            ! - \int_{\Omega} T\nabla\Phi + \int_{\partial\Omega} f\Phi.
            ! The idea is to store in:
            !   nabla_SF(1,:),
            !   nabla_SF(NDIM,:),
            !   nabla_TF(1,1,:),
            !   nabla_TF(1,NDIM,:),
            !   nabla_TF(NDIM,1,:),
            !   nabla_TF(NDIM,NDIM,:)
            ! the values of the approximated integrals:
            !  $\int \partial_xT\Phi_x$,
            !  $\int \partial_zT\Phi_z$,
            !  $\int \partial_xVx\Phi_x$,
            !  etc.
            ! and "desintegrate" those values by multiplying the obtained vector by the inverse mass matrix. As explained before, the integrals are computed by using the divergence theorem and taking into account the surface terms.
            
            ! Compute inner contributions.
            if(swSF) then
              dSF_dXZ_ctrbXi(i,j,:)   = Jac_WzWx_L(1) * (DXiEta_L(1,:) * SF(iglob))
              dSF_dXZ_ctrbEta(i,j,:)  = Jac_WzWx_L(2) * (DXiEta_L(2,:) * SF(iglob))
            endif
            if(swTF) then
              do k=1, NDIM
                dTF_dXZ_ctrbXi(k,i,j,:)  = Jac_WzWx_L(1) * (DXiEta_L(1,:) * TF(k,iglob))
                dTF_dXZ_ctrbEta(k,i,j,:) = Jac_WzWx_L(2) * (DXiEta_L(2,:) * TF(k,iglob))
              enddo
            endif
          else
            ! In that case, we want to compute \int_{\Omega} (\nabla f)\Phi directly as it.
            ! The idea is to store in:
            !   nabla_SF(1,:),
            !   nabla_SF(NDIM,:),
            !   nabla_TF(1,1,:),
            !   nabla_TF(1,NDIM,:),
            !   nabla_TF(NDIM,1,:),
            !   nabla_TF(NDIM,NDIM,:)
            ! the actual values of the quantities:
            !   \partial_xT,
            !   \partial_zT,
            !   \partial_xVx,
            !   etc.
            ! which is immediate through the SEM formulation.
            
            if(swSF) then
              dSF_dXiEta    = ZEROcr
            endif
            if(swTF) then
              dTF_dXi       = ZEROcr ! Vector operation.
              dTF_dEta      = ZEROcr ! Vector operation.
            endif
            
            ! Compute derivatives in unit element \Lambda.
            ! Note: we can merge the two loops because NGLLX=NGLLZ.
            do k = 1, NGLLX
              if(swSF) then
                dSF_dXiEta(1)  = dSF_dXiEta(1)  + SF(ibool_DG(k,j,ispec)) * real(hprime_xx(i,k), kind=CUSTOM_REAL)
                dSF_dXiEta(2)  = dSF_dXiEta(2)  + SF(ibool_DG(i,k,ispec)) * real(hprime_zz(j,k), kind=CUSTOM_REAL)
              endif
              if(swTF) then
                dTF_dXi(1)     = dTF_dXi(1)     + TF(1,ibool_DG(k,j,ispec)) * real(hprime_xx(i,k), kind=CUSTOM_REAL)
                dTF_dEta(1)    = dTF_dEta(1)    + TF(1,ibool_DG(i,k,ispec)) * real(hprime_zz(j,k), kind=CUSTOM_REAL)
                dTF_dXi(NDIM)  = dTF_dXi(NDIM)  + TF(NDIM,ibool_DG(k,j,ispec)) * real(hprime_xx(i,k), kind=CUSTOM_REAL)
                dTF_dEta(NDIM) = dTF_dEta(NDIM) + TF(NDIM,ibool_DG(i,k,ispec)) * real(hprime_zz(j,k), kind=CUSTOM_REAL)
              endif
            enddo ! Enddo on k.

            ! Compute derivatives in element \Omega using the chain rule.
            if(swSF) then
              nabla_SF(:,iglob) = MATMUL(TRANSPOSE(DXiEta_L), dSF_dXiEta)
            endif
            if(swTF) then
              do k=1,NDIM ! Re-using index k for spatial dimension.
                nabla_TF(k,1,iglob)    = dTF_dXi(k) * DXiEta_L(1,1) + dTF_dEta(k) * DXiEta_L(2,1)
                nabla_TF(k,NDIM,iglob) = dTF_dXi(k) * DXiEta_L(1,2) + dTF_dEta(k) * DXiEta_L(2,2)
              enddo ! Enddo on k.
            endif ! Endif on swTF.
          endif ! Endif on swMETHOD.
        enddo ! Enddo on i.
      enddo ! Enddo on j.
      
      if(swMETHOD) then
        ! "Desintegrate".
        do j = 1, NGLLZ
          do i = 1, NGLLX
            iglob = ibool_DG(i,j,ispec)
            ! Assemble the contributions using the inner contributions.
            do k = 1, NGLLX
              if(swSF) then
                nabla_SF(1,iglob)    = nabla_SF(1,iglob) - &
                                       (dSF_dXZ_ctrbXi(k,j,1) * real(hprimewgll_xx(k,i), kind=CUSTOM_REAL) + &
                                       dSF_dXZ_ctrbEta(i,k,1) * real(hprimewgll_zz(k,j), kind=CUSTOM_REAL))
                nabla_SF(NDIM,iglob) = nabla_SF(NDIM,iglob) - &
                                       (dSF_dXZ_ctrbXi(k,j,NDIM) * real(hprimewgll_xx(k,i), kind=CUSTOM_REAL) + &
                                       dSF_dXZ_ctrbEta(i,k,NDIM) * real(hprimewgll_zz(k,j), kind=CUSTOM_REAL))
              endif
              if(swTF) then
                nabla_TF(1,1,iglob)       = nabla_TF(1,1,iglob) - &
                                            (dTF_dXZ_ctrbXi(1,k,j,1) * real(hprimewgll_xx(k,i), kind=CUSTOM_REAL) + &
                                            dTF_dXZ_ctrbEta(1,i,k,1) * real(hprimewgll_zz(k,j), kind=CUSTOM_REAL))
                nabla_TF(1,NDIM,iglob)    = nabla_TF(1,NDIM,iglob) - &
                                            (dTF_dXZ_ctrbXi(1,k,j,NDIM) * real(hprimewgll_xx(k,i), kind=CUSTOM_REAL) + &
                                            dTF_dXZ_ctrbEta(1,i,k,NDIM) * real(hprimewgll_zz(k,j), kind=CUSTOM_REAL))
                nabla_TF(NDIM,1,iglob)    = nabla_TF(NDIM,1,iglob) - &
                                            (dTF_dXZ_ctrbXi(NDIM,k,j,1) * real(hprimewgll_xx(k,i), kind=CUSTOM_REAL) + &
                                            dTF_dXZ_ctrbEta(NDIM,i,k,1) * real(hprimewgll_zz(k,j), kind=CUSTOM_REAL))
                nabla_TF(NDIM,NDIM,iglob) = nabla_TF(NDIM,NDIM,iglob) - &
                                            (dTF_dXZ_ctrbXi(NDIM,k,j,NDIM) * real(hprimewgll_xx(k,i), kind=CUSTOM_REAL) + &
                                            dTF_dXZ_ctrbEta(NDIM,i,k,NDIM) * real(hprimewgll_zz(k,j), kind=CUSTOM_REAL))
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
            iglob = ibool_DG(i,j,ispec)
            if(swSF) then
              SF_P         = ZEROcr
            endif
            if(swTF) then
              TF_P         = ZEROcr
            endif
            
            ! Step 1: prepare the normals' parameters (n_out(1), n_out(NDIM), weight, etc.).
            ! Interior point
            n_out(1)    = nx_iface(iface, ispec)
            n_out(NDIM) = nz_iface(iface, ispec)
            halfWeight = weight_iface(iface1,iface, ispec)*HALFcr
            neighbor = -1
            if(neighbor_DG_iface(iface1, iface, ispec, 3) > -1) then
              iface1_neighbor = neighbor_DG_iface(iface1, iface, ispec, 1)
              iface_neighbor  = neighbor_DG_iface(iface1, iface, ispec, 2)
              ispec_neighbor  = neighbor_DG_iface(iface1, iface, ispec, 3)
              neighbor(1:2)   = link_iface_ijispec(iface1_neighbor, iface_neighbor, ispec_neighbor,1:2)
              neighbor(3)     = ispec_neighbor
              neighbour_type = 1 ! <=> a neighbouring LNS DG element was found in the same partition
            else
              neighbour_type = 10
            endif
            
            ! Step 2: knowing the normals' parameters, compute now the fluxes.
            iglobP = 1
            if(neighbor(1) > -1) then
              iglobP = ibool_DG(neighbor(1), neighbor(2), neighbor(3))
            endif
            exact_interface_flux = .false. ! Reset this variable to .false.: by default, the fluxes have to be computed (jump!=0). In some specific cases (assigned during the call to LNS_get_interfaces_unknowns), the flux can be exact (jump==0).
            call LNS_get_interfaces_unknowns(i, j, ispec, iface1, iface, neighbor, neighbour_type, timelocal, & ! Point identifier (input).
                  LNS_drho(iglob), LNS_rho0dv(:,iglob), & ! Input constitutive variables, "M" side.
                  LNS_drho(iglobP), LNS_rho0dv(:,iglobP), LNS_dE(iglobP), & ! Input constitutive variables, "P" side.
                  LNS_dp(iglob), & ! Input other variable, "M" side.
                  n_out, & ! Normal vector (input).
                  exact_interface_flux, & ! Switch to disable jump in some cases (output).
                  LNS_dummy_1d(1), LNS_dummy_2d(:,1), LNS_dummy_1d(1), & ! Output constitutive variables.
                  LNS_dummy_2d(:,1), LNS_dummy_1d(1), TF_P, & ! Output other variables.
                  .false., LNS_dummy_2d(:,1), LNS_dummy_1d(1:3), & ! Output other variables: viscous.
                  .true., SF_P) ! Output other variables.
            ! We know that the scalar field will always be temperature, and that the tensor field will always be the velocity, so we use the already-built LNS_get_interfaces_unknowns routine to get them. The switch is here to prevent unecessary quantites to be computed during that specific call. If other fields need to be computed, one will have to use a dedicated routine.
            call check_neighbour_type(neighbour_type) ! Safeguard: crash the program if neighbour_type outside of possible values.
            
            ! For the real stretching boundary conditions , \Sigma for each constitutive variable becomes \Ya\Sigma. It is heavy to change each and every expression where \Sigma arises. Rather, we make use of what is multiplying \Sigma.
            ! Here, in the surface integrations, easiest is the normals. But we have to do it after the call to the 'LNS_get_interfaces_unknowns' routine and the affectation of lambda.
            ! If you change anything to this, remember to do it also in the 'compute_gradient_TFSF' subroutine ('compute_forces_acoustic_LNS_calling_routine.f90').
            if(ABC_STRETCH .and. stretching_buffer(ibool_before_perio(i,j,ispec))>0) then
              do SPCDM = 1, NDIM
                n_out(SPCDM) = n_out(SPCDM) * stretching_ya(SPCDM, ibool_before_perio(i,j,ispec))
              enddo
            endif
            
            if(swTF) then
              if(abs(timelocal)<TINYVAL) then
                ! At timelocal==0, we want to compute \nabla v_0, and thus we need to add v_0 to TF_P. TF already contains v_0 (in the call itself).
                ! At timelocal==0, we do not care about \nabla T.
                ! For timelocal>0, we know we are computing \nabla v' or \nabla T', and thus we do not need to add anything.
                TF_P = TF_P + LNS_v0(:,iglob)
              endif
            endif
            ! Dot products.
            do i=1,NDIM
              if(swSF) then
                nabla_SF(i,iglob) = nabla_SF(i,iglob) + halfWeight*(SF(iglob)+SF_P)*n_out(i)
              endif
              if(swTF) then
                do k=1,NDIM
                  nabla_TF(k,i,iglob) = nabla_TF(k,i,iglob) + halfWeight*(TF(k,iglob)+TF_P(k))*n_out(i)
                enddo
              endif
            enddo
          enddo ! Enddo on iface1.
        enddo ! Enddo on iface.
      endif ! Endif on swMETHOD.
    endif ! End of test if acoustic element.
  enddo ! Enddo on ispec.
  
  if(swMETHOD) then
    ! "Desintegrate".
    do i=1,NDIM
      if(swSF) then
        nabla_SF(i,:)  = nabla_SF(i,:) * rmass_inverse_acoustic_DG(:)
      endif
      !nabla_SF(NDIM,:)  = nabla_SF(NDIM,:) * rmass_inverse_acoustic_DG(:)
      if(swTF) then
        do j=1,NDIM
          nabla_TF(i,j,:) = nabla_TF(i,j,:) * rmass_inverse_acoustic_DG(:)
        enddo
      endif
    enddo
  endif
end subroutine compute_gradient_TFSF


#if 0
! This routine is unused.
! ------------------------------------------------------------ !
! compute_p                                                    !
! ------------------------------------------------------------ !
! Computes pressure from constitutive variables.

subroutine compute_p(in_rho, in_v, in_E, out_p)
  use constants, only: CUSTOM_REAL, NDIM
  use specfem_par, only: gammaext_DG, nglob_DG
  use specfem_par_LNS, only: norm2
  implicit none
  ! Input/Output.
  real(kind=CUSTOM_REAL), dimension(nglob_DG), intent(in) :: in_rho, in_E
  real(kind=CUSTOM_REAL), dimension(NDIM, nglob_DG), intent(in) :: in_v
  real(kind=CUSTOM_REAL), dimension(nglob_DG), intent(out) :: out_p
  ! Local.
  real(kind=CUSTOM_REAL), parameter :: ONEcr = 1._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: HALFcr = 0.5_CUSTOM_REAL
  out_p =   (gammaext_DG - ONEcr) &
          !* (in_E - HALFcr * in_rho * (  in_v(1,:)**2 &
          !                             + in_v(NDIM,:)**2))
          * (in_E - HALFcr*in_rho*norm2(in_v))
end subroutine compute_p
#endif


! ------------------------------------------------------------ !
! compute_dp                                                   !
! ------------------------------------------------------------ !
! Computes pressure perturbation from constitutive variables.
! Note: this could have been done nearly inline by using the subroutine compute_p, but defining this function enables one to use less RAM.

subroutine compute_dp(in_rho, in_v, in_E, out_dp)
  use constants, only: CUSTOM_REAL, NDIM
  use specfem_par, only: gammaext_DG, nglob_DG
  use specfem_par_LNS, only: LNS_p0, norm2
  implicit none
  ! Input/Output.
  real(kind=CUSTOM_REAL), dimension(nglob_DG), intent(in) :: in_rho, in_E
  real(kind=CUSTOM_REAL), dimension(NDIM, nglob_DG), intent(in) :: in_v
  real(kind=CUSTOM_REAL), dimension(nglob_DG), intent(out) :: out_dp
  ! Local.
  real(kind=CUSTOM_REAL), parameter :: ONEcr = 1._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: HALFcr = 0.5_CUSTOM_REAL
  out_dp =   (gammaext_DG - ONEcr) &
           !* (in_E - HALFcr * in_rho * (  in_v(1,:)**2 &
           !                             + in_v(NDIM,:)**2)) &
           * (in_E - HALFcr*in_rho*norm2(in_v)) &
           - LNS_p0
end subroutine compute_dp


! ------------------------------------------------------------ !
! compute_dp_i                                                 !
! ------------------------------------------------------------ !
! Same as compute_dp, but point by point (unvectorised).

subroutine compute_dp_i(in_rho, in_v, in_E, out_p, iglob)
  use constants, only: CUSTOM_REAL, NDIM
  use specfem_par, only: gammaext_DG
  use specfem_par_LNS, only: LNS_p0, norm2r1
  implicit none
  ! Input/Output.
  real(kind=CUSTOM_REAL), intent(in) :: in_rho, in_E
  real(kind=CUSTOM_REAL), dimension(NDIM), intent(in) :: in_v
  integer, intent(in) :: iglob
  real(kind=CUSTOM_REAL), intent(out) :: out_p
  ! Local.
  real(kind=CUSTOM_REAL), parameter :: ONEcr = 1._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: HALFcr = 0.5_CUSTOM_REAL
  out_p =   (gammaext_DG(iglob) - ONEcr) &
          !* (in_E - HALFcr*in_rho*(in_v(1)**2 + in_v(NDIM)**2)) &
          !* (in_E - HALFcr*in_rho*minval(norm2(reshape(in_v,(/2,1/))))) & ! Could not make our "norm2" function work both with rank 1 arrays and rank 2 arrays, so used a trick: reshape the 2-sized vector (rank 1) to a (2,1) matrix (rank 2), send it to norm2, retrieve a 1-sized vector (rank 1), use minval to send it back to a scalar value (rank 0) as needed. It proved very unefficient in terms of computation time, so we defined another dedicated "norm2" function.
          * (in_E - HALFcr*in_rho*norm2r1(in_v)) & ! Performance seems comparable to explicit formulation above.
          - LNS_p0(iglob)
end subroutine compute_dp_i


! ------------------------------------------------------------ !
! compute_T                                                    !
! ------------------------------------------------------------ !
! Computes temperature from constitutive variables.

subroutine compute_T(in_rho, in_v, in_E, out_T)
  use constants, only: CUSTOM_REAL, NDIM, TINYVAL
  use specfem_par, only: c_V, nglob_DG
  use specfem_par_LNS, only: norm2
  implicit none
  ! Input/Output.
  real(kind=CUSTOM_REAL), dimension(nglob_DG), intent(in) :: in_rho, in_E
  real(kind=CUSTOM_REAL), dimension(NDIM, nglob_DG), intent(in) :: in_v
  real(kind=CUSTOM_REAL), dimension(nglob_DG), intent(out) :: out_T
  out_T = 0._CUSTOM_REAL
  where(in_rho>TINYVAL) out_T = (in_E/in_rho - 0.5_CUSTOM_REAL*norm2(in_v))/c_V
end subroutine compute_T


#if 0
! This routine is unused.
! ------------------------------------------------------------ !
! compute_T_i                                                  !
! ------------------------------------------------------------ !
! Same as compute_T, but point by point (unvectorised).
! Unused as of 190124.

subroutine compute_T_i(in_rho, in_v, in_E, out_T)
  use constants, only: CUSTOM_REAL, NDIM
  use specfem_par, only: c_V
  use specfem_par_LNS, only: norm2r1
  implicit none
  ! Input/Output.
  real(kind=CUSTOM_REAL), intent(in) :: in_rho, in_E
  real(kind=CUSTOM_REAL), dimension(NDIM), intent(in) :: in_v
  real(kind=CUSTOM_REAL), intent(out) :: out_T
  ! Local.
  real(kind=CUSTOM_REAL), parameter :: HALFcr = 0.5_CUSTOM_REAL
  !out_T=(in_E/in_rho - HALFcr*(in_v(1)**2+in_v(NDIM)**2))/c_V
  !out_T=(in_E/in_rho - HALFcr*minval(norm2(reshape(in_v,(/2,1/)))))/c_V ! Could not make our "norm2" function work both with rank 1 arrays and rank 2 arrays, so used a trick: reshape the 2-sized vector (rank 1) to a (2,1) matrix (rank 2), send it to norm2, retrieve a 1-sized vector (rank 1), use minval to send it back to a scalar value (rank 0) as needed. It proved very unefficient in terms of computation time, so we defined another dedicated "norm2" function.
  out_T=(in_E/in_rho - HALFcr*norm2r1(in_v))/c_V ! Performance seems comparable to explicit formulation above.
end subroutine compute_T_i
#endif


! ------------------------------------------------------------ !
! compute_dT                                                   !
! ------------------------------------------------------------ !
! Computes temperature perturbation from constitutive variables.
! Note: this could have been done nearly inline by using the subroutine compute_T, but defining this function enables one to use less RAM.

subroutine compute_dT(in_rho, in_v, in_E, out_dT)
  use constants, only: CUSTOM_REAL, NDIM, TINYVAL
  use specfem_par, only: c_V, nglob_DG
  use specfem_par_LNS, only: LNS_T0, norm2
  implicit none
  ! Input/Output.
  real(kind=CUSTOM_REAL), dimension(nglob_DG), intent(in) :: in_rho, in_E
  real(kind=CUSTOM_REAL), dimension(NDIM, nglob_DG), intent(in) :: in_v
  real(kind=CUSTOM_REAL), dimension(nglob_DG), intent(out) :: out_dT
  out_dT = 0._CUSTOM_REAL
  where(in_rho>TINYVAL) out_dT =   (in_E/in_rho - 0.5_CUSTOM_REAL*norm2(in_v))/c_V &
                                 - LNS_T0
end subroutine compute_dT


! ------------------------------------------------------------ !
! compute_dT_i                                                 !
! ------------------------------------------------------------ !
! Same as compute_dT, but point by point (unvectorised).

subroutine compute_dT_i(in_rho, in_v, in_E, out_T, iglob)
  use constants, only: CUSTOM_REAL, NDIM
  use specfem_par, only: c_V
  use specfem_par_LNS, only: LNS_T0, norm2r1
  implicit none
  ! Input/Output.
  real(kind=CUSTOM_REAL), intent(in) :: in_rho, in_E
  real(kind=CUSTOM_REAL), dimension(NDIM), intent(in) :: in_v
  real(kind=CUSTOM_REAL), intent(out) :: out_T
  integer, intent(in) :: iglob
  ! Local.
  real(kind=CUSTOM_REAL), parameter :: HALFcr = 0.5_CUSTOM_REAL
  out_T =   (in_E/in_rho - HALFcr*norm2r1(in_v))/c_V &
          - LNS_T0(iglob)
end subroutine compute_dT_i


! ------------------------------------------------------------ !
! compute_E                                                    !
! ------------------------------------------------------------ !
! Computes energy from constitutive variables.

subroutine compute_E(in_rho, in_v, in_p, out_E)
  use constants, only: CUSTOM_REAL, NDIM
  use specfem_par, only: gammaext_DG, nglob_DG
  use specfem_par_LNS, only: norm2
  implicit none
  ! Input/Output.
  real(kind=CUSTOM_REAL), dimension(nglob_DG), intent(in) :: in_rho, in_p
  real(kind=CUSTOM_REAL), dimension(NDIM, nglob_DG), intent(in) :: in_v
  real(kind=CUSTOM_REAL), dimension(nglob_DG), intent(out) :: out_E
  ! Local.
  real(kind=CUSTOM_REAL), parameter :: ONEcr = 1._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: HALFcr = 0.5_CUSTOM_REAL
  out_E =   in_p/(gammaext_DG - ONEcr) &
          + in_rho*HALFcr*norm2(in_v)
end subroutine compute_E


! ------------------------------------------------------------ !
! compute_E_i                                                  !
! ------------------------------------------------------------ !
! Same as compute_E, but point by point (unvectorised).

subroutine compute_E_i(in_rho, in_v, in_p, out_E, iglob)
  use constants, only: CUSTOM_REAL, NDIM
  use specfem_par, only: gammaext_DG
  use specfem_par_LNS, only: norm2r1
  implicit none
  ! Input/Output.
  real(kind=CUSTOM_REAL), intent(in) :: in_rho, in_p
  real(kind=CUSTOM_REAL), dimension(NDIM), intent(in) :: in_v
  integer, intent(in) :: iglob
  real(kind=CUSTOM_REAL), intent(out) :: out_E
  ! Local.
  real(kind=CUSTOM_REAL), parameter :: ONEcr = 1._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: HALFcr = 0.5_CUSTOM_REAL
  out_E =   in_p/(gammaext_DG(iglob) - ONEcr) &
          + in_rho*HALFcr*norm2r1(in_v)
end subroutine compute_E_i


! ------------------------------------------------------------ !
! compute_dE_i                                                 !
! ------------------------------------------------------------ !
! Computes energy perturbation from constitutive variables, but point by point (unvectorised).

subroutine compute_dE_i(in_rho, in_v, in_p, out_E, iglob)
  use constants, only: CUSTOM_REAL, NDIM
  use specfem_par, only: gammaext_DG
  use specfem_par_LNS, only: LNS_E0, norm2r1
  implicit none
  ! Input/Output.
  real(kind=CUSTOM_REAL), intent(in) :: in_rho, in_p
  real(kind=CUSTOM_REAL), dimension(NDIM), intent(in) :: in_v
  integer, intent(in) :: iglob
  real(kind=CUSTOM_REAL), intent(out) :: out_E
  ! Local.
  real(kind=CUSTOM_REAL), parameter :: ONEcr = 1._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: HALFcr = 0.5_CUSTOM_REAL
  out_E =   in_p/(gammaext_DG(iglob) - ONEcr) &
          + in_rho*HALFcr*norm2r1(in_v) &
          - LNS_E0(iglob)
end subroutine compute_dE_i


! ------------------------------------------------------------ !
! stressBuilder_addInviscidFluid                               !
! ------------------------------------------------------------ !
! Routine which builds on top of an already-existing stress tensor, adding the inviscid part of the fluid stress.

subroutine stressBuilder_addInviscidFluid(rho, v1, v2, pressure, out_sigma)
  use constants, only: CUSTOM_REAL, NDIM
  
  implicit none
  
  ! Input/Output.
  real(kind=CUSTOM_REAL), intent(in) :: rho, pressure
  real(kind=CUSTOM_REAL), dimension(NDIM), intent(in) :: v1, v2
  real(kind=CUSTOM_REAL), dimension(NDIM, NDIM), intent(inout) :: out_sigma
  
  ! Local.
  integer :: i, j
  
  do i=1, NDIM
    do j=1, NDIM
      out_sigma(i, j) = out_sigma(i, j) + rho*v1(i)*v2(j)
      if(i==j) then
        out_sigma(i, j) = out_sigma(i, j) + pressure
      endif
    enddo
  enddo
end subroutine stressBuilder_addInviscidFluid


! ------------------------------------------------------------ !
! stressBuilder_addViscousFluid                                !
! ------------------------------------------------------------ !
! Routine which builds on top of an already-existing stress tensor, adding the viscous part of the fluid stress.

subroutine stressBuilder_addViscousFluid(viscous_tensor_local, out_sigma)
  use constants, only: CUSTOM_REAL, NDIM
  use specfem_par_lns, only: NVALSIGMA
  
  implicit none
  
  ! Input/Output.
  real(kind=CUSTOM_REAL), dimension(NVALSIGMA), intent(in) :: viscous_tensor_local
  real(kind=CUSTOM_REAL), dimension(NDIM, NDIM), intent(inout) :: out_sigma
  
  ! Local.
  ! N./A.
  
  out_sigma(1, 1) = out_sigma(1, 1) + viscous_tensor_local(1)
  out_sigma(1, 2) = out_sigma(1, 2) + viscous_tensor_local(2)
  out_sigma(2, 1) = out_sigma(2, 1) + viscous_tensor_local(2)
  out_sigma(2, 2) = out_sigma(2, 2) + viscous_tensor_local(3)
end subroutine stressBuilder_addViscousFluid


! ------------------------------------------------------------ !
! LNS_prevent_nonsense                                         !
! ------------------------------------------------------------ !
! Attempts to detect nonsense in LNS variables, and stop program if they are found.

subroutine LNS_prevent_nonsense()
  use constants, only: TINYVAL
  use specfem_par, only: myrank
  use specfem_par_LNS, only: LNS_rho0, LNS_E0, LNS_p0, LNS_mu, LNS_eta, LNS_kappa, LNS_viscous
  implicit none
  ! Input/Output.
  ! N./A.
  ! Local.
  ! N./A.
  
  ! Initial state.
  if(minval(LNS_rho0)<TINYVAL) then
    call exit_MPI(myrank, "LNS_rho0 is non-positive (<= 0) somewhere.")
  endif
  if(minval(LNS_E0)<TINYVAL) then
    call exit_MPI(myrank, "LNS_E0 is non-positive (<= 0) somewhere.")
  endif
  
  ! Auxiliary initial variables.
  if(minval(LNS_p0)<TINYVAL) then
    call exit_MPI(myrank, "LNS_p0 is non-positive (<= 0) somewhere.")
  endif
  
  ! Physical parameters.
  if(LNS_viscous) then ! Check if viscosity exists whatsoever.
    if(minval(LNS_mu)<-TINYVAL) then
      call exit_MPI(myrank, "LNS_mu is negative (< 0) somewhere.")
    endif
    if(minval(LNS_eta)<-TINYVAL) then
      call exit_MPI(myrank, "LNS_eta is negative (< 0) somewhere.")
    endif
    if(minval(LNS_kappa)<-TINYVAL) then
      call exit_MPI(myrank, "LNS_kappa is negative (< 0) somewhere.")
    endif
  endif
end subroutine LNS_prevent_nonsense


! ------------------------------------------------------------ !
! LNS_warn_nonsense                                            !
! ------------------------------------------------------------ !
! Attempts to detect nonsense in LNS variables, and warn user if they are found.

subroutine LNS_warn_nonsense()
  use constants, only: CUSTOM_REAL, NDIM, NGLLX, NGLLZ, TINYVAL
  use specfem_par, only: nspec, coord, myrank, ispec_is_acoustic_DG, ibool_DG, ibool_before_perio
  use specfem_par_LNS, only: LNS_drho, LNS_dv, LNS_dE, LNS_rho0, LNS_v0, LNS_E0
  implicit none
  ! Input/Output.
  ! N./A.
  ! Local.
  real(kind=CUSTOM_REAL), parameter :: thresholdRatioUp = real(1.0d4, kind=CUSTOM_REAL)
  integer :: i,j,ispec,iglob,d,broken
  
  broken = 0
  
  outer:do ispec=1,nspec
    if(ispec_is_acoustic_DG(ispec)) then
      do j=1,NGLLZ
        do i=1,NGLLX
          iglob=ibool_DG(i,j,ispec)
          if((LNS_drho(iglob)/LNS_rho0(iglob))>thresholdRatioUp) then
            broken = 1
            !write(*,*) i,j,ispec,ispec_is_acoustic_DG(ispec) ! DEBUG
            exit outer
          endif
          do d=1,NDIM
            if((LNS_v0(d,iglob))>TINYVAL) then
              if((LNS_dv(d,iglob)/LNS_v0(d,iglob))>thresholdRatioUp) then
                broken = 2
                exit outer
              endif
            endif
          enddo
          if((LNS_dE(iglob)/LNS_E0(iglob))>thresholdRatioUp) then
            broken = 3
            exit outer
          endif
        enddo
      enddo
    endif
  enddo outer
  
  if(broken/=0) then
    write(*,*) "********************************"
    write(*,*) "*           WARNING            *"
    write(*,*) "********************************"
    write(*,*) "* A ratio is high somewhere.   *"
    write(*,*) "* Linear hypothesis might      *"
    write(*,*) "* break.                       *"
    write(*,*) "********************************"
    write(*,*) "* X         | ", coord(:, ibool_before_perio(i, j, ispec))
    select case (broken)
      case (1)
        write(*,*) "* ratio     | LNS_drho/LNS_rho0     *"
        write(*,*) "* max(drho) | ", maxval(LNS_drho)
        write(*,*) "* min(rho0) | ", minval(LNS_rho0)
      case (2)
        write(*,*) "* ratio     | LNS_dv/LNS_v0         *"
        write(*,*) "* d         | ", d
        write(*,*) "* max(dv)   | ", maxval(LNS_dv)
        write(*,*) "* min(v0)   | ", minval(LNS_v0)
      case (3)
        write(*,*) "* ratio     | LNS_dE/LNS_E0         *"
        write(*,*) "* max(dE)   | ", maxval(LNS_dE)
        write(*,*) "* min(E0)   | ", minval(LNS_E0)
    end select
    write(*,*) "* CPU       | ", myrank
    write(*,*) "* threshold | ", thresholdRatioUp
    write(*,*) "********************************"
    write(*,*) "* If run crashes, we advise    *"
    write(*,*) "* the user to try using FNS.   *"
    write(*,*) "********************************"
  endif
end subroutine LNS_warn_nonsense


! ------------------------------------------------------------ !
! check_neighbour_type                                         !
! ------------------------------------------------------------ !
! Checks the value taken by the variable neighbour_type, used in the 'compute_forces_acoustic_LNS' routine
! ('compute_forces_acoustic_LNS.f90').
! Note: the acceptable values for neighbour_type are:
!   neighbour_type = 1  <=> A neighbouring LNS DG element was found in the same partition.
!   neighbour_type = 10 <=> Temporary value waiting for precision.
!   neighbour_type = 11 <=> A neighbouring LNS DG element was found in another partition.
!   neighbour_type = 21 <=> Other material: viscoelastic.
!   neighbour_type = 22 <=> Other material: potential acoustic (not implemented).
!   neighbour_type = 99 <=> Outside boundary.

subroutine check_neighbour_type(neighbour_type)
  use specfem_par, only: myrank
  integer, intent(in) :: neighbour_type
  if(.not.(neighbour_type==1 .or. &
           neighbour_type==11 .or. &
           neighbour_type==21 .or. &
           neighbour_type==99)) then
    write(*,*) 'neighbour_type = ', neighbour_type
    call exit_MPI(myrank,'neighbour_type should not have this value.')
  endif
end subroutine check_neighbour_type


! ------------------------------------------------------------ !
! initialise_VALIDATION_MMS                                    !
! ------------------------------------------------------------ !
! Initialises the constants used for the validation of the fluid equations using the method of manufactured solutions.

subroutine initialise_VALIDATION_MMS()
  use constants, only: CUSTOM_REAL
  use specfem_par, only: myrank
  use specfem_par_LNS, only: USE_LNS, LNS_viscous, &
                             VALIDATION_MMS, VALIDATION_MMS_IV, VALIDATION_MMS_KA, VALIDATION_MMS_MU, &
                             MMS_dRHO_cst, MMS_dVX_cst, MMS_dVZ_cst, MMS_dE_cst, &
                             MMS_dRHO_x, MMS_dRHO_z, &
                             MMS_dVX_x, MMS_dVX_z, MMS_dVZ_x, MMS_dVZ_z, &
                             MMS_dE_x, MMS_dE_z
  
  implicit none
  
  ! Input/Output.
  ! N.A.
  
  ! Local.
  ! N.A.
  
  if(.not. USE_LNS) call exit_MPI(myrank, "CANNOT TEST MMS WITHOUT LNS")
  if(.not. VALIDATION_MMS) call exit_MPI(myrank, "THIS ROUTINE SHOULD NOT BE CALLED IF NOT VALIDATION_MMS")
  if(VALIDATION_MMS_IV) then
    if(LNS_viscous) call exit_MPI(myrank, "CANNOT TEST MMS INVISCID IF VISCOSITY")
  endif
  if(VALIDATION_MMS_KA) then
    if(.not. LNS_viscous) call exit_MPI(myrank, "CANNOT TEST MMS KAPPA IF NO VISCOSITY")
  endif
  if(VALIDATION_MMS_MU) then
    if(.not. LNS_viscous) call exit_MPI(myrank, "CANNOT TEST MMS MU IF NO VISCOSITY")
  endif
  if(VALIDATION_MMS_IV) then
    MMS_dRHO_cst = 0.001_CUSTOM_REAL
    MMS_dVX_cst  = 0._CUSTOM_REAL
    MMS_dE_cst   = 0.05_CUSTOM_REAL
  endif
  if(VALIDATION_MMS_KA) then
    MMS_dRHO_cst = 0._CUSTOM_REAL
    MMS_dVX_cst  = 0._CUSTOM_REAL
    MMS_dE_cst   = 0.05_CUSTOM_REAL
  endif
  if(VALIDATION_MMS_MU) then
    MMS_dRHO_cst = 0._CUSTOM_REAL
    MMS_dVX_cst  = 0.001_CUSTOM_REAL
    MMS_dE_cst   = 0._CUSTOM_REAL
  endif
  MMS_dVZ_cst = 0._CUSTOM_REAL
  MMS_dRHO_x  = 1._CUSTOM_REAL
  MMS_dRHO_z  = 2._CUSTOM_REAL
  MMS_dVX_x   = 2._CUSTOM_REAL
  MMS_dVX_z   = 0._CUSTOM_REAL
  MMS_dVZ_x   = 0._CUSTOM_REAL
  MMS_dVZ_z   = 0._CUSTOM_REAL
  MMS_dE_x    = 3._CUSTOM_REAL
  MMS_dE_z    = 4._CUSTOM_REAL
end subroutine initialise_VALIDATION_MMS


! ------------------------------------------------------------ !
! VALIDATION_MMS_source_terms                                  !
! ------------------------------------------------------------ !
! Produces the source terms for the validation of the fluid equations using the method of manufactured solutions.

subroutine VALIDATION_MMS_source_terms(outrhs_drho, outrhs_rho0dv, outrhs_dE)
  use constants, only: CUSTOM_REAL, NGLLX, NGLLZ, NDIM, TINYVAL, PI
  use specfem_par, only: jacobian, wxgll, wzgll, &
                         nspec, nglob_DG, ibool_DG, ibool_before_perio, &
                         gammaext_dg, myrank, &
                         coord, ibool_before_perio, c_V, sound_velocity
  use specfem_par_LNS, only: USE_LNS, NVALSIGMA, &
                             LNS_rho0, LNS_v0, &
                             LNS_kappa, LNS_mu, &
                             VALIDATION_MMS, VALIDATION_MMS_IV, VALIDATION_MMS_KA, VALIDATION_MMS_MU, &
                             MMS_dRHO_cst, MMS_dVX_cst, MMS_dE_cst, &!MMS_dVZ_cst, &
                             MMS_dRHO_x, MMS_dRHO_z, &
                             MMS_dVX_x, &!MMS_dVX_z, &
                             !MMS_dVZ_x, MMS_dVZ_z, &
                             MMS_dE_x, MMS_dE_z, &
                             LNS_verbose, LNS_modprint
  
  implicit none
  
  ! Input/Output.
  real(kind=CUSTOM_REAL), dimension(nglob_DG), intent(out) :: outrhs_drho, outrhs_dE
  real(kind=CUSTOM_REAL), dimension(NDIM, nglob_DG), intent(out) :: outrhs_rho0dv
  
  ! Local.
  real(kind=CUSTOM_REAL), parameter :: EIGHT_THIRDScr  = (8._CUSTOM_REAL/3._CUSTOM_REAL)
  real(kind=CUSTOM_REAL), parameter :: ONEcr  = 1._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: HALFcr  = 0.5_CUSTOM_REAL
  integer :: ispec, i, j, iglob
  real(kind=CUSTOM_REAL), dimension(NDIM) :: X
  real(kind=CUSTOM_REAL) :: GAM, w
  if(.not. USE_LNS) call exit_MPI(myrank, "CANNOT TEST MMS WITHOUT LNS")
  if(.not. VALIDATION_MMS) call exit_MPI(myrank, "THIS ROUTINE SHOULD NOT BE CALLED IF NOT VALIDATION_MMS")
  do ispec = 1, nspec
    do j = 1, NGLLZ
      do i = 1, NGLLX
        iglob = ibool_DG(i, j, ispec)
        X = coord(:, ibool_before_perio(i, j, ispec))
        GAM = gammaext_DG(iglob)
        w = LNS_v0(1, iglob)
        ! Mass conservation.
        if(VALIDATION_MMS_IV) then
            outrhs_drho(iglob) = MMS_dRHO_cst*MMS_dRHO_x*PI*cos(MMS_dRHO_x*PI*X(1))*w
        endif
        if(VALIDATION_MMS_KA) then
            outrhs_drho(iglob) = MMS_dRHO_cst*MMS_dRHO_x*PI*cos(MMS_dRHO_x*PI*X(1))*w
        endif
        if(VALIDATION_MMS_MU) then
            outrhs_drho(iglob) = LNS_rho0(iglob)*MMS_dVX_cst*MMS_dVX_x*PI*cos(MMS_dVX_x*PI*X(1))
        endif
        outrhs_drho(iglob) = outrhs_drho(iglob) * wxgll(i)*wzgll(j)*jacobian(i,j,ispec)
        ! x-momentum.
        if(VALIDATION_MMS_IV) then
          outrhs_rho0dv(1, iglob) = (GAM-ONEcr)*PI &
                                         *(  MMS_dE_cst*MMS_dE_x*cos(MMS_dE_x*PI*X(1)) &
                                           - HALFcr*MMS_dRHO_cst*MMS_dRHO_x*cos(MMS_dRHO_x*PI*X(1))*w**2 )
        endif
        if(VALIDATION_MMS_KA) then
          outrhs_rho0dv(1, iglob) = (GAM-ONEcr)*PI &
                                         *(  MMS_dE_cst*MMS_dE_x*cos(MMS_dE_x*PI*X(1)) &
                                           - HALFcr*MMS_dRHO_cst*MMS_dRHO_x*cos(MMS_dRHO_x*PI*X(1))*w**2 )
        endif
        if(VALIDATION_MMS_MU) then
          outrhs_rho0dv(1, iglob) = ( &
                                         LNS_rho0(iglob)*MMS_dVX_cst*MMS_dVX_x*PI*cos(MMS_dVX_x*PI*X(1)) &
                                         *(w-(GAM-ONEcr)*(w+MMS_dVX_cst*sin(MMS_dVX_x*PI*X(1)))) &
                                       + &
                                         PI**2*MMS_dVX_cst*MMS_dVX_x**2*EIGHT_THIRDScr*LNS_mu(iglob)*sin(MMS_dVX_x*PI*X(1)) &
                                       )
        endif
        outrhs_rho0dv(1, iglob) = outrhs_rho0dv(1, iglob) * wxgll(i)*wzgll(j)*jacobian(i,j,ispec)
        ! z-momentum.
        if(VALIDATION_MMS_IV) then
          outrhs_rho0dv(2, iglob) = (GAM-ONEcr)*PI &
                                         *(  MMS_dE_cst*MMS_dE_z*cos(MMS_dE_z*PI*X(2)) &
                                           - HALFcr*MMS_dRHO_cst*MMS_dRHO_z*cos(MMS_dRHO_z*PI*X(2))*w**2 )
        endif
        if(VALIDATION_MMS_KA) then
          outrhs_rho0dv(2, iglob) = (GAM-ONEcr)*PI &
                                         *(  MMS_dE_cst*MMS_dE_z*cos(MMS_dE_z*PI*X(2)) &
                                           - HALFcr*MMS_dRHO_cst*MMS_dRHO_z*cos(MMS_dRHO_z*PI*X(2))*w**2 )
        endif
        if(VALIDATION_MMS_MU) then
          outrhs_rho0dv(2, iglob) = 0._CUSTOM_REAL
        endif
        outrhs_rho0dv(2, iglob) = outrhs_rho0dv(2, iglob) * wxgll(i)*wzgll(j)*jacobian(i,j,ispec)
        ! Energy.
        if(VALIDATION_MMS_IV) then
          outrhs_dE(iglob) = w*PI &
                                  *(  GAM*MMS_dE_cst*MMS_dE_x*cos(MMS_dE_x*PI*X(1)) &
                                    - HALFcr*(GAM-ONEcr)*MMS_dRHO_cst*MMS_dRHO_x*cos(MMS_dRHO_x*PI*X(1))*w**2)
        endif
        if(VALIDATION_MMS_KA) then
          outrhs_dE(iglob) = (   w*PI*GAM*MMS_dE_cst*MMS_dE_x*cos(MMS_dE_x*PI*X(1)) &
                                  + ((MMS_dE_cst*LNS_kappa(iglob)*PI**2)/(LNS_rho0(iglob)*c_V)) &
                                    *(sin(MMS_dE_x*PI*X(1))*MMS_dE_x**2 + sin(MMS_dE_z*PI*X(2))*MMS_dE_z**2) )
        endif
        if(VALIDATION_MMS_MU) then
          outrhs_dE(iglob) = ( &
                                  (LNS_rho0(iglob)*MMS_dVX_cst*MMS_dVX_x*PI*cos(MMS_dVX_x*PI*X(1))/(GAM-ONEcr)) &
                                  *(   HALFcr*(GAM-ONEcr)*w**2 &
                                     + sound_velocity**2 &
                                     - (GAM-ONEcr)**2*w*(w+MMS_dVX_cst*sin(MMS_dVX_x*PI*X(1))) &
                                   ) &
                                + &
                                  PI**2*MMS_dVX_cst*MMS_dVX_x**2 &
                                  *(   EIGHT_THIRDScr*LNS_mu(iglob)*w*sin(MMS_dVX_x*PI*X(1)) &
                                     + (LNS_kappa(iglob)/c_V) &
                                       *( &
                                            MMS_dVX_cst*sin(2._CUSTOM_REAL*MMS_dVX_x*PI*X(1)) &
                                          + w*sin(MMS_dVX_x*PI*X(1)) &
                                        ) &
                                  ) &
                                )
        endif
        outrhs_dE(iglob) = outrhs_dE(iglob) * wxgll(i)*wzgll(j)*jacobian(i,j,ispec)
      enddo
    enddo
  enddo
end subroutine VALIDATION_MMS_source_terms


! ------------------------------------------------------------ !
! VALIDATION_MMS_boundary_terms                                !
! ------------------------------------------------------------ !
! Produces the boundary terms for the validation of the fluid equations using the method of manufactured solutions.

subroutine VALIDATION_MMS_boundary_terms(iglob, iglobM, &
                                         exact_interface_flux, out_drho_P, out_dv_P, out_dE_P, &
                                         out_dp_P, out_rho0dv_P, &
                                         swCompVisc, out_nabla_dT_P, out_sigma_dv_P, &
                                         swCompdT, out_dT_P)
  use constants, only: CUSTOM_REAL, NDIM, PI
  use specfem_par, only: coord, myrank
  use specfem_par_LNS, only: USE_LNS, NVALSIGMA, &
                             LNS_rho0, LNS_v0, LNS_E0, nabla_dT, sigma_dv, &
                             VALIDATION_MMS, &
                             MMS_dRHO_cst, MMS_dVX_cst, MMS_dVZ_cst, MMS_dE_cst, &
                             MMS_dRHO_x, MMS_dRHO_z, &
                             MMS_dVX_x, MMS_dVX_z, MMS_dVZ_x, MMS_dVZ_z, &
                             MMS_dE_x, MMS_dE_z
  
  implicit none
  
  ! Input/Output.
  integer, intent(in) :: iglob, iglobM
  logical, intent(in) :: swCompVisc, swCompdT!, swCompv ! Do not unnecessarily compute some quantities.
  logical, intent(out) :: exact_interface_flux ! Output switch.
  real(kind=CUSTOM_REAL), intent(out) :: out_drho_P, out_dE_P ! Output constitutive variables.
  real(kind=CUSTOM_REAL), dimension(NDIM), intent(out) :: out_rho0dv_P ! Output constitutive variable.
  real(kind=CUSTOM_REAL), dimension(NDIM), intent(out) :: out_dv_P, out_nabla_dT_P ! Output other variables.
  real(kind=CUSTOM_REAL), intent(out) :: out_dp_P ! Output other variables.
  real(kind=CUSTOM_REAL), dimension(NVALSIGMA), intent(out) :: out_sigma_dv_P ! Output other variables.
  real(kind=CUSTOM_REAL), intent(out) :: out_dT_P ! In compute_gradient_TFSF (desintegration method), we need temperature on the other side of the boundary in order to compute the flux.
  
  ! Local.
  real(kind=CUSTOM_REAL), parameter :: ZEROcr = 0._CUSTOM_REAL
  if(.not. USE_LNS) call exit_MPI(myrank, "CANNOT TEST MMS WITHOUT LNS")
  if(.not. VALIDATION_MMS) call exit_MPI(myrank, "THIS ROUTINE SHOULD NOT BE CALLED IF NOT VALIDATION_MMS")
  exact_interface_flux = .false. ! Do not force jump to zero.     
  out_drho_P     = MMS_dRHO_cst*(  sin(MMS_dRHO_x*PI*coord(1, iglob)) &
                                 + sin(MMS_dRHO_z*PI*coord(2, iglob)))
  out_dv_P(1)    = MMS_dVX_cst *(  sin(MMS_dVX_x*PI*coord(1, iglob)) &
                                 + sin(MMS_dVX_z*PI*coord(2, iglob)))
  out_dv_P(NDIM) = MMS_dVZ_cst *(  sin(MMS_dVZ_x*PI*coord(1, iglob)) &
                                 + sin(MMS_dVZ_z*PI*coord(2, iglob)))
  out_dE_P       = MMS_dE_cst  *(  sin(MMS_dE_x*PI*coord(1, iglob)) &
                                 + sin(MMS_dE_z*PI*coord(2, iglob)))
  call compute_dp_i(LNS_rho0(iglobM)+out_drho_P, LNS_v0(:,iglobM)+out_dv_P, LNS_E0(iglobM)+out_dE_P, out_dp_P, iglobM)
  out_rho0dv_P = LNS_rho0(iglobM)*out_dv_P
  if(swCompVisc) then
    out_nabla_dT_P = nabla_dT(:,iglobM)
    out_sigma_dv_P = sigma_dv(:,iglobM)
  else
   out_nabla_dT_P = ZEROcr
   out_sigma_dv_P = ZEROcr
  endif
  if(swCompdT) then
    call compute_dT_i(LNS_rho0(iglobM)+out_drho_P, LNS_v0(:, iglobM)+out_dv_P, LNS_E0(iglobM)+out_dE_P, out_dT_P, iglobM)
  else
    out_dT_P = ZEROcr
  endif
end subroutine VALIDATION_MMS_boundary_terms


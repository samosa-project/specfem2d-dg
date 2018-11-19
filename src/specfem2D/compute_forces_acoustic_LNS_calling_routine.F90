! ------------------------------------------------------------ !
! compute_forces_acoustic_LNS_main                              !
! ------------------------------------------------------------ !
! TODO: Description.

subroutine compute_forces_acoustic_LNS_main()
  ! TODO: select variables to use.
  use constants
  use specfem_par
  use specfem_par_LNS
  !use constants, only: &
  !  CUSTOM_REAL, &
  !  rk4a_d, rk4b_d, rk4c_d, &
  !  ls33rk_a, ls33rk_b, ls33rk_c, &
  !  HALF, ONE, NGLLX, NGLLZ
  !use specfem_par, only:&! kmato,&
  !  deltat, nglob_DG,stage_time_scheme,&
  !  !gravityext, muext, etaext, kappa_DG,&
  !  !tau_epsilon, tau_sigma, &
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
  !  NDIM,&
  !  LNS_dummy_1d, LNS_dummy_2d!, LNS_dummy_3d

  implicit none

  ! Local.
  real(kind=CUSTOM_REAL), parameter :: ZEROcr = 0._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: ONEcr = 1._CUSTOM_REAL
  real(kind=CUSTOM_REAL) :: timelocal
  real(kind=CUSTOM_REAL), dimension(stage_time_scheme) :: scheme_A, scheme_B, scheme_C
  !real(kind=CUSTOM_REAL), dimension(NDIM,nglob_DG) :: LNS_dv
  real(kind=CUSTOM_REAL), dimension(NDIM,NDIM,nglob_DG) :: nabla_dv
  !real(kind=CUSTOM_REAL), dimension(3,nglob_DG) :: sigma_dv
  !real(kind=CUSTOM_REAL), dimension(NDIM,nglob_DG) :: nabla_dT
  integer :: ier, i_aux,i,j,ispec,iglob,ispec_PML
  logical check_linearHypothesis!, check_linearHypothesis_ON_ALL_PROCS, check_linearHypothesis_FIND_POINT
  
  ! PMLs.
  real(kind=CUSTOM_REAL), dimension(NDIM) :: pml_alpha
  
  ! Checks if anything has to be done.
  if (.not. any_acoustic_DG) then
    return
  endif
  
  ! Debug switches. They are computationaly very heavy and should not be used on every simulation.
  check_linearHypothesis=.true.
  
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
    
    ! Initialise state registers. Note: since constitutive variables are perturbations, they are necessarily zero at start.
    LNS_drho   = ZEROcr
    LNS_rho0dv = ZEROcr
    LNS_dE     = ZEROcr
    ! Initialise state auxiliary quantities.
    !LNS_T0     = ZEROcr
    !LNS_p0     = ZEROcr
    !nabla_v0   = ZEROcr
    !sigma_v_0  = ZEROcr
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
    
    if(PML_BOUNDARY_CONDITIONS) then
      !stop "PML WITH LNS ARE NOT FULLY IMPLEMENTED YET."
      aux_PML_drho = ZEROcr
      aux_PML_rho0dv = ZEROcr
      aux_PML_dE = ZEROcr
      RHS_PML_drho = ZEROcr
      RHS_PML_rho0dv = ZEROcr
      RHS_PML_dE = ZEROcr
      LNS_PML_drho = ZEROcr
      LNS_PML_rho0dv = ZEROcr
      LNS_PML_dE = ZEROcr
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
      stop "Error allocating 'ispec_is_acoustic_coupling_ac' arrays (see 'compute_forces_acoustic_LNS_calling_routine.F90')."
    endif
    ispec_is_acoustic_coupling_ac = -1
    if(.not. only_DG_acoustic) then
      call find_DG_acoustic_coupling()
    endif
  
    call LNS_prevent_nonsense() ! Check initial conditions.
  endif ! Endif on (it == 1) and (i_stage == 1).
  
  if(LNS_VERBOSE>=1 .and. myrank == 0 .AND. mod(it, LNS_MODPRINT)==0) then
    WRITE(*,*) "****************************************************************"
    WRITE(*,"(a,i9,a,i1,a,e23.16,a)") "Iteration", it, ", stage ", i_stage, ", local time", timelocal, "."
  endif
  
  ! Precompute momentum perturbation and velocity perturbation.
  LNS_dv=ZEROcr
  do i_aux=1,NDIM
    LNS_dm(i_aux,:) = LNS_rho0dv(i_aux,:)+LNS_drho*LNS_v0(i_aux,:)
    where(LNS_rho0/=ZEROcr) LNS_dv(i_aux,:)=LNS_rho0dv(i_aux,:)/LNS_rho0 ! 'where(...)' as safeguard, as rho0=0 should not happen.
  enddo
  
  ! Recompute temperature and pressure perturbation.
  call compute_dT(LNS_rho0+LNS_drho, LNS_v0+LNS_dv, LNS_E0+LNS_dE, LNS_dT)
  call compute_dp(LNS_rho0+LNS_drho, LNS_v0+LNS_dv, LNS_E0+LNS_dE, LNS_dp)
  
  !write(*,*)it, i_stage, "rho0", LNS_rho0(1), "E0", LNS_E0(1), "T0", LNS_T0(1), "p0", LNS_p0(1), &
  !             "gamma", gammaext_DG(1), "c_V", c_V, &
  !             "dp", LNS_dp(1)! DEBUG
  !stop
  
  ! Precompute gradients.
  if(LNS_viscous) then ! Check if viscosity exists whatsoever.
    call compute_gradient_TFSF(LNS_dv, LNS_dT, LNS_viscous, LNS_viscous, LNS_switch_gradient, nabla_dv, nabla_dT, timelocal) ! Note: we compute nabla_dv only if viscosity is activated.
    ! Precompute \Sigma_v'.
    !write(*,*) "computing sigma_dv"
    call LNS_compute_viscous_stress_tensor(nabla_dv, sigma_dv)
    !write(*,*) "computing sigma_dv", maxval(sigma_dv), maxval(nabla_dv)
  endif

#ifdef USE_MPI
  ! If MPI, assemble at each iteration.
  if (NPROC > 1 .and. ninterface_acoustic > 0) then
    ! Assemble state buffers.
    call assemble_MPI_vector_DG(LNS_drho, buffer_LNS_drho_P)
    do i_aux=1,NDIM
      call assemble_MPI_vector_DG(LNS_rho0dv(i_aux,:), buffer_LNS_rho0dv_P(i_aux,:,:))
    enddo
    call assemble_MPI_vector_DG(LNS_dE, buffer_LNS_dE_P)
    
    ! Assemble viscous buffers.
    if(LNS_viscous) then ! Check if viscosity exists whatsoever.
      do i_aux=1,NDIM
        call assemble_MPI_vector_DG(nabla_dT(i_aux, :), buffer_LNS_nabla_dT(i_aux,:,:))
      enddo
      do i_aux=1,NVALSIGMA
        call assemble_MPI_vector_DG(sigma_dv(i_aux, :), buffer_LNS_sigma_dv(i_aux,:,:))
      enddo
    endif
  endif
#endif
  
  !do ispec=1,nspec
  !do i=1,NGLLX
  !do j=1,NGLLZ
  !if(      abs(coord(1,ibool_before_perio(i,j,ispec))-10.)<2. &
  !   .and. abs(coord(2,ibool_before_perio(i,j,ispec))-10.)<2.) then
  !  write(*,*) it, i_stage, coord(:,ibool_before_perio(i,j,ispec)), nabla_dv(:,:,ibool_DG(i,j,ispec))
  !endif
  !enddo
  !enddo
  !enddo
  !if(it==6 .and. i_stage==1) stop 'kekest'
  
  ! Compute RHS.
  ! Note: if there is PML BCs, additionnal terms will be queried directly inside the call to compute_forces_acoustic_LNS, since auxiliary variables and coefficients are global variables.
  call compute_forces_acoustic_LNS(LNS_drho, LNS_rho0dv, LNS_dE, & ! Constitutive variables.
                                   LNS_dm, LNS_dp, nabla_dT, sigma_dv, & ! Precomputed quantities. sigma_dv is sent even if viscosity is deactivated, but in that case it should be zero and unused in the subroutine.
                                   RHS_drho, RHS_rho0dv, RHS_dE, & ! Output.
                                   timelocal) ! Time.
  
  if (time_stepping_scheme == 3 .or. time_stepping_scheme == 4) then
    ! Inverse mass matrix multiplication, in order to obtain actual RHS.
    RHS_drho(:)            = RHS_drho(:)*rmass_inverse_acoustic_DG(:) ! RHS = A^{-1}*b
    RHS_dE(:)              = RHS_dE(:)  *rmass_inverse_acoustic_DG(:)
    ! Update the auxiliary register. Note: no need to zero it beforehand if scheme_A(1) is 0.
    aux_drho               = scheme_A(i_stage)*aux_drho + deltat*RHS_drho ! U_{i} = a_{i}*U_{i-1} + dt*RHS
    aux_dE                 = scheme_A(i_stage)*aux_dE   + deltat*RHS_dE
    ! Update the state register.
    LNS_drho               = LNS_drho + scheme_B(i_stage)*aux_drho ! Y^{n+1} = Y^{n+1} + b_i*U_{i}
    LNS_dE                 = LNS_dE   + scheme_B(i_stage)*aux_dE
    ! Group momentum treatment in loop on dimension.
    do i_aux=1,NDIM
      RHS_rho0dv(i_aux,:) = RHS_rho0dv(i_aux,:)*rmass_inverse_acoustic_DG(:)
      aux_rho0dv(i_aux,:) = scheme_A(i_stage)*aux_rho0dv(i_aux,:) + deltat*RHS_rho0dv(i_aux,:)
      LNS_rho0dv(i_aux,:) = LNS_rho0dv(i_aux,:) + scheme_B(i_stage)*aux_rho0dv(i_aux,:)
    enddo
    ! If one wants to do all at once (method which should be working but seems to be somewhat more unstable), comment the previous lines and uncomment the following lines.
    !LNS_drho               =   LNS_drho &
    !                         + scheme_B(i_stage)*(  scheme_A(i_stage)*aux_drho &
    !                                              + deltat*(RHS_drho(:)*rmass_inverse_acoustic_DG(:))) ! Y^{n+1} = Y^{n+1} + b_i*U_{i}
    !LNS_dE                 =   LNS_dE &
    !                         + scheme_B(i_stage)*(  scheme_A(i_stage)*aux_dE &
    !                                              + deltat*(RHS_dE(:)*rmass_inverse_acoustic_DG(:)))
    !do i_aux=1,NDIM
    !  LNS_rho0dv(i_aux,:) =   LNS_rho0dv(i_aux,:) &
    !                        + scheme_B(i_stage)*(  scheme_A(i_stage)*aux_rho0dv(i_aux,:) &
    !                                             + deltat*(RHS_rho0dv(i_aux,:)*rmass_inverse_acoustic_DG(:)))
    !enddo
    
    ! PML, iterate ADEs.
    if(PML_BOUNDARY_CONDITIONS) then
      !do i_aux=1,NDIM ! Loop on ADEs.
      
        ! Prepare PML ADE RHS, vectorially. TODO: a dedicated routine, or include it in some other loop (typically, using the compute_forces_acoustic_LNS call would be a good idea).
        do ispec=1,nspec; if(ispec_is_PML(ispec)) then; ispec_PML=spec_to_PML(ispec); do j=1,NGLLZ; do i=1,NGLLX
          iglob=ibool_DG(i, j, ispec)
          !Write(*,*) "iglob", iglob, "is_acoustic_dg", ispec_is_acoustic_DG(ispec)
          !WRITE(*,*) "LNS_PML_b", LNS_PML_b(i_aux,i,j,ispec_PML)
          !WRITE(*,*) "LNS_PML_drho", LNS_PML_drho(i_aux,i,j,ispec_PML)
          !WRITE(*,*) "LNS_drho", LNS_drho(iglob)
          pml_alpha(1)    = alpha_x_store(i,j,ispec_PML) + 0.001_CUSTOM_REAL ! TODO: make adding the t0 is systematic, or check if this can work without
          pml_alpha(NDIM) = alpha_z_store(i,j,ispec_PML) + 0.002_CUSTOM_REAL ! TODO: make adding the t0 is systematic, or check if this can work without
          RHS_PML_drho(:,i,j,ispec_PML) =   pml_alpha*LNS_PML_drho(:,i,j,ispec_PML) &
                                          - LNS_drho(iglob) ! A minus sign might have to be used here.
          RHS_PML_dE(:,i,j,ispec_PML) =   pml_alpha*LNS_PML_dE(:,i,j,ispec_PML) &
                                        - LNS_dE(iglob) ! A minus sign might have to be used here.
          do i_aux=1,NDIM ! Loop on momenta.
            RHS_PML_rho0dv(:,i_aux,i,j,ispec_PML) =   pml_alpha*LNS_PML_rho0dv(:,i_aux,i,j,ispec_PML) &
                                                    - LNS_rho0dv(i_aux, iglob) ! A minus sign might have to be used here.
          enddo
        enddo; enddo; endif; enddo
        
        ! Update PML ADE, vectorially.
        aux_PML_drho(:,:,:,:) = scheme_A(i_stage)*aux_PML_drho(:,:,:,:) + deltat*RHS_PML_drho(:,:,:,:)
        aux_PML_dE(:,:,:,:)   = scheme_A(i_stage)*aux_PML_dE(:,:,:,:)   + deltat*RHS_PML_dE(:,:,:,:)
        LNS_PML_drho(:,:,:,:) = LNS_PML_drho(:,:,:,:) + scheme_B(i_stage)*aux_PML_drho(:,:,:,:)
        LNS_PML_dE(:,:,:,:)   = LNS_PML_dE(:,:,:,:)   + scheme_B(i_stage)*aux_PML_dE(:,:,:,:)
        aux_PML_rho0dv(:,:,:,:,:) =   scheme_A(i_stage)*aux_PML_rho0dv(:,:,:,:,:) &
                                    + deltat*RHS_PML_rho0dv(:,:,:,:,:)
        LNS_PML_rho0dv(:,:,:,:,:) =   LNS_PML_rho0dv(:,:,:,:,:) &
                                    + scheme_B(i_stage)*aux_PML_rho0dv(:,:,:,:,:)
        !do i_aux=1,NDIM ! Loop on momenta.
        !  aux_PML_rho0dv(:,i_aux,:,:,:) =   scheme_A(i_stage)*aux_PML_rho0dv(:,i_aux,:,:,:) &
        !                                  + deltat*RHS_PML_rho0dv(:,i_aux,:,:,:)
        !  LNS_PML_rho0dv(:,i_aux,:,:,:) =   LNS_PML_rho0dv(:,i_aux,:,:,:) &
        !                                  + scheme_B(i_stage)*aux_PML_rho0dv(:,i_aux,:,:,:)
        !enddo
      !enddo
    endif
    
    ! test
    !if(PML_BOUNDARY_CONDITIONS) then
    !  if(myrank==1) then
    !    write(*,*) timelocal, minval(aux_PML_drho), maxval(aux_PML_drho)
    !  endif
    !endif
    !write(*,*) maxval(abs(LNS_dT)), minval(abs(LNS_dT)), maxval(abs(nabla_dT)), minval(abs(nabla_dT))
    !stop
    !write(*,*) "soudnspeed lol", minval(sqrt(gammaext_DG*LNS_p0/LNS_rho0)), maxval(sqrt(gammaext_DG*LNS_p0/LNS_rho0))
    !write(*,*) allocated(gravityext), allocated(muext), allocated(etaext), &
    !           allocated(kappa_DG), allocated(tau_epsilon), allocated(tau_sigma)
    !kek = reshape(LNS_rho0dv(:,5),(/2,1/))
    !write(*,*) 'kek', kek
    !write(*,*) 'v1', norm2(LNS_rho0dv(:,5:6))
    !write(*,*) 'v2', norm2(reshape(LNS_rho0dv(:,5),(/2,1/)))
    !write(*,*) '|v1-v2|', abs(norm2(LNS_rho0dv(:,5:8))-(LNS_rho0dv(1,5:8)**2 + LNS_rho0dv(NDIM,5:8)**2))
    !write(*,*) "kek"
    !write(*,*) size(ibool_before_perio) ! 0 to nspec
    !write(*,*) size(ibool) ! 0 to nspec
    !write(*,*) size(ibool_DG) ! 0 to nglob_dg
    !write(*,*) size(coord(2,:))
    !write(*,*) ibool
    !write(*,*) size(pack(ibool(:,:,:), ibool(:,:,:)<4))
    !write(*,*) size(LNS_rho0dv(2,:))
    !write(*,*) pack(LNS_rho0dv(2,:),abs(coord(2,:))<1.)
    !stop
    
    ! Eventually check non-positivity.
    !LNS_dv(1,:)=1e4*LNS_dv(1,:)
    !LNS_v0=0.*LNS_v0+1d-5
    
    if(check_linearHypothesis) then
      call LNS_warn_nonsense()
      !if(     check_linearHypothesis_ON_ALL_PROCS &
      !   .or. ((.not. check_linearHypothesis_ON_ALL_PROCS) .and. myrank==0) &
      !  ) then
      !  ! If:    we check on all procs,
      !  !     or we check only on proc 0 and we are on proc 0.
      !  if(minval(LNS_drho) < 10d-14) then
      !    WRITE(*,*) "***************************************************************"
      !    WRITE(*,*) "* CAREFUL, VERY SMALL DENSITY: ", minval(LNS_drho), "    *"
      !    if(check_linearHypothesis_FIND_POINT) then
      !      ! Find where density is low.
      !      do ispec = 1,nspec
      !        do j = 1,NGLLZ
      !          do i = 1,NGLLX
      !            !write(*,*) LNS_drho(ibool_DG(ispec, i, j))
      !            !write(*,*) ispec, i, j, ibool_DG(ispec, i, j)
      !            if(LNS_drho(ibool_DG(i, j, ispec))==minval(LNS_drho)) then
      !              WRITE(*, *) "* Element", ispec, ", GLL", i, j, ".         *"
      !              write(*, *) "* Coords", coord(1, ibool_before_perio(i, j, ispec)), &
      !                          coord(NDIM, ibool_before_perio(i, j, ispec)), ".*"
      !              !call virtual_stretch_prime(i, j, ispec, coef_stretch_x_ij_prime, coef_stretch_z_ij_prime)
      !              !write(*, *) coef_stretch_x_ij_prime, coef_stretch_z_ij_prime
      !            endif
      !          enddo
      !        enddo
      !      enddo
      !    endif ! Endif on check_linearHypothesis_FIND_POINT.
      !    WRITE(*,*) "***************************************************************"
      !  endif ! Endif on minval(LNS_drho).
      !endif ! Endif on check_linearHypothesis_ON_ALL_PROCS.
    endif ! Endif on check_linearHypothesis.
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
  !    veloc_x = LNS_rho0dv(NDIM,:) - LNS_v0(NDIM,:)
  !    call SlopeLimit1(veloc_x, timelocal, 1)
  !    LNS_rho0dv(NDIM,:) = veloc_x + LNS_v0(NDIM,:)
  !  else
  !    call SlopeLimit1(LNS_rho0dv(NDIM,:), timelocal, 3)
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
  use constants, only: CUSTOM_REAL, NGLLX, NGLLZ, TINYVAL
  use specfem_par, only: ibool_DG, nspec, gravityext,muext, etaext, kappa_DG, tau_epsilon, &
                         tau_sigma, NPROC, buffer_DG_gamma_P, gammaext_DG, ninterface_acoustic!, &
                         !PML_BOUNDARY_CONDITIONS !, ispec_is_acoustic_DG, nspec
  use specfem_par_LNS, only: LNS_E0, LNS_p0, LNS_rho0, LNS_v0, LNS_T0, LNS_mu, &
                             LNS_eta, LNS_kappa, sigma_dv, LNS_dummy_1d, LNS_dummy_2d, &
                             LNS_viscous, sigma_v_0, nabla_v0, buffer_LNS_nabla_dT, &
                             buffer_LNS_sigma_dv, LNS_switch_gradient, LNS_g

  implicit none
  
  ! Input/Output.
  ! N/A.
  
  ! Local.
  real(kind=CUSTOM_REAL), parameter :: ZEROcr = 0._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: ONEcr = 1._CUSTOM_REAL
  integer :: ispec, iglob, i, j
  !real(kind=CUSTOM_REAL) :: rho_P, veloc_x_P, veloc_z_P, E_P
  
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
    !if(ispec_is_acoustic_DG(ispec)) then
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
    !endif
  enddo
  
  deallocate(gravityext, muext, etaext, kappa_DG, tau_epsilon, tau_sigma) ! Not really optimal, but less invasive.
  
  ! Ugly patch to make it work.
  ! It seems nglob_DG includes EVERY POINT on the mesh, including those in elastic materials. See nglob_DG definition in 'createnum_slow.f90'.
  ! Thus, since DG-related variables are badly initialised in other materials (since they make little to no sense in those), arithmetic errors can arise.
  ! However, we leave it as is and use the following patches, in fear of breaking the FNS implementation.
  ! TODO: something else, maybe involving correcting the FNS implementation.
  where(LNS_rho0 < TINYVAL) LNS_rho0 = ONEcr
  where(LNS_p0 < TINYVAL) LNS_p0 = ONEcr
  where(LNS_E0 < TINYVAL) LNS_E0 = ONEcr
  
  ! Initialise T_0.
  call compute_T(LNS_rho0, LNS_v0, LNS_E0, LNS_T0)
  !LNS_T0=LNS_p0/(LNS_rho0*287.05684504212125397) ! Approximation assuming air molar mass is constant at 0.028964.
  
  ! Detect if viscosity exists somewhere on this CPU.
  if(     maxval(LNS_mu) > TINYVAL &
     .OR. maxval(LNS_eta) > TINYVAL &
     .OR. maxval(LNS_kappa) > TINYVAL) then
    LNS_viscous=.true.
    !write(*,*) "LNS: min(mu,eta,kappa)>0, computation will be viscous."
  else
    LNS_viscous=.false.
    !write(*,*) "LNS: max(mu,eta,kappa)=0, computation will be inviscid."
    deallocate(LNS_mu, LNS_eta, LNS_kappa) ! Ambitious deallocate to free memory in the inviscid case.
#ifdef USE_MPI
    if(NPROC > 1) then
      deallocate(buffer_LNS_nabla_dT, buffer_LNS_sigma_dv) ! Even more ambitious deallocate.
    endif
#endif
    deallocate(sigma_dv) ! Even more ambitious deallocate. nabla_dv may be added here if it was declared in specfem_par_lns.
  endif
  
  ! Initialise \nabla v_0.
  call compute_gradient_TFSF(LNS_v0, LNS_dummy_1d, .true., .false., LNS_switch_gradient, nabla_v0, LNS_dummy_2d, ZEROcr) ! Dummy variables are not optimal, but prevent from duplicating subroutines.
  if(LNS_viscous) then ! Check if viscosity exists whatsoever.
    ! Initialise \Sigma_v_0.
    call LNS_compute_viscous_stress_tensor(nabla_v0, sigma_v_0)
  endif
  
#ifdef USE_MPI
  if (NPROC > 1 .and. ninterface_acoustic > 0) then
    call assemble_MPI_vector_DG(gammaext_DG, buffer_DG_gamma_P)
  endif
#endif
end subroutine initial_state_LNS

! ------------------------------------------------------------ !
! LNS_PML_init_coefs                                           !
! ------------------------------------------------------------ !
! Initialises coefficients needed for PML implementation.
! Note: base coefficients ("K_*_store", "d_*_store", and "alpha_*_store" are initialised in pml_init.F90).

subroutine LNS_PML_init_coefs()
  ! TODO: select variables to use.
  use constants
  use specfem_par
  use specfem_par_LNS
  
  implicit none
  
  ! Input/Output.
  ! N/A.

  ! Local.
  real(kind=CUSTOM_REAL), parameter :: ZEROcr = 0._CUSTOM_REAL
  integer :: i,j,ispec,ispec_PML
  real(kind=CUSTOM_REAL), dimension(NDIM) :: pmlk, pmld, pmla
  
  ! Safeguards.
  if(     (.not. allocated(LNS_PML_a0)) &
     .or. (.not. allocated(LNS_PML_b))) then
    write(*,*) "********************************"
    write(*,*) "*            ERROR             *"
    write(*,*) "********************************"
    write(*,*) "* Some PML coefficients arrays *"
    write(*,*) "* are not allocated but should *"
    write(*,*) "* be.                          *"
    write(*,*) "********************************"
    stop
  endif
  if(.not. PML_BOUNDARY_CONDITIONS) then
    write(*,*) "********************************"
    write(*,*) "*            ERROR             *"
    write(*,*) "********************************"
    write(*,*) "* PML are not activated and    *"
    write(*,*) "* you are trying to initialise *"
    write(*,*) "* the parameters.              *"
    write(*,*) "********************************"
    stop
  endif
  
  ! These arrays were allocated in prepare_timerun_pml.f90.
  
  ! Initialise.
  LNS_PML_a0 = ZEROcr
  LNS_PML_b  = ZEROcr
  
  ! Value.
  do ispec=1,nspec
    if(ispec_is_PML(ispec)) then
      ispec_PML=spec_to_PML(ispec)
      do j=1,NGLLZ
        do i=1,NGLLX
          pmlk(1)=K_x_store(i,j,ispec_PML) ! Decrease performance, but increases readability. Since this routine only runs once, we decide it's okay.
          pmlk(2)=K_z_store(i,j,ispec_PML)
          pmld(1)=d_x_store(i,j,ispec_PML)
          pmld(2)=d_z_store(i,j,ispec_PML)
          pmla(1)=alpha_x_store(i,j,ispec_PML) + 0.001_CUSTOM_REAL ! TODO: make adding the t0 is systematic, or check if this can work without
          pmla(2)=alpha_z_store(i,j,ispec_PML) + 0.002_CUSTOM_REAL ! TODO: make adding the t0 is systematic, or check if this can work without
          
          
          LNS_PML_a0(i,j,ispec_PML) = pmlk(1)*pmld(2) + pmlk(2)*pmld(1)
          
          if(abs(pmla(2)-pmla(1)) < TINYVAL) then
            write(*,*) "|a2-a1| is very smol at ", coord(:, ibool_before_perio(i, j, ispec)), ": a1,a2=", pmla ! DEBUG
          else
            LNS_PML_b(1,i,j,ispec_PML) = - pmla(1)*pmld(1) * (   pmlk(2) &
                                                               - pmld(2)/(pmla(1)-pmla(2)) )
            LNS_PML_b(2,i,j,ispec_PML) = - pmla(2)*pmld(2) * (   pmlk(1) &
                                                               + pmld(1)/(pmla(1)-pmla(2)) )
          endif
          
          if(abs(coord(1, ibool_before_perio(i, j, ispec)))<TINYVAL) & ! Monitor slice at x==0.
            write(*,*) coord(2, ibool_before_perio(i, j, ispec)), ":", pmlk, pmld, pmla, &
                       LNS_PML_a0(i,j,ispec_PML), LNS_PML_b(:,i,j,ispec_PML) ! DEBUG
        enddo
      enddo
    endif
  enddo
end subroutine LNS_PML_init_coefs

! ------------------------------------------------------------ !
! background_physical_parameters                               !
! ------------------------------------------------------------ !
! Affects values of background state. May thus be used as initialiser (if time is 0), for far-field boundary conditions, or for bottom forcings.

subroutine background_physical_parameters(i, j, ispec, timelocal, out_rho, swComputeV, out_v, swComputeE, out_E, swComputeP, out_p)
  use constants, only: CUSTOM_REAL, TINYVAL, NDIM
  use specfem_par, only: assign_external_model, gammaext_DG, ibool_DG, pext_dg, rhoext, &
SCALE_HEIGHT, sound_velocity, surface_density, TYPE_FORCING, USE_ISOTHERMAL_MODEL, wind, windxext
  use specfem_par_LNS, only: LNS_g

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
  real(kind=CUSTOM_REAL) :: z, H!, G!, A
  
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
    stop
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
    stop
  endif
  
  iglob = ibool_DG(i, j, ispec)
  
  !write(*,*) "assign_external_model", assign_external_model
  
  if(assign_external_model) then
    ! If an external model data file is given for initial conditions, read from it.
    ! > If initialisation (condition on timelocal), set gravity, viscosity coefficients (mu, eta, and kappa).
    !write(*,*) "abs(timelocal)", abs(timelocal)
    !if(abs(timelocal)<TINYVAL) then
    !  LNS_g(iglob) = gravityext(i, j, ispec)
    !  LNS_mu(iglob) = muext(i, j, ispec)
    !  LNS_eta(iglob) = etaext(i, j, ispec)
    !  LNS_kappa(iglob) = kappa_DG(i, j, ispec)
    !endif
    
    ! > Set density.
    out_rho = rhoext(i, j, ispec)
    
    ! > Set pressure.
    if(swComputeP) then
      out_p = pext_DG(i, j, ispec)
    endif
    
    ! > Set wind.
    if(swComputeV) then
      out_v(1) = windxext(i, j, ispec)
      ! One might want to set vertical wind here, too.
    endif
  else
    ! If no external model data file is given (no initial conditions were specified), build model.
    
    ! > If initialisation (condition on timelocal), set gravity, viscosity coefficients (mu, eta, and kappa), and gamma.
    !write(*,*) "kek" ! DEBUG
    !if(abs(timelocal)<TINYVAL) then
    !  ! We are at t=0.
    !  if(USE_ISOTHERMAL_MODEL) then
    !    ! > Isothermal case.
    !    LNS_g(iglob) = real(gravity_cte_DG, kind=CUSTOM_REAL)
    !  else
    !    ! > Isobaric case. Since we need to stay hydrostatic, the gravity field needs to stay 0.
    !    LNS_g(iglob) = ZEROcr
    !  endif
    !  gammaext_DG(iglob) = cp/c_V
    !  LNS_mu(iglob) = dynamic_viscosity_cte_DG
    !  LNS_eta(iglob) = (4./3.)*dynamic_viscosity_cte_DG
    !  LNS_kappa(iglob) = thermal_conductivity_cte_DG
    !endif
    
    if(USE_ISOTHERMAL_MODEL) then
      ! > Set density.
      H = SCALE_HEIGHT ! Also for pressure, below.
      out_rho = surface_density*exp(-z/H)
      ! > Set pressure.
      if(swComputeP) then
        out_p = out_rho*LNS_g(iglob)*H
      endif
    else
      ! > Set density.
      out_rho = surface_density
      ! > Set pressure.
      !write(*,*) sound_velocity ! DEBUG
      !write(*,*) out_rho ! DEBUG
      !write(*,*) gammaext_DG(ibool_DG(i, j, ispec)) ! DEBUG
      if(swComputeP) then
        out_p = (sound_velocity**2)*out_rho/gammaext_DG(iglob) ! Acoustic only (under ideal gas hypothesis): p = c^2 * \rho / \gamma.
      endif
    endif
    
    ! > Set wind.
    if(swComputeV) then
      out_v(1) = wind ! Re-read horizontal wind from the scalar value read from parfile.
      ! One might want to set vertical wind here, too.
    endif
  endif ! Endif on assign_external_model.
  
  ! Eventually set vertical wind.
  if(swComputeV) then
    if(TYPE_FORCING/=0) then
      ! If forcing exists, apply it.
      call forcing_DG(i, j, ispec, timelocal, out_v(NDIM))
    else
      ! Else, force it to be zero.
      out_v(NDIM) = ZEROcr
      ! One might want to set vertical wind here, too.
    endif
  endif
  
  !! > Set gravity potentials.
  !if(timelocal == ZEROcr) then
  !  potential_dphi_dx_DG(ibool(i, j, ispec)) = ZEROcr
  !  potential_dphi_dz_DG(ibool(i, j, ispec)) = gravityext(i, j, ispec)
  !endif
  
  ! Set energy based on pressure.
  if(swComputeE) then
    call compute_E_i(out_rho, out_v, out_p, out_E, iglob)
    !out_E =   out_p/(gammaext_DG(iglob) - ONEcr) &
    !        + out_rho*HALFcr*( out_v(1)**2 + out_v(NDIM)**2 )
  endif
end subroutine background_physical_parameters

! ------------------------------------------------------------ !
! set_fluid_properties                                         !
! ------------------------------------------------------------ !
! Set fluid properties.

subroutine set_fluid_properties(i, j, ispec)
  use constants, only: CUSTOM_REAL, TINYVAL, NDIM
  use specfem_par, only: assign_external_model, cp, c_v, dynamic_viscosity_cte_DG, etaext, &
gammaext_DG, gravityext, gravity_cte_DG, ibool_DG, kappa_DG, muext, thermal_conductivity_cte_DG, USE_ISOTHERMAL_MODEL
  use specfem_par_LNS, only: LNS_eta, LNS_kappa, LNS_g, LNS_mu

  implicit none
  
  ! Input/Output.
  integer, intent(in) :: i, j, ispec
  
  ! Local.
  integer :: iglob
  real(kind=CUSTOM_REAL), parameter :: ZEROcr = 0._CUSTOM_REAL
  
  iglob = ibool_DG(i, j, ispec)
  
  if(assign_external_model) then
    ! If an external model data file is given for initial conditions, read from it.
    LNS_g(iglob) = gravityext(i, j, ispec)
    ! gammaext_DG unchanged, since already defined.
    LNS_mu(iglob) = muext(i, j, ispec)
    LNS_eta(iglob) = etaext(i, j, ispec)
    LNS_kappa(iglob) = kappa_DG(i, j, ispec)
    !LNS_v0(1, iglob) = windxext(i, j, ispec)
    ! One might want to initialise vertical wind here, too.
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
    LNS_eta(iglob) = (4._CUSTOM_REAL/3._CUSTOM_REAL)*dynamic_viscosity_cte_DG
    LNS_kappa(iglob) = thermal_conductivity_cte_DG
    !LNS_v0(1, iglob) = wind ! Read horizontal wind from the scalar value read from parfile.
    ! One might want to initialise vertical wind here, too.
  endif ! Endif on assign_external_model.
end subroutine set_fluid_properties

















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
  
  LNS_lambda = LNS_eta-(2._CUSTOM_REAL/3._CUSTOM_REAL)*LNS_mu
  
  ! Explicit.
  !sigma_v(1,:) = TWOcr*LNS_mu*nabla_v(1,1,:) + LNS_lambda*(nabla_v(1,1,:)+nabla_v(2,2,:))
  !sigma_v(2,:) = LNS_mu*(nabla_v(1,2,:)+nabla_v(2,1,:))
  !sigma_v(3,:) = TWOcr*LNS_mu*nabla_v(2,2,:) + LNS_lambda*(nabla_v(1,1,:)+nabla_v(2,2,:))
  
  ! Compact.
  sigma_v(1,:) = (TWOcr*LNS_mu+LNS_lambda)*nabla_v(1,1,:) + LNS_lambda*nabla_v(NDIM,2,:)
  sigma_v(2,:) = LNS_mu*(nabla_v(1,2,:)+nabla_v(NDIM,1,:))
  sigma_v(3,:) = (TWOcr*LNS_mu+LNS_lambda)*nabla_v(NDIM,2,:) + LNS_lambda*nabla_v(1,1,:)
  
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
   
subroutine compute_gradient_TFSF(TF, SF, swTF, swSF, swMETHOD, nabla_TF, nabla_SF, timelocal)
  use constants, only: CUSTOM_REAL, NGLLX, NGLLZ, TINYVAL, NDIM
  use specfem_par, only: gammax, gammaz, hprime_xx, hprime_zz, hprimewgll_xx, hprimewgll_zz,&
ibool_DG, ispec_is_acoustic_DG, jacobian, link_iface_ijispec, neighbor_dg_iface, nglob_DG, nspec, nx_iface,&
nz_iface, rmass_inverse_acoustic_DG, weight_iface, wxgll, wzgll, xix, xiz
  use specfem_par_LNS, only: LNS_dE, LNS_drho, LNS_dummy_1d, LNS_dummy_2d, LNS_rho0dv, LNS_v0
  
  implicit none
  
  ! Input/Output.
  real(kind=CUSTOM_REAL), dimension(nglob_DG) :: SF
  real(kind=CUSTOM_REAL), dimension(NDIM, nglob_DG) :: TF
  logical, intent(in) :: swTF, swSF, swMETHOD
  real(kind=CUSTOM_REAL), dimension(NDIM, nglob_DG), intent(out) :: nabla_SF
  real(kind=CUSTOM_REAL), dimension(NDIM, NDIM, nglob_DG), intent(out) :: nabla_TF
  real(kind=CUSTOM_REAL), intent(in) :: timelocal
  
  ! Local.
  real(kind=CUSTOM_REAL), parameter :: ZEROcr = 0._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: HALFcr = 0.5_CUSTOM_REAL
  real(kind=CUSTOM_REAL) :: SF_P ! When swMETHOD==.true., this variable is used to store the value of the scalar field SF across the element's boundary, in order to compute the flux.
  real(kind=CUSTOM_REAL), dimension(NDIM) :: TF_P, n_out ! When swMETHOD==.true., those variables are used to store the value of the tensor field TF across the element's boundary, in order to compute the flux.
  integer :: ispec,i,j,k,iglob, iglobM, iglobP!, iglob_unique
  real(kind=CUSTOM_REAL) :: weight!, flux_n, flux_x, flux_z,  !nx, nz, !rho_DG_P, rhovx_DG_P, rhovz_DG_P, &
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
  !       nabla_SF(1,:), nabla_SF(NDIM,:), nabla_TF(1,1,:), nabla_TF(NDIM,NDIM,:), nabla_TF(1,NDIM,:), nabla_TF(NDIM,1,:)
  
  !real(kind=CUSTOM_REAL), dimension(NDIM,nglob_DG) :: locgrad_SF
  !real(kind=CUSTOM_REAL), dimension(NDIM,NDIM,nglob_DG) :: locgrad_TF
  
  !real(kind=CUSTOM_REAL) :: dxiTF(1), dgamTF(1), dxiTF(NDIM), dgamTF(NDIM), dSF_dxi, dSF_dgamma ! Derivatives in \Lambda.
  real(kind=CUSTOM_REAL) :: dSF_dxi, dSF_dgamma ! Derivatives in \Lambda.
  real(kind=CUSTOM_REAL), dimension(NDIM) :: dxiTF, dgamTF ! Derivatives in \Lambda.
  !real(kind=CUSTOM_REAL) :: dTFx_dx, dTFx_dz, dTFz_dx, dTFz_dz, dSF_dx, dSF_dz
  real(kind=CUSTOM_REAL) :: wxl, wzl ! Quadrature weights.
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: temp_SFx_1, temp_SFx_2, &
        temp_SFz_1, temp_SFz_2, temp_TFxx_1, temp_TFxx_2, &
        temp_TFxz_1, temp_TFxz_2, temp_TFzx_1, temp_TFzx_2, temp_TFzz_1, temp_TFzz_2
  !real(kind=CUSTOM_REAL) :: vx_init, vz_init
  !real(kind=CUSTOM_REAL) :: ya_x_l, ya_z_l
  
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
    !veloc_x_DG = rhovx_DG/rho_DG
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
          
          wzl = wzgll(j)
          wxl = wxgll(i)
          
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
              temp_TFzx_1(i,j) = wzl * jacLoc * (xixl * TF(NDIM,iglob))
              temp_TFzz_1(i,j) = wzl * jacLoc * (xizl * TF(NDIM,iglob))
              temp_TFxx_2(i,j) = wxl * jacLoc * (gammaxl * TF(1,iglob))
              temp_TFxz_2(i,j) = wxl * jacLoc * (gammazl * TF(1,iglob))
              temp_TFzx_2(i,j) = wxl * jacLoc * (gammaxl * TF(NDIM,iglob))
              temp_TFzz_2(i,j) = wxl * jacLoc * (gammazl * TF(NDIM,iglob))
            endif
            !else
            !  vx_init = rhovx_init(iglob)/rho_init(iglob)
            !  vz_init = rhovz_init(iglob)/rho_init(iglob)
            !  
            !  temp_SFx_1(i,j)  = wzl * jacLoc * (xixl * (SF(iglob) - T_init(iglob)))
            !  temp_SFz_1(i,j)  = wzl * jacLoc * (xizl * (SF(iglob) - T_init(iglob)))
            !  temp_TFxx_1(i,j) = wzl * jacLoc * (xixl * (TF(1,iglob) - vx_init))
            !  temp_TFxz_1(i,j) = wzl * jacLoc * (xizl * (TF(1,iglob))) ! Some hypothesis on initial velocity is used here.
            !  temp_TFzx_1(i,j) = wzl * jacLoc * (xixl * (TF(NDIM,iglob) - vz_init))
            !  temp_TFzz_1(i,j) = wzl * jacLoc * (xizl * (TF(NDIM,iglob) - vz_init))
            !  
            !  temp_SFx_2(i,j)  = wxl * jacLoc * (gammaxl * (SF(iglob) - T_init(iglob)))
            !  temp_SFz_2(i,j)  = wxl * jacLoc * (gammazl * (SF(iglob) - T_init(iglob)))
            !  temp_TFxx_2(i,j) = wxl * jacLoc * (gammaxl * (TF(1,iglob) - vx_init))
            !  temp_TFxz_2(i,j) = wxl * jacLoc * (gammazl * (TF(1,iglob))) ! Some hypothesis on initial velocity is used here.
            !  temp_TFzx_2(i,j) = wxl * jacLoc * (gammaxl * (TF(NDIM,iglob) - vz_init))
            !  temp_TFzz_2(i,j) = wxl * jacLoc * (gammazl * (TF(NDIM,iglob) - vz_init))
            !endif
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
                dxiTF(NDIM)  = dxiTF(NDIM)  + TF(NDIM,ibool_DG(k,j,ispec)) * real(hprime_xx(i,k), kind=CUSTOM_REAL)
                dgamTF(NDIM) = dgamTF(NDIM) + TF(NDIM,ibool_DG(i,k,ispec)) * real(hprime_zz(j,k), kind=CUSTOM_REAL)
              endif
              !else
              !  vx_init = rhovx_init(ibool_DG(k,j,ispec))/rho_init(ibool_DG(k,j,ispec))
              !  dxiTF(1) = dxiTF(1) + (TF(1,ibool_DG(k,j,ispec)) - vx_init) &
              !                      * real(hprime_xx(i,k), kind=CUSTOM_REAL)
              !  vx_init = rhovx_init(ibool_DG(i,k,ispec))/rho_init(ibool_DG(i,k,ispec))       
              !  dgamTF(1) = dgamTF(1) + (TF(1,ibool_DG(i,k,ispec)) - vx_init) &
              !                            * real(hprime_zz(j,k), kind=CUSTOM_REAL)
              !  vz_init = rhovz_init(ibool_DG(k,j,ispec))/rho_init(ibool_DG(k,j,ispec))
              !  dxiTF(NDIM) = dxiTF(NDIM) + (TF(NDIM,ibool_DG(k,j,ispec)) - vz_init) &
              !                      * real(hprime_xx(i,k), kind=CUSTOM_REAL)
              !  vz_init = rhovz_init(ibool_DG(i,k,ispec))/rho_init(ibool_DG(i,k,ispec))
              !  dgamTF(NDIM) = dgamTF(NDIM) + (TF(NDIM,ibool_DG(i,k,ispec)) - vz_init) &
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
              !nabla_SF(NDIM,iglob) = dSF_dz
              nabla_SF(1,iglob) = dSF_dxi * xixl + dSF_dgamma * gammaxl
              nabla_SF(NDIM,iglob) = dSF_dxi * xizl + dSF_dgamma * gammazl
            endif
            if(swTF) then
              !dTFx_dx   = dxiTF(1) * xixl + dgamTF(1) * gammaxl
              !dTFx_dz   = dxiTF(1) * xizl + dgamTF(1) * gammazl
              !dTFz_dx   = dxiTF(NDIM) * xixl + dgamTF(NDIM) * gammaxl
              !dTFz_dz   = dxiTF(NDIM) * xizl + dgamTF(NDIM) * gammazl
              !temp_TFxx = dTFx_dx * jacLoc
              !temp_TFxz = dTFx_dz * jacLoc
              !temp_TFzx = dTFz_dx * jacLoc
              !temp_TFzz = dTFz_dz * jacLoc
              !nabla_TF(1,1,iglob) = dTFx_dx
              !nabla_TF(1,NDIM,iglob) = dTFx_dz
              !nabla_TF(NDIM,1,iglob) = dTFz_dx
              !nabla_TF(NDIM,NDIM,iglob) = dTFz_dz
              do k=1,NDIM ! Re-using index k for spatial dimension.
                nabla_TF(k,1,iglob) = dxiTF(k) * xixl + dgamTF(k) * gammaxl
                nabla_TF(k,NDIM,iglob) = dxiTF(k) * xizl + dgamTF(k) * gammazl
                !nabla_TF(1,1,iglob) = dxiTF(1) * xixl + dgamTF(1) * gammaxl
                !nabla_TF(1,NDIM,iglob) = dxiTF(1) * xizl + dgamTF(1) * gammazl
                !nabla_TF(NDIM,1,iglob) = dxiTF(NDIM) * xixl + dgamTF(NDIM) * gammaxl
                !nabla_TF(NDIM,NDIM,iglob) = dxiTF(NDIM) * xizl + dgamTF(NDIM) * gammazl
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
              nabla_SF(NDIM,iglob) = nabla_SF(NDIM,iglob) - &
                               (temp_SFz_1(k,j) * real(hprimewgll_xx(k,i), kind=CUSTOM_REAL) + &
                                temp_SFz_2(i,k) * real(hprimewgll_zz(k,j), kind=CUSTOM_REAL))
              endif
              if(swTF) then
              nabla_TF(1,1,iglob) = nabla_TF(1,1,iglob) - &
                                (temp_TFxx_1(k,j) * real(hprimewgll_xx(k,i), kind=CUSTOM_REAL) + &
                                 temp_TFxx_2(i,k) * real(hprimewgll_zz(k,j), kind=CUSTOM_REAL))
              nabla_TF(1,NDIM,iglob) = nabla_TF(1,NDIM,iglob) - &
                                (temp_TFxz_1(k,j) * real(hprimewgll_xx(k,i), kind=CUSTOM_REAL) + &
                                 temp_TFxz_2(i,k) * real(hprimewgll_zz(k,j), kind=CUSTOM_REAL))
              nabla_TF(NDIM,1,iglob) = nabla_TF(NDIM,1,iglob) - &
                                (temp_TFzx_1(k,j) * real(hprimewgll_xx(k,i), kind=CUSTOM_REAL) + &
                                 temp_TFzx_2(i,k) * real(hprimewgll_zz(k,j), kind=CUSTOM_REAL))
              nabla_TF(NDIM,NDIM,iglob) = nabla_TF(NDIM,NDIM,iglob) - &
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
            if(swSF) then
              SF_P         = ZEROcr
            endif
            if(swTF) then
              TF_P         = ZEROcr
            endif
            !p_DG_P       = ZEROcr
            !Vxx_DG_P     = ZEROcr
            !Vzz_DG_P     = ZEROcr
            !Vzx_DG_P     = ZEROcr
            !Vxz_DG_P     = ZEROcr
            n_out(1) = nx_iface(iface, ispec)
            n_out(2) = nz_iface(iface, ispec)
            
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
          !call compute_interface_unknowns(i,j,ispec, drho_P, rho0dv_P(1), &
          !        rho0dv_P(2), dE_P, TF_P(1), TF_P(NDIM), in_dp_P, T_P, &
          !        Tx_DG_P, Tz_DG_P, Vxx_DG_P, Vzz_DG_P, Vzx_DG_P, Vxz_DG_P, gamma_P,&
          !        neighbor, &
          !        exact_interface_flux, &
          !        cv_drho(iglobM), cv_dE(iglobM), cv_rho0dv(1,iglobM), cv_rho0dv(NDIM,iglobM), &
          !        V_DG(:,:,iglobM), T_DG(:,iglobM), &
          !        cv_drho(iglobP), cv_dE(iglobP), cv_rho0dv(1,iglobP), cv_rho0dv(NDIM,iglobP), &
          !        V_DG(:,:,iglobP), T_DG(:,iglobP), &
          !        nx, nz, weight, currentTime, iface1, iface)
          !        !TEST STRETCH
          !        !nx_unit, nz_unit, weight, currentTime, iface1, iface)
          !call LNS_get_interfaces_unknowns(i, j, ispec, &
          !        drho_P, rho0dv_P(1), &
          !        rho0dv_P(2), dE_P, TF_P(1), TF_P(NDIM), in_dp_P, T_P, &
          !        Tx_DG_P, Tz_DG_P, Vxx_DG_P, Vzz_DG_P, Vzx_DG_P, Vxz_DG_P, gamma_P,&
          !        neighbor, &
          !        exact_interface_flux, &
          !        cv_drho(iglobM), cv_dE(iglobM), cv_rho0dv(1,iglobM), cv_rho0dv(NDIM,iglobM), &
          !        V_DG(:,:,iglobM), T_DG(:,iglobM), &
          !        cv_drho(iglobP), cv_dE(iglobP), cv_rho0dv(1,iglobP), cv_rho0dv(NDIM,iglobP), &
          !        V_DG(:,:,iglobP), T_DG(:,iglobP), &
          !        nx, nz, weight, currentTime, iface1, iface)
            
            call LNS_get_interfaces_unknowns(i, j, ispec, iface1, iface, neighbor, timelocal, & ! Point identifier (input).
                  LNS_drho(iglobM), LNS_rho0dv(:,iglobM), & ! Input constitutive variables, "M" side.
                  LNS_drho(iglobP), LNS_rho0dv(:,iglobP), LNS_dE(iglobP), & ! Input constitutive variables, "P" side.
                  LNS_dummy_1d(1), & ! Input other variable, "M" side.
                  !V_DG(:,:,iglobM), T_DG(:,iglobM), & ! Input derivatives, "M" side. MIGHT NEED.
                  !V_DG(:,:,iglobP), T_DG(:,iglobP), & ! Input derivatives, "M" side. MIGHT NEED.
                  n_out, & ! Normal vector (input).
                  exact_interface_flux, & ! Switch to disable jump in some cases (output).
                  LNS_dummy_1d(1), LNS_dummy_2d(:,1), LNS_dummy_1d(1), & ! Output constitutive variables.
                  !Tx_DG_P, Tz_DG_P, Vxx_DG_P, Vzz_DG_P, Vzx_DG_P, Vxz_DG_P,& ! Output derivatives. MIGHT NEED.
                  LNS_dummy_2d(:,1), LNS_dummy_1d(1), TF_P, & ! Output other variables.
                  .false., LNS_dummy_2d(:,1), LNS_dummy_1d(1:3), & ! Output other variables: viscous.
                  .true., SF_P) ! Output other variables.
            ! We know that the scalar field will always be temperature, and that the tensor field will always be the velocity, so we use the already-built LNS_get_interfaces_unknowns routine to get them. The switch is here to prevent unecessary quantites to be computed during that specific call. If other fields need to be computed, one will have to use a dedicated routine.

            !vx_init = rhovx_init(iglobM)/rho_init(iglobM)
            !vz_init = rhovz_init(iglobM)/rho_init(iglobM)
            
            if(swTF) then
              if(abs(timelocal)<TINYVAL) then
                ! At timelocal==0, we want to compute \nabla v_0, and thus we need to add v_0 to TF_P. At timelocal==0, we do not care about \nabla T.
                ! For timelocal>0, we know we are computing \nabla v' or \nabla T', and thus we do not need to add anything.
                TF_P=TF_P+LNS_v0(:,iglobM)
              endif
            endif
            ! Dot products.
            do i=1,NDIM
              if(swSF) then
                nabla_SF(i,iglobM) = nabla_SF(i,iglobM) + weight*(SF(iglobM)+SF_P)*n_out(i)*HALFcr
              endif
              if(swTF) then
                do j=1,NDIM
                  nabla_TF(i,j,iglobM) = nabla_TF(i,j,iglobM) + weight*(TF(i,iglobM)+TF_P(i))*n_out(j)*HALFcr
                enddo
              endif
            enddo
            !if(swSF) then
            !  !!flux_x = SF(iglobM) + SF_P ! Once multiplied with HALFcr below, will represent flux along x of the scalar field SF.
            !  !!!if(CONSTRAIN_HYDROSTATIC) flux_x = flux_x - 2*T_init(iglobM)
            !  !!flux_n = flux_x*n_out(1) ! Recall: n_out(1)=n_x.
            !  !flux_n = (SF(iglobM)+SF_P)*n_out(1) ! Recall: n_out(1)=n_x.
            !  !nabla_SF(1,iglobM) = nabla_SF(1,iglobM) + weight*flux_n*HALFcr
            !  !!flux_z = SF(iglobM) + SF_P ! Once multiplied with HALFcr below, will represent flux along z of the scalar field SF.
            !  !!!if(CONSTRAIN_HYDROSTATIC) flux_z = flux_z - 2*T_init(iglobM)
            !  !!flux_n = flux_z*n_out(NDIM) ! Recall: n_out(NDIM)=n_z.
            !  !flux_n = (SF(iglobM)+SF_P)*n_out(NDIM) ! Recall: n_out(NDIM)=n_z.
            !  !nabla_SF(NDIM,iglobM) = nabla_SF(NDIM,iglobM) + weight*flux_n*HALFcr
            !  do i=1,NDIM
            !    !flux_n = (SF(iglobM)+SF_P)*n_out(i) ! Recall: n_out(1)=n_x.
            !    !nabla_SF(i,iglobM) = nabla_SF(i,iglobM) + weight*flux_n*HALFcr
            !    nabla_SF(i,iglobM) = nabla_SF(i,iglobM) + weight*(SF(iglobM)+SF_P)*n_out(i)*HALFcr
            !  enddo
            !endif
            !if(swTF) then
            !  !!flux_x = TF(1,iglobM) + TF_P(1) ! Once multiplied with HALFcr below, will represent flux along x of the x-component of the tensor field TF.
            !  !!!if(CONSTRAIN_HYDROSTATIC) flux_x = flux_x - 2*vx_init
            !  !!flux_n = flux_x*n_out(1) ! Recall: n_out(1)=n_x.
            !  !flux_n = (TF(1,iglobM)+TF_P(1))*n_out(1) ! Recall: n_out(1)=n_x.
            !  !nabla_TF(1,1,iglobM) = nabla_TF(1,1,iglobM) + weight*flux_n*HALFcr
            
            !  !flux_z = TF(1,iglobM) + TF_P(1) ! Once multiplied with HALFcr below, will represent flux along z of the x-component of the tensor field TF.
            !  !!if(CONSTRAIN_HYDROSTATIC) flux_z = flux_z - 2*vx_init
            !  !flux_n = flux_z*n_out(NDIM) ! Recall: n_out(NDIM)=n_z.
            !  !nabla_TF(1,NDIM,iglobM) = nabla_TF(1,NDIM,iglobM) + weight*flux_n*HALFcr
            
            !  !flux_x = TF(NDIM,iglobM) + TF_P(NDIM) ! Once multiplied with HALFcr below, will represent flux along x of the z-component of the tensor field TF.
            !  !!if(CONSTRAIN_HYDROSTATIC) flux_x = flux_x - 2*vz_init
            !  !flux_n = flux_x*n_out(1) ! Recall: n_out(1)=n_x.
            !  !nabla_TF(NDIM,1,iglobM) = nabla_TF(NDIM,1,iglobM) + weight*flux_n*HALFcr
            
            !  !flux_z = TF(NDIM,iglobM) + TF_P(NDIM) ! Once multiplied with HALFcr below, will represent flux along z of the z-component of the tensor field TF.
            !  !!if(CONSTRAIN_HYDROSTATIC) flux_z = flux_z - 2*vz_init
            !  !flux_n = flux_z*n_out(NDIM) ! Recall: n_out(NDIM)=n_z.
            !  !nabla_TF(NDIM,NDIM,iglobM) = nabla_TF(NDIM,NDIM,iglobM) + weight*flux_n*HALFcr
            !  do i=1,NDIM
            !    do j=1,NDIM
            !      nabla_TF(i,j,iglobM) = nabla_TF(i,j,iglobM) + weight*(TF(i,iglobM)+TF_P(i))*n_out(j)*HALFcr
            !    enddo
            !  enddo
            !endif
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
          !nabla_TF(i,NDIM,:) = nabla_TF(i,NDIM,:) * rmass_inverse_acoustic_DG(:)
          !nabla_TF(1,1,:) = nabla_TF(1,1,:) * rmass_inverse_acoustic_DG(:)
          !nabla_TF(1,NDIM,:) = nabla_TF(1,NDIM,:) * rmass_inverse_acoustic_DG(:)
          !nabla_TF(NDIM,1,:) = nabla_TF(NDIM,1,:) * rmass_inverse_acoustic_DG(:)
          !nabla_TF(NDIM,NDIM,:) = nabla_TF(NDIM,NDIM,:) * rmass_inverse_acoustic_DG(:)
        enddo
      endif
    enddo
  endif
  
  !! Store in variables intended for output.
  !if(swSF) then
  !  nabla_SF(1, :) = nabla_SF(1,:)
  !  nabla_SF(2, :) = nabla_SF(NDIM,:)
  !endif
  !if(swTF) then
  !  nabla_TF(1, 1, :) = nabla_TF(1,1,:) ! dxTFx
  !  nabla_TF(1, NDIM, :) = nabla_TF(1,NDIM,:) ! dzTFx
  !  nabla_TF(NDIM, 1, :) = nabla_TF(NDIM,1,:) ! dxTFz
  !  nabla_TF(NDIM, NDIM, :) = nabla_TF(NDIM,NDIM,:) ! dzTFz
  !endif
end subroutine compute_gradient_TFSF




















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
  use constants, only: CUSTOM_REAL, NDIM
  use specfem_par, only: c_V, nglob_DG
  use specfem_par_LNS, only: norm2
  implicit none
  ! Input/Output.
  real(kind=CUSTOM_REAL), dimension(nglob_DG), intent(in) :: in_rho, in_E
  real(kind=CUSTOM_REAL), dimension(NDIM, nglob_DG), intent(in) :: in_v
  real(kind=CUSTOM_REAL), dimension(nglob_DG), intent(out) :: out_T
  ! Local.
  real(kind=CUSTOM_REAL), parameter :: HALFcr = 0.5_CUSTOM_REAL
  !out_T=(in_E/in_rho - HALFcr*(in_v(1,:)**2+in_v(NDIM,:)**2))/c_V
  out_T=(in_E/in_rho - HALFcr*norm2(in_v))/c_V
end subroutine compute_T
! ------------------------------------------------------------ !
! compute_T_i                                                  !
! ------------------------------------------------------------ !
! Same as compute_T, but point by point (unvectorised).
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
! ------------------------------------------------------------ !
! compute_dT                                                   !
! ------------------------------------------------------------ !
! Computes temperature perturbation from constitutive variables.
! Note: this could have been done nearly inline by using the subroutine compute_T, but defining this function enables one to use less RAM.
subroutine compute_dT(in_rho, in_v, in_E, out_dT)
  use constants, only: CUSTOM_REAL, NDIM
  use specfem_par, only: c_V, nglob_DG
  use specfem_par_LNS, only: LNS_T0, norm2
  implicit none
  ! Input/Output.
  real(kind=CUSTOM_REAL), dimension(nglob_DG), intent(in) :: in_rho, in_E
  real(kind=CUSTOM_REAL), dimension(NDIM, nglob_DG), intent(in) :: in_v
  real(kind=CUSTOM_REAL), dimension(nglob_DG), intent(out) :: out_dT
  ! Local.
  real(kind=CUSTOM_REAL), parameter :: HALFcr = 0.5_CUSTOM_REAL
  !out_dT=  (in_E/in_rho - HALFcr*(in_v(1,:)**2+in_v(NDIM,:)**2))/c_V &
  out_dT=  (in_E/in_rho - HALFcr*norm2(in_v))/c_V &
         - LNS_T0
end subroutine compute_dT
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
  !out_T=(in_E/in_rho - HALFcr*(in_v(1)**2+in_v(NDIM)**2))/c_V
  !out_T=(in_E/in_rho - HALFcr*minval(norm2(reshape(in_v,(/2,1/)))))/c_V ! Could not make our "norm2" function work both with rank 1 arrays and rank 2 arrays, so used a trick: reshape the 2-sized vector (rank 1) to a (2,1) matrix (rank 2), send it to norm2, retrieve a 1-sized vector (rank 1), use minval to send it back to a scalar value (rank 0) as needed. It proved very unefficient in terms of computation time, so we defined another dedicated "norm2" function.
  out_T=(in_E/in_rho - HALFcr*norm2r1(in_v))/c_V &
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
          !+ in_rho*HALFcr*( in_v(1,:)**2 + in_v(NDIM,:)**2 )
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
          !+ in_rho*HALFcr*( in_v(1,:)**2 + in_v(NDIM,:)**2 )
          + in_rho*HALFcr*norm2r1(in_v)
end subroutine compute_E_i
! ------------------------------------------------------------ !
! compute_dE_i                                                 !
! ------------------------------------------------------------ !
! Same as compute_dE, but point by point (unvectorised).
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
          !+ in_rho*HALFcr*( in_v(1,:)**2 + in_v(NDIM,:)**2 )
          + in_rho*HALFcr*norm2r1(in_v) &
          - LNS_E0(iglob)
end subroutine compute_dE_i


























! ------------------------------------------------------------ !
! LNS_prevent_nonsense                                         !
! ------------------------------------------------------------ !
! Attempts to detect nonsense in LNS variables, and stop program if they are found.
subroutine LNS_prevent_nonsense()
  ! TODO: select variables to use.
  use constants
  use specfem_par
  use specfem_par_LNS
  implicit none
  ! Input/Output.
  ! N./A.
  ! Local.
  ! N./A.
  
  ! Initial state.
  if(minval(LNS_rho0)<TINYVAL) then
    stop "LNS_rho0 is non-positive (<= 0) somewhere."
  endif
  if(minval(LNS_E0)<TINYVAL) then
    stop "LNS_E0 is non-positive (<= 0) somewhere."
  endif
  
  ! Auxiliary initial variables.
  if(minval(LNS_p0)<TINYVAL) then
    stop "LNS_p0 is non-positive (<= 0) somewhere."
  endif
  
  ! Physical parameters.
  if(LNS_viscous) then ! Check if viscosity exists whatsoever.
    if(minval(LNS_mu)<-TINYVAL) then
      stop "LNS_mu is negative (< 0) somewhere."
    endif
    if(minval(LNS_eta)<-TINYVAL) then
      stop "LNS_eta is negative (< 0) somewhere."
    endif
    if(minval(LNS_kappa)<-TINYVAL) then
      stop "LNS_kappa is negative (< 0) somewhere."
    endif
  endif
end subroutine LNS_prevent_nonsense

! ------------------------------------------------------------ !
! LNS_warn_nonsense                                            !
! ------------------------------------------------------------ !
! Attempts to detect nonsense in LNS variables, and warn user if they are found.
subroutine LNS_warn_nonsense()
  ! TODO: select variables to use.
  use constants
  use specfem_par
  use specfem_par_LNS
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
  
  !write(*,*) i,j,ispec,ispec_is_acoustic_DG(ispec) ! DEBUG
  
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










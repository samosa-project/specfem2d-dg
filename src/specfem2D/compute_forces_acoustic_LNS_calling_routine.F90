! ------------------------------------------------------------ !
! compute_forces_acoustic_LNS_main                             !
! ------------------------------------------------------------ !
! TODO: Description.
! 
! Initialisation (routine initial_state_LNS()) was done way before, in
! prepare_timerun().
! 
! This routine is called by iterate_time().

subroutine compute_forces_acoustic_LNS_main()
  ! TODO: select variables to use.
  use constants
  use specfem_par
  use specfem_par_LNS

  implicit none

  ! Local.
  real(kind=CUSTOM_REAL), parameter :: ZEROcr = 0._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: ONEcr = 1._CUSTOM_REAL
  real(kind=CUSTOM_REAL) :: timelocal
  real(kind=CUSTOM_REAL), dimension(NDIM,NDIM,nglob_DG) :: nabla_dv
  integer :: ier, i_aux!, i, j, ispec, iglob
  logical check_linearHypothesis
  
  ! Checks if anything has to be done.
  if (.not. any_acoustic_DG) then
    return
  endif
  
  ! Debug switches. They are computationaly very heavy and should not be used on every simulation.
  check_linearHypothesis=.true.
  
  ! Compute current time.
  timelocal = (it-1)*deltat + LNS_scheme_C(i_stage)*deltat
  
  ! Intialisation.
  if(it == 1 .and. i_stage == 1) then
    !if(USE_SLOPE_LIMITER) then
    !  ! The Vandermonde matrices are only used when the slope limiter is.
    !  call setUpVandermonde()
    !endif
    
    !if(USE_DISCONTINUOUS_METHOD .and. USE_LNS) then
    !  ! The subroutine 'initial_state_LNS' needs to have stretching initialised to compute \nabla\bm{v}_0. This is why it is put here.
    !  call initial_state_LNS() ! This routine can be found in compute_forces_acoustic_LNS.F90.
    !endif
    
    ! Initialise state registers. Note: since constitutive variables are perturbations, they are necessarily zero at start.
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
    WRITE(*, *) "****************************************************************"
    WRITE(*, "(a,i9,a,i1,a,e23.16,a)") "Iteration", it, ", stage ", i_stage, ", local time", timelocal, "."
  endif
  
  ! Precompute momentum perturbation, velocity perturbation,
  ! temperature perturbation, and pressure perturbation.
  LNS_dv = ZEROcr
  do i_aux = 1, NDIM
    LNS_dm(i_aux, :) = LNS_rho0dv(i_aux, :) + LNS_drho*LNS_v0(i_aux, :)
    where(LNS_rho0/=ZEROcr) LNS_dv(i_aux, :) = LNS_rho0dv(i_aux, :)/LNS_rho0 ! 'where(...)' as safeguard, as rho0=0 should not happen.
  enddo
  call compute_dp(LNS_rho0+LNS_drho, LNS_v0+LNS_dv, LNS_E0+LNS_dE, LNS_dp)
#if 0
  ! TEST MANUFACTURED SOLUTIONS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Exact formulas, for manufactured solutions' validation check.
  !write(*,*) "minval(ibool_DG), maxval(ibool_DG), nglob_DG", minval(ibool_DG), maxval(ibool_DG), nglob_DG
  do ispec = 1, nspec; do j = 1, NGLLZ; do i = 1, NGLLX ! V1: get on all GLL points
    !if(     abs(coord(1,ibool_before_perio(i,j,ispec))-0.)<TINYVAL &
    !   .or. abs(coord(1,ibool_before_perio(i,j,ispec))-1.)<TINYVAL &
    !   .or. abs(coord(2,ibool_before_perio(i,j,ispec))-0.)<TINYVAL &
    !   .or. abs(coord(2,ibool_before_perio(i,j,ispec))-1.)<TINYVAL) then
      iglob = ibool_DG(i,j,ispec)
      LNS_drho(iglob) = 0.
#define aa 2.
#define cc 3.
#define ee 5.
#define gg 7.
#define bb 4.
#define dd 4.
#define ff 4.
#define hh 4.
      LNS_dv(1,iglob) = aa*coord(1,ibool_before_perio(i,j,ispec))**bb + cc*coord(2,ibool_before_perio(i,j,ispec))**dd
      LNS_dv(2,iglob) = ee*coord(1,ibool_before_perio(i,j,ispec))**ff + gg*coord(2,ibool_before_perio(i,j,ispec))**hh
      !LNS_dv(1,iglob) = coord(1,ibool_before_perio(i,j,ispec))+coord(2,ibool_before_perio(i,j,ispec))
      !LNS_dv(2,iglob) = coord(1,ibool_before_perio(i,j,ispec))+coord(2,ibool_before_perio(i,j,ispec))
      LNS_rho0dv(:,iglob) = LNS_rho0(iglob)*LNS_dv(:,iglob)
      LNS_dE(iglob) = 0.
      !LNS_dE(iglob) = coord(2,ibool_before_perio(i,j,ispec)) - 0.5
      !LNS_dp(iglob) = (gammaext_DG(iglob)-1.)*(  sin(1.5*PI*coord(1,ibool_before_perio(i,j,ispec))) &
      !                                         + sin(1.5*PI*coord(2,ibool_before_perio(i,j,ispec))))
      !LNS_dp(iglob) = (gammaext_DG(iglob)-1.)*LNS_dE(iglob)
      LNS_dp(iglob) = 0.
    !endif
  enddo; enddo; enddo
#if 0
    do ispec = 1, nspec; do j = 1, NGLLZ; do i = 1, NGLLX ! V1: get on all GLL points
      if(abs(coord(2,ibool_before_perio(i,j,ispec))-2.5)<TINYVAL) then
        write(*,*) '(x z) (vx vz)', coord(:,ibool_before_perio(i,j,ispec)), LNS_dv(:,ibool_DG(i,j,ispec))
      endif
    enddo; enddo; enddo
#endif
  ! TEST MANUFACTURED SOLUTIONS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#endif
  !call compute_dT(LNS_rho0+LNS_drho, LNS_p0+LNS_dp, LNS_dT)
  call compute_dT(LNS_rho0+LNS_drho, LNS_v0+LNS_dv, LNS_E0+LNS_dE, LNS_dT)
  
  ! Precompute gradients, only if viscosity exists whatsoever.
  if(LNS_viscous) then
    ! Note: this will compute nabla_dv only if viscosity is activated.
    call compute_gradient_TFSF(LNS_dv, LNS_dT, LNS_viscous, LNS_viscous, LNS_switch_gradient, nabla_dv, nabla_dT, timelocal)
    !call compute_viscous_tensors(nabla_dT, nabla_dv, LNS_dummy_1d, LNS_dummy_1d, LNS_rho0, &
    !                             LNS_rho0dv(1,:), LNS_rho0dv(2,:), LNS_dE, timelocal) ! Litteraly no difference.
    !call compute_gradient_TFSF(LNS_dummy_2d, LNS_dv(1,:), .false., .true., LNS_switch_gradient, &
    !                           LNS_dummy_2d, nabla_dv(1,:,:), timelocal)
    !call compute_gradient_TFSF(LNS_dummy_2d, LNS_dv(2,:), .false., .true., LNS_switch_gradient, &
    !                           LNS_dummy_2d, nabla_dv(2,:,:), timelocal) ! Litteraly no difference.
#if 0
    do ispec = 1, nspec; do j = 1, NGLLZ; do i = 1, NGLLX ! V1: get on all GLL points
      !if(abs(coord(1,ibool_before_perio(i,j,ispec))-0.)<TINYVAL .and. abs(coord(2,ibool_before_perio(i,j,ispec))-2.5)<TINYVAL) then
      if(abs(nabla_dv(1, 1, ibool_DG(i,j,ispec))-aa*bb*coord(1,ibool_before_perio(i,j,ispec))**(bb-1))>TINYVAL .or. &
         abs(nabla_dv(1, 2, ibool_DG(i,j,ispec))-cc*dd*coord(2,ibool_before_perio(i,j,ispec))**(dd-1))>TINYVAL .or. &
         abs(nabla_dv(2, 1, ibool_DG(i,j,ispec))-ee*ff*coord(1,ibool_before_perio(i,j,ispec))**(ff-1))>TINYVAL .or. &
         abs(nabla_dv(2, 2, ibool_DG(i,j,ispec))-gg*hh*coord(2,ibool_before_perio(i,j,ispec))**(hh-1))>TINYVAL ) then
        write(*,*) '(x z) (dxvx dzvx dxvz dzvz)',&
                                            coord(:,ibool_before_perio(i,j,ispec)), &
                                            nabla_dv(1, 1, ibool_DG(i,j,ispec)), nabla_dv(1, 2, ibool_DG(i,j,ispec)), & ! TEST
                                            nabla_dv(2, 1, ibool_DG(i,j,ispec)), nabla_dv(2, 2, ibool_DG(i,j,ispec)) ! TEST
        write(*,*) 'th =                       ', &
                                                  coord(:,ibool_before_perio(i,j,ispec)), &
                                                  aa*bb*coord(1,ibool_before_perio(i,j,ispec))**(bb-1), &
                                                  cc*dd*coord(2,ibool_before_perio(i,j,ispec))**(dd-1), &
                                                  ee*ff*coord(1,ibool_before_perio(i,j,ispec))**(ff-1), &
                                                  gg*hh*coord(2,ibool_before_perio(i,j,ispec))**(hh-1)
        endif
      !endif
    enddo; enddo; enddo
#endif
    ! Precompute \Sigma_v' from \nabla v'.
    call LNS_compute_viscous_stress_tensor(nabla_dv, sigma_dv)
#if 0
    do ispec = 1, nspec; do j = 1, NGLLZ; do i = 1, NGLLX ! V1: get on all GLL points
      if(isNotClose(  sigma_dv(1, ibool_DG(i,j,ispec)) &
                    , (8./3.)*LNS_mu(ibool_DG(i,j,ispec))*(  aa*bb*coord(1,ibool_before_perio(i,j,ispec))**(bb-1) &
                                                           + 0.25*gg*hh*coord(2,ibool_before_perio(i,j,ispec))**(hh-1)) &
                    , TINYVAL) .or. &
         isNotClose(  sigma_dv(2, ibool_DG(i,j,ispec)) &
                    , LNS_mu(ibool_DG(i,j,ispec))*(  cc*dd*coord(2,ibool_before_perio(i,j,ispec))**(dd-1) &
                                                   + ee*ff*coord(1,ibool_before_perio(i,j,ispec))**(ff-1)) &
                    , TINYVAL) .or. &
         isNotClose(  sigma_dv(3, ibool_DG(i,j,ispec)) &
                    , (8./3.)*LNS_mu(ibool_DG(i,j,ispec))*(  gg*hh*coord(2,ibool_before_perio(i,j,ispec))**(hh-1) &
                                                           + 0.25*aa*bb*coord(1,ibool_before_perio(i,j,ispec))**(bb-1)) &
                    , TINYVAL)) then
        write(*,*) '(x z) (Sv_11 Sv_12 Sv_22)', coord(:,ibool_before_perio(i,j,ispec)), sigma_dv(:, ibool_DG(i,j,ispec)) ! TEST
        write(*,*) 'th =                     ', coord(:,ibool_before_perio(i,j,ispec)), &
                   (8./3.)*LNS_mu(ibool_DG(i,j,ispec))*(  aa*bb*coord(1,ibool_before_perio(i,j,ispec))**(bb-1) &
                                                        + 0.25*gg*hh*coord(2,ibool_before_perio(i,j,ispec))**(hh-1)), &
                   LNS_mu(ibool_DG(i,j,ispec))*(  cc*dd*coord(2,ibool_before_perio(i,j,ispec))**(dd-1) &
                                                + ee*ff*coord(1,ibool_before_perio(i,j,ispec))**(ff-1)), &
                   (8./3.)*LNS_mu(ibool_DG(i,j,ispec))*(  gg*hh*coord(2,ibool_before_perio(i,j,ispec))**(hh-1) &
                                                        + 0.25*aa*bb*coord(1,ibool_before_perio(i,j,ispec))**(bb-1))
      endif
    enddo; enddo; enddo
#endif
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
  
  !write(*,*) "minval(LNS_dv), maxval(LNS_dv)", minval(LNS_dv), maxval(LNS_dv)! DEBUG
  
  ! Compute RHS.
  ! Note: if there is PML BCs, additionnal terms will be queried directly inside the call to compute_forces_acoustic_LNS, since auxiliary variables and coefficients are global variables.
  call compute_forces_acoustic_LNS(LNS_drho, LNS_rho0dv, LNS_dE, & ! Constitutive variables.
                                   LNS_dm, LNS_dp, nabla_dT, sigma_dv, & ! Precomputed quantities. sigma_dv is sent even if viscosity is deactivated, but in that case it should be zero and unused in the subroutine.
                                   RHS_drho, RHS_rho0dv, RHS_dE, & ! Output.
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
  ! Group momentum treatment in loop on dimension.
  do i_aux=1,NDIM
    RHS_rho0dv(i_aux,:) = RHS_rho0dv(i_aux,:) * rmass_inverse_acoustic_DG(:)
    aux_rho0dv(i_aux,:) = LNS_scheme_A(i_stage)*aux_rho0dv(i_aux,:) + deltat*RHS_rho0dv(i_aux,:)
    LNS_rho0dv(i_aux,:) = LNS_rho0dv(i_aux,:) + LNS_scheme_B(i_stage)*aux_rho0dv(i_aux,:)
  enddo
  ! If one wants to do all at once (method which should be working but seems to be somewhat more unstable), comment the previous lines and uncomment the following lines.
  !LNS_drho               =   LNS_drho &
  !                         + LNS_scheme_B(i_stage)*(  LNS_scheme_A(i_stage)*aux_drho &
  !                                              + deltat*(RHS_drho(:)*rmass_inverse_acoustic_DG(:))) ! Y^{n+1} = Y^{n+1} + b_i*U_{i}
  !LNS_dE                 =   LNS_dE &
  !                         + LNS_scheme_B(i_stage)*(  LNS_scheme_A(i_stage)*aux_dE &
  !                                              + deltat*(RHS_dE(:)*rmass_inverse_acoustic_DG(:)))
  !do i_aux=1,NDIM
  !  LNS_rho0dv(i_aux,:) =   LNS_rho0dv(i_aux,:) &
  !                        + LNS_scheme_B(i_stage)*(  LNS_scheme_A(i_stage)*aux_rho0dv(i_aux,:) &
  !                                             + deltat*(RHS_rho0dv(i_aux,:)*rmass_inverse_acoustic_DG(:)))
  !enddo
  
#if 0
  ! PML, iterate ADEs.
  if(PML_BOUNDARY_CONDITIONS) then
    LNS_PML_aux(:,:,:) = LNS_scheme_A(i_stage)*LNS_PML_aux(:,:,:) + deltat*LNS_PML_RHS(:,:,:)
    LNS_PML(:,:,:)     = LNS_PML(:,:,:) + LNS_scheme_B(i_stage)*LNS_PML_aux(:,:,:)
    LNS_PML(1,5:6,:)   = ZEROcr ! for q=rho', G=0, thus no convolution, thus no need for aux vars. hack, TODO something better
  endif
#endif
  
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
  endif ! Endif on check_linearHypothesis.
  
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
  
  ! TEST POSTERIORI DAMPING
  !if(.false.) then
  !if(ABC_STRETCH) then
  !  if(it==1 .and. i_stage==1 .and. myrank==0) then
  !      write(*,*) "********************************"
  !      write(*,*) "*           WARNING            *"
  !      write(*,*) "********************************"
  !      write(*,*) "* A posteriori damping of the  *"
  !      write(*,*) "* solution activated. Solution *"
  !      write(*,*) "* is being damped in the       *"
  !      write(*,*) "* absorbing buffers.           *"
  !      write(*,*) "********************************"
  !  endif
  !  call damp_solution_LNS(LNS_drho, LNS_rho0dv, LNS_dE) ! See 'prepare_stretching.f90'.
  !endif
  !endif
end subroutine compute_forces_acoustic_LNS_main










! ------------------------------------------------------------ !
! damp_solution_LNS                                            !
! ------------------------------------------------------------ !
! TODO: Description.

subroutine damp_solution_LNS(drho, rho0dv, dE)

  use specfem_par, only: nspec, coord, ibool_DG, nglob_DG, &
        ibool_before_perio,&
        ABC_STRETCH_TOP_LBUF, ABC_STRETCH_LEFT_LBUF, ABC_STRETCH_BOTTOM_LBUF, ABC_STRETCH_RIGHT_LBUF,&
        !ABC_STRETCH_LEFT, ABC_STRETCH_RIGHT, ABC_STRETCH_TOP, ABC_STRETCH_BOTTOM,&
        ispec_is_acoustic_DG,stretching_buffer,&
        mesh_xmin, mesh_xmax, mesh_zmin, mesh_zmax
  
  use constants,only: CUSTOM_REAL,NGLLX,NGLLZ,NDIM

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
            !if(stretching_buffer(ibool_before_perio(i,j,ispec))>0) then
              sigma = ONE ! In case something bad happens.
              ! Load damping value.
              if(ibits(stretching_buffer(ibool_before_perio(i,j,ispec)),0,1)==1) then
              !if(ABC_STRETCH_TOP) then
                r_l = (z - mesh_zmax)/ABC_STRETCH_TOP_LBUF + ONE
                if(r_l>ZERO .and. r_l<=ONE) call damp_function(r_l, sigma)
              endif
              if(ibits(stretching_buffer(ibool_before_perio(i,j,ispec)),1,1)==1) then
              !if(ABC_STRETCH_LEFT) then
                r_l = ONE - (x - mesh_xmin)/ABC_STRETCH_LEFT_LBUF
                if(r_l>ZERO .and. r_l<=ONE) call damp_function(r_l, sigma)
              endif
              if(ibits(stretching_buffer(ibool_before_perio(i,j,ispec)),2,1)==1) then
              !if(ABC_STRETCH_BOTTOM) then
                r_l = ONE - (z - mesh_zmin)/ABC_STRETCH_BOTTOM_LBUF
                if(r_l>ZERO .and. r_l<=ONE) call damp_function(r_l, sigma)
              endif
              if(ibits(stretching_buffer(ibool_before_perio(i,j,ispec)),3,1)==1) then
              !if(ABC_STRETCH_RIGHT) then
                r_l = (x - mesh_xmax)/ABC_STRETCH_RIGHT_LBUF + ONE
                if(r_l>ZERO .and. r_l<=ONE) call damp_function(r_l, sigma)
              endif
              ! Quiet state is zero.
              ! Damp perturbation.
              drho(iglob)     = sigma*drho(iglob)
              rho0dv(:,iglob) = sigma*rho0dv(:,iglob)
              dE(iglob)       = sigma*dE(iglob)
            !endif
          enddo
        enddo
      endif
    endif
  enddo
end subroutine damp_solution_LNS









#if 0
subroutine LNS_PML_buildRHS(idQ, idR, iglobPML, beta, q)
  ! TODO: select variables to use.
  use constants
  use specfem_par
  use specfem_par_LNS
  implicit none  
  ! Input/Output.
  integer, intent(in) :: idQ, idR, iglobPML
  real(kind=CUSTOM_REAL), intent(in) :: beta
  real(kind=CUSTOM_REAL), intent(in) :: q
  ! Local.
  ! N./A.
  !dt(auxvar) = beta * auxvar - q
  LNS_PML_RHS(idQ, idR, iglobPML) = beta*LNS_PML(idQ, idR, iglobPML) - q
end subroutine LNS_PML_buildRHS

subroutine LNS_PML_updateD0(d0cntrb, q, idQ, a1, a0, boa, b, Jac_L, iglobPML)
  ! TODO: select variables to use.
  use constants
  use specfem_par
  use specfem_par_LNS
  implicit none  
  ! Input/Output.
  real(kind=CUSTOM_REAL), intent(inout) :: d0cntrb
  real(kind=CUSTOM_REAL), intent(in) :: q, a1, a0, Jac_L
  real(kind=CUSTOM_REAL), dimension(NDIM), intent(in) :: boa, b
  integer, intent(in) :: idQ, iglobPML
  ! Local.
  ! N./A.
  d0cntrb = a1 * d0cntrb
  if(idQ==1) then
    ! Nothing to add on G for rho' (idQ==1)
    d0cntrb = d0cntrb - (   a0*q &
                          + b(1)  *LNS_PML_RHS(1,1,iglobPML) & ! YU1
                          + b(2)  *LNS_PML_RHS(1,2,iglobPML) & ! YU2
                        )*Jac_L
  else
    d0cntrb = d0cntrb - (   a0*q &
                          + boa(1)*LNS_PML_RHS(idQ,5,iglobPML) & ! G1
                          + boa(2)*LNS_PML_RHS(idQ,6,iglobPML) & ! G2
                          + b(1)  *LNS_PML_RHS(idQ,1,iglobPML) & ! YU1
                          + b(2)  *LNS_PML_RHS(idQ,2,iglobPML) & ! YU2
                        )*Jac_L
  endif
end subroutine LNS_PML_updateD0
#endif






! ------------------------------------------------------------ !
! initial_state_LNS                                            !
! ------------------------------------------------------------ !
! Computes initial state.
! 
! This routine is called by prepare_timerun(), well before
! iterate_time() and thus compute_forces_acoustic_LNS_main().

subroutine initial_state_LNS()
  use constants, only: CUSTOM_REAL, NGLLX, NGLLZ, TINYVAL
  use specfem_par, only: ibool_DG, nspec, &
                         gravityext, muext, kappa_DG, tau_epsilon, tau_sigma, &!etaext, &
                         NPROC, buffer_DG_gamma_P, gammaext_DG, ninterface_acoustic!, &
                         !PML_BOUNDARY_CONDITIONS !, ispec_is_acoustic_DG, nspec
  use specfem_par_LNS, only: LNS_E0, LNS_p0, LNS_rho0, LNS_v0, LNS_T0, LNS_mu, &
                             LNS_eta, LNS_kappa, sigma_dv, LNS_dummy_1d, LNS_dummy_2d, &
                             LNS_viscous, sigma_v_0, nabla_v0, buffer_LNS_nabla_dT, &
                             buffer_LNS_sigma_dv, LNS_switch_gradient, LNS_g,LNS_c0

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
  
  ! Not really optimal, but less invasive.
  if(allocated(gravityext)) deallocate(gravityext)
  if(allocated(muext)) deallocate(muext)
  !if(allocated(etaext)) deallocate(etaext) ! Somehow this now (at least as of commit 5d9152885007a8c8dfc42a60f203ec12afe43b1e, maybe earlier) crashes the code when using external models. I do not have the time to look more into it.
  if(allocated(kappa_DG)) deallocate(kappa_DG)
  if(allocated(tau_epsilon)) deallocate(tau_epsilon)
  if(allocated(tau_sigma)) deallocate(tau_sigma)
  
  ! Ugly patch to make it work.
  ! It seems nglob_DG includes EVERY POINT on the mesh, including those in elastic materials. See nglob_DG definition in 'createnum_slow.f90'.
  ! Thus, since DG-related variables are badly initialised in other materials (since they make little to no sense in those), arithmetic errors can arise.
  ! However, we leave it as is and use the following patches, in fear of breaking the FNS implementation.
  ! TODO: something else, maybe involving correcting the FNS implementation.
  where(LNS_rho0 < TINYVAL) LNS_rho0 = ONEcr
  where(LNS_p0 < TINYVAL) LNS_p0 = ONEcr
  where(LNS_E0 < TINYVAL) LNS_E0 = ONEcr
  
  ! Save globally c0.
  LNS_c0 = sqrt(gammaext_DG*LNS_p0/LNS_rho0)
  
  ! Initialise T_0.
  call compute_T(LNS_rho0, LNS_v0, LNS_E0, LNS_T0)
  !call compute_T(LNS_rho0, LNS_p0, LNS_T0)
  !LNS_T0=LNS_p0/(LNS_rho0*287.05684504212125397) ! Approximation assuming air molar mass is constant at 0.028964.
  
  ! Detect if viscosity exists somewhere on this CPU.
  if(     maxval(LNS_mu) > TINYVAL &
     .OR. maxval(LNS_eta) > TINYVAL &
     .OR. maxval(LNS_kappa) > TINYVAL) then
    LNS_viscous=.true.
  else
    LNS_viscous=.false.
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
  
!  write(*,*) 'min max v0x', minval(LNS_v0(1,:)), maxval(LNS_v0(1,:)) ! DEBUG
!  stop
end subroutine initial_state_LNS

#if 0
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
  real(kind=CUSTOM_REAL), parameter :: ONEcr = 1._CUSTOM_REAL
  integer :: i,j,ispec,ispec_PML
  real(kind=CUSTOM_REAL), dimension(NDIM) :: pmlk, pmld, pmla
  
  ! Safeguards.
  if(     (.not. allocated(LNS_PML_kapp)) &
     .or. (.not. allocated(LNS_PML_alpha)) &
     .or. (.not. allocated(LNS_PML_a0)) &
     .or. (.not. allocated(LNS_PML_b)) &
     .or. (.not. allocated(LNS_PML_d))) then
    write(*,*) "********************************"
    write(*,*) "*            ERROR             *"
    write(*,*) "********************************"
    write(*,*) "* Some PML coefficients arrays *"
    write(*,*) "* are not allocated but should *"
    write(*,*) "* be. Warning occured in       *"
    write(*,*) "* 'compute_forces_acoustic_LNS_calling_routine.F90',"
    write(*,*) "* allocation should happen in  *"
    write(*,*) "* 'prepare_timerun_pml.f90'.   *"
    write(*,*) "********************************"
    ! Should happen in 'prepare_timerun_pml.f90'.
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
  
  ! Value.
  do ispec=1,nspec
    if(ispec_is_PML(ispec)) then
      ispec_PML=spec_to_PML(ispec)
      do j=1,NGLLZ
        do i=1,NGLLX
          ! Note: K_x_store, K_z_store, d_x_store, d_z_store, alpha_x_store, alpha_z_store are initialised in 'pml_init.F90'.
          
          LNS_PML_alpha(1,i,j,ispec_PML) = alpha_x_store(i,j,ispec_PML)
          LNS_PML_alpha(2,i,j,ispec_PML) = alpha_z_store(i,j,ispec_PML)
          !LNS_PML_alpha(1,i,j,ispec_PML) = LNS_PML_alpha(1,i,j,ispec_PML)
          LNS_PML_alpha(2,i,j,ispec_PML) = LNS_PML_alpha(2,i,j,ispec_PML) + 0.001_CUSTOM_REAL ! TODO: check if this is necessary (must be in angles)
          
          LNS_PML_kapp(1,i,j,ispec_PML) = K_x_store(i,j,ispec_PML)
          LNS_PML_kapp(2,i,j,ispec_PML) = K_z_store(i,j,ispec_PML)
          
          ! If 1<=LNS_PML_kapp<=2 linearly, we can transform it.
          !write(*,*) 'minmax LNS_PML_kapp', minval(LNS_PML_kapp), maxval(LNS_PML_kapp) ! DEBUG
          !LNS_PML_kapp(:,i,j,ispec_PML) = ONEcr - (ONEcr - 0.2_CUSTOM_REAL) &
          !                     * (ONEcr - (ONEcr - (LNS_PML_kapp(:,i,j,ispec_PML)-ONEcr))**3.25_CUSTOM_REAL)**6._CUSTOM_REAL ! Arina's stretching
          !LNS_PML_kapp(:,i,j,ispec_PML) = ONEcr/LNS_PML_kapp(:,i,j,ispec_PML)
          
          pmlk=LNS_PML_kapp(:,i,j,ispec_PML) ! Decrease performance, but increases readability.
          
          LNS_PML_d(1,i,j,ispec_PML)=d_x_store(i,j,ispec_PML)
          LNS_PML_d(2,i,j,ispec_PML)=d_z_store(i,j,ispec_PML)
          pmld=LNS_PML_d(:,i,j,ispec_PML)
          !pmld=0._CUSTOM_REAL ! test pure stretching
          
          pmla=LNS_PML_alpha(:,i,j,ispec_PML) ! Decrease performance, but increases readability.
          
          LNS_PML_a0(i,j,ispec_PML) = pmlk(1)*pmld(2) + pmlk(2)*pmld(1)
          
          if(abs(pmla(2)-pmla(1)) < TINYVAL) then
            write(*,*) "|a2-a1| is very smol at ", coord(:, ibool_before_perio(i, j, ispec)), ": a1,a2=", pmla ! DEBUG
          else
            LNS_PML_b(1,i,j,ispec_PML) = - pmla(1)*pmld(1) * (   pmlk(2) &
                                                               - pmld(2)/(pmla(1)-pmla(2)) )
            LNS_PML_b(2,i,j,ispec_PML) = - pmla(2)*pmld(2) * (   pmlk(1) &
                                                               + pmld(1)/(pmla(1)-pmla(2)) )
          endif
          
          if(abs(coord(1, ibool_before_perio(i, j, ispec)))<TINYVAL) & ! Monitor slice at x==0. ! DEBUG
            write(*,*) coord(2, ibool_before_perio(i, j, ispec)), ": kap=", pmlk, "d=", pmld, "alpha=", pmla, & ! DEBUG
                       "a0=", LNS_PML_a0(i,j,ispec_PML), "b=",LNS_PML_b(:,i,j,ispec_PML) ! DEBUG
        enddo
      enddo
    endif
  enddo
end subroutine LNS_PML_init_coefs
#endif

! ------------------------------------------------------------ !
! background_physical_parameters                               !
! ------------------------------------------------------------ !
! Affects values of background state. May thus be used as initialiser (if time is 0), for far-field boundary conditions, or for bottom forcings.
! Note: This model-building routine builds essentially the same model as the 'boundary_condition_DG' (in 'boundary_terms_DG.f90') routine does.

subroutine background_physical_parameters(i, j, ispec, timelocal, out_rho, swComputeV, out_v, swComputeE, out_E, swComputeP, out_p)
  use constants, only: CUSTOM_REAL, TINYVAL, NDIM
  use specfem_par, only: assign_external_model, gammaext_DG, ibool_DG, pext_dg, rhoext, &
SCALE_HEIGHT, sound_velocity, surface_density, TYPE_FORCING, USE_ISOTHERMAL_MODEL, wind, windxext!, &
!ABC_STRETCH, ibool_before_perio, stretching_buffer, stretching_ya
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
  
  if(assign_external_model) then
    ! If an external model data file is given for initial conditions, read from it.
    
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
  
  !! Test for stretching
  !if(ABC_STRETCH .and. stretching_buffer(ibool_before_perio(i,j,ispec))>0) then
  !  out_v(1) = out_v(1) + (1.-stretching_ya(1, ibool_before_perio(i,j,ispec)))*340.*0.1
  !endif
  
  ! Set energy based on pressure.
  if(swComputeE) then
    call compute_E_i(out_rho, out_v, out_p, out_E, iglob)
  endif
end subroutine background_physical_parameters

! ------------------------------------------------------------ !
! set_fluid_properties                                         !
! ------------------------------------------------------------ !
! Set fluid properties.
! Note: This property-setting routine sets essentially the same values as the 'boundary_condition_DG' (in 'boundary_terms_DG.f90') routine does.

subroutine set_fluid_properties(i, j, ispec)
  use constants, only: CUSTOM_REAL, TINYVAL, NDIM, FOUR_THIRDS
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
!  use constants, only: CUSTOM_REAL, NGLLX, NGLLZ, TINYVAL, NDIM
!  use specfem_par, only: gammax, gammaz, hprime_xx, hprime_zz, hprimewgll_xx, hprimewgll_zz,&
!ibool_DG, ispec_is_acoustic_DG, jacobian, link_iface_ijispec, neighbor_dg_iface, nglob_DG, nspec, nx_iface,&
!nz_iface, rmass_inverse_acoustic_DG, weight_iface, wxgll, wzgll, xix, xiz
!  use specfem_par_LNS, only: LNS_dE, LNS_drho, LNS_dummy_1d, LNS_dummy_2d, LNS_rho0dv, LNS_v0
  ! TODO: select variables to use.
  use constants!, only: CUSTOM_REAL,NGLLX,NGLLZ,NDIM
  use specfem_par!, only: nglob_DG,nspec, ispec_is_acoustic_DG,&
                         !xix,xiz,gammax,gammaz,jacobian, &
                         !hprimewgll_xx, &
                         !hprimewgll_zz,wxgll,wzgll, &
                         !ibool_DG, &
                         !it,potential_dphi_dx_DG, potential_dphi_dz_DG, ibool, &
                         !DIR_RIGHT, DIR_LEFT, DIR_UP, DIR_DOWN, &
                         !myrank, &
                         !i_stage, p_DG_init, gammaext_DG, muext, etaext, kappa_DG,tau_epsilon, tau_sigma, &
                         !!rhovx_init, rhovz_init, E_init, &
                         !rho_init, &
                         !CONSTRAIN_HYDROSTATIC, TYPE_SOURCE_DG, &
                         !link_iface_ijispec, nx_iface, nz_iface, weight_iface, neighbor_DG_iface,&
                         !!mesh_xmin, mesh_xmax, mesh_zmin, mesh_zmax,&
                         !!coord, &
                         !ibool_before_perio,stretching_buffer!,c_V
  use specfem_par_LNS
  
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
  real(kind=CUSTOM_REAL) :: halfWeight!, flux_n, flux_x, flux_z,  !nx, nz, !rho_DG_P, rhovx_DG_P, rhovz_DG_P, &
        !E_DG_P&!, p_DG_P, &
        !Tx_DG_P, Tz_DG_P, Vxx_DG_P, Vzz_DG_P, Vzx_DG_P, Vxz_DG_P, &
        !timelocal,gamma_P
  logical :: exact_interface_flux
  !integer, dimension(nglob_DG) :: MPI_iglob
  integer, dimension(3) :: neighbor
  integer :: iface1, iface, iface1_neighbor, iface_neighbor, ispec_neighbor ,SPCDM
  real(kind=CUSTOM_REAL) :: dXi_dX_L,dXi_dZ_L,dEta_dX_L,dEta_dZ_L,Jac_L ! Jacobian matrix and determinant
  !real(kind=CUSTOM_REAL) :: temp_SFx, temp_SFz, temp_TFxx, temp_TFzx, temp_TFxz, temp_TFzz
  !real(kind=CUSTOM_REAL), dimension(nglob_DG) :: &
  !       rho_DG, rhovx_DG, rhovz_DG, E_DG, veloc_x_DG, veloc_z_DG, T, &
  !       nabla_SF(1,:), nabla_SF(NDIM,:), nabla_TF(1,1,:), nabla_TF(NDIM,NDIM,:), nabla_TF(1,NDIM,:), nabla_TF(NDIM,1,:)
  
  !real(kind=CUSTOM_REAL), dimension(NDIM,nglob_DG) :: locgrad_SF
  !real(kind=CUSTOM_REAL), dimension(NDIM,NDIM,nglob_DG) :: locgrad_TF
  
  !real(kind=CUSTOM_REAL) :: dTF_dXi(1), dTF_dEta(1), dTF_dXi(NDIM), dTF_dEta(NDIM), dSF_dXi, dSF_dEta ! Derivatives in \Lambda.
  real(kind=CUSTOM_REAL) :: dSF_dXi, dSF_dEta ! Derivatives in \Lambda.
  real(kind=CUSTOM_REAL), dimension(NDIM) :: dTF_dXi, dTF_dEta ! Derivatives in \Lambda.
  !real(kind=CUSTOM_REAL) :: dTFx_dx, dTFx_dz, dTFz_dx, dTFz_dz, dSF_dx, dSF_dz
  real(kind=CUSTOM_REAL) :: wxl, wzl ! Quadrature weights.
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: dSF_dX_tmp1, dSF_dX_tmp2, &
        dSF_dZ_tmp1, dSF_dZ_tmp2, dTF1_dX_tmp1, dTF1_dX_tmp2, &
        dTF1_dZ_tmp1, dTF1_dZ_tmp2, dTF2_dX_tmp1, dTF2_dX_tmp2, dTF2_dZ_tmp1, dTF2_dZ_tmp2
  !real(kind=CUSTOM_REAL) :: vx_init, vz_init
  !real(kind=CUSTOM_REAL) :: ya_x_l, ya_z_l
  real(kind=CUSTOM_REAL), dimension(NDIM) :: Ya_SABC ! Stretching absorbing boundary conditions.
  
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
  
  ! PASSING OF TF IS OKAY
#if 0
  if(timelocal>ZEROcr) then
    do ispec = 1, nspec; do j = 1, NGLLZ; do i = 1, NGLLX ! V1: get on all GLL points
      if(abs(TF(1,ibool_DG(i,j,ispec)) - (  aa*coord(1,ibool_before_perio(i,j,ispec))**bb &
                            + cc*coord(2,ibool_before_perio(i,j,ispec))**dd))>TINYVAL) then
        write(*,*) coord(:,ibool_before_perio(i,j,ispec)), TF(:,ibool_DG(i,j,ispec))
        stop "passing was wrong"
      endif
      if(abs(TF(2,ibool_DG(i,j,ispec)) - (   ee*coord(1,ibool_before_perio(i,j,ispec))**ff &
                             + gg*coord(2,ibool_before_perio(i,j,ispec))**hh))>TINYVAL) then
        stop "passing was wrong"
      endif
    enddo; enddo; enddo
  endif
#endif
  
  if(swMETHOD) then
    ! "Desintegrate."
    !veloc_x_DG = rhovx_DG/rho_DG
    !veloc_z_DG = rhovz_DG/rho_DG
    !T = (E_DG/rho_DG - HALFcr*(veloc_x_DG**2 + veloc_z_DG**2))/c_V
    n_out = ZEROcr
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
          iglob     = ibool_DG(i,j,ispec)
          Jac_L     = jacobian(i,j,ispec)
          dXi_dX_L  = xix(i,j,ispec)
          dXi_dZ_L  = xiz(i,j,ispec)
          dEta_dX_L = gammax(i,j,ispec)
          dEta_dZ_L = gammaz(i,j,ispec)
          
          if(ABC_STRETCH .and. stretching_buffer(ibool_before_perio(i,j,ispec))>0) then
            ! See beginning of subroutine compute_forces_acoustic_LNS for detailed explanations.
            !iglob_unique = ibool_before_perio(i, j, ispec)
            Ya_SABC  = stretching_ya(:, ibool_before_perio(i, j, ispec))
            !ya_z_l  = stretching_ya(2, iglob_unique)
            ! If you change anything, remember to do it also in the 'compute_forces_acoustic_LNS' subroutine.
            ! Multiply x component by ya_x.
            dXi_dX_L    = Ya_SABC(1) * dXi_dX_L
            ! Multiply z component by ya_z.
            dXi_dZ_L    = Ya_SABC(2) * dXi_dZ_L
            ! Multiply x component by ya_x.
            dEta_dX_L   = Ya_SABC(1) * dEta_dX_L
            ! Multiply z component by ya_z.
            dEta_dZ_L   = Ya_SABC(2) * dEta_dZ_L
            ! TODO: something on jacobian??
            !Jac_L = ya_x_l*ya_z_l*Jac_L
          endif
          
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
              dSF_dX_tmp1(i,j)  = wzl * Jac_L * (dXi_dX_L * SF(iglob))
              dSF_dZ_tmp1(i,j)  = wzl * Jac_L * (dXi_dZ_L * SF(iglob))
              dSF_dX_tmp2(i,j)  = wxl * Jac_L * (dEta_dX_L * SF(iglob))
              dSF_dZ_tmp2(i,j)  = wxl * Jac_L * (dEta_dZ_L * SF(iglob))
            endif
            if(swTF) then
              dTF1_dX_tmp1(i,j) = wzl * Jac_L * (dXi_dX_L * TF(1,iglob))
              dTF1_dZ_tmp1(i,j) = wzl * Jac_L * (dXi_dZ_L * TF(1,iglob))
              dTF2_dX_tmp1(i,j) = wzl * Jac_L * (dXi_dX_L * TF(NDIM,iglob))
              dTF2_dZ_tmp1(i,j) = wzl * Jac_L * (dXi_dZ_L * TF(NDIM,iglob))
              dTF1_dX_tmp2(i,j) = wxl * Jac_L * (dEta_dX_L * TF(1,iglob))
              dTF1_dZ_tmp2(i,j) = wxl * Jac_L * (dEta_dZ_L * TF(1,iglob))
              dTF2_dX_tmp2(i,j) = wxl * Jac_L * (dEta_dX_L * TF(NDIM,iglob))
              dTF2_dZ_tmp2(i,j) = wxl * Jac_L * (dEta_dZ_L * TF(NDIM,iglob))
            endif
            !else
            !  vx_init = rhovx_init(iglob)/rho_init(iglob)
            !  vz_init = rhovz_init(iglob)/rho_init(iglob)
            !  
            !  dSF_dX_tmp1(i,j)  = wzl * Jac_L * (dXi_dX_L * (SF(iglob) - T_init(iglob)))
            !  dSF_dZ_tmp1(i,j)  = wzl * Jac_L * (dXi_dZ_L * (SF(iglob) - T_init(iglob)))
            !  dTF1_dX_tmp1(i,j) = wzl * Jac_L * (dXi_dX_L * (TF(1,iglob) - vx_init))
            !  dTF1_dZ_tmp1(i,j) = wzl * Jac_L * (dXi_dZ_L * (TF(1,iglob))) ! Some hypothesis on initial velocity is used here.
            !  dTF2_dX_tmp1(i,j) = wzl * Jac_L * (dXi_dX_L * (TF(NDIM,iglob) - vz_init))
            !  dTF2_dZ_tmp1(i,j) = wzl * Jac_L * (dXi_dZ_L * (TF(NDIM,iglob) - vz_init))
            !  
            !  dSF_dX_tmp2(i,j)  = wxl * Jac_L * (dEta_dX_L * (SF(iglob) - T_init(iglob)))
            !  dSF_dZ_tmp2(i,j)  = wxl * Jac_L * (dEta_dZ_L * (SF(iglob) - T_init(iglob)))
            !  dTF1_dX_tmp2(i,j) = wxl * Jac_L * (dEta_dX_L * (TF(1,iglob) - vx_init))
            !  dTF1_dZ_tmp2(i,j) = wxl * Jac_L * (dEta_dZ_L * (TF(1,iglob))) ! Some hypothesis on initial velocity is used here.
            !  dTF2_dX_tmp2(i,j) = wxl * Jac_L * (dEta_dX_L * (TF(NDIM,iglob) - vz_init))
            !  dTF2_dZ_tmp2(i,j) = wxl * Jac_L * (dEta_dZ_L * (TF(NDIM,iglob) - vz_init))
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
              dSF_dXi       = ZEROcr
              dSF_dEta      = ZEROcr
            endif
            if(swTF) then
              dTF_dXi       = ZEROcr ! Vector operation.
              dTF_dEta      = ZEROcr ! Vector operation.
              !dTF_dXi(1)     = ZEROcr ! Safeguard: do scalar instead.
              !dTF_dXi(2)     = ZEROcr ! Safeguard: do scalar instead.
              !dTF_dEta(1)    = ZEROcr ! Safeguard: do scalar instead.
              !dTF_dEta(2)    = ZEROcr ! Safeguard: do scalar instead.
            endif
            
            ! Compute derivatives in unit element \Lambda.
            ! Note: we can merge the two loops because NGLLX=NGLLZ.
            do k = 1, NGLLX
              !if(.not. CONSTRAIN_HYDROSTATIC) then
              if(swSF) then
                dSF_dXi   = dSF_dXi   + SF(ibool_DG(k,j,ispec)) * real(hprime_xx(i,k), kind=CUSTOM_REAL)
                dSF_dEta  = dSF_dEta  + SF(ibool_DG(i,k,ispec)) * real(hprime_zz(j,k), kind=CUSTOM_REAL)
              endif
              if(swTF) then
                dTF_dXi(1)     = dTF_dXi(1)     + TF(1,ibool_DG(k,j,ispec)) * real(hprime_xx(i,k), kind=CUSTOM_REAL)
                dTF_dEta(1)    = dTF_dEta(1)    + TF(1,ibool_DG(i,k,ispec)) * real(hprime_zz(j,k), kind=CUSTOM_REAL)
                dTF_dXi(NDIM)  = dTF_dXi(NDIM)  + TF(NDIM,ibool_DG(k,j,ispec)) * real(hprime_xx(i,k), kind=CUSTOM_REAL)
                dTF_dEta(NDIM) = dTF_dEta(NDIM) + TF(NDIM,ibool_DG(i,k,ispec)) * real(hprime_zz(j,k), kind=CUSTOM_REAL)
              endif
              !else
              !  vx_init = rhovx_init(ibool_DG(k,j,ispec))/rho_init(ibool_DG(k,j,ispec))
              !  dTF_dXi(1) = dTF_dXi(1) + (TF(1,ibool_DG(k,j,ispec)) - vx_init) &
              !                      * real(hprime_xx(i,k), kind=CUSTOM_REAL)
              !  vx_init = rhovx_init(ibool_DG(i,k,ispec))/rho_init(ibool_DG(i,k,ispec))       
              !  dTF_dEta(1) = dTF_dEta(1) + (TF(1,ibool_DG(i,k,ispec)) - vx_init) &
              !                            * real(hprime_zz(j,k), kind=CUSTOM_REAL)
              !  vz_init = rhovz_init(ibool_DG(k,j,ispec))/rho_init(ibool_DG(k,j,ispec))
              !  dTF_dXi(NDIM) = dTF_dXi(NDIM) + (TF(NDIM,ibool_DG(k,j,ispec)) - vz_init) &
              !                      * real(hprime_xx(i,k), kind=CUSTOM_REAL)
              !  vz_init = rhovz_init(ibool_DG(i,k,ispec))/rho_init(ibool_DG(i,k,ispec))
              !  dTF_dEta(NDIM) = dTF_dEta(NDIM) + (TF(NDIM,ibool_DG(i,k,ispec)) - vz_init) &
              !                            * real(hprime_zz(j,k), kind=CUSTOM_REAL)
              !  dSF_dXi = dSF_dXi + (SF(ibool_DG(k,j,ispec)) - T_init(ibool_DG(k,j,ispec))) &
              !                    * real(hprime_xx(i,k), kind=CUSTOM_REAL)
              !  dSF_dEta = dSF_dEta + (SF(ibool_DG(i,k,ispec)) - T_init(ibool_DG(i,k,ispec))) &
              !                           * real(hprime_zz(j,k), kind=CUSTOM_REAL)
              !endif
#if 0
              if(abs(coord(1,ibool_before_perio(i,j,ispec))-2.)<TINYVAL .and. &
                 abs(coord(2,ibool_before_perio(i,j,ispec))-3.)<TINYVAL) then
                  write(*,*) 'k, dTF_dXi(2)', k, dTF_dXi(2), &
                             'added', TF(2,ibool_DG(k,j,ispec)), real(hprime_xx(i,k), kind=CUSTOM_REAL)
              endif
#endif
            enddo ! Enddo on k.

            ! Compute derivatives in element \Omega using the chain rule.
            if(swSF) then
              !dSF_dx   = dSF_dXi * dXi_dX_L + dSF_dEta * dEta_dX_L
              !dSF_dz   = dSF_dXi * dXi_dZ_L + dSF_dEta * dEta_dZ_L
              !temp_SFx = dSF_dx * Jac_L
              !temp_SFz = dSF_dz * Jac_L
              !nabla_SF(1,iglob) = dSF_dx
              !nabla_SF(NDIM,iglob) = dSF_dz
              nabla_SF(1,iglob)    = dSF_dXi * dXi_dX_L + dSF_dEta * dEta_dX_L
              nabla_SF(NDIM,iglob) = dSF_dXi * dXi_dZ_L + dSF_dEta * dEta_dZ_L
            endif
            if(swTF) then
              !dTFx_dx   = dTF_dXi(1) * dXi_dX_L + dTF_dEta(1) * dEta_dX_L
              !dTFx_dz   = dTF_dXi(1) * dXi_dZ_L + dTF_dEta(1) * dEta_dZ_L
              !dTFz_dx   = dTF_dXi(NDIM) * dXi_dX_L + dTF_dEta(NDIM) * dEta_dX_L
              !dTFz_dz   = dTF_dXi(NDIM) * dXi_dZ_L + dTF_dEta(NDIM) * dEta_dZ_L
              !temp_TFxx = dTFx_dx * Jac_L
              !temp_TFxz = dTFx_dz * Jac_L
              !temp_TFzx = dTFz_dx * Jac_L
              !temp_TFzz = dTFz_dz * Jac_L
              !nabla_TF(1,1,iglob) = dTFx_dx
              !nabla_TF(1,NDIM,iglob) = dTFx_dz
              !nabla_TF(NDIM,1,iglob) = dTFz_dx
              !nabla_TF(NDIM,NDIM,iglob) = dTFz_dz
              
              do k=1,NDIM ! Re-using index k for spatial dimension.
                nabla_TF(k,1,iglob)    = dTF_dXi(k) * dXi_dX_L + dTF_dEta(k) * dEta_dX_L
                nabla_TF(k,NDIM,iglob) = dTF_dXi(k) * dXi_dZ_L + dTF_dEta(k) * dEta_dZ_L
              enddo ! Enddo on k.
#if 0
              if(abs(coord(1,ibool_before_perio(i,j,ispec))-2.)<TINYVAL .and. &
                 abs(coord(2,ibool_before_perio(i,j,ispec))-3.)<TINYVAL) then
                write(*,*) '(x z) (dxvx dzvx dxvz dzvz)',&
                            !nabla_dv(1, 1, ibool_DG(i,j,ispec)), nabla_dv(1, 2, ibool_DG(i,j,ispec)), & ! TEST
                            nabla_TF(2, :, ibool_DG(i,j,ispec))! TEST
                write(*,*) 'components',&
                            dTF_dXi(2), dXi_dX_L, dTF_dEta(2), dEta_dX_L, & ! dTF_dXi at fault
                            dTF_dXi(2), dXi_dZ_L, dTF_dEta(2), dEta_dZ_L ! dTF_dEta at fault (for high degrees only it seems)
                write(*,*) 'th =                       ', &
                            !aa*bb*coord(1,ibool_before_perio(i,j,ispec))**(bb-1), &
                            !cc*dd*coord(2,ibool_before_perio(i,j,ispec))**(dd-1), &
                            ee*ff*coord(1,ibool_before_perio(i,j,ispec))**(ee-1), &
                            gg*hh*coord(2,ibool_before_perio(i,j,ispec))**(hh-1)
             endif
#endif
              
             ! nabla_TF(1,1,iglob) = dTF_dXi(1) * dXi_dX_L + dTF_dEta(1) * dEta_dX_L ! dTFx_dx = d[TF(1)]/dXi * dXi/dX + d[TF(1)]/dEta * dEta/dX
             ! nabla_TF(1,2,iglob) = dTF_dXi(1) * dXi_dZ_L + dTF_dEta(1) * dEta_dZ_L ! dTFx_dz = d[TF(1)]/dXi * dXi/dZ + d[TF(1)]/dEta * dEta/dZ
             ! nabla_TF(2,1,iglob) = dTF_dXi(2) * dXi_dX_L + dTF_dEta(2) * dEta_dX_L ! dTFz_dx = d[TF(2)]/dXi * dXi/dX + d[TF(2)]/dEta * dEta/dX
             ! nabla_TF(2,2,iglob) = dTF_dXi(2) * dXi_dZ_L + dTF_dEta(2) * dEta_dZ_L ! dTFz_dz = d[TF(2)]/dXi * dXi/dZ + d[TF(2)]/dEta * dEta/dZ
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
                                       (dSF_dX_tmp1(k,j) * real(hprimewgll_xx(k,i), kind=CUSTOM_REAL) + &
                                        dSF_dX_tmp2(i,k) * real(hprimewgll_zz(k,j), kind=CUSTOM_REAL))
                nabla_SF(NDIM,iglob) = nabla_SF(NDIM,iglob) - &
                                       (dSF_dZ_tmp1(k,j) * real(hprimewgll_xx(k,i), kind=CUSTOM_REAL) + &
                                        dSF_dZ_tmp2(i,k) * real(hprimewgll_zz(k,j), kind=CUSTOM_REAL))
              endif
              if(swTF) then
                nabla_TF(1,1,iglob)       = nabla_TF(1,1,iglob) - &
                                            (dTF1_dX_tmp1(k,j) * real(hprimewgll_xx(k,i), kind=CUSTOM_REAL) + &
                                             dTF1_dX_tmp2(i,k) * real(hprimewgll_zz(k,j), kind=CUSTOM_REAL))
                nabla_TF(1,NDIM,iglob)    = nabla_TF(1,NDIM,iglob) - &
                                            (dTF1_dZ_tmp1(k,j) * real(hprimewgll_xx(k,i), kind=CUSTOM_REAL) + &
                                             dTF1_dZ_tmp2(i,k) * real(hprimewgll_zz(k,j), kind=CUSTOM_REAL))
                nabla_TF(NDIM,1,iglob)    = nabla_TF(NDIM,1,iglob) - &
                                            (dTF2_dX_tmp1(k,j) * real(hprimewgll_xx(k,i), kind=CUSTOM_REAL) + &
                                             dTF2_dX_tmp2(i,k) * real(hprimewgll_zz(k,j), kind=CUSTOM_REAL))
                nabla_TF(NDIM,NDIM,iglob) = nabla_TF(NDIM,NDIM,iglob) - &
                                            (dTF2_dZ_tmp1(k,j) * real(hprimewgll_xx(k,i), kind=CUSTOM_REAL) + &
                                             dTF2_dZ_tmp2(i,k) * real(hprimewgll_zz(k,j), kind=CUSTOM_REAL))
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
            endif
            
            ! Step 2: knowing the normals' parameters, compute now the fluxes.
            iglobP = 1
            if(neighbor(1) > -1) then
              iglobP = ibool_DG(neighbor(1), neighbor(2), neighbor(3))
            endif
            exact_interface_flux = .false. ! Reset this variable to .false.: by default, the fluxes have to be computed (jump!=0). In some specific cases (assigned during the call to LNS_get_interfaces_unknowns), the flux can be exact (jump==0).
            call LNS_get_interfaces_unknowns(i, j, ispec, iface1, iface, neighbor, timelocal, & ! Point identifier (input).
                  LNS_drho(iglob), LNS_rho0dv(:,iglob), & ! Input constitutive variables, "M" side.
                  LNS_drho(iglobP), LNS_rho0dv(:,iglobP), LNS_dE(iglobP), & ! Input constitutive variables, "P" side.
                  LNS_dp(iglob), & ! Input other variable, "M" side.
                  !V_DG(:,:,iglob), T_DG(:,iglob), & ! Input derivatives, "M" side. MIGHT NEED.
                  !V_DG(:,:,iglobP), T_DG(:,iglobP), & ! Input derivatives, "M" side. MIGHT NEED.
                  n_out, & ! Normal vector (input).
                  exact_interface_flux, & ! Switch to disable jump in some cases (output).
                  LNS_dummy_1d(1), LNS_dummy_2d(:,1), LNS_dummy_1d(1), & ! Output constitutive variables.
                  !Tx_DG_P, Tz_DG_P, Vxx_DG_P, Vzz_DG_P, Vzx_DG_P, Vxz_DG_P,& ! Output derivatives. MIGHT NEED.
                  LNS_dummy_2d(:,1), LNS_dummy_1d(1), TF_P, & ! Output other variables.
                  .false., LNS_dummy_2d(:,1), LNS_dummy_1d(1:3), & ! Output other variables: viscous.
                  .true., SF_P) ! Output other variables.
            ! We know that the scalar field will always be temperature, and that the tensor field will always be the velocity, so we use the already-built LNS_get_interfaces_unknowns routine to get them. The switch is here to prevent unecessary quantites to be computed during that specific call. If other fields need to be computed, one will have to use a dedicated routine.
            
            ! For real stretching, \Sigma for each constitutive variable becomes \Ya\Sigma. It is heavy to change each and every expression where \Sigma arises. Rather, we make use of what is multiplying \Sigma.
            ! Here, in the surface integrations, easiest is the normals. But we have to do it after the call to LNS_get_interfaces_unknowns. If you change anything, remember to do it also in the 'compute_forces_acoustic_LNS' subroutine.
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
            !if(swSF) then
            !  !!flux_x = SF(iglob) + SF_P ! Once multiplied with HALFcr below, will represent flux along x of the scalar field SF.
            !  !!!if(CONSTRAIN_HYDROSTATIC) flux_x = flux_x - 2*T_init(iglob)
            !  !!flux_n = flux_x*n_out(1) ! Recall: n_out(1)=n_x.
            !  !flux_n = (SF(iglob)+SF_P)*n_out(1) ! Recall: n_out(1)=n_x.
            !  !nabla_SF(1,iglob) = nabla_SF(1,iglob) + halfWeight*flux_n
            !  !!flux_z = SF(iglob) + SF_P ! Once multiplied with HALFcr below, will represent flux along z of the scalar field SF.
            !  !!!if(CONSTRAIN_HYDROSTATIC) flux_z = flux_z - 2*T_init(iglob)
            !  !!flux_n = flux_z*n_out(NDIM) ! Recall: n_out(NDIM)=n_z.
            !  !flux_n = (SF(iglob)+SF_P)*n_out(NDIM) ! Recall: n_out(NDIM)=n_z.
            !  !nabla_SF(NDIM,iglob) = nabla_SF(NDIM,iglob) + halfWeight*flux_n
            !  do i=1,NDIM
            !    !flux_n = (SF(iglob)+SF_P)*n_out(i) ! Recall: n_out(1)=n_x.
            !    !nabla_SF(i,iglob) = nabla_SF(i,iglob) + halfWeight*flux_n
            !    nabla_SF(i,iglob) = nabla_SF(i,iglob) + halfWeight*(SF(iglob)+SF_P)*n_out(i)
            !  enddo
            !endif
            !if(swTF) then
            !  !!flux_x = TF(1,iglob) + TF_P(1) ! Once multiplied with HALFcr below, will represent flux along x of the x-component of the tensor field TF.
            !  !!!if(CONSTRAIN_HYDROSTATIC) flux_x = flux_x - 2*vx_init
            !  !!flux_n = flux_x*n_out(1) ! Recall: n_out(1)=n_x.
            !  !flux_n = (TF(1,iglob)+TF_P(1))*n_out(1) ! Recall: n_out(1)=n_x.
            !  !nabla_TF(1,1,iglob) = nabla_TF(1,1,iglob) + halfWeight*flux_n
            
            !  !flux_z = TF(1,iglob) + TF_P(1) ! Once multiplied with HALFcr below, will represent flux along z of the x-component of the tensor field TF.
            !  !!if(CONSTRAIN_HYDROSTATIC) flux_z = flux_z - 2*vx_init
            !  !flux_n = flux_z*n_out(NDIM) ! Recall: n_out(NDIM)=n_z.
            !  !nabla_TF(1,NDIM,iglob) = nabla_TF(1,NDIM,iglob) + halfWeight*flux_n
            
            !  !flux_x = TF(NDIM,iglob) + TF_P(NDIM) ! Once multiplied with HALFcr below, will represent flux along x of the z-component of the tensor field TF.
            !  !!if(CONSTRAIN_HYDROSTATIC) flux_x = flux_x - 2*vz_init
            !  !flux_n = flux_x*n_out(1) ! Recall: n_out(1)=n_x.
            !  !nabla_TF(NDIM,1,iglob) = nabla_TF(NDIM,1,iglob) + halfWeight*flux_n
            
            !  !flux_z = TF(NDIM,iglob) + TF_P(NDIM) ! Once multiplied with HALFcr below, will represent flux along z of the z-component of the tensor field TF.
            !  !!if(CONSTRAIN_HYDROSTATIC) flux_z = flux_z - 2*vz_init
            !  !flux_n = flux_z*n_out(NDIM) ! Recall: n_out(NDIM)=n_z.
            !  !nabla_TF(NDIM,NDIM,iglob) = nabla_TF(NDIM,NDIM,iglob) + halfWeight*flux_n
            !  do i=1,NDIM
            !    do j=1,NDIM
            !      nabla_TF(i,j,iglob) = nabla_TF(i,j,iglob) + halfWeight*(TF(i,iglob)+TF_P(i))*n_out(j)
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
  
#if 0
    do ispec = 1, nspec; do j = 1, NGLLZ; do i = 1, NGLLX ! V1: get on all GLL points
      !if(abs(coord(1,ibool_before_perio(i,j,ispec))-0.)<TINYVAL .and. abs(coord(2,ibool_before_perio(i,j,ispec))-2.5)<TINYVAL) then
      !if(abs(nabla_dv(1, 1, ibool_DG(i,j,ispec))-aa*bb*coord(1,ibool_before_perio(i,j,ispec))**(bb-1))>TINYVAL .or. &
      !   abs(nabla_dv(1, 2, ibool_DG(i,j,ispec))-cc*dd*coord(2,ibool_before_perio(i,j,ispec))**(dd-1))>TINYVAL) then! .or. &
      if(abs(nabla_TF(2, 1, ibool_DG(i,j,ispec))-ee*ff*coord(1,ibool_before_perio(i,j,ispec))**(ff-1))>TINYVAL .or. &
         abs(nabla_TF(2, 2, ibool_DG(i,j,ispec))-gg*hh*coord(2,ibool_before_perio(i,j,ispec))**(hh-1))>TINYVAL ) then
        write(*,*) '(x z) (dxvx dzvx dxvz dzvz)',&
                                            coord(:,ibool_before_perio(i,j,ispec)), &
                                            !nabla_dv(1, 1, ibool_DG(i,j,ispec)), nabla_dv(1, 2, ibool_DG(i,j,ispec)), & ! TEST
                                            nabla_TF(2, 1, ibool_DG(i,j,ispec)), nabla_TF(2, 2, ibool_DG(i,j,ispec)) ! TEST
        write(*,*) 'th =                       ', &
                                                  coord(:,ibool_before_perio(i,j,ispec)), &
                                                  !aa*bb*coord(1,ibool_before_perio(i,j,ispec))**(bb-1), &
                                                  !cc*dd*coord(2,ibool_before_perio(i,j,ispec))**(dd-1), &
                                                  ee*ff*coord(1,ibool_before_perio(i,j,ispec))**(ff-1.), &
                                                  gg*hh*coord(2,ibool_before_perio(i,j,ispec))**(hh-1.)
        endif
      !endif
    enddo; enddo; enddo
#endif
end subroutine compute_gradient_TFSF




















#if 0
! UNUSED
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
  !out_T = (in_E/in_rho - HALFcr*(in_v(1,:)**2+in_v(NDIM,:)**2))/c_V
  out_T = (in_E/in_rho - HALFcr*norm2(in_v))/c_V
end subroutine compute_T
!subroutine compute_T(in_rho, in_p, out_T)
!  use constants, only: CUSTOM_REAL, NDIM
!  use specfem_par, only: nglob_DG
!  use specfem_par_LNS, only: R_ADIAB
!  implicit none
!  ! Input/Output.
!  real(kind=CUSTOM_REAL), dimension(nglob_DG), intent(in) :: in_rho, in_p
!  real(kind=CUSTOM_REAL), dimension(nglob_DG), intent(out) :: out_T
!  out_T = in_p / (R_ADIAB * in_rho)
!end subroutine compute_T
! ------------------------------------------------------------ !
! compute_T_i                                                  !
! ------------------------------------------------------------ !
! Same as compute_T, but point by point (unvectorised).
! Unused as of 190124.
!subroutine compute_T_i(in_rho, in_v, in_E, out_T)
!  use constants, only: CUSTOM_REAL, NDIM
!  use specfem_par, only: c_V
!  use specfem_par_LNS, only: norm2r1
!  implicit none
!  ! Input/Output.
!  real(kind=CUSTOM_REAL), intent(in) :: in_rho, in_E
!  real(kind=CUSTOM_REAL), dimension(NDIM), intent(in) :: in_v
!  real(kind=CUSTOM_REAL), intent(out) :: out_T
!  ! Local.
!  real(kind=CUSTOM_REAL), parameter :: HALFcr = 0.5_CUSTOM_REAL
!  !out_T=(in_E/in_rho - HALFcr*(in_v(1)**2+in_v(NDIM)**2))/c_V
!  !out_T=(in_E/in_rho - HALFcr*minval(norm2(reshape(in_v,(/2,1/)))))/c_V ! Could not make our "norm2" function work both with rank 1 arrays and rank 2 arrays, so used a trick: reshape the 2-sized vector (rank 1) to a (2,1) matrix (rank 2), send it to norm2, retrieve a 1-sized vector (rank 1), use minval to send it back to a scalar value (rank 0) as needed. It proved very unefficient in terms of computation time, so we defined another dedicated "norm2" function.
!  out_T=(in_E/in_rho - HALFcr*norm2r1(in_v))/c_V ! Performance seems comparable to explicit formulation above.
!end subroutine compute_T_i
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
  !out_dT =   (in_E/in_rho - HALFcr*(in_v(1,:)**2+in_v(NDIM,:)**2))/c_V &
  out_dT =   (in_E/in_rho - HALFcr*norm2(in_v))/c_V &
           - LNS_T0
end subroutine compute_dT
!subroutine compute_dT(in_rho, in_p, out_dT)
!  use constants, only: CUSTOM_REAL, NDIM
!  use specfem_par, only: nglob_DG!,c_V
!  use specfem_par_LNS, only: LNS_T0,R_ADIAB
!  implicit none
!  ! Input/Output.
!  real(kind=CUSTOM_REAL), dimension(nglob_DG), intent(in) :: in_rho, in_p
!  real(kind=CUSTOM_REAL), dimension(nglob_DG), intent(out) :: out_dT
!  out_dT = in_p / (R_ADIAB * in_rho) - LNS_T0
!end subroutine compute_dT
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
  out_T =   (in_E/in_rho - HALFcr*norm2r1(in_v))/c_V &
          - LNS_T0(iglob)
end subroutine compute_dT_i
!subroutine compute_dT_i(in_rho, in_p, out_dTi, iglob)
!  use constants, only: CUSTOM_REAL, NDIM
!  !use specfem_par, only: c_V
!  use specfem_par_LNS, only: LNS_T0, R_ADIAB
!  implicit none
!  ! Input/Output.
!  real(kind=CUSTOM_REAL), intent(in) :: in_rho, in_p
!  real(kind=CUSTOM_REAL), intent(out) :: out_dTi
!  integer, intent(in) :: iglob
!  out_dTi = in_p / (R_ADIAB * in_rho) - LNS_T0(iglob)
!end subroutine compute_dT_i
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
! stressBuilder_addInviscidFluid                               !
! ------------------------------------------------------------ !
! Routine building on top of a stress tensor, adding the inviscid part of the fluid stress.
subroutine stressBuilder_addInviscidFluid(rho, v1, v2, pressure, out_sigma)
  use constants, only: CUSTOM_REAL, NDIM
  
  implicit none
  
  ! Input/Output.
  real(kind=CUSTOM_REAL), intent(in) :: rho, pressure
  real(kind=CUSTOM_REAL), dimension(NDIM), intent(in) :: v1, v2
  real(kind=CUSTOM_REAL), dimension(NDIM, NDIM), intent(out) :: out_sigma
  
  ! Local.
  integer :: i, j
  
  !out_sigma(1, 1) = out_sigma(1, 1) + rho*v1(1)*v2(1) + pressure
  !out_sigma(1, 2) = out_sigma(1, 2) + rho*v1(1)*v2(2)
  !out_sigma(2, 1) = out_sigma(2, 1) + rho*v1(2)*v2(1)
  !out_sigma(2, 2) = out_sigma(2, 2) + rho*v1(2)*v2(2) + pressure
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
! Routine building on top of a stress tensor, adding the viscous part of the fluid stress.
subroutine stressBuilder_addViscousFluid(viscous_tensor_local, out_sigma)
  use constants, only: CUSTOM_REAL, NDIM
  use specfem_par_lns, only: NVALSIGMA
  
  implicit none
  
  ! Input/Output.
  real(kind=CUSTOM_REAL), dimension(NVALSIGMA), intent(in) :: viscous_tensor_local
  real(kind=CUSTOM_REAL), dimension(NDIM, NDIM), intent(out) :: out_sigma
  
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










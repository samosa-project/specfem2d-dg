! ------------------------------------------------------------ !
! compute_forces_acoustic_LNS_main                              !
! ------------------------------------------------------------ !
! TODO: Description.

subroutine compute_forces_acoustic_LNS_main()

  use constants, only: &
    CUSTOM_REAL, &
    rk4a_d, rk4b_d, rk4c_d, &
    ls33rk_a, ls33rk_b, ls33rk_c, &
    HALF, ONE, NGLLX, NGLLZ
  use specfem_par, only: &
    deltat, nglob_DG,stage_time_scheme,&
    gravityext, muext, etaext, kappa_DG,&
    tau_epsilon, tau_sigma, & ! Temporary hack until boundary_condition_DG is modified.
    ispec_is_acoustic_coupling_ac,&
    rmass_inverse_acoustic_DG,&
    assign_external_model, any_acoustic_DG, only_DG_acoustic, &
    c_V,gammaext_DG,&
    it, i_stage,&
    myrank, nspec, nproc, ibool_DG, ibool_before_perio, &
    ninterface_acoustic, &
    time_stepping_scheme, &
    CONSTRAIN_HYDROSTATIC, coord, &
    T_DG, V_DG ! Use those to store respectively \nabla T' and \nabla v'.
  use specfem_par_LNS, only: &
    LNS_g, LNS_mu, LNS_eta, LNS_kappa,&
    LNS_rho0, LNS_v0, LNS_E0,&
    LNS_drho, LNS_rho0dv, LNS_dE,&
    aux_drho, aux_rho0dvx, aux_rho0dvz, aux_dE,&
    RHS_drho, RHS_rho0dvx, RHS_rho0dvz, RHS_dE,&
    nabla_v0, LNS_p0, LNS_T0, sigma_v_0,&
    LNS_dm, LNS_dp, LNS_dT,&
    buffer_LNS_drho_P, buffer_LNS_rho0dvx_P, buffer_LNS_rho0dvz_P, buffer_LNS_dE_P,&
    LNS_VERBOSE,&
    LNS_dummy_1d, LNS_dummy_2d!, LNS_dummy_3d

  implicit none

  ! Local.
  real(kind=CUSTOM_REAL), parameter :: ZEROcr = 0._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: threshold = 0.0000001_CUSTOM_REAL
  real(kind=CUSTOM_REAL) :: timelocal
  real(kind=CUSTOM_REAL), dimension(stage_time_scheme) :: scheme_A, scheme_B, scheme_C
  real(kind=CUSTOM_REAL), dimension(2,nglob_DG) :: LNS_dv
  real(kind=CUSTOM_REAL), dimension(2,2,nglob_DG) :: nabla_dv
  real(kind=CUSTOM_REAL), dimension(3,nglob_DG) :: sigma_dv
  real(kind=CUSTOM_REAL), dimension(2,nglob_DG) :: nabla_dT
  integer :: i,j,ispec
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
    
    ! Zero registers.
    aux_drho    = ZEROcr
    aux_rho0dvx = ZEROcr
    aux_rho0dvz = ZEROcr
    aux_dE      = ZEROcr
    RHS_drho    = ZEROcr
    RHS_rho0dvx = ZEROcr
    RHS_rho0dvz = ZEROcr
    RHS_dE      = ZEROcr
    T_DG        = ZEROcr
    V_DG        = ZEROcr
    
    ! Physical parameters.
    if(.not. assign_external_model) then
      deallocate(gravityext, muext, etaext, kappa_DG)
      allocate(gravityext(NGLLX, NGLLZ, nspec), &
               etaext(NGLLX, NGLLZ, nspec), &
               muext(NGLLX, NGLLZ, nspec), &
               kappa_DG(NGLLX, NGLLZ, nspec))
      deallocate(tau_epsilon, tau_sigma) ! Temporary hack until boundary_condition_DG is modified.
      allocate(tau_epsilon(NGLLX, NGLLZ, nspec), tau_sigma(NGLLX, NGLLZ, nspec)) ! Temporary hack until boundary_condition_DG is modified.
    endif
    
    ! Prepare MPI buffers.
    call prepare_MPI_DG()
    
    ! Allocate acoustic coupling array.
    allocate(ispec_is_acoustic_coupling_ac(nglob_DG))
    ispec_is_acoustic_coupling_ac = -1
    if(.not. only_DG_acoustic) then
      call find_DG_acoustic_coupling()
    endif
    
    ! Initialise state registers.
    ! Note: since constitutive variables are perturbations, they are zero at start.
    LNS_drho   = ZEROcr
    LNS_rho0dv = ZEROcr
    LNS_dE     = ZEROcr

    ! Save initial state.
    !write(*,*) "DEBUG LNS: initial_state_LNS started."
    call initial_state_LNS(LNS_rho0, LNS_v0, LNS_E0)
    !write(*,*) "DEBUG LNS: initial_state_LNS finished."
    
    ! Send the (:,:,:) classic registers depending on (i,j,ispec) into one (:) register depending on (iglob), in order to be able to use vector calculus.
    ! This is only done here because we use the subroutine "boundary_condition_DG" in "initial_state_LNS", which needs the classic registers.
    allocate(LNS_g(nglob_DG), &
             LNS_eta(nglob_DG), &
             LNS_mu(nglob_DG), &
             LNS_kappa(nglob_DG))
    ! TODO: Somehow send the (:,:,:) classic registers depending on (i,j,ispec) into one (:) register depending on (iglob), in order to be able to use vector calculus.
    ! Not really optimal, but less invasive.
    deallocate(gravityext, muext, etaext, kappa_DG)
    
    ! Compute p_0 and T_0.
    LNS_p0 = (gammaext_DG - ONE) &
              * (LNS_E0 - HALF*LNS_rho0*(LNS_v0(1,:)**2+LNS_v0(2,:)**2))
    LNS_T0 = (LNS_E0/LNS_rho0 - HALF*(LNS_v0(1,:)**2+LNS_v0(2,:)**2))/c_V
    where(LNS_p0 <= 0._CUSTOM_REAL) LNS_p0 = 1._CUSTOM_REAL ! When elastic-DG simulations, p_DG = 0 in elastic elements and 1/p_DG will not be properly defined. Use a hack.
    
    ! Compute \nabla v_0.
    !write(*,*) "DEBUG LNS: compute_gradient_TFSF started."
    call compute_gradient_TFSF(LNS_v0, LNS_dummy_1d, .true., .false., switch_gradient, nabla_v0, LNS_dummy_2d) ! Dummy variables are not optimal, but prevent from duplicating subroutines.
    !write(*,*) "DEBUG LNS: compute_gradient_TFSF finished."
    
    ! Compute \Sigma_v_0.
    !write(*,*) "DEBUG LNS: LNS_compute_viscous_stress_tensor started."
    call LNS_compute_viscous_stress_tensor(nabla_v0, sigma_v_0)
    !write(*,*) "DEBUG LNS: LNS_compute_viscous_stress_tensor finished."
    
#ifdef USE_MPI
    if (NPROC > 1 .and. ninterface_acoustic > 0) then
      !call assemble_MPI_vector_DG(gammaext_DG, buffer_LNS_gamma_P)
    endif
#endif
    
    !write(*,*) "DEBUG LNS: initialisation finished."
  endif ! Endif on (it == 1) and (i_stage == 1).
  
  if(LNS_VERBOSE>1 .and. myrank == 0 .AND. mod(it, 100)==0) then
    WRITE(*,*) "****************************************************************"
    WRITE(*,"(a,i9,a,i1,a,e23.16,a,i3,a)") "Iteration", it, ", stage ", i_stage, ", local time", timelocal, &
    ". Informations for process number ", myrank, "."
  endif
  
#ifdef USE_MPI
  if (NPROC > 1 .and. ninterface_acoustic > 0) then
    call assemble_MPI_vector_DG(LNS_drho, buffer_LNS_drho_P)
    call assemble_MPI_vector_DG(LNS_rho0dv(1,:), buffer_LNS_rho0dvx_P)
    call assemble_MPI_vector_DG(LNS_rho0dv(2,:), buffer_LNS_rho0dvz_P)
    call assemble_MPI_vector_DG(LNS_dE, buffer_LNS_dE_P)
  endif
#endif

  ! Local Discontinuous Galerkin for viscous fluxes.
  if((maxval(muext) > 0 .OR. maxval(etaext) > 0 .OR. maxval(kappa_DG) > 0) .OR. CONSTRAIN_HYDROSTATIC) then
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
    !call compute_viscous_tensors(T_DG, V_DG, LNS_drho, LNS_rho0dv(1,:), LNS_rho0dv(2,:), LNS_dE, timelocal)
    ! TODO: If CONSTRAIN_HYDROSTATIC=.true. and muext=etaext=kappa_DG=0, only T_DG is needed. Consider computing only T_DG in that case.
  endif
  
  ! Compute momentum perturbation.
  LNS_dm(1,:)=LNS_rho0dv(1,:)+LNS_drho*LNS_v0(1,:)
  LNS_dm(2,:)=LNS_rho0dv(2,:)+LNS_drho*LNS_v0(2,:)
  
  ! Compute temperature and pressure perturbation.
  LNS_dp = (gammaext_DG - ONE) &
            * (LNS_E0+LNS_dE - HALF*(LNS_rho0+LNS_drho) &
                                   *(  (LNS_v0(1,:)+LNS_rho0dv(1,:)/LNS_rho0)**2 &
                                     + (LNS_v0(2,:)+LNS_rho0dv(2,:)/LNS_rho0)**2)) &
           - LNS_p0
  LNS_dT = (  (LNS_E0+LNS_dE)/(LNS_rho0+LNS_drho) &
            - HALF*(  (LNS_v0(1,:)+LNS_rho0dv(1,:)/LNS_rho0)**2 &
                    + (LNS_v0(2,:)+LNS_rho0dv(2,:)/LNS_rho0)**2)  )/c_V &
           - LNS_T0
  
  ! Compute gradients.
  LNS_dv=ZEROcr
  where(LNS_rho0/=ZEROcr) LNS_dv(1,:)=LNS_rho0dv(1,:)/LNS_rho0 ! 'where(...)' as safeguard, as rho0=0 should not happen.
  where(LNS_rho0/=ZEROcr) LNS_dv(2,:)=LNS_rho0dv(1,:)/LNS_rho0 ! 'where(...)' as safeguard, as rho0=0 should not happen.
  call compute_gradient_TFSF(LNS_dv, LNS_dT, .true., .true., switch_gradient, nabla_dv, nabla_dT) ! Dummy variables are not optimal, but prevent from duplicating
  ! Compute \Sigma_v'.
  call LNS_compute_viscous_stress_tensor(nabla_dv, sigma_dv)
  
  ! Compute RHS.
  !call compute_forces_acoustic_DG(LNS_drho, LNS_rho0dv(1,:), LNS_rho0dv(2,:), LNS_dE, &
  !                                T_DG, V_DG, e1_DG, &
  !                                RHS_drho, RHS_rho0dvx, RHS_rho0dvz, RHS_dE, dot_e1, &
  !                                timelocal)
  
  if (time_stepping_scheme == 3 .or. time_stepping_scheme == 4) then
    ! Inverse mass matrix multiplication, in order to obtain actual RHS.
    RHS_drho(:)    = RHS_drho(:)    * rmass_inverse_acoustic_DG(:)
    RHS_rho0dvx(:) = RHS_rho0dvx(:) * rmass_inverse_acoustic_DG(:)
    RHS_rho0dvz(:) = RHS_rho0dvz(:) * rmass_inverse_acoustic_DG(:)
    RHS_dE(:)      = RHS_dE(:)      * rmass_inverse_acoustic_DG(:)
    ! Update the auxiliary register.
    ! Note: no need to zero it beforehand if scheme_A(1) is 0.
    aux_drho    = scheme_A(i_stage)*aux_drho    + deltat*RHS_drho
    aux_rho0dvx = scheme_A(i_stage)*aux_rho0dvx + deltat*RHS_rho0dvx
    aux_rho0dvz = scheme_A(i_stage)*aux_rho0dvz + deltat*RHS_rho0dvz
    aux_dE      = scheme_A(i_stage)*aux_dE      + deltat*RHS_dE
    ! Update the state register.
    LNS_drho        = LNS_drho        + scheme_B(i_stage)*aux_drho
    LNS_rho0dv(1,:) = LNS_rho0dv(1,:) + scheme_B(i_stage)*aux_rho0dvx
    LNS_rho0dv(2,:) = LNS_rho0dv(2,:) + scheme_B(i_stage)*aux_rho0dvz
    LNS_dE          = LNS_dE          + scheme_B(i_stage)*aux_dE
    
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
                    write(*, *) "* Coords", coord(1, ibool_before_perio(i, j, ispec)), coord(2, ibool_before_perio(i, j, ispec)), &
                                ".*"
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
  !    veloc_x = LNS_rho0dv(2,:) - LNS_v0(2,:)
  !    call SlopeLimit1(veloc_x, timelocal, 1)
  !    LNS_rho0dv(2,:) = veloc_x + LNS_v0(2,:)
  !  else
  !    call SlopeLimit1(LNS_rho0dv(2,:), timelocal, 3)
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

subroutine initial_state_LNS(rho0, v0, E0)
  use constants,only: CUSTOM_REAL,NGLLX,NGLLZ,PI
  use specfem_par, only: ibool_DG, nglob_DG, nspec
  use specfem_par_LNS, only: LNS_dummy_1d

  implicit none
  
  ! Input/Output.
  real(kind=CUSTOM_REAL), dimension(nglob_DG), intent(out) :: rho0, E0
  real(kind=CUSTOM_REAL), dimension(2, nglob_DG), intent(out) :: v0
  
  ! Local.
  real(kind=CUSTOM_REAL), parameter :: ZEROcr = 0._CUSTOM_REAL
  !real(kind=CUSTOM_REAL), parameter :: ONEl  = 1._CUSTOM_REAL
  !real(kind=CUSTOM_REAL), parameter :: TWOl  = 2._CUSTOM_REAL
  !real(kind=CUSTOM_REAL), parameter :: SIXl  = 6._CUSTOM_REAL
  !real(kind=CUSTOM_REAL), parameter :: HALFl = 0.5_CUSTOM_REAL
  integer :: ispec, iglob, i, j
  real(kind=CUSTOM_REAL) :: rho_P, veloc_x_P, veloc_z_P, E_P
  
  do ispec = 1, nspec
    do j = 1, NGLLZ
      do i = 1, NGLLX
       iglob = ibool_DG(i, j, ispec)
       call boundary_condition_DG(i, j, ispec, ZEROcr, rho_P, LNS_dummy_1d, LNS_dummy_1d, E_P, &
                                  veloc_x_P, veloc_z_P, LNS_dummy_1d, LNS_dummy_1d)
       ! TODO: Ultimately, implement a dedicated subroutine.
       rho0(iglob) = rho_P
       v0(1,iglob) = veloc_x_P
       v0(2,iglob) = veloc_z_P
       E0(iglob)   = E_P
      enddo
    enddo
  enddo
end subroutine initial_state_LNS

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
  real(kind=CUSTOM_REAL), dimension(2, 2, nglob_DG), intent(in) :: nabla_v
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
  sigma_v(1,:) = (TWOcr*LNS_mu+LNS_lambda)*nabla_v(1,1,:) + LNS_lambda*nabla_v(2,2,:)
  sigma_v(2,:) = LNS_mu*(nabla_v(1,2,:)+nabla_v(2,1,:))
  sigma_v(3,:) = (TWOcr*LNS_mu+LNS_lambda)*nabla_v(2,2,:) + LNS_lambda*nabla_v(1,1,:)
  
end subroutine LNS_compute_viscous_stress_tensor

! ------------------------------------------------------------ !
! compute_viscous_tensors                                      !
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
  use specfem_par_LNS ! TODO: select variables to use.
  
  implicit none
  
  ! Input/Output.
  real(kind=CUSTOM_REAL), dimension(nglob_DG), intent(out) :: SF
  real(kind=CUSTOM_REAL), dimension(2, nglob_DG), intent(out) :: TF
  logical, intent(in) :: swTF, swSF, swMETHOD
  real(kind=CUSTOM_REAL), dimension(2, nglob_DG), intent(out) :: nabla_SF
  real(kind=CUSTOM_REAL), dimension(2, 2, nglob_DG), intent(out) :: nabla_TF
  
  ! Local.
  real(kind=CUSTOM_REAL), parameter :: ZEROcr = 0._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: HALFcr = 0.5_CUSTOM_REAL
  integer :: ispec,i,j,k,iglob, iglobM, iglobP!, iglob_unique
  real(kind=CUSTOM_REAL) :: rho_DG_P, rhovx_DG_P, rhovz_DG_P, &
        E_DG_P, veloc_x_DG_P, veloc_z_DG_P, p_DG_P, T_P, &
        Tx_DG_P, Tz_DG_P, Vxx_DG_P, Vzz_DG_P, Vzx_DG_P, Vxz_DG_P, &
        flux_n, flux_x, flux_z, nx, nz, timelocal, weight, gamma_P
  logical :: exact_interface_flux
  !integer, dimension(nglob_DG) :: MPI_iglob
  integer, dimension(3) :: neighbor
  integer :: iface1, iface, iface1_neighbor, iface_neighbor, ispec_neighbor
  real(kind=CUSTOM_REAL) :: xixl,xizl,gammaxl,gammazl,jacobianl ! Jacobian matrix and determinant
  !real(kind=CUSTOM_REAL) :: temp_SFx, temp_SFz, temp_TFxx, temp_TFzx, temp_TFxz, temp_TFzz
  real(kind=CUSTOM_REAL), dimension(nglob_DG) :: &
  !       rho_DG, rhovx_DG, rhovz_DG, E_DG, veloc_x_DG, veloc_z_DG, T, &
         grad_SF_x, grad_SF_z, grad_TFx_x, grad_TFz_z, grad_TFx_z, grad_TFz_x
  real(kind=CUSTOM_REAL) :: dTFx_dxi, dTFx_dgamma, dTFz_dxi, dTFz_dgamma, dSF_dxi, dSF_dgamma ! Derivatives in \Lambda.
  !real(kind=CUSTOM_REAL) :: dTFx_dx, dTFx_dz, dTFz_dx, dTFz_dz, dSF_dx, dSF_dz
  real(kind=CUSTOM_REAL) :: wxl, wzl ! Quadrature weights.
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: temp_SFx_1, temp_SFx_2, &
        temp_SFz_1, temp_SFz_2, temp_TFxx_1, temp_TFxx_2, &
        temp_TFxz_1, temp_TFxz_2, temp_TFzx_1, temp_TFzx_2, temp_TFzz_1, temp_TFzz_2
  !real(kind=CUSTOM_REAL) :: vx_init, vz_init
  !real(kind=CUSTOM_REAL) :: ya_x_l, ya_z_l
  
  !veloc_x_DG = rhovx_DG/rho_DG
  !veloc_z_DG = rhovz_DG/rho_DG
  !T = (E_DG/rho_DG - HALFcr*(veloc_x_DG**2 + veloc_z_DG**2))/c_V
  
  if(swSF) then
    grad_SF_x  = ZEROcr
    grad_SF_z  = ZEROcr
  endif
  if(swTF) then
    grad_TFx_x = ZEROcr
    grad_TFz_z = ZEROcr
    grad_TFx_z = ZEROcr
    grad_TFz_x = ZEROcr
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
          jacobianl = jacobian(i,j,ispec)
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
          !  jacobianl = ya_x_l*ya_z_l*jacobianl
          !endif
          
          wzl = wzgll(j)
          wxl = wxgll(i)
          
          if(swMETHOD) then
            ! In that case, we want to compute \int_{\Omega} (\nabla f)\Phi as:
            ! - \int_{\Omega} T\nabla\Phi + \int_{\partial\Omega} f\Phi.
            ! The idea is to store in:
            !   grad_SF_x,
            !   grad_SF_z,
            !   grad_TFx_x,
            !   grad_TFx_z,
            !   grad_TFz_x,
            !   grad_TFz_z
            ! the values of the approximated integrals:
            !  $\int \partial_xT\Phi_x$,
            !  $\int \partial_zT\Phi_z$,
            !  $\int \partial_xVx\Phi_x$,
            !  etc.
            ! and "desintegrate" those values by multiplying the obtained vector by the inverse mass matrix. As explained before, the integrals are computed by using the divergence theorem and taking into account the surface terms.
            
            ! Compute inner contributions.
            !if(.not. CONSTRAIN_HYDROSTATIC) then
            if(swSF) then
              temp_SFx_1(i,j)  = wzl * jacobianl * (xixl * SF(iglob))
              temp_SFz_1(i,j)  = wzl * jacobianl * (xizl * SF(iglob))
              temp_SFx_2(i,j)  = wxl * jacobianl * (gammaxl * SF(iglob))
              temp_SFz_2(i,j)  = wxl * jacobianl * (gammazl * SF(iglob))
            endif
            if(swTF) then
              temp_TFxx_1(i,j) = wzl * jacobianl * (xixl * TF(1,iglob))
              temp_TFxz_1(i,j) = wzl * jacobianl * (xizl * TF(1,iglob))
              temp_TFzx_1(i,j) = wzl * jacobianl * (xixl * TF(2,iglob))
              temp_TFzz_1(i,j) = wzl * jacobianl * (xizl * TF(2,iglob))
              temp_TFxx_2(i,j) = wxl * jacobianl * (gammaxl * TF(1,iglob))
              temp_TFxz_2(i,j) = wxl * jacobianl * (gammazl * TF(1,iglob))
              temp_TFzx_2(i,j) = wxl * jacobianl * (gammaxl * TF(2,iglob))
              temp_TFzz_2(i,j) = wxl * jacobianl * (gammazl * TF(2,iglob))
            endif
            !else
            !  vx_init = rhovx_init(iglob)/rho_init(iglob)
            !  vz_init = rhovz_init(iglob)/rho_init(iglob)
            !  
            !  temp_SFx_1(i,j)  = wzl * jacobianl * (xixl * (SF(iglob) - T_init(iglob)))
            !  temp_SFz_1(i,j)  = wzl * jacobianl * (xizl * (SF(iglob) - T_init(iglob)))
            !  temp_TFxx_1(i,j) = wzl * jacobianl * (xixl * (TF(1,iglob) - vx_init))
            !  temp_TFxz_1(i,j) = wzl * jacobianl * (xizl * (TF(1,iglob))) ! Some hypothesis on initial velocity is used here.
            !  temp_TFzx_1(i,j) = wzl * jacobianl * (xixl * (TF(2,iglob) - vz_init))
            !  temp_TFzz_1(i,j) = wzl * jacobianl * (xizl * (TF(2,iglob) - vz_init))
            !  
            !  temp_SFx_2(i,j)  = wxl * jacobianl * (gammaxl * (SF(iglob) - T_init(iglob)))
            !  temp_SFz_2(i,j)  = wxl * jacobianl * (gammazl * (SF(iglob) - T_init(iglob)))
            !  temp_TFxx_2(i,j) = wxl * jacobianl * (gammaxl * (TF(1,iglob) - vx_init))
            !  temp_TFxz_2(i,j) = wxl * jacobianl * (gammazl * (TF(1,iglob))) ! Some hypothesis on initial velocity is used here.
            !  temp_TFzx_2(i,j) = wxl * jacobianl * (gammaxl * (TF(2,iglob) - vz_init))
            !  temp_TFzz_2(i,j) = wxl * jacobianl * (gammazl * (TF(2,iglob) - vz_init))
            !endif
          else
            ! In that case, we want to compute \int_{\Omega} (\nabla f)\Phi directly as it.
            ! The idea is to store in:
            !   grad_SF_x,
            !   grad_SF_z,
            !   grad_TFx_x,
            !   grad_TFx_z,
            !   grad_TFz_x,
            !   grad_TFz_z
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
              dTFx_dxi    = ZEROcr
              dTFx_dgamma = ZEROcr
              dTFz_dxi    = ZEROcr
              dTFz_dgamma = ZEROcr
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
                dTFx_dxi    = dTFx_dxi    + TF(1,ibool_DG(k,j,ispec)) * real(hprime_xx(i,k), kind=CUSTOM_REAL)
                dTFx_dgamma = dTFx_dgamma + TF(1,ibool_DG(i,k,ispec)) * real(hprime_zz(j,k), kind=CUSTOM_REAL)
                dTFz_dxi    = dTFz_dxi    + TF(2,ibool_DG(k,j,ispec)) * real(hprime_xx(i,k), kind=CUSTOM_REAL)
                dTFz_dgamma = dTFz_dgamma + TF(2,ibool_DG(i,k,ispec)) * real(hprime_zz(j,k), kind=CUSTOM_REAL)
              endif
              !else
              !  vx_init = rhovx_init(ibool_DG(k,j,ispec))/rho_init(ibool_DG(k,j,ispec))
              !  dTFx_dxi = dTFx_dxi + (TF(1,ibool_DG(k,j,ispec)) - vx_init) &
              !                      * real(hprime_xx(i,k), kind=CUSTOM_REAL)
              !  vx_init = rhovx_init(ibool_DG(i,k,ispec))/rho_init(ibool_DG(i,k,ispec))       
              !  dTFx_dgamma = dTFx_dgamma + (TF(1,ibool_DG(i,k,ispec)) - vx_init) &
              !                            * real(hprime_zz(j,k), kind=CUSTOM_REAL)
              !  vz_init = rhovz_init(ibool_DG(k,j,ispec))/rho_init(ibool_DG(k,j,ispec))
              !  dTFz_dxi = dTFz_dxi + (TF(2,ibool_DG(k,j,ispec)) - vz_init) &
              !                      * real(hprime_xx(i,k), kind=CUSTOM_REAL)
              !  vz_init = rhovz_init(ibool_DG(i,k,ispec))/rho_init(ibool_DG(i,k,ispec))
              !  dTFz_dgamma = dTFz_dgamma + (TF(2,ibool_DG(i,k,ispec)) - vz_init) &
              !                            * real(hprime_zz(j,k), kind=CUSTOM_REAL)
              !  dSF_dxi = dSF_dxi + (SF(ibool_DG(k,j,ispec)) - T_init(ibool_DG(k,j,ispec))) &
              !                    * real(hprime_xx(i,k), kind=CUSTOM_REAL)
              !  dSF_dgamma = dSF_dgamma + (SF(ibool_DG(i,k,ispec)) - T_init(ibool_DG(i,k,ispec))) &
              !                           * real(hprime_zz(j,k), kind=CUSTOM_REAL)
              !endif
            enddo

            ! Compute derivatives in element \Omega using the chain rule.
            if(swSF) then
              !dSF_dx   = dSF_dxi * xixl + dSF_dgamma * gammaxl
              !dSF_dz   = dSF_dxi * xizl + dSF_dgamma * gammazl
              !temp_SFx = dSF_dx * jacobianl
              !temp_SFz = dSF_dz * jacobianl
              !grad_SF_x(iglob) = dSF_dx
              !grad_SF_z(iglob) = dSF_dz
              grad_SF_x(iglob) = dSF_dxi * xixl + dSF_dgamma * gammaxl
              grad_SF_z(iglob) = dSF_dxi * xizl + dSF_dgamma * gammazl
            endif
            if(swTF) then
              !dTFx_dx   = dTFx_dxi * xixl + dTFx_dgamma * gammaxl
              !dTFx_dz   = dTFx_dxi * xizl + dTFx_dgamma * gammazl
              !dTFz_dx   = dTFz_dxi * xixl + dTFz_dgamma * gammaxl
              !dTFz_dz   = dTFz_dxi * xizl + dTFz_dgamma * gammazl
              !temp_TFxx = dTFx_dx * jacobianl
              !temp_TFxz = dTFx_dz * jacobianl
              !temp_TFzx = dTFz_dx * jacobianl
              !temp_TFzz = dTFz_dz * jacobianl
              !grad_TFx_x(iglob) = dTFx_dx
              !grad_TFx_z(iglob) = dTFx_dz
              !grad_TFz_x(iglob) = dTFz_dx
              !grad_TFz_z(iglob) = dTFz_dz
              grad_TFx_x(iglob) = dTFx_dxi * xixl + dTFx_dgamma * gammaxl
              grad_TFx_z(iglob) = dTFx_dxi * xizl + dTFx_dgamma * gammazl
              grad_TFz_x(iglob) = dTFz_dxi * xixl + dTFz_dgamma * gammaxl
              grad_TFz_z(iglob) = dTFz_dxi * xizl + dTFz_dgamma * gammazl
            endif
          endif
        enddo
      enddo
      
      if(swMETHOD) then
        do j = 1, NGLLZ
          do i = 1, NGLLX
            iglob = ibool_DG(i,j,ispec)
            ! Assemble the contributions using the inner contributions.
            do k = 1, NGLLX
              if(swSF) then
              grad_SF_x(iglob) = grad_SF_x(iglob) - &
                               (temp_SFx_1(k,j) * real(hprimewgll_xx(k,i), kind=CUSTOM_REAL) + &
                                temp_SFx_2(i,k) * real(hprimewgll_zz(k,j), kind=CUSTOM_REAL))
              grad_SF_z(iglob) = grad_SF_z(iglob) - &
                               (temp_SFz_1(k,j) * real(hprimewgll_xx(k,i), kind=CUSTOM_REAL) + &
                                temp_SFz_2(i,k) * real(hprimewgll_zz(k,j), kind=CUSTOM_REAL))
              endif
              if(swTF) then
              grad_TFx_x(iglob) = grad_TFx_x(iglob) - &
                                (temp_TFxx_1(k,j) * real(hprimewgll_xx(k,i), kind=CUSTOM_REAL) + &
                                 temp_TFxx_2(i,k) * real(hprimewgll_zz(k,j), kind=CUSTOM_REAL))
              grad_TFx_z(iglob) = grad_TFx_z(iglob) - &
                                (temp_TFxz_1(k,j) * real(hprimewgll_xx(k,i), kind=CUSTOM_REAL) + &
                                 temp_TFxz_2(i,k) * real(hprimewgll_zz(k,j), kind=CUSTOM_REAL))
              grad_TFz_x(iglob) = grad_TFz_x(iglob) - &
                                (temp_TFzx_1(k,j) * real(hprimewgll_xx(k,i), kind=CUSTOM_REAL) + &
                                 temp_TFzx_2(i,k) * real(hprimewgll_zz(k,j), kind=CUSTOM_REAL))
              grad_TFz_z(iglob) = grad_TFz_z(iglob) - &
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
            rho_DG_P     = ZEROcr
            rhovx_DG_P   = ZEROcr
            rhovz_DG_P   = ZEROcr
            E_DG_P       = ZEROcr
            veloc_x_DG_P = ZEROcr
            veloc_z_DG_P = ZEROcr
            p_DG_P       = ZEROcr
            T_P          = ZEROcr
            Vxx_DG_P     = ZEROcr
            Vzz_DG_P     = ZEROcr
            Vzx_DG_P     = ZEROcr
            Vxz_DG_P     = ZEROcr
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
            
            exact_interface_flux = .false. ! Reset this variable to .false.: by default, the fluxes have to be computed (jump!=0). In some specific cases (assigned during the call to compute_interface_unknowns), the flux can be exact (jump==0).
            call compute_interface_unknowns(i,j,ispec, rho_DG_P, rhovx_DG_P, &
                    rhovz_DG_P, E_DG_P, veloc_x_DG_P, veloc_z_DG_P, p_DG_P, T_P, &
                    Tx_DG_P, Tz_DG_P, Vxx_DG_P, Vzz_DG_P, Vzx_DG_P, Vxz_DG_P, gamma_P, &
                    neighbor,&
                    exact_interface_flux, &
                    rho_DG(iglobM), E_DG(iglobM), rhovx_DG(iglobM), rhovz_DG(iglobM), &
                    nabla_TF(:,:,iglobM), nabla_SF(:,iglobM), &
                    rho_DG(iglobP), E_DG(iglobP), rhovx_DG(iglobP), rhovz_DG(iglobP), &
                    nabla_TF(:,:,iglobP), nabla_SF(:,iglobP), &
                    nx, nz, weight, timelocal,iface1, iface)

            !vx_init = rhovx_init(iglobM)/rho_init(iglobM)
            !vz_init = rhovz_init(iglobM)/rho_init(iglobM)

            ! Dot products.
            if(swSF) then
              flux_x = SF(iglobM) + T_P
              !if(CONSTRAIN_HYDROSTATIC) flux_x = flux_x - 2*T_init(iglobM)
              flux_n = flux_x*nx
              grad_SF_x(iglobM) = grad_SF_x(iglobM) + weight*flux_n*HALFcr
              flux_z = SF(iglobM) + T_P
              !if(CONSTRAIN_HYDROSTATIC) flux_z = flux_z - 2*T_init(iglobM)
              flux_n = flux_z*nz
              grad_SF_z(iglobM) = grad_SF_z(iglobM) + weight*flux_n*HALFcr
            endif
            if(swTF) then
              flux_x = TF(1,iglobM) + veloc_x_DG_P
              !if(CONSTRAIN_HYDROSTATIC) flux_x = flux_x - 2*vx_init
              flux_n = flux_x*nx
              grad_TFx_x(iglobM) = grad_TFx_x(iglobM) + weight*flux_n*HALFcr
              flux_z = TF(1,iglobM) + veloc_x_DG_P
              !if(CONSTRAIN_HYDROSTATIC) flux_z = flux_z - 2*vx_init
              flux_n = flux_z*nz
              grad_TFx_z(iglobM) = grad_TFx_z(iglobM) + weight*flux_n*HALFcr
              flux_x = TF(2,iglobM) + veloc_z_DG_P
              !if(CONSTRAIN_HYDROSTATIC) flux_x = flux_x - 2*vz_init
              flux_n = flux_x*nx
              grad_TFz_x(iglobM) = grad_TFz_x(iglobM) + weight*flux_n*HALFcr
              flux_z = TF(2,iglobM) + veloc_z_DG_P
              !if(CONSTRAIN_HYDROSTATIC) flux_z = flux_z - 2*vz_init
              flux_n = flux_z*nz
              grad_TFz_z(iglobM) = grad_TFz_z(iglobM) + weight*flux_n*HALFcr
            endif
          enddo
        enddo
      endif ! Endif on swMETHOD.
    endif ! End of test if acoustic element.
  enddo ! End of loop over elements.
  
  if(swMETHOD) then
    ! "Desintegrate".
    if(swSF) then
      grad_SF_x(:)  = grad_SF_x(:) * rmass_inverse_acoustic_DG(:)
      grad_SF_z(:)  = grad_SF_z(:) * rmass_inverse_acoustic_DG(:)
    endif
    if(swTF) then
      grad_TFx_x(:) = grad_TFx_x(:) * rmass_inverse_acoustic_DG(:)
      grad_TFx_z(:) = grad_TFx_z(:) * rmass_inverse_acoustic_DG(:)
      grad_TFz_x(:) = grad_TFz_x(:) * rmass_inverse_acoustic_DG(:)
      grad_TFz_z(:) = grad_TFz_z(:) * rmass_inverse_acoustic_DG(:)
    endif
  endif
  
  ! Store in variables intended for output.
  if(swSF) then
    nabla_SF(1, :) = grad_SF_x
    nabla_SF(2, :) = grad_SF_z
  endif
  if(swTF) then
    nabla_TF(1, 1, :) = grad_TFx_x ! dxTFx
    nabla_TF(1, 2, :) = grad_TFx_z ! dzTFx
    nabla_TF(2, 1, :) = grad_TFz_x ! dxTFz
    nabla_TF(2, 2, :) = grad_TFz_z ! dzTFz
  endif
end subroutine compute_gradient_TFSF

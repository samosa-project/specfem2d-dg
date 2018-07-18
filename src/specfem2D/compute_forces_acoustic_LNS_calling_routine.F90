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
    LNS_rho0, LNS_v0x,     LNS_v0z,     LNS_E0,&
    LNS_drho, LNS_rho0dvx, LNS_rho0dvz, LNS_dE,&
    aux_drho, aux_rho0dvx, aux_rho0dvz, aux_dE,&
    RHS_drho, RHS_rho0dvx, RHS_rho0dvz, RHS_dE,&
    nabla_v0, LNS_p0, LNS_T0, sigma_v_0,&
    buffer_LNS_drho_P, buffer_LNS_rho0dvx_P, buffer_LNS_rho0dvz_P, buffer_LNS_dE_P,&
    LNS_VERBOSE

  implicit none

  ! Local variables.
  real(kind=CUSTOM_REAL), parameter :: ZEROl = 0._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: threshold = 0.0000001_CUSTOM_REAL
  real(kind=CUSTOM_REAL) :: timelocal
  real(kind=CUSTOM_REAL), dimension(stage_time_scheme) :: scheme_A, scheme_B, scheme_C
  integer :: i,j,ispec
  
  logical CHECK_NONPOSITIVITY, CHECK_NONPOSITIVITY_ON_ALL_PROCS, CHECK_NONPOSITIVITY_FIND_POINT
  
  ! Checks if anything has to be done.
  if (.not. any_acoustic_DG) then
    return
  endif
  
  ! Those are debug switches. They are computationaly very heavy and should not be used on every simulation.
  ! CHECK_NONPOSITIVITY_ON_ALL_PROCS=.true. and CHECK_NONPOSITIVITY_FIND_POINT=.true. are particularly heavy.
  CHECK_NONPOSITIVITY=.false.              ! Set to .true. to enable nonpositivity checking.
  CHECK_NONPOSITIVITY_ON_ALL_PROCS=.false. ! Only used if CHECK_NONPOSITIVITY==.true.. Set to .false. for checking only on proc 0. Set to .true. for checking on all procs.
  CHECK_NONPOSITIVITY_FIND_POINT=.false.   ! Only used if CHECK_NONPOSITIVITY==.true.. Set to .true. to find where nonpositivity was encountered.
  
  ! Intialisation.
  if(it == 1 .AND. i_stage == 1) then
    !if(USE_SLOPE_LIMITER) then
    !  ! The Vandermonde matrices are only used when the slope limiter is.
    !  call setUpVandermonde()
    !endif
    
    ! Zero registers.
    aux_drho    = ZEROl
    aux_rho0dvx = ZEROl
    aux_rho0dvz = ZEROl
    aux_dE      = ZEROl
    RHS_drho    = ZEROl
    RHS_rho0dvx = ZEROl
    RHS_rho0dvz = ZEROl
    RHS_dE      = ZEROl
    T_DG        = ZEROl
    V_DG        = ZEROl
    
    ! Physical parameters.
    if(.not. assign_external_model) then
      deallocate(gravityext, muext, etaext, kappa_DG)
      allocate(gravityext(NGLLX, NGLLZ, nspec), &
               etaext(NGLLX, NGLLZ, nspec), &
               muext(NGLLX, NGLLZ, nspec), &
               kappa_DG(NGLLX, NGLLZ, nspec)) 
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
    LNS_drho    = ZEROl
    LNS_rho0dvx = ZEROl
    LNS_rho0dvz = ZEROl
    LNS_dE      = ZEROl

    ! Save initial state.
    allocate(LNS_rho0(nglob_DG), &
             LNS_v0x(nglob_DG), &
             LNS_v0z(nglob_DG), &
             LNS_E0(nglob_DG))
    ! TODO: correctly initialise initial state (use somthing like the "initial_condition_DG" subroutine).
    LNS_rho0   = LNS_drho ! Not true.
    LNS_v0x    = LNS_rho0dvx ! Not true.
    LNS_v0z    = LNS_rho0dvz ! Not true.
    LNS_E0     = LNS_dE ! Not true.
    
    ! Compute \nabla v_0.
    ! TODO: compute nabla_v0.
    write(*,*) nabla_v0 ! Temporary hack to shut compilation errors down.
    
    ! Compute T_0 and p_0.
    LNS_T0 = (LNS_E0/LNS_rho0 - HALF*(LNS_v0x**2+LNS_v0z**2))/(c_V)
    LNS_p0 = (gammaext_DG - ONE) &
              * (LNS_E0 - HALF*LNS_rho0*(LNS_v0x**2+LNS_v0z**2))
    where(LNS_p0 <= 0._CUSTOM_REAL) LNS_p0 = 1._CUSTOM_REAL ! When elastic-DG simulations, p_DG = 0 in elastic elements and 1/p_DG will not be properly defined. Use a hack.
    
    ! Compute \Sigma_v_0.
    ! TODO: compute sigma_v_0
    write(*,*) sigma_v_0 ! Temporary hack to shut compilation errors down.
    
#ifdef USE_MPI
    if (NPROC > 1 .and. ninterface_acoustic > 0) then
      !call assemble_MPI_vector_DG(gammaext_DG, buffer_LNS_gamma_P)
    endif
#endif

  endif ! Endif on (it == 1) and (i_stage == 1).
  
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
  else
    if(myrank==0) then
      write(*,*) "********************************"
      write(*,*) "*            ERROR             *"
      write(*,*) "********************************"
      write(*,*) "* This time_stepping_scheme is *"
      write(*,*) "* not implemented for DG yet.  *"
      write(*,*) "* time_stepping_scheme ", time_stepping_scheme
      write(*,*) "********************************"
      stop
    endif
  endif
  
  timelocal = (it-1)*deltat + scheme_C(i_stage)*deltat ! Compute current time.
  
  if(LNS_VERBOSE>1 .and. myrank == 0 .AND. mod(it, 100)==0) then
    WRITE(*,*) "****************************************************************"
    WRITE(*,"(a,i9,a,i1,a,e23.16,a,i3,a)") "Iteration", it, ", stage ", i_stage, ", local time", timelocal, &
    ". Informations for process number ", myrank, "."
  endif
  
#ifdef USE_MPI
  if (NPROC > 1 .and. ninterface_acoustic > 0) then
    call assemble_MPI_vector_DG(LNS_drho, buffer_LNS_drho_P)
    call assemble_MPI_vector_DG(LNS_rho0dvx, buffer_LNS_rho0dvx_P)
    call assemble_MPI_vector_DG(LNS_rho0dvz, buffer_LNS_rho0dvz_P)
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
    !call compute_viscous_tensors(T_DG, V_DG, LNS_drho, LNS_rho0dvx, LNS_rho0dvz, LNS_dE, timelocal)
    ! TODO: If CONSTRAIN_HYDROSTATIC=.true. and muext=etaext=kappa_DG=0, only T_DG is needed. Consider computing only T_DG in that case.
  endif
  
  ! Compute RHS.
  !call compute_forces_acoustic_DG(LNS_drho, LNS_rho0dvx, LNS_rho0dvz, LNS_dE, &
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
    LNS_drho    = LNS_drho    + scheme_B(i_stage)*aux_drho
    LNS_rho0dvx = LNS_rho0dvx + scheme_B(i_stage)*aux_rho0dvx
    LNS_rho0dvz = LNS_rho0dvz + scheme_B(i_stage)*aux_rho0dvz
    LNS_dE      = LNS_dE      + scheme_B(i_stage)*aux_dE
    
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
  else
    if(myrank==0) then
      ! Safeguard only, as normally the error is prompted before (see above).
      write(*,*) "********************************"
      write(*,*) "*            ERROR             *"
      write(*,*) "********************************"
      write(*,*) "* This time_stepping_scheme is *"
      write(*,*) "* not implemented for DG yet.  *"
      write(*,*) "* time_stepping_scheme ", time_stepping_scheme
      write(*,*) "********************************"
      stop
    endif
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
  !    veloc_x = LNS_rho0dvx - LNS_v0x
  !    call SlopeLimit1(veloc_x, timelocal, 1)
  !    LNS_rho0dvx = veloc_x + LNS_v0x
  !  else
  !    call SlopeLimit1(LNS_rho0dvx, timelocal, 2)
  !  endif
  !  ! rhovz.
  !  if(CONSTRAIN_HYDROSTATIC) then
  !    veloc_x = LNS_rho0dvz - LNS_v0z
  !    call SlopeLimit1(veloc_x, timelocal, 1)
  !    LNS_rho0dvz = veloc_x + LNS_v0z
  !  else
  !    call SlopeLimit1(LNS_rho0dvz, timelocal, 3)
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
! prepare_MPI_DG                                               !
! ------------------------------------------------------------ !
! See compute_forces_acoustic_DG_calling_routine.

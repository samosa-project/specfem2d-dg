! ------------------------------------------------------------ !
! compute_forces_acoustic_LNS                                  !
! ------------------------------------------------------------ !
! TODO: Description.

subroutine compute_forces_acoustic_LNS(cv_drho, cv_rho0dv, cv_dE, & ! Constitutive variables.
                                       in_dm, in_dp, in_nabla_dT, in_sigma_dv, & ! Precomputed quantities.
                                       outrhs_drho, outrhs_rho0dv, outrhs_dE, & ! Output (RHS for each constitutive variable).
                                       currentTime) ! Time.
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
  real(kind=CUSTOM_REAL), dimension(nglob_DG), intent(in) :: cv_drho, cv_dE, in_dp!, in_dT
  real(kind=CUSTOM_REAL), dimension(NDIM, nglob_DG), intent(in) :: cv_rho0dv, in_dm, in_nabla_dT
  real(kind=CUSTOM_REAL), dimension(NVALSIGMA, nglob_DG), intent(in) :: in_sigma_dv
  real(kind=CUSTOM_REAL), dimension(nglob_DG), intent(out) :: outrhs_drho, outrhs_dE
  real(kind=CUSTOM_REAL), dimension(NDIM, nglob_DG), intent(out) :: outrhs_rho0dv
  real(kind=CUSTOM_REAL), intent(in) :: currentTime
  
  ! Local.
  real(kind=CUSTOM_REAL), parameter :: ZEROcr = 0._CUSTOM_REAL
  !real(kind=CUSTOM_REAL), parameter :: ONE  = 1._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: TWOcr  = 2._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: HALFcr = 0.5_CUSTOM_REAL
  
  ! MESSY
  !real(kind=CUSTOM_REAL), dimension(nglob_DG) :: rho_DG, rhovx_DG, rhovz_DG, E_DG
  !real(kind=CUSTOM_REAL), dimension(2, nglob_DG) :: T_DG
  !real(kind=CUSTOM_REAL), dimension(2, 2, nglob_DG) :: V_DG
  integer :: ispec, i,j, k,iglob,SPCDM, iglob_unique
  real(kind=CUSTOM_REAL), dimension(NGLLX, NGLLZ, NDIM) :: cntrb_drho, cntrb_dE
  real(kind=CUSTOM_REAL), dimension(NDIM, NGLLX, NGLLZ, NDIM) :: cntrb_rho0dv
  ! Variables "cntrb_*" are aimed at assembling the different contributions to constitutive variables.
  ! For drho and dE, "cntrb_*(:,:,1)" contains the contribution along xi, and "cntrb_*(:,:,2)" the contribution along gamma.
  ! For rho0dv, "cntrb_*(1,:,:,1)" contains the contribution to rho0dvx along xi, "cntrb_*(1,:,:,2)" the contribution to rho0dvx along gamma, "cntrb_*(2,:,:,1)" the contribution to rhodvz along xi, and "cntrb_*(2,:,:,2)" the contribution to rhodvz along gamma.
  real(kind=CUSTOM_REAL), dimension(NGLLX, NGLLZ) :: d0cntrb_dE
  real(kind=CUSTOM_REAL), dimension(NDIM, NGLLX,NGLLZ) :: d0cntrb_rho0dv
  ! Variable "d0cntrb_*" are aimed at assembling contributions of zero-th degree (outside divergence operator).
  
  
  real(kind=CUSTOM_REAL) :: jacLoc ! Jacobian matrix determinant.
  real(kind=CUSTOM_REAL) :: lambda, halfWeight, jump!, flux_n
  
  real(kind=CUSTOM_REAL), dimension(NDIM) :: locSigma!, locG ! Local tensor, local gravity.
  real(kind=CUSTOM_REAL), dimension(NDIM,NDIM) :: xigammal
  real(kind=CUSTOM_REAL), dimension(NDIM) :: wzlwxljacLoc
  real(kind=CUSTOM_REAL), dimension(NDIM) :: n_out
  real(kind=CUSTOM_REAL), dimension(NDIM) :: rho0dv_P
  integer :: iglobM, iglobP
  integer, dimension(3) :: neighbor
  
!  ! Variables specifically for PML.
! 
!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  ! IMPORTANT INFORMATION:
!  ! The following block of code is commented out to prevent confusion.
!  ! But it might be useful in the future, when someone actually tries to finish the full PML implementation.
!  ! Please do not remove those lines, even if they are commented out for now.
!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
!  real(kind=CUSTOM_REAL), dimension(NGLLX, NGLLZ) :: d0cntrb_drho
!  integer :: ispec_PML, iglobPML
!  real(kind=CUSTOM_REAL) :: pml_a0,pml_a1
!  real(kind=CUSTOM_REAL), dimension(NDIM) :: pml_b, pml_boa, pml_alp!, pml_ade_rho0dv ! Those two are of dimension NDIM, but not for the same reason: for rho0dv it's because it's actually of dimension NDIM, while for pml_b it's because we happen to have as many ADEs as space dimensions.
  
  ! Variables specifically for LNS_get_interfaces_unknowns.
  real(kind=CUSTOM_REAL), dimension(NDIM) :: dv_P, dm_P, nabla_dT_P
  real(kind=CUSTOM_REAL) :: drho_P, dp_P, dE_P
  real(kind=CUSTOM_REAL), dimension(NVALSIGMA) :: sigma_dv_P
  logical :: viscousComputation
  real(kind=CUSTOM_REAL) :: wxlwzl!,wxljacLoc, wzljacLoc ! Integration weigths.
  logical :: exact_interface_flux
  integer :: iface1, iface, iface1_neighbor, iface_neighbor, ispec_neighbor
  real(kind=CUSTOM_REAL) :: ya_x_l, ya_z_l ! Stretching absorbing boundary conditions.
  
!  ! Variables specifically for manufactured solutions.
!  real(kind=CUSTOM_REAL), dimension(NDIM) :: X
!  real(kind=CUSTOM_REAL) :: GAM
  
  ! Initialisation of the RHS.
  outrhs_drho    = ZEROcr
  outrhs_rho0dv  = ZEROcr
  outrhs_dE      = ZEROcr
  
  ! Start by adding source terms.
!  if(.false.) then ! TEST MANUFACTURED SOLUTIONS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  select case (TYPE_SOURCE_DG)
    case (1)
      call LNS_mass_source(outrhs_drho, outrhs_rho0dv, outrhs_dE, it, i_stage)
    !case (2)
    !  do SPCDM = 1, NDIM
    !    call compute_add_sources_acoustic_DG_spread(outrhs_rho0dv(SPCDM,:), it, i_stage)
    !  enddo
    case (3)
      call compute_add_sources_acoustic_DG_spread(outrhs_dE, it, i_stage)
    case default
      stop "TYPE_SOURCE_DG not implemented."
  end select
!  endif ! TEST MANUFACTURED SOLUTIONS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  ! TEST MANUFACTURED SOLUTIONS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  do ispec = 1, nspec; do j = 1, NGLLZ; do i = 1, NGLLX ! V1: get on all GLL points
!  !do ispec = 1, nspec; do j = 2,4; do i = 2,4 ! V2: get on center GLL point only
!    iglob = ibool_DG(i,j,ispec)
!#define RHO_cst 0.001
!#define dRHO_x 1.
!#define dRHO_z 2.
!#define E_cst 0.05
!#define dE_x 3.
!#define dE_z 4.
!    X = coord(:,ibool_before_perio(i,j,ispec))
!    GAM = gammaext_DG(iglob)
!    outrhs_drho(iglob) = &!outrhs_drho(iglob) + &
!                           RHO_cst*dRHO_x*PI*cos(dRHO_x*PI*X(1))*LNS_v0(1,iglob) &
!                         * wxgll(i)*wzgll(j)*jacobian(i,j,ispec)
!    
!    !outrhs_rho0dv(1, iglob) = outrhs_rho0dv(1, iglob) + 1.
!    outrhs_rho0dv(1, iglob) = &!outrhs_rho0dv(1, iglob) + &
!                                (GAM-1.) &
!                                  *(  E_cst*dE_x*PI*cos(dE_x*PI*X(1)) &
!                                    - 0.5*RHO_cst*dRHO_x*PI*cos(dRHO_x*PI*X(1))*LNS_v0(1,iglob)**2 ) &
!                              * wxgll(i)*wzgll(j)*jacobian(i,j,ispec)
!    
!    outrhs_rho0dv(2, iglob) = &!outrhs_rho0dv(2, iglob) + &
!                                (GAM-1.) &
!                                  *(  E_cst*dE_z*PI*cos(dE_z*PI*X(2)) &
!                                    - 0.5*RHO_cst*dRHO_z*PI*cos(dRHO_z*PI*X(2))*LNS_v0(1,iglob)**2 ) &
!                              * wxgll(i)*wzgll(j)*jacobian(i,j,ispec)
!    
!    outrhs_dE(iglob) = &!outrhs_dE(iglob) + &
!                         LNS_v0(1,iglob)*PI &
!                           *(  GAM*E_cst*dE_x*cos(dE_x*PI*X(1)) &
!                             - 0.5*(GAM-1.)*RHO_cst*dRHO_x*cos(dRHO_x*PI*X(1))*LNS_v0(1,iglob)**2) &
!                       * wxgll(i)*wzgll(j)*jacobian(i,j,ispec)
!  enddo; enddo; enddo
!  !! DO NOT FORGET TO UPDATE BOUNDARY TERMS BELOW, lines ~1100.
!  ! TEST MANUFACTURED SOLUTIONS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  ! TEST MANUFACTURED SOLUTIONS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  ! Look what's inside source
!  do ispec = 1, nspec; do j = 1, NGLLZ; do i = 1, NGLLX
!    iglob = ibool_DG(i,j,ispec)
!    if(coord(1,ibool_before_perio(i,j,ispec))>0.21 .and. coord(1,ibool_before_perio(i,j,ispec))<0.22) then
!      if(coord(2,ibool_before_perio(i,j,ispec))>0.16 .and. coord(2,ibool_before_perio(i,j,ispec))<0.21) then
!        write(*,*) 'x z outrhs_rho0dv(1, iglob)', coord(1,ibool_before_perio(i,j,ispec)), &
!                   coord(2,ibool_before_perio(i,j,ispec)), outrhs_rho0dv(1, iglob)
!      endif
!    endif
!  enddo; enddo; enddo
!  stop
!  ! TEST MANUFACTURED SOLUTIONS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  if(myrank == 0 .and. LNS_VERBOSE>=51 .and. mod(it, LNS_MODPRINT) == 0) then
    write(*,"(a,i6,a)") "Informations for process number ", myrank, "."
    write(*,"(a)")                   " quantity [                 max    ,                  min    ]"
    WRITE(*,"(a,e24.16,a,e24.16,a)") " drho     [", maxval(cv_drho), ", ", minval(cv_drho), "]"
    do SPCDM = 1, NDIM
      WRITE(*,"(a,e24.16,a,e24.16,a,i1)") " rho0dv_i [", maxval(cv_rho0dv(SPCDM,:)), ", ", minval(cv_rho0dv(SPCDM,:)), "], i=",SPCDM
    enddo
    WRITE(*,"(a,e24.16,a,e24.16,a)") " dE       [", maxval(cv_dE), ", ", minval(cv_dE), "]"
    !WRITE(*,"(a,e23.16,a)")        "Ratio |p-p_{init}|/p_{init}:", maxval(abs((p_DG-p_DG_init)/p_DG_init)), "."
  endif
  
  do ispec = 1, nspec ! Loop over elements.
    if (ispec_is_acoustic_DG(ispec)) then ! Only do something for DG elements.
! 
!          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!          ! IMPORTANT INFORMATION:
!          ! The following block of code is commented out to prevent confusion.
!          ! But it might be useful in the future, when someone actually tries to finish the full PML implementation.
!          ! Please do not remove those lines, even if they are commented out for now.
!          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
!      if(PML_BOUNDARY_CONDITIONS .and. ispec_is_PML(ispec)) then
!        ! If PML, save ispec_PML now (instead of recalling it for each (i,j)).
!        ispec_PML=spec_to_PML(ispec)
!      endif
      
      ! --------------------------- !
      ! First set of loops: compute !
      ! volumic contributions.      !
      ! --------------------------- !
      ! First, zero the temporary contribution fields.
      ! This might not be needed, since they are strictly overwritten below, but this is safer.
      cntrb_drho     = ZEROcr
      cntrb_rho0dv   = ZEROcr
      cntrb_dE       = ZEROcr
    !  d0cntrb_drho   = ZEROcr
      d0cntrb_rho0dv = ZEROcr
      d0cntrb_dE     = ZEROcr
      
      ! Then, loop on the GLL points to set each of the (i, j) component to the temporary contribution fields.
      do j = 1, NGLLZ
        do i = 1, NGLLX
          iglob = ibool_DG(i,j,ispec)
! 
!          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!          ! IMPORTANT INFORMATION:
!          ! The following block of code is commented out to prevent confusion.
!          ! But it might be useful in the future, when someone actually tries to finish the full PML implementation.
!          ! Please do not remove those lines, even if they are commented out for now.
!          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
!          if(PML_BOUNDARY_CONDITIONS .and. anyabs .and. ispec_is_PML(ispec)) then
!            ! If PML, save iglobPML now.
!            iglobPML = ibool_LNS_PML(i, j, ispec_PML)
!          endif
          jacLoc = jacobian(i,j,ispec)
          xigammal(1,    1)    = xix(i,j,ispec) ! = \partial_x\xi from report
          xigammal(1,    NDIM) = xiz(i,j,ispec) ! = \partial_z\xi from report
          xigammal(NDIM, 1)    = gammax(i,j,ispec) ! = \partial_x\eta from report
          xigammal(NDIM, NDIM) = gammaz(i,j,ispec) ! = \partial_z\eta from report
          
          if(ABC_STRETCH .and. stretching_buffer(ibool_before_perio(i,j,ispec))>0) then
            ! For real stretching, \Sigma for each constitutive variable becomes \Ya\Sigma. It is heavy to change each and every expression where \Sigma arises. Rather, we make use of what is multiplying \Sigma.
            ! Here, in the volume integrations, easiest is the local derivatives involved in the change of reference.
            ! \partial_x becomes \ya_x\partial_x, and since \partial_x=(\partial_x\xi)\partial_\xi+(\partial_x\eta)\partial_\eta, only updating \partial_x\xi and \partial_x\eta is enough. Idem for \partial_z. Hence, only updating xix to \ya_x * xix, xiz to \ya_z * xiz, etc. is enough to update the operator.
            iglob_unique      = ibool_before_perio(i,j,ispec)
            ya_x_l            = stretching_ya(1, iglob_unique)
            ya_z_l            = stretching_ya(2, iglob_unique)
            xigammal(:, 1)    = ya_x_l * xigammal(:, 1)    ! Multiply x component by ya_x. ! If you change anything, remember to do it also in the 'compute_gradient_TFSF' subroutine.
            xigammal(:, NDIM) = ya_z_l * xigammal(:, NDIM) ! Multiply x component by ya_z. ! If you change anything, remember to do it also in the 'compute_gradient_TFSF' subroutine.
            !jacLoc            = ya_x_l * ya_z_l * jacLoc     ! TODO: something on jacobian?? ! If you change anything, remember to do it also in the 'compute_gradient_TFSF' subroutine.
          endif
          
! 
!          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!          ! IMPORTANT INFORMATION:
!          ! The following block of code is commented out to prevent confusion.
!          ! But it might be useful in the future, when someone actually tries to finish the full PML implementation.
!          ! Please do not remove those lines, even if they are commented out for now.
!          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
!          if(PML_BOUNDARY_CONDITIONS .and. ispec_is_PML(ispec)) then
!            ! Need to include stretching on base stress tensor.
!            ! We do it as follows. While not very readable, it gets the job done in a somewhat efficient way.
!            ! \Sigma_x needs to become \kappa_2\Sigma_x, and \Sigma_z needs to become \kappa_1\Sigma_z.
!            ! Since all \Sigma_i contributions are added through the dot product with (xix,xiz,gammax,gammaz), we include kappa_i in those directly.
!            xigammal(:, 1)    = LNS_PML_kapp(NDIM, i,j,ispec_PML)*xigammal(:, 1)    ! Multiply x component by kappa_2.
!            xigammal(:, NDIM) = LNS_PML_kapp(1,    i,j,ispec_PML)*xigammal(:, NDIM) ! Multiply z component by kappa_1.
!          endif
          
          wzlwxljacLoc(1)    = wzgll(j)*jacLoc ! Notice how 1 is z and 2 is x.
          wzlwxljacLoc(NDIM) = wxgll(i)*jacLoc
          
          if(LNS_viscous) then ! Check if viscosity exists whatsoever.
            ! Activate/deactivate, for this particular point (iglob), computation of quantities only needed when viscosity is present.
            if(     LNS_mu(iglob) > TINYVAL &
               .OR. LNS_eta(iglob) > TINYVAL &
               .OR. LNS_kappa(iglob) > TINYVAL) then
              viscousComputation = .true.
            else
              viscousComputation = .false.
            endif
          else
            viscousComputation = .false. ! If viscosity is globally disabled, deactivate it for this element.
          endif
          
          ! Note:
          ! All contributions to dimension 'SPCDM' are built using the pattern:
          !   wzlwxljacLoc(SPCDM)*DOT_PRODUCT(xigammal(SPCDM,:), locSigma)
          ! in order to store successively contributions to xi (\partial_x\xi * \Sigma_x + \partial_z\xi * \Sigma_z) and to gamma (\partial_x\gamma * \Sigma_x + \partial_z\gamma * \Sigma_z).
          
          ! 1) Inviscid stress tensor's contributions.
          ! 1.1-2) Group mass conservation and x-Momentum in loop for efficiency.
          locSigma(1)    = cv_rho0dv(1,iglob)*LNS_v0(1,iglob) + in_dp(iglob) ! x-Momentum \Sigma_1.
          locSigma(NDIM) = cv_rho0dv(NDIM,iglob)*LNS_v0(1,iglob) ! x-Momentum \Sigma_2.
          do SPCDM = 1, NDIM
            cntrb_drho(i,j,SPCDM)     = wzlwxljacLoc(SPCDM)*DOT_PRODUCT(xigammal(SPCDM,:), in_dm(:,iglob))
            cntrb_rho0dv(1,i,j,SPCDM) = wzlwxljacLoc(SPCDM)*DOT_PRODUCT(xigammal(SPCDM,:), locSigma)
          enddo
          ! 1.2) z-Momentum.
          locSigma(1)    = cv_rho0dv(1,iglob)*LNS_v0(NDIM,iglob)
          locSigma(NDIM) = cv_rho0dv(NDIM,iglob)*LNS_v0(NDIM,iglob) + in_dp(iglob)
          do SPCDM = 1, NDIM
            cntrb_rho0dv(NDIM,i,j,SPCDM) =   wzlwxljacLoc(SPCDM)&
                                           * DOT_PRODUCT(xigammal(SPCDM,:), locSigma)
          enddo
          
          ! 2) Add viscous stress tensor's contributions.
          if(viscousComputation) then
            ! 2.1) Mass conservation: no viscous contribution.
            ! 2.2) x-Momentum.
            locSigma(1)    = -in_sigma_dv(1,iglob) ! Recall, this corresponds to \Sigma_v(v')_{1,1}.
            locSigma(NDIM) = -in_sigma_dv(2,iglob) ! Recall, this corresponds to \Sigma_v(v')_{1,2} (and \Sigma_v(v')_{2,1}).
            do SPCDM = 1, NDIM
              cntrb_rho0dv(1,i,j,SPCDM) =   cntrb_rho0dv(1,i,j,SPCDM) &
                                          + wzlwxljacLoc(SPCDM)*DOT_PRODUCT(xigammal(SPCDM,:), locSigma)
            enddo
            ! 2.2) z-Momentum.
            locSigma(1)    = -in_sigma_dv(2,iglob) ! Recall, this corresponds to \Sigma_v(v')_{2,1} (and \Sigma_v(v')_{1,2}).
            locSigma(NDIM) = -in_sigma_dv(3,iglob) ! Recall, this corresponds to \Sigma_v(v')_{2,2}.
            do SPCDM = 1, NDIM
              cntrb_rho0dv(NDIM,i,j,SPCDM) =   cntrb_rho0dv(NDIM,i,j,SPCDM) &
                                             + wzlwxljacLoc(SPCDM)*DOT_PRODUCT(xigammal(SPCDM,:), locSigma)
            enddo
          endif ! Endif on viscousComputation.
          
          ! 3) Special case: energy.
          ! Indeed, if they exist, viscous contributions can be grouped to inviscid ones.
          if(viscousComputation) then
            locSigma(1) =   LNS_dv(1,iglob)*(LNS_E0(iglob) + LNS_p0(iglob) - sigma_v_0(1,iglob)) &
                          + LNS_v0(1,iglob)*(cv_dE(iglob) + in_dp(iglob) - in_sigma_dv(1,iglob)) &
                          - LNS_dv(NDIM,iglob)*sigma_v_0(2,iglob) &
                          - LNS_v0(NDIM,iglob)*in_sigma_dv(2,iglob) &
                          - LNS_kappa(iglob)*in_nabla_dT(1,iglob)
            locSigma(NDIM) =   LNS_dv(NDIM,iglob)*(LNS_E0(iglob) + LNS_p0(iglob) - sigma_v_0(3,iglob)) &
                             + LNS_v0(NDIM,iglob)*(cv_dE(iglob) + in_dp(iglob) - in_sigma_dv(3,iglob)) &
                             - LNS_dv(1,iglob)*sigma_v_0(2,iglob) &
                             - LNS_v0(1,iglob)*in_sigma_dv(2,iglob) &
                             - LNS_kappa(iglob)*in_nabla_dT(NDIM,iglob)
          else ! Else on viscousComputation.
            do SPCDM = 1, NDIM
              locSigma(SPCDM) =   LNS_dv(SPCDM,iglob)*(LNS_E0(iglob) + LNS_p0(iglob)) &
                                + LNS_v0(SPCDM,iglob)*(cv_dE(iglob) + in_dp(iglob))
            enddo
          endif ! Endif on viscousComputation.
          do SPCDM = 1, NDIM
            cntrb_dE(i,j,SPCDM) = wzlwxljacLoc(SPCDM)*DOT_PRODUCT(xigammal(SPCDM,:), locSigma)
          enddo
          
          ! 4) Zero-th degree contributions.
          !locG(1) = potential_dphi_dx_DG(ibool(i,j,ispec)); ! TODO: something more general than this hack
          !locG(1) = ZEROcr; ! TODO: something more general than this hack
          !locG(NDIM) = - LNS_g(iglob)
          ! Notice the "-" sign, needed to make sure these terms are added on the right "side" of the RHS.
          ! 4.1) Mass conservation: none.
          !d0cntrb_drho(i,j)     = ZERO
          ! 4.2) Momenta.
          ! 4.2.1) Version 1: most general.
          !do SPCDM = 1, NDIM
          !  d0cntrb_rho0dv(SPCDM,i,j) = - (   cv_drho(iglob)*locG(SPCDM) & ! \rho'g_i
          !                                  - DOT_PRODUCT(in_dm(:,iglob), nabla_v0(SPCDM,:,iglob))) ! {\delta_m}\cdot\nabla v_{0,i}
          !enddo
          ! 4.2.2) Version 2: [gravity is vertical (locG(1)=0, locG(2)=-LNS_g)].
          !do SPCDM = 1, NDIM
          d0cntrb_rho0dv(1,i,j)    = DOT_PRODUCT(in_dm(:,iglob), nabla_v0(SPCDM,:,iglob)) ! {\delta_m}\cdot\nabla v_{0,i}
          d0cntrb_rho0dv(NDIM,i,j) =   cv_drho(iglob)*LNS_g(iglob) & ! \rho'g_i
                                     + DOT_PRODUCT(in_dm(:,iglob), nabla_v0(SPCDM,:,iglob)) ! {\delta_m}\cdot\nabla v_{0,i}
          !enddo
          ! OLD4.2.2) Version 2: under HV0 (v_{0,z}=0), and SM (d_xv_0=0), dm part of the zero-th degree RHS is simplified.
          !d0cntrb_rho0dv(1,i,j) = - (  cv_drho(iglob)*potential_dphi_dx_DG(ibool(i,j,ispec)) & ! \rho'g_x
          !                           + cv_rho0dv(NDIM,iglob)*nabla_v0(1,NDIM,iglob)) * jacLoc ! \rho_0v'_z\partial_zv_{0,x}
          !d0cntrb_rho0dv(NDIM,i,j) = - cv_drho(iglob)*potential_dphi_dz_DG(ibool(i,j,ispec)) * jacLoc ! \rho'g_z
          ! OLD4.2.3) Version 3: under HV0 (v_{0,z}=0), potential_dphi_dx_DG=0, and potential_dphi_dz_DG=LNS_g
          !d0cntrb_rho0dv(1,i,j) = - cv_rho0dv(NDIM,iglob)*nabla_v0(1,NDIM,iglob) * jacLoc ! \rho_0v'_z\partial_zv_{0,x}
          !d0cntrb_rho0dv(NDIM,i,j) = - cv_drho(iglob)*LNS_g(iglob) * jacLoc ! \rho'g_z
          ! OLD4.2.4) Version 4: under HV0 (v_{0,z}=0), SM (d_xv_0=0), potential_dphi_dx_DG=0, and potential_dphi_dz_DG=LNS_g
          !d0cntrb_rho0dv(1,i,j)    = - cv_rho0dv(NDIM,iglob)*nabla_v0(1,NDIM,iglob) * jacLoc !  \rho_0v'_z\partial_zv_{0,x}
          !d0cntrb_rho0dv(NDIM,i,j) = - cv_drho(iglob)*LNS_g(iglob) * jacLoc ! -\rho'g_z, but LNS_g=9.81=-g_z
          d0cntrb_rho0dv(:,i,j) = d0cntrb_rho0dv(:,i,j) * jacLoc ! small optimisation
          ! 4.3) Energy.
          ! 4.2.1) Version 1: most general.
          !d0cntrb_dE(i,j) = - DOT_PRODUCT(locG, in_dm(:,iglob)) * jacLoc
          ! 4.2.2) Version 2: [gravity is vertical (locG(1)=0, locG(2)=-LNS_g)].
          d0cntrb_dE(i,j) = LNS_g(iglob)*in_dm(NDIM,iglob) * jacLoc
          
!          ! 5) Eventually, PML.
! 
!          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!          ! IMPORTANT INFORMATION:
!          ! These blocks of code are commented out to prevent confusion.
!          ! But they may be useful in the future, when someone actually tries to finish the full PML implementation.
!          ! Please do not remove those lines, even if they are commented out for now.
!          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
!          if(PML_BOUNDARY_CONDITIONS .and. ispec_is_PML(ispec)) then
!            ! Do two things here:
!            ! 1) Add auxiliary variables to main RHS.
!            ! 2) Prepare auxiliary variables RHS.
!            
!            pml_a1 = product(LNS_PML_kapp(:,i,j,ispec_PML))
!            pml_a0 = LNS_PML_a0(i,j,ispec_PML) ! LNS_PML_a0 initialised in 'compute_forces_acoustic_LNS_calling_routine.F90'.
!            pml_alp = LNS_PML_alpha(:,i,j,ispec_PML)
!            pml_b  = LNS_PML_b(:,i,j,ispec_PML) ! LNS_PML_b initialised in 'compute_forces_acoustic_LNS_calling_routine.F90'.
!            pml_boa = ZEROcr
!            where(pml_b/=ZEROcr) pml_boa = pml_b / pml_alp ! b_i/alpha_i, only where b_i!=0
!            !write(*,*) "pml_a1, pml_a0, pml_alp, pml_b, pml_boa", pml_a1, pml_a0, pml_alp, pml_b, pml_boa! debuG
!            !stop 'kek'
!            
!            ! 1) Add auxiliary variables to main RHS.
!            ! 1.1) Update degree 0 contributions (set G=a1*G ***FIRST***), add a_0q, add G auxiliary variables, and add YU's auxiliary variables.
!            call LNS_PML_updateD0(d0cntrb_drho(i,j), cv_drho(iglob), 1, &
!                                  pml_a1, pml_a0, pml_boa, pml_b, jacLoc, iglobPML)
!            do k = 1, NDIM
!              call LNS_PML_updateD0(d0cntrb_rho0dv(k,i,j), cv_rho0dv(k,iglob), 1+k, &
!                                    pml_a1, pml_a0, pml_boa, pml_b, jacLoc, iglobPML)
!            enddo
!            call LNS_PML_updateD0(d0cntrb_dE(i,j), cv_dE(iglob), 2+NDIM, &
!                                  pml_a1, pml_a0, pml_boa, pml_b, jacLoc, iglobPML)
!            ! 1.2) Update degree 1 contributions.
!            ! 1.2.1) First, reset xigammal.
!            xigammal(1,    1)    = xix(i,j,ispec) ! = \partial_x\xi from report
!            xigammal(1,    NDIM) = xiz(i,j,ispec) ! = \partial_z\xi from report
!            xigammal(NDIM, 1)    = gammax(i,j,ispec) ! = \partial_x\eta from report
!            xigammal(NDIM, NDIM) = gammaz(i,j,ispec) ! = \partial_z\eta from report
!            ! 1.2.2) Need to include d_i on auxiliary stress tensor \Sigma=(R^{alpha_2}_{\sigma_1}, R^{alpha_1}_{\sigma_2}). We do it as follows. While not very readable, it gets the job done in a somewhat efficient way. \Sigma_1 needs to become d_2\Sigma_1, and \Sigma_2 needs to become d_1\Sigma_2. Since all \Sigma_i contributions are added through the dot product with (xix,xiz,gammax,gammaz), we include d_i in those directly.
!            xigammal(:, 1)    = LNS_PML_d(NDIM, i,j,ispec_PML)*xigammal(:, 1)    ! Multiply x component by d_2.
!            xigammal(:, NDIM) = LNS_PML_d(1,    i,j,ispec_PML)*xigammal(:, NDIM) ! Multiply z component by d_1.
!            ! 1.2.3)
!            do SPCDM = 1, NDIM ! rho'
!              cntrb_drho(i,j,SPCDM) =   cntrb_drho(i,j,SPCDM) &
!                                      + wzlwxljacLoc(SPCDM)*DOT_PRODUCT(xigammal(SPCDM,:), LNS_PML(1,3:4,iglobPML))
!            enddo
!            do k = 1, NDIM ! rho0v'
!              do SPCDM = 1, NDIM
!                cntrb_rho0dv(k,i,j,SPCDM) =   cntrb_rho0dv(k,i,j,SPCDM) &
!                                            + wzlwxljacLoc(SPCDM)*DOT_PRODUCT(xigammal(SPCDM,:), LNS_PML(1+k,3:4,iglobPML))
!              enddo
!            enddo
!            do SPCDM = 1, NDIM ! E'
!              cntrb_dE(i,j,SPCDM) =   cntrb_dE(i,j,SPCDM) &
!                                    + wzlwxljacLoc(SPCDM)*DOT_PRODUCT(xigammal(SPCDM,:), LNS_PML(2+NDIM,3:4,iglobPML))
!            enddo
!            
!            ! 2) Prepare auxiliary variables RHS. Recall pattern for R^beta_q: rhs = beta*R^beta_q - q.
!            !if(.false.) then
!            ! 2.1) YU*q (:,1:2,:)
!            ! 2.1.1) rho' (1,:,:)
!            call LNS_PML_buildRHS(1, 1, iglobPML, pml_alp(1), cv_drho(iglob))
!            call LNS_PML_buildRHS(1, 2, iglobPML, pml_alp(2), cv_drho(iglob))
!            ! 2.1.2) rho0v' (2:3,:,:)
!            do k=1,NDIM
!              call LNS_PML_buildRHS(1+k, 1, iglobPML, pml_alp(1), cv_rho0dv(k, iglob))
!              call LNS_PML_buildRHS(1+k, 2, iglobPML, pml_alp(2), cv_rho0dv(k, iglob))
!            enddo
!            ! 2.1.3) E' (2+NDIM,:,:)
!            call LNS_PML_buildRHS(2+NDIM, 1, iglobPML, pml_alp(1), cv_dE(iglob))
!            call LNS_PML_buildRHS(2+NDIM, 2, iglobPML, pml_alp(2), cv_dE(iglob))
!            ! 2.2) Sigma (:,3:4,:) : care for indices flipping (alpha2 with sigma1 and alpha1 with sigma2) ! TODO: TAKE CARE OF VISCOUS TENSOR
!            ! 2.2.1) rho' (1,:,:)
!            call LNS_PML_buildRHS(1, 3, iglobPML, pml_alp(2), in_dm(1,iglob))
!            call LNS_PML_buildRHS(1, 4, iglobPML, pml_alp(1), in_dm(2,iglob))
!            ! 2.2.2) rho0v' (2:3,:,:)
!            call LNS_PML_buildRHS(2, 3, iglobPML, pml_alp(2), cv_rho0dv(1,iglob)   *LNS_v0(1,iglob)   +in_dp(iglob)) ! vx sig1
!            call LNS_PML_buildRHS(2, 4, iglobPML, pml_alp(1), cv_rho0dv(NDIM,iglob)*LNS_v0(1,iglob)) ! vx sig2
!            call LNS_PML_buildRHS(3, 3, iglobPML, pml_alp(2), cv_rho0dv(1,iglob)   *LNS_v0(NDIM,iglob)) ! vz sig1
!            call LNS_PML_buildRHS(3, 4, iglobPML, pml_alp(1), cv_rho0dv(NDIM,iglob)*LNS_v0(NDIM,iglob)+in_dp(iglob)) ! vz sig2
!            ! 2.2.3) E' (2+NDIM,:,:)
!            call LNS_PML_buildRHS(2+NDIM, 3, iglobPML, pml_alp(2), &
!                                     LNS_dv(1,iglob)*(LNS_E0(iglob)+LNS_p0(iglob)) &
!                                   + LNS_v0(1,iglob)*(cv_dE(iglob)+in_dp(iglob))) ! E sig1
!            call LNS_PML_buildRHS(2+NDIM, 4, iglobPML, pml_alp(1), &
!                                     LNS_dv(2,iglob)*(LNS_E0(iglob)+LNS_p0(iglob)) &
!                                   + LNS_v0(2,iglob)*(cv_dE(iglob)+in_dp(iglob))) ! E sig2
!            ! 2.3) G (:,5:6,:) : zeroth degree contribution
!            ! 2.3.1) rho' (1,:,:): G=0
!            !call LNS_PML_buildRHS(1, 5, iglobPML, pml_alp(1), ZEROcr) ! rho a1G, do nothing
!            !call LNS_PML_buildRHS(1, 6, iglobPML, pml_alp(2), ZEROcr) ! rho a2G, do nothing
!            LNS_PML_RHS(1,5:6,iglobPML) = ZEROcr ! force to zero
!            ! 2.3.2) rho0v' (2:3,:,:): G=...
!            call LNS_PML_buildRHS(2, 5, iglobPML, pml_alp(1), cv_rho0dv(NDIM,iglob)*nabla_v0(1,NDIM,iglob)) ! vx a1G
!            call LNS_PML_buildRHS(2, 6, iglobPML, pml_alp(2), cv_rho0dv(NDIM,iglob)*nabla_v0(1,NDIM,iglob)) ! vx a2G
!            call LNS_PML_buildRHS(3, 5, iglobPML, pml_alp(1), - cv_drho(iglob)*LNS_g(iglob)) ! vz a1G
!            call LNS_PML_buildRHS(3, 6, iglobPML, pml_alp(2), - cv_drho(iglob)*LNS_g(iglob)) ! vz a2G
!            ! 2.3.3) E' (2+NDIM,:,:): G=g*dm_z
!            call LNS_PML_buildRHS(2+NDIM, 5, iglobPML, pml_alp(1), - LNS_g(iglob)*in_dm(NDIM,iglob)) ! vx a1G
!            call LNS_PML_buildRHS(2+NDIM, 6, iglobPML, pml_alp(2), - LNS_g(iglob)*in_dm(NDIM,iglob)) ! vx a2G
!            !endif
!          else
!            d0cntrb_drho = ZEROcr ! Safeguard: make sure d0cntrb_drho is zero when not in PMLs or when not using PMLs at all.
!          endif ! Endif on PML_BOUNDARY_CONDITIONS.
        enddo ! Enddo on i.
      enddo ! Enddo on j.
      
      ! Assemble the contributions previously computed, and add gravity's contribution.
      ! The integration by quadrature on the GLL points leads to three sums. See in particular Komatitsch (Méthodes spectrales et éléments spectraux pour l'équation de l'élastodynamique 2D et 3D en milieu hétérogène), Annexe 3.A.
      do j = 1, NGLLZ
        do i = 1, NGLLX
          iglob = ibool_DG(i,j,ispec)
          do k = 1, NGLLX
            outrhs_drho(iglob) =   outrhs_drho(iglob) &
                                + (  cntrb_drho(k,j,1) * real(hprimewgll_xx(k,i), kind=CUSTOM_REAL) &
                                   + cntrb_drho(i,k,2) * real(hprimewgll_zz(k,j), kind=CUSTOM_REAL))
            do SPCDM = 1, NDIM
              outrhs_rho0dv(SPCDM,iglob) =   outrhs_rho0dv(SPCDM,iglob) &
                                          + (  cntrb_rho0dv(SPCDM,k,j,1) * real(hprimewgll_xx(k,i), kind=CUSTOM_REAL) &
                                             + cntrb_rho0dv(SPCDM,i,k,2) * real(hprimewgll_zz(k,j), kind=CUSTOM_REAL))
            enddo
            outrhs_dE(iglob) =   outrhs_dE(iglob) &
                              + (  cntrb_dE(k,j,1) * real(hprimewgll_xx(k,i), kind=CUSTOM_REAL) &
                                 + cntrb_dE(i,k,2) * real(hprimewgll_zz(k,j), kind=CUSTOM_REAL))
          enddo ! Enddo on k.
          
          wxlwzl = real(wxgll(i)*wzgll(j), kind=CUSTOM_REAL)
          
          ! Add zero-th order terms.
! 
!          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!          ! IMPORTANT INFORMATION:
!          ! The following block of code is commented out to prevent confusion.
!          ! But it might be useful in the future, when someone actually tries to finish the full PML implementation.
!          ! Please do not remove those lines, even if they are commented out for now.
!          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
!          if(PML_BOUNDARY_CONDITIONS .and. anyabs .and. ispec_is_PML(ispec)) then ! Degree 0 contribution to rho' only exists if there are PMLs. TODO: a better condition.
!            outrhs_drho(iglob) = outrhs_drho(iglob) + d0cntrb_drho(i,j) * wxlwzl
!          endif
          do SPCDM = 1, NDIM
            outrhs_rho0dv(SPCDM,iglob) = outrhs_rho0dv(SPCDM,iglob) + d0cntrb_rho0dv(SPCDM,i,j) * wxlwzl
          enddo
          outrhs_dE(iglob) = outrhs_dE(iglob) + d0cntrb_dE(i,j) * wxlwzl
        enddo
      enddo
      
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
          iglobM = ibool_DG(i,j,ispec)
          
          if(LNS_viscous) then ! Check if viscosity exists whatsoever.
            ! Activate/deactivate, for this particular point (iglob), computation of quantities only needed when viscosity is present.
            if(     LNS_mu(iglobM) > TINYVAL &
               .OR. LNS_eta(iglobM) > TINYVAL &
               .OR. LNS_kappa(iglobM) > TINYVAL) then
              viscousComputation=.true.
            else
              viscousComputation=.false.
            endif
          else
            viscousComputation=.false. ! If viscosity is globally disabled, deactivate it for this element.
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
            ispec_neighbor = neighbor_DG_iface(iface1, iface, ispec, 3)
            neighbor(1:2) = link_iface_ijispec(iface1_neighbor, iface_neighbor, ispec_neighbor,1:2)
            neighbor(3) = ispec_neighbor
          endif
          
          ! Step 2: knowing the normals' parameters, compute now the fluxes.
          iglobP = 1
          if(neighbor(1) > -1) then
            iglobP = ibool_DG(neighbor(1), neighbor(2), neighbor(3))
          endif
          exact_interface_flux = .false. ! Reset this variable to .false.: by default, the fluxes have to be computed (jump!=0). In some specific cases (assigned during the call to LNS_get_interfaces_unknowns), the flux can be exact (jump==0).
          call LNS_get_interfaces_unknowns(i, j, ispec, iface1, iface, neighbor, currentTime, & ! Point identifier (input).
                  cv_drho(iglobM), cv_rho0dv(:,iglobM), & ! Input constitutive variables, "M" side.
                  cv_drho(iglobP), cv_rho0dv(:,iglobP), cv_dE(iglobP), & ! Input constitutive variables, "P" side. Note they make no sense if neighbor(1)<=-1.
                  in_dp(iglobM), & ! Input other variable, "M" side.
                  !V_DG(:,:,iglobM), T_DG(:,iglobM), & ! Input derivatives, "M" side. MIGHT NEED.
                  !V_DG(:,:,iglobP), T_DG(:,iglobP), & ! Input derivatives, "M" side. MIGHT NEED.
                  n_out, & ! Normal vector (input).
                  exact_interface_flux, & ! Switch to disable jump in some cases (output).
                  drho_P, rho0dv_P, dE_P, & ! Output constitutive variables.
                  !Tx_DG_P, Tz_DG_P, Vxx_DG_P, Vzz_DG_P, Vzx_DG_P, Vxz_DG_P,& ! Output derivatives. MIGHT NEED.
                  dm_P, dp_P, dv_P,& ! Output other variables.
                  viscousComputation, nabla_dT_P, sigma_dv_P, & ! Output other variables: viscous.
                  .false., LNS_dummy_1d(1)) ! Use dummies and set the switch to false not to compute unecessary quantities.
          
          jump   = ZERO ! Safeguard.
          
          ! Recover an approximate local maximum linearized acoustic wave speed. See for example Hesthaven (doi.org/10.1007/9780387720678), page 208.
          lambda = ZERO
          !lambda = max(  abs(dot_product(n_out, LNS_v0(:,iglobM))) & ! v_-\cdot n
          !             + LNS_c0(iglobM) & ! Local sound speed, side "M".
          !             , abs(dot_product(n_out, LNS_v0(:,iglobP))) & ! v_+\cdot n
          !             + LNS_c0(iglobP) & ! Local sound speed, side "P".
          !             )
          ! Exact lambda (i.e. with exact full v=v0+v' and c=c0+c').
          ! This may be important for manufactured solutions' validation with non-zero (non-negligible) dv field.
          lambda = max(  abs(dot_product(n_out, LNS_v0(:,iglobM)+LNS_dv(:,iglobM))) &
                       + sqrt(gammaext_DG(iglobM)*(LNS_p0(iglobM)+in_dp(iglobM))/(LNS_rho0(iglobM)+cv_drho(iglobM))) &
                       , abs(dot_product(n_out, LNS_v0(:,iglobP)+LNS_dv(:,iglobP))) &
                       + sqrt(gammaext_DG(iglobP)*(LNS_p0(iglobP)+dp_P)/(LNS_rho0(iglobP)+drho_P)) &
                       )
          
          ! For real stretching, \Sigma for each constitutive variable becomes \Ya\Sigma. It is heavy to change each and every expression where \Sigma arises. Rather, we make use of what is multiplying \Sigma.
          ! Here, in the surface integrations, easiest is the normals. But we have to do it after the call to LNS_get_interfaces_unknowns and the affectation of lambda. If you change anything, remember to do it also in the 'compute_gradient_TFSF' subroutine.
          if(ABC_STRETCH .and. stretching_buffer(ibool_before_perio(i,j,ispec))>0) then
            iglob_unique = ibool_before_perio(i,j,ispec)
            do SPCDM = 1, NDIM
              n_out(SPCDM) = n_out(SPCDM) * stretching_ya(SPCDM, iglob_unique)
            enddo
          endif
          
          ! 1) Inviscid contributions.
          ! 1.1) Mass conservation (fully inviscid).
          !flux_n = DOT_PRODUCT(n_out, in_dm(:,iglobM)+dm_P)
          if(exact_interface_flux) then
            jump = ZERO
          else
            jump = cv_drho(iglobM) - drho_P
          endif
          !outrhs_drho(iglobM) = outrhs_drho(iglobM) - halfWeight*(flux_n + lambda*jump) ! Add flux' contribution.
          outrhs_drho(iglobM) =   outrhs_drho(iglobM) &
                                - halfWeight*(DOT_PRODUCT(n_out, in_dm(:,iglobM)+dm_P) + lambda*jump)
          ! 1.2) x-Momentum inviscid contributions.
          locSigma(1)    =   cv_rho0dv(1,iglobM)*LNS_v0(1,iglobM) + in_dp(iglobM) & ! "M" side.
                           + rho0dv_P(1)        *LNS_v0(1,iglobP) + dp_P ! "P" side.
          locSigma(NDIM) =   cv_rho0dv(NDIM,iglobM)*LNS_v0(1,iglobM) & ! "M" side.
                           + rho0dv_P(NDIM)        *LNS_v0(1,iglobP) ! "P" side.
          !flux_n = DOT_PRODUCT(n_out, locSigma)
          if(exact_interface_flux) then
            jump = ZERO
          else
            jump = cv_rho0dv(1,iglobM) - rho0dv_P(1)
          endif
          !outrhs_rho0dv(1,iglobM) = outrhs_rho0dv(1,iglobM) - halfWeight*(flux_n + lambda*jump) ! Add flux' contribution.
          outrhs_rho0dv(1,iglobM) =   outrhs_rho0dv(1,iglobM) &
                                    - halfWeight*(DOT_PRODUCT(n_out, locSigma) + lambda*jump)
          ! 1.2) z-Momentum inviscid contributions.
          locSigma(1)    =   cv_rho0dv(1,iglobM)*LNS_v0(NDIM,iglobM) & ! "M" side.
                           + rho0dv_P(1)        *LNS_v0(NDIM,iglobP) ! "P" side.
          locSigma(NDIM) =   cv_rho0dv(NDIM,iglobM)*LNS_v0(NDIM,iglobM) + in_dp(iglobM) & ! "M" side.
                           + rho0dv_P(NDIM)        *LNS_v0(NDIM,iglobP) + dp_P ! "P" side.
          !flux_n = DOT_PRODUCT(n_out, locSigma)
          if(exact_interface_flux) then
            jump = ZERO
          else
            jump = cv_rho0dv(NDIM,iglobM) - rho0dv_P(NDIM)
          endif
          !outrhs_rho0dv(NDIM,iglobM) = outrhs_rho0dv(NDIM,iglobM) - halfWeight*(flux_n + lambda*jump) ! Add flux' contribution.
          outrhs_rho0dv(NDIM,iglobM) =   outrhs_rho0dv(NDIM,iglobM) &
                                       - halfWeight*(DOT_PRODUCT(n_out, locSigma) + lambda*jump)
          ! 1.3) Energy inviscid contributions.
          locSigma =   LNS_dv(:,iglobM)*(LNS_E0(iglobM) + LNS_p0(iglobM)) & ! "M" side, part.
                     + LNS_v0(:,iglobM)*(cv_dE(iglobM)  + in_dp(iglobM)) & ! "M" side, part.
                     + dv_P(:)         *(LNS_E0(iglobP) + LNS_p0(iglobP)) & ! "P" side, part.
                     + LNS_v0(:,iglobP)*(dE_P           + dp_P) ! "P" side, part.
          !flux_n = DOT_PRODUCT(n_out, locSigma)
          if(exact_interface_flux) then
            jump = ZERO
          else
            jump = cv_dE(iglobM) - dE_P
          endif
          !outrhs_dE(iglobM) = outrhs_dE(iglobM) - halfWeight*(flux_n + lambda*jump) ! Add flux' contribution.
          outrhs_dE(iglobM) =   outrhs_dE(iglobM) &
                              - halfWeight*(DOT_PRODUCT(n_out, locSigma) + lambda*jump)
          
          ! 2) Viscous contributions.
          if(viscousComputation) then
            ! Note: since we consider here viscous terms, that is purely diffusive phenomena, the acoustic wave speed makes no sense, and thus the jump in the Lax-Friedrich approximation does not appear in the formulations.
            ! 2.1) Mass conservation: nothing.
            ! 2.2) x-Momentum viscous contributions. The vector [locSigma_x, locSigma_z] represents the mean average flux at the boundary of the x-momentum (based on the 1st line of \Sigma_v).
            locSigma(1)    = -in_sigma_dv(1,iglob) -sigma_dv_P(1) ! Recall, this corresponds to \Sigma_v(v')_{1,1}.
            locSigma(NDIM) = -in_sigma_dv(2,iglob) -sigma_dv_P(2) ! Recall, this corresponds to \Sigma_v(v')_{1,2} (and \Sigma_v(v')_{2,1}).
            outrhs_rho0dv(1,iglobM) = outrhs_rho0dv(1,iglobM) - halfWeight*DOT_PRODUCT(n_out, locSigma) ! Add flux' contribution, with dot product \Sigma\cdot n.
            ! 2.2) z-Momentum viscous contributions. The vector [locSigma_x, locSigma_z] represents the mean average flux at the boundary of the z-momentum (based on the 2nd line of \Sigma_v).
            locSigma(1)    = -in_sigma_dv(2,iglob) -sigma_dv_P(2) ! Recall, this corresponds to \Sigma_v(v')_{2,1} (and \Sigma_v(v')_{1,2}).
            locSigma(NDIM) = -in_sigma_dv(3,iglob) -sigma_dv_P(3) ! Recall, this corresponds to \Sigma_v(v')_{2,2}.
            outrhs_rho0dv(NDIM,iglobM) = outrhs_rho0dv(NDIM,iglobM) - halfWeight*DOT_PRODUCT(n_out, locSigma) ! Add flux' contribution, with dot product \Sigma\cdot n.
            ! 2.3) Energy viscous contributions.
            locSigma(1)    = - (  LNS_dv(1,iglobM)*sigma_v_0(1,iglobM) & ! "M" side, part.
                                + LNS_v0(1,iglobM)*in_sigma_dv(1,iglobM) & ! "M" side, part.
                                + LNS_dv(NDIM,iglobM)*sigma_v_0(2,iglobM) & ! "M" side, part.
                                + LNS_v0(NDIM,iglobM)*in_sigma_dv(2,iglobM) & ! "M" side, part.
                                + LNS_kappa(iglobM)*in_nabla_dT(1,iglobM)) & ! "M" side, part.
                             - (  dv_P(1)*sigma_v_0(1,iglobP) & ! "P" side, part.
                                + LNS_v0(1,iglobM)*sigma_dv_P(1) & ! "P" side, part.
                                + dv_P(NDIM)*sigma_v_0(2,iglobP) & ! "P" side, part.
                                + LNS_v0(NDIM,iglobP)*sigma_dv_P(2) & ! "P" side, part.
                                + LNS_kappa(iglobP)*nabla_dT_P(1)) ! "P" side, part.
            locSigma(NDIM) = - (  LNS_dv(NDIM,iglobM)*sigma_v_0(3,iglobM) & ! "M" side, part.
                                + LNS_v0(NDIM,iglobM)*in_sigma_dv(3,iglobM) & ! "M" side, part.
                                + LNS_dv(1,iglobM)*sigma_v_0(2,iglobM) & ! "M" side, part.
                                + LNS_v0(1,iglobM)*in_sigma_dv(2,iglobM) & ! "M" side, part.
                                + LNS_kappa(iglobM)*in_nabla_dT(2,iglobM)) & ! "M" side, part.
                             - (  dv_P(NDIM)*sigma_v_0(3,iglobP) & ! "P" side, part.
                                + LNS_v0(NDIM,iglobP)*sigma_dv_P(3) & ! "P" side, part.
                                + dv_P(1)*sigma_v_0(2,iglobP) & ! "P" side, part.
                                + LNS_v0(1,iglobP)*sigma_dv_P(2) & ! "P" side, part.
                                + LNS_kappa(iglobP)*nabla_dT_P(NDIM)) ! "P" side, part.
            outrhs_dE(iglobM) = outrhs_dE(iglobM) + halfWeight*DOT_PRODUCT(n_out, locSigma)
          endif ! Endif on viscousComputation.
        enddo ! Enddo on iface.
      enddo ! Enddo on iface1.
    endif ! End of test if acoustic element
  enddo ! End of loop on elements.
  
end subroutine compute_forces_acoustic_LNS




















! ------------------------------------------------------------ !
! LNS_get_interfaces_unknowns                                  !
! ------------------------------------------------------------ !
! From the coordinates of a GLL point in an element (local coordinates (i, j) and global element number ispec), and its neighbour's identifier (neighbor), compute the values of the constitutive variables at the neighbour.
! Variables ending in "out_*_P" (for "plus") are output-intended and are the sought values, or exterior values.
! Variables ending in "inp_*_M" (for "minus") are input-intended and should correspond to the interior values.
! Variables ending in "inp_*_P" are input-intended, and correspond to the exterior values, if known. Remark that those variables are only used if neighbor(3) != -1.
! n_out: outward-pointing normal vector.
! exact_interface_flux: switch to disable jump in some cases.
! drho_P: \rho', for "P" side.
! rho0dv_P: \rho_0v', for "P" side.
! dE_P: E', for "P" side.
! dm_P: \rho_0v'+\rho'v_0, for "P" side. It is needed only for the mass conservation equation.
! dv_P: v', for "P" side. It is needed only for the energy equation, both for inviscid and viscous contributions.
! nabla_dT_P: \nabla T', for "P" side. It is needed only for the energy equation viscous contribution.
! sigma_dv_P: \Sigma_v(v'), for "P" side. It is (obviously) needed only for viscous contributions, both for momenta and energy equations.

!subroutine LNS_get_interfaces_unknowns(i, j, ispec, neighbor, & ! Point identifier.
!                  veloc_x_DG_P, veloc_z_DG_P, in_dp_P, T_P, & ! Precomputed quantities.
!                  exact_interface_flux, &
!                  inp_drho_M, inp_rho0dv_M, inp_dE_M, & ! Input constitutive variables, "M" side.
!                  inp_drho_P, inp_dE_P, inp_rho0dv_P, & ! Input constitutive variables, "P" side.
!                  !V_DG_iM, T_DG_iM, & ! Input derivatives, "M" side. MIGHT NEED.
!                  !V_DG_iP, T_DG_iP, & ! Input derivatives, "P" side. MIGHT NEED.
!                  nx, nz, weight, currentTime, iface1, iface&
!                  
!                  out_drho_P, out_rho0dv_P, out_dE_P& ! Output constitutive variables.
!                  !Tx_DG_P, Tz_DG_P, Vxx_DG_P, Vzz_DG_P, Vzx_DG_P, Vxz_DG_P,& ! Output derivatives. MIGHT NEED.
!                  out_dv_P) ! Output other variables.
subroutine LNS_get_interfaces_unknowns(i, j, ispec, iface1, iface, neighbor, timelocal, & ! Point identifier (input).
             inp_drho_M, inp_rho0dv_M, & ! Input constitutive variables, "M" side.
             inp_drho_P, inp_rho0dv_P, inp_dE_P, & ! Input constitutive variables, "P" side.
             inp_dp_M, & ! Input other variable, "M" side.
             !V_DG_iM, T_DG_iM, V_DG_iP, T_DG_iP, & ! Input derivatives. MIGHT NEED.
             n_out, & ! Normal vector (input).
             exact_interface_flux, & ! Switch to disable jump in some cases (output).
             out_drho_P, out_rho0dv_P, out_dE_P, & ! Output constitutive variables.
             !Tx_DG_P, Tz_DG_P, Vxx_DG_P, Vzz_DG_P, Vzx_DG_P, Vxz_DG_P,& ! Output derivatives. MIGHT NEED.
             out_dm_P, out_dp_P, out_dv_P, & ! Output other variables.
             swCompVisc, out_nabla_dT_P, out_sigma_dv_P, & ! Output other variables: viscous.
             swCompdT, out_dT_P) ! Output other variables.
  
  use constants ! TODO: select variables to use.,only: CUSTOM_REAL,NGLLX,NGLLZ,gamma_euler
  use specfem_par ! TODO: select variables to use.,only:  ibool_DG, &
                         !ispec_is_acoustic_forcing, &
                         !ACOUSTIC_FORCING, &
                         ! ispec_is_acoustic_coupling_el, ispec_is_acoustic_coupling_ac, potential_dot_dot_acoustic, veloc_elastic,&
                         !DIR_RIGHT, DIR_LEFT, DIR_UP, DIR_DOWN, &
                         !buffer_DG_rho_P, buffer_DG_rhovx_P, buffer_DG_rhovz_P, buffer_DG_E_P, NPROC, &
                         !buffer_DG_Vxx_P, buffer_DG_Vzz_P, buffer_DG_Vxz_P, buffer_DG_Vzx_P, buffer_DG_Tz_P, buffer_DG_Tx_P, &
                         !LNS_p0, gammaext_DG, muext, etaext, kappa_DG, ibool, c_V, &
                         !buffer_DG_gamma_P, coord,  &
                         !rho_init, rhovx_init, rhovz_init, E_init, &
                         !potential_dot_dot_acoustic, &
                         !veloc_vector_acoustic_DG_coupling, MPI_transfer_iface,&
                         !ibool_before_perio, ABC_STRETCH, ABC_STRETCH_LEFT,&
                         !ABC_STRETCH_RIGHT, ABC_STRETCH_LEFT_LBUF, ABC_STRETCH_RIGHT_LBUF,&
                         !mesh_xmin, mesh_xmax
  use specfem_par_LNS ! TODO: select variables to use.
  
  implicit none
  
  ! Input/Output.
  integer, intent(in) :: i, j, ispec, iface1, iface
  integer, dimension(3), intent(in) :: neighbor
  real(kind=CUSTOM_REAL), intent(in) :: timelocal
  real(kind=CUSTOM_REAL), intent(in) :: inp_drho_M, inp_drho_P, inp_dE_P ! Input constitutive variables.
  real(kind=CUSTOM_REAL), dimension(NDIM), intent(in) :: inp_rho0dv_M, inp_rho0dv_P ! Input constitutive variables.
  real(kind=CUSTOM_REAL), intent(in) :: inp_dp_M ! Input other variables.
  !real(kind=CUSTOM_REAL), dimension(NDIM), intent(in) :: T_DG_iP, T_DG_iM ! Input derivatives. MIGHT NEED.
  !real(kind=CUSTOM_REAL), dimension(NDIM, NDIM), intent(in) :: V_DG_iP, V_DG_iM ! Input derivatives. MIGHT NEED.
  real(kind=CUSTOM_REAL), dimension(NDIM), intent(in) :: n_out
  logical, intent(out) :: exact_interface_flux ! Output switch.
  real(kind=CUSTOM_REAL), intent(out) :: out_drho_P, out_dE_P ! Output constitutive variables.
  !real(kind=CUSTOM_REAL), intent(out) :: Tx_DG_P, Tz_DG_P, Vxx_DG_P, Vzz_DG_P, Vzx_DG_P, Vxz_DG_P ! Output derivatives. MIGHT NEED.
  real(kind=CUSTOM_REAL), dimension(NDIM), intent(out) :: out_rho0dv_P ! Output constitutive variable.
  real(kind=CUSTOM_REAL), dimension(NDIM), intent(out) :: out_dv_P, out_dm_P, out_nabla_dT_P ! Output other variables.
  real(kind=CUSTOM_REAL), intent(out) :: out_dp_P ! Output other variables.
  real(kind=CUSTOM_REAL), dimension(NVALSIGMA), intent(out) :: out_sigma_dv_P ! Output other variables.
  real(kind=CUSTOM_REAL), intent(out) :: out_dT_P ! In compute_gradient_TFSF (desintegration method), we need temperature on the other side of the boundary in order to compute the flux.
  !real(kind=CUSTOM_REAL), dimension(NDIM), intent(out) :: velocity_P ! In compute_gradient_TFSF (desintegration method), we need velocity on the other side of the boundary in order to compute the flux.
  logical, intent(in) :: swCompVisc, swCompdT!, swCompv ! Do not unnecessarily compute some quantities.
  
  ! Local.
  integer :: SPCDM
  real(kind=CUSTOM_REAL), dimension(NDIM) :: velocity_P
  real(kind=CUSTOM_REAL), dimension(NDIM) :: tang ! Tangential vector.
  real(kind=CUSTOM_REAL) :: normal_v, tangential_v!, &
                            !veloc_x, veloc_z!, gamma_P!, &
                            !dv_M(1), dv_M(NDIM), inp_dp_M, gamma_P, e1_DG_P
  real(kind=CUSTOM_REAL), dimension(NDIM) :: veloc_P
  real(kind=CUSTOM_REAL), parameter :: ZEROcr = 0._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: ONEcr  = 1._CUSTOM_REAL
  !real(kind=CUSTOM_REAL), parameter :: TWO  = 2._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: HALFcr = 0.5_CUSTOM_REAL
  integer :: iglobM, i_el, j_el, ispec_el, i_ac, j_ac, ispec_ac, iglob, iglobP, ipoin, num_interface
  real(kind=CUSTOM_REAL), dimension(NDIM, NDIM) :: trans_boundary
  !real(kind=CUSTOM_REAL) :: veloc_x_dg_p, veloc_z_dg_p!,x ! For coupling deactivation in buffers.
  ! Characteristic based BC
  !real(kind=CUSTOM_REAL) :: s_b, rho_inf, v_b_x, v_b_z, un_b, &!, rho_b, p_b, c_b, &
  !                          rlambda_max, rlambda_min, un_in, un_inf!, &!c_in, c_inf, &
                            !deltaZ1, deltaZ2star, p_n!, a_n, alpha0
  
  
  ! Initialise output variables to default values.
  exact_interface_flux=.false.
  out_drho_P=ZEROcr
  out_rho0dv_P=ZEROcr
  out_dE_P=ZEROcr
  out_dv_P=ZEROcr
  out_dm_P=ZEROcr
  out_dp_P=ZEROcr
  out_nabla_dT_P=ZEROcr
  out_sigma_dv_P=ZEROcr
  out_dT_P=ZEROcr
  velocity_P=ZEROcr
  
  ! Check each and every output variable is correctly set, and this whatever the case below.
  ! Set exact_interface_flux.
  ! Set out_drho_P.
  ! Set out_rho0dv_P.
  ! Set out_dE_P.
  ! Set out_dv_P.
  ! Set out_dm_P: see bottom of routine.
  ! Set out_dp_P.
  ! Set out_nabla_dT_P.
  ! Set out_sigma_dv_P.
  ! Set out_dT_P.
  
  ! Extract calling point iglob.
  iglobM = ibool_DG(i, j, ispec)
  
  !exact_interface_flux = .false. ! Unless otherwise specified, use the Lax-Friedrich approximation.
  
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
      
      if(ispec_is_acoustic_coupling_ac(ibool_DG(i, j, ispec)) >= 0) then
        ! --------------------------- !
        ! MPI acoustic potential      !
        ! neighbour.                  !
        ! --------------------------- !
        write(*,*) "********************************"
        write(*,*) "*            ERROR             *"
        write(*,*) "********************************"
        write(*,*) "* Interfaces between LNS and   *"
        write(*,*) "* acoustic potential elements  *"
        write(*,*) "* are not implemented yet.     *"
        write(*,*) "********************************"
        stop
        ! Set exact_interface_flux.
        exact_interface_flux = .true.
        ! Coordinates of elastic element
        iglob = ibool(i, j, ispec)
        ! Set out_drho_P.
        call background_physical_parameters(i, j, ispec, timelocal, out_drho_P, &
                                            .true., velocity_P(1), & ! "(1)" here is to please Intel compilers. It should not interfere too much since the concerned array is nicely allocated.
                                            .false., LNS_dummy_1d(1), &
                                            .false., LNS_dummy_1d(1)) ! Get needed background parameters. Use dummies for values we're not interested in.
        ! Set velocity_P: free slip and normal velocity continuity.
        !call boundary_condition_DG(i, j, ispec, timelocal, out_drho_P, out_rho0dv_P(1), out_rho0dv_P(NDIM), out_dE_P, &
        !                           out_dv_P(1), out_dv_P(NDIM), out_dp_P, e1_DG_P) ! Warning, expression of out_dp_P might not be exact.
        veloc_P = veloc_vector_acoustic_DG_coupling(iglob, :)
        !veloc_x = veloc_vector_acoustic_DG_coupling(iglob, 1)
        !veloc_z = veloc_vector_acoustic_DG_coupling(iglob, 2)
        call build_trans_boundary(n_out, tang, trans_boundary)
        ! Tangential vector
        ! Since only bottom topography n_out(NDIM) > 0
        !tang(1) = -n_out(NDIM) ! Recall: n_out(NDIM)=n_z.
        !tang(NDIM) = n_out(1) ! Recall: n_out(1)=n_x.
        ! Normal velocity of the solid perturbation and tangential velocity of the fluid flow
        !normal_v     = veloc_P(1)*n_out(1) + veloc_P(1)*n_out(NDIM) ! Recall: n_out(NDIM)=n_z.
        !tangential_v = veloc_P(1)*tang(1) + veloc_P(1)*tang(NDIM)
        normal_v     = DOT_PRODUCT(n_out, veloc_P)
        tangential_v = DOT_PRODUCT(tang, veloc_P)
        
        do SPCDM = 1, NDIM
          velocity_P(SPCDM) = trans_boundary(SPCDM, 1)*normal_v + trans_boundary(SPCDM, 2)*tangential_v!veloc_elastic(1,iglob)
          !velocity_P(1) = trans_boundary(1, 1)*normal_v + trans_boundary(1, 2)*tangential_v!veloc_elastic(1,iglob)
          !velocity_P(NDIM) = trans_boundary(2, 1)*normal_v + trans_boundary(2, 2)*tangential_v
        enddo
        
        ! TODO: ! QUICK HACK: DEACTIVATE COUPLING IN BUFFER ZONES. See boundary_terms_DG.f90.
        
        ! Set out_dv_P.
        out_dv_P = velocity_P - LNS_v0(:, iglobM) ! Requesting iglobM might be technically inexact, but on elements' boundaries points should overlap. Plus, iglobP does not exist on outer computational domain boundaries.
        ! Set out_dp_P: traction continuity.
        out_dp_P = LNS_p0(iglobM) - potential_dot_dot_acoustic(iglob) ! Warning, expression of out_dp_P might not be exact.
        ! Set out_dE_P.
        call compute_dE_i(LNS_rho0(iglobM)+out_drho_P, velocity_P, LNS_p0(iglobM)+out_dp_P, out_dE_P, iglobM) ! Warning, expression of out_dp_P might not be exact.
        !out_dE_P       = out_dp_P/(gammaext_DG(iglobM) - 1.) + HALFcr*out_drho_P*( out_dv_P(1)**2 + out_dv_P(NDIM)**2 )
        ! Set out_rho0dv_P.
        do SPCDM = 1, NDIM
          out_rho0dv_P(SPCDM) = out_drho_P*out_dv_P(SPCDM) ! Safe version, just in case element-wise mutiplication fails somehow.
        enddo
        !out_rho0dv_P=out_drho_P*out_dv_P
        ! Set out_dm_P: see bottom of routine.
        if(swCompVisc) then
          ! Set out_nabla_dT_P: same as other side, that is a Neumann condition.
          out_nabla_dT_P = nabla_dT(:,iglobM)
          ! Set out_sigma_dv_P: same as other side, that is a Neumann condition.
          out_sigma_dv_P = sigma_dv(:,iglobM)
        endif
        ! Set out_dT_P.
        if(swCompdT) then
          !out_dT_P = (inp_dE_M/inp_drho_M - 0.5*(dv_M(1)**2 + dv_M(NDIM)**2))/c_V
          call compute_dT_i(LNS_rho0(iglobM)+out_drho_P, velocity_P, LNS_E0(iglobM)+out_dE_P, out_dT_P, iglobM)
          !call compute_dT_i(LNS_rho0(iglobM)+out_drho_P, LNS_p0(iglobM)+out_dp_P, out_dT_P, iglobM)
        endif
        
      else
        ! --------------------------- !
        ! MPI acoustic DG neighbour   !
        ! (not acoustic potential     !
        ! neighbour).                 !
        ! --------------------------- !
        
        ! Set exact_interface_flux.
        exact_interface_flux = .false.
        
        ! Set out_drho_P.
        out_drho_P     = buffer_LNS_drho_P(ipoin, num_interface)
        
        ! Set out_rho0dv_P.
        out_rho0dv_P   = buffer_LNS_rho0dv_P(:, ipoin, num_interface)
        
        ! Set out_dE_P.
        out_dE_P       = buffer_LNS_dE_P(ipoin, num_interface)
        
        !gamma_P        = buffer_DG_gamma_P(ipoin,num_interface)
        
        ! Set out_dv_P.
        out_dv_P = out_rho0dv_P/LNS_rho0(iglobM)
        !out_dv_P(1) = out_rho0dv_P(1)/out_drho_P
        !out_dv_P(NDIM) = out_rho0dv_P(NDIM)/out_drho_P
        
        ! Set out_dp_P.
        call compute_dp_i(LNS_rho0(iglobM)+out_drho_P, LNS_v0(:,iglobM)+out_dv_P, LNS_E0(iglobM)+out_dE_P, out_dp_P, iglobM)
        !out_dp_P = (gamma_P - ONEcr)*( out_dE_P & ! Warning, expression of out_dp_P might not be exact.
        !         - (HALFcr)*out_drho_P*( out_dv_P(1)**2 + out_dv_P(NDIM)**2 ) )
        
        ! Set out_dm_P: see bottom of routine.
        
        if(swCompVisc) then
          ! Set out_nabla_dT_P: get the values from the MPI buffers.
          out_nabla_dT_P = buffer_LNS_nabla_dT(:, ipoin, num_interface)
          
          ! Set out_sigma_dv_P: get the values from the MPI buffers.
          out_sigma_dv_P = buffer_LNS_sigma_dv(:, ipoin, num_interface)
        endif
        
        ! Set out_dT_P.
        if(swCompdT) then
          call compute_dT_i(LNS_rho0(iglobM)+out_drho_P, LNS_v0(:, iglobM)+out_dv_P, LNS_E0(iglobM)+out_dE_P, out_dT_P, iglobM)
          !call compute_dT_i(LNS_rho0(iglobM)+out_drho_P, LNS_p0(iglobM)+out_dp_P, out_dT_P, iglobM)
        endif
        !out_dT_P = (out_dE_P/out_drho_P - 0.5*((out_rho0dv_P(1)/out_drho_P)**2 + (out_rho0dv_P(NDIM)/out_drho_P)**2))/c_V
        
      endif
      
    elseif(ACOUSTIC_FORCING .AND. ispec_is_acoustic_forcing(i, j, ispec)) then
      ! --------------------------- !
      ! ipoin == -1                 !
      !   and acoustic forcing      !
      ! --------------------------- !
      write(*,*) "********************************"
      write(*,*) "*            ERROR             *"
      write(*,*) "********************************"
      write(*,*) "* ACOUSTIC_FORCING is obsolete *"
      write(*,*) "* for DG simulations.          *"
      write(*,*) "********************************"
      stop
         
    elseif(ispec_is_acoustic_coupling_el(i, j, ispec, 3) >= 0) then
    !elseif(ispec_is_acoustic_coupling_el(i, j, ispec, 3) >= 0 .AND. abs(n_out(1)) < 1.) then ! The condition '|n_x|<1' is there to remove coupling on outer-pointing edges (those having |n_x|>1). However, while this works with flat horizontal topography, this also removes coupling at the vertical boundaries within the domain of interest.
      ! --------------------------- !
      ! ipoin == -1                 !
      !   and elastic coupling      !
      ! --------------------------- !
      
      ! Set exact_interface_flux.
      exact_interface_flux = .false.
      
      ! Coordinates of elastic element
      i_el     = ispec_is_acoustic_coupling_el(i, j, ispec, 1)
      j_el     = ispec_is_acoustic_coupling_el(i, j, ispec, 2)
      ispec_el = ispec_is_acoustic_coupling_el(i, j, ispec, 3)
      !iglob    = ibool(i_el, j_el, ispec_el)
      
      ! Set out_drho_P: same as other side, that is a Neumann condition.
      out_drho_P = inp_drho_M
      
      ! Set velocity_P.
      !! VERSION 1 : normal velocity continuity and slip condition for tangential
      !!veloc_P = veloc_elastic(:, iglob)
      !call build_trans_boundary(n_out, tang, trans_boundary)
      !!normal_v     = veloc_elastic(1, iglob)*n_out(1)  + veloc_elastic(NDIM, iglob)*n_out(NDIM)
      !normal_v     = DOT_PRODUCT(n_out, veloc_elastic(:,ibool(i_el, j_el, ispec_el)))
      !!tangential_v = veloc_x_DG_P*tang(1) + veloc_z_DG_P*tang(NDIM)
      !tangential_v = DOT_PRODUCT(tang, LNS_v0(:,iglobM)+LNS_dv(:,iglobM))
      !do SPCDM = 1, NDIM
      !  velocity_P(SPCDM) = trans_boundary(SPCDM, 1)*normal_v + trans_boundary(SPCDM, 2)*tangential_v
      !enddo
      ! VERSION 2: TERRANA.
      ! If you comment/uncomment this, don't forget to uncomment/comment the setting of out_dp_P below.
      call S2F_Terrana_coupling(n_out, &
                                LNS_rho0(iglobM)+inp_drho_M, &
                                LNS_v0(:,iglobM)+LNS_dv(:,iglobM), &
                                !LNS_p0(iglobM)+inp_dp_M, &
                                inp_dp_M, &
                                LNS_c0(iglobM), & ! either that, or recomputing using sqrt(gammaext_DG(iglobM)*(p0+dp)/rho_fluid)
                                iglobM, &
                                veloc_elastic(:,ibool(i_el, j_el, ispec_el)), &
                                sigma_elastic(:,:,ibool(i_el, j_el, ispec_el)), &
                                i_el, j_el, ispec_el, &
                                velocity_P, out_dp_P)
      
      ! Set out_dv_P.
      out_dv_P = velocity_P - LNS_v0(:, iglobM) ! Requesting iglobM might be technically inexact, but on elements' boundaries points should overlap. Plus, iglobP does not exist on outer computational domain boundaries.
      
      ! Set out_dp_P.
      ! VERSION 1: no stress continuity.
      !out_dp_P = inp_dp_M
      ! VERSION 2: TERRANA.
      ! done above in the call to S2F_Terrana_coupling

      ! Set out_dE_P.
      call compute_dE_i(LNS_rho0(iglobM)+out_drho_P, velocity_P, LNS_p0(iglobM)+out_dp_P, out_dE_P, iglobM)

      ! Set out_rho0dv_P.
      !do SPCDM = 1, NDIM
      !  out_rho0dv_P(SPCDM) = LNS_rho0(iglobM)*out_dv_P(SPCDM) ! Safe version, just in case element-wise mutiplication fails somehow.
      !enddo
      out_rho0dv_P(:)=LNS_rho0(iglobM)*out_dv_P(:)
      
      ! Set out_dm_P: see bottom of routine.
        
      if(swCompVisc) then
        ! Set out_nabla_dT_P: same as other side, that is a Neumann condition.
        out_nabla_dT_P = nabla_dT(:,iglobM)
        
        ! Set out_sigma_dv_P: same as other side, that is a Neumann condition.
        out_sigma_dv_P = sigma_dv(:,iglobM)
      endif
      
      ! Set out_dT_P.
      if(swCompdT) then
        call compute_dT_i(LNS_rho0(iglobM)+out_drho_P, velocity_P, LNS_E0(iglobM)+out_dE_P, out_dT_P, iglobM)
        !call compute_dT_i(LNS_rho0(iglobM)+out_drho_P, LNS_p0(iglobM)+out_dp_P, out_dT_P, iglobM)
      endif
      
      !if(      coord(2,ibool_before_perio(i,j,ispec))<1. & ! DEBUG
      !   .and. coord(2,ibool_before_perio(i,j,ispec))>=ZEROcr & ! DEBUG
      !   .and. abs(coord(1,ibool_before_perio(i,j,ispec)))<2.) then ! DEBUG
      !  write(*,*) timelocal, coord(:,ibool_before_perio(i,j,ispec)), & ! DEBUG
      !             !LNS_p0(iglobM)+inp_dp_M, LNS_p0(iglobM)+out_dp_P ! DEBUG
      !             !out_rho0dv_P ! DEBUG
      !             !out_dv_P ! DEBUG
      !             !trans_boundary ! DEBUG
      !endif ! DEBUG

    else
      ! --------------------------- !
      ! ipoin == -1                 !
      !   classical outer boundary  !
      !   conditions                !
      ! --------------------------- !
!      if(PML_BOUNDARY_CONDITIONS .and. anyabs .and. ispec_is_PML(ispec)) then
        !if(abs(coord(2, ibool_before_perio(i, j, ispec))-25.)<TINYVAL) & ! DEBUG
        !  write(*,*) "TOP PML OUTER BOUNDARY CONDITION IS ATTAINED" ! DEBUG
        ! --------------------------- !
        ! Outer PML.                  !
        ! --------------------------- !
        ! Classically (see the 'pml_boundary_acoustic' routine in 'pml_compute.f90'), condition at the PML outer boundary is Neumann for the PML variables which are potential variables because Dirichlet conditons are needed for displacement.
        ! Here, the PML variables are directly related to the actual variables, thus we apply the same outer condition.
        ! Set exact_interface_flux.
!        exact_interface_flux = .true.
!        out_drho_P = inp_drho_M ! Set out_drho_P: same as other side, that is a Neumann condition.
!        out_dp_P = inp_dp_M ! Set out_dp_P: same as other side, that is a Neumann condition.
!        out_dv_P = inp_rho0dv_M/LNS_rho0(iglobM) ! Set out_dv_P: same as other side, that is a Neumann condition.
!        call compute_dE_i(LNS_rho0(iglobM)+out_drho_P, LNS_v0(:,iglobM)+out_dv_P, LNS_p0(iglobM)+out_dp_P, out_dE_P, iglobM) ! Set out_dE_P: same as other side, that is a Neumann condition.
!        out_rho0dv_P = inp_rho0dv_M ! Set out_rho0dv_P: same as other side, that is a Neumann condition.
        ! Set out_dm_P: see bottom of routine.
!        if(swCompVisc) then
!          out_nabla_dT_P = nabla_dT(:,iglobM) ! Set out_nabla_dT_P: same as other side, that is a Neumann condition.
!          out_sigma_dv_P = sigma_dv(:,iglobM) ! Set out_sigma_dv_P: same as other side, that is a Neumann condition.
!        endif
!        if(swCompdT) then
!           call compute_dT_i(LNS_rho0(iglobM)+out_drho_P, LNS_v0(:,iglobM)+out_dv_P, LNS_E0(iglobM)+out_dE_P, out_dT_P, iglobM) ! Set out_dT_P: same as other side, that is a Neumann condition.
!         !call compute_dT_i(LNS_rho0(iglobM)+out_drho_P, LNS_p0(iglobM)+out_dp_P, out_dT_P, iglobM)
!        endif
!      else ! Else on PML.
        ! --------------------------- !
        ! Not outer PML.              !
        ! --------------------------- !
        ! Set exact_interface_flux.
        exact_interface_flux = .false.
        ! Get out_drho_P and out_dp_P (actually, this call gets rho=(rho0+drho) and p=(p0+dp), and this needs to be correctly right after).
        call background_physical_parameters(i, j, ispec, timelocal, out_drho_P, &
                                            .true., velocity_P(1), & ! We get the far-field v field here. "(1)" here is to please Intel compilers. It should not interfere too much since the concerned array is nicely allocated.
                                            .false., LNS_dummy_1d(1), &
                                            .true., out_dp_P) ! Get needed background parameters. Use dummies for values we're not interested in.
        out_drho_P = out_drho_P - LNS_rho0(iglobM) ! Warning: out_drho_P contains rho=(rho0 + drho) here. Correct this.
        out_dp_P = out_dp_P - LNS_p0(iglobM) ! Warning: out_dp_P contains p=(p0 + dp) here. Correct this.
        ! Set velocity_P. Convert the velocity components from mesh coordinates to normal/tangential coordinates. This is more practical to set the boundary conditions. The setting of the boundary conditions is thus made here.
        call build_trans_boundary(n_out, tang, trans_boundary)
        normal_v     = DOT_PRODUCT(n_out, velocity_P)
        tangential_v = DOT_PRODUCT(tang, velocity_P)
        ! Notes:
        !   At that point, normal_v and tangential_v are set to their "background" (or "far-field", or "unperturbed") values.
        !   Treatment of boundary conditions based on normal/tangential velocities should be done here. The "free slip" condition and the "normal velocity continuity" conditions can be set here.
        !   These outer boundary conditions agree with the current state of the classical DG (FNS) outer boundary conditions. It is recommended to check wether or not the two outer boundary conditions agree. To do so, compare these lines to the corresponding ones (search 'outer boundary conditions' or something like that) in 'boudnary_terms_DG.f90'. TODO: maybe a unified subroutine between FNS and LNS for far-field BC.
        ! Convert (back) the velocity components from normal/tangential coordinates to mesh coordinates.
        do SPCDM = 1, NDIM
          velocity_P(SPCDM) = trans_boundary(SPCDM, 1)*normal_v + trans_boundary(SPCDM, 2)*tangential_v
        enddo
        ! Set out_dv_P.
        out_dv_P = velocity_P - LNS_v0(:, iglobM) ! Requesting iglobM might be technically inexact, but on elements' boundaries points should overlap. Plus, iglobP does not exist on outer computational domain boundaries.
        ! Set out_dE_P.
        call compute_dE_i(LNS_rho0(iglobM)+out_drho_P, velocity_P, LNS_p0(iglobM)+out_dp_P, out_dE_P, iglobM)
        ! Set out_rho0dv_P.
        !do SPCDM = 1, NDIM
        !  out_rho0dv_P(SPCDM) = out_drho_P*out_dv_P(SPCDM) ! Safe version, just in case element-wise mutiplication fails somehow.
        !enddo
        out_rho0dv_P(:) = LNS_rho0(iglobM)*out_dv_P(:)
        ! Set out_dm_P: see bottom of routine.
        if(swCompVisc) then
          out_nabla_dT_P = nabla_dT(:,iglobM) ! Set out_nabla_dT_P: same as other side, that is a Neumann condition.
          out_sigma_dv_P = sigma_dv(:,iglobM) ! Set out_sigma_dv_P: same as other side, that is a Neumann condition.
        endif
        ! Set out_dT_P.
        if(swCompdT) then
          call compute_dT_i(LNS_rho0(iglobM)+out_drho_P, velocity_P, LNS_E0(iglobM)+out_dE_P, out_dT_P, iglobM)
          !call compute_dT_i(LNS_rho0(iglobM)+out_drho_P, LNS_p0(iglobM)+out_dp_P, out_dT_P, iglobM)
        endif
        
        
!        ! TEST MANUFACTURED SOLUTIONS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        exact_interface_flux = .false. ! Do not force jump to zero.
!        
!        out_drho_P = RHO_cst*(  sin(dRHO_x*PI*coord(1,ibool_before_perio(i,j,ispec))) &
!                              + sin(dRHO_z*PI*coord(2,ibool_before_perio(i,j,ispec))))
!        out_dv_P = 0.*out_dv_P
!        out_dE_P = E_cst*(  sin(dE_x*PI*coord(1,ibool_before_perio(i,j,ispec))) &
!                          + sin(dE_z*PI*coord(2,ibool_before_perio(i,j,ispec))))
!        call compute_dp_i(LNS_rho0(iglobM)+out_drho_P, LNS_v0(:,iglobM)+out_dv_P, LNS_E0(iglobM)+out_dE_P, out_dp_P, iglobM)
!        
!        out_rho0dv_P = LNS_rho0(iglobM)*out_dv_P
!        if(swCompVisc) then
!          out_nabla_dT_P = nabla_dT(:,iglobM) ! Set out_nabla_dT_P: same as other side, that is a Neumann condition.
!          out_sigma_dv_P = sigma_dv(:,iglobM) ! Set out_sigma_dv_P: same as other side, that is a Neumann condition.
!        endif
!        if(swCompdT) then
!          call compute_dT_i(LNS_rho0(iglobM)+out_drho_P, LNS_v0(:, iglobM)+out_dv_P, LNS_E0(iglobM)+out_dE_P, out_dT_P, iglobM)
!          !call compute_dT_i(LNS_rho0(iglobM)+out_drho_P, LNS_p0(iglobM)+out_dp_P, out_dT_P, iglobM)
!        endif
!        !! DO NOT FORGET TO UPDATE SOURCE TERM ABOVE, lines ~110.
!        ! TEST MANUFACTURED SOLUTIONS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
!      endif ! Endif on PML.
    endif ! Endif on ipoin.
    
  elseif(ispec_is_acoustic_coupling_ac(ibool_DG(i, j, ispec)) >= 0) then
    ! --------------------------- !
    ! neighbor(3) /= -1           !
    !   and acoustic coupling.    !
    ! --------------------------- !
    write(*,*) "********************************"
    write(*,*) "*            ERROR             *"
    write(*,*) "********************************"
    write(*,*) "* Interfaces between LNS and   *"
    write(*,*) "* acoustic potential elements  *"
    write(*,*) "* are not implemented yet.     *"
    write(*,*) "********************************"
    stop
    ! Set exact_interface_flux.
    exact_interface_flux = .true.
    ! Coordinates of elastic element
    i_ac     = neighbor(1)
    j_ac     = neighbor(2)
    ispec_ac = neighbor(3)
    iglob = ibool(i_ac,j_ac,ispec_ac)!ibool(i, j, ispec)
    !if(.false.) then ! DEBUG
    !  WRITE(*,*) ibool(i_ac,j_ac,ispec_ac), xixl, xizl, gammaxl, gammazl, duz_dxi, duz_dgamma, k
    !endif
    ! Set out_drho_P.
    call background_physical_parameters(i_ac, j_ac, ispec_ac, timelocal, out_drho_P, &
                                        .false., LNS_dummy_1d(1), &
                                        .false., LNS_dummy_1d(1), &
                                        .false., LNS_dummy_1d(1)) ! Get needed background parameters. Use dummies for values we're not interested in.
    !!out_drho_P = inp_drho_M
    !!duz_dxi    = ZEROcr
    !!duz_dgamma = ZEROcr
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
    veloc_P = veloc_vector_acoustic_DG_coupling(iglob, :)
    !veloc_x = veloc_vector_acoustic_DG_coupling(iglob, 1)
    !veloc_z = veloc_vector_acoustic_DG_coupling(iglob, 2)
    call build_trans_boundary(n_out, tang, trans_boundary)
    ! Tangential vector, since only bottom topography n_out(NDIM) > 0.
    !tang(1) = -n_out(NDIM) ! Recall: n_out(NDIM)=n_z.
    !tang(NDIM) =  n_out(1) ! Recall: n_out(1)=n_x.
    ! Normal velocity of the solid perturbation and tangential velocity of the fluid flow
    !normal_v     = veloc_P(1)*n_out(1) + veloc_P(1)*n_out(NDIM) ! Recall: n_out(NDIM)=n_z.
    !tangential_v = veloc_P(1)*tang(1) + veloc_P(1)*tang(NDIM)
    normal_v     = DOT_PRODUCT(n_out, veloc_P)
    tangential_v = DOT_PRODUCT(tang, veloc_P)
    ! Set the matrix for the transformation from normal/tangential coordinates to mesh coordinates.
    !trans_boundary(1, 1) =  tang(NDIM)
    !trans_boundary(1, 2) = -n_out(NDIM) ! Recall: n_out(NDIM)=n_z.
    !trans_boundary(2, 1) = -tang(1)
    !trans_boundary(2, 2) =  n_out(1) ! Recall: n_out(1)=n_x.
    !trans_boundary = trans_boundary/(n_out(1) * tang(NDIM) - tang(1) * n_out(NDIM))
    ! Set velocity_P: free slip and normal velocity continuity.
    do SPCDM = 1, NDIM
      velocity_P(SPCDM) = trans_boundary(SPCDM, 1)*normal_v + trans_boundary(SPCDM, 2)*tangential_v!veloc_elastic(1,iglob)
      !velocity_P(1) = trans_boundary(1, 1)*normal_v + trans_boundary(1, 2)*tangential_v!veloc_elastic(1,iglob)
      !velocity_P(NDIM) = trans_boundary(2, 1)*normal_v + trans_boundary(2, 2)*tangential_v
    enddo
    ! Set out_dv_P.
    out_dv_P = velocity_P - LNS_v0(:,iglobM) ! Requesting iglobM might be technically inexact, but on elements' boundaries points should overlap. Plus, iglobP does not exist on outer computational domain boundaries.
    ! Set out_dp_P: traction continuity.
    out_dp_P = LNS_p0(iglobM) - potential_dot_dot_acoustic(iglob) ! Warning, expression of out_dp_P might not be exact.
    ! Set out_dE_P.
    call compute_dE_i(LNS_rho0(iglobM)+out_drho_P, LNS_v0(:,iglobM)+out_dv_P, LNS_p0(iglobM)+out_dp_P, out_dE_P, iglobM) ! Warning, expression of out_dp_P might not be exact.
    !out_dE_P = out_dp_P/(gammaext_DG(iglobM) - 1.) + HALFcr*out_drho_P*( out_dv_P(1)**2 + out_dv_P(NDIM)**2 ) ! Warning, expression of out_dp_P might not be exact.
    ! Set out_rho0dv_P.
    do SPCDM = 1, NDIM
      out_rho0dv_P(SPCDM) = out_drho_P*out_dv_P(SPCDM) ! Safe version, just in case element-wise mutiplication fails somehow.
    enddo
    !out_rho0dv_P=out_drho_P*out_dv_P
    ! Set gradients, and sigma
    !out_nabla_dT_P
    !call compute_gradient_TFSF(LNS_dv, LNS_dT, .true., .true., switch_gradient, nabla_dv, nabla_dT)
    !out_dT_P = (inp_dE_M/inp_drho_M - 0.5*(dv_M(1)**2 + dv_M(NDIM)**2))/c_V
    ! Set out_dm_P: see bottom of routine.
    if(swCompVisc) then
      ! Set out_nabla_dT_P: same as other side, that is a Neumann condition.
      out_nabla_dT_P = nabla_dT(:,iglobM)
      ! Set out_sigma_dv_P: same as other side, that is a Neumann condition.
      out_sigma_dv_P = sigma_dv(:,iglobM)
    endif
    ! Set out_dT_P.
    if(swCompdT) then
      call compute_dT_i(LNS_rho0(iglobM)+out_drho_P, velocity_P, LNS_E0(iglobM)+out_dE_P, out_dT_P, iglobM)
      !call compute_dT_i(LNS_rho0(iglobM)+out_drho_P, LNS_p0(iglobM)+out_dp_P, out_dT_P, iglobM)
    endif
    
  else
    ! --------------------------- !
    ! neighbor(3) /= -1           !
    !   and not an 'outside' edge !
    !   and not acoustic          !
    !           potential         !
    ! --------------------------- !
    ! Here, values are known.
    
    iglobP = ibool_DG(neighbor(1), neighbor(2), neighbor(3))
    
    ! Set exact_interface_flux.
    exact_interface_flux = .false.
    !exact_interface_flux = .true. ! TEST
    
    ! Set out_drho_P.
    out_drho_P = inp_drho_P
    
    ! Set out_dE_P.
    out_dE_P = inp_dE_P
    
    ! Set out_rho0dv_P.
    out_rho0dv_P = inp_rho0dv_P
    
    ! Set out_dv_P.
    out_dv_P(:) = out_rho0dv_P(:)/LNS_rho0(iglobP)
    
    ! Set out_dp_P.
    ! iglobP exists and has meaning.
    call compute_dp_i(LNS_rho0(iglobP)+out_drho_P, LNS_v0(:,iglobP)+out_dv_P, LNS_E0(iglobP)+out_dE_P, out_dp_P, iglobP)
    
    ! Set out_dm_P: see bottom of routine.
    if(swCompVisc) then
    ! iglobP exists and has meaning.
      out_nabla_dT_P = nabla_dT(:,iglobP) ! Set out_nabla_dT_P: same as other side, that is a Neumann condition.
      out_sigma_dv_P = sigma_dv(:,iglobP) ! Set out_sigma_dv_P: same as other side, that is a Neumann condition.
    endif
    
    ! Set out_dT_P.
    if(swCompdT) then
      ! iglobP exists and has meaning.
      call compute_dT_i(LNS_rho0(iglobP)+out_drho_P, LNS_v0(:, iglobP)+out_dv_P, LNS_E0(iglobP)+out_dE_P, out_dT_P, iglobP)
      !call compute_dT_i(LNS_rho0(iglobP)+out_drho_P, LNS_p0(iglobP)+out_dp_P, out_dT_P, iglobP)
    endif
    !out_dT_P = (out_dE_P/out_drho_P - 0.5*((out_rho0dv_P(1)/out_drho_P)**2 + (out_rho0dv_P(NDIM)/out_drho_P)**2))/c_V
    
  endif
  
  ! Set out_dm_P.
  out_dm_P=out_rho0dv_P+out_drho_P*LNS_v0(:,iglobM)
  
end subroutine LNS_get_interfaces_unknowns



! ------------------------------------------------------------ !
! build_trans_boundary                                         !
! ------------------------------------------------------------ !
subroutine build_trans_boundary(normal, tangential, transf_matrix)
  use constants, only: CUSTOM_REAL, NDIM
  
  implicit none
  
  ! Input/Output.
  real(kind=CUSTOM_REAL), dimension(NDIM), intent(in) :: normal
  real(kind=CUSTOM_REAL), dimension(NDIM), intent(out) :: tangential
  real(kind=CUSTOM_REAL), dimension(NDIM, NDIM), intent(out) :: transf_matrix
  
  ! Local.
  ! N./A.
  
  ! For the tangential vector, it is assumed that we will only have the elastic media under the DG medium, hence we always have n_out(NDIM)>=0.
  tangential(1)    = -normal(NDIM) ! Recall: normal(NDIM)=n_z.
  tangential(NDIM) =  normal(1) ! Recall: normal(1)=n_x.
  transf_matrix(1, 1)       =   tangential(NDIM)
  transf_matrix(1, NDIM)    = - normal(NDIM) ! Recall: normal(NDIM)=n_z.
  transf_matrix(NDIM,    1) = - tangential(1)
  transf_matrix(NDIM, NDIM) =   normal(1) ! Recall: normal(1)=n_x.
  transf_matrix = transf_matrix/(normal(1)*tangential(NDIM) - tangential(1)*normal(NDIM))
end subroutine






! ------------------------------------------------------------ !
! S2F_Terrana_coupling                                         !
! ------------------------------------------------------------ !
! Implements [Terrana et al., 2018]'s (56).
! Terrana, S., Vilotte, J. P., & Guillot, L. (2017). A spectral hybridizable discontinuous Galerkin method for elastic–acoustic wave propagation. Geophysical Journal International, 213(1), 574-602.
subroutine S2F_Terrana_coupling(normal, rho_fluid, v_fluid, dp_fluid, soundspeed, iglob_DG, &
                                v_solid, sigma_el_local, i_el, j_el, ispec_el, &
                                v_hat, dp_hat)
  use constants, only: CUSTOM_REAL, NDIM
  use specfem_par_lns, only: inverse_2x2, LNS_rho0, LNS_v0, LNS_dv!, LNS_viscous, sigma_dv
  
  implicit none
  
  ! Input/Output.
  real(kind=CUSTOM_REAL), dimension(NDIM), intent(in) :: normal, v_fluid, v_solid
  real(kind=CUSTOM_REAL), intent(in) :: rho_fluid, dp_fluid, soundspeed
  integer, intent(in) :: iglob_DG, i_el, j_el, ispec_el
  real(kind=CUSTOM_REAL), dimension(NDIM, NDIM), intent(in) :: sigma_el_local
  real(kind=CUSTOM_REAL), dimension(NDIM), intent(out) :: v_hat
  real(kind=CUSTOM_REAL), intent(out) :: dp_hat
  
  ! Local.
  real(kind=CUSTOM_REAL), dimension(NDIM, NDIM) :: TAU_F, TAU_S, sigma_fluid
  
  ! Build tensors \tau.
  call build_tau_s(normal, i_el, j_el, ispec_el, TAU_S) ! TAU_S could be pre-computed and stored, since it doesn't change with time. ! TODO, maybe: an optimisation.
  call build_tau_f(normal, rho_fluid, soundspeed, TAU_F) ! TAU_F couldn't be pre-computed, since some values change with time.
  ! Build inverse.
!  INV_SUMTAU = inverse_2x2(TAU_S + TAU_F)
  
  ! Prepare fluid stress tensor.
  sigma_fluid = 0. ! Initialise stress.
  call stressBuilder_addInviscidFluid(LNS_rho0(iglob_DG), LNS_v0(:,iglob_DG), &
                                      LNS_dv(:,iglob_DG), dp_fluid, sigma_fluid)
  !if(LNS_viscous) then ! Check if viscosity exists whatsoever.
  !  call stressBuilder_addViscousFluid(-sigma_dv(:,iglob_DG), sigma_fluid) ! Send viscous tensor.
  !endif
  ! Build actual velocity from [Terrana et al., 2018]'s (56).
  !v_hat = matmul(TAU_F,v_fluid) + matmul(TAU_S,v_solid) + matmul(sigma_el_local,normal) - dp_fluid*normal
  !v_hat = matmul(TAU_F, v_fluid) + matmul(TAU_S, v_solid) + matmul(sigma_el_local - sigma_fluid, normal)
  !v_hat = matmul(inverse_2x2(TAU_S+TAU_F), v_hat)
  ! Do everything inline. Four 2x2 little matrix-vector products won't crash the stack.
  v_hat = matmul(inverse_2x2(TAU_S+TAU_F), &
                 matmul(TAU_F, v_fluid) + matmul(TAU_S, v_solid) + matmul(sigma_el_local - sigma_fluid, normal))
  
  ! Build actual pressure perturbation from [Terrana et al., 2018]'s (51).
  dp_hat = dp_fluid + DOT_PRODUCT(matmul(TAU_F, v_fluid - v_solid), normal)
end subroutine S2F_Terrana_coupling

! ------------------------------------------------------------ !
! build_tau                                                    !
! ------------------------------------------------------------ !
! Generic tensor \tau, from [Terrana et al., 2018]'s (37).
! Tau is a symmetric matrix. It is a bit heavy to store all components, maybe consider storing only the upper coefficients.
! Terrana, S., Vilotte, J. P., & Guillot, L. (2017). A spectral hybridizable discontinuous Galerkin method for elastic–acoustic wave propagation. Geophysical Journal International, 213(1), 574-602.
subroutine build_tau_s(normal, i_el, j_el, ispec_el, TAU_S)
  use constants, only: CUSTOM_REAL, NDIM
  use specfem_par, only: assign_external_model, density, kmato, poroelastcoef, &
                         rhoext, vpext, vsext
  
  implicit none
  
  ! Input/Output.
  real(kind=CUSTOM_REAL), dimension(NDIM), intent(in) :: normal
  integer, intent(in) :: i_el, j_el, ispec_el
  real(kind=CUSTOM_REAL), dimension(NDIM, NDIM), intent(out) :: TAU_S
  
  ! Local.
  real(kind=CUSTOM_REAL) :: mu_elastic_unrelaxed, rho, vp, vs
  
  ! Get elastic parameters.
  mu_elastic_unrelaxed = poroelastcoef(2, 1, kmato(ispec_el))
  rho = density(1, kmato(ispec_el))
  vp = sqrt((poroelastcoef(1, 1, kmato(ispec_el)) + 2.*mu_elastic_unrelaxed)/rho) ! vp = ((lambda+2mu)/rho)^0.5 http://www.subsurfwiki.org/wiki/P-wave_modulus, poroelastcoef(1,1,kmato(ispec_el)) = lambda
  vs = sqrt(mu_elastic_unrelaxed/rho) ! vs = (mu/rho)^0.5 http://www.subsurfwiki.org/wiki/P-wave_modulus
  if(assign_external_model) then
    rho = rhoext(i_el, j_el, ispec_el)
    vp = vpext(i_el, j_el, ispec_el)
    vs = vsext(i_el, j_el, ispec_el)
  endif
  
  TAU_S(1,1) = vp*normal(1)**2 + vs*normal(NDIM)**2
  TAU_S(1,2) = (vp - vs)*normal(1)*normal(NDIM)
  TAU_S(2,1) = TAU_S(1,2)
  TAU_S(2,2) = vs*normal(1)**2 + vp*normal(NDIM)**2
  TAU_S = rho*TAU_S
end subroutine build_tau_s
! ------------------------------------------------------------ !
! build_tau_f                                                  !
! ------------------------------------------------------------ !
! Acoustic tensor \tau_{ac}, from [Terrana et al., 2018]'s (52).
! Tau is a symmetric matrix. It is a bit heavy to store all components, maybe consider storing only the upper coefficients.
! Terrana, S., Vilotte, J. P., & Guillot, L. (2017). A spectral hybridizable discontinuous Galerkin method for elastic–acoustic wave propagation. Geophysical Journal International, 213(1), 574-602.
subroutine build_tau_f(normal, rho, c, TAU_F)
  use constants, only: CUSTOM_REAL, NDIM
  
  implicit none
  
  ! Input/Output.
  real(kind=CUSTOM_REAL), dimension(NDIM), intent(in) :: normal
  real(kind=CUSTOM_REAL), intent(in) :: rho, c
  real(kind=CUSTOM_REAL), dimension(NDIM, NDIM), intent(out) :: TAU_F
  
  ! Local.
  ! N./A.
  
  ! Expansive version.
  TAU_F(1,1) = normal(1)**2
  TAU_F(1,2) = normal(1)*normal(NDIM)
  TAU_F(2,1) = TAU_F(1,2)
  TAU_F(2,2) = normal(NDIM)**2
  TAU_F = rho*c*TAU_F
end subroutine build_tau_f





! ------------------------------------------------------------ !
! LNS_mass_source                                              !
! ------------------------------------------------------------ !
! Implements a source on the mass conservation equation:
! 1) Add said source on the mass equation.
! 2) Ensure compatibility by adding it also to the other equations.
! Based off the compute_add_sources_acoustic_DG_mass routine in compute_add_sources_acoustic_DG.f90. If you make changes here, make sure to have them also there.

subroutine LNS_mass_source(d_drho, d_rho0dv, d_dE, it, i_stage)

  use constants,only: CUSTOM_REAL, NGLLX, NGLLZ, PI, HUGEVAL,NDIM
  use specfem_par, only: nglob_DG,ispec_is_acoustic_DG,&
                         NSOURCES,myrank,&
                         source_time_function,&
                         ibool_DG,gammaext_DG,&
                         jacobian, wxgll, wzgll, ibool_before_perio, &
                         USE_SPREAD_SSF, nspec, source_spatial_function_DG
  use specfem_par_LNS ! TODO: select variables to use.
  implicit none

  ! Input/output.
  real(kind=CUSTOM_REAL), dimension(nglob_DG), intent(inout) :: d_drho
  real(kind=CUSTOM_REAL), dimension(NDIM, nglob_DG), intent(inout) :: d_rho0dv
  real(kind=CUSTOM_REAL), dimension(nglob_DG), intent(inout) :: d_dE
  integer, intent(in) :: it, i_stage
  
  ! Local variables.
  real(kind=CUSTOM_REAL) :: temp_sourcewxlwzljacobianl
  integer :: i, j, i_source, ispec, iglob, iglob_unique, SPCDM
  real(kind=CUSTOM_REAL) :: stf ! In order to store the source time function at current timestep outside the many loops.
  do i_source = 1, NSOURCES ! Loop on sources.
    stf = source_time_function(i_source, it, i_stage) ! Store the source time function outside the many loops.
    if(USE_SPREAD_SSF) then
      do ispec = 1, nspec
        if(ispec_is_acoustic_DG(ispec)) then
          do j = 1, NGLLZ
            do i = 1, NGLLX
              iglob_unique = ibool_before_perio(i, j, ispec)
              iglob = ibool_DG(i, j, ispec)
              ! See "prepare_source_spatial_function.f90" for the subroutine initialising the vector "source_spatial_function_DG".
              temp_sourcewxlwzljacobianl =   stf &
                                           * source_spatial_function_DG(i_source, iglob_unique) &
                                           * real(wxgll(i), kind=CUSTOM_REAL) &
                                           * real(wzgll(j), kind=CUSTOM_REAL) &
                                           * jacobian(i, j, ispec)
              
              d_drho(iglob) = d_drho(iglob) + temp_sourcewxlwzljacobianl
              
              do SPCDM = 1, NDIM
                d_rho0dv(SPCDM, iglob) = d_rho0dv(SPCDM, iglob) + LNS_dv(SPCDM, iglob) * temp_sourcewxlwzljacobianl
              enddo
              
              ! Approx. c^2=c_0^2 & v=v_0. Lowest CPU use, 18% error at Mach .3, 7% error at Mach .3.
              d_dE(iglob) =   d_dE(iglob) &
                           + (   LNS_c0(iglob)**2/(gammaext_DG(iglob)-1.) &
                               + 0.5*norm2r1(LNS_v0(:, iglob)) &
                             ) * temp_sourcewxlwzljacobianl
              ! Approx. c^2=c_0^2. TBTested.
              !d_dE(iglob) =   d_dE(iglob) &
              !             + (   LNS_c0(iglob)**2/(gammaext_DG(iglob)-1.) &
              !                 + 0.5*norm2r1(LNS_v0(:, iglob)+LNS_dv(:, iglob)) &
              !               ) * temp_sourcewxlwzljacobianl
              ! Full. Highest CPU use, but technically exactly what the equations are, same results as "Approx. c^2=c_0^2 & v=v_0".
              !d_dE(iglob) =   d_dE(iglob) &
              !             + (   gammaext_DG(iglob) * (LNS_p0(iglob)+LNS_dp(iglob)) &
              !                   / ((LNS_rho0(iglob)+LNS_drho(iglob))*(gammaext_DG(iglob)-1.)) &
              !                 + 0.5*norm2r1(LNS_v0(:, iglob)+LNS_dv(:, iglob)) &
              !               ) * temp_sourcewxlwzljacobianl
            enddo
          enddo
        endif
      enddo
    else
      if(myrank==0) then
        write(*,*) "********************************"
        write(*,*) "*            ERROR             *"
        write(*,*) "********************************"
        write(*,*) "* Mass source is not yet       *"
        write(*,*) "* implemented with             *"
        write(*,*) "* USE_SPREAD_SSF=.false.. See  *"
        write(*,*) "* compute_forces_acoustic_LNS.f90."
        write(*,*) "********************************"
        stop
      endif
    endif ! Endif on SIGMA_SSF.
  enddo ! Enddo on i_source.
end subroutine LNS_mass_source














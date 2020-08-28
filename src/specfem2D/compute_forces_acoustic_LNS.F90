! ------------------------------------------------------------ !
! compute_forces_acoustic_LNS                                  !
! ------------------------------------------------------------ !
! Computes the right-hand side of the differential LNS (Linear Navier-Stokes) system

subroutine compute_forces_acoustic_LNS(cv_drho, cv_rho0dv, cv_dE, cv_e1, & ! Constitutive variables.
                                       in_dm, in_dp, in_nabla_dT, in_nabla_dv, in_sigma_dv, & ! Precomputed quantities.
                                       outrhs_drho, outrhs_rho0dv, outrhs_dE, outrhs_e1, & ! Output (RHS for each constitutive variable).
                                       currentTime) ! Time.
  use constants, only: CUSTOM_REAL, NGLLX, NGLLZ, NDIM, TINYVAL, PI
  use specfem_par, only: myrank, &
                         jacobian, xix, xiz, gammax, gammaz, wxgll, wzgll, hprimewgll_xx, hprimewgll_zz, &
                         link_iface_ijispec, nx_iface, nz_iface, weight_iface, neighbor_dg_iface, &
                         nspec, nglob_DG, ibool_DG, ibool_before_perio, ispec_is_acoustic_dg, it, i_stage, gammaext_dg, &
                         TYPE_SOURCE_DG, ABC_STRETCH, stretching_ya, stretching_buffer
  use specfem_par_LNS, only: USE_LNS, NVALSIGMA, LNS_viscous, &
                             LNS_dummy_1d, &
                             LNS_rho0, &
                             LNS_v0, LNS_p0, LNS_E0, LNS_dv, &
                             nabla_v0, sigma_v_0, &
                             LNS_kappa, LNS_mu, LNS_eta, LNS_g, &
                             VALIDATION_MMS, &
                             LNS_avib, LNS_avib_taueps, LNS_avib_tausig, &
                             LNS_verbose, LNS_modprint
  
  implicit none
  
  ! Input/Output.
  real(kind=CUSTOM_REAL), dimension(nglob_DG), intent(in) :: cv_drho, cv_dE, cv_e1, in_dp!, in_dT
  real(kind=CUSTOM_REAL), dimension(NDIM, nglob_DG), intent(in) :: cv_rho0dv, in_dm, in_nabla_dT
  real(kind=CUSTOM_REAL), dimension(NDIM, NDIM, nglob_DG), intent(in) :: in_nabla_dv
  real(kind=CUSTOM_REAL), dimension(NVALSIGMA, nglob_DG), intent(in) :: in_sigma_dv
  real(kind=CUSTOM_REAL), dimension(nglob_DG), intent(out) :: outrhs_drho, outrhs_dE, outrhs_e1
  real(kind=CUSTOM_REAL), dimension(NDIM, nglob_DG), intent(out) :: outrhs_rho0dv
  real(kind=CUSTOM_REAL), intent(in) :: currentTime
  
  ! Local.
  real(kind=CUSTOM_REAL), parameter :: ZEROcr = 0._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: ONEcr  = 1._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: TWOcr  = 2._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: HALFcr = 0.5_CUSTOM_REAL
  integer :: ispec, i, j, k, iglob, SPCDM, iglob_unique
  
  ! > Variables "cntrb_*" are aimed at assembling the different contributions to constitutive variables.
  !   For drho and dE, "cntrb_*(:,:,1)" contains the contribution along xi, and "cntrb_*(:,:,2)" the contribution along gamma.
  !   For rho0dv, "cntrb_*(1,:,:,1)" contains the contribution to rho0dvx along xi, "cntrb_*(1,:,:,2)" the contribution to rho0dvx along gamma, "cntrb_*(2,:,:,1)" the contribution to rhodvz along xi, and "cntrb_*(2,:,:,2)" the contribution to rhodvz along gamma.
  ! > Variables "d0cntrb_*" are aimed at assembling contributions of zero-th degree (outside divergence operator).
  real(kind=CUSTOM_REAL), dimension(NGLLX, NGLLZ, NDIM) :: cntrb_drho, cntrb_dE
  real(kind=CUSTOM_REAL), dimension(NDIM, NGLLX, NGLLZ, NDIM) :: cntrb_rho0dv
  real(kind=CUSTOM_REAL), dimension(NGLLX, NGLLZ) :: d0cntrb_dE
  real(kind=CUSTOM_REAL), dimension(NDIM, NGLLX,NGLLZ) :: d0cntrb_rho0dv
  
  real(kind=CUSTOM_REAL) :: Jac_L ! Jacobian matrix determinant.
  real(kind=CUSTOM_REAL) :: lambda, halfWeight, jump
  real(kind=CUSTOM_REAL), dimension(NDIM) :: Sigma_L ! Local stress tensor.
  real(kind=CUSTOM_REAL), dimension(NDIM,NDIM) :: DXiEta_L
  real(kind=CUSTOM_REAL), dimension(NDIM) :: Jac_WzWx_L
  real(kind=CUSTOM_REAL), dimension(NDIM) :: n_out
  real(kind=CUSTOM_REAL), dimension(NDIM) :: rho0dv_P
  integer :: iglobP
  integer, dimension(3) :: neighbor
  integer :: neighbour_type
  
  ! Variables specifically for LNS_get_interfaces_unknowns.
  real(kind=CUSTOM_REAL), dimension(NDIM) :: v0_P, dv_P, dm_P, nabla_dT_P, rho_MP
  real(kind=CUSTOM_REAL) :: rho0_P, drho_P, p0_P, dp_P, E0_P, dE_P
  real(kind=CUSTOM_REAL) :: gam_P, kap_P
  real(kind=CUSTOM_REAL), dimension(NVALSIGMA) :: sigma_dv_P, sigma_v0_p
  real(kind=CUSTOM_REAL) :: wxlwzl
  logical :: viscousComputation, exact_interface_flux
  integer :: iface1, iface, iface1_neighbor, iface_neighbor, ispec_neighbor
  
  if(.not. USE_LNS) then
    stop "THIS ROUTINE SHOULD NOT BE CALLED IF USE_LNS=.false."
  endif
  
  ! Initialisation of the RHS.
  outrhs_drho    = ZEROcr
  outrhs_rho0dv  = ZEROcr
  outrhs_dE      = ZEROcr
  if(LNS_avib) outrhs_e1 = ZEROcr
  
  ! Start by adding source terms.
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
  
  IF(VALIDATION_MMS) call VALIDATION_MMS_source_terms(outrhs_drho, outrhs_rho0dv, outrhs_dE)
  
  if(myrank == 0 .and. LNS_VERBOSE>=51 .and. mod(it, LNS_MODPRINT) == 0) then
    write(*,"(a,i6,a)") "Informations for process number ", myrank, "."
    write(*,"(a)")                   " quantity [                 max    ,                  min    ]"
    write(*,"(a,e24.16,a,e24.16,a)") " drho     [", maxval(cv_drho), ", ", minval(cv_drho), "]"
    do SPCDM = 1, NDIM
      write(*,"(a,e24.16,a,e24.16,a,i1)") " rho0dv_i [", maxval(cv_rho0dv(SPCDM,:)), ", ", minval(cv_rho0dv(SPCDM,:)), "], i=",SPCDM
    enddo
    write(*,"(a,e24.16,a,e24.16,a)") " dE       [", maxval(cv_dE), ", ", minval(cv_dE), "]"
  endif
  
  do ispec = 1, nspec ! Loop over elements.
    if (ispec_is_acoustic_DG(ispec)) then ! Only do something for DG elements.
      
      ! --------------------------- !
      ! First set of loops: compute !
      ! volumic contributions.      !
      ! --------------------------- !
      ! First, zero the temporary contribution fields. This might not be needed, since they are strictly overwritten below, but this is safer.
      cntrb_drho     = ZEROcr
      cntrb_rho0dv   = ZEROcr
      cntrb_dE       = ZEROcr
      d0cntrb_rho0dv = ZEROcr
      d0cntrb_dE     = ZEROcr
      
      ! Then, loop on the GLL points to set each of the (i, j) component to the temporary contribution fields.
      do j = 1, NGLLZ
        do i = 1, NGLLX
          iglob = ibool_DG(i,j,ispec)
          
          Jac_L                = jacobian(i,j,ispec)
          DXiEta_L(1,    1)    = xix(i,j,ispec) ! = \partial_x\xi
          DXiEta_L(1,    NDIM) = xiz(i,j,ispec) ! = \partial_z\xi
          DXiEta_L(NDIM, 1)    = gammax(i,j,ispec) ! = \partial_x\eta
          DXiEta_L(NDIM, NDIM) = gammaz(i,j,ispec) ! = \partial_z\eta
          
          if(ABC_STRETCH .and. stretching_buffer(ibool_before_perio(i,j,ispec))>0) then
            ! For real stretching, \Sigma for each constitutive variable becomes \Ya\Sigma. It is heavy to change each and every expression where \Sigma arises. Rather, we make use of what is multiplying \Sigma.
            ! Here, in the volume integrations, easiest is the local derivatives involved in the change of reference.
            ! \partial_x becomes \ya_x\partial_x, and since \partial_x=(\partial_x\xi)\partial_\xi+(\partial_x\eta)\partial_\eta, only updating \partial_x\xi and \partial_x\eta is enough. Idem for \partial_z. Hence, only updating xix to \ya_x * xix, xiz to \ya_z * xiz, etc. is enough to update the operator.
            iglob_unique      = ibool_before_perio(i,j,ispec)
            DXiEta_L(:, 1)    = stretching_ya(1, iglob_unique)*DXiEta_L(:, 1)    ! Multiply x component by ya_x. ! If you change anything, remember to do it also in the 'compute_gradient_TFSF' subroutine.
            DXiEta_L(:, NDIM) = stretching_ya(2, iglob_unique)*DXiEta_L(:, NDIM) ! Multiply x component by ya_z. ! If you change anything, remember to do it also in the 'compute_gradient_TFSF' subroutine.
          endif
          
          Jac_WzWx_L(1)    = wzgll(j)*Jac_L
          Jac_WzWx_L(NDIM) = wxgll(i)*Jac_L
          ! In Jac_WzWx_L, notice how 1 is z and 2 is x. This comes from the use of the chain rule on the bilinear form, and reorganisation. See Martire's PhD thesis manuscript.
          
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
          !   Jac_WzWx_L(SPCDM)*DOT_PRODUCT(DXiEta_L(SPCDM,:), Sigma_L)
          ! in order to store successively contributions to xi (\partial_x\xi * \Sigma_x + \partial_z\xi * \Sigma_z) and to gamma (\partial_x\gamma * \Sigma_x + \partial_z\gamma * \Sigma_z).
          
          ! 1) Inviscid stress tensor's contributions.
          ! 1.1-2) Group mass conservation and x-Momentum in loop for efficiency.
          Sigma_L(1)    = cv_rho0dv(1,iglob)*LNS_v0(1,iglob) + in_dp(iglob) ! x-Momentum \Sigma_1.
          Sigma_L(NDIM) = cv_rho0dv(NDIM,iglob)*LNS_v0(1,iglob) ! x-Momentum \Sigma_2.
          do SPCDM = 1, NDIM
            cntrb_drho(i,j,SPCDM)     = Jac_WzWx_L(SPCDM)*DOT_PRODUCT(DXiEta_L(SPCDM,:), in_dm(:,iglob))
            cntrb_rho0dv(1,i,j,SPCDM) = Jac_WzWx_L(SPCDM)*DOT_PRODUCT(DXiEta_L(SPCDM,:), Sigma_L)
          enddo
          ! 1.2) z-Momentum.
          Sigma_L(1)    = cv_rho0dv(1,iglob)*LNS_v0(NDIM,iglob)
          Sigma_L(NDIM) = cv_rho0dv(NDIM,iglob)*LNS_v0(NDIM,iglob) + in_dp(iglob)
          do SPCDM = 1, NDIM
            cntrb_rho0dv(NDIM,i,j,SPCDM) =   Jac_WzWx_L(SPCDM)&
                                           * DOT_PRODUCT(DXiEta_L(SPCDM,:), Sigma_L)
          enddo
          
          ! 2) Add viscous stress tensor's contributions.
          if(viscousComputation) then
            ! 2.1) Mass conservation: no viscous contribution.
            ! 2.2) x-Momentum.
            Sigma_L(1)    = -in_sigma_dv(1,iglob) ! Recall, this corresponds to \Sigma_v(v')_{1,1}.
            Sigma_L(NDIM) = -in_sigma_dv(2,iglob) ! Recall, this corresponds to \Sigma_v(v')_{1,2} (and \Sigma_v(v')_{2,1}).
            do SPCDM = 1, NDIM
              cntrb_rho0dv(1,i,j,SPCDM) =   cntrb_rho0dv(1,i,j,SPCDM) &
                                          + Jac_WzWx_L(SPCDM)*DOT_PRODUCT(DXiEta_L(SPCDM,:), Sigma_L)
            enddo
            ! 2.2) z-Momentum.
            Sigma_L(1)    = -in_sigma_dv(2,iglob) ! Recall, this corresponds to \Sigma_v(v')_{2,1} (and \Sigma_v(v')_{1,2}).
            Sigma_L(NDIM) = -in_sigma_dv(3,iglob) ! Recall, this corresponds to \Sigma_v(v')_{2,2}.
            do SPCDM = 1, NDIM
              cntrb_rho0dv(NDIM,i,j,SPCDM) =   cntrb_rho0dv(NDIM,i,j,SPCDM) &
                                             + Jac_WzWx_L(SPCDM)*DOT_PRODUCT(DXiEta_L(SPCDM,:), Sigma_L)
            enddo
          endif ! Endif on viscousComputation.
          
          ! 3) Special case: energy.
          ! Indeed, if they exist, viscous contributions can be grouped to inviscid ones.
          if(viscousComputation) then
            Sigma_L(1) =   LNS_dv(1,iglob)*(LNS_E0(iglob) + LNS_p0(iglob) - sigma_v_0(1,iglob)) &
                          + LNS_v0(1,iglob)*(cv_dE(iglob) + in_dp(iglob) - in_sigma_dv(1,iglob)) &
                          - LNS_dv(NDIM,iglob)*sigma_v_0(2,iglob) &
                          - LNS_v0(NDIM,iglob)*in_sigma_dv(2,iglob) &
                          - LNS_kappa(iglob)*in_nabla_dT(1,iglob)
            Sigma_L(NDIM) =   LNS_dv(NDIM,iglob)*(LNS_E0(iglob) + LNS_p0(iglob) - sigma_v_0(3,iglob)) &
                             + LNS_v0(NDIM,iglob)*(cv_dE(iglob) + in_dp(iglob) - in_sigma_dv(3,iglob)) &
                             - LNS_dv(1,iglob)*sigma_v_0(2,iglob) &
                             - LNS_v0(1,iglob)*in_sigma_dv(2,iglob) &
                             - LNS_kappa(iglob)*in_nabla_dT(NDIM,iglob)
          else ! Else on viscousComputation.
            do SPCDM = 1, NDIM
              Sigma_L(SPCDM) =   LNS_dv(SPCDM,iglob)*(LNS_E0(iglob) + LNS_p0(iglob)) &
                               + LNS_v0(SPCDM,iglob)*(cv_dE(iglob)  + in_dp(iglob))
            enddo
          endif ! Endif on viscousComputation.
          do SPCDM = 1, NDIM
            cntrb_dE(i,j,SPCDM) = Jac_WzWx_L(SPCDM)*DOT_PRODUCT(DXiEta_L(SPCDM,:), Sigma_L)
          enddo
          
          ! 4) Zero-th degree contributions.
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
          d0cntrb_rho0dv(1,i,j)    = DOT_PRODUCT(in_dm(:,iglob), nabla_v0(1,:,iglob)) ! \rho'g_1=0; + {\delta_m}\cdot[\nabla v_{0}]_{1,:}
          d0cntrb_rho0dv(NDIM,i,j) = - cv_drho(iglob)*LNS_g(iglob) & ! \rho'g_i!=0;
                                     + DOT_PRODUCT(in_dm(:,iglob), nabla_v0(NDIM,:,iglob)) ! {\delta_m}\cdot[\nabla v_{0}]_{d,:}
          d0cntrb_rho0dv(:,i,j) = d0cntrb_rho0dv(:,i,j) * Jac_L ! Factor Jacobian (small optimisation).
          ! 4.3) Energy.
          ! 4.2.1) Version 1: most general.
          !d0cntrb_dE(i,j) = - DOT_PRODUCT(locG, in_dm(:,iglob)) * Jac_L
          ! 4.2.2) Version 2: [gravity is vertical (locG(1)=0, locG(2)=-LNS_g)].
          d0cntrb_dE(i,j) = - LNS_g(iglob)*in_dm(NDIM,iglob) * Jac_L
          
          if(LNS_avib) then
            ! Add vibrational attenuation memory variable contribution.
            ! 10.1007/s11214-016-0324-6, equation (19), describes where in the pressure evolution equation the memory variable would be used.
            d0cntrb_dE(i,j) =   d0cntrb_dE(i,j) &
                              - Jac_L &
                                * ((LNS_p0(iglob)*gammaext_DG(iglob))/(gammaext_DG(iglob)-ONEcr)) &
                                * (   (LNS_avib_taueps(iglob)/LNS_avib_tausig(iglob)-ONEcr) &
                                      * (in_nabla_dv(1,1,iglob) + in_nabla_dv(NDIM,NDIM,iglob) - cv_e1(iglob)) &
                                    + cv_e1(iglob))
          endif
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
          
          ! Vibrational attenuation.
          if(LNS_avib) then
            outrhs_e1(iglob) =   outrhs_e1(iglob) &
                               - (ONEcr/LNS_avib_tausig(iglob)) &
                                 * (   (ONEcr - (LNS_avib_tausig(iglob)/LNS_avib_taueps(iglob)))              &
                                     * (in_nabla_dv(1,1,iglob)+in_nabla_dv(NDIM,NDIM,iglob)) + cv_e1(iglob) ) &
                                 * wxlwzl*jacobian(i,j,ispec)
          endif
          
          ! Add zero-th order terms.
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
          iglob = ibool_DG(i,j,ispec)
          
          if(LNS_viscous) then ! Check if viscosity exists whatsoever.
            ! Activate/deactivate, for this particular point (iglob), computation of quantities only needed when viscosity is present.
            if(     LNS_mu(iglob) > TINYVAL &
               .OR. LNS_eta(iglob) > TINYVAL &
               .OR. LNS_kappa(iglob) > TINYVAL) then
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
            ispec_neighbor  = neighbor_DG_iface(iface1, iface, ispec, 3)
            neighbor(1:2)   = link_iface_ijispec(iface1_neighbor, iface_neighbor, ispec_neighbor,1:2)
            neighbor(3)     = ispec_neighbor
            neighbour_type = 1 ! <=> a neighbouring LNS DG element was found in the same partition
          else
            neighbour_type = 10
          endif
          ! At this point, either:
          !   neighbor(1:3)=[-1, -1, -1]
          !   neighbour_type=10
          !     (because neighbor_DG_iface(iface1, iface, ispec, 3)>-1)
          !     which means the neighbouring "element" is either
          !       *) of type acoustic DG but in another partition (neighbour_type=11)
          !       *) of another material type which for now means elastic only (neighbour_type=21)
          !       *) the outside boundary (neighbour_type=99)
          ! or:
          !   neighbor(1:3)=[link_iface_ijispec(...), ispec_neighbor]
          !   neighbour_type = 1
          !     (because because neighbor_DG_iface(iface1, iface, ispec, 3)<=-1)
          !     which means 1) a neighbouring element was found in the same partition, and
          !                 2) the (i,j,ispec) coordinates are exactly (neighbor(1), neighbor(2), neighbor(3))
          
          ! Step 2: knowing the normals' parameters, compute now the fluxes.
          !iglobP = -1 ! Set iglobP like this in order to crash the program when this variable is badly used.
          ! Note: iglobP should NOT be used except when we are sure it makes sense, i.e. when iglobP does indeed point toward an element within th same CPU partition.
          !       Otherwise, i.e. if neighbor(1)<=-1, iglobP remains 1, and all quantities queried are from the 1-th index in this partition, which makes NO sense.
          !if(neighbor(1) > -1) then
          !  iglobP = ibool_DG(neighbor(1), neighbor(2), neighbor(3))
          !endif
          if(neighbour_type==1) then
            ! A neighbouring LNS DG element was found in the same partition.
            iglobP = ibool_DG(neighbor(1), neighbor(2), neighbor(3))
          else
            iglobP = -1000 ! Set iglobP like this in order to deliberately crash the program whenever this variable is badly used.
          endif
          
          exact_interface_flux = .false. ! Reset this variable to .false.: by default, the fluxes have to be computed (jump!=0). In some specific cases (assigned during the call to LNS_get_interfaces_unknowns), the flux can be exact (jump==0).
          call LNS_get_interfaces_unknowns(i, j, ispec, iface1, iface, neighbor, neighbour_type, currentTime, & ! Point identifier (input).
                  cv_drho(iglob), cv_rho0dv(:,iglob), & ! Input constitutive variables, "M" side.
                  cv_drho(iglobP), cv_rho0dv(:,iglobP), cv_dE(iglobP), & ! Input constitutive variables, "P" side. Note they make no sense if neighbor(1)<=-1.
                  in_dp(iglob), & ! Input other variable, "M" side.
                  n_out, & ! Normal vector (input).
                  exact_interface_flux, & ! Switch to disable jump in some cases (output).
                  drho_P, rho0dv_P, dE_P, & ! Output constitutive variables.
                  dm_P, dp_P, dv_P,& ! Output other variables.
                  viscousComputation, nabla_dT_P, sigma_dv_P, & ! Output other variables: viscous.
                  .false., LNS_dummy_1d(1)) ! Use dummies and set the switch to false not to compute unecessary quantities.
          call check_neighbour_type(neighbour_type) ! Safeguard: crash the program if neighbour_type outside of possible values.
          
          call LNS_get_bgqts_at_interface(i, j, ispec, iface1, iface, neighbor, neighbour_type, iglobP, &
                                          rho0_P, v0_P, E0_P, p0_P, sigma_v0_P, gam_P, kap_P)
          
          jump   = ZEROcr ! Safeguard.
          
          ! Recover an approximate local maximum linearized acoustic wave speed. See for example Hesthaven (doi.org/10.1007/9780387720678), page 208.
          lambda = ZEROcr
#if 0
          ! Approximate lambda (i.e., using the initial speed of sound c0).
          lambda = max(  abs(dot_product(n_out, LNS_v0(:,iglob))) & ! v_-\cdot n
                       + LNS_c0(iglob) & ! Local sound speed, side "M".
                       , abs(dot_product(n_out, LNS_v0(:,iglobP))) & ! v_+\cdot n
                       + LNS_c0(iglobP) & ! Local sound speed, side "P".
                       )
#endif          
#if 1
          ! Exact lambda (i.e., with the exact full v=v0+v' and c=c0+c').
          rho_MP(1) = LNS_rho0(iglob)+cv_drho(iglob)
          rho_MP(2) = rho0_P+drho_P
          if(rho_MP(1)>TINYVAL .and. rho_MP(2)>TINYVAL) then
            lambda = max(  abs(dot_product(n_out, LNS_v0(:,iglob)+LNS_dv(:,iglob))) &
                         + sqrt(gammaext_DG(iglob)*(LNS_p0(iglob)+in_dp(iglob))/rho_MP(1)) &
                         , abs(dot_product(n_out, v0_P(:)+dv_P)) &
                         + sqrt(gam_P*(p0_P+dp_P)/rho_MP(2)) &
                         )
          else
            ! Possibly, in elastic elements, rho_MP is zero (since both LNS_rho0 and cv_drho and drho_P are zero there). Correct lambda over there to prevent errors.
            lambda = ZEROcr
          endif
#endif
          
          ! For the real stretching boundary conditions , \Sigma for each constitutive variable becomes \Ya\Sigma. It is heavy to change each and every expression where \Sigma arises. Rather, we make use of what is multiplying \Sigma.
          ! Here, in the surface integrations, easiest is the normals. But we have to do it after the call to the 'LNS_get_interfaces_unknowns' routine and the affectation of lambda.
          ! If you change anything to this, remember to do it also in the 'compute_gradient_TFSF' subroutine ('compute_forces_acoustic_LNS_calling_routine.f90').
          if(ABC_STRETCH .and. stretching_buffer(ibool_before_perio(i,j,ispec))>0) then
            iglob_unique = ibool_before_perio(i,j,ispec)
            do SPCDM = 1, NDIM
              n_out(SPCDM) = n_out(SPCDM) * stretching_ya(SPCDM, iglob_unique)
            enddo
          endif
          
          ! 1) Inviscid contributions.
          ! 1.1) Mass conservation (fully inviscid).
          !flux_n = DOT_PRODUCT(n_out, in_dm(:,iglob)+dm_P)
          if(exact_interface_flux) then
            jump = ZEROcr
          else
            jump = cv_drho(iglob) - drho_P
          endif
          outrhs_drho(iglob) =   outrhs_drho(iglob) &
                                - halfWeight*(DOT_PRODUCT(n_out, in_dm(:,iglob)+dm_P) + lambda*jump)
          ! 1.2) x-Momentum inviscid contributions.
          Sigma_L(1)    =   cv_rho0dv(1,iglob)*LNS_v0(1,iglob) + in_dp(iglob) & ! "M" side.
                          + rho0dv_P(1)       *v0_P(1) + dp_P ! "P" side.
          Sigma_L(NDIM) =   cv_rho0dv(NDIM,iglob)*LNS_v0(1,iglob) & ! "M" side.
                          + rho0dv_P(NDIM)       *v0_P(1) ! "P" side.
          if(exact_interface_flux) then
            jump = ZEROcr
          else
            jump = cv_rho0dv(1,iglob) - rho0dv_P(1)
          endif
          outrhs_rho0dv(1,iglob) =   outrhs_rho0dv(1,iglob) &
                                   - halfWeight*(DOT_PRODUCT(n_out, Sigma_L) + lambda*jump)
          ! 1.2) z-Momentum inviscid contributions.
          Sigma_L(1)    =   cv_rho0dv(1,iglob)*LNS_v0(NDIM,iglob) & ! "M" side.
                          + rho0dv_P(1)       *v0_P(NDIM) ! "P" side.
          Sigma_L(NDIM) =   cv_rho0dv(NDIM,iglob)*LNS_v0(NDIM,iglob) + in_dp(iglob) & ! "M" side.
                          + rho0dv_P(NDIM)       *v0_P(NDIM) + dp_P ! "P" side.
          if(exact_interface_flux) then
            jump = ZEROcr
          else
            jump = cv_rho0dv(NDIM,iglob) - rho0dv_P(NDIM)
          endif
          outrhs_rho0dv(NDIM,iglob) =   outrhs_rho0dv(NDIM,iglob) &
                                      - halfWeight*(DOT_PRODUCT(n_out, Sigma_L) + lambda*jump)
          ! 1.3) Energy inviscid contributions.
          Sigma_L =   LNS_dv(:,iglob)*(LNS_E0(iglob) + LNS_p0(iglob)) & ! "M" side, part.
                    + LNS_v0(:,iglob)*(cv_dE(iglob)  + in_dp(iglob)) & ! "M" side, part.
                    + dv_P(:)        *(E0_P          + p0_P) & ! "P" side, part.
                    + v0_P(:)        *(dE_P          + dp_P) ! "P" side, part.
          if(exact_interface_flux) then
            jump = ZEROcr
          else
            jump = cv_dE(iglob) - dE_P
          endif
          outrhs_dE(iglob) =   outrhs_dE(iglob) &
                             - halfWeight*(DOT_PRODUCT(n_out, Sigma_L) + lambda*jump)
          
          ! 2) Viscous contributions.
          if(viscousComputation) then
            ! Note: since we consider here viscous terms, that is purely diffusive phenomena, the acoustic wave speed makes no sense, and thus the jump in the Lax-Friedrich approximation does not appear in the formulations.
            ! 2.1) Mass conservation: nothing.
            ! 2.2) x-Momentum viscous contributions.
            outrhs_rho0dv(1,iglob) =   outrhs_rho0dv(1,iglob) &
                                      - halfWeight*DOT_PRODUCT(n_out, -(in_sigma_dv(1:2,iglob)+sigma_dv_P(1:2)))
            ! 2.2) z-Momentum viscous contributions.
            outrhs_rho0dv(NDIM,iglob) =   outrhs_rho0dv(NDIM,iglob) &
                                         - halfWeight*DOT_PRODUCT(n_out, -(in_sigma_dv(2:3,iglob)+sigma_dv_P(2:3)))
            ! 2.3) Energy viscous contributions.
            Sigma_L(1)    = - (  LNS_dv(1,iglob)*sigma_v_0(1,iglob) & ! "M" side, part.
                               + LNS_v0(1,iglob)*in_sigma_dv(1,iglob) & ! "M" side, part.
                               + LNS_dv(NDIM,iglob)*sigma_v_0(2,iglob) & ! "M" side, part.
                               + LNS_v0(NDIM,iglob)*in_sigma_dv(2,iglob) & ! "M" side, part.
                               + LNS_kappa(iglob)*in_nabla_dT(1,iglob)) & ! "M" side, part.
                            - (  dv_P(1)*sigma_v0_P(1) & ! "P" side, part.
                               + LNS_v0(1,iglob)*sigma_dv_P(1) & ! "P" side, part.
                               + dv_P(NDIM)*sigma_v0_P(2) & ! "P" side, part.
                               + v0_P(NDIM)*sigma_dv_P(2) & ! "P" side, part.
                               + kap_P*nabla_dT_P(1)) ! "P" side, part.
            Sigma_L(NDIM) = - (  LNS_dv(NDIM,iglob)*sigma_v_0(3,iglob) & ! "M" side, part.
                               + LNS_v0(NDIM,iglob)*in_sigma_dv(3,iglob) & ! "M" side, part.
                               + LNS_dv(1,iglob)*sigma_v_0(2,iglob) & ! "M" side, part.
                               + LNS_v0(1,iglob)*in_sigma_dv(2,iglob) & ! "M" side, part.
                               + LNS_kappa(iglob)*in_nabla_dT(2,iglob)) & ! "M" side, part.
                            - (  dv_P(NDIM)*sigma_v0_P(3) & ! "P" side, part.
                               + v0_P(NDIM)*sigma_dv_P(3) & ! "P" side, part.
                               + dv_P(1)*sigma_v0_P(2) & ! "P" side, part.
                               + v0_P(1)*sigma_dv_P(2) & ! "P" side, part.
                               + kap_P*nabla_dT_P(NDIM)) ! "P" side, part.
            outrhs_dE(iglob) = outrhs_dE(iglob) - halfWeight*DOT_PRODUCT(n_out, Sigma_L)
          endif ! Endif on viscousComputation.
        enddo ! Enddo on iface.
      enddo ! Enddo on iface1.
    endif ! End of test if acoustic element
  enddo ! End of loop on elements.
  
end subroutine compute_forces_acoustic_LNS


! ------------------------------------------------------------ !
! LNS_get_bgqts_at_interface                                   !
! ------------------------------------------------------------ !
! This routine complements the 'LNS_get_interfaces_unknowns' routine below.
! This routine returns the background state quantities at the queried point, while 'LNS_get_interfaces_unknowns' returns the constitutive variables.

subroutine LNS_get_bgqts_at_interface(i, j, ispec, iface1, iface, neighbor, neighbour_type, &
                                      iglob_P, &
                                      rho0_P, v0_P, E0_P, p0_P, &
                                      sigma_v0_P, gam_P, kap_P)
  use constants, only: CUSTOM_REAL, NDIM
  use specfem_par, only: myrank, NPROC, MPI_transfer_iface, ibool_DG, &
                         gammaext_DG, buffer_DG_gamma_P
  use specfem_par_lns, only: LNS_viscous, NVALSIGMA, &
                             LNS_rho0, LNS_v0, LNS_E0, LNS_p0, LNS_kappa, sigma_v_0, &
                             buffer_LNS_rho0, buffer_LNS_E0, buffer_LNS_p0, buffer_LNS_kappa, buffer_LNS_v0, buffer_sigma_v0_P
  ! Input/Output.
  integer, intent(in) :: i, j, ispec, iglob_P, iface1, iface
  integer, dimension(3), intent(in) :: neighbor
  integer, intent(in) :: neighbour_type
  real(kind=CUSTOM_REAL), dimension(NDIM), intent(out) :: v0_P
  real(kind=CUSTOM_REAL), intent(out) :: rho0_P, p0_P, E0_P
  real(kind=CUSTOM_REAL), dimension(NVALSIGMA), intent(out) :: sigma_v0_P
  real(kind=CUSTOM_REAL), intent(out) :: gam_P, kap_P
  ! Local.
  real(kind=CUSTOM_REAL), parameter :: ZEROcr = 0._CUSTOM_REAL
  integer :: iglob_M, ipoin, num_interface
  ! Initialisations.
  iglob_M = ibool_DG(i, j, ispec)
  rho0_P = ZEROcr
  v0_P = ZEROcr
  E0_P = ZEROcr
  p0_P = ZEROcr
  gam_P = ZEROcr
  kap_P = ZEROcr
  ipoin = -1
  num_interface = -1
  if(neighbour_type==1) then
    ! A neighbouring LNS DG element was found in the same partition.
    ! Simply grab values from iglob_P.
    rho0_P = LNS_rho0(iglob_P)
    v0_P = LNS_v0(:, iglob_P)
    E0_P = LNS_E0(iglob_P)
    p0_P = LNS_p0(iglob_P)
    gam_P = gammaext_DG(iglob_P)
    if(LNS_viscous) then
      kap_P = LNS_kappa(iglob_P)
      sigma_v0_P = sigma_v_0(:, iglob_P)
    endif
    
  elseif(neighbour_type==11) then
    if(NPROC <= 1) then
      call exit_MPI(myrank, 'Should not have neighbour_type=11 with NPROC<=1.')
    else
      ipoin         = MPI_transfer_iface(iface1, iface, ispec, 1)
      num_interface = MPI_transfer_iface(iface1, iface, ispec, 2)
    endif
    ! A neighbouring LNS DG element was found in another partition.
    ! Use buffers.
    rho0_P = buffer_LNS_rho0(ipoin, num_interface)
    v0_P = buffer_LNS_v0(:, ipoin, num_interface)
    E0_P = buffer_LNS_E0(ipoin, num_interface)
    p0_P = buffer_LNS_p0(ipoin, num_interface)
    gam_P = buffer_DG_gamma_P(ipoin, num_interface)
    if(LNS_viscous) then
      kap_P = buffer_LNS_kappa(ipoin, num_interface)
      sigma_v0_P = buffer_sigma_v0_P(:, ipoin, num_interface)
    endif
    
  elseif(neighbour_type==21 .or. neighbour_type==99) then
    ! Other material viscoelastic or outside boundary.
    ! Neumann condition.
    rho0_P = LNS_rho0(iglob_M)
    v0_P = LNS_v0(:, iglob_M)
    E0_P = LNS_E0(iglob_M)
    p0_P = LNS_p0(iglob_M)
    gam_P = gammaext_DG(iglob_M)
    if(LNS_viscous) then
      kap_P = LNS_kappa(iglob_M)
      sigma_v0_P = sigma_v_0(:, iglob_M)
    endif
    
  else
    call exit_MPI(myrank, 'Unexpected value for neighbour_type.')
    
  endif
end subroutine LNS_get_bgqts_at_interface


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

subroutine LNS_get_interfaces_unknowns(i, j, ispec, iface1, iface, neighbor, neighbour_type, timelocal, & ! Point identifier (input).
             inp_drho_M, inp_rho0dv_M, & ! Input constitutive variables, "M" side.
             inp_drho_P, inp_rho0dv_P, inp_dE_P, & ! Input constitutive variables, "P" side. They make no sense if neighbor(1)<=-1.
             inp_dp_M, & ! Input other variable, "M" side.
             n_out, & ! Normal vector (input).
             exact_interface_flux, & ! Switch to disable jump in some cases (output).
             out_drho_P, out_rho0dv_P, out_dE_P, & ! Output constitutive variables.
             out_dm_P, out_dp_P, out_dv_P, & ! Output other variables.
             swCompVisc, out_nabla_dT_P, out_sigma_dv_P, & ! Output other variables: viscous.
             swCompdT, out_dT_P) ! Output other variables.
  
  use constants, only: CUSTOM_REAL, NDIM, PI
  use specfem_par, only: NPROC, ibool, ibool_DG, acoustic_forcing, myrank, &!coord, &
                         veloc_elastic, sigma_elastic, &
                         mpi_transfer_iface, &
                         ibool_before_perio, & ! For MMS validation.
                         ispec_is_acoustic_coupling_el, ispec_is_acoustic_forcing
  use specfem_par_LNS, only: NVALSIGMA, LNS_dummy_1d, &
                             LNS_rho0, LNS_v0, LNS_E0, LNS_p0, LNS_c0, LNS_dv, nabla_dT, sigma_dv, &
                             VALIDATION_MMS, &
                             buffer_LNS_drho_P, buffer_LNS_rho0dv_P, buffer_LNS_dE_P, buffer_LNS_nabla_dT, buffer_LNS_sigma_dv
  
  implicit none
  
  ! Input/Output.
  integer, intent(in) :: i, j, ispec, iface1, iface
  integer, dimension(3), intent(in) :: neighbor
  integer, intent(inout) :: neighbour_type ! May be modified (see below).
  real(kind=CUSTOM_REAL), intent(in) :: timelocal
  real(kind=CUSTOM_REAL), intent(in) :: inp_drho_M, inp_drho_P, inp_dE_P ! Input constitutive variables. They make no sense if neighbor(1)<=-1.
  real(kind=CUSTOM_REAL), dimension(NDIM), intent(in) :: inp_rho0dv_M, inp_rho0dv_P ! Input constitutive variables.
  real(kind=CUSTOM_REAL), intent(in) :: inp_dp_M ! Input other variables.
  real(kind=CUSTOM_REAL), dimension(NDIM), intent(in) :: n_out
  logical, intent(out) :: exact_interface_flux ! Output switch.
  real(kind=CUSTOM_REAL), intent(out) :: out_drho_P, out_dE_P ! Output constitutive variables.
  real(kind=CUSTOM_REAL), dimension(NDIM), intent(out) :: out_rho0dv_P ! Output constitutive variable.
  real(kind=CUSTOM_REAL), dimension(NDIM), intent(out) :: out_dv_P, out_dm_P, out_nabla_dT_P ! Output other variables.
  real(kind=CUSTOM_REAL), intent(out) :: out_dp_P ! Output other variables.
  real(kind=CUSTOM_REAL), dimension(NVALSIGMA), intent(out) :: out_sigma_dv_P ! Output other variables.
  real(kind=CUSTOM_REAL), intent(out) :: out_dT_P ! In compute_gradient_TFSF (desintegration method), we need temperature on the other side of the boundary in order to compute the flux.
  logical, intent(in) :: swCompVisc, swCompdT!, swCompv ! Do not unnecessarily compute some quantities.
  
  ! Local.
  real(kind=CUSTOM_REAL), parameter :: ZEROcr = 0._CUSTOM_REAL
  integer :: SPCDM
  real(kind=CUSTOM_REAL), dimension(NDIM) :: velocity_P
  real(kind=CUSTOM_REAL), dimension(NDIM) :: tang ! Tangential vector.
  real(kind=CUSTOM_REAL) :: normal_v, tangential_v
  integer :: iglobM, i_el, j_el, ispec_el, iglobP, ipoin, num_interface
  real(kind=CUSTOM_REAL), dimension(NDIM, NDIM) :: trans_boundary
  
  !if(.false.) write(*,*) coord, ibool_before_perio ! Only here because compiler is going to râler because of imports necessary to MMS validation.
  
  ! Initialise output variables to default values.
  exact_interface_flux = .false.
  out_drho_P     = ZEROcr
  out_rho0dv_P   = ZEROcr
  out_dE_P       = ZEROcr
  out_dv_P       = ZEROcr
  out_dm_P       = ZEROcr
  out_dp_P       = ZEROcr
  out_nabla_dT_P = ZEROcr
  out_sigma_dv_P = ZEROcr
  out_dT_P       = ZEROcr
  velocity_P     = ZEROcr
  
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
  
  iglobM = ibool_DG(i, j, ispec) ! Extract calling point iglob on the "minus" side.
  
  ! Characterise the type of neighbour.
  ! The definition of each type of neighbour is detailed alongside the 'check_neighbour_type' routine ('compute_forces_acoustic_LNS_calling_routine.f90').
  if(neighbour_type >= 10) then
    ipoin         = -1 ! By default, specify that for the triplet (iface1,iface,ispec), the values should not be sought in another partition. That is, either the edge is a coupling (with another material) edge, or an edge on the outer boundary of the computational domain.
    num_interface = -1 ! Initialised to this value, but should not be used anywhere as is.
    if(NPROC > 1) then
      ipoin         = MPI_transfer_iface(iface1, iface, ispec, 1)
      num_interface = MPI_transfer_iface(iface1, iface, ispec, 2)
    endif
    if(ipoin > -1) then
      neighbour_type = 11 ! A neighbouring LNS DG element was found in another partition.
    else
      if(ispec_is_acoustic_coupling_el(i, j, ispec, 3) >= 0) then
        neighbour_type = 21 ! Other material: viscoelastic.
      elseif(ACOUSTIC_FORCING .AND. ispec_is_acoustic_forcing(i, j, ispec)) then
        neighbour_type = 22 ! Other material: potential acoustic (not implemented).
      else
        neighbour_type = 99 ! Outside boundary.
      endif
    endif
  else
    ! Unchanged: (neighbour_type < 10) <=> (neighbour_type == 1) since no types defined in the 2:10 range.
    neighbour_type = 1 ! A neighbouring LNS DG element was found in the same partition.
  endif
  
  ! Depending on the type of neighbour, set the boundary conditionz.
  !--------------------------------------------------------------
  if(neighbour_type==11) then
    exact_interface_flux = .false. ! Set exact_interface_flux.
    out_drho_P     = buffer_LNS_drho_P(ipoin, num_interface) ! Set out_drho_P.
    out_rho0dv_P   = buffer_LNS_rho0dv_P(:, ipoin, num_interface) ! Set out_rho0dv_P.
    out_dE_P       = buffer_LNS_dE_P(ipoin, num_interface) ! Set out_dE_P.
    out_dv_P(:) = out_rho0dv_P(:)/LNS_rho0(iglobM) ! Set out_dv_P.
    call compute_dp_i(LNS_rho0(iglobM)+out_drho_P, LNS_v0(:,iglobM)+out_dv_P, LNS_E0(iglobM)+out_dE_P, out_dp_P, iglobM) ! Set out_dp_P.
    ! Set out_dm_P: see bottom of routine.
    if(swCompVisc) then
      out_nabla_dT_P = buffer_LNS_nabla_dT(:, ipoin, num_interface) ! Set out_nabla_dT_P: get the values from the MPI buffers.
      out_sigma_dv_P = buffer_LNS_sigma_dv(:, ipoin, num_interface) ! Set out_sigma_dv_P: get the values from the MPI buffers.
    endif
    if(swCompdT) then
      call compute_dT_i(LNS_rho0(iglobM)+out_drho_P, LNS_v0(:, iglobM)+out_dv_P, LNS_E0(iglobM)+out_dE_P, out_dT_P, iglobM) ! Set out_dT_P.
    endif
  
  !--------------------------------------------------------------  
  elseif(neighbour_type==22) then
    write(*,*) "********************************"
    write(*,*) "*            ERROR             *"
    write(*,*) "********************************"
    write(*,*) "* ACOUSTIC_FORCING is obsolete *"
    write(*,*) "* for DG simulations.          *"
    write(*,*) "********************************"
    stop
  
  !--------------------------------------------------------------
  elseif(neighbour_type==21) then
    exact_interface_flux = .false. ! Set exact_interface_flux.
    ! Coordinates of elastic element
    i_el     = ispec_is_acoustic_coupling_el(i, j, ispec, 1)
    j_el     = ispec_is_acoustic_coupling_el(i, j, ispec, 2)
    ispec_el = ispec_is_acoustic_coupling_el(i, j, ispec, 3)
    
    out_drho_P = inp_drho_M ! Set out_drho_P: same as other side, that is a Neumann condition.
    
    ! Set velocity_P.
    call build_trans_boundary(n_out, tang, trans_boundary)
!#define USECLASSICALCOUPLING 1
#define USECLASSICALCOUPLING 0
#if USECLASSICALCOUPLING
    ! VERSION 1: normal velocity from elastic velocity, and slip condition for tangential.
    normal_v     = DOT_PRODUCT(n_out, veloc_elastic(:,ibool(i_el, j_el, ispec_el)))
    tangential_v = DOT_PRODUCT(tang, LNS_v0(:,iglobM)+LNS_dv(:,iglobM))
    do SPCDM = 1, NDIM
      velocity_P(SPCDM) = trans_boundary(SPCDM, 1)*normal_v + trans_boundary(SPCDM, 2)*tangential_v
    enddo
#else
    ! VERSION 2: normal velocity from Terrana velocity, and slip condition for tangential.
    call S2F_Terrana_coupling(n_out, &
                              LNS_rho0(iglobM)+inp_drho_M, &
                              LNS_v0(:,iglobM)+LNS_dv(:,iglobM), &
                              inp_dp_M, &
                              LNS_c0(iglobM), &
                              iglobM, &
                              veloc_elastic(:,ibool(i_el, j_el, ispec_el)), &
                              sigma_elastic(:,:,ibool(i_el, j_el, ispec_el)), &
                              i_el, j_el, ispec_el, &
                              velocity_P, out_dp_P)
    ! Set out_dv_P: do as in FNS, put Terrana's velocity in normal component, and leave tangential untouched (slip condition).
    call build_trans_boundary(n_out, tang, trans_boundary)
    normal_v     = DOT_PRODUCT(n_out, velocity_P)
    tangential_v = DOT_PRODUCT(tang, LNS_v0(:,iglobM)+LNS_dv(:,iglobM))
    do SPCDM = 1, NDIM
      velocity_P(SPCDM) = trans_boundary(SPCDM, 1)*normal_v + trans_boundary(SPCDM, 2)*tangential_v
    enddo
#endif
    ! Set out_dv_P.
    out_dv_P = velocity_P - LNS_v0(:, iglobM) ! Requesting iglobM might be technically inexact, but on elements' boundaries points should overlap. Plus, iglobP does not exist on outer computational domain boundaries.
    
    ! Set out_dp_P.
#if USECLASSICALCOUPLING
    ! VERSION 1: no stress continuity.
    out_dp_P = inp_dp_M
#else
    ! VERSION 2: TERRANA.
    ! done above in the call to S2F_Terrana_coupling
#endif

    ! Set out_dE_P.
    call compute_dE_i(LNS_rho0(iglobM)+out_drho_P, velocity_P, LNS_p0(iglobM)+out_dp_P, out_dE_P, iglobM)

    ! Set out_rho0dv_P.
    out_rho0dv_P(:) = LNS_rho0(iglobM)*out_dv_P(:)
    
    ! Set out_dm_P: see bottom of routine.
      
    if(swCompVisc) then
      out_nabla_dT_P = nabla_dT(:,iglobM) ! Set out_nabla_dT_P: same as other side, that is a Neumann condition.
      out_sigma_dv_P = sigma_dv(:,iglobM) ! Set out_sigma_dv_P: same as other side, that is a Neumann condition.
    endif
    if(swCompdT) then
      call compute_dT_i(LNS_rho0(iglobM)+out_drho_P, velocity_P, LNS_E0(iglobM)+out_dE_P, out_dT_P, iglobM) ! Set out_dT_P.
    endif
    
  !--------------------------------------------------------------
  elseif(neighbour_type==99) then
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
    
    if(VALIDATION_MMS) call VALIDATION_MMS_boundary_terms(ibool_before_perio(i, j, ispec), iglobM, &
                                                          exact_interface_flux, out_drho_P, out_dv_P, out_dE_P, &
                                                          out_dp_P, out_rho0dv_P, &
                                                          swCompVisc, out_nabla_dT_P, out_sigma_dv_P, &
                                                          swCompdT, out_dT_P)
!      endif ! Endif on PML.
!    endif ! Endif on ipoin.
  
#if 0
! The case below will stop the program anyhow afterwards. Best remove it altogether for now.
! Moreover, ispec_is_acoustic_coupling_ac is allocated in a bad place (in compute_forces_acoustic_DG_calling_routine.F90), and so testing it causes a segfault. So, best remove it altogether for now.
! You may have to bring it back onboard if you plan on implementing acoustic_DG/acoustic coupling.  
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
    !  write(*,*) ibool(i_ac,j_ac,ispec_ac), xixl, xizl, gammaxl, gammazl, duz_dxi, duz_dgamma, k
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
#endif
    
!  else
  elseif(neighbour_type==1) then
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
      out_nabla_dT_P = nabla_dT(:,iglobP) ! Set out_nabla_dT_P.
      out_sigma_dv_P = sigma_dv(:,iglobP) ! Set out_sigma_dv_P.
      !out_nabla_dT_P = nabla_dT(:,iglobM) ! Set out_nabla_dT_P.
      !out_sigma_dv_P = sigma_dv(:,iglobM) ! Set out_sigma_dv_P.
    endif
    
    ! Set out_dT_P.
    if(swCompdT) then
      ! iglobP exists and has meaning.
      call compute_dT_i(LNS_rho0(iglobP)+out_drho_P, LNS_v0(:, iglobP)+out_dv_P, LNS_E0(iglobP)+out_dE_P, out_dT_P, iglobP)
      !call compute_dT_i(LNS_rho0(iglobP)+out_drho_P, LNS_p0(iglobP)+out_dp_P, out_dT_P, iglobP)
    endif
    !out_dT_P = (out_dE_P/out_drho_P - 0.5*((out_rho0dv_P(1)/out_drho_P)**2 + (out_rho0dv_P(NDIM)/out_drho_P)**2))/c_V
  
  else
    ! Safeguard.
    call exit_MPI(myrank, 'Unexpected value for neighbour_type.')
    
  endif ! Endif on neighbour_type.
  
  ! Set out_dm_P.
  out_dm_P = out_rho0dv_P + out_drho_P*LNS_v0(:,iglobM)
#if 0
! DEBUG
          if(      abs(coord(1,ibool_before_perio(i,j,ispec))+200.)<=5.&
             .and. (     abs(coord(2,ibool_before_perio(i,j,ispec))-200.)<=3. &
                    .or. abs(coord(2,ibool_before_perio(i,j,ispec))-233.)<=3.) &
             .and. abs(timelocal-.9689999999999)<=0.00001 &
             .and. myrank==0) then
          write(*,*) myrank, 'X', coord(:,ibool_before_perio(i,j,ispec)), &
                     'dm other side', out_dm_P
        endif
#endif
  
end subroutine LNS_get_interfaces_unknowns



! ------------------------------------------------------------ !
! build_trans_boundary                                         !
! ------------------------------------------------------------ !
! From an input normal vector, compute the tangential vector and the matrix for passin from one reference frame to the other.
! The matrix transf_matrix converts a vector from the normal reference frame (n, t) to the canonic reference frame (x, z).
! 
subroutine build_trans_boundary(normal, tangential, transf_matrix)
  use constants, only: CUSTOM_REAL, NDIM
  
  implicit none
  
  ! Input/Output.
  real(kind=CUSTOM_REAL), dimension(NDIM), intent(in) :: normal
  real(kind=CUSTOM_REAL), dimension(NDIM), intent(out) :: tangential
  real(kind=CUSTOM_REAL), dimension(NDIM, NDIM), intent(out) :: transf_matrix
  
  ! Local.
  ! N./A.
  
  ! The convention for the tangential vector is as follows.
  !             ^
  !             | normal
  !             |
  ! .....<------+............ interface
  ! tangential
  ! Though non-standard, it poses no issue and is coherent with DG implementation.
  tangential(1)    = -normal(NDIM) ! Recall: normal(NDIM)=n_z.
  tangential(NDIM) =  normal(1) ! Recall: normal(1)=n_x.
  
  ! Build the transformation matrix.
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
  !use specfem_par_lns, only: inverse_2x2, LNS_v0, LNS_dv!, LNS_p0!, LNS_viscous, sigma_dv
  
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
  call stressBuilder_addInviscidFluid(LNS_rho0(iglob_DG), LNS_v0(:,iglob_DG), & ! ORIGINAL LINE
  !call stressBuilder_addInviscidFluid(0., LNS_v0(:,iglob_DG), & ! Typically, this builds simply the pressure.
                                      LNS_dv(:,iglob_DG), dp_fluid, sigma_fluid) ! ORIGINAL LINE
  !                                    LNS_dv(:,iglob_DG), LNS_p0(iglob_DG) + dp_fluid, sigma_fluid) ! Take full pressure.
  !sigma_fluid(1,2) = 0.; sigma_fluid(2,1) = 0.; ! Force tangent stress to zero.
  !sigma_fluid = 0. ! Force whole stress to zero (test purposes).
  !sigma_fluid(1,1) = 
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
#if 0
  write(*,*) 'v_solid', v_solid, 'v_fluid', v_fluid ! DEBUG
  write(*,*) 'TAU_F', TAU_F
  write(*,*) 'matmul(TAU_F, v_fluid)', matmul(TAU_F, v_fluid), 'matmul(TAU_S, v_solid)', matmul(TAU_S, v_solid)
  write(*,*) 'sigma_el_local - sigma_fluid', (sigma_el_local - sigma_fluid)
  write(*,*) 'v_hat', v_hat
#endif
  
  ! Build actual pressure perturbation from [Terrana et al., 2018]'s (51).
  dp_hat = dp_fluid + DOT_PRODUCT(matmul(TAU_F, v_fluid - v_solid), normal)
  !dp_hat = dp_fluid + DOT_PRODUCT(matmul(TAU_F, v_fluid - v_hat), normal)
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
  vp = sqrt((poroelastcoef(1, 1, kmato(ispec_el)) + 2._CUSTOM_REAL*mu_elastic_unrelaxed)/rho) ! vp = ((lambda+2mu)/rho)^0.5 http://www.subsurfwiki.org/wiki/P-wave_modulus, poroelastcoef(1,1,kmato(ispec_el)) = lambda
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
              !d_dE(iglob) =   d_dE(iglob) &
              !             + (   LNS_c0(iglob)**2/(gammaext_DG(iglob)-1._CUSTOM_REAL) &
              !                 + 0.5_CUSTOM_REAL*norm2r1(LNS_v0(:, iglob)) &
              !               ) * temp_sourcewxlwzljacobianl
              ! Approx. c^2=c_0^2. TBTested.
              !d_dE(iglob) =   d_dE(iglob) &
              !             + (   LNS_c0(iglob)**2/(gammaext_DG(iglob)-1.) &
              !                 + 0.5*norm2r1(LNS_v0(:, iglob)+LNS_dv(:, iglob)) &
              !               ) * temp_sourcewxlwzljacobianl
              ! Full. Highest CPU use, but technically exactly what the equations are, same results as "Approx. c^2=c_0^2 & v=v_0".
              d_dE(iglob) =   d_dE(iglob) &
                           + (   gammaext_DG(iglob) * (LNS_p0(iglob)+LNS_dp(iglob)) &
                                 / ((LNS_rho0(iglob)+LNS_drho(iglob))*(gammaext_DG(iglob)-1._CUSTOM_REAL)) &
                               + 0.5_CUSTOM_REAL*norm2r1(LNS_v0(:, iglob)+LNS_dv(:, iglob)) &
                             ) * temp_sourcewxlwzljacobianl
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











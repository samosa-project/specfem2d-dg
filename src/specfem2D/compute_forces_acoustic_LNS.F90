! ------------------------------------------------------------ !
! compute_forces_acoustic_LNS                                  !
! ------------------------------------------------------------ !
! TODO: Description.

subroutine compute_forces_acoustic_LNS(cv_drho, cv_rho0dv, cv_dE, & ! Constitutive variables.
                                       in_dm, in_dp, in_dT, in_nabla_dT, in_sigma_dv, & ! Precomputed quantities.
                                       outrhs_drho, outrhs_rho0dv, outrhs_dE, & ! Output (RHS for each constitutive variable).
                                       currentTime) ! Time.
  ! TODO: select variables to use.
  use constants!, only: CUSTOM_REAL,NGLLX,NGLLZ,SPACEDIM
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
  real(kind=CUSTOM_REAL), dimension(nglob_DG), intent(in) :: cv_drho, cv_dE, in_dp, in_dT
  real(kind=CUSTOM_REAL), dimension(SPACEDIM, nglob_DG), intent(in) :: cv_rho0dv, in_dm, in_nabla_dT
  real(kind=CUSTOM_REAL), dimension(3, nglob_DG), intent(in) :: in_sigma_dv
  real(kind=CUSTOM_REAL), dimension(nglob_DG), intent(out) :: outrhs_drho, outrhs_dE
  real(kind=CUSTOM_REAL), dimension(SPACEDIM, nglob_DG), intent(out) :: outrhs_rho0dv
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
  integer :: ispec, i,j, k,iglob,SPCDM!, iglob_unique
  real(kind=CUSTOM_REAL), dimension(NGLLX, NGLLZ, 2) :: cntrb_drho, cntrb_dE
  real(kind=CUSTOM_REAL), dimension(SPACEDIM, NGLLX, NGLLZ, 2) :: cntrb_rho0dv
  ! Variables "cntrb_*" are aimed at assembling the different contributions to constitutive variables.
  ! For drho and dE, "cntrb_*(:,:,1)" contains the contribution along xi, and "cntrb_*(:,:,2)" the contribution along gamma.
  ! For rho0dv, "cntrb_*(1,:,:,1)" contains the contribution to rho0dvx along xi, "cntrb_*(1,:,:,2)" the contribution to rho0dvx along gamma, "cntrb_*(2,:,:,1)" the contribution to rhodvz along xi, and "cntrb_*(2,:,:,2)" the contribution to rhodvz along gamma.
  real(kind=CUSTOM_REAL), dimension(NGLLX, NGLLZ) :: d0cntrb_dE
  real(kind=CUSTOM_REAL), dimension(SPACEDIM, NGLLX,NGLLZ) :: d0cntrb_rho0dv
  ! Variable "d0cntrb_*" are aimed at assembling contributions of zero-th degree (outside divergence operator).
  
  
  real(kind=CUSTOM_REAL) :: xixl,xizl,gammaxl,gammazl,jacLoc ! Jacobian matrix and determinant.
  real(kind=CUSTOM_REAL) :: lambda, weight, &
                            !tmp_unknown_x_M, tmp_unknown_x_P, tmp_unknown_z_M, tmp_unknown_z_P, &
                            tmp_unknown_x, tmp_unknown_z, &
                            flux_n, jump!, &
                            !veloc_x_DG_P, veloc_z_DG_P!, &
                            !in_dp_P!, &!currentTime,T_P,  &
                            !Tx_DG_P, Tz_DG_P, Vxx_DG_P, Vzz_DG_P, Vxz_DG_P, Vzx_DG_P, & ! MIGHT BE NEEDED
                            ! TEST
                            !gamma_P!,e1_DG_P
  
  real(kind=CUSTOM_REAL), dimension(SPACEDIM) :: n_out
  real(kind=CUSTOM_REAL), dimension(SPACEDIM) :: rho0dv_P
  !real(kind=CUSTOM_REAL) :: dT_dx, dT_dz
  !real(kind=CUSTOM_REAL) :: veloc_n_M, veloc_n_P
  integer :: iglobM, iglobP
  integer, dimension(3) :: neighbor
  !real(kind=CUSTOM_REAL), dimension(nglob_DG) :: veloc_x_DG, veloc_z_DG, in_dp
  !real(kind=CUSTOM_REAL), dimension(SPACEDIM, nglob_DG) :: LNS_dv
  
  ! Variables specifically for LNS_get_interfaces_unknowns.
  real(kind=CUSTOM_REAL), dimension(SPACEDIM) :: dv_P, dm_P, nabla_dT_P
  real(kind=CUSTOM_REAL) :: drho_P, dp_P, dE_P
  real(kind=CUSTOM_REAL), dimension(NVALSIGMA) :: sigma_dv_P
  logical :: viscousComputation
  !real(kind=CUSTOM_REAL) :: dux_dx, dux_dz, duz_dx, duz_dz ! Viscosity
  real(kind=CUSTOM_REAL) :: wxl, wzl ! Integration weigths.
  logical :: exact_interface_flux
  !integer :: cnst_hdrsttc ! For computationnaly better CONSTRAIN_HYDROSTATIC switches.
  integer :: iface1, iface, iface1_neighbor, iface_neighbor, ispec_neighbor
  !real(kind=CUSTOM_REAL) :: ya_x_l, ya_z_l ! Stretching absorbing boundary conditions.
  
  !T_DG = T_DG_main
  !V_DG = V_DG_main
  
  ! Initialise auxiliary unknowns from constitutive variables.
  ! WAS ALREADY PRECOMPUTED IN CALLING ROUTINE
  !LNS_dv = 0.
  !do SPCDM=1,SPACEDIM
  !  LNS_dv(SPCDM,:) = cv_rho0dv(SPCDM,:)/LNS_rho0(:) ! Element-wise division should occur naturally.
  !enddo
  !p_DG       = (gammaext_DG - ONE)*( cv_dE &
  !             - (HALFcr)*cv_drho*( veloc_x_DG**2 + veloc_z_DG**2 ) )
  
  ! Initialisation of the RHS.
  outrhs_drho    = ZEROcr
  outrhs_rho0dv  = ZEROcr
  outrhs_dE      = ZEROcr
  d0cntrb_rho0dv = ZEROcr
  d0cntrb_dE     = ZEROcr
  
  ! Start by adding source terms.
  select case (TYPE_SOURCE_DG)
    case (1)
      call compute_add_sources_acoustic_DG_spread(outrhs_drho, it, i_stage)
    case (2)
      do SPCDM=1,SPACEDIM
        call compute_add_sources_acoustic_DG_spread(outrhs_rho0dv(SPCDM,:), it, i_stage)
      enddo
    case (3)
      call compute_add_sources_acoustic_DG_spread(outrhs_dE, it, i_stage)
    case default
      stop "TYPE_SOURCE_DG not implemented."
  end select
  
  if(myrank == 0 .and. LNS_VERBOSE>=51 .and. mod(it, LNS_MODPRINT) == 0) then
    write(*,"(a,i6,a)") "Informations for process number ", myrank, "."
    write(*,"(a)")                   " quantity [                 max    ,                  min    ]"
    WRITE(*,"(a,e24.16,a,e24.16,a)") " drho     [", maxval(cv_drho), ", ", minval(cv_drho), "]"
    do SPCDM=1,SPACEDIM
      WRITE(*,"(a,e24.16,a,e24.16,a,i1)") " rho0dv_i [", maxval(cv_rho0dv(SPCDM,:)), ", ", minval(cv_rho0dv(SPCDM,:)), "], i=",SPCDM
    enddo
    WRITE(*,"(a,e24.16,a,e24.16,a)") " dE       [", maxval(cv_dE), ", ", minval(cv_dE), "]"
    !WRITE(*,"(a,e23.16,a)")        "Ratio |p-p_{init}|/p_{init}:", maxval(abs((p_DG-p_DG_init)/p_DG_init)), "."
  endif
  
  do ispec = 1, nspec ! Loop over elements.
    if (ispec_is_acoustic_DG(ispec)) then ! Only do something for DG elements.
      ! --------------------------- !
      ! First set of loops: compute !
      ! volumic contributions.      !
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
          
          if(LNS_viscous) then ! Check if viscosity exists whatsoever.
            ! Activate/deactivate, for this particular point (iglob), computation of quantities only needed when viscosity is present.
            if(     LNS_mu(iglob) > 0. &
               .OR. LNS_eta(iglob) > 0. &
               .OR. LNS_kappa(iglob) > 0.) then
              viscousComputation=.true.
            else
              viscousComputation=.false.
            endif
          else
            viscousComputation=.false. ! If viscosity is globally disabled, deactivate it for this element.
          endif
          
          ! Inviscid stress tensor's contributions.
          !   Mass conservation.
          tmp_unknown_x = in_dm(1,iglob)
          tmp_unknown_z = in_dm(SPACEDIM,iglob)
          cntrb_drho(i,j,1) = wzl * jacLoc * (xixl * tmp_unknown_x + xizl * tmp_unknown_z) ! Contribution along xi.
          cntrb_drho(i,j,2) = wxl * jacLoc * (gammaxl * tmp_unknown_x + gammazl * tmp_unknown_z) ! Contribution along gamma.
          !   x-Momentum.
          tmp_unknown_x = cv_rho0dv(1,iglob)*LNS_v0(1,iglob) + in_dp(iglob)
          tmp_unknown_z = cv_rho0dv(SPACEDIM,iglob)*LNS_v0(1,iglob)
          cntrb_rho0dv(1,i,j,1) = wzl * jacLoc * (xixl * tmp_unknown_x + xizl * tmp_unknown_z) ! Contribution along xi.
          cntrb_rho0dv(1,i,j,2) = wxl * jacLoc * (gammaxl * tmp_unknown_x + gammazl * tmp_unknown_z) ! Contribution along gamma.
          !   z-Momentum.
          tmp_unknown_x = cv_rho0dv(1,iglob)*LNS_v0(SPACEDIM,iglob)
          tmp_unknown_z = cv_rho0dv(SPACEDIM,iglob)*LNS_v0(SPACEDIM,iglob) + in_dp(iglob)
          cntrb_rho0dv(SPACEDIM,i,j,1) = wzl * jacLoc * (xixl * tmp_unknown_x + xizl * tmp_unknown_z) ! Contribution along xi.
          cntrb_rho0dv(SPACEDIM,i,j,2) = wxl * jacLoc * (gammaxl * tmp_unknown_x + gammazl * tmp_unknown_z) ! Contribution along gamma.
          
          ! Add viscous stress tensor's contributions.
          if(viscousComputation) then
            !   Mass conservation: no viscous contribution.
            !   x-Momentum.
            tmp_unknown_x = -in_sigma_dv(1,iglob) ! Recall, this corresponds to \Sigma_v(v')_{1,1}.
            tmp_unknown_z = -in_sigma_dv(2,iglob) ! Recall, this corresponds to \Sigma_v(v')_{1,2} (and \Sigma_v(v')_{2,1}).
            cntrb_rho0dv(1,i,j,1) = cntrb_rho0dv(1,i,j,1) + wzl * jacLoc * (xixl * tmp_unknown_x + xizl * tmp_unknown_z) ! Contribution along xi.
            cntrb_rho0dv(1,i,j,2) = cntrb_rho0dv(1,i,j,2) + wxl * jacLoc * (gammaxl * tmp_unknown_x + gammazl * tmp_unknown_z) ! Contribution along gamma.
            !   z-Momentum.
            tmp_unknown_x = -in_sigma_dv(2,iglob) ! Recall, this corresponds to \Sigma_v(v')_{2,1} (and \Sigma_v(v')_{1,2}).
            tmp_unknown_z = -in_sigma_dv(3,iglob) ! Recall, this corresponds to \Sigma_v(v')_{2,2}.
            cntrb_rho0dv(SPACEDIM,i,j,1) =   cntrb_rho0dv(SPACEDIM,i,j,1) &
                                           + wzl * jacLoc * (xixl * tmp_unknown_x + xizl * tmp_unknown_z) ! Contribution along xi.
            cntrb_rho0dv(SPACEDIM,i,j,2) =   cntrb_rho0dv(SPACEDIM,i,j,2) &
                                           + wxl * jacLoc * (gammaxl * tmp_unknown_x + gammazl * tmp_unknown_z) ! Contribution along gamma.
          endif ! Endif on viscousComputation.
          
          ! Special case: energy. Indeed, if they exist, viscous contributions can be grouped to inviscid ones.
          if(viscousComputation) then
            tmp_unknown_x =   LNS_dv(1,iglob)*(LNS_E0(iglob) + LNS_p0(iglob) - sigma_v_0(1,iglob)) &
                            + LNS_v0(1,iglob)*(cv_dE(iglob) + in_dp(iglob) - in_sigma_dv(1,iglob)) &
                            - LNS_dv(SPACEDIM,iglob)*sigma_v_0(2,iglob) &
                            - LNS_v0(SPACEDIM,iglob)*in_sigma_dv(2,iglob) &
                            + LNS_kappa(iglob)*in_nabla_dT(1,iglob)
            tmp_unknown_z =   LNS_dv(SPACEDIM,iglob)*(LNS_E0(iglob) + LNS_p0(iglob) - sigma_v_0(3,iglob)) &
                            + LNS_v0(SPACEDIM,iglob)*(cv_dE(iglob) + in_dp(iglob) - in_sigma_dv(3,iglob)) &
                            - LNS_dv(1,iglob)*sigma_v_0(2,iglob) &
                            - LNS_v0(1,iglob)*in_sigma_dv(2,iglob) &
                            + LNS_kappa(iglob)*in_nabla_dT(SPACEDIM,iglob)
          else ! Else on viscousComputation.
            tmp_unknown_x =   LNS_dv(1,iglob)*(LNS_E0(iglob) + LNS_p0(iglob)) &
                            + LNS_v0(1,iglob)*(cv_dE(iglob) + in_dp(iglob))
            tmp_unknown_z =   LNS_dv(SPACEDIM,iglob)*(LNS_E0(iglob) + LNS_p0(iglob)) &
                            + LNS_v0(SPACEDIM,iglob)*(cv_dE(iglob) + in_dp(iglob))
          endif ! Endif on viscousComputation.
          cntrb_dE(i,j,1) = wzl * jacLoc * (xixl * tmp_unknown_x + xizl * tmp_unknown_z) ! Contribution along xi.
          cntrb_dE(i,j,2) = wxl * jacLoc * (gammaxl * tmp_unknown_x + gammazl * tmp_unknown_z) ! Contribution along gamma.
          
          ! Zero-th degree contributions.
          ! > Mass conservation: none.
          !d0cntrb_drho(i,j)     = ZERO
          ! > Momenta.
          !   Version 1: most general.
          !d0cntrb_rho0dv(1,i,j) = (  cv_drho(iglob)*potential_dphi_dx_DG(ibool(i,j,ispec)) & ! \rho'g_x
          !                         + in_dm(1,iglob)*nabla_v0(1,1,iglob) & ! {\delta_m}_x\partial_xv_{0,x}
          !                         + in_dm(SPACEDIM,iglob)*nabla_v0(1,SPACEDIM,iglob)) * jacLoc ! {\delta_m}_z\partial_zv_{0,x}
          !d0cntrb_rho0dv(SPACEDIM,i,j) = (  cv_drho(iglob)*potential_dphi_dz_DG(ibool(i,j,ispec)) & ! \rho'g_z
          !                                + in_dm(1,iglob)*nabla_v0(SPACEDIM,1,iglob) & ! {\delta_m}_x\partial_xv_{0,z}
          !                                + in_dm(SPACEDIM,iglob)*nabla_v0(SPACEDIM,SPACEDIM,iglob)) * jacLoc ! {\delta_m}_z\partial_zv_{0,z}
          !   Version 2: under HV0 (v_{0,z}=0), and SM (d_xv_0=0), dm part of the zero-th degree RHS is simplified.
          !d0cntrb_rho0dv(1,i,j) = (  cv_drho(iglob)*potential_dphi_dx_DG(ibool(i,j,ispec)) & ! \rho'g_x
          !                         + cv_rho0dv(SPACEDIM,iglob)*nabla_v0(1,SPACEDIM,iglob)) * jacLoc ! \rho_0v'_z\partial_zv_{0,x}
          !d0cntrb_rho0dv(SPACEDIM,i,j) = cv_drho(iglob)*potential_dphi_dz_DG(ibool(i,j,ispec)) * jacLoc ! \rho'g_z
          !   Version 3: under HV0 (v_{0,z}=0), potential_dphi_dx_DG=0, and potential_dphi_dz_DG=LNS_g
          d0cntrb_rho0dv(1,i,j) = cv_rho0dv(SPACEDIM,iglob)*nabla_v0(1,SPACEDIM,iglob) * jacLoc ! \rho_0v'_z\partial_zv_{0,x}
          d0cntrb_rho0dv(SPACEDIM,i,j) = cv_drho(iglob)*LNS_g(iglob) * jacLoc ! \rho'g_z
          !   Version 4: under HV0 (v_{0,z}=0), SM (d_xv_0=0), potential_dphi_dx_DG=0, and potential_dphi_dz_DG=LNS_g
          !d0cntrb_rho0dv(1,i,j) = ZEROcr
          d0cntrb_rho0dv(SPACEDIM,i,j) = cv_drho(iglob)*LNS_g(iglob) * jacLoc ! \rho'g_z
          ! > Energy.
          !d0cntrb_dE(i,j)       = (  potential_dphi_dx_DG(ibool(i,j,ispec))*in_dm(1,iglob) &
          !                         + potential_dphi_dz_DG(ibool(i,j,ispec))*in_dm(SPACEDIM,iglob)) * jacLoc
          ! Recall potential_dphi_dx_DG=0 and potential_dphi_dz_DG=LNS_g
          d0cntrb_dE(i,j)       = LNS_g(iglob)*in_dm(SPACEDIM,iglob) * jacLoc
        enddo
      enddo
      
      ! Assemble the contributions previously computed, and add gravity's contribution.
      ! The integration by quadrature on the GLL points leads to three sums. See in particular Komatitsch (Méthodes spectrales et éléments spectraux pour l'équation de l'élastodynamique 2D et 3D en milieu hétérogène), Annexe 3.A.
      do j = 1, NGLLZ
        do i = 1, NGLLX
          iglob = ibool_DG(i,j,ispec)
          do k = 1, NGLLX
            outrhs_drho(iglob) =   outrhs_drho(iglob) &
                                + (  cntrb_drho(k,j,1) * real(hprimewgll_xx(k,i), kind=CUSTOM_REAL) &
                                   + cntrb_drho(i,k,2) * real(hprimewgll_zz(k,j), kind=CUSTOM_REAL))
            do SPCDM=1,SPACEDIM
              outrhs_rho0dv(SPCDM,iglob) =   outrhs_rho0dv(SPCDM,iglob) &
                                          + (  cntrb_rho0dv(SPCDM,k,j,1) * real(hprimewgll_xx(k,i), kind=CUSTOM_REAL) &
                                             + cntrb_rho0dv(SPCDM,i,k,2) * real(hprimewgll_zz(k,j), kind=CUSTOM_REAL))
            enddo
            !outrhs_rho0dv(SPACEDIM,iglob) =   outrhs_rho0dv(SPACEDIM,iglob) &
            !                      + (  cntrb_rho0dv(SPACEDIM,k,j,1) * real(hprimewgll_xx(k,i), kind=CUSTOM_REAL) &
            !                         + cntrb_rho0dv(SPACEDIM,i,k,2) * real(hprimewgll_zz(k,j), kind=CUSTOM_REAL))
            outrhs_dE(iglob) =   outrhs_dE(iglob) &
                              + (  cntrb_dE(k,j,1) * real(hprimewgll_xx(k,i), kind=CUSTOM_REAL) &
                                 + cntrb_dE(i,k,2) * real(hprimewgll_zz(k,j), kind=CUSTOM_REAL))
          enddo ! Enddo on k.
          
          wzl = real(wzgll(j), kind=CUSTOM_REAL)
          wxl = real(wxgll(i), kind=CUSTOM_REAL)
          
          ! Add zero-th order terms.
          !outrhs_drho(iglob)     = outrhs_drho(iglob)     + d0cntrb_drho(i,j) * wxl * wzl
          do SPCDM=1,SPACEDIM
            outrhs_rho0dv(SPCDM,iglob) = outrhs_rho0dv(SPCDM,iglob) + d0cntrb_rho0dv(SPCDM,i,j) * wxl * wzl
          enddo
          outrhs_dE(iglob) = outrhs_dE(iglob) + d0cntrb_dE(i,j) * wxl * wzl
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
          
          ! Step 1: prepare the normals' parameters (n_out(1), n_out(SPACEDIM), weight, etc.).
          ! Interior point
          iglobM = ibool_DG(i,j,ispec)
          
          if(LNS_viscous) then ! Check if viscosity exists whatsoever.
            ! Activate/deactivate, for this particular point (iglob), computation of quantities only needed when viscosity is present.
            if(     LNS_mu(iglob) > 0. &
               .OR. LNS_eta(iglob) > 0. &
               .OR. LNS_kappa(iglob) > 0.) then
              viscousComputation=.true.
            else
              viscousComputation=.false.
            endif
          else
            viscousComputation=.false. ! If viscosity is globally disabled, deactivate it for this element.
          endif
          
          ! TEST WITH IFACE FORMULATION
          n_out(1)     = nx_iface(iface, ispec)
          n_out(SPACEDIM)     = nz_iface(iface, ispec)
          weight = weight_iface(iface1,iface, ispec)
          neighbor = -1
          if(neighbor_DG_iface(iface1, iface, ispec, 3) > -1) then
            iface1_neighbor = neighbor_DG_iface(iface1, iface, ispec, 1)
            iface_neighbor  = neighbor_DG_iface(iface1, iface, ispec, 2)
            ispec_neighbor = neighbor_DG_iface(iface1, iface, ispec, 3)
            neighbor(1) = link_iface_ijispec(iface1_neighbor, iface_neighbor, ispec_neighbor,1)
            neighbor(2) = link_iface_ijispec(iface1_neighbor, iface_neighbor, ispec_neighbor,2)
            neighbor(3) = ispec_neighbor
          endif
          
          ! Step 2: knowing the normals' parameters, compute now the fluxes.
          !drho_P     = ZERO
          !rho0dv_P   = ZERO
          !dE_P       = ZERO
          !veloc_x_DG_P = ZERO
          !veloc_z_DG_P = ZERO
          !in_dp_P       = ZERO
          
          iglobP = 1
          if(neighbor(1) > -1) then
            iglobP = ibool_DG(neighbor(1), neighbor(2), neighbor(3))
          endif
          
          exact_interface_flux = .false. ! Reset this variable to .false.: by default, the fluxes have to be computed (jump!=0). In some specific cases (assigned during the call to LNS_get_interfaces_unknowns), the flux can be exact (jump==0).
          !call compute_interface_unknowns(i,j,ispec, drho_P, rho0dv_P(1), &
          !        rho0dv_P(2), dE_P, veloc_x_DG_P, veloc_z_DG_P, in_dp_P, T_P, &
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
          
          call LNS_get_interfaces_unknowns(i, j, ispec, iface1, iface, neighbor, currentTime, & ! Point identifier (input).
                  cv_drho(iglobM), cv_rho0dv(:,iglobM), cv_dE(iglobM), & ! Input constitutive variables, "M" side.
                  cv_drho(iglobP), cv_rho0dv(:,iglobP), cv_dE(iglobP), & ! Input constitutive variables, "P" side.
                  in_dp(iglobM), & ! Input other variable, "M" side.
                  !V_DG(:,:,iglobM), T_DG(:,iglobM), & ! Input derivatives, "M" side. MIGHT NEED.
                  !V_DG(:,:,iglobP), T_DG(:,iglobP), & ! Input derivatives, "M" side. MIGHT NEED.
                  n_out, & ! Normal vector (input).
                  exact_interface_flux, & ! Switch to disable jump in some cases (output).
                  drho_P, rho0dv_P, dE_P, & ! Output constitutive variables.
                  !Tx_DG_P, Tz_DG_P, Vxx_DG_P, Vzz_DG_P, Vzx_DG_P, Vxz_DG_P,& ! Output derivatives. MIGHT NEED.
                  dm_P, dp_P, dv_P,& ! Output other variables.
                  viscousComputation, nabla_dT_P, sigma_dv_P, & ! Output other variables: viscous.
                  .false., LNS_dummy_1d(1), LNS_dummy_2d(:,1)) ! Use dummies and set the switch to false not to compute unecessary quantities.
          
          ! Recover an approximate local maximum linearized acoustic wave speed. See for example Hesthaven (doi.org/10.1007/9780387720678), page 208.
          lambda = ZERO
          jump   = ZERO
          !veloc_n_M = abs(LNS_dv(1,iglobM)*n_out(1) + LNS_dv(SPACEDIM,iglobM)*n_out(SPACEDIM))
          !veloc_n_P = abs(veloc_x_DG_P*n_out(1) + veloc_z_DG_P*n_out(SPACEDIM))
          ! Save some allocations.
          !lambda = max( veloc_n_M + sqrt(abs(gammaext_DG(iglobM)*in_dp(iglobM)/cv_drho(iglobM))), &
          !              veloc_n_P + sqrt(abs(gamma_P*in_dp_P/drho_P)) )
          !lambda = max(  abs(  (LNS_v0(1,iglobM)+LNS_dv(1,iglobM))*n_out(1) & ! [v_-\cdot n]_1
          !                   + (LNS_v0(SPACEDIM,iglobM)+LNS_dv(SPACEDIM,iglobM))*n_out(SPACEDIM) & ! [v_-\cdot n]_d
          !                   + sqrt(abs(  gammaext_DG(iglobM) &
          !                              * (LNS_p0(iglobM)+in_dp(iglobM)) &
          !                              / (LNS_rho0(iglobM)+cv_drho(iglobM)))) & ! Local sound speed.
          !                  ) & ! Local sound speed, side "M".
          !             , abs(  (LNS_v0(1,iglobP)+dv_P(1))*n_out(1) & ! [v_+\cdot n]_1
          !                   + (LNS_v0(SPACEDIM,iglobP)+dv_P(SPACEDIM))*n_out(SPACEDIM) & ! [v_+\cdot n]_d
          !                   + sqrt(abs(  gammaext_DG(iglobP) &
          !                              * (LNS_p0(iglobP)+dp_P) &
          !                              / (LNS_rho0(iglobP)+drho_P))) & ! Local sound speed.
          !                  ) & ! Local sound speed, side "P".
          !             )
          lambda = max(  abs(  dot_product(n_out, LNS_v0(:,iglobM)+LNS_dv(:,iglobM)) & ! v_-\cdot n
                             + sqrt(abs(  gammaext_DG(iglobM) &
                                        * (LNS_p0(iglobM)+in_dp(iglobM)) &
                                        / (LNS_rho0(iglobM)+cv_drho(iglobM)))) & ! Local sound speed.
                            ) & ! Local sound speed, side "M".
                       , abs(  dot_product(n_out, LNS_v0(:,iglobP)+LNS_dv(:,iglobP)) & ! v_+\cdot n
                             + sqrt(abs(  gammaext_DG(iglobP) &
                                        * (LNS_p0(iglobP)+dp_P) &
                                        / (LNS_rho0(iglobP)+drho_P))) & ! Local sound speed.
                            ) & ! Local sound speed, side "P".
                       )
          !if(      coord(2,ibool_before_perio(i,j,ispec))<1. & ! DEBUG
          !   .and. coord(2,ibool_before_perio(i,j,ispec))>=0. & ! DEBUG
          !   .and. abs(coord(1,ibool_before_perio(i,j,ispec)))<2.) then ! DEBUG
          !  write(*,*) currentTime, coord(:,ibool_before_perio(i,j,ispec)), & ! DEBUG
          !             !LNS_p0(iglobM)+inp_dp_M, LNS_p0(iglobM)+out_dp_P ! DEBUG
          !             !out_rho0dv_P ! DEBUG
          !             !out_dv_P ! DEBUG
          !             !trans_boundary ! DEBUG
          !             lambda
          !endif ! DEBUG
          
          ! Viscous stress tensor's contributions (already under the form of the average mean flux).
          !dux_dx = ZERO
          !dux_dz = ZERO
          !duz_dx = ZERO
          !duz_dz = ZERO
          !dT_dx = ZERO
          !dT_dz = ZERO
          !if(LNS_mu(iglob) > 0 .OR. &
          !   LNS_eta(iglob) > 0 .OR. &
          !   LNS_kappa(iglob)  > 0) then
          !  dux_dx = HALFcr*(V_DG(1, 1, iglobM) + Vxx_DG_P)
          !  dux_dz = HALFcr*(V_DG(1, 2, iglobM) + Vxz_DG_P)
          !  duz_dx = HALFcr*(V_DG(2, 1, iglobM) + Vzx_DG_P)
          !  duz_dz = HALFcr*(V_DG(2, 2, iglobM) + Vzz_DG_P)
          !  dT_dx = HALFcr*(T_DG(1, iglobM) + Tx_DG_P)
          !  dT_dz = HALFcr*(T_DG(2, iglobM) + Tz_DG_P)
          !endif
          
          ! Mass conservation (fully inviscid).
          flux_n = (in_dm(1,iglobM)+dm_P(1))*n_out(1) + (in_dm(SPACEDIM,iglobM)+dm_P(SPACEDIM))*n_out(SPACEDIM) ! Dot product: \Sigma\cdot n.
          !flux_n = (tmp_unknown_x_M+tmp_unknown_x_P)*n_out(1) + (tmp_unknown_z_M+tmp_unknown_z_P)*n_out(SPACEDIM) ! Recall: n_out(SPACEDIM)=n_z.
          if(exact_interface_flux) then
            jump = ZERO
          else
            jump = cv_drho(iglobM) - drho_P
          endif
          outrhs_drho(iglobM) = outrhs_drho(iglobM) - weight*(flux_n + lambda*jump)*HALFcr ! Add flux' contribution.
          ! x-Momentum inviscid contributions.
          tmp_unknown_x =   cv_rho0dv(1,iglobM)*LNS_v0(1,iglobM) + in_dp(iglobM) & ! "M" side.
                          + rho0dv_P(1)*LNS_v0(1,iglobP) + dp_P ! "P" side.
          tmp_unknown_z =   cv_rho0dv(SPACEDIM,iglobM)*LNS_v0(1,iglobM) & ! "M" side.
                          + rho0dv_P(SPACEDIM)*LNS_v0(1,iglobP) ! "P" side.
          flux_n = tmp_unknown_x*n_out(1) + tmp_unknown_z*n_out(SPACEDIM) ! Dot product: \Sigma\cdot n.
          !tmp_unknown_x_M = cv_rho0dv(1,iglobM)*LNS_v0(1,iglobM) + in_dp(iglobM)
          !tmp_unknown_z_M = cv_rho0dv(SPACEDIM,iglobM)*LNS_v0(1,iglobM)
          !tmp_unknown_x_P = rho0dv_P(1)*LNS_v0(1,iglobP) + dp_P
          !tmp_unknown_z_P = rho0dv_P(SPACEDIM)*LNS_v0(1,iglobP)
          !flux_n = (tmp_unknown_x_M+tmp_unknown_x_P)*n_out(1) + (tmp_unknown_z_M+tmp_unknown_z_P)*n_out(SPACEDIM) ! Dot product: \Sigma\cdot n.
          if(exact_interface_flux) then
            jump = ZERO
          else
            jump = cv_rho0dv(1,iglobM) - rho0dv_P(1)
          endif
          outrhs_rho0dv(1,iglobM) = outrhs_rho0dv(1,iglobM) - weight*(flux_n + lambda*jump)*HALFcr ! Add flux' contribution.
          ! z-Momentum inviscid contributions.
          tmp_unknown_x =   cv_rho0dv(1,iglobM)*LNS_v0(SPACEDIM,iglobM) & ! "M" side.
                          + rho0dv_P(1)*LNS_v0(SPACEDIM,iglobP) ! "P" side.
          tmp_unknown_z =   cv_rho0dv(SPACEDIM,iglobM)*LNS_v0(SPACEDIM,iglobM) + in_dp(iglobM) & ! "M" side.
                          + rho0dv_P(SPACEDIM)*LNS_v0(SPACEDIM,iglobP) + dp_P ! "P" side.
          flux_n = tmp_unknown_x*n_out(1) + tmp_unknown_z*n_out(SPACEDIM) ! Dot product: \Sigma\cdot n.
          !tmp_unknown_x_M = cv_rho0dv(1,iglobM)*LNS_v0(SPACEDIM,iglobM)
          !tmp_unknown_z_M = cv_rho0dv(SPACEDIM,iglobM)*LNS_v0(SPACEDIM,iglobM) + in_dp(iglobM)
          !tmp_unknown_x_P = rho0dv_P(1)*LNS_v0(SPACEDIM,iglobP)
          !tmp_unknown_z_P = rho0dv_P(SPACEDIM)*LNS_v0(SPACEDIM,iglobP) + dp_P
          !flux_n = (tmp_unknown_x_M+tmp_unknown_x_P)*n_out(1) + (tmp_unknown_z_M+tmp_unknown_z_P)*n_out(SPACEDIM) ! Dot product: \Sigma\cdot n.
          if(exact_interface_flux) then
            jump = ZERO
          else
            jump = cv_rho0dv(SPACEDIM,iglobM) - rho0dv_P(2)
          endif
          outrhs_rho0dv(SPACEDIM,iglobM) = outrhs_rho0dv(SPACEDIM,iglobM) - weight*(flux_n + lambda*jump)*HALFcr ! Add flux' contribution.
          ! Energy inviscid contributions.
          tmp_unknown_x =   LNS_dv(1,iglobM)*(LNS_E0(iglobM) + LNS_p0(iglobM)) & ! "M" side, part.
                          + LNS_v0(1,iglobM)*(cv_dE(iglobM) + in_dp(iglobM)) & ! "M" side, part.
                          + dv_P(1)*(LNS_E0(iglobP) + LNS_p0(iglobP)) & ! "P" side, part.
                          + LNS_v0(1,iglobP)*(dE_P + dp_P) ! "P" side, part.
          tmp_unknown_z =   LNS_dv(SPACEDIM,iglobM)*(LNS_E0(iglobM) + LNS_p0(iglobM)) & ! "M" side, part.
                          + LNS_v0(SPACEDIM,iglobM)*(cv_dE(iglobM) + in_dp(iglobM)) & ! "M" side, part.
                          + dv_P(SPACEDIM)*(LNS_E0(iglobP) + LNS_p0(iglobP)) & ! "P" side, part.
                          + LNS_v0(SPACEDIM,iglobP)*(dE_P + dp_P) ! "P" side, part.
          flux_n = tmp_unknown_x*n_out(1) + tmp_unknown_z*n_out(SPACEDIM) ! Dot product: \Sigma\cdot n.
          !tmp_unknown_x_M = LNS_dv(1,iglobM)*(cv_dE(iglobM) + in_dp(iglobM))
          !tmp_unknown_x_P = veloc_x_DG_P*(dE_P + in_dp_P)
          !tmp_unknown_z_M = LNS_dv(SPACEDIM,iglobM)*(cv_dE(iglobM) + in_dp(iglobM))
          !tmp_unknown_z_P = veloc_z_DG_P*(dE_P + in_dp_P)
          !flux_n = (tmp_unknown_x_M+tmp_unknown_x_P)*n_out(1) + (tmp_unknown_z_M+tmp_unknown_z_P)*n_out(SPACEDIM) ! Dot product: \Sigma\cdot n.
          if(exact_interface_flux) then
            jump = ZERO
          else
            jump = cv_dE(iglobM) - dE_P
          endif
          outrhs_dE(iglobM) = outrhs_dE(iglobM) - weight*(flux_n + lambda*jump)*HALFcr ! Add flux' contribution.
          
          ! Viscous stress tensor's contributions.
          if(viscousComputation) then
            ! Note: since we consider here viscous terms, that is purely diffusive phenomena, the acoustic wave speed makes no sense, and thus the jump in the Lax-Friedrich approximation does not appear in the formulations.
            ! x-Momentum viscous contributions. The vector [tmp_unknown_x, tmp_unknown_z] represents the mean average flux at the boundary of the x-momentum (based on the 1st line of \Sigma_v).
            tmp_unknown_x = -in_sigma_dv(1,iglob) -sigma_dv_P(1) ! Recall, this corresponds to \Sigma_v(v')_{1,1}.
            tmp_unknown_z = -in_sigma_dv(2,iglob) -sigma_dv_P(2) ! Recall, this corresponds to \Sigma_v(v')_{1,2} (and \Sigma_v(v')_{2,1}).
            !tmp_unknown_x = LNS_mu(iglob)*TWOcr*dux_dx + (LNS_eta(iglob) - (TWOcr/3.)*LNS_mu(iglob))*(dux_dx + duz_dz) 
            !tmp_unknown_z = LNS_mu(iglob)*( dux_dz + duz_dx )
            outrhs_rho0dv(1,iglobM) = outrhs_rho0dv(1,iglobM) - weight*(tmp_unknown_x*n_out(1)+tmp_unknown_z*n_out(SPACEDIM))*HALFcr ! Add flux' contribution, with dot product \Sigma\cdot n.
            ! z-Momentum viscous contributions. The vector [tmp_unknown_x, tmp_unknown_z] represents the mean average flux at the boundary of the z-momentum (based on the 2nd line of \Sigma_v).
            tmp_unknown_x = -in_sigma_dv(2,iglob) -sigma_dv_P(2) ! Recall, this corresponds to \Sigma_v(v')_{2,1} (and \Sigma_v(v')_{1,2}).
            tmp_unknown_z = -in_sigma_dv(3,iglob) -sigma_dv_P(3) ! Recall, this corresponds to \Sigma_v(v')_{2,2}.
            !tmp_unknown_x = LNS_mu(iglob)*( dux_dz + duz_dx )
            !tmp_unknown_z = LNS_mu(iglob)*TWOcr*duz_dz + (LNS_eta(iglob) - (TWOcr/3.)*LNS_mu(iglob))*(dux_dx + duz_dz)
            outrhs_rho0dv(SPACEDIM,iglobM) =   outrhs_rho0dv(SPACEDIM,iglobM) &
                                             - weight*(tmp_unknown_x*n_out(1)+tmp_unknown_z*n_out(SPACEDIM))*HALFcr ! Add flux' contribution, with dot product \Sigma\cdot n.
            ! Energy viscous contributions.
            tmp_unknown_x = - (  LNS_dv(1,iglobM)*sigma_v_0(1,iglobM) & ! "M" side, part.
                               + LNS_v0(1,iglobM)*in_sigma_dv(1,iglobM) & ! "M" side, part.
                               + LNS_dv(SPACEDIM,iglobM)*sigma_v_0(2,iglobM) & ! "M" side, part.
                               + LNS_v0(SPACEDIM,iglobM)*in_sigma_dv(2,iglobM)) & ! "M" side, part.
                            + LNS_kappa(iglobM)*in_nabla_dT(1,iglobM) & ! "M" side, part.
                            - (  dv_P(1)*sigma_v_0(1,iglobP) & ! "P" side, part.
                               + LNS_v0(1,iglobM)*sigma_dv_P(1) & ! "P" side, part.
                               + dv_P(SPACEDIM)*sigma_v_0(2,iglobP) & ! "P" side, part.
                               + LNS_v0(SPACEDIM,iglobP)*sigma_dv_P(2)) & ! "P" side, part.
                            + LNS_kappa(iglobP)*nabla_dT_P(1) ! "P" side, part.
            tmp_unknown_z = - (  LNS_dv(SPACEDIM,iglobM)*sigma_v_0(3,iglobM) & ! "M" side, part.
                               + LNS_v0(SPACEDIM,iglobM)*in_sigma_dv(3,iglobM) & ! "M" side, part.
                               + LNS_dv(1,iglobM)*sigma_v_0(2,iglobM) & ! "M" side, part.
                               + LNS_v0(1,iglobM)*in_sigma_dv(2,iglobM)) & ! "M" side, part.
                            + LNS_kappa(iglobM)*in_nabla_dT(2,iglobM) & ! "M" side, part.
                            - (  dv_P(SPACEDIM)*sigma_v_0(3,iglobP) & ! "P" side, part.
                               + LNS_v0(SPACEDIM,iglobP)*sigma_dv_P(3) & ! "P" side, part.
                               + dv_P(1)*sigma_v_0(2,iglobP) & ! "P" side, part.
                               + LNS_v0(1,iglobP)*sigma_dv_P(2)) & ! "P" side, part.
                            + LNS_kappa(iglobP)*nabla_dT_P(SPACEDIM) ! "P" side, part.
            outrhs_dE(iglobM) = outrhs_dE(iglobM) + weight*(tmp_unknown_x*n_out(1) + tmp_unknown_z*n_out(SPACEDIM))*HALFcr
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
! Variables ending in "_P" (for "plus") are output-intended and are the sought values, or exterior values.
! Variables ending in "_iM" (for "minus") are input-intended and should correspond to the interior values.
! Variables ending in "_iP" are input-intended, and correspond to the exterior values, if known. Remark that those variables are only used if neighbor(3) != -1.
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
             inp_drho_M, inp_rho0dv_M, inp_dE_M, & ! Input constitutive variables, "M" side.
             inp_drho_P, inp_rho0dv_P, inp_dE_P, & ! Input constitutive variables, "P" side.
             inp_dp_M, & ! Input other variable, "M" side.
             !V_DG_iM, T_DG_iM, V_DG_iP, T_DG_iP, & ! Input derivatives. MIGHT NEED.
             n_out, & ! Normal vector (input).
             exact_interface_flux, & ! Switch to disable jump in some cases (output).
             out_drho_P, out_rho0dv_P, out_dE_P, & ! Output constitutive variables.
             !Tx_DG_P, Tz_DG_P, Vxx_DG_P, Vzz_DG_P, Vzx_DG_P, Vxz_DG_P,& ! Output derivatives. MIGHT NEED.
             out_dm_P, out_dp_P, out_dv_P, & ! Output other variables.
             swCompVisc, out_nabla_dT_P, out_sigma_dv_P, & ! Output other variables: viscous.
             swCompT, out_T_P, out_v_P) ! Output other variables.
  
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
  real(kind=CUSTOM_REAL), intent(in) :: inp_drho_M, inp_dE_M, inp_drho_P, inp_dE_P ! Input constitutive variables.
  real(kind=CUSTOM_REAL), dimension(SPACEDIM), intent(in) :: inp_rho0dv_M, inp_rho0dv_P ! Input constitutive variables.
  real(kind=CUSTOM_REAL), intent(in) :: inp_dp_M ! Input other variables.
  !real(kind=CUSTOM_REAL), dimension(SPACEDIM), intent(in) :: T_DG_iP, T_DG_iM ! Input derivatives. MIGHT NEED.
  !real(kind=CUSTOM_REAL), dimension(SPACEDIM, SPACEDIM), intent(in) :: V_DG_iP, V_DG_iM ! Input derivatives. MIGHT NEED.
  real(kind=CUSTOM_REAL), dimension(SPACEDIM), intent(in) :: n_out
  logical, intent(out) :: exact_interface_flux ! Output switch.
  real(kind=CUSTOM_REAL), intent(out) :: out_drho_P, out_dE_P ! Output constitutive variables.
  !real(kind=CUSTOM_REAL), intent(out) :: Tx_DG_P, Tz_DG_P, Vxx_DG_P, Vzz_DG_P, Vzx_DG_P, Vxz_DG_P ! Output derivatives. MIGHT NEED.
  real(kind=CUSTOM_REAL), dimension(SPACEDIM), intent(out) :: out_rho0dv_P ! Output constitutive variable.
  real(kind=CUSTOM_REAL), dimension(SPACEDIM), intent(out) :: out_dv_P, out_dm_P, out_nabla_dT_P ! Output other variables.
  real(kind=CUSTOM_REAL), intent(out) :: out_dp_P ! Output other variables.
  real(kind=CUSTOM_REAL), dimension(NVALSIGMA), intent(out) :: out_sigma_dv_P ! Output other variables.
  real(kind=CUSTOM_REAL), intent(out) :: out_T_P ! In compute_gradient_TFSF (desintegration method), we need temperature on the other side of the boundary in order to compute the flux.
  real(kind=CUSTOM_REAL), dimension(SPACEDIM), intent(out) :: out_v_P ! In compute_gradient_TFSF (desintegration method), we need velocity on the other side of the boundary in order to compute the flux.
  logical, intent(in) :: swCompVisc, swCompT ! Do not unnecessarily compute some quantities.
  
  ! Local.
  integer :: SPCDM
  real(kind=CUSTOM_REAL), dimension(SPACEDIM) :: tang ! Tangential vector.
  real(kind=CUSTOM_REAL) :: normal_v, tangential_v!, &
                            !veloc_x, veloc_z!, gamma_P!, &
                            !dv_M(1), dv_M(SPACEDIM), inp_dp_M, gamma_P, e1_DG_P
  real(kind=CUSTOM_REAL), dimension(SPACEDIM) :: veloc_P
  real(kind=CUSTOM_REAL), parameter :: ZEROcr = 0._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: ONEcr  = 1._CUSTOM_REAL
  !real(kind=CUSTOM_REAL), parameter :: TWO  = 2._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: HALFcr = 0.5_CUSTOM_REAL
  integer :: iglobM, i_el, j_el, ispec_el, i_ac, j_ac, ispec_ac, iglob, iglobP, ipoin, num_interface
  real(kind=CUSTOM_REAL), dimension(SPACEDIM, SPACEDIM) :: trans_boundary
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
  out_T_P=ZEROcr
  out_v_P=ZEROcr
  
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
  ! Set out_T_P.
  ! Set out_v_P.
  
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
                                            .false., out_v_P, &
                                            .false., LNS_dummy_1d(1), &
                                            .false., LNS_dummy_1d(1)) ! Get needed background parameters. Use dummies for values we're not interested in.
        ! Set out_v_P: free slip and normal velocity continuity.
        !call boundary_condition_DG(i, j, ispec, timelocal, out_drho_P, out_rho0dv_P(1), out_rho0dv_P(SPACEDIM), out_dE_P, &
        !                           out_dv_P(1), out_dv_P(SPACEDIM), out_dp_P, e1_DG_P) ! Warning, expression of out_dp_P might not be exact.
        veloc_P = veloc_vector_acoustic_DG_coupling(iglob, :)
        !veloc_x = veloc_vector_acoustic_DG_coupling(iglob, 1)
        !veloc_z = veloc_vector_acoustic_DG_coupling(iglob, 2)
        call build_trans_boundary(n_out, tang, trans_boundary)
        ! Tangential vector
        ! Since only bottom topography n_out(SPACEDIM) > 0
        !tang(1) = -n_out(SPACEDIM) ! Recall: n_out(SPACEDIM)=n_z.
        !tang(SPACEDIM) = n_out(1) ! Recall: n_out(1)=n_x.
        ! Normal velocity of the solid perturbation and tangential velocity of the fluid flow
        !normal_v     = veloc_P(1)*n_out(1) + veloc_P(1)*n_out(SPACEDIM) ! Recall: n_out(SPACEDIM)=n_z.
        !tangential_v = veloc_P(1)*tang(1) + veloc_P(1)*tang(SPACEDIM)
        normal_v     = DOT_PRODUCT(n_out, veloc_P)
        tangential_v = DOT_PRODUCT(tang, veloc_P)
        
        do SPCDM=1,SPACEDIM
          out_v_P(SPCDM) = trans_boundary(SPCDM, 1)*normal_v + trans_boundary(SPCDM, 2)*tangential_v!veloc_elastic(1,iglob)
          !out_v_P(1) = trans_boundary(1, 1)*normal_v + trans_boundary(1, 2)*tangential_v!veloc_elastic(1,iglob)
          !out_v_P(SPACEDIM) = trans_boundary(2, 1)*normal_v + trans_boundary(2, 2)*tangential_v
        enddo
        ! Set out_dv_P.
        out_dv_P=out_v_P-LNS_v0(:,iglobM) ! Requesting iglobM might be technically inexact, but on elements' boundaries points should overlap. Plus, iglobP does not exist on outer computational domain boundaries.
        ! Set out_dp_P: traction continuity.
        out_dp_P = LNS_p0(iglobM) - potential_dot_dot_acoustic(iglob) ! Warning, expression of out_dp_P might not be exact.
        ! Set out_dE_P.
        call compute_dE_i(LNS_rho0(iglobM)+out_drho_P, out_v_P, LNS_p0(iglobM)+out_dp_P, out_dE_P, iglobM) ! Warning, expression of out_dp_P might not be exact.
        !out_dE_P       = out_dp_P/(gammaext_DG(iglobM) - 1.) + HALFcr*out_drho_P*( out_dv_P(1)**2 + out_dv_P(SPACEDIM)**2 )
        ! Set out_rho0dv_P.
        do SPCDM=1,SPACEDIM
          out_rho0dv_P(SPCDM) = out_drho_P*out_dv_P(SPCDM) ! Safe version, just in case element-wise mutiplication fails somehow.
        enddo
        !out_rho0dv_P=out_drho_P*out_dv_P
        ! Set out_dm_P: see bottom of routine.
        if(swCompVisc) then
          ! Set out_nabla_dT_P: same as other side.
          out_nabla_dT_P = nabla_dT(:,iglobM)
          ! Set out_sigma_dv_P: same as other side.
          out_sigma_dv_P = sigma_dv(:,iglobM)
        endif
        ! Set out_T_P.
        if(swCompT) then
          !out_dT_P = (inp_dE_M/inp_drho_M - 0.5*(dv_M(1)**2 + dv_M(SPACEDIM)**2))/c_V
          call compute_T_i(LNS_rho0(iglobM)+out_drho_P, out_v_P, LNS_E0(iglobM)+out_dE_P, out_T_P)
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
        out_dv_P=out_rho0dv_P/LNS_rho0(iglobM)
        !out_dv_P(1) = out_rho0dv_P(1)/out_drho_P
        !out_dv_P(SPACEDIM) = out_rho0dv_P(SPACEDIM)/out_drho_P
        
        ! Set out_dp_P.
        call compute_dp_i(LNS_rho0(iglobM)+out_drho_P, LNS_v0(:,iglobM)+out_dv_P, LNS_E0(iglobM)+out_dE_P, out_dp_P, iglobM)
        !out_dp_P = (gamma_P - ONEcr)*( out_dE_P & ! Warning, expression of out_dp_P might not be exact.
        !         - (HALFcr)*out_drho_P*( out_dv_P(1)**2 + out_dv_P(SPACEDIM)**2 ) )
        
        ! Set out_dm_P: see bottom of routine.
        
        if(swCompVisc) then
          ! Set out_nabla_dT_P: get the values from the MPI buffers.
          out_nabla_dT_P = buffer_LNS_nabla_dT(:, ipoin, num_interface)
          
          ! Set out_sigma_dv_P: get the values from the MPI buffers.
          out_sigma_dv_P = buffer_LNS_sigma_dv(:, ipoin, num_interface)
        endif
        
        ! Set out_v_P.
        out_v_P = out_dv_P + LNS_v0(:,iglobM)
        
        ! Set out_T_P.
        if(swCompT) then
          call compute_T_i(LNS_rho0(iglobM)+out_drho_P, out_v_P, LNS_E0(iglobM)+out_dE_P, out_T_P)
        endif
        !out_dT_P = (out_dE_P/out_drho_P - 0.5*((out_rho0dv_P(1)/out_drho_P)**2 + (out_rho0dv_P(SPACEDIM)/out_drho_P)**2))/c_V
        
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
      
      ! Set out_drho_P: same as other side.
      out_drho_P = inp_drho_M
      
      ! Set out_v_P: get elastic medium velocities.
      !veloc_P = veloc_elastic(:, iglob)
      call build_trans_boundary(n_out, tang, trans_boundary)
      !normal_v     = veloc_elastic(1, iglob)*n_out(1)  + veloc_elastic(SPACEDIM, iglob)*n_out(SPACEDIM)
      normal_v     = DOT_PRODUCT(n_out, veloc_elastic(:,ibool(i_el, j_el, ispec_el)))
      !tangential_v = veloc_x_DG_P*tang(1) + veloc_z_DG_P*tang(SPACEDIM)
      tangential_v = DOT_PRODUCT(tang, LNS_v0(:,iglobM)+LNS_dv(:,iglobM))
      do SPCDM=1,SPACEDIM
        out_v_P(SPCDM) = trans_boundary(SPCDM, 1)*normal_v + trans_boundary(SPCDM, 2)*tangential_v
      enddo
      
      ! Set out_dv_P.
      out_dv_P=out_v_P-LNS_v0(:,iglobM) ! Requesting iglobM might be technically inexact, but on elements' boundaries points should overlap. Plus, iglobP does not exist on outer computational domain boundaries.
      
      ! Set out_dp_P: no stress continuity.
      out_dp_P = inp_dp_M

      ! Set out_dE_P.
      call compute_dE_i(LNS_rho0(iglobM)+out_drho_P, out_v_P, LNS_p0(iglobM)+out_dp_P, out_dE_P, iglobM)

      ! Set out_rho0dv_P.
      !do SPCDM=1,SPACEDIM
      !  out_rho0dv_P(SPCDM) = LNS_rho0(iglobM)*out_dv_P(SPCDM) ! Safe version, just in case element-wise mutiplication fails somehow.
      !enddo
      out_rho0dv_P(:)=LNS_rho0(iglobM)*out_dv_P(:)
      
      ! Set out_dm_P: see bottom of routine.
        
      if(swCompVisc) then
        ! Set out_nabla_dT_P: same as other side.
        out_nabla_dT_P = nabla_dT(:,iglobM)
        
        ! Set out_sigma_dv_P: same as other side.
        out_sigma_dv_P = sigma_dv(:,iglobM)
      endif
      
      ! Set out_T_P.
      if(swCompT) then
        call compute_T_i(LNS_rho0(iglobM)+out_drho_P, out_v_P, LNS_E0(iglobM)+out_dE_P, out_T_P)
      endif
      
      !if(      coord(2,ibool_before_perio(i,j,ispec))<1. & ! DEBUG
      !   .and. coord(2,ibool_before_perio(i,j,ispec))>=0. & ! DEBUG
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
      
      ! Set exact_interface_flux.
      exact_interface_flux = .false.
      
      ! Set out_drho_P and out_dp_P.
      call background_physical_parameters(i, j, ispec, timelocal, out_drho_P, &
                                          .true., out_v_P, & ! We get the far-field v field here.
                                          .false., LNS_dummy_1d(1), &
                                          .true., out_dp_P) ! Get needed background parameters. Use dummies for values we're not interested in.
      ! Warning: out_drho_P contains rho=(rho0 + drho) here. Correct this.
      out_drho_P = out_drho_P - LNS_rho0(iglobM)
      ! Warning: out_dp_P contains p=(p0 + dp) here. Correct this.
      out_dp_P = out_dp_P - LNS_p0(iglobM)
      
      ! Set out_v_P.
      ! Convert the velocity components from mesh coordinates to normal/tangential coordinates.
      ! This is more practical to set the boundary conditions.
      ! The setting of the boundary conditions is thus made here.
      call build_trans_boundary(n_out, tang, trans_boundary)
      !tang(1) = -n_out(SPACEDIM) ! Recall: n_out(SPACEDIM)=n_z.
      !tang(SPACEDIM) =  n_out(1) ! Recall: n_out(1)=n_x.
      !normal_v = (out_v_P(1)*n_out(1) + out_v_P(SPACEDIM)*n_out(SPACEDIM))
      !tangential_v = (out_v_P(1)*tang(1) + out_v_P(SPACEDIM)*tang(SPACEDIM))
      normal_v     = DOT_PRODUCT(n_out, out_v_P)
      tangential_v = DOT_PRODUCT(tang, out_v_P)
      ! Notes:
      !   At that point, normal_v and tangential_v are set to their "background" (or "far-field", or "unperturbed") values.
      !   Treatment of boundary conditions based on normal/tangential velocities should be done here. The "free slip" condition and the "normal velocity continuity" conditions can be set here.
      ! Convert (back) the velocity components from normal/tangential coordinates to mesh coordinates.
      do SPCDM=1,SPACEDIM
        out_v_P(SPCDM) = trans_boundary(SPCDM, 1)*normal_v + trans_boundary(SPCDM, 2)*tangential_v!veloc_elastic(1,iglob)
        !out_v_P(1) = trans_boundary(1, 1)*normal_v + trans_boundary(1, 2)*tangential_v
        !out_v_P(SPACEDIM) = trans_boundary(2, 1)*normal_v + trans_boundary(2, 2)*tangential_v
      enddo
      
      ! Set out_dv_P.
      out_dv_P=out_v_P-LNS_v0(:,iglobM) ! Requesting iglobM might be technically inexact, but on elements' boundaries points should overlap. Plus, iglobP does not exist on outer computational domain boundaries.
      
      ! Set out_dE_P.
      call compute_dE_i(LNS_rho0(iglobM)+out_drho_P, out_v_P, LNS_p0(iglobM)+out_dp_P, out_dE_P, iglobM) ! TODO: Warning, expression of out_dp_P might not be exact.
      !out_dE_P = out_dp_P/(gammaext_DG(iglobM) - ONEcr) &
      !         + out_drho_P*HALFcr*( out_dv_P(1)**2 + out_dv_P(SPACEDIM)**2 )
      
      ! Set out_rho0dv_P.
      do SPCDM=1,SPACEDIM
        out_rho0dv_P(SPCDM) = out_drho_P*out_dv_P(SPCDM) ! Safe version, just in case element-wise mutiplication fails somehow.
      enddo
      !out_rho0dv_P=out_drho_P*out_dv_P
      
      ! Set out_dm_P: see bottom of routine.
        
      if(swCompVisc) then
        ! Set out_nabla_dT_P: same as other side.
        out_nabla_dT_P = nabla_dT(:,iglobM)
        
        ! Set out_sigma_dv_P: same as other side.
        out_sigma_dv_P = sigma_dv(:,iglobM)
      endif
      
      ! Set out_T_P.
      if(swCompT) then
        call compute_T_i(LNS_rho0(iglobM)+out_drho_P, out_v_P, LNS_E0(iglobM)+out_dE_P, out_T_P)
      endif
      
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
    ! Tangential vector, since only bottom topography n_out(SPACEDIM) > 0.
    !tang(1) = -n_out(SPACEDIM) ! Recall: n_out(SPACEDIM)=n_z.
    !tang(SPACEDIM) =  n_out(1) ! Recall: n_out(1)=n_x.
    ! Normal velocity of the solid perturbation and tangential velocity of the fluid flow
    !normal_v     = veloc_P(1)*n_out(1) + veloc_P(1)*n_out(SPACEDIM) ! Recall: n_out(SPACEDIM)=n_z.
    !tangential_v = veloc_P(1)*tang(1) + veloc_P(1)*tang(SPACEDIM)
    normal_v     = DOT_PRODUCT(n_out, veloc_P)
    tangential_v = DOT_PRODUCT(tang, veloc_P)
    ! Set the matrix for the transformation from normal/tangential coordinates to mesh coordinates.
    !trans_boundary(1, 1) =  tang(SPACEDIM)
    !trans_boundary(1, 2) = -n_out(SPACEDIM) ! Recall: n_out(SPACEDIM)=n_z.
    !trans_boundary(2, 1) = -tang(1)
    !trans_boundary(2, 2) =  n_out(1) ! Recall: n_out(1)=n_x.
    !trans_boundary = trans_boundary/(n_out(1) * tang(SPACEDIM) - tang(1) * n_out(SPACEDIM))
    ! Set out_v_P: free slip and normal velocity continuity.
    do SPCDM=1,SPACEDIM
      out_v_P(SPCDM) = trans_boundary(SPCDM, 1)*normal_v + trans_boundary(SPCDM, 2)*tangential_v!veloc_elastic(1,iglob)
      !out_v_P(1) = trans_boundary(1, 1)*normal_v + trans_boundary(1, 2)*tangential_v!veloc_elastic(1,iglob)
      !out_v_P(SPACEDIM) = trans_boundary(2, 1)*normal_v + trans_boundary(2, 2)*tangential_v
    enddo
    ! Set out_dv_P.
    out_dv_P=out_v_P-LNS_v0(:,iglobM) ! Requesting iglobM might be technically inexact, but on elements' boundaries points should overlap. Plus, iglobP does not exist on outer computational domain boundaries.
    ! Set out_dp_P: traction continuity.
    out_dp_P = LNS_p0(iglobM) - potential_dot_dot_acoustic(iglob) ! Warning, expression of out_dp_P might not be exact.
    ! Set out_dE_P.
    call compute_dE_i(LNS_rho0(iglobM)+out_drho_P, LNS_v0(:,iglobM)+out_dv_P, LNS_p0(iglobM)+out_dp_P, out_dE_P, iglobM) ! Warning, expression of out_dp_P might not be exact.
    !out_dE_P = out_dp_P/(gammaext_DG(iglobM) - 1.) + HALFcr*out_drho_P*( out_dv_P(1)**2 + out_dv_P(SPACEDIM)**2 ) ! Warning, expression of out_dp_P might not be exact.
    ! Set out_rho0dv_P.
    do SPCDM=1,SPACEDIM
      out_rho0dv_P(SPCDM) = out_drho_P*out_dv_P(SPCDM) ! Safe version, just in case element-wise mutiplication fails somehow.
    enddo
    !out_rho0dv_P=out_drho_P*out_dv_P
    ! Set gradients, and sigma
    !out_nabla_dT_P
    !call compute_gradient_TFSF(LNS_dv, LNS_dT, .true., .true., switch_gradient, nabla_dv, nabla_dT)
    !out_dT_P = (inp_dE_M/inp_drho_M - 0.5*(dv_M(1)**2 + dv_M(SPACEDIM)**2))/c_V
    ! Set out_dm_P: see bottom of routine.
    if(swCompVisc) then
      ! Set out_nabla_dT_P: same as other side.
      out_nabla_dT_P = nabla_dT(:,iglobM)
      ! Set out_sigma_dv_P: same as other side.
      out_sigma_dv_P = sigma_dv(:,iglobM)
    endif
    ! Set out_T_P.
    if(swCompT) then
      call compute_T_i(LNS_rho0(iglobM)+out_drho_P, out_v_P, LNS_E0(iglobM)+out_dE_P, out_T_P)
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
    
    ! Set out_drho_P.
    out_drho_P = inp_drho_P
    
    ! Set out_dE_P.
    out_dE_P = inp_dE_P
    
    ! Set out_rho0dv_P.
    out_rho0dv_P = inp_rho0dv_P
    
    ! Set out_dv_P.
    out_dv_P(:) = out_rho0dv_P(:)/LNS_rho0(iglobP)
    
    ! Set out_dp_P.
    call compute_dp_i(LNS_rho0(iglobM)+out_drho_P, LNS_v0(:,iglobM)+out_dv_P, LNS_E0(iglobM)+out_dE_P, out_dp_P, iglobM)
    !out_dp_P       = (gamma_P - ONEcr)*( out_dE_P & ! Warning, expression of out_dp_P might not be exact.
    !               - (HALFcr)*out_drho_P*( out_dv_P(1)**2 + out_dv_P(SPACEDIM)**2 ) )
    
    ! Set out_dm_P: see bottom of routine.
    if(swCompVisc) then
      ! Set out_nabla_dT_P: same as other side.
      out_nabla_dT_P = nabla_dT(:,iglobM)
      ! Set out_sigma_dv_P: same as other side.
      out_sigma_dv_P = sigma_dv(:,iglobM)
    endif
    
    ! Set out_v_P.
    out_v_P = out_dv_P + LNS_v0(:,iglobM)
    
    ! Set out_T_P.
    if(swCompT) then
      call compute_T_i(LNS_rho0(iglobM)+out_drho_P, out_v_P, LNS_E0(iglobM)+out_dE_P, out_T_P)
    endif
    !out_dT_P = (out_dE_P/out_drho_P - 0.5*((out_rho0dv_P(1)/out_drho_P)**2 + (out_rho0dv_P(SPACEDIM)/out_drho_P)**2))/c_V
    
  endif
  
  ! Set out_dm_P.
  out_dm_P=out_rho0dv_P+out_drho_P*LNS_rho0(iglobM)
  
end subroutine LNS_get_interfaces_unknowns



! ------------------------------------------------------------ !
! build_trans_boundary                                         !
! ------------------------------------------------------------ !
subroutine build_trans_boundary(normal, tangential, transf_matrix)
  use constants, only: CUSTOM_REAL
  use specfem_par_LNS, only: SPACEDIM
  
  implicit none
  
  ! Input/Output.
  real(kind=CUSTOM_REAL), dimension(SPACEDIM), intent(in) :: normal
  real(kind=CUSTOM_REAL), dimension(SPACEDIM), intent(out) :: tangential
  real(kind=CUSTOM_REAL), dimension(SPACEDIM, SPACEDIM), intent(out) :: transf_matrix
  
  ! Local.
  ! N./A.
  
  ! For the tangential vector, it is assumed that we will only have the elastic media under the DG medium, hence we always have n_out(SPACEDIM)>=0.
  tangential(1)                     = -normal(SPACEDIM) ! Recall: normal(SPACEDIM)=n_z.
  tangential(SPACEDIM)              =  normal(1) ! Recall: normal(1)=n_x.
  transf_matrix(1,        1)        =   tangential(SPACEDIM)
  transf_matrix(1,        SPACEDIM) = - normal(SPACEDIM) ! Recall: normal(SPACEDIM)=n_z.
  transf_matrix(SPACEDIM, 1)        = - tangential(1)
  transf_matrix(SPACEDIM, SPACEDIM) =   normal(1) ! Recall: normal(1)=n_x.
  transf_matrix = transf_matrix/(normal(1)*tangential(SPACEDIM) - tangential(1)*normal(SPACEDIM))
end subroutine 















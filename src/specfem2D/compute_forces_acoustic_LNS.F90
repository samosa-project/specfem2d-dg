! ------------------------------------------------------------ !
! compute_forces_acoustic_LNS                                  !
! ------------------------------------------------------------ !
! TODO: Description.
subroutine compute_forces_acoustic_LNS(cv_drho, cv_rho0dv, cv_dE, & ! Constitutive variables.
                                       in_dm, in_dp, in_dT, in_nabla_dT, in_sigma_dv, & ! Precomputed quantities.
                                       outrhs_drho, outrhs_rho0dv, outrhs_dE, & ! Output.
                                       currentTime) ! Time.
  use constants ! TODO: select variables to use. , only: CUSTOM_REAL,NGLLX,NGLLZ,SPACEDIM
  use specfem_par ! TODO: select variables to use. , only: nglob_DG,nspec, ispec_is_acoustic_DG,&
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
  use specfem_par_LNS ! TODO: select variables to use.
  
  implicit none
  
  ! Input/Output.
  real(kind=CUSTOM_REAL), dimension(nglob_DG), intent(in) :: cv_drho, cv_dE, in_dp, in_dT
  real(kind=CUSTOM_REAL), dimension(SPACEDIM, nglob_DG), intent(in) :: cv_rho0dv, in_dm, in_nabla_dT
  real(kind=CUSTOM_REAL), dimension(3, nglob_DG), intent(in) :: in_sigma_dv
  real(kind=CUSTOM_REAL), dimension(nglob_DG), intent(out) :: outrhs_drho, outrhs_dE
  real(kind=CUSTOM_REAL), dimension(SPACEDIM, nglob_DG), intent(out) :: outrhs_rho0dv
  real(kind=CUSTOM_REAL), dimension(nglob_DG), intent(in) :: currentTime
  
  ! Local.
  real(kind=CUSTOM_REAL), parameter :: ZEROcr = 0._CUSTOM_REAL
  !real(kind=CUSTOM_REAL), parameter :: ONE  = 1._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: TWOcr  = 2._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: HALFcr = 0.5_CUSTOM_REAL
  
  ! MESSY
  !real(kind=CUSTOM_REAL), dimension(nglob_DG) :: rho_DG, rhovx_DG, rhovz_DG, E_DG
  !real(kind=CUSTOM_REAL), dimension(2, nglob_DG) :: T_DG
  !real(kind=CUSTOM_REAL), dimension(2, 2, nglob_DG) :: V_DG
  integer :: ispec, i,j, k,iglob!, iglob_unique
  real(kind=CUSTOM_REAL), dimension(NGLLX, NGLLZ, 2) :: cntrb_drho, cntrb_dE
  real(kind=CUSTOM_REAL), dimension(SPACEDIM, NGLLX, NGLLZ, 2) :: cntrb_rho0dv
  ! Variables "cntrb_*" are aimed at assembling the different contributions to constitutive variables.
  ! For drho and dE, "cntrb_*(:,:,1)" contains the contribution along xi, and "cntrb_*(:,:,2)" the contribution along gamma.
  ! For rho0dv, "cntrb_*(1,:,:,1)" contains the contribution to rho0dvx along xi, "cntrb_*(1,:,:,2)" the contribution to rho0dvx along gamma, "cntrb_*(2,:,:,1)" the contribution to rhodvz along xi, and "cntrb_*(2,:,:,2)" the contribution to rhodvz along gamma.
  real(kind=CUSTOM_REAL), dimension(NGLLX, NGLLZ) :: d0cntrb_dE
  real(kind=CUSTOM_REAL), dimension(SPACEDIM, NGLLX,NGLLZ) :: d0cntrb_rho0dv
  ! Variable "d0cntrb_*" are aimed at assembling contributions of zero-th degree (outside divergence operator).
  
  
  real(kind=CUSTOM_REAL) :: xixl,xizl,gammaxl,gammazl,jacLoc ! Jacobian matrix and determinant.
  real(kind=CUSTOM_REAL) :: lambda, nx, nz, weight, &
                            tmp_unknown_x_M, tmp_unknown_x_P, tmp_unknown_z_M, tmp_unknown_z_P, &
                            tmp_unknown_x, tmp_unknown_z, &
                            flux_n, jump, &
                            drho_P, veloc_x_DG_P, veloc_z_DG_P, &
                            dE_P, in_dp_P, &!currentTime, &
                            Tx_DG_P, Tz_DG_P, Vxx_DG_P, Vzz_DG_P, Vxz_DG_P, Vzx_DG_P, &!T_P, &
                            ! TEST
                            gamma_P!,e1_DG_P
  real(kind=CUSTOM_REAL), dimension(2) :: rho0dv_P
  real(kind=CUSTOM_REAL) :: dT_dx, dT_dz
  !real(kind=CUSTOM_REAL) :: veloc_n_M, veloc_n_P
  integer :: iglobM, iglobP
  integer, dimension(3) :: neighbor
  !real(kind=CUSTOM_REAL), dimension(nglob_DG) :: veloc_x_DG, veloc_z_DG, in_dp
  real(kind=CUSTOM_REAL), dimension(2, nglob_DG) :: LNS_dv
  real(kind=CUSTOM_REAL) :: dux_dx, dux_dz, duz_dx, duz_dz ! Viscosity
  real(kind=CUSTOM_REAL) :: wxl, wzl ! Integration weigths.
  logical :: exact_interface_flux
  !integer :: cnst_hdrsttc ! For computationnaly better CONSTRAIN_HYDROSTATIC switches.
  integer :: iface1, iface, iface1_neighbor, iface_neighbor, ispec_neighbor
  !real(kind=CUSTOM_REAL) :: ya_x_l, ya_z_l ! Stretching absorbing boundary conditions.
  
  ! TESTS
  !real(kind=CUSTOM_REAL) :: x,z ! Artifical advection and a priori damping.
  !real(kind=CUSTOM_REAL) :: maxval_rho,maxval_rhovx,maxval_rhovz,maxval_E ! ABSORB
  !logical :: ABSORB_BC ! ABSORB
  
  ! For more convinient CONSTRAIN_HYDROSTATIC switches.
  ! TODO: Replace the CONSTRAIN_HYDROSTATIC switches using this variable.
  !if(CONSTRAIN_HYDROSTATIC) then
  !  cnst_hdrsttc=ONE
  !else
  !  cnst_hdrsttc=ZERO
  !endif
  
  !T_DG = T_DG_main
  !V_DG = V_DG_main
  
  ! Initialise auxiliary unknowns from constitutive variables.
  LNS_dv = 0.
  !LNS_dv = cv_rho0dv(1,:)/LNS_rho0 ! TODO: Element-wise division.
  !p_DG       = (gammaext_DG - ONE)*( cv_dE &
  !             - (HALFcr)*cv_drho*( veloc_x_DG**2 + veloc_z_DG**2 ) )
  
  ! Initialisation of the RHS.
  outrhs_drho   = ZEROcr
  outrhs_rho0dv = ZEROcr
  outrhs_dE     = ZEROcr
  
  ! Start by adding source terms.
  ! TODO: dedicated routines
  if(TYPE_SOURCE_DG == 1) then
    call compute_add_sources_acoustic_DG_spread(outrhs_drho, it, i_stage)   
  elseif(TYPE_SOURCE_DG == 2) then
    call compute_add_sources_acoustic_DG_spread(outrhs_rho0dv(1,:), it, i_stage)
    call compute_add_sources_acoustic_DG_spread(outrhs_rho0dv(2,:), it, i_stage)
  elseif(TYPE_SOURCE_DG == 3) then
    call compute_add_sources_acoustic_DG_spread(outrhs_dE, it, i_stage)
  endif
  
  if(myrank == 0 .and. LNS_VERBOSE>30 .and. mod(it, 100) == 0) then
    write(*,"(a)")                   " quantity: [max                     , min                     ]"
    WRITE(*,"(a,e24.16,a,e24.16,a)") " drho:     [", maxval(cv_drho), ", ", minval(cv_drho), "]"
    WRITE(*,"(a,e24.16,a,e24.16,a)") " rho0dvx:  [", maxval(cv_rho0dv(1,:)), ", ", minval(cv_rho0dv(1,:)), "]"
    WRITE(*,"(a,e24.16,a,e24.16,a)") " rho0dvz:  [", maxval(cv_rho0dv(1,:)), ", ", minval(cv_rho0dv(1,:)), "]"
    WRITE(*,"(a,e24.16,a,e24.16,a)") " dE        [", maxval(cv_dE), ", ", minval(cv_dE), "]"
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
          
          ! Inviscid stress tensor's contributions.
          !   Mass conservation.
          tmp_unknown_x = in_dm(1,iglob)
          tmp_unknown_z = in_dm(2,iglob)
          cntrb_drho(i,j,1) = wzl * jacLoc * (xixl * tmp_unknown_x + xizl * tmp_unknown_z) ! Contribution along xi.
          cntrb_drho(i,j,2) = wxl * jacLoc * (gammaxl * tmp_unknown_x + gammazl * tmp_unknown_z) ! Contribution along gamma.
          !   x-Momentum.
          tmp_unknown_x = cv_rho0dv(1,iglob)*LNS_v0(1,iglob) + in_dp(iglob)
          tmp_unknown_z = cv_rho0dv(2,iglob)*LNS_v0(1,iglob)
          cntrb_rho0dv(1,i,j,1) = wzl * jacLoc * (xixl * tmp_unknown_x + xizl * tmp_unknown_z) ! Contribution along xi.
          cntrb_rho0dv(1,i,j,2) = wxl * jacLoc * (gammaxl * tmp_unknown_x + gammazl * tmp_unknown_z) ! Contribution along gamma.
          !   z-Momentum.
          tmp_unknown_x = cv_rho0dv(1,iglob)*LNS_v0(2,iglob)
          tmp_unknown_z = cv_rho0dv(2,iglob)*LNS_v0(2,iglob) + in_dp(iglob)
          cntrb_rho0dv(2,i,j,1) = wzl * jacLoc * (xixl * tmp_unknown_x + xizl * tmp_unknown_z) ! Contribution along xi.
          cntrb_rho0dv(2,i,j,2) = wxl * jacLoc * (gammaxl * tmp_unknown_x + gammazl * tmp_unknown_z) ! Contribution along gamma.
          
          ! Add viscous stress tensor's contributions.
          if(LNS_mu(iglob) > 0 .OR. &
             LNS_eta(iglob) > 0 .OR. &
             LNS_kappa(iglob) > 0) then
            !   Mass conservation: no viscous contribution.
            !   x-Momentum.
            tmp_unknown_x = -in_sigma_dv(1,iglob) ! Recall, this corresponds to \Sigma_v(v')_{1,1}.
            tmp_unknown_z = -in_sigma_dv(2,iglob) ! Recall, this corresponds to \Sigma_v(v')_{1,2} (and \Sigma_v(v')_{2,1}).
            cntrb_rho0dv(1,i,j,1) = cntrb_rho0dv(1,i,j,1) + wzl * jacLoc * (xixl * tmp_unknown_x + xizl * tmp_unknown_z) ! Contribution along xi.
            cntrb_rho0dv(1,i,j,2) = cntrb_rho0dv(1,i,j,2) + wxl * jacLoc * (gammaxl * tmp_unknown_x + gammazl * tmp_unknown_z) ! Contribution along gamma.
            !   z-Momentum.
            tmp_unknown_x = -in_sigma_dv(2,iglob) ! Recall, this corresponds to \Sigma_v(v')_{2,1} (and \Sigma_v(v')_{1,2}).
            tmp_unknown_z = -in_sigma_dv(3,iglob) ! Recall, this corresponds to \Sigma_v(v')_{2,2}.
            cntrb_rho0dv(2,i,j,1) = cntrb_rho0dv(2,i,j,1) + wzl * jacLoc * (xixl * tmp_unknown_x + xizl * tmp_unknown_z) ! Contribution along xi.
            cntrb_rho0dv(2,i,j,2) = cntrb_rho0dv(2,i,j,2) + wxl * jacLoc * (gammaxl * tmp_unknown_x + gammazl * tmp_unknown_z) ! Contribution along gamma.
          endif
          
          ! Special case: energy.
          ! Indeed, if they exist, viscous contributions can be grouped to inviscid ones.
          if(LNS_mu(iglob) > 0 .OR. &
             LNS_eta(iglob) > 0 .OR. &
             LNS_kappa(iglob) > 0) then
            tmp_unknown_x =   LNS_dv(1,iglob)*(LNS_E0(iglob) + LNS_p0(iglob) - sigma_v_0(1,iglob)) &
                            + LNS_v0(1,iglob)*(cv_dE(iglob) + in_dp(iglob) - in_sigma_dv(1,iglob)) &
                            - LNS_dv(2,iglob)*sigma_v_0(2,iglob) &
                            - LNS_v0(2,iglob)*in_sigma_dv(2,iglob) &
                            + LNS_kappa(iglob)*in_nabla_dT(1,iglob)
            tmp_unknown_z =   LNS_dv(2,iglob)*(LNS_E0(iglob) + LNS_p0(iglob) - sigma_v_0(3,iglob)) &
                            + LNS_v0(2,iglob)*(cv_dE(iglob) + in_dp(iglob) - in_sigma_dv(3,iglob)) &
                            - LNS_dv(1,iglob)*sigma_v_0(2,iglob) &
                            - LNS_v0(1,iglob)*in_sigma_dv(2,iglob) &
                            + LNS_kappa(iglob)*in_nabla_dT(2,iglob)
          else
            tmp_unknown_x =   LNS_dv(1,iglob)*(LNS_E0(iglob) + LNS_p0(iglob)) &
                            + LNS_v0(1,iglob)*(cv_dE(iglob) + in_dp(iglob))
            tmp_unknown_z =   LNS_dv(2,iglob)*(LNS_E0(iglob) + LNS_p0(iglob)) &
                            + LNS_v0(2,iglob)*(cv_dE(iglob) + in_dp(iglob))
          endif
          cntrb_dE(i,j,1) = wzl * jacLoc * (xixl * tmp_unknown_x + xizl * tmp_unknown_z) ! Contribution along xi.
          cntrb_dE(i,j,2) = wxl * jacLoc * (gammaxl * tmp_unknown_x + gammazl * tmp_unknown_z) ! Contribution along gamma.
          
          ! Zero-th degree contributions.
          !d0cntrb_drho(i,j)     = ZERO
          d0cntrb_rho0dv(1,i,j) = (  cv_drho(iglob)*potential_dphi_dx_DG(ibool(i,j,ispec)) & ! \rho'g_x
                                   + in_dm(1,iglob)*nabla_v0(1,1,iglob) & ! {\delta_m}_x\partial_xv_{0,x}
                                   + in_dm(2,iglob)*nabla_v0(1,2,iglob)) * jacLoc ! {\delta_m}_z\partial_zv_{0,x}
          d0cntrb_rho0dv(2,i,j) = (  cv_drho(iglob)*potential_dphi_dz_DG(ibool(i,j,ispec)) & ! \rho'g_z
                                   + in_dm(1,iglob)*nabla_v0(2,1,iglob) & ! {\delta_m}_x\partial_xv_{0,z}
                                   + in_dm(2,iglob)*nabla_v0(2,2,iglob)) * jacLoc ! {\delta_m}_z\partial_zv_{0,z}
          d0cntrb_dE(i,j)       = (  potential_dphi_dx_DG(ibool(i,j,ispec))*in_dm(1,iglob) &
                                   + potential_dphi_dz_DG(ibool(i,j,ispec))*in_dm(2,iglob)) * jacLoc
        enddo
      enddo
      
      ! Assemble the contributions previously computed, and add gravity's contribution.
      ! The integration by quadrature on the GLL points leads to three sums. See in particular Komatitsch (Méthodes spectrales et éléments spectraux pour l'équation de l'élastodynamique 2D et 3D en milieu hétérogène), Annexe 3.A.
      do j = 1, NGLLZ
        do i = 1, NGLLX
          iglob = ibool_DG(i,j,ispec)
          do k = 1, NGLLX
            outrhs_drho(iglob)     =   outrhs_drho(iglob) &
                                  + (  cntrb_drho(k,j,1) * real(hprimewgll_xx(k,i), kind=CUSTOM_REAL) &
                                     + cntrb_drho(i,k,2) * real(hprimewgll_zz(k,j), kind=CUSTOM_REAL))
            outrhs_rho0dv(1,iglob) =   outrhs_rho0dv(1,iglob) &
                                  + (  cntrb_rho0dv(1,k,j,1) * real(hprimewgll_xx(k,i), kind=CUSTOM_REAL) &
                                     + cntrb_rho0dv(1,i,k,2) * real(hprimewgll_zz(k,j), kind=CUSTOM_REAL))
            outrhs_rho0dv(2,iglob) =   outrhs_rho0dv(2,iglob) &
                                  + (  cntrb_rho0dv(2,k,j,1) * real(hprimewgll_xx(k,i), kind=CUSTOM_REAL) &
                                     + cntrb_rho0dv(2,i,k,2) * real(hprimewgll_zz(k,j), kind=CUSTOM_REAL))
            outrhs_dE(iglob)       =   outrhs_dE(iglob) &
                                  + (  cntrb_dE(k,j,1) * real(hprimewgll_xx(k,i), kind=CUSTOM_REAL) &
                                     + cntrb_dE(i,k,2) * real(hprimewgll_zz(k,j), kind=CUSTOM_REAL))
          enddo ! Enddo on k.
          
          wzl = real(wzgll(j), kind=CUSTOM_REAL)
          wxl = real(wxgll(i), kind=CUSTOM_REAL)
          
          ! Add zero-th order terms.
          !outrhs_drho(iglob)     = outrhs_drho(iglob)     + d0cntrb_drho(i,j) * wxl * wzl
          outrhs_rho0dv(1,iglob) = outrhs_rho0dv(1,iglob) + d0cntrb_rho0dv(1,i,j) * wxl * wzl
          outrhs_rho0dv(2,iglob) = outrhs_rho0dv(2,iglob) + d0cntrb_rho0dv(2,i,j) * wxl * wzl
          outrhs_dE(iglob)       = outrhs_dE(iglob)       + d0cntrb_dE(i,j) * wxl * wzl
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
          
          ! Step 1: prepare the normals' parameters (nx, nz, weight, etc.).
          ! Interior point
          iglobM = ibool_DG(i,j,ispec)
          
          ! TEST WITH IFACE FORMULATION
          nx     = nx_iface(iface, ispec)
          nz     = nz_iface(iface, ispec)
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
          drho_P     = ZERO
          rho0dv_P   = ZERO
          dE_P       = ZERO
          veloc_x_DG_P = ZERO
          veloc_z_DG_P = ZERO
          in_dp_P       = ZERO
          
          iglobP = 1
          if(neighbor(1) > -1) then
            iglobP = ibool_DG(neighbor(1), neighbor(2), neighbor(3))
          endif
          
          exact_interface_flux = .false. ! Reset this variable to .false.: by default, the fluxes have to be computed (jump!=0). In some specific cases (assigned during the call to compute_interface_unknowns), the flux can be exact (jump==0).
          ! TODO: dedicated routine.
          !call compute_interface_unknowns(i,j,ispec, drho_P, rho0dv_P(1), &
          !        rho0dv_P(2), dE_P, veloc_x_DG_P, veloc_z_DG_P, in_dp_P, T_P, &
          !        Tx_DG_P, Tz_DG_P, Vxx_DG_P, Vzz_DG_P, Vzx_DG_P, Vxz_DG_P, gamma_P,&
          !        neighbor, &
          !        exact_interface_flux, &
          !        cv_drho(iglobM), cv_dE(iglobM), cv_rho0dv(1,iglobM), cv_rho0dv(2,iglobM), &
          !        V_DG(:,:,iglobM), T_DG(:,iglobM), &
          !        cv_drho(iglobP), cv_dE(iglobP), cv_rho0dv(1,iglobP), cv_rho0dv(2,iglobP), &
          !        V_DG(:,:,iglobP), T_DG(:,iglobP), &
          !        nx, nz, weight, currentTime, iface1, iface)
          !        !TEST STRETCH
          !        !nx_unit, nz_unit, weight, currentTime, iface1, iface)
          
          ! Recover an approximate local maximum linearized acoustic wave speed. See for example Hesthaven (doi.org/10.1007/9780387720678), page 208.
          lambda = ZERO
          jump   = ZERO
          !veloc_n_M = abs(LNS_dv(1,iglobM)*nx + LNS_dv(2,iglobM)*nz)
          !veloc_n_P = abs(veloc_x_DG_P*nx + veloc_z_DG_P*nz)
          ! Save some operations.
          !lambda = max( veloc_n_M + sqrt(abs(gammaext_DG(iglobM)*in_dp(iglobM)/cv_drho(iglobM))), &
          !              veloc_n_P + sqrt(abs(gamma_P*in_dp_P/drho_P)) )
          lambda = max(abs(LNS_dv(1,iglobM)*nx + LNS_dv(2,iglobM)*nz) &
                       + sqrt(abs(gammaext_DG(iglobM)*in_dp(iglobM)/cv_drho(iglobM))), &
                       abs(veloc_x_DG_P*nx + veloc_z_DG_P*nz) &
                       + sqrt(abs(gamma_P*in_dp_P/drho_P)))
          
          ! Viscous stress tensor's contributions (already under the form of the average mean flux).
          dux_dx = ZERO
          dux_dz = ZERO
          duz_dx = ZERO
          duz_dz = ZERO
          dT_dx = ZERO
          dT_dz = ZERO
          if(LNS_mu(iglob) > 0 .OR. &
             LNS_eta(iglob) > 0 .OR. &
             LNS_kappa(iglob)  > 0) then
            dux_dx = HALFcr*(V_DG(1, 1, iglobM) + Vxx_DG_P)
            dux_dz = HALFcr*(V_DG(1, 2, iglobM) + Vxz_DG_P)
            duz_dx = HALFcr*(V_DG(2, 1, iglobM) + Vzx_DG_P)
            duz_dz = HALFcr*(V_DG(2, 2, iglobM) + Vzz_DG_P)
            dT_dx = HALFcr*(T_DG(1, iglobM) + Tx_DG_P)
            dT_dz = HALFcr*(T_DG(2, iglobM) + Tz_DG_P)
          endif
          
          ! Mass conservation equation's contributions (which are fully inviscid).
          !tmp_unknown_x_M  = cv_rho0dv(1,iglobM)
          !tmp_unknown_x_P  = rho0dv_P
          !tmp_unknown_z_M = cv_rho0dv(2,iglobM)
          !tmp_unknown_z_P = rho0dv_P
          ! Dot product.
          flux_n = (cv_rho0dv(1,iglobM)+rho0dv_P(1))*nx + (cv_rho0dv(2,iglobM)+rho0dv_P(2))*nz
          !flux_n = (tmp_unknown_x_M+tmp_unknown_x_P)*nx + (tmp_unknown_z_M+tmp_unknown_z_P)*nz ! [5 operations + 1 affectation], instead of [5 operations + 3 affectations].
          jump   = cv_drho(iglobM) - drho_P
          ! Add flux' contribution.
          if(exact_interface_flux) then
            jump = ZERO
          endif
          outrhs_drho(iglobM) = outrhs_drho(iglobM) - weight*(flux_n + lambda*jump)*HALFcr
          
          ! x-Momentum equation's inviscid contributions.
          tmp_unknown_x_M = cv_drho(iglobM)*LNS_dv(1,iglobM)**2 + in_dp(iglobM)
          tmp_unknown_x_P = drho_P*veloc_x_DG_P**2 + in_dp_P
          tmp_unknown_z_M = cv_drho(iglobM)*LNS_dv(1,iglobM)*LNS_dv(2,iglobM)
          tmp_unknown_z_P = drho_P*veloc_x_DG_P*veloc_z_DG_P
          ! Dot product.
          flux_n = (tmp_unknown_x_M+tmp_unknown_x_P)*nx + (tmp_unknown_z_M+tmp_unknown_z_P)*nz ! [5 operations + 1 affectation], instead of [5 operations + 3 affectations].
          jump   = cv_rho0dv(1,iglobM) - rho0dv_P(1)
          ! Add flux' contribution.
          if(exact_interface_flux) then
            jump = ZERO
          endif
          outrhs_rho0dv(1,iglobM) = outrhs_rho0dv(1,iglobM) - weight*(flux_n + lambda*jump)*HALFcr

          ! x-Momentum equation's viscous contributions.
          ! 1st line of \Sigma_v.
          ! The vector [tmp_unknown_x, tmp_unknown_z] represents the mean average flux at the boundary of the x-momentum.
          ! Recall: dux_dx, duz_dx, dux_dz, and duz_dz already contain the 0.5 factor to put the flux under mean average form.
          tmp_unknown_x = LNS_mu(iglob)*TWOcr*dux_dx + (LNS_eta(iglob) - (TWOcr/3.)*LNS_mu(iglob))*(dux_dx + duz_dz) 
          tmp_unknown_z = LNS_mu(iglob)*( dux_dz + duz_dx )
          ! Dot product.
          !flux_n = tmp_unknown_x*nx + tmp_unknown_z*nz ! [3 operations + 1 affectation], instead of [3 operations + 3 affectations].
          !outrhs_rho0dv(1,iglobM) = outrhs_rho0dv(1,iglobM) + weight*flux_n
          outrhs_rho0dv(1,iglobM) = outrhs_rho0dv(1,iglobM) + weight*(tmp_unknown_x*nx+tmp_unknown_z*nz) ! [1 affectation], instead of [2 affections].
          ! The computed values contained in the variables tmp_unknown_x and tmp_unknown_z can be used to compute the energy's x component of the mean average flux at the boundary. Thus, we add this contribution here.
          outrhs_dE(iglobM)     = outrhs_dE(iglobM) &
                              + weight * HALFcr * (  (LNS_dv(1,iglobM) + veloc_x_DG_P) * tmp_unknown_x &
                                                  +(LNS_dv(2,iglobM) + veloc_z_DG_P) * tmp_unknown_z )*nx
          
          ! z-Momentum equation's inviscid contributions.
          tmp_unknown_x_M  = cv_drho(iglobM)*LNS_dv(1,iglobM)*LNS_dv(2,iglobM)
          tmp_unknown_x_P  = drho_P*veloc_x_DG_P*veloc_z_DG_P
          !if(.not. CONSTRAIN_HYDROSTATIC) then        
            tmp_unknown_z_M = cv_drho(iglobM)*LNS_dv(2,iglobM)**2 + in_dp(iglobM)
            tmp_unknown_z_P = drho_P*veloc_z_DG_P**2 + in_dp_P
          !else
          !  tmp_unknown_z_M = cv_drho(iglobM)*LNS_dv(2,iglobM)**2 + (in_dp(iglobM) - in_dp_init(iglobM))
          !  tmp_unknown_z_P = drho_P*veloc_z_DG_P**2 + (in_dp_P - in_dp_init(iglobM))
          !endif
          ! Dot product.
          flux_n = (tmp_unknown_x_M+tmp_unknown_x_P)*nx + (tmp_unknown_z_M+tmp_unknown_z_P)*nz ! [5 operations + 1 affectation], instead of [5 operations + 3 affectations].
          jump   = cv_rho0dv(2,iglobM) - rho0dv_P(2)
          ! Add flux' contribution.
          if(exact_interface_flux) then
            jump = ZERO
          endif
          outrhs_rho0dv(2,iglobM) = outrhs_rho0dv(2,iglobM) - weight*(flux_n + lambda*jump)*HALFcr
          
          ! z-Momentum equation's viscous contributions.
          ! 2nd line of \Sigma_v.
          ! The vector [tmp_unknown_x, tmp_unknown_z] represents the mean average flux at the boundary of the z-momentum.
          ! Recall: dux_dx, duz_dx, dux_dz, and duz_dz already contain the 0.5 factor to put the flux under mean average form.
          tmp_unknown_x = LNS_mu(iglob)*( dux_dz + duz_dx )
          tmp_unknown_z = LNS_mu(iglob)*TWOcr*duz_dz + (LNS_eta(iglob) - (TWOcr/3.)*LNS_mu(iglob))*(dux_dx + duz_dz)
          ! Dot product.
          !flux_n = tmp_unknown_x*nx + tmp_unknown_z*nz ! [3 operations + 1 affectation], instead of [3 operations + 3 affectations].
          !outrhs_rho0dv(2,iglobM) = outrhs_rho0dv(2,iglobM) + weight*flux_n
          outrhs_rho0dv(2,iglobM) = outrhs_rho0dv(2,iglobM) + weight*(tmp_unknown_x*nx+tmp_unknown_z*nz) ! [1 affectation], instead of [2 affections].
          ! The computed values contained in the variables tmp_unknown_x and tmp_unknown_z can be used to compute the energy's z component of the mean average flux at the boundary. Thus, we add this contribution here.
          outrhs_dE(iglobM)     = outrhs_dE(iglobM) &
                              + weight * HALFcr * (  (LNS_dv(1,iglobM) + veloc_x_DG_P) * tmp_unknown_x &
                                                  +(LNS_dv(2,iglobM) + veloc_z_DG_P) * tmp_unknown_z )*nz

          ! Energy equation's fully inviscid contributions.
          tmp_unknown_x_M = LNS_dv(1,iglobM)*(cv_dE(iglobM) + in_dp(iglobM))
          tmp_unknown_x_P = veloc_x_DG_P*(dE_P + in_dp_P)
          tmp_unknown_z_M = LNS_dv(2,iglobM)*(cv_dE(iglobM) + in_dp(iglobM))
          tmp_unknown_z_P = veloc_z_DG_P*(dE_P + in_dp_P)
          ! Dot product.
          flux_n = (tmp_unknown_x_M+tmp_unknown_x_P)*nx + (tmp_unknown_z_M+tmp_unknown_z_P)*nz ! [5 operations + 1 affectation], instead of [5 operations + 3 affectations].
          jump   = cv_dE(iglobM) - dE_P
          ! Add flux' contribution.
          if(exact_interface_flux) then
            jump = ZERO
          endif
          outrhs_dE(iglobM) = outrhs_dE(iglobM) - weight*(flux_n + lambda*jump)*HALFcr
          
          ! Energy equation's heat flux' contribution (last remaining term, viscous).
          ! Recall: dT_dx already contains the 0.5 factor to put the flux under mean average form.
          outrhs_dE(iglobM) = outrhs_dE(iglobM) &
                          + weight*( LNS_kappa(iglob)*( dT_dx*nx + dT_dz*nz ) )
        enddo
      enddo
    endif ! End of test if acoustic element
  enddo ! End of loop on elements.
  
end subroutine compute_forces_acoustic_LNS

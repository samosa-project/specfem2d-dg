!========================================================================
!
!                   S P E C F E M 2 D  Version 7 . 0
!                   --------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!                and CNRS / University of Marseille, France
!                 (there are currently many more authors!)
! (c) Princeton University and CNRS / University of Marseille, April 2014
!
! This software is a computer program whose purpose is to solve
! the two-dimensional viscoelastic anisotropic or poroelastic wave equation
! using a spectral-element method (SEM).
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
! The full text of the license is available in file "LICENSE".
!
!========================================================================

  ! Compute initial condition
  subroutine initial_condition_DG()

! compute forces in the acoustic elements in forward simulation and in adjoint simulation in adjoint inversion

  use constants,only: CUSTOM_REAL,NGLLX,NGLLZ,PI

  use specfem_par, only: nglob_DG, nglob, ibool_DG, &
        rho_DG, rhovx_DG, rhovz_DG, E_DG, e1_DG, nspec, &
        potential_dphi_dx_DG, potential_dphi_dz_DG, &
        hprime_xx, hprime_zz, jacobian, xix, xiz, gammax, gammaz, ibool, &
        gammaext_DG, rhoext, windxext, coord, myrank, pext_DG

  implicit none
  
  ! Parameters
  real(kind=CUSTOM_REAL), parameter :: ZEROl = 0._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: ONEl  = 1._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: TWOl  = 2._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: SIXl  = 6._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: HALFl = 0.5_CUSTOM_REAL
  
  integer :: ispec, iglob, i, j
  
  real(kind=CUSTOM_REAL) :: rho_DG_P, rhovx_DG_P, rhovz_DG_P, E_DG_P, &
      veloc_x_DG_P, veloc_z_DG_P, p_DG_P, e1_DG_P
  
    
  do ispec = 1,nspec

        ! first double loop over GLL points to compute and store gradients
        do j = 1,NGLLZ
          do i = 1,NGLLX
          
           iglob = ibool_DG(i,j,ispec)
  
           call boundary_condition_DG(i, j, ispec, ZEROl, rho_DG_P, rhovx_DG_P, rhovz_DG_P, E_DG_P, &
                veloc_x_DG_P, veloc_z_DG_P, p_DG_P, e1_DG_P)
  
           rho_DG(iglob)     = rho_DG_P
           rhovx_DG(iglob)   = rhovx_DG_P
           rhovz_DG(iglob)   = rhovz_DG_P
           E_DG(iglob)       = E_DG_P
           e1_DG(iglob)      = e1_DG_P
           
          enddo
       enddo
       
  enddo  
  
  if(.false.) call recompute_density()
  
  end subroutine initial_condition_DG
  
  subroutine recompute_density()

  use constants,only: CUSTOM_REAL, NGLLX, NGLLZ

  use specfem_par, only: xix, xiz, gammax, gammaz, coord, ibool_before_perio, &
        jacobian, hprime_xx, hprime_zz, ibool, ibool_DG, nglob, gravityext, &
        E_DG, rhovx_DG, rhovz_DG, rho_DG, gammaext_DG, pext_DG, vpext, windxext, myrank, nspec, rhoext, &
        nglob_DG, neighbor_DG, neighbor_DG_corner, V_DG, T_DG, rmass_inverse_acoustic_DG,elastic_tensor, &
        ispec_is_acoustic, normal_DG, wzgll, wxgll, hprimewgll_zz, normal_DG_corner, weight_DG_corner, &
        hprimewgll_xx, weight_DG, dir_normal_DG_corner, dir_normal_DG, is_corner , potential_dphi_dz_DG, &
        hprime_xx, hprime_zz    

  implicit none
  
  ! local parameters
  integer :: ispec,i,j,k,iglob
  integer :: ifirstelem,ilastelem

  ! Jacobian matrix and determinant
  real(kind=CUSTOM_REAL) :: xixl,xizl,gammaxl,gammazl
  
  ! Viscosity
  real(kind=CUSTOM_REAL) :: duz_dxi, duz_dgamma, duz_dz
  
  ! Parameters
  real(kind=CUSTOM_REAL), parameter :: ZERO = 0._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: ONE  = 1._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: TWO  = 2._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: HALF = 0.5_CUSTOM_REAL
  
  real(kind=CUSTOM_REAL), dimension(nglob_DG) :: p_DG

  p_DG       = (gammaext_DG - ONE)*( E_DG &
        - (HALF/rho_DG)*( rhovx_DG**2 + rhovz_DG**2 ) )
  
  ifirstelem = 1
  ilastelem = nspec

   do ispec = ifirstelem,ilastelem

    ! acoustic spectral element
    if (ispec_is_acoustic(ispec)) then

      ! first double loop over GLL points to compute and store gradients
      do j = 1,NGLLZ
        do i = 1,NGLLX
   
          iglob = ibool_DG(i,j,ispec)
   
          duz_dxi    = ZERO
          duz_dgamma = ZERO
          
          xixl = xix(i,j,ispec)
          xizl = xiz(i,j,ispec)
          gammaxl = gammax(i,j,ispec)
          gammazl = gammaz(i,j,ispec)
          
          ! first double loop over GLL points to compute and store gradients
          ! we can merge the two loops because NGLLX == NGLLZ
          do k = 1,NGLLX
                duz_dxi    = duz_dxi    + log(p_DG(ibool_DG(k,j,ispec))) * real(hprime_xx(i,k), kind=CUSTOM_REAL)
                duz_dgamma = duz_dgamma + log(p_DG(ibool_DG(i,k,ispec))) * real(hprime_zz(j,k), kind=CUSTOM_REAL)
          enddo
        
          ! derivatives of velocities
          duz_dz = ( duz_dxi * xizl + duz_dgamma * gammazl )
         
          rhoext(i,j,ispec) = -duz_dz*p_DG(iglob)/gravityext(i,j,ispec)
          rho_DG(iglob)     = rhoext(i,j,ispec)
          
          E_DG(iglob)         = pext_DG(i,j,ispec)/(gammaext_DG(iglob) - 1.) &
                                + 0.5*rho_DG(iglob)*( windxext(i,j,ispec)**2 ) 
          
          rhovx_DG(iglob)    = rho_DG(iglob)*windxext(i,j,ispec)
         
         enddo
      enddo
     
     endif
     
   enddo
   
   
  end subroutine recompute_density
  
  ! -----------------------------------------------------------------------------------------------
  
  ! Compute initial condition
  subroutine boundary_condition_DG(i, j, ispec, timelocal, rho_DG_P, rhovx_DG_P, rhovz_DG_P, E_DG_P, &
      veloc_x_DG_P, veloc_z_DG_P, p_DG_P, e1_DG_P)

! compute forces in the acoustic elements in forward simulation and in adjoint simulation in adjoint inversion

  use constants,only: CUSTOM_REAL,gamma_euler,PI

  use specfem_par, only: ibool_before_perio, ibool_DG, coord, MODEL, &
        rhoext, windxext, pext_DG, gravityext, gammaext_DG, &
        etaext, muext, coord_interface, kappa_DG, cp, cnu, T_init, &
        tau_epsilon, tau_sigma, myrank, &
        rhovx_init, rhovz_init, E_init, rho_init, vpext, &
        ispec_is_acoustic, p_DG_init, T_init, &
        gravity_cte_DG, dynamic_viscosity_cte_DG, thermal_conductivity_cte_DG, tau_eps_cte_DG, tau_sig_cte_DG, SCALE_HEIGHT, &
        USE_ISOTHERMAL_MODEL, potential_dphi_dx_DG, potential_dphi_dz_DG, ibool

  implicit none
  
  integer :: i, j, ispec
  
  real(kind=CUSTOM_REAL) :: rho_DG_P, rhovx_DG_P, rhovz_DG_P, E_DG_P, &
      veloc_x_DG_P, veloc_z_DG_P, p_DG_P, e1_DG_P, timelocal
  
  ! Parameters
  real(kind=CUSTOM_REAL), parameter :: ZEROl = 0._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: ONEl  = 1._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: TWOl  = 2._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: SIXl  = 6._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: HALFl = 0.5_CUSTOM_REAL
  
  real(kind=CUSTOM_REAL) :: x, z, H, G!, A
  !real(kind=CUSTOM_REAL) :: x0, z0, r, beta
  ! Forcing
  real(kind=CUSTOM_REAL) :: to, perio, lambdo, xo
  
  real(kind=CUSTOM_REAL) :: mu_visco, eta_visco
  ! Hydrostatic solution
  !real(kind=CUSTOM_REAL) :: RR, p0, rho0
  ! Density current
  !real(kind=CUSTOM_REAL) :: cp, cnu, exner, RR, p0, rho0, rs, theta, theta0
  ! Linear mountains
  !real(kind=CUSTOM_REAL) :: cp, cnu, exner, RR, p0, rho0, rs, theta, theta0, Nsq
  !real(kind=CUSTOM_REAL) :: Tl, p0,  RR
  ! Thermal bubble
  real(kind=CUSTOM_REAL) :: Tl, Tu, rho0, p0, RR
  !real(kind=CUSTOM_REAL) :: Tl, Tu, p0, rho0, RR, rs!, theta, thetaprime, pibar
  !real(kind=CUSTOM_REAL) :: Tl, Tu, p0, RR, rs, theta0, Nsq, rho0
  !real(kind=CUSTOM_REAL) :: cp, cnu, R, P0, N, theta0, exner, theta
  !real(kind=CUSTOM_REAL) :: x0, z0, r, rs, p0, RR, cnu, cp, theta0, exner
  !real(kind=CUSTOM_REAL) :: x0, z0, r, theta, RR, p0!, exner, theta0
  !real(kind=CUSTOM_REAL) :: &
  !      lz, vs, vs_x, vs_z, Ms, p_1, rho_1, veloc_x_1, veloc_z_1, p_2, rho_2, veloc_x_2, veloc_z_2 
  
  ! Coordinate of the solid-fluid interface
  coord_interface = 304.!250.
  
  x = real(coord(1,ibool_before_perio(i,j,ispec)), kind=CUSTOM_REAL)
  z = real(coord(2,ibool_before_perio(i,j,ispec)), kind=CUSTOM_REAL)
  z = z - coord_interface
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! If an external data file is given for initial conditions
  if(trim(MODEL) == 'external') then
  
   ! Overwrite external model
   rho_DG_P     = rhoext(i,j,ispec)
   p_DG_P       = pext_DG(i,j,ispec)!rho_DG_P*gravityext(i,j,ispec)*Htabext_DG(ibool_DG(i,j,ispec))!pext_DG(i,j,ispec)
   veloc_x_DG_P = windxext(i,j,ispec)
   veloc_z_DG_P = 0.
   
   if(timelocal == 0) then
   !gravityext(i,j,ispec) = G!9.831
   !gammaext_DG(ibool_DG(i,j,ispec))   = 1.4!(rho_DG_P/p_DG_P)*(vpext(i,j,ispec)**2)!gamma_euler!
   !muext(i,j,ispec)  = 1.0d-05
   !etaext(i,j,ispec)  = (4/3)*muext(i,j,ispec)
   !kappa_DG(i,j,ispec) = 0.025!kappa_DG(i,j,ispec) / 10.!0.2
   cp = 7/2
   cnu = 5/2
   
   tau_epsilon(i,j,ispec) = 1.
   tau_sigma(i,j,ispec)   = 1.!0.4013
   
   potential_dphi_dx_DG(ibool(i,j,ispec)) = 0.
   potential_dphi_dz_DG(ibool(i,j,ispec)) = gravityext(i,j,ispec)
   endif
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! If not an external model => means that no initial conditions specified
  else
  
  ! rely on Mars model (large) published in SSR et 1km altitude
  ! AGW_FD_2016_eos/DATA_FD/atmos_model/MARS_MCD_AGW_model_INSIGHT_LS_240_LT_12h.txt
  ! and MARS_MCD_AGW_model_INSIGHT_LS_90_LT_12h_large.dat
  ! At initial time if no external model we create an uniform gravity fields
  if(timelocal == 0) then
        
        gravityext(i,j,ispec) = real(gravity_cte_DG, kind=CUSTOM_REAL)!9.81!gravity_cte_DG!0.0 !3.7247 real(coord(1,ibool_before_perio(i,j,ispec)), kind=CUSTOM_REAL)
        !WRITE(*,*) ">>>", gravity_cte_DG, gravityext(i,j,ispec), abs(real(gravity_cte_DG, kind=CUSTOM_REAL)- gravityext(i,j,ispec))
        !cp = 29 !7!1010!7/2
        !cnu = 20.8 !5!718!5/2
        gammaext_DG(ibool_DG(i,j,ispec))   = cp/cnu!1.33076167! 1.4!gamma_euler
        mu_visco  = dynamic_viscosity_cte_DG!1.092656e-05
        eta_visco = 2.5*dynamic_viscosity_cte_DG!2.67e-05
        muext(i,j,ispec)  = mu_visco
        etaext(i,j,ispec) = eta_visco
!        kappa_DG(i,j,ispec) = 0.0
        kappa_DG(i,j,ispec) = thermal_conductivity_cte_DG!0.!0.05!4.79046750E-04
       tau_epsilon(i,j,ispec) = tau_eps_cte_DG!1 !1.5!2.!1.5
        tau_sigma(i,j,ispec)   = tau_sig_cte_DG!1 !1.!1./(8*PI**2)!0.025!1
        
        potential_dphi_dx_DG(ibool(i,j,ispec)) = 0.
        potential_dphi_dz_DG(ibool(i,j,ispec)) = gravityext(i,j,ispec)
  endif
  
  e1_DG_P = 0.
  
   ! Acoustic only
   rho_DG_P = (1.164)
   ! Gravi
   if(USE_ISOTHERMAL_MODEL) then
           H = SCALE_HEIGHT!9350!10000
           G = gravityext(i,j,ispec)
           rho_DG_P = (1.164)*exp(-z/H)
   endif
   
   ! Acoustic only
   p_DG_P = (234.0151**2)*rho_DG_P/gammaext_DG(ibool_DG(i,j,ispec))
   ! Gravi
   if(USE_ISOTHERMAL_MODEL) p_DG_P = rho_DG_P*G*H
   
   veloc_x_DG_P = ZEROl
   veloc_z_DG_P = ZEROl
   
   endif
   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Rayleigh-Taylor instability (case 5.3 paper)
   if(.false.) then
   Tl = 2.!TWOl
   Tu = 1.
   p0   = ONEl
   rho0 = ONEl
   RR   = 1.!8.3145_CUSTOM_REAL
   gammaext_DG(ibool_DG(i,j,ispec))   = 1.4
   gravityext(i,j,ispec) = 1
   !if(z <= 1._CUSTOM_REAL) &
   p_DG_P = p0*exp(-(z-ONEl/(TWOl))/(RR*Tl))
   if(z > ONEl/TWOl) &
   p_DG_P = p0*exp(-(z-ONEl/(TWOl))/(RR*Tu))
  !!if(z <= 1._CUSTOM_REAL) &
   rho_DG_P = p_DG_p/(RR*Tl)
   if(z > ONEl/TWOl) &
   rho_DG_P = p_DG_p/(RR*Tu)
   veloc_x_DG_P = ZEROl
   veloc_z_DG_P = ZEROl
   endif
   
   !call boundary_forcing_DG(timelocal, rho_DG_P, veloc_x_DG_P, veloc_z_DG_P, E_DG_P)
   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Acoustic plane wave forcing
   to = 100!5.!/5.
   perio = 100!4.!/5.
   if(z == ZEROl .AND. .false.) &     
   veloc_z_DG_P = 0.01*(&
                  - (2d0/(perio/4d0))*((timelocal-(to-perio/4d0))/(perio/4d0))* &
                           (exp(-((timelocal-(to-perio/4d0))/(perio/4d0))**2)) &
                  + (2d0/(perio/4d0))*((timelocal-(to+perio/4d0))/(perio/4d0))* &
                           (exp(-((timelocal-(to+perio/4d0))/(perio/4d0))**2)) ) 
   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Acoustic sinus plane wave forcing
   to = 5!5.!/5.
   perio = 150!4.!/5.
   if(z == ZEROl .AND. .false.) then
   veloc_z_DG_P = 0.
   if(timelocal - to >= 0. .AND. timelocal - to <= perio ) &
   veloc_z_DG_P = 0.1*sin((timelocal-to)*2*PI/perio)
   endif
   
   !!!!!!!!!!!!!!!!!!!!!!!!
   ! Gravi wave forcing
   lambdo = 100000
   perio  = 40!200
   xo = 300000
   to = 50!35!250.
   if(z == ZEROl .AND. .false.) &     
   veloc_z_DG_P = 0.001*(&
                  - (2d0/(perio/4d0))*((timelocal-(to-perio/4d0))/(perio/4d0))* &
                           (exp(-((timelocal-(to-perio/4d0))/(perio/4d0))**2)) &
                  + (2d0/(perio/4d0))*((timelocal-(to+perio/4d0))/(perio/4d0))* &
                           (exp(-((timelocal-(to+perio/4d0))/(perio/4d0))**2)) ) &
                   * ( exp(-((x-(xo-lambdo/4))/(lambdo/4))**2) - &
                          exp(-((x-(xo+lambdo/4))/(lambdo/4))**2) ) 
   
   E_DG_P       = p_DG_P/(gammaext_DG(ibool_DG(i,j,ispec)) - ONEl) &
                        + rho_DG_P*HALFl*( veloc_x_DG_P**2 + veloc_z_DG_P**2 ) 
   
   rhovx_DG_P   = rho_DG_P*veloc_x_DG_P
   rhovz_DG_P   = rho_DG_P*veloc_z_DG_P
   
  end subroutine boundary_condition_DG
 
  subroutine compute_interface_unknowns(i, j, ispec, rho_DG_P, rhovx_DG_P, &
                rhovz_DG_P, E_DG_P, veloc_x_DG_P, veloc_z_DG_P, p_DG_P, T_P, &
                Tx_DG_P, Tz_DG_P, Vxx_DG_P, Vzz_DG_P, Vzx_DG_P, Vxz_DG_P, gamma_P, &
                neighbor, MPI_change_cpt, &
                exact_interface_flux, &
                rho_DG_iM, E_DG_iM, rhovx_DG_iM, rhovz_DG_iM, &
                V_DG_iM, T_DG_iM, &
                rho_DG_iP, E_DG_iP, rhovx_DG_iP, rhovz_DG_iP, &
                V_DG_iP, T_DG_iP, &
                MPI_iglob, chosen_nxnz_forMPI, dir_normal, nx, nz, weight, timelocal, &
                elastic_tensor)
  
! compute forces in the acoustic elements in forward simulation and in adjoint simulation in adjoint inversion
  use constants,only: CUSTOM_REAL,NGLLX,NGLLZ,gamma_euler

  use specfem_par,only: nglob_DG,nspec, &!ispec_is_acoustic, &
                         normal_DG, normal_DG_corner, ibool_DG, weight_DG, weight_DG_corner, &
                         ispec_is_acoustic_forcing, &
                         ACOUSTIC_FORCING, is_corner, &
                          ispec_is_acoustic_coupling_el, veloc_elastic,&
                         dir_normal_DG, dir_normal_DG_corner, &
                         DIR_RIGHT, DIR_LEFT, DIR_UP, DIR_DOWN, &
                         buffer_DG_rho_P, buffer_DG_rhovx_P, buffer_DG_rhovz_P, buffer_DG_E_P, NPROC, &
                         buffer_DG_Vxx_P, buffer_DG_Vzz_P, buffer_DG_Vxz_P, buffer_DG_Vzx_P, buffer_DG_Tz_P, buffer_DG_Tx_P, &
                         MPI_transfer, p_DG_init, gammaext_DG, muext, etaext, kappa_DG, ibool, cnu, &
                         ! TEST
                         buffer_DG_gamma_P, myrank, coord, it, i_stage, &
                         rho_init, rhovx_init, rhovz_init, E_init!, ibool
                         
  implicit none
  
  real(kind=CUSTOM_REAL), dimension(NGLLX, NGLLZ, nspec, 4) :: elastic_tensor
  
  integer, intent(in) :: i, j, ispec, chosen_nxnz_forMPI
  
  integer, dimension(3), intent(in) :: neighbor
  
  real(kind=CUSTOM_REAL), intent(out) :: rho_DG_P, rhovx_DG_P, rhovz_DG_P, &
        E_DG_P, veloc_x_DG_P, veloc_z_DG_P, p_DG_P, T_P, &
        Tx_DG_P, Tz_DG_P, Vxx_DG_P, Vzz_DG_P, Vzx_DG_P, Vxz_DG_P
        
  real(kind=CUSTOM_REAL), intent(in) :: timelocal
        
  logical, intent(inout) :: exact_interface_flux, MPI_change_cpt
  integer, dimension(nglob_DG), intent(inout) :: MPI_iglob
  integer, intent(inout) :: dir_normal
  
  real(kind=CUSTOM_REAL) :: tx, tz, normal_v, tangential_v, &
        nx, nz, veloc_x, veloc_z, weight, &
        veloc_x_DG_iM, veloc_z_DG_iM, p_DG_iM, gamma_P, e1_DG_P
        
  real(kind=CUSTOM_REAL), intent(in) :: rho_DG_iM, E_DG_iM, rhovx_DG_iM, rhovz_DG_iM, &
                rho_DG_iP, E_DG_iP, rhovx_DG_iP, rhovz_DG_iP
  real(kind=CUSTOM_REAL), dimension(2), intent(in) :: T_DG_iP, T_DG_iM
  real(kind=CUSTOM_REAL), dimension(2,2), intent(in) :: V_DG_iP, V_DG_iM
        
  ! Local variables     
  integer :: iglobM, i_el, j_el, ispec_el, iglob, iglobP, ipoin, num_interface
  real(kind=CUSTOM_REAL), dimension(2,2) :: trans_boundary, tensor_temp
  
  ! Parameters
  real(kind=CUSTOM_REAL), parameter :: ZERO = 0._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: ONE  = 1._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: TWO  = 2._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: HALF = 0.5_CUSTOM_REAL
  
  ! Characteristic based BC
  real(kind=CUSTOM_REAL) :: p_b, rho_b, s_b, rho_inf, v_b_x, v_b_z, c_b, un_b, &
  rlambda_max, rlambda_min, c_in, c_inf, un_in, un_inf, &
  deltaZ1, deltaZ2star, p_n, a_n, alpha0
  
  !character(len=100) file_name
  
  iglobM = ibool_DG(i,j,ispec)
  
  veloc_x_DG_iM = rhovx_DG_iM/rho_DG_iM
  veloc_z_DG_iM = rhovz_DG_iM/rho_DG_iM
  p_DG_iM       = (gammaext_DG(iglobM) - ONE)*( E_DG_iM &
        - (HALF)*rho_DG_iM*( veloc_x_DG_iM**2 + veloc_z_DG_iM**2 ) )
  
  Vxx_DG_P = -V_DG_iM(1,1)
  Vzz_DG_P = -V_DG_iM(2,2)
  Vxz_DG_P = -V_DG_iM(1,2)
  Vzx_DG_P = -V_DG_iM(2,1)
  Tx_DG_P = -T_DG_iM(1)
  Tz_DG_P = -T_DG_iM(2)
  
  gamma_P = gammaext_DG(iglobM)
  
  e1_DG_P = 0
  
  exact_interface_flux = .false.
  MPI_change_cpt = .false.
  if(neighbor(3) == -1) then
        
        ! MPI neighbor
        ipoin         = -1
        num_interface = -1
        if(NPROC > 1) then
        ipoin         = MPI_transfer(iglobM,MPI_iglob(iglobM),1)
        num_interface = MPI_transfer(iglobM,MPI_iglob(iglobM),2)
        endif                  
        ! If MPI neighbor, but not diagonal corner element and not outside element
        if(ipoin > -1) then
                        
                        ! Check for mistakes at corners
                        if(is_corner(i,j) .AND. &
                                ( (MPI_transfer(iglobM,MPI_iglob(iglobM),3) == i .AND. &
                                        (dir_normal == DIR_LEFT .OR. dir_normal == DIR_RIGHT) ) .OR. &
                                  (MPI_transfer(iglobM,MPI_iglob(iglobM),4) == j .AND. &
                                        (dir_normal == DIR_UP .OR. dir_normal == DIR_DOWN) ) ) ) then
                        
                                MPI_change_cpt = .true.
                        
                                if(chosen_nxnz_forMPI == 1) then
                                nx     = normal_DG(i,j,ispec, 1)
                                nz     = normal_DG(i,j,ispec, 2)
                                weight = weight_DG(i, j, ispec)
                                dir_normal = dir_normal_DG(i,j,ispec)
                                elseif(chosen_nxnz_forMPI == 0) then
                                nx     = normal_DG_corner(i,j,ispec, 1)
                                nz     = normal_DG_corner(i,j,ispec, 2)
                                weight = weight_DG_corner(i, j, ispec)
                                dir_normal = dir_normal_DG_corner(i,j,ispec)
                                endif
                                
                        endif
                        
                        MPI_iglob(iglobM) = MPI_iglob(iglobM) + 1
                        
                        rho_DG_P     = buffer_DG_rho_P(ipoin,num_interface)
                        E_DG_P       = buffer_DG_E_P(ipoin,num_interface)
                        rhovx_DG_P   = buffer_DG_rhovx_P(ipoin,num_interface)
                        rhovz_DG_P   = buffer_DG_rhovz_P(ipoin,num_interface)
                        veloc_x_DG_P = rhovx_DG_P/rho_DG_P
                        veloc_z_DG_P = rhovz_DG_P/rho_DG_P
                        
                        gamma_P   = buffer_DG_gamma_P(ipoin,num_interface)
                        
                        p_DG_P       = (gamma_P - ONE)*( E_DG_P &
                                - (HALF)*rho_DG_P*( veloc_x_DG_P**2 + veloc_z_DG_P**2 ) )
                              
                        if(muext(i,j,ispec) > 0 .OR. etaext(i,j,ispec) > 0 &
                        .OR. kappa_DG(i,j,ispec) > 0) then
                        ! Viscosity  
                        Vxx_DG_P = buffer_DG_Vxx_P(ipoin,num_interface)
                        Vzz_DG_P = buffer_DG_Vzz_P(ipoin,num_interface)
                        Vxz_DG_P = buffer_DG_Vxz_P(ipoin,num_interface)
                        Vzx_DG_P = buffer_DG_Vzx_P(ipoin,num_interface)
                        Tx_DG_P = buffer_DG_Tx_P(ipoin,num_interface)
                        Tz_DG_P = buffer_DG_Tz_P(ipoin,num_interface)
                        endif
                        
                        T_P = (E_DG_P/rho_DG_P - 0.5*((rhovx_DG_P/rho_DG_P)**2 + (rhovz_DG_P/rho_DG_P)**2))/(cnu)
                        
        elseif(ACOUSTIC_FORCING .AND. ispec_is_acoustic_forcing(i,j,ispec)) then
        
                stop 'ACOUSTIC_FORCING obsolete for DG simulations'
                
        ! Elastic coupling        
        elseif(ispec_is_acoustic_coupling_el(i,j,ispec,3) >= 0) then
        
               ! If we already know the "real" flux at boundary
               exact_interface_flux = .true.
        
               ! Coordinates of elastic element
               i_el     = ispec_is_acoustic_coupling_el(i,j,ispec,1)
               j_el     = ispec_is_acoustic_coupling_el(i,j,ispec,2)
               ispec_el = ispec_is_acoustic_coupling_el(i,j,ispec,3)
               
               iglob = ibool(i_el,j_el,ispec_el)
        
               ! Only for density
               call boundary_condition_DG(i, j, ispec, timelocal, rho_DG_P, rhovx_DG_P, rhovz_DG_P, E_DG_P, &
                veloc_x_DG_P, veloc_z_DG_P, p_DG_P, e1_DG_P)
        
               !rho_DG_P = rho_DG_iM
               
               ! Elastic velocities
               veloc_x = veloc_elastic(1,iglob)
               veloc_z = veloc_elastic(2,iglob)
               
               ! Tangential vector
               ! Since only bottom topography nz > 0
               tx = -nz
               tz = nx
               
               ! Normal velocity of the solid perturbation and tangential velocity of the fluid flow
               normal_v     = (veloc_x*nx + veloc_z*nz) 
               !normal_v     = -(veloc_x_DG_iM*nx + veloc_z_DG_iM*nz) + 2*(veloc_x*nx + veloc_z*nz) 
               tangential_v = veloc_x_DG_P*tx + veloc_z_DG_P*tz
               !tangential_v = veloc_x_DG_iM*tx + veloc_z_DG_iM*tz
               !tangential_v    = (veloc_x*tx + veloc_z*tz) 
               
               ! Transformation matrix between mesh coordinates and normal/tangential coordinates
               trans_boundary(1,1) = tz
               trans_boundary(1,2) = -nz
               trans_boundary(2,1) = -tx
               trans_boundary(2,2) = nx
               trans_boundary = trans_boundary/(nx*tz - tx*nz)
               
               ! From free slip and normal velocity continuity
               veloc_x_DG_P = trans_boundary(1,1)*normal_v + trans_boundary(1,2)*tangential_v!veloc_elastic(1,iglob)
               veloc_z_DG_P = trans_boundary(2,1)*normal_v + trans_boundary(2,2)*tangential_v
               
               !veloc_x_DG_P = veloc_x 
               !veloc_z_DG_P = veloc_z 
               
               !veloc_z_DG_P = veloc_z
               !if(abs(veloc_z_DG_P) > 0) &
               !veloc_x_DG_P = -veloc_x_DG_iM!-elastic_tensor(i_el,j_el,ispec_el,2)/(rho_DG_P*veloc_z_DG_P)
                        !( -(elastic_tensor(i_el,j_el,ispec_el,1)*nx + elastic_tensor(i_el,j_el,ispec_el,2)*nz)*tx &
                        !-(elastic_tensor(i_el,j_el,ispec_el,3)*nx + elastic_tensor(i_el,j_el,ispec_el,4)*nz)*tz ) &
                        !/(rho_DG_P*veloc_z_DG_P*nz*tx)
               
               !tensor_temp = 0
               ! Maybe we should remove the initial velocity field to veloc_x_DG_P/veloc_z_DG_P in the following
               ! Since we should only have the perturbated non-linear stress tensor
               tensor_temp(1,1) = -elastic_tensor(i_el,j_el,ispec_el,1) - 0*rho_DG_P*veloc_x_DG_P**2 
               tensor_temp(1,2) = -elastic_tensor(i_el,j_el,ispec_el,2) - 0*rho_DG_P*veloc_x_DG_P*veloc_z_DG_P
               tensor_temp(2,1) = -elastic_tensor(i_el,j_el,ispec_el,3) - 0*rho_DG_P*veloc_x_DG_P*veloc_z_DG_P
               tensor_temp(2,2) = -elastic_tensor(i_el,j_el,ispec_el,4) - 0*rho_DG_P*veloc_z_DG_P**2 
               
               ! From traction continuity
               p_DG_P = p_DG_init(iglobM) - (&
                  nx*( nx*tensor_temp(1,1) + nz*tensor_temp(1,2) ) &
                + nz*( nx*tensor_temp(2,1) + nz*tensor_temp(2,2) ) )!* 2 - (p_DG_init(iglobM) - p_DG_iM)
               
               E_DG_P       = p_DG_P/(gammaext_DG(iglobM) - 1.) + HALF*rho_DG_P*( veloc_x_DG_P**2 + veloc_z_DG_P**2 )
               
               !p_DG_P = 2*p_DG_P - p_DG_iM
               !veloc_x_DG_P = 2*veloc_x -veloc_x_DG_iM
               !veloc_z_DG_P = 2*veloc_z -veloc_z_DG_iM
               !E_DG_P = 2*E_DG_P - E_DG_iM
               
               if( abs(p_DG_P) > 1.2*abs(p_DG_init(iglobM)) .OR. abs(p_DG_P) < 0.8*abs(p_DG_init(iglobM)) ) &
               WRITE(*,*) it,i,j,ispec,">>", p_DG_P, p_DG_init(iglobM), (&
                  nx*( nx*tensor_temp(1,1) + nz*tensor_temp(1,2) ) &
                + nz*( nx*tensor_temp(2,1) + nz*tensor_temp(2,2) ) ), rho_DG_P*veloc_x_DG_P**2 , rho_DG_P*veloc_z_DG_P**2, &
                        rho_DG_P*veloc_x_DG_P*veloc_z_DG_P, tensor_temp(1,1), tensor_temp(1,2), tensor_temp(2,1), tensor_temp(2,2)
               !if(abs(elastic_tensor(i_el,j_el,ispec_el,2)*nz- elastic_tensor(i_el,j_el,ispec_el,4)*nz) > 0. ) &
               !WRITE(*,*) ">>>", i,j,ispec, elastic_tensor(i_el,j_el,ispec_el,2)*nz, elastic_tensor(i_el,j_el,ispec_el,4)*nz
                
!               p_DG_P = p_DG_init(iglobM) - (nz*tensor_temp(1,2) + nz*tensor_temp(2,2))/((nx+nz))!p_DG_iM
               
               !E_DG_P       = p_DG_P/(gammaext_DG(iglobM) - 1.) + HALF*rho_DG_P*( veloc_x_DG_P**2 + veloc_z_DG_P**2 )
               
               !veloc_x_DG_P = 2*veloc_x_DG_P - veloc_x_DG_iM
               !veloc_z_DG_P = 2*veloc_z_DG_P - veloc_z_DG_iM
               !p_DG_P       = 2*p_DG_P - p_DG_iM
               !E_DG_P       = 2*E_DG_P - E_DG_iM
               
               rhovx_DG_P   = rho_DG_P*veloc_x_DG_P
               rhovz_DG_P   = rho_DG_P*veloc_z_DG_P
               
               T_P = (E_DG_iM/rho_DG_iM - 0.5*(veloc_x_DG_iM**2 + veloc_z_DG_iM**2))/(cnu)
               !T_P = (E_DG_P/rho_DG_P - 0.5*((rhovx_DG_P/rho_DG_P)**2 + (rhovz_DG_P/rho_DG_P)**2))/(cnu)
        ! Classical boundary conditions
        else
        
                exact_interface_flux = .false.
        
                call boundary_condition_DG(i, j, ispec, timelocal, rho_DG_P, rhovx_DG_P, rhovz_DG_P, E_DG_P, &
                veloc_x_DG_P, veloc_z_DG_P, p_DG_P, e1_DG_P)
                
                !rho_DG_P = rho_DG_iM
                
                tx = -nz
                tz = nx
               ! Normal velocity of the solid perturbation and tangential velocity of the fluid flow
                normal_v     = veloc_x_DG_P*nx + veloc_z_DG_P*nz
                !normal_v     = 2*normal_v-(veloc_x_DG_iM*nx + veloc_z_DG_iM*nz)
                !tangential_v = -(veloc_x_DG_iM*tx + veloc_z_DG_iM*tz)
                tangential_v = (veloc_x_DG_P*tx + veloc_z_DG_P*tz)
                
               ! Transformation matrix between mesh coordinates and normal/tangential coordinates
                trans_boundary(1,1) = tz
                trans_boundary(1,2) = -nz
                trans_boundary(2,1) = -tx
                trans_boundary(2,2) = nx
                trans_boundary = trans_boundary/(nx*tz - tx*nz)
               
               ! From free slip and normal velocity continuity
               veloc_x_DG_P = trans_boundary(1,1)*normal_v + trans_boundary(1,2)*tangential_v!veloc_elastic(1,iglob)
               veloc_z_DG_P = trans_boundary(2,1)*normal_v + trans_boundary(2,2)*tangential_v
                
                !if(coord(2,ibool(i,j,ispec)) == 0) &
                !WRITE(*,*) ">>>>", veloc_x_DG_P, veloc_z_DG_P, nx,nz, coord(:,ibool(i,j,ispec)), veloc_x_DG_P, veloc_z_DG_P, &
                !        veloc_x_DG_iM, veloc_z_DG_iM
                
                E_DG_P       = p_DG_P/(gammaext_DG(iglobM) - ONE) &
                        + rho_DG_P*HALF*( veloc_x_DG_P**2 + veloc_z_DG_P**2 ) 
                !E_DG_P       = p_DG_iM/(gammaext_DG(iglobM) - ONE)! &
                        !+ rho_DG_P*HALF*( veloc_x_DG_P**2 + veloc_z_DG_P**2 ) 
                !E_DG_P = E_DG_iM
                
                rhovx_DG_P   = rho_DG_P*veloc_x_DG_P
                rhovz_DG_P   = rho_DG_P*veloc_z_DG_P
                
                if(coord(2,ibool(i,j,ispec)) == 20000. .AND. .false.) then
                
                ! Characteristic based BC
                un_in       = veloc_x_DG_iM*nx + veloc_z_DG_iM*nz
                un_inf      = veloc_x_DG_P*nx  + veloc_z_DG_P*nz
                c_in        = sqrt(gammaext_DG(iglobM)*p_DG_iM/rho_DG_iM)
                c_inf       = sqrt(gammaext_DG(iglobM)*p_DG_P/rho_DG_P)
                
                rlambda_max = un_in  - 2.*c_in/(gammaext_DG(iglobM) - 1.)
                rlambda_min = un_inf + 2.*c_inf/(gammaext_DG(iglobM) - 1.)
                
                un_b = (rlambda_min + rlambda_max)/2.
                c_b  = (gammaext_DG(iglobM) - 1.)*(rlambda_min - rlambda_max)/4.
                
                ! Boudnary values for velocity
                v_b_x  = veloc_x_DG_P + ( un_b - un_inf )*nx
                v_b_z  = veloc_z_DG_P + ( un_b - un_inf )*nz
                
                rho_inf = rho_DG_P
                s_b = (c_in**2)/( gammaext_DG(iglobM)*(rho_DG_iM**(gammaext_DG(iglobM) - 1.)) )
                
                ! Boundary values for density and pressure
                rho_b = ((c_b**2)/(gammaext_DG(iglobM)*s_b))**(1/(gammaext_DG(iglobM)-1.))
                p_b    = rho_b*(c_b**2)/gammaext_DG(iglobM)
                
                !veloc_x_DG_P = v_b_x
                !veloc_z_DG_P = v_b_z
                !rho_DG_P     = rho_b
                !p_DG_P       = p_b
                
                !WRITE(*,*) it, "TEST", i, j, ispec, ">>", rho_b, v_b_x, v_b_z, p_b, rho_inf**(gammaext_DG(iglobM) - 1.), &
                !        c_inf, c_b, rlambda_min, rlambda_max, un_in, un_inf, c_in, c_inf, &
                !        veloc_x_DG_iM, veloc_z_DG_iM, veloc_x_DG_P, veloc_z_DG_P,s_b
                !WRITE(*,*) it,i_stage, veloc_x_DG_iM, veloc_z_DG_iM, veloc_x_DG_P, veloc_z_DG_P
                !stop 'TOTO'
                
                !E_DG_P       = p_DG_P/(gammaext_DG(iglobM) - ONE) &
                 !       + rho_DG_P*HALF*( veloc_x_DG_P**2 + veloc_z_DG_P**2 )
                
                !else
                
                rho_DG_P      = -rho_DG_iM + rho_b
                veloc_x_DG_P  = v_b_x!veloc_x_DG_P -2.*( normal_v + veloc_x_DG_P*nx + veloc_z_DG_P*nz )*nx
                veloc_z_DG_P  = v_b_z!veloc_z_DG_P -2.*( normal_v + veloc_x_DG_P*nx + veloc_z_DG_P*nz )*nz
                p_DG_P = p_b
                E_DG_P = p_DG_P/(gammaext_DG(iglobM) - ONE) &
                        + rho_DG_P*HALF*( veloc_x_DG_P**2 + veloc_z_DG_P**2 )
                
                rhovx_DG_P   = -rhovx_DG_iM + 2.*rho_DG_P*veloc_x_DG_P
                rhovz_DG_P   = -rhovz_DG_iM + 2.*rho_DG_P*veloc_z_DG_P
                rho_DG_P     = -rho_DG_iM + 2.*rho_DG_P
                !veloc_x_DG_P = -veloc_x_DG_iM + 2.*veloc_x_DG_P
                !veloc_z_DG_P = -veloc_z_DG_iM + 2.*veloc_z_DG_P     
                E_DG_P       = -E_DG_iM + 2.*E_DG_P
                   
                
                endif
                
                if(coord(2,ibool(i,j,ispec)) == 0. .AND. .false.) then
                        
                        p_n = (gammaext_DG(iglobM) - ONE)*( E_init(iglobM) &
                                - (HALF/rho_init(iglobM))*( rhovx_init(iglobM)**2 + rhovz_init(iglobM)**2 ) )
                        a_n = sqrt(gammaext_DG(iglobM)*p_n/rho_init(iglobM))
                        deltaZ1 = 2*p_DG_P - p_DG_iM - p_n &
                                - rho_init(iglobM)*a_n&
                                        *((veloc_x_DG_iM - rhovx_init(iglobM)/rho_init(iglobM))*nx &
                                                +(veloc_z_DG_iM-rhovz_init(iglobM)/rho_init(iglobM))*nz)
                
                        deltaZ2star = (p_DG_iM - p_n) &
                           + rho_init(iglobM)*a_n*((veloc_x_DG_iM - rhovx_init(iglobM)/rho_init(iglobM))*nx &
                                                 + (veloc_z_DG_iM-rhovz_init(iglobM)/rho_init(iglobM))*nz)
                
                       ! rho_DG_P = rho_DG_iM
                        
                        alpha0 = 0.0001          
                        
                        p_DG_P = p_n - (alpha0*deltaZ1 + deltaZ2star)/2
                        
                        normal_v     = ((rhovx_init(iglobM)/rho_init(iglobM))*nx+(rhovz_init(iglobM)/rho_init(iglobM))*nz) &
                                - (alpha0*deltaZ1 - deltaZ2star)/(2.*a_n*rho_init(iglobM))
                        tangential_v = 0.
                
                        ! Transformation matrix between mesh coordinates and normal/tangential coordinates
                        trans_boundary(1,1) = tz
                        trans_boundary(1,2) = -nz
                        trans_boundary(2,1) = -tx
                        trans_boundary(2,2) = nx
                        trans_boundary = trans_boundary/(nx*tz - tx*nz)
               
                        ! From free slip and normal velocity continuity
                        veloc_x_DG_P = trans_boundary(1,1)*normal_v + trans_boundary(1,2)*tangential_v!veloc_elastic(1,iglob)
                        veloc_z_DG_P = trans_boundary(2,1)*normal_v + trans_boundary(2,2)*tangential_v
                        
                        
                        
                endif
                
                !E_DG_P = p_DG_P/(gammaext_DG(iglobM) - ONE) &
                !        + rho_DG_P*HALF*( veloc_x_DG_P**2 + veloc_z_DG_P**2 )
                
                !T_P = (E_DG_iM/rho_DG_iM - 0.5*(veloc_x_DG_iM**2 + veloc_z_DG_iM**2))/(cnu)
                T_P = (E_DG_P/rho_DG_P - 0.5*(veloc_x_DG_P**2 + veloc_z_DG_P**2))/(cnu)
                !Tx_DG_P = -T_DG_iM(1)
                !Tz_DG_P = -T_DG_iM(2)
                
        endif
        
   ! Not an outside edge
   else
                iglobP       = ibool_DG(neighbor(1),neighbor(2),neighbor(3))
                
                gamma_P = gammaext_DG(iglobP)
                
                rho_DG_P     = rho_DG_iP
                E_DG_P       = E_DG_iP
                rhovx_DG_P   = rhovx_DG_iP
                rhovz_DG_P   = rhovz_DG_iP
                
                veloc_x_DG_P = rhovx_DG_P/rho_DG_P
                veloc_z_DG_P = rhovz_DG_P/rho_DG_P
                p_DG_P       = (gamma_P - ONE)*( E_DG_P &
                - (HALF)*rho_DG_P*( veloc_x_DG_P**2 + veloc_z_DG_P**2 ) )
                
                if(muext(i,j,ispec) > 0 .OR. etaext(i,j,ispec) > 0 &
                        .OR. kappa_DG(i,j,ispec) > 0) then
                        ! Viscosity  
                        Vxx_DG_P = V_DG_iP(1,1)!,iglobP)
                        Vzz_DG_P = V_DG_iP(2,2)!,iglobP)
                        Vxz_DG_P = V_DG_iP(1,2)!,iglobP)
                        Vzx_DG_P = V_DG_iP(2,1)!,iglobP)
                        Tx_DG_P = T_DG_iP(1)!,iglobP)
                        Tz_DG_P = T_DG_iP(2)!,iglobP)
                endif

                T_P = (E_DG_P/rho_DG_P - 0.5*((rhovx_DG_P/rho_DG_P)**2 + (rhovz_DG_P/rho_DG_P)**2))/(cnu)

  endif
  
  ! Temperature computation
  !T_P = (E_DG_P/rho_DG_P - 0.5*((rhovx_DG_P/rho_DG_P)**2 + (rhovz_DG_P/rho_DG_P)**2))/(cnu)
  ! (E_DG_iM/rho_DG_iM - 0.5*(veloc_x_DG_iM**2 + veloc_z_DG_iM**2))/(cnu)
  
  end subroutine compute_interface_unknowns
  

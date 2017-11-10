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

! ------------------------------------------------------------ !
! compute_forces_acoustic_DG                                   !
! ------------------------------------------------------------ !

  subroutine compute_forces_acoustic_DG(rho_DG_main, rhovx_DG_main, rhovz_DG_main, E_DG_main, &
        T_DG_main, V_DG_main, e1_DG, &
        dot_rho, dot_rhovx, dot_rhovz, dot_E, dot_e1, timelocal)


! compute forces in the acoustic elements in forward simulation and in adjoint simulation in adjoint inversion

  use constants,only: CUSTOM_REAL,NGLLX,NGLLZ,gamma_euler

  use specfem_par, only: nglob_DG,nspec, ispec_is_acoustic_DG,&!ispec_is_acoustic,MINMOD_FACTOR, USE_SLOPE_LIMITER
                         xix,xiz,gammax,gammaz,jacobian, &
                         hprimewgll_xx, &
                         hprimewgll_zz,wxgll,wzgll, &
                         normal_DG, normal_DG_corner, ibool_DG, weight_DG, weight_DG_corner, &
                         neighbor_DG, &
                         neighbor_DG_corner, is_corner, &
                         it,potential_dphi_dx_DG, potential_dphi_dz_DG, ibool, &
                         elastic_tensor, &
                         dir_normal_DG, dir_normal_DG_corner, &
                         DIR_RIGHT, DIR_LEFT, DIR_UP, DIR_DOWN, &
                         myrank, &
                         i_stage, p_DG_init, gammaext_DG, muext, etaext, kappa_DG,tau_epsilon, tau_sigma, &
                         !hprime_xx, hprime_zz,  cnu, &
                         !ibool_before_perio, coord, &
                         rhovx_init, rhovz_init, E_init, rho_init, &!T_init
                         CONSTRAIN_HYDROSTATIC, TYPE_SOURCE_DG
                         
  implicit none

  real(kind=CUSTOM_REAL), dimension(nglob_DG) :: rho_DG_main, rhovx_DG_main, rhovz_DG_main, E_DG_main, e1_DG
  real(kind=CUSTOM_REAL), dimension(2,nglob_DG) :: T_DG_main
  real(kind=CUSTOM_REAL), dimension(2,2,nglob_DG) :: V_DG_main
  real(kind=CUSTOM_REAL), dimension(nglob_DG) :: rho_DG, rhovx_DG, rhovz_DG, E_DG
  real(kind=CUSTOM_REAL), dimension(2,nglob_DG) :: T_DG
  real(kind=CUSTOM_REAL), dimension(2,2,nglob_DG) :: V_DG
  real(kind=CUSTOM_REAL), dimension(nglob_DG) :: dot_rho, dot_rhovx, dot_rhovz, dot_E, dot_e1
  
  ! local parameters
  integer :: ispec,i,j,k,iglob, it_corner
  integer :: ifirstelem,ilastelem

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: temp_rho_1, temp_rho_2, &
        temp_rhovx_1, temp_rhovx_2, temp_rhovz_1, temp_rhovz_2, &
        temp_E_1, temp_E_2, &
        temp_rho_gravi, temp_rhovx_gravi, temp_rhovz_gravi, temp_E_gravi!, temp_e1

  ! Jacobian matrix and determinant
  real(kind=CUSTOM_REAL) :: xixl,xizl,gammaxl,gammazl,jacobianl

  real(kind=CUSTOM_REAL) :: lambda, nx, nz, weight, &
        temp_unknown_M, temp_unknown_P, temp_unknown2_M, temp_unknown2_P, &
        temp_unknown, temp_unknown2, &
        flux_x, flux_z, flux_n, jump, &
        rho_DG_P, veloc_x_DG_P, veloc_z_DG_P, &
        E_DG_P, p_DG_P, rhovx_DG_P, rhovz_DG_P, timelocal, &
        Tx_DG_P, Tz_DG_P, Vxx_DG_P, Vzz_DG_P, Vxz_DG_P, Vzx_DG_P, T_P, &
        ! TEST
        gamma_P!, &
        !templ_rhovx_gravi, templ_rhovz_gravi, templ_rho_gravi, templ_E_gravi
    
  real(kind=CUSTOM_REAL) :: dT_dx, dT_dz, veloc_n_M, veloc_n_P
        
  integer :: iglobM, iglobP
  
  integer, dimension(3) :: neighbor
  
  ! Local
  real(kind=CUSTOM_REAL), dimension(nglob_DG) :: veloc_x_DG, veloc_z_DG, p_DG
  
  ! Viscosity
  !real(kind=CUSTOM_REAL) :: dux_dxi, dux_dgamma, duz_dxi, duz_dgamma
  real(kind=CUSTOM_REAL) :: dux_dx, dux_dz, duz_dx, duz_dz
  real(kind=CUSTOM_REAL) :: wxl, wzl
  
  ! Parameters
  real(kind=CUSTOM_REAL), parameter :: ZERO = 0._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: ONE  = 1._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: TWO  = 2._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: HALF = 0.5_CUSTOM_REAL
  
  real(kind=CUSTOM_REAL), parameter :: threshold = 0.000000001_CUSTOM_REAL
  
  ! Temporary
  !real(kind=CUSTOM_REAL), dimension(nglob_DG) :: rho0_DG, E0_DG, p0_DG
  logical, parameter :: ONLY_PERTURBATION = .false.
  
  logical :: neighbor_exists
  integer :: dir_normal, i_ex, j_ex, ispec_ex
  
  ! TEmporary test
  integer :: chosen_nxnz_forMPI
  !character(len=100) file_name
  
  logical :: MPI_change_cpt, exact_interface_flux
  
  integer, dimension(nglob_DG) :: MPI_iglob
  
  real(kind=CUSTOM_REAL) :: flux_rho, flux_rhovx, flux_rhovz, flux_E
  
  ! TEST ABSORB
  real(kind=CUSTOM_REAL) :: maxval_rho,maxval_rhovx,maxval_rhovz,maxval_E
  logical :: ABSORB_BC
  
  ! for MPI transfer
  MPI_iglob = 1
  
  !found = .false.
  
  ifirstelem = 1
  ilastelem = nspec
  
  rho_DG   = rho_DG_main
  rhovx_DG = rhovx_DG_main
  rhovz_DG = rhovz_DG_main
  E_DG     = E_DG_main
  
  ABSORB_BC = .false.
  if(ABSORB_BC) then
  maxval_rho   = maxval(rho_DG-rho_init)
  maxval_rhovx = maxval(rhovx_DG-rhovx_init)
  maxval_rhovz = maxval(rhovz_DG-rhovz_init)
  maxval_E     = maxval(E_DG-E_init)
  endif
  
  T_DG = T_DG_main
  V_DG = V_DG_main
  
  ! Initialise auxiliary unknwons.
  veloc_x_DG = rhovx_DG/rho_DG
  veloc_z_DG = rhovz_DG/rho_DG
  p_DG       = (gammaext_DG - ONE)*( E_DG &
        - (HALF)*rho_DG*( veloc_x_DG**2 + veloc_z_DG**2 ) )
        
  ! Initialisation.
  dot_rho   = ZERO
  dot_rhovx = ZERO
  dot_rhovz = ZERO
  dot_E     = ZERO
  dot_e1    = ZERO
  
  ! Add force source.
  if(TYPE_SOURCE_DG == 1) then
    call compute_add_sources_acoustic_DG_spread(dot_rho, it, i_stage)   
  elseif(TYPE_SOURCE_DG == 2) then
    call compute_add_sources_acoustic_DG_spread(dot_rhovx, it, i_stage)
    call compute_add_sources_acoustic_DG_spread(dot_rhovz, it, i_stage)
  elseif(TYPE_SOURCE_DG == 3) then
    call compute_add_sources_acoustic_DG_spread(dot_E, it, i_stage)
  endif
  
  ! TODO: introduce a verbosity parameter in order to prevent unwanted flooding of the terminal.
  if(myrank == 0 .AND. mod(it, 50) == 0) then
    write(*,"(a)") "               | max                     | min"
    WRITE(*,"(a,e24.16,a,e24.16)") " rho           |", maxval(rho_DG), " |", minval(rho_DG)
    WRITE(*,"(a,e24.16,a,e24.16)") " rhovx         |", maxval(rhovx_DG), " |", minval(rhovx_DG)
    WRITE(*,"(a,e24.16,a,e24.16)") " rhovz         |", maxval(rhovz_DG), " |", minval(rhovz_DG)
    WRITE(*,"(a,e24.16,a,e24.16)") " E             |", maxval(E_DG), " |", minval(E_DG)
    WRITE(*,"(a,e23.16,a)") "Ratio |p-p_{init}|/p_{init}:", maxval(abs((p_DG-p_DG_init)/p_DG_init)), "."
  endif
  
  do ispec = ifirstelem, ilastelem ! Loop over elements.
    ! acoustic spectral element
    if (ispec_is_acoustic_DG(ispec)) then !if (ispec_is_acoustic(ispec)) then
      ! --------------------------- !
      ! First set of loops: compute !
      ! volumic contributions.      !
      ! --------------------------- !
      do j = 1, NGLLZ
        do i = 1, NGLLX
          iglob = ibool_DG(i, j, ispec)
          
          jacobianl = jacobian(i, j, ispec)
        
          xixl = xix(i, j, ispec)
          xizl = xiz(i, j, ispec)
          gammaxl = gammax(i, j, ispec)
          gammazl = gammaz(i, j, ispec)
          
          wzl = wzgll(j)
          wxl = wxgll(i)
          
          ! Inviscid stress tensor's contributions.
          temp_unknown = rhovx_DG(iglob)
          temp_unknown2 = rhovz_DG(iglob)
          
          temp_rho_1(i,j) = wzl * jacobianl * (xixl * temp_unknown + xizl * temp_unknown2) 
          temp_rho_2(i,j) = wxl * jacobianl * (gammaxl * temp_unknown + gammazl * temp_unknown2) 
          
          if(.not. CONSTRAIN_HYDROSTATIC) then
            temp_unknown = rho_DG(iglob)*veloc_x_DG(iglob)**2 + p_DG(iglob)
          else
            temp_unknown = rho_DG(iglob)*veloc_x_DG(iglob)**2 + (p_DG(iglob) - p_DG_init(iglob))
          endif
          temp_unknown2 = rho_DG(iglob)*veloc_x_DG(iglob)*veloc_z_DG(iglob)
          
          temp_rhovx_1(i,j) = wzl * jacobianl * (xixl * temp_unknown + xizl * temp_unknown2) 
          temp_rhovx_2(i,j) = wxl * jacobianl * (gammaxl * temp_unknown + gammazl * temp_unknown2) 
          
          temp_unknown = rho_DG(iglob)*veloc_x_DG(iglob)*veloc_z_DG(iglob)
          if(.not. CONSTRAIN_HYDROSTATIC) then
            temp_unknown2 = rho_DG(iglob)*veloc_z_DG(iglob)**2 + p_DG(iglob)
          else
            temp_unknown2 = rho_DG(iglob)*veloc_z_DG(iglob)**2 + (p_DG(iglob) - p_DG_init(iglob))
          endif
          
          temp_rhovz_1(i,j) = wzl * jacobianl * (xixl * temp_unknown + xizl * temp_unknown2) 
          temp_rhovz_2(i,j) = wxl * jacobianl * (gammaxl * temp_unknown + gammazl * temp_unknown2) 
          
          if(.not. CONSTRAIN_HYDROSTATIC) then
            temp_unknown = veloc_x_DG(iglob)*(E_DG(iglob) + p_DG(iglob))
            temp_unknown2 = veloc_z_DG(iglob)*(E_DG(iglob) + p_DG(iglob))
          else          
            temp_unknown = veloc_x_DG(iglob)*(E_DG(iglob) + (p_DG(iglob) - p_DG_init(iglob)))
            temp_unknown2 = veloc_z_DG(iglob)*(E_DG(iglob) + (p_DG(iglob) - p_DG_init(iglob)))
          endif
          
          temp_E_1(i,j) = wzl * jacobianl * (xixl * temp_unknown + xizl * temp_unknown2) 
          temp_E_2(i,j) = wxl * jacobianl * (gammaxl * temp_unknown + gammazl * temp_unknown2) 
          
          ! Add the viscous stress tensor's contributions.
          dux_dx = V_DG(1, 1, iglob)
          dux_dz = V_DG(1, 2, iglob)
          duz_dx = V_DG(2, 1, iglob)
          duz_dz = V_DG(2, 2, iglob)
          if(muext(i, j, ispec) > 0 .OR. &
             etaext(i, j, ispec) > 0 .OR. &
             kappa_DG(i, j, ispec) > 0) then
          
            dT_dx  = T_DG(1, iglob)
            dT_dz  = T_DG(2, iglob)
            
            temp_unknown = muext(i, j, ispec)*TWO*dux_dx + (etaext(i, j, ispec) - (TWO/3.)*muext(i, j, ispec))*(dux_dx + duz_dz) 
            temp_unknown2 = muext(i, j, ispec)*( dux_dz + duz_dx )
            temp_rhovx_1(i,j) = temp_rhovx_1(i,j) - wzl * jacobianl * (xixl * temp_unknown + xizl * temp_unknown2) 
            temp_rhovx_2(i,j) = temp_rhovx_2(i,j) - wxl * jacobianl * (gammaxl * temp_unknown + gammazl * temp_unknown2) 
            
            temp_unknown  = veloc_x_DG(iglob)*temp_unknown + veloc_z_DG(iglob)*temp_unknown2
            temp_E_1(i,j) = temp_E_1(i,j) - wzl * jacobianl * (xixl * temp_unknown) 
            temp_E_2(i,j) = temp_E_2(i,j) - wxl * jacobianl * (gammaxl * temp_unknown) 
            
            temp_unknown = muext(i, j, ispec)*( dux_dz + duz_dx )
            temp_unknown2 = muext(i, j, ispec)*TWO*duz_dz + (etaext(i, j, ispec) - (TWO/3.)*muext(i, j, ispec))*(dux_dx + duz_dz) 
            temp_rhovz_1(i,j) = temp_rhovz_1(i,j) - wzl * jacobianl * (xixl * temp_unknown + xizl * temp_unknown2) 
            temp_rhovz_2(i,j) = temp_rhovz_2(i,j) - wxl * jacobianl * (gammaxl * temp_unknown + gammazl * temp_unknown2) 
            
            temp_unknown2  = veloc_x_DG(iglob)*temp_unknown + veloc_z_DG(iglob)*temp_unknown2
            
            temp_E_1(i,j) = temp_E_1(i,j) - wzl * jacobianl * (xizl * temp_unknown2) 
            temp_E_2(i,j) = temp_E_2(i,j) - wxl * jacobianl * (gammazl * temp_unknown2) 
            
            ! Add heat's contributions.
            temp_unknown  = kappa_DG(i, j, ispec)*dT_dx
            temp_unknown2 = kappa_DG(i, j, ispec)*dT_dz
            temp_E_1(i,j) = temp_E_1(i,j) - wzl * jacobianl * (xixl * temp_unknown + xizl * temp_unknown2) 
            temp_E_2(i,j) = temp_E_2(i,j) - wxl * jacobianl * (gammaxl * temp_unknown + gammazl * temp_unknown2) 
          endif
          
          
          ! Gravity contributions (separated from the rest).
          temp_rho_gravi(i,j) = 0.
          
          temp_rhovx_gravi(i,j) = -rho_DG(iglob)*potential_dphi_dx_DG(ibool(i, j, ispec))* jacobianl!
          
          if(.not. CONSTRAIN_HYDROSTATIC) then
            temp_rhovz_gravi(i,j) = -rho_DG(iglob)*potential_dphi_dz_DG(ibool(i, j, ispec))* jacobianl
          else
            temp_rhovz_gravi(i,j) = -(rho_DG(iglob) - rho_init(iglob)) * potential_dphi_dz_DG(ibool(i, j, ispec)) * jacobianl 
          endif
          
          if(.not. CONSTRAIN_HYDROSTATIC) then
            temp_E_gravi(i,j) = -rho_DG(iglob)*(veloc_x_DG(iglob)*potential_dphi_dx_DG(ibool(i, j, ispec)) + &
                                veloc_z_DG(iglob)*potential_dphi_dz_DG(ibool(i, j, ispec)))* jacobianl
          else
            temp_E_gravi(i,j) = &
                                -(rho_DG(iglob) - rho_init(iglob))*(veloc_x_DG(iglob)*potential_dphi_dx_DG(ibool(i, j, ispec)) + &
                                veloc_z_DG(iglob)*potential_dphi_dz_DG(ibool(i, j, ispec)))* jacobianl   
                        
            temp_E_gravi(i,j) = temp_E_gravi(i,j) - p_DG_init(iglob)*(dux_dx + duz_dz)* jacobianl       
          endif
                
          temp_E_gravi(i,j) = temp_E_gravi(i,j) - jacobianl * (p_DG_init(iglob)*gammaext_DG(iglob)) &
                              * ( (tau_epsilon(i, j, ispec)/tau_sigma(i, j, ispec)) - 1. ) &
                              * ( dux_dx + duz_dz - e1_DG(iglob))/(gammaext_DG(iglob) - ONE)
          
          
          dot_e1(iglob) = dot_e1(iglob) - (1/tau_sigma(i, j, ispec))*( &
                (1 - (tau_sigma(i, j, ispec)/tau_epsilon(i, j, ispec)))*(dux_dx + duz_dz) + e1_DG(iglob) )
          
        enddo
      enddo
      
      ! Assemble the contributions previously computed, and add gravity's contribution.
      do j = 1, NGLLZ
        do i = 1, NGLLX
          iglob = ibool_DG(i, j, ispec)
          ! along x direction and z direction
          do k = 1, NGLLX
            dot_rho(iglob) = dot_rho(iglob) + &
                             (temp_rho_1(k,j) * real(hprimewgll_xx(k,i), kind=CUSTOM_REAL) + &
                              temp_rho_2(i,k) * real(hprimewgll_zz(k,j), kind=CUSTOM_REAL))
                     
            dot_rhovx(iglob) = dot_rhovx(iglob) + &
                               (temp_rhovx_1(k,j) * real(hprimewgll_xx(k,i), kind=CUSTOM_REAL) + &
                                temp_rhovx_2(i,k) * real(hprimewgll_zz(k,j), kind=CUSTOM_REAL))
            
            dot_rhovz(iglob) = dot_rhovz(iglob) + &
                               (temp_rhovz_1(k,j) * real(hprimewgll_xx(k,i), kind=CUSTOM_REAL) + &
                                temp_rhovz_2(i,k) * real(hprimewgll_zz(k,j), kind=CUSTOM_REAL))
            
            dot_E(iglob) = dot_E(iglob) + &
                           (temp_E_1(k,j) * real(hprimewgll_xx(k,i), kind=CUSTOM_REAL) + &
                            temp_E_2(i,k) * real(hprimewgll_zz(k,j), kind=CUSTOM_REAL))
          enddo
          
          wzl = real(wzgll(j), kind=CUSTOM_REAL)
          wxl = real(wxgll(i), kind=CUSTOM_REAL)
          
          dot_rho(iglob)   = dot_rho(iglob)   + temp_rho_gravi(i, j) * wxl * wzl
          dot_rhovx(iglob) = dot_rhovx(iglob) + temp_rhovx_gravi(i, j) * wxl * wzl
          dot_rhovz(iglob) = dot_rhovz(iglob) + temp_rhovz_gravi(i, j) * wxl * wzl
          dot_E(iglob)     = dot_E(iglob)     + temp_E_gravi(i, j) * wxl * wzl
        enddo
      enddo
      
      it_corner = 0
      MPI_change_cpt = .false.
      
      ! --------------------------- !
      ! Second set of loops: add    !
      ! fluxes between elements,    !
      ! but only on exterior        !
      ! points.                     !
      ! --------------------------- !
      do  i = 1, NGLLX
        j = 1
        do while (j <= NGLLZ)
          ! Skip interior points.
          if(i > 1 .AND. i < NGLLX .AND. j > 1 .AND. j < NGLLZ) then
            j = j + 1
            cycle
          endif
          
          ! Step 1: prepare the normals' parameters (nx, nz, weight, dir_normal, etc.).
          
          ! Recover neighbor location
          neighbor = neighbor_DG(i, j, ispec, :)
          
          chosen_nxnz_forMPI = -1
        
          ! Reinit boolean to know if neighbor exists.
          neighbor_exists = .false.
          
          ! TODO: Maybe, put this in the else of the if(it_corner == 1) block below, in order to avoid double affectation.
          !nx = normal_DG(ispec, ind, 1)
          !nz = normal_DG(ispec, ind, 2)
          nx     = normal_DG(i, j, ispec, 1)
          nz     = normal_DG(i, j, ispec, 2)
          weight = weight_DG(i, j, ispec)
          dir_normal = dir_normal_DG(i, j, ispec)
          chosen_nxnz_forMPI = 0
          
          ! Needs x2 points at corners to correctly map edges
          ! => 2 step on the same point
          if(it_corner == 1) then
            neighbor = neighbor_DG_corner(i, j, ispec, :)
            nx     = normal_DG_corner(i, j, ispec, 1)
            nz     = normal_DG_corner(i, j, ispec, 2)
            weight = weight_DG_corner(i, j, ispec)
            dir_normal = dir_normal_DG_corner(i, j, ispec)
            chosen_nxnz_forMPI = 1
            it_corner = 2
          endif
          
          if(neighbor_DG(i, j, ispec, 3) > -1 .OR. &
             neighbor_DG_corner(i, j, ispec, 3) > -1) then
            neighbor_exists = .true.
          endif
          
          ! If not outer boundary check for corresponding neighbor normal
          if(is_corner(i, j)) then
            ! If at least one neighbor exists
            if(neighbor_exists) then
              i_ex = neighbor(1)
              j_ex = neighbor(2)
              ispec_ex = neighbor(3)
              ! If corner of an outside edge
              if(it_corner == 2 .AND. neighbor(3) == -1) then
                i_ex = neighbor_DG(i, j, ispec, 1)
                j_ex = neighbor_DG(i, j, ispec, 2)
                ispec_ex = neighbor_DG(i, j, ispec, 3)
              elseif(it_corner < 2 .AND. neighbor(3) == -1) then
                i_ex = neighbor_DG_corner(i, j, ispec, 1)
                j_ex = neighbor_DG_corner(i, j, ispec, 2)
                ispec_ex = neighbor_DG_corner(i, j, ispec, 3)
              endif
              
              ! Cross product to verify if the normal corresponds to the normal
              !normal_DG(i_ex,j_ex,ispec_ex, 1) 
              !normal_DG(i_ex,j_ex,ispec_ex, 2)
              if(dir_normal /= -dir_normal_DG(i_ex,j_ex,ispec_ex) .AND. &
                 dir_normal /= -dir_normal_DG_corner(i_ex,j_ex,ispec_ex) ) then
                ! Only change normal if inner element
                if(neighbor(3) > -1 .AND. it_corner < 2) then
                  nx     = normal_DG_corner(i, j, ispec, 1)
                  nz     = normal_DG_corner(i, j, ispec, 2)
                  weight = weight_DG_corner(i, j, ispec)
                  dir_normal = dir_normal_DG_corner(i, j, ispec)
                  ! MODIF for MPI
                  chosen_nxnz_forMPI = 1
                elseif(neighbor(3) > -1 .AND. it_corner == 2) then
                  nx     = normal_DG(i, j, ispec, 1)
                  nz     = normal_DG(i, j, ispec, 2)
                  weight = weight_DG(i, j, ispec)
                  dir_normal = dir_normal_DG(i, j, ispec)
                  ! MODIF for MPI
                  chosen_nxnz_forMPI = 0
                endif
              ! If outside element, if the normal corresponds to the one computed here
              ! it means that we should take the other one
              elseif(neighbor(3) == -1) then
                ! Only change normal if inner element
                if(it_corner < 2) then
                  nx     = normal_DG_corner(i, j, ispec, 1)
                  nz     = normal_DG_corner(i, j, ispec, 2)
                  weight = weight_DG_corner(i, j, ispec)
                  dir_normal = dir_normal_DG_corner(i, j, ispec)
                  ! MODIF for MPI
                  chosen_nxnz_forMPI = 1
                elseif(it_corner == 2) then
                  nx     = normal_DG(i, j, ispec, 1)
                  nz     = normal_DG(i, j, ispec, 2)
                  weight = weight_DG(i, j, ispec)
                  dir_normal = dir_normal_DG(i, j, ispec)
                  ! MODIF for MPI
                  chosen_nxnz_forMPI = 0
                endif
              endif
            endif
          endif ! End of if(is_corner(i, j)).
          
          ! TODO: Maybe, put all that follows (before "Step 2") inside the else block of the previous if(is_corner(i, j)) block.
          
          ! Interior point
          iglobM = ibool_DG(i, j, ispec)
          
          ! If a MPI surface node has been ill referenced and we need to switch between normal_DG and normal_DG_corner.
          if(MPI_change_cpt) then
            if(chosen_nxnz_forMPI == 1) then
              nx     = normal_DG(i, j, ispec, 1)
              nz     = normal_DG(i, j, ispec, 2)
              weight = weight_DG(i, j, ispec)
              dir_normal = dir_normal_DG(i, j, ispec)
            elseif(chosen_nxnz_forMPI == 0) then
              nx     = normal_DG_corner(i, j, ispec, 1)
              nz     = normal_DG_corner(i, j, ispec, 2)
              weight = weight_DG_corner(i, j, ispec)
              dir_normal = dir_normal_DG_corner(i, j, ispec)
            endif
          endif !if(is_MPI_interface_DG(iglobM) .AND. NPROC > 1)
          
          ! If at corner, then notify that we will need to go one more time.
          if(is_corner(i, j) .AND. it_corner == 0) then
            it_corner = 1
          endif
          
          ! Step 2: knowing the normals' parameters, compute now the fluxes. We use the Lax-Friedrich flux.
          
          rho_DG_P     = ZERO
          rhovx_DG_P   = ZERO
          rhovz_DG_P   = ZERO
          E_DG_P       = ZERO
          veloc_x_DG_P = ZERO
          veloc_z_DG_P = ZERO
          p_DG_P       = ZERO
          
          iglobP = 1
          if(neighbor(1) > -1) then
            iglobP = ibool_DG(neighbor(1), neighbor(2), neighbor(3))
          endif
          
          exact_interface_flux = .false. ! Reset this variable to .false. in case it does not get assigned during the call to compute_interface_unknowns.
          call compute_interface_unknowns(i, j, ispec, rho_DG_P, rhovx_DG_P, &
                  rhovz_DG_P, E_DG_P, veloc_x_DG_P, veloc_z_DG_P, p_DG_P, T_P, &
                  Tx_DG_P, Tz_DG_P, Vxx_DG_P, Vzz_DG_P, Vzx_DG_P, Vxz_DG_P, gamma_P,&
                  neighbor, MPI_change_cpt, &
                  exact_interface_flux, &
                  rho_DG(iglobM), E_DG(iglobM), rhovx_DG(iglobM), rhovz_DG(iglobM), &
                  V_DG(:,:,iglobM), T_DG(:,iglobM), &
                  rho_DG(iglobP), E_DG(iglobP), rhovx_DG(iglobP), rhovz_DG(iglobP), &
                  V_DG(:,:,iglobP), T_DG(:,iglobP), &
                  MPI_iglob, chosen_nxnz_forMPI, dir_normal, nx, nz, weight, timelocal, elastic_tensor)
          
          flux_rho   = rho_DG_P*veloc_x_DG_P*nx + rho_DG_P*veloc_z_DG_P*nz
          flux_rhovx = (rho_DG_P*veloc_x_DG_P**2 + p_DG_P)*nx + rho_DG_P*veloc_x_DG_P*veloc_z_DG_P*nz
          flux_rhovz = (rho_DG_P*veloc_z_DG_P**2 + p_DG_P)*nz + rho_DG_P*veloc_x_DG_P*veloc_z_DG_P*nx
          flux_E     = veloc_x_DG_P*(E_DG_P + p_DG_P)*nx + veloc_z_DG_P*(E_DG_P + p_DG_P)*nz
          
          ! Recover an approximate local maximum linearized acoustic wave speed. See for example Hesthaven (doi.org/10.1007/9780387720678), page 208.
          lambda = 0.
          jump   = 0.
          !gamma_P = gammaext_DG(iglobM) ! DEBUG
          veloc_n_M = sqrt(veloc_x_DG(iglobM)**2 + veloc_z_DG(iglobM)**2)
          veloc_n_P = sqrt(veloc_x_DG_P**2 + veloc_z_DG_P**2)
          lambda = max( veloc_n_M + sqrt(abs(gammaext_DG(iglobM)*p_DG(iglobM)/rho_DG(iglobM))), &
                        veloc_n_P + sqrt(abs(gamma_P*p_DG_P/rho_DG_P)) )
          
          ! Viscous stress tensor's contributions (already under the form of the average mean flux).
          dux_dx = ZERO
          dux_dz = ZERO
          duz_dx = ZERO
          duz_dz = ZERO
          dT_dx = ZERO
          dT_dz = ZERO
          if(muext(i, j, ispec) > 0 .OR. &
             etaext(i, j, ispec) > 0 .OR. &
             kappa_DG(i, j, ispec)  > 0) then
            dux_dx = 0.5*(V_DG(1, 1, iglobM) + Vxx_DG_P)
            dux_dz = 0.5*(V_DG(1, 2, iglobM) + Vxz_DG_P)
            duz_dx = 0.5*(V_DG(2, 1, iglobM) + Vzx_DG_P)
            duz_dz = 0.5*(V_DG(2, 2, iglobM) + Vzz_DG_P)
            dT_dx = 0.5*(T_DG(1, iglobM) + Tx_DG_P)
            dT_dz = 0.5*(T_DG(2, iglobM) + Tz_DG_P)
          endif
          
          ! Mass conservation equation's contributions (fully inviscid).
          temp_unknown_M  = rhovx_DG(iglobM)
          temp_unknown_P  = rhovx_DG_P
          temp_unknown2_M = rhovz_DG(iglobM)
          temp_unknown2_P = rhovz_DG_P
          ! Dot product.
          flux_x = temp_unknown_M + temp_unknown_P
          flux_z = temp_unknown2_M + temp_unknown2_P
          flux_n = flux_x*nx + flux_z*nz
          jump   = rho_DG(iglobM) - rho_DG_P
          ! Add flux' contribution.
          if(exact_interface_flux) then
            jump = 0.
            !flux_n = flux_rho*2.
          endif
          dot_rho(iglobM) = dot_rho(iglobM) - weight*(flux_n + lambda*jump)*HALF
          
          ! x-Momentum equation's inviscid contributions.
          if(.not. CONSTRAIN_HYDROSTATIC) then
            temp_unknown_M = rho_DG(iglobM)*veloc_x_DG(iglobM)**2 + p_DG(iglobM)
            temp_unknown_P = rho_DG_P*veloc_x_DG_P**2 + p_DG_P
            temp_unknown2_M = rho_DG(iglobM)*veloc_x_DG(iglobM)*veloc_z_DG(iglobM)
            temp_unknown2_P = rho_DG_P*veloc_x_DG_P*veloc_z_DG_P
          else        
            temp_unknown_M = rho_DG(iglobM)*veloc_x_DG(iglobM)**2 + (p_DG(iglobM) - p_DG_init(iglobM))
            temp_unknown_P = rho_DG_P*veloc_x_DG_P**2 + (p_DG_P - p_DG_init(iglobM))
            temp_unknown2_M = rho_DG(iglobM)*veloc_x_DG(iglobM)*veloc_z_DG(iglobM)
            temp_unknown2_P = rho_DG_P*veloc_x_DG_P*veloc_z_DG_P
          endif
          ! Dot product.
          flux_x = temp_unknown_M + temp_unknown_P
          flux_z = temp_unknown2_M + temp_unknown2_P
          flux_n = flux_x*nx + flux_z*nz
          jump   = rhovx_DG(iglobM) - rhovx_DG_P
          ! Add flux' contribution.
          if(exact_interface_flux) then
            jump = 0.
          endif
          dot_rhovx(iglobM) = dot_rhovx(iglobM) - weight*(flux_n + lambda*jump)*HALF

          ! x-Momentum equation's viscous contributions.
          ! TODO: Check the expressions of the viscous momentum flux' terms (temp_unknown and temp_unknown2).
          ! Some of the energy equation's flux' terms are computed and added here. TODO: Explain how.
          ! Recall: dux_dx, duz_dx, dux_dz, and duz_dz already contain the 0.5 factor to put the flux under mean average form.
          temp_unknown = muext(i, j, ispec)*TWO*dux_dx + (etaext(i, j, ispec) - (TWO/3.)*muext(i, j, ispec))*(dux_dx + duz_dz) 
          temp_unknown2 = muext(i, j, ispec)*( dux_dz + duz_dx )
          ! Dot product.
          !flux_x = temp_unknown
          !flux_z = temp_unknown2
          !flux_n = flux_x*nx + flux_z*nz
          flux_n = temp_unknown*nx + temp_unknown2*nz ! [3 operations + 1 affectation], instead of [3 operations + 3 affectations]. Keep the lines above for comprehension.
          dot_rhovx(iglobM) = dot_rhovx(iglobM) + weight*flux_n
          dot_E(iglobM)     = dot_E(iglobM) &
                              + weight * HALF * (  (veloc_x_DG(iglobM) + veloc_x_DG_P) * temp_unknown &
                                                  +(veloc_z_DG(iglobM) + veloc_z_DG_P) * temp_unknown2 )*nx
          
          ! z-Momentum equation's inviscid contributions.
          temp_unknown_M  = rho_DG(iglobM)*veloc_x_DG(iglobM)*veloc_z_DG(iglobM)
          temp_unknown_P  = rho_DG_P*veloc_x_DG_P*veloc_z_DG_P
          if(.not. CONSTRAIN_HYDROSTATIC) then        
            temp_unknown2_M = rho_DG(iglobM)*veloc_z_DG(iglobM)**2 + p_DG(iglobM)
            temp_unknown2_P = rho_DG_P*veloc_z_DG_P**2 + p_DG_P
          else
            temp_unknown2_M = rho_DG(iglobM)*veloc_z_DG(iglobM)**2 + (p_DG(iglobM) - p_DG_init(iglobM))
            temp_unknown2_P = rho_DG_P*veloc_z_DG_P**2 + (p_DG_P - p_DG_init(iglobM))
          endif
          ! Dot product.
          flux_x = temp_unknown_M + temp_unknown_P
          flux_z = temp_unknown2_M + temp_unknown2_P
          flux_n = flux_x*nx + flux_z*nz
          jump   = rhovz_DG(iglobM) - rhovz_DG_P
          ! Add flux' contribution.
          if(exact_interface_flux) then
            jump = 0.
          endif
          dot_rhovz(iglobM) = dot_rhovz(iglobM) - weight*(flux_n + lambda*jump)*HALF
          
          ! z-Momentum equation's viscous contributions.
          ! TODO: Check the expressions of the viscous momentum flux' terms (temp_unknown and temp_unknown2).
          ! Some of the energy equation's flux' terms are computed and added here. TODO: Explain how.
          ! Recall: dux_dx, duz_dx, dux_dz, and duz_dz already contain the 0.5 factor to put the flux under mean average form.
          temp_unknown = muext(i, j, ispec)*( dux_dz + duz_dx )
          temp_unknown2 = muext(i, j, ispec)*TWO*duz_dz + (etaext(i, j, ispec) - (TWO/3.)*muext(i, j, ispec))*(dux_dx + duz_dz)
          ! Dot product.
          !flux_x = temp_unknown
          !flux_z = temp_unknown2
          !flux_n = flux_x*nx + flux_z*nz
          flux_n = temp_unknown*nx + temp_unknown2*nz ! [3 operations + 1 affectation], instead of [3 operations + 3 affectations]. Keep the lines above for comprehension.
          dot_rhovz(iglobM) = dot_rhovz(iglobM) + weight*flux_n
          dot_E(iglobM)     = dot_E(iglobM) &
                              + weight * HALF * (  (veloc_x_DG(iglobM) + veloc_x_DG_P) * temp_unknown &
                                                  +(veloc_z_DG(iglobM) + veloc_z_DG_P) * temp_unknown2 )*nz
          
          ! Energy equation's fully inviscid contributions.
          if(.not. CONSTRAIN_HYDROSTATIC) then
            temp_unknown_M = veloc_x_DG(iglobM)*(E_DG(iglobM) + p_DG(iglobM))
            temp_unknown_P = veloc_x_DG_P*(E_DG_P + p_DG_P)
            temp_unknown2_M = veloc_z_DG(iglobM)*(E_DG(iglobM) + p_DG(iglobM))
            temp_unknown2_P = veloc_z_DG_P*(E_DG_P + p_DG_P)
          else        
            temp_unknown_M = veloc_x_DG(iglobM)*(E_DG(iglobM) + (p_DG(iglobM)- p_DG_init(iglobM)))
            temp_unknown_P = veloc_x_DG_P*(E_DG_P + (p_DG_P - p_DG_init(iglobM)))
            temp_unknown2_M = veloc_z_DG(iglobM)*(E_DG(iglobM) + (p_DG(iglobM)- p_DG_init(iglobM)))
            temp_unknown2_P = veloc_z_DG_P*(E_DG_P + (p_DG_P - p_DG_init(iglobM)))
          endif        
          ! Dot product.
          flux_x = temp_unknown_M + temp_unknown_P
          flux_z = temp_unknown2_M + temp_unknown2_P
          flux_n = flux_x*nx + flux_z*nz
          jump   = E_DG(iglobM) - E_DG_P
          ! Add flux' contribution.
          if(exact_interface_flux) then
            jump = 0.
          endif
          dot_E(iglobM) = dot_E(iglobM) - weight*(flux_n + lambda*jump)*HALF
          
          ! Energy equation's heat flux' contribution (last remaining term, viscous).
          ! Recall: dT_dx, and dT_dx already contain the 0.5 factor to put the flux under mean average form.
          dot_E(iglobM) = dot_E(iglobM) &
                          + weight*( kappa_DG(i, j, ispec)*( dT_dx*nx + dT_dx*nz ) )
          
          ! Increment NGLLZ counter.
          j = j + 1
          
          if(it_corner == 1) then
            ! If we are at a corner and this was the first step, then we need to go one more time.
            j = j - 1
          endif
          if(it_corner == 2) then
            ! If we are at a corner and this was the second step, then we need to reset the corner notification.
            it_corner = 0
          endif

        enddo
      enddo
    endif ! End of test if acoustic element
  enddo ! End of loop on elements.
  
  end subroutine compute_forces_acoustic_DG
  
! ------------------------------------------------------------ !
! compute_viscous_tensors                                      !
! ------------------------------------------------------------ !
   
  subroutine compute_viscous_tensors(T_DG, V_DG, rho_DG, rhovx_DG, rhovz_DG, E_DG, timelocal)

! compute forces in the acoustic elements in forward simulation and in adjoint simulation in adjoint inversion
  
  use constants,only: CUSTOM_REAL,NGLLX,NGLLZ,gamma_euler

  use specfem_par, only: nglob_DG,nspec, ispec_is_acoustic_DG,&!ispec_is_acoustic
                         xix,xiz,gammax,gammaz,jacobian, &
                         wxgll,wzgll, ibool_DG, &
                         hprimewgll_zz, hprimewgll_xx, &
                         hprime_xx, hprime_zz, rmass_inverse_acoustic_DG, &
                         neighbor_DG, neighbor_DG_corner, normal_DG, normal_DG_corner, &
                         weight_DG, weight_DG_corner, dir_normal_DG, dir_normal_DG_corner, is_corner, cnu, &
                         elastic_tensor, &
                         rhovx_init, rhovz_init, rho_init, T_init, CONSTRAIN_HYDROSTATIC
                         
  implicit none
  
  ! local parameters
  integer :: ispec,i,j,k,iglob, iglobM, iglobP, it_corner
  !integer :: ifirstelem,ilastelem
  integer :: i_ex, j_ex, ispec_ex, chosen_nxnz_forMPI, dir_normal
  real(kind=CUSTOM_REAL) :: rho_DG_P, rhovx_DG_P, rhovz_DG_P, &
        E_DG_P, veloc_x_DG_P, veloc_z_DG_P, p_DG_P, T_P, &
        Tx_DG_P, Tz_DG_P, Vxx_DG_P, Vzz_DG_P, Vzx_DG_P, Vxz_DG_P, &
        flux_n, flux_x, flux_z, nx, nz, timelocal, weight, gamma_P
  logical :: exact_interface_flux, MPI_change_cpt
  integer, dimension(nglob_DG) :: MPI_iglob
  integer, dimension(3) :: neighbor
  logical :: neighbor_exists

  ! Jacobian matrix and determinant
  real(kind=CUSTOM_REAL) :: xixl,xizl,gammaxl,gammazl,jacobianl

  real(kind=CUSTOM_REAL) :: temp_Tx, temp_Tz, temp_Vxx, temp_Vzx, temp_Vxz, temp_Vzz
  
  ! Local
  real(kind=CUSTOM_REAL), dimension(2, 2, nglob_DG) :: V_DG
  real(kind=CUSTOM_REAL), dimension(2, nglob_DG) :: T_DG
  real(kind=CUSTOM_REAL), dimension(nglob_DG) :: &
         rho_DG, rhovx_DG, rhovz_DG, E_DG, veloc_x_DG, veloc_z_DG, T, &
         grad_Tx, grad_Tz, grad_Vxx, grad_Vzz, grad_Vxz, grad_Vzx
  
  ! Viscosity
  real(kind=CUSTOM_REAL) :: dux_dxi, dux_dgamma, duz_dxi, duz_dgamma, dT_dxi, dT_dgamma
  real(kind=CUSTOM_REAL) :: dux_dx, dux_dz, duz_dx, duz_dz, dT_dx, dT_dz
  real(kind=CUSTOM_REAL) :: wxl, wzl
  
  ! Parameters
  real(kind=CUSTOM_REAL), parameter :: ZERO = 0._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: ONE  = 1._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: TWO  = 2._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: HALF = 0.5_CUSTOM_REAL
  
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: temp_Tx_1, temp_Tx_2, &
        temp_Tz_1, temp_Tz_2, temp_Vxx_1, temp_Vxx_2, &
        temp_Vxz_1, temp_Vxz_2, temp_Vzx_1, temp_Vzx_2, temp_Vzz_1, temp_Vzz_2
  
  logical :: ADD_SURFACE_TERMS
  real(kind=CUSTOM_REAL) :: vx_init, vz_init
  
  integer :: coef_surface
  
  ADD_SURFACE_TERMS = .true. ! TODO: Decide to set a value (.true. or .false.), or to introduce a parameter in the parfile.
  !T_init = (E_DG/rho_DG - 0.5*((rhovx_DG/rho_DG)**2 + (rhovz_DG/rho_DG)**2))/(cnu)
  
  MPI_iglob = 1
  
  veloc_x_DG = rhovx_DG/rho_DG
  veloc_z_DG = rhovz_DG/rho_DG
  T = (E_DG/rho_DG - 0.5*(veloc_x_DG**2 + veloc_z_DG**2))/(cnu)
  
  grad_Tx = ZERO
  grad_Tz = ZERO
  grad_Vxx = ZERO
  grad_Vzz = ZERO
  grad_Vxz = ZERO
  grad_Vzx = ZERO

  do ispec = 1, nspec ! Loop over elements.
    ! acoustic spectral element
    !if (ispec_is_acoustic(ispec)) then
    if (ispec_is_acoustic_DG(ispec)) then
      
      ! --------------------------- !
      ! Volume terms.               !
      ! --------------------------- !
      
      do j = 1, NGLLZ
        do i = 1, NGLLX
          
          iglob = ibool_DG(i, j, ispec)

          jacobianl = jacobian(i, j, ispec)
        
          xixl = xix(i, j, ispec)
          xizl = xiz(i, j, ispec)
          gammaxl = gammax(i, j, ispec)
          gammazl = gammaz(i, j, ispec)
          
          wzl = wzgll(j)
          wxl = wxgll(i)
          
          if(ADD_SURFACE_TERMS) then
            ! In that case, compute \int_{\Omega^k} \mathcal{T}\Phi d\Omega^k = \int_{\Omega^k} (\nabla T)\Phi d\Omega^k as:
            ! - \int_{\Omega^k} T\nabla\Phi d\Omega^k + \int_{\partial\Omega^k} T\Phi n d\Gamma^k.
            ! Here, prepare the components which will be later used to compute the "- \int_{\Omega^k} T\nabla\Phi d\Omega^k" part.
            
            ! Viscous stress tensors
            if(.not. CONSTRAIN_HYDROSTATIC) then
              temp_Tx_1(i,j)  = wzl * jacobianl * (xixl * T(iglob)) 
              temp_Tz_1(i,j)  = wzl * jacobianl * (xizl * T(iglob)) 
              temp_Vxx_1(i,j) = wzl * jacobianl * (xixl * veloc_x_DG(iglob) )
              temp_Vxz_1(i,j) = wzl * jacobianl * (xizl * veloc_x_DG(iglob) )
              temp_Vzx_1(i,j) = wzl * jacobianl * (xixl * veloc_z_DG(iglob)) 
              temp_Vzz_1(i,j) = wzl * jacobianl * (xizl * veloc_z_DG(iglob)) 
              
              temp_Tx_2(i,j)  = wxl * jacobianl * (gammaxl * T(iglob)) 
              temp_Tz_2(i,j)  = wxl * jacobianl * (gammazl * T(iglob)) 
              temp_Vxx_2(i,j) = wxl * jacobianl * (gammaxl * veloc_x_DG(iglob)) 
              temp_Vxz_2(i,j) = wxl * jacobianl * (gammazl * veloc_x_DG(iglob)) 
              temp_Vzx_2(i,j) = wxl * jacobianl * (gammaxl * veloc_z_DG(iglob)) 
              temp_Vzz_2(i,j) = wxl * jacobianl * (gammazl * veloc_z_DG(iglob)) 
            else
              vx_init = rhovx_init(iglob)/rho_init(iglob)
              vz_init = rhovz_init(iglob)/rho_init(iglob)
              
              temp_Tx_1(i,j)  = wzl * jacobianl * (xixl * (T(iglob) - T_init(iglob)) ) 
              temp_Tz_1(i,j)  = wzl * jacobianl * (xizl * (T(iglob) - T_init(iglob)) ) 
              temp_Vxx_1(i,j) = wzl * jacobianl * (xixl * (veloc_x_DG(iglob) - vx_init) )
              temp_Vxz_1(i,j) = wzl * jacobianl * (xizl * (veloc_x_DG(iglob) - vx_init) )
              temp_Vzx_1(i,j) = wzl * jacobianl * (xixl * (veloc_z_DG(iglob) - vz_init) ) 
              temp_Vzz_1(i,j) = wzl * jacobianl * (xizl * (veloc_z_DG(iglob) - vz_init) ) 
              
              temp_Tx_2(i,j)  = wxl * jacobianl * (gammaxl * (T(iglob) - T_init(iglob)) ) 
              temp_Tz_2(i,j)  = wxl * jacobianl * (gammazl * (T(iglob) - T_init(iglob)) ) 
              temp_Vxx_2(i,j) = wxl * jacobianl * (gammaxl * (veloc_x_DG(iglob) - vx_init) ) 
              temp_Vxz_2(i,j) = wxl * jacobianl * (gammazl * (veloc_x_DG(iglob) - vx_init) ) 
              temp_Vzx_2(i,j) = wxl * jacobianl * (gammaxl * (veloc_z_DG(iglob) - vz_init) ) 
              temp_Vzz_2(i,j) = wxl * jacobianl * (gammazl * (veloc_z_DG(iglob) - vz_init) ) 
            endif
          else
            ! In that case, compute \int_{\Omega^k} \mathcal{T}\Phi d\Omega^k = \int_{\Omega^k} (\nabla T)\Phi d\Omega^k directly as it.
            dux_dxi    = ZERO
            dux_dgamma = ZERO
            duz_dxi    = ZERO
            duz_dgamma = ZERO
            dT_dxi     = ZERO
            dT_dgamma  = ZERO
            
            ! first double loop over GLL points to compute and store gradients
            ! we can merge the two loops because NGLLX == NGLLZ
            do k = 1, NGLLX
              if(.not. CONSTRAIN_HYDROSTATIC) then
                dux_dxi    = dux_dxi    + veloc_x_DG(ibool_DG(k,j,ispec)) * real(hprime_xx(i,k), kind=CUSTOM_REAL)
                dux_dgamma = dux_dgamma + veloc_x_DG(ibool_DG(i,k,ispec)) * real(hprime_zz(j,k), kind=CUSTOM_REAL)
                duz_dxi    = duz_dxi    + veloc_z_DG(ibool_DG(k,j,ispec)) * real(hprime_xx(i,k), kind=CUSTOM_REAL)
                duz_dgamma = duz_dgamma + veloc_z_DG(ibool_DG(i,k,ispec)) * real(hprime_zz(j,k), kind=CUSTOM_REAL)
                dT_dxi     = dT_dxi    + T(ibool_DG(k,j,ispec)) * real(hprime_xx(i,k), kind=CUSTOM_REAL)
                dT_dgamma  = dT_dgamma + T(ibool_DG(i,k,ispec)) * real(hprime_zz(j,k), kind=CUSTOM_REAL)
              else
                vx_init = rhovx_init(ibool_DG(k,j,ispec))/rho_init(ibool_DG(k,j,ispec))
                dux_dxi = dux_dxi    + (veloc_x_DG(ibool_DG(k,j,ispec)) - vx_init) &
                          * real(hprime_xx(i,k), kind=CUSTOM_REAL)
                vx_init = rhovx_init(ibool_DG(i,k,ispec))/rho_init(ibool_DG(i,k,ispec))       
                dux_dgamma = dux_dgamma + (veloc_x_DG(ibool_DG(i,k,ispec)) - vx_init) &
                             * real(hprime_zz(j,k), kind=CUSTOM_REAL)
                vz_init = rhovz_init(ibool_DG(k,j,ispec))/rho_init(ibool_DG(k,j,ispec))
                duz_dxi = duz_dxi    + (veloc_z_DG(ibool_DG(k,j,ispec)) - vz_init) &
                          * real(hprime_xx(i,k), kind=CUSTOM_REAL)
                vz_init = rhovz_init(ibool_DG(i,k,ispec))/rho_init(ibool_DG(i,k,ispec))
                duz_dgamma = duz_dgamma + (veloc_z_DG(ibool_DG(i,k,ispec)) - vz_init) &
                             * real(hprime_zz(j,k), kind=CUSTOM_REAL)
                dT_dxi = dT_dxi    + (T(ibool_DG(k,j,ispec)) - T_init(ibool_DG(k,j,ispec))) &
                         * real(hprime_xx(i,k), kind=CUSTOM_REAL)
                dT_dgamma = dT_dgamma + (T(ibool_DG(i,k,ispec)) - T_init(ibool_DG(i,k,ispec))) &
                            * real(hprime_zz(j,k), kind=CUSTOM_REAL)
              endif
            enddo

            ! Derivatives of velocities.
            dux_dx = ( dux_dxi * xixl + dux_dgamma * gammaxl )
            dux_dz = ( dux_dxi * xizl + dux_dgamma * gammazl )
            duz_dx = ( duz_dxi * xixl + duz_dgamma * gammaxl )
            duz_dz = ( duz_dxi * xizl + duz_dgamma * gammazl )
            ! Derivatives of temperature.
            dT_dx = ( dT_dxi * xixl + dT_dgamma * gammaxl )
            dT_dz = ( dT_dxi * xizl + dT_dgamma * gammazl )
            
            temp_Tx  = dT_dx * jacobianl
            temp_Tz  = dT_dz * jacobianl
            temp_Vxx = dux_dx * jacobianl
            temp_Vxz = dux_dz * jacobianl
            temp_Vzx = duz_dx * jacobianl
            temp_Vzz = duz_dz * jacobianl
            
            grad_Tx(iglob)  = dT_dx
            grad_Tz(iglob)  = dT_dz
            grad_Vxx(iglob) = dux_dx
            grad_Vxz(iglob) = dux_dz
            grad_Vzx(iglob) = duz_dx
            grad_Vzz(iglob) = duz_dz
          endif
        enddo
      enddo
      
      coef_surface = 0
      if(ADD_SURFACE_TERMS) then
        coef_surface = NGLLZ
      endif
      do j = 1, coef_surface
        do i = 1, NGLLX
          iglob = ibool_DG(i, j, ispec)
          ! along x direction and z direction
          ! and assemble the contributions
          do k = 1, NGLLX
            grad_Tx(iglob) = grad_Tx(iglob) - &
                             (temp_Tx_1(k,j) * real(hprimewgll_xx(k,i), kind=CUSTOM_REAL) + &
                             temp_Tx_2(i,k) * real(hprimewgll_zz(k,j), kind=CUSTOM_REAL))
            grad_Tz(iglob) = grad_Tz(iglob) - &
                             (temp_Tz_1(k,j) * real(hprimewgll_xx(k,i), kind=CUSTOM_REAL) + &
                             temp_Tz_2(i,k) * real(hprimewgll_zz(k,j), kind=CUSTOM_REAL))
            grad_Vxx(iglob) = grad_Vxx(iglob) - &
                              (temp_Vxx_1(k,j) * real(hprimewgll_xx(k,i), kind=CUSTOM_REAL) + &
                              temp_Vxx_2(i,k) * real(hprimewgll_zz(k,j), kind=CUSTOM_REAL))
            grad_Vxz(iglob) = grad_Vxz(iglob) - &
                              (temp_Vxz_1(k,j) * real(hprimewgll_xx(k,i), kind=CUSTOM_REAL) + &
                              temp_Vxz_2(i,k) * real(hprimewgll_zz(k,j), kind=CUSTOM_REAL))
            grad_Vzx(iglob) = grad_Vzx(iglob) - &
                              (temp_Vzx_1(k,j) * real(hprimewgll_xx(k,i), kind=CUSTOM_REAL) + &
                              temp_Vzx_2(i,k) * real(hprimewgll_zz(k,j), kind=CUSTOM_REAL))
            grad_Vzz(iglob) = grad_Vzz(iglob) - &
                              (temp_Vzz_1(k,j) * real(hprimewgll_xx(k,i), kind=CUSTOM_REAL) + &
                              temp_Vzz_2(i,k) * real(hprimewgll_zz(k,j), kind=CUSTOM_REAL))
          enddo
        enddo
      enddo
        
      ! --------------------------- !
      ! Interface terms.            !
      ! --------------------------- !
        
      coef_surface = 0
      if(ADD_SURFACE_TERMS) then
        coef_surface = NGLLZ
      endif
      it_corner = 0
      MPI_change_cpt = .false.
      do i = 1, coef_surface
        j = 1
        do while (j <= NGLLZ)
          ! We skip interior points
          if(i > 1 .AND. i < NGLLX .AND. j > 1 .AND. j < NGLLZ) then
            j = j + 1
            cycle
          endif
          
          ! Recover neighbor location
          neighbor = neighbor_DG(i, j, ispec, :)
          
          chosen_nxnz_forMPI = -1
        
          ! Reinit boolean to know if neighbor exists
          neighbor_exists = .false.
          
          !nx = normal_DG(ispec, ind, 1)
          !nz = normal_DG(ispec, ind, 2)
          nx     = normal_DG(i, j, ispec, 1)
          nz     = normal_DG(i, j, ispec, 2)
          weight = weight_DG(i, j, ispec)
          dir_normal = dir_normal_DG(i, j, ispec)
          chosen_nxnz_forMPI = 0
          
          ! Needs x2 points at corners to correctly map edges
          ! => 2 step on the same point
          if(it_corner == 1) then
            neighbor = neighbor_DG_corner(i, j, ispec,:)
            nx     = normal_DG_corner(i, j, ispec, 1)
            nz     = normal_DG_corner(i, j, ispec, 2)
            weight = weight_DG_corner(i, j, ispec)
            dir_normal = dir_normal_DG_corner(i, j, ispec)
            chosen_nxnz_forMPI = 1
            it_corner = 2
          endif
          
          if(neighbor_DG(i, j, ispec, 3) > -1 .OR. &
                  neighbor_DG_corner(i, j, ispec, 3) > -1) neighbor_exists = .true.
          
          ! If not outer boundary check for corresponding neighbor normal
          if(is_corner(i, j)) then
            ! If at least one neighbor exists
            if(neighbor_exists) then
              i_ex = neighbor(1)
              j_ex = neighbor(2)
              ispec_ex = neighbor(3)
              ! If corner of an outside edge
              if(it_corner == 2 .AND. neighbor(3) == -1) then
                i_ex = neighbor_DG(i, j, ispec,1)
                j_ex = neighbor_DG(i, j, ispec,2)
                ispec_ex = neighbor_DG(i, j, ispec,3)
              elseif(it_corner < 2 .AND. neighbor(3) == -1) then
                i_ex = neighbor_DG_corner(i, j, ispec,1)
                j_ex = neighbor_DG_corner(i, j, ispec,2)
                ispec_ex = neighbor_DG_corner(i, j, ispec,3)
              endif
              
              ! Cross product to verify if the normal corresponds to the normal
              !normal_DG(i_ex,j_ex,ispec_ex, 1) 
              !normal_DG(i_ex,j_ex,ispec_ex, 2)
              if( dir_normal /= -dir_normal_DG(i_ex,j_ex,ispec_ex) .AND. &
                dir_normal /= -dir_normal_DG_corner(i_ex,j_ex,ispec_ex) ) then
                ! Only change normal if inner element
                if(neighbor(3) > -1 .AND. it_corner < 2) then
                  nx     = normal_DG_corner(i, j, ispec, 1)
                  nz     = normal_DG_corner(i, j, ispec, 2)
                  weight = weight_DG_corner(i, j, ispec)
                  dir_normal = dir_normal_DG_corner(i, j, ispec)
                  ! MODIF for MPI
                  chosen_nxnz_forMPI = 1
                elseif(neighbor(3) > -1 .AND. it_corner == 2) then
                  nx     = normal_DG(i, j, ispec, 1)
                  nz     = normal_DG(i, j, ispec, 2)
                  weight = weight_DG(i, j, ispec)
                  dir_normal = dir_normal_DG(i, j, ispec)
                  ! MODIF for MPI
                  chosen_nxnz_forMPI = 0
                endif
              ! If outside element, if the normal corresponds to the one computed here
              ! it means that we should take the other one
              elseif(neighbor(3) == -1) then
                ! Only change normal if inner element
                if(it_corner < 2) then
                  nx     = normal_DG_corner(i, j, ispec, 1)
                  nz     = normal_DG_corner(i, j, ispec, 2)
                  weight = weight_DG_corner(i, j, ispec)
                  dir_normal = dir_normal_DG_corner(i, j, ispec)
                  ! MODIF for MPI
                  chosen_nxnz_forMPI = 1
                elseif(it_corner == 2) then
                  nx     = normal_DG(i, j, ispec, 1)
                  nz     = normal_DG(i, j, ispec, 2)
                  weight = weight_DG(i, j, ispec)
                  dir_normal = dir_normal_DG(i, j, ispec)
                  ! MODIF for MPI
                  chosen_nxnz_forMPI = 0
                endif
              endif
            endif
          endif
          
          ! Interior point
          iglobM = ibool_DG(i, j, ispec)
          
          ! If a MPI surface node has been ill referenced and we need to witch between
          ! normal_DG and normal_DG_corner
          if(MPI_change_cpt) then
            if(chosen_nxnz_forMPI == 1) then
              nx     = normal_DG(i, j, ispec, 1)
              nz     = normal_DG(i, j, ispec, 2)
              weight = weight_DG(i, j, ispec)
              dir_normal = dir_normal_DG(i, j, ispec)
            elseif(chosen_nxnz_forMPI == 0) then
              nx     = normal_DG_corner(i, j, ispec, 1)
              nz     = normal_DG_corner(i, j, ispec, 2)
              weight = weight_DG_corner(i, j, ispec)
              dir_normal = dir_normal_DG_corner(i, j, ispec)
            endif
          endif
          
          ! If at corner notify that we will need to go again 
          if(is_corner(i, j) .AND. it_corner == 0) it_corner = 1
          
          rho_DG_P     = ZERO
          rhovx_DG_P   = ZERO
          rhovz_DG_P   = ZERO
          E_DG_P       = ZERO
          veloc_x_DG_P = ZERO
          veloc_z_DG_P = ZERO
          p_DG_P       = ZERO
          T_P          = ZERO
          Vxx_DG_P     = ZERO
          Vzz_DG_P     = ZERO
          Vzx_DG_P     = ZERO
          Vxz_DG_P     = ZERO
          
          exact_interface_flux = .false.
          
          iglobP = 1
          if(neighbor(1) > -1) &
          iglobP = ibool_DG(neighbor(1),neighbor(2),neighbor(3))
          
          ! Get the interface unknowns in order to compute the fluxes.
          call compute_interface_unknowns(i, j, ispec, rho_DG_P, rhovx_DG_P, &
                  rhovz_DG_P, E_DG_P, veloc_x_DG_P, veloc_z_DG_P, p_DG_P, T_P, &
                  Tx_DG_P, Tz_DG_P, Vxx_DG_P, Vzz_DG_P, Vzx_DG_P, Vxz_DG_P, gamma_P, &
                  neighbor, MPI_change_cpt, &
                  exact_interface_flux, &
                  rho_DG(iglobM), E_DG(iglobM), rhovx_DG(iglobM), rhovz_DG(iglobM), &
                  V_DG(:,:,iglobM), T_DG(:,iglobM), &
                  rho_DG(iglobP), E_DG(iglobP), rhovx_DG(iglobP), rhovz_DG(iglobP), &
                  V_DG(:,:,iglobP), T_DG(:,iglobP), &
                  MPI_iglob, chosen_nxnz_forMPI, dir_normal, nx, nz, weight, timelocal,elastic_tensor)

          vx_init = rhovx_init(iglobM)/rho_init(iglobM)
          vz_init = rhovz_init(iglobM)/rho_init(iglobM)

          ! Dot product.
          flux_x = T(iglobM) + T_P
          if(CONSTRAIN_HYDROSTATIC) flux_x = flux_x - 2*T_init(iglobM)
          flux_n = flux_x*nx
          grad_Tx(iglobM) = grad_Tx(iglobM) + weight*flux_n*HALF
          
          ! Dot product.
          flux_z = T(iglobM) + T_P 
          if(CONSTRAIN_HYDROSTATIC) flux_z = flux_z - 2*T_init(iglobM)
          flux_n = flux_z*nz
          grad_Tz(iglobM) = grad_Tz(iglobM) + weight*flux_n*HALF
          
          ! Dot product.
          flux_x = veloc_x_DG(iglobM) + veloc_x_DG_P 
          if(CONSTRAIN_HYDROSTATIC) flux_x = flux_x - 2*vx_init
          flux_n = flux_x*nx
          grad_Vxx(iglobM) = grad_Vxx(iglobM) + weight*flux_n*HALF
          
          ! Dot product.
          flux_z = veloc_x_DG(iglobM) + veloc_x_DG_P 
          if(CONSTRAIN_HYDROSTATIC) flux_z = flux_z - 2*vx_init
          flux_n = flux_z*nz
          grad_Vxz(iglobM) = grad_Vxz(iglobM) + weight*flux_n*HALF
          
          ! Dot product.
          flux_x = veloc_z_DG(iglobM) + veloc_z_DG_P 
          if(CONSTRAIN_HYDROSTATIC) flux_x = flux_x - 2*vz_init
          flux_n = flux_x*nx
          grad_Vzx(iglobM) = grad_Vzx(iglobM) + weight*flux_n*HALF
          
          ! Dot product.
          flux_z = veloc_z_DG(iglobM) + veloc_z_DG_P 
          if(CONSTRAIN_HYDROSTATIC) flux_z = flux_z - 2*vz_init
          flux_n = flux_z*nz
          grad_Vzz(iglobM) = grad_Vzz(iglobM) + weight*flux_n*HALF
          
          ! Increment NGLLZ counter
          j = j + 1
          
          ! If at corner and first step => go again 
          if(it_corner == 1) j = j - 1
          ! Reset corner notification
          if(it_corner == 2) it_corner = 0

        enddo
      enddo
    endif ! End of test if acoustic element.
  enddo ! End of loop over elements.
  
  if(ADD_SURFACE_TERMS) then
    ! TODO: Explain why only if ADD_SURFACE_TERMS == .true..
    grad_Tx(:)  = grad_Tx(:) * rmass_inverse_acoustic_DG(:)
    grad_Tz(:)  = grad_Tz(:) * rmass_inverse_acoustic_DG(:)
    grad_Vxx(:) = grad_Vxx(:) * rmass_inverse_acoustic_DG(:)
    grad_Vxz(:) = grad_Vxz(:) * rmass_inverse_acoustic_DG(:)
    grad_Vzx(:) = grad_Vzx(:) * rmass_inverse_acoustic_DG(:)
    grad_Vzz(:) = grad_Vzz(:) * rmass_inverse_acoustic_DG(:)
  endif
  
  ! Store in variables inteded for output.
  T_DG(1, :) = grad_Tx
  T_DG(2, :) = grad_Tz
  V_DG(1, 1, :) = grad_Vxx
  V_DG(1, 2, :) = grad_Vxz
  V_DG(2, 1, :) = grad_Vzx
  V_DG(2, 2, :) = grad_Vzz
  
  end subroutine compute_viscous_tensors
  

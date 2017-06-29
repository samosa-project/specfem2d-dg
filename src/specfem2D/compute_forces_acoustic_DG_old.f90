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


  subroutine compute_forces_acoustic_DG(rho_DG_main, rhovx_DG_main, rhovz_DG_main, E_DG_main, &
        T_DG_main, V_DG_main, e1_DG, &
        dot_rho, dot_rhovx, dot_rhovz, dot_E, dot_e1, timelocal)


! compute forces in the acoustic elements in forward simulation and in adjoint simulation in adjoint inversion

  use constants,only: CUSTOM_REAL,NGLLX,NGLLZ,gamma_euler

  use specfem_par, only: nglob_DG,nspec,ispec_is_acoustic, &
                         xix,xiz,gammax,gammaz,jacobian, &
                         hprimewgll_xx, &
                         hprimewgll_zz,wxgll,wzgll, &
                         normal_DG, normal_DG_corner, ibool_DG, weight_DG, weight_DG_corner, &
                         neighbor_DG, &
                         neighbor_DG_corner, is_corner, &!ACOUSTIC_FORCING
                         !ispec_is_acoustic_surface, ispec_is_acoustic_surface_corner,  &
                         it,potential_dphi_dx_DG, potential_dphi_dz_DG, ibool, &
                         !hprime_xx, hprime_zz, coord,ibool_before_perio, ispec_is_acoustic_forcing, &
                         elastic_tensor, &!ispec_is_acoustic_coupling_el, veloc_elastic,&
                         dir_normal_DG, dir_normal_DG_corner, &
                         DIR_RIGHT, DIR_LEFT, DIR_UP, DIR_DOWN, &
                         !is_MPI_interface_DG, &! ninterface_acoustic_DG, MPI_transfer, &
                         !      inum_interfaces_acoustic_DG, nibool_interfaces_acoustic_DG, ibool_interfaces_acoustic_dg,&
                         myrank, &
                         !buffer_DG_rho_P, buffer_DG_rhovx_P, buffer_DG_rhovz_P, buffer_DG_E_P, NPROC, &
                         !buffer_DG_Vxx_P, buffer_DG_Vzz_P, buffer_DG_Vxz_P, buffer_DG_Vzx_P, buffer_DG_Tz_P, buffer_DG_Tx_P, &
                         i_stage, p_DG_init, gammaext_DG, muext, etaext, T_init, kappa_DG, cnu, tau_epsilon, tau_sigma, &
                         hprime_xx, hprime_zz!, diag_MPI!, &
                         !max_ibool_interfaces_size_ac,ninterface_acoustic
                         !rho0_DG, E0_DG!, &!this_ibool_is_a_periodic_edge, &
                         !ibool, coord, ibool_before_perio!, ADD_PERIODIC_CONDITIONS, rmass_inverse_acoustic_DG! &
                         !,coord, ibool!, link_DG_CG
                         
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
  
  real(kind=CUSTOM_REAL) :: dux_dxi, dux_dgamma, duz_dxi, duz_dgamma
  
  ! for MPI transfer
  MPI_iglob = 1
  
  !found = .false.
  
  ifirstelem = 1
  ilastelem = nspec
  
  rho_DG   = rho_DG_main
  rhovx_DG = rhovx_DG_main
  rhovz_DG = rhovz_DG_main
  E_DG     = E_DG_main
  
  T_DG = T_DG_main
  V_DG = V_DG_main
  
  ! Init auxiliary unknwons
  veloc_x_DG = rhovx_DG/rho_DG
  veloc_z_DG = rhovz_DG/rho_DG
  p_DG       = (gammaext_DG - ONE)*( E_DG &
        - (HALF)*rho_DG*( veloc_x_DG**2 + veloc_z_DG**2 ) )
        
  if(timelocal == 0) then
        p_DG_init = p_DG
        T_init = (E_DG/rho_DG - 0.5*((rhovx_DG/rho_DG)**2 + (rhovz_DG/rho_DG)**2))/(cnu)
  endif
  !T = (E_DG/rho_DG - 0.5*((rhovx_DG/rho_DG)**2 + (rhovz_DG/rho_DG)**2))/(5/2)
        
  ! Initialization
  dot_rho   = ZERO
  dot_rhovx = ZERO
  dot_rhovz = ZERO
  dot_E     = ZERO
  dot_e1    = ZERO
  
  ! add force source
  !call compute_add_sources_acoustic_DG(dot_rhovx,it,i_stage)
  !call compute_add_sources_acoustic_DG(dot_rhovz,it,i_stage)
  call compute_add_sources_acoustic_DG(dot_E,it,i_stage)
  !call compute_add_sources_acoustic_DG(dot_rho,it,i_stage)   
  
  if(myrank == 0) then
  WRITE(*,*) it,"MAXVAL ", maxval(rho_DG), minval(rho_DG)
  WRITE(*,*) it,"MAXVAL ", maxval(rhovx_DG), minval(rhovx_DG)
  WRITE(*,*) it,"MAXVAL ", maxval(rhovz_DG), minval(rhovz_DG)
  WRITE(*,*) it,"MAXVAL ", maxval(E_DG), minval(E_DG)
  WRITE(*,*) it,"RATIO PRESSURE", maxval(abs((p_DG-p_DG_init)/p_DG_init))
  WRITE(*,*) "*****************"
  endif
  
! loop over spectral elements
  do ispec = ifirstelem,ilastelem

    ! acoustic spectral element
    if (ispec_is_acoustic(ispec)) then

      ! first double loop over GLL points to compute and store gradients
      do j = 1,NGLLZ
        do i = 1,NGLLX
        
          iglob = ibool_DG(i,j,ispec)
          
          jacobianl = jacobian(i,j,ispec)
        
          xixl = xix(i,j,ispec)
          xizl = xiz(i,j,ispec)
          gammaxl = gammax(i,j,ispec)
          gammazl = gammaz(i,j,ispec)
          
          wzl = wzgll(j)
          wxl = wxgll(i)
        
          !!!!!!!!!!!!!!!!!!!!!!!!
          ! Inviscid stress tensor
          temp_unknown = rhovx_DG(iglob)
          temp_unknown2 = rhovz_DG(iglob)
          temp_rho_1(i,j) = wzl * jacobianl * (xixl * temp_unknown + xizl * temp_unknown2) 
          temp_rho_2(i,j) = wxl * jacobianl * (gammaxl * temp_unknown + gammazl * temp_unknown2) 
          
          temp_unknown = rho_DG(iglob)*veloc_x_DG(iglob)**2 + p_DG(iglob)
          temp_unknown2 = rho_DG(iglob)*veloc_x_DG(iglob)*veloc_z_DG(iglob)
          temp_rhovx_1(i,j) = wzl * jacobianl * (xixl * temp_unknown + xizl * temp_unknown2) 
          temp_rhovx_2(i,j) = wxl * jacobianl * (gammaxl * temp_unknown + gammazl * temp_unknown2) 
          
          temp_unknown = rho_DG(iglob)*veloc_x_DG(iglob)*veloc_z_DG(iglob)
          temp_unknown2 = rho_DG(iglob)*veloc_z_DG(iglob)**2 + p_DG(iglob)
          temp_rhovz_1(i,j) = wzl * jacobianl * (xixl * temp_unknown + xizl * temp_unknown2) 
          temp_rhovz_2(i,j) = wxl * jacobianl * (gammaxl * temp_unknown + gammazl * temp_unknown2) 
          
          temp_unknown = veloc_x_DG(iglob)*(E_DG(iglob) + p_DG(iglob))
          temp_unknown2 = veloc_z_DG(iglob)*(E_DG(iglob) + p_DG(iglob))
          temp_E_1(i,j) = wzl * jacobianl * (xixl * temp_unknown + xizl * temp_unknown2) 
          temp_E_2(i,j) = wxl * jacobianl * (gammaxl * temp_unknown + gammazl * temp_unknown2) 
          
          !!!!!!!!!!!!!!!!!!!!!!!
          ! Viscous stress tensor
          if(muext(i,j,ispec) > 0 .OR. etaext(i,j,ispec) > 0 .OR. kappa_DG(i,j,ispec) > 0) then
          
          dux_dx = V_DG(1,1,iglob)
          dux_dz = V_DG(1,2,iglob)
          duz_dx = V_DG(2,1,iglob)
          duz_dz = V_DG(2,2,iglob)
          dT_dx  = T_DG(1,iglob)
          dT_dz  = T_DG(2,iglob)
          
          temp_unknown = muext(i,j,ispec)*TWO*dux_dx + (etaext(i,j,ispec) - (TWO/3.)*muext(i,j,ispec))*(dux_dx + duz_dz) 
          temp_unknown2 = muext(i,j,ispec)*( dux_dz + duz_dx )
          temp_rhovx_1(i,j) = temp_rhovx_1(i,j) - wzl * jacobianl * (xixl * temp_unknown + xizl * temp_unknown2) 
          temp_rhovx_2(i,j) = temp_rhovx_2(i,j) - wxl * jacobianl * (gammaxl * temp_unknown + gammazl * temp_unknown2) 
          
          temp_unknown  = veloc_x_DG(iglob)*temp_unknown + veloc_z_DG(iglob)*temp_unknown2
          temp_E_1(i,j) = temp_E_1(i,j) - wzl * jacobianl * (xixl * temp_unknown) 
          temp_E_2(i,j) = temp_E_2(i,j) - wxl * jacobianl * (gammaxl * temp_unknown) 
          
          temp_unknown = muext(i,j,ispec)*( dux_dz + duz_dx )
          temp_unknown2 = muext(i,j,ispec)*TWO*duz_dz + (etaext(i,j,ispec) - (TWO/3.)*muext(i,j,ispec))*(dux_dx + duz_dz) 
          temp_rhovz_1(i,j) = temp_rhovz_1(i,j) - wzl * jacobianl * (xixl * temp_unknown + xizl * temp_unknown2) 
          temp_rhovz_2(i,j) = temp_rhovz_2(i,j) - wxl * jacobianl * (gammaxl * temp_unknown + gammazl * temp_unknown2) 
          
          temp_unknown2  = veloc_x_DG(iglob)*temp_unknown + veloc_z_DG(iglob)*temp_unknown2
          
          temp_E_1(i,j) = temp_E_1(i,j) - wzl * jacobianl * (xizl * temp_unknown2) 
          temp_E_2(i,j) = temp_E_2(i,j) - wxl * jacobianl * (gammazl * temp_unknown2) 
          
          !!!!!!!!!!!!!
          ! Heat flux
          temp_unknown  = kappa_DG(i,j,ispec)*dT_dx
          temp_unknown2 = kappa_DG(i,j,ispec)*dT_dz
          temp_E_1(i,j) = temp_E_1(i,j) - wzl * jacobianl * (xixl * temp_unknown + xizl * temp_unknown2) 
          temp_E_2(i,j) = temp_E_2(i,j) - wxl * jacobianl * (gammaxl * temp_unknown + gammazl * temp_unknown2) 
           
          endif !if(muext(i,j,ispec) > 0 .OR. etaext(i,j,ispec) > 0) then
                
          !!!!!!!!!!!!!!!!!!!!!
          ! Gravity potentials
          temp_rho_gravi(i,j)   = ZERO
          temp_rhovx_gravi(i,j) = -rho_DG(iglob)*potential_dphi_dx_DG(ibool(i,j,ispec))* jacobianl
          temp_rhovz_gravi(i,j) = -rho_DG(iglob)*potential_dphi_dz_DG(ibool(i,j,ispec))* jacobianl
          temp_E_gravi(i,j)     = -rho_DG(iglob)*(veloc_x_DG(iglob)*potential_dphi_dx_DG(ibool(i,j,ispec)) + &
                veloc_z_DG(iglob)*potential_dphi_dz_DG(ibool(i,j,ispec)))* jacobianl &
                - jacobianl * e1_DG(iglob)/(gammaext_DG(iglob) - ONE)
          
          
          dux_dxi    = ZERO
          dux_dgamma = ZERO
          duz_dxi    = ZERO
          duz_dgamma = ZERO
          
          !if(abs(V_DG(2,2,iglob)) > 0) WRITE(*,*) "V_DG(2,2,iglob): ", V_DG(2,2,iglob) 
          wzl = real(wzgll(j), kind=CUSTOM_REAL)
          wxl = real(wxgll(i), kind=CUSTOM_REAL)
          dot_e1(iglob) = dot_e1(iglob) - (1/tau_sigma(i,j,ispec))*( &
                (1 - (tau_sigma(i,j,ispec)/tau_epsilon(i,j,ispec)))*(dux_dx &
                        + duz_dz)*( gammaext_DG(iglob)*p_DG_init(iglob) ) + e1_DG(iglob) )
          
        enddo
      enddo

!    
! second double-loop over GLL to compute all the terms
!
      do j = 1,NGLLZ
        do i = 1,NGLLX
          iglob = ibool_DG(i,j,ispec)
          ! along x direction and z direction
          ! and assemble the contributions
            do k = 1,NGLLX
            
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
            
            dot_rho(iglob)   = dot_rho(iglob)   + temp_rho_gravi(i,j) * wxl * wzl
            dot_rhovx(iglob) = dot_rhovx(iglob) + temp_rhovx_gravi(i,j) * wxl * wzl
            dot_rhovz(iglob) = dot_rhovz(iglob) + temp_rhovz_gravi(i,j) * wxl * wzl
            dot_E(iglob)     = dot_E(iglob)     + temp_E_gravi(i,j) * wxl * wzl
            
        enddo ! second loop over the GLL points
      enddo
      
      !WRITE(*,*) "1 - dot_rhovx2", maxval(dot_rhovx), minval(dot_rhovx)
      
    it_corner = 0
    MPI_change_cpt = .false.
    
    do  i = 1,NGLLX
      
      j = 1
      do while (j <= NGLLZ)
      
        ! We skip interior points
        if(i > 1 .AND. i < NGLLX .AND. j > 1 .AND. j < NGLLZ) then
                j = j + 1
                cycle
        endif
        
        ! Recover neighbor location
        neighbor = neighbor_DG(i,j,ispec,:)
        
        chosen_nxnz_forMPI = -1
      
        ! Reinit boolean to know if neighbor exists
        neighbor_exists = .false.
        
        !nx = normal_DG(ispec, ind, 1)
        !nz = normal_DG(ispec, ind, 2)
        nx     = normal_DG(i,j,ispec, 1)
        nz     = normal_DG(i,j,ispec, 2)
        weight = weight_DG(i, j, ispec)
        dir_normal = dir_normal_DG(i,j,ispec)
        chosen_nxnz_forMPI = 0
        
        ! Needs x2 points at corners to correctly map edges
        ! => 2 step on the same point
        if(it_corner == 1) then
                neighbor = neighbor_DG_corner(i,j,ispec,:)
                nx     = normal_DG_corner(i,j,ispec, 1)
                nz     = normal_DG_corner(i,j,ispec, 2)
                weight = weight_DG_corner(i, j, ispec)
                dir_normal = dir_normal_DG_corner(i,j,ispec)
                chosen_nxnz_forMPI = 1
                it_corner = 2
        endif
        
        if(neighbor_DG(i,j,ispec,3) > -1 .OR. &
                neighbor_DG_corner(i,j,ispec,3) > -1) neighbor_exists = .true.
        
        ! If not outer boundary check for corresponding neighbor normal
        if(is_corner(i,j)) then
                
                ! If at least one neighbor exists
                if(neighbor_exists) then
                
                i_ex = neighbor(1)
                j_ex = neighbor(2)
                ispec_ex = neighbor(3)
                ! If corner of an outside edge
                if(it_corner == 2 .AND. neighbor(3) == -1) then
                        i_ex = neighbor_DG(i,j,ispec,1)
                        j_ex = neighbor_DG(i,j,ispec,2)
                        ispec_ex = neighbor_DG(i,j,ispec,3)
                elseif(it_corner < 2 .AND. neighbor(3) == -1) then
                        i_ex = neighbor_DG_corner(i,j,ispec,1)
                        j_ex = neighbor_DG_corner(i,j,ispec,2)
                        ispec_ex = neighbor_DG_corner(i,j,ispec,3)
                endif
                
                ! Cross product to verify if the normal corresponds to the normal
                !normal_DG(i_ex,j_ex,ispec_ex, 1) 
                !normal_DG(i_ex,j_ex,ispec_ex, 2)
                if( dir_normal /= -dir_normal_DG(i_ex,j_ex,ispec_ex) .AND. &
                     dir_normal /= -dir_normal_DG_corner(i_ex,j_ex,ispec_ex) ) then
                        ! Only change normal if inner element
                        if(neighbor(3) > -1 .AND. it_corner < 2) then
                                nx     = normal_DG_corner(i,j,ispec, 1)
                                nz     = normal_DG_corner(i,j,ispec, 2)
                                weight = weight_DG_corner(i, j, ispec)
                                dir_normal = dir_normal_DG_corner(i,j,ispec)
                                ! MODIF for MPI
                                chosen_nxnz_forMPI = 1
                        elseif(neighbor(3) > -1 .AND. it_corner == 2) then
                                nx     = normal_DG(i,j,ispec, 1)
                                nz     = normal_DG(i,j,ispec, 2)
                                weight = weight_DG(i, j, ispec)
                                dir_normal = dir_normal_DG(i,j,ispec)
                                ! MODIF for MPI
                                chosen_nxnz_forMPI = 0
                        endif
                ! If outside element, if the normal corresponds to the one computed here
                ! it means that we should take the other one
                elseif(neighbor(3) == -1) then
                        ! Only change normal if inner element
                        if(it_corner < 2) then
                                nx     = normal_DG_corner(i,j,ispec, 1)
                                nz     = normal_DG_corner(i,j,ispec, 2)
                                weight = weight_DG_corner(i, j, ispec)
                                dir_normal = dir_normal_DG_corner(i,j,ispec)
                                ! MODIF for MPI
                                chosen_nxnz_forMPI = 1
                        elseif(it_corner == 2) then
                                nx     = normal_DG(i,j,ispec, 1)
                                nz     = normal_DG(i,j,ispec, 2)
                                weight = weight_DG(i, j, ispec)
                                dir_normal = dir_normal_DG(i,j,ispec)
                                ! MODIF for MPI
                                chosen_nxnz_forMPI = 0
                        endif
                endif
                
                endif
                
        endif
        
        ! Interior point
        iglobM = ibool_DG(i,j,ispec)
        
        ! If a MPI surface node has been ill referenced and we need to witch between
        ! normal_DG and normal_DG_corner
        if(MPI_change_cpt) then
        
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
                
        endif !if(is_MPI_interface_DG(iglobM) .AND. NPROC > 1)
        
        ! If at corner notify that we will need to go again 
        if( is_corner(i,j) .AND. it_corner == 0) it_corner = 1
        
        rho_DG_P     = ZERO
        rhovx_DG_P   = ZERO
        rhovz_DG_P   = ZERO
        E_DG_P       = ZERO
        veloc_x_DG_P = ZERO
        veloc_z_DG_P = ZERO
        p_DG_P       = ZERO
        
        iglobP = 1
        if(neighbor(1) > -1) &
        iglobP = ibool_DG(neighbor(1),neighbor(2),neighbor(3))
        
        exact_interface_flux = .false.
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
        
                flux_rhovx = (rho_DG_P*veloc_x_DG_P**2 + p_DG_P)*nx + rho_DG_P*veloc_x_DG_P*veloc_z_DG_P*nz
                flux_rhovz = (rho_DG_P*veloc_z_DG_P**2 + p_DG_P)*nz + rho_DG_P*veloc_x_DG_P*veloc_z_DG_P*nx
                
                flux_rho = rho_DG_P*veloc_x_DG_P*nx + rho_DG_P*veloc_z_DG_P*nz
                flux_E   = veloc_x_DG_P*(E_DG_P + p_DG_P)*nx + veloc_z_DG_P*(E_DG_P + p_DG_P)*nz
        
        ! Approximate local maximum linearized acoustic wave speed
        ! Lax-Friedrich
        lambda = 0.
        jump   = 0.
        if(.not. exact_interface_flux) then
        
        ! TEST
        !gamma_P = gammaext_DG(iglobM)
        veloc_n_M = sqrt(veloc_x_DG(iglobM)**2+veloc_z_DG(iglobM)**2)
        veloc_n_P = sqrt(veloc_x_DG_P**2+veloc_z_DG_P**2)
        lambda = max( veloc_n_M + sqrt(abs(gammaext_DG(iglobM)*p_DG(iglobM)/rho_DG(iglobM))), &
                veloc_n_P + sqrt(abs(gamma_P*p_DG_P/rho_DG_P)) )
       
        endif

        !!!!!!!!!!!!!!!!!!!!!!!
        ! Viscous stress tensor
        dux_dx = ZERO
        dux_dz = ZERO
        duz_dx = ZERO
        duz_dz = ZERO
        dT_dx = ZERO
        dT_dz = ZERO
        if(muext(i,j,ispec) > 0 .OR. etaext(i,j,ispec) > 0 .OR. kappa_DG(i,j,ispec)  > 0) then

        dux_dx = 0.5*(V_DG(1,1,iglobM) + Vxx_DG_P)
        dux_dz = 0.5*(V_DG(1,2,iglobM) + Vxz_DG_P)
        duz_dx = 0.5*(V_DG(2,1,iglobM) + Vzx_DG_P)
        duz_dz = 0.5*(V_DG(2,2,iglobM) + Vzz_DG_P)
        dT_dx = 0.5*(T_DG(1,iglobM) + Tx_DG_P)
        dT_dz = 0.5*(T_DG(2,iglobM) + Tz_DG_P)
        
        endif ! if(muext(i,j,ispec) > 0)
      
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Mass conservation equation
        if(exact_interface_flux) then
                flux_n = flux_rho
                flux_n = TWO*flux_n
        else
                temp_unknown_M  = rhovx_DG(iglobM)
                temp_unknown_P  = rhovx_DG_P
                
                temp_unknown2_M = rhovz_DG(iglobM)
                temp_unknown2_P = rhovz_DG_P
                
                ! compute dot product
                flux_x = temp_unknown_M + temp_unknown_P
                flux_z = temp_unknown2_M + temp_unknown2_P
                flux_n = flux_x*nx + flux_z*nz
                jump   = rho_DG(iglobM) - rho_DG_P
        endif
        
        dot_rho(iglobM) = dot_rho(iglobM) - weight*(flux_n + lambda*jump)*HALF
        
        !!!!!!!!!!!!!!!!!!!!!
        ! x-Momentum equation
        if(exact_interface_flux) then
                flux_n = flux_rhovx
                flux_n = TWO*flux_n
        else
                temp_unknown_M = rho_DG(iglobM)*veloc_x_DG(iglobM)**2 + p_DG(iglobM)
                temp_unknown_P = rho_DG_P*veloc_x_DG_P**2 + p_DG_P
                
                temp_unknown2_M = rho_DG(iglobM)*veloc_x_DG(iglobM)*veloc_z_DG(iglobM)
                temp_unknown2_P = rho_DG_P*veloc_x_DG_P*veloc_z_DG_P
                
                ! compute dot product
                flux_x = temp_unknown_M + temp_unknown_P
                flux_z = temp_unknown2_M + temp_unknown2_P
                flux_n = flux_x*nx + flux_z*nz
                jump   = rhovx_DG(iglobM) - rhovx_DG_P
        endif
        
        dot_rhovx(iglobM) = dot_rhovx(iglobM) - weight*(flux_n + lambda*jump)*HALF

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Viscosity x-Momentum contribution
        !if(.not. exact_interface_flux) then
                
                temp_unknown = muext(i,j,ispec)*TWO*dux_dx + (etaext(i,j,ispec) - (TWO/3.)*muext(i,j,ispec))*(dux_dx + duz_dz) 
                temp_unknown2 = muext(i,j,ispec)*( dux_dz + duz_dx )
                
                ! compute dot product
                flux_x = temp_unknown
                flux_z = temp_unknown2
                flux_n = flux_x*nx + flux_z*nz
                
        !endif
        
        dot_rhovx(iglobM) = dot_rhovx(iglobM) + weight*flux_n
        dot_E(iglobM)     = dot_E(iglobM) &
                + weight*( 0.5*(veloc_x_DG(iglobM) + veloc_x_DG_P) * temp_unknown &
                + 0.5*(veloc_z_DG(iglobM) + veloc_z_DG_P) * temp_unknown2 )*nx
                ! weight*( veloc_x_DG(iglobM) * temp_unknown + veloc_z_DG(iglobM) * temp_unknown2 )*nx

        !!!!!!!!!!!!!!!!!!!!!
        ! z-Momentum equation
        if(exact_interface_flux) then
                flux_n = flux_rhovz
                flux_n = TWO*flux_n
        else
                temp_unknown_M  = rho_DG(iglobM)*veloc_x_DG(iglobM)*veloc_z_DG(iglobM)
                temp_unknown_P  = rho_DG_P*veloc_x_DG_P*veloc_z_DG_P
                
                temp_unknown2_M = rho_DG(iglobM)*veloc_z_DG(iglobM)**2 + p_DG(iglobM)
                temp_unknown2_P = rho_DG_P*veloc_z_DG_P**2 + p_DG_P
                  
                ! compute dot product
                flux_x = temp_unknown_M + temp_unknown_P
                flux_z = temp_unknown2_M + temp_unknown2_P
                flux_n = flux_x*nx + flux_z*nz
                jump   = rhovz_DG(iglobM) - rhovz_DG_P
                
        endif

        dot_rhovz(iglobM) = dot_rhovz(iglobM) - weight*(flux_n + lambda*jump)*HALF
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Viscosity z-Momentum contribution
        !if(.not. exact_interface_flux) then
                
                temp_unknown = muext(i,j,ispec)*( dux_dz + duz_dx )
                temp_unknown2 = muext(i,j,ispec)*TWO*duz_dz + (etaext(i,j,ispec) - (TWO/3.)*muext(i,j,ispec))*(dux_dx + duz_dz) 
                
                ! compute dot product
                flux_x = temp_unknown
                flux_z = temp_unknown2
                flux_n = flux_x*nx + flux_z*nz
                
        !endif
        
        dot_rhovz(iglobM) = dot_rhovz(iglobM) + weight*flux_n
        dot_E(iglobM)     = dot_E(iglobM) &
                + weight*( 0.5*(veloc_x_DG(iglobM) + veloc_x_DG_P) * temp_unknown &
                + 0.5*(veloc_z_DG(iglobM) + veloc_z_DG_P) * temp_unknown2 )*nz
         
        !!!!!!!!!!!!!!!!!
        ! Energy equation
        if(exact_interface_flux) then
                flux_n = flux_E
                flux_n = TWO*flux_n
        else
                temp_unknown_M = veloc_x_DG(iglobM)*(E_DG(iglobM) + p_DG(iglobM))
                temp_unknown_P = veloc_x_DG_P*(E_DG_P + p_DG_P)
                
                temp_unknown2_M = veloc_z_DG(iglobM)*(E_DG(iglobM) + p_DG(iglobM))
                temp_unknown2_P = veloc_z_DG_P*(E_DG_P + p_DG_P)
                
                ! compute dot product
                flux_x = temp_unknown_M + temp_unknown_P
                flux_z = temp_unknown2_M + temp_unknown2_P
                flux_n = flux_x*nx + flux_z*nz
                jump   = E_DG(iglobM) - E_DG_P
        endif
        dot_E(iglobM) = dot_E(iglobM) - weight*(flux_n + lambda*jump)*HALF
        
        !!!!!!!!!!!!!
        ! Heat flux
        dot_E(iglobM)     = dot_E(iglobM) &
                + weight*( kappa_DG(i,j,ispec)*( dT_dx*nx + dT_dz*nz ) )
        
        ! Increment NGLLZ counter
        j = j + 1
        
        ! If at corner and first step => go again 
        if(it_corner == 1) j = j - 1
        ! Reset corner notification
        if(it_corner == 2) it_corner = 0

      enddo
    enddo
    
    endif ! end of test if acoustic element
    
  enddo
  
  end subroutine compute_forces_acoustic_DG
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  ! Compute initial condition
  subroutine initial_condition_DG()

! compute forces in the acoustic elements in forward simulation and in adjoint simulation in adjoint inversion

  use constants,only: CUSTOM_REAL,NGLLX,NGLLZ,PI

  use specfem_par, only: nglob_DG, ibool_DG, &
        rho_DG, rhovx_DG, rhovz_DG, E_DG, e1_DG, nspec, &
        potential_dphi_dx_DG, potential_dphi_dz_DG, &
        hprime_xx, hprime_zz, jacobian, xix, xiz, gammax, gammaz, ibool, &
        gammaext_DG, rhoext, windxext!, &!gravityext, coord &

  implicit none
  
  real(kind=CUSTOM_REAL), dimension(nglob_DG) :: veloc_x_DG, veloc_z_DG, p_DG
  
  ! Parameters
  real(kind=CUSTOM_REAL), parameter :: ZEROl = 0._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: ONEl  = 1._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: TWOl  = 2._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: SIXl  = 6._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: HALFl = 0.5_CUSTOM_REAL
  
  integer :: ispec, iglob, i, j!, k, iglob_DG
  
  real(kind=CUSTOM_REAL) :: rho_DG_P, rhovx_DG_P, rhovz_DG_P, E_DG_P, &
      veloc_x_DG_P, veloc_z_DG_P, p_DG_P, e1_DG_P
  
  !real(kind=CUSTOM_REAL) :: dphi_dxi, dphi_dgamma, dphi_dzl, &
  !      xixl, xizl, gammaxl, gammazl, &
  !      jacobianl
  
  do ispec = 1,nspec

        ! first double loop over GLL points to compute and store gradients
        do j = 1,NGLLZ
          do i = 1,NGLLX
          
           iglob = ibool_DG(i,j,ispec)
  
           call boundary_condition_DG(i, j, ispec, ZEROl, rho_DG_P, rhovx_DG_P, rhovz_DG_P, E_DG_P, &
                veloc_x_DG_P, veloc_z_DG_P, p_DG_P, e1_DG_P)
  
           rho_DG(iglob)     = rho_DG_P
           veloc_x_DG(iglob) = veloc_x_DG_P
           veloc_z_DG(iglob) = veloc_z_DG_P
           p_DG(iglob)       = p_DG_P
           rhovx_DG(iglob)   = rhovx_DG_P
           rhovz_DG(iglob)   = rhovz_DG_P
           E_DG(iglob)       = E_DG_P
           e1_DG(iglob)      = e1_DG_P
           
          enddo
       enddo
       
       call setup_gravity_potential(potential_dphi_dx_DG, potential_dphi_dz_DG, ispec)
       
  enddo  
  
  end subroutine initial_condition_DG
  
  ! Compute gravity potential derivatives
  subroutine setup_gravity_potential(potential_dphi_dx_DG,potential_dphi_dz_DG,ispec)

! compute forces in the acoustic elements in forward simulation and in adjoint simulation in adjoint inversion

  use constants,only: CUSTOM_REAL, NGLLX, NGLLZ

  use specfem_par, only: xix, xiz, gammax, gammaz, coord, ibool_before_perio, &
        jacobian, hprime_xx, hprime_zz, ibool, nglob, gravityext, coord_interface

  implicit none
  
  ! Parameters
  real(kind=CUSTOM_REAL), parameter :: ZEROl = 0._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: ONEl  = 1._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: TWOl  = 2._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: SIXl  = 6._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: HALFl = 0.5_CUSTOM_REAL
  
  integer :: ispec, i, j, k
  
  real(kind=CUSTOM_REAL) :: dphi_dxi, dphi_dgamma, dphi_dxl, dphi_dzl, &
        xixl, xizl, gammaxl, gammazl, &
        jacobianl
        
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: potential_phi_DG
        
  real(kind=CUSTOM_REAL), dimension(nglob) :: potential_dphi_dx_DG, potential_dphi_dz_DG
  
    do i=1,NGLLX
        do j=1,NGLLZ
        potential_phi_DG(i,j) = &
                gravityext(i,j,ispec)*( real(coord(2,ibool_before_perio(i,j,ispec)),kind=CUSTOM_REAL) - coord_interface )
        enddo
    enddo
    
    do i=1,NGLLX
        do j=1,NGLLZ
  
          ! Gravity potential derivatives
          jacobianl = jacobian(i,j,ispec)
        
          xixl = xix(i,j,ispec)
          xizl = xiz(i,j,ispec)
          gammaxl = gammax(i,j,ispec)
          gammazl = gammaz(i,j,ispec)
           
          dphi_dxi    = ZEROl
          dphi_dgamma = ZEROl

          ! first double loop over GLL points to compute and store gradients
          ! we can merge the two loops because NGLLX == NGLLZ
          do k = 1,NGLLX
            dphi_dxi    = dphi_dxi    + potential_phi_DG(k,j) * hprime_xx(i,k)
            dphi_dgamma = dphi_dgamma + potential_phi_DG(i,k) * hprime_zz(j,k)
          enddo

          ! derivatives of potential
          dphi_dxl = dphi_dxi * xixl + dphi_dgamma * gammaxl
          dphi_dzl = dphi_dxi * xizl + dphi_dgamma * gammazl 
          potential_dphi_dx_DG(ibool(i,j,ispec)) = dphi_dxl
          if(abs( potential_dphi_dx_DG(ibool(i,j,ispec))) < 0.0000000001) potential_dphi_dx_DG(ibool(i,j,ispec)) = 0.
          
          potential_dphi_dz_DG(ibool(i,j,ispec)) = gravityext(i,j,ispec)!dphi_dzl
          
     enddo
  enddo   
          
  end subroutine setup_gravity_potential
  
  ! Compute initial condition
  subroutine boundary_condition_DG(i, j, ispec, timelocal, rho_DG_P, rhovx_DG_P, rhovz_DG_P, E_DG_P, &
      veloc_x_DG_P, veloc_z_DG_P, p_DG_P, e1_DG_P)

! compute forces in the acoustic elements in forward simulation and in adjoint simulation in adjoint inversion

  use constants,only: CUSTOM_REAL,gamma_euler,PI

  use specfem_par, only: ibool_before_perio, ibool_DG, coord, MODEL, &
        rhoext, windxext, pext_DG, gravityext, gammaext_DG, &
        etaext, muext, coord_interface, kappa_DG, cp, cnu, T_init, &
        tau_epsilon, tau_sigma, myrank

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
  !real(kind=CUSTOM_REAL) :: Tl, Tu, rho0!, p0, RR
  !real(kind=CUSTOM_REAL) :: Tl, Tu, p0, rho0, RR, rs!, theta, thetaprime, pibar
  !real(kind=CUSTOM_REAL) :: Tl, Tu, p0, RR, rs, theta0, Nsq, rho0
  !real(kind=CUSTOM_REAL) :: cp, cnu, R, P0, N, theta0, exner, theta
  !real(kind=CUSTOM_REAL) :: x0, z0, r, rs, p0, RR, cnu, cp, theta0, exner
  !real(kind=CUSTOM_REAL) :: x0, z0, r, theta, RR, p0!, exner, theta0
  !real(kind=CUSTOM_REAL) :: &
  !      lz, vs, vs_x, vs_z, Ms, p_1, rho_1, veloc_x_1, veloc_z_1, p_2, rho_2, veloc_x_2, veloc_z_2 
  
  ! Coordinate of the solid-fluid interface
  coord_interface = 5000.
  
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
   
   H = 15000!30964!30964.2932
   G = 9.81
   rho_DG_P = ONEl*exp(-z/H)
   p_DG_P   = rho_DG_P*G*H
   
   if(timelocal == 0) then
   gravityext(i,j,ispec) = G!9.831
   gammaext_DG(ibool_DG(i,j,ispec))   = gamma_euler!(rho_DG_P/p_DG_P)*(vpext(i,j,ispec)**2)!gamma_euler!
   muext(i,j,ispec)  = 2e-05
   etaext(i,j,ispec)  = 0!(4/3)*muext(i,j,ispec)
   kappa_DG(i,j,ispec) = 0.2
   cp = 7/2
   cnu = 5/2
   endif
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! If not an external model => means that no initial conditions specified
  else
  
  ! rely on Mars model (large) published in SSR et 1km altitude
  ! AGW_FD_2016_eos/DATA_FD/atmos_model/MARS_MCD_AGW_model_INSIGHT_LS_240_LT_12h.txt
  ! and MARS_MCD_AGW_model_INSIGHT_LS_90_LT_12h_large.dat
  ! At initial time if no external model we create an uniform gravity fields
  if(timelocal == 0) then
        gravityext(i,j,ispec) = 0.0 !3.7247 
        gammaext_DG(ibool_DG(i,j,ispec))   = 1.33076167! 1.4!gamma_euler
        mu_visco  = 1.092656e-05
!        eta_visco = mu_visco*(4/3)!4!(4/3)*mu_visco!1!4e-05!1!4e-05!
        eta_visco = 2.45033366e-05
        muext(i,j,ispec)  = mu_visco
        etaext(i,j,ispec) = eta_visco
!        kappa_DG(i,j,ispec) = 0.0
        kappa_DG(i,j,ispec) = 4.79046750E-04
       tau_epsilon(i,j,ispec) = 9.22879479e-04 !1.5!2.!1.5
        tau_sigma(i,j,ispec)   = 8.99908325e-04 !1.!1./(8*PI**2)!0.025!1
  endif
  
  e1_DG_P = 0.
  
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! BASIC FORCING
   H = 10000
   G = gravityext(i,j,ispec)
   ! http://catalog.conveyorspneumatic.com/Asset/FLS%20Specific%20Heat%20Capacities%20of%20Gases.pdf
   cp = 33.4300880 !7!1010!7/2
   cnu = 25.1155891 !5!718!5/2
   rho_DG_P     = 1.579412d-2 !exp(-z/H)
   
!   veloc_x_DG_P = ZEROl!20
   veloc_x_DG_P = -1.06393795e+01
   veloc_z_DG_P = ZEROl
   
   ! Acoustic only
   p_DG_P = (234.0151**2)*rho_DG_P/gammaext_DG(ibool_DG(i,j,ispec))
   ! Gravi
   !p_DG_P = rho_DG_P*G*H
   
   endif
   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Acoustic plane wave forcing
   to = 2./5.
   perio = 1./5.
   if(z == ZEROl .AND. .false.) &     
   veloc_z_DG_P = 0.001*(&
                  - (2d0/(perio/4d0))*((timelocal-(to-perio/4d0))/(perio/4d0))* &
                           (exp(-((timelocal-(to-perio/4d0))/(perio/4d0))**2)) &
                  + (2d0/(perio/4d0))*((timelocal-(to+perio/4d0))/(perio/4d0))* &
                           (exp(-((timelocal-(to+perio/4d0))/(perio/4d0))**2)) ) 
   
   !!!!!!!!!!!!!!!!!!!!!!!!
   ! Gravi wave forcing
   lambdo = 20000
   perio  = 400!200
   xo = 75000
   to = 350!250.
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
  
  ! SLOPE LIMITER
  
  subroutine SlopeLimit1(Q, timelocal, type_var)
  
    use constants,only: CUSTOM_REAL,NGLLX,NGLLZ

    use specfem_par, only: nspec, nglob_DG, ibool_DG, Vandermonde, invVandermonde, &
        neighbor_DG, neighbor_DG_corner,is_corner, NGLLX, NGLLZ, &
        max_interface_size, ninterface, NPROC, MPI_transfer!, myrank!, ibool, ibool_before_perio!, coord, i_stage&
        !ispec_is_acoustic_surface, ispec_is_acoustic_surface_corner!,ibool_before_perio,

    implicit none 
  
    integer type_var
    real(kind=CUSTOM_REAL) timelocal
    real(kind=CUSTOM_REAL) :: rho_DG_P, rhovx_DG_P, rhovz_DG_P, E_DG_P, &
                veloc_x_DG_P, veloc_z_DG_P, p_DG_P, e1_DG_P
    
    integer :: iglob, iglob2, ispec, j, k, l, h, prod_ind
    integer, parameter :: NGLL = NGLLX*NGLLZ
    real(kind=CUSTOM_REAL) :: vkm1_x,vkp1_x,vkm1_z,vkp1_z, vkm1_x_save
    real(kind=CUSTOM_REAL), dimension(nglob_DG) :: Q, ulimit
    real(kind=CUSTOM_REAL), dimension(NGLL) :: uh, uavg, uhl, ulimit_temp
    real(kind=CUSTOM_REAL), dimension(1:nspec,1:NGLL) :: ul
    real(kind=CUSTOM_REAL), dimension(1:nspec) :: v!, vk
    real(kind=CUSTOM_REAL), dimension(nglob_DG) :: v_MPI
    real(kind=CUSTOM_REAL), dimension(NGLLX*max_interface_size,ninterface) :: buffer_v
    
    integer, dimension(nglob_DG) :: MPI_iglob
    integer :: ipoin,num_interface
    logical, dimension(NGLLX, NGLLZ) :: ispec_bound
    logical ispec_bound_MPI
    !real(kind=CUSTOM_REAL), dimension(1:NGLLX,1:NGLLZ) :: val_in_elmnt
    
    logical :: activate
    
    ! Parameters
    real(kind=CUSTOM_REAL), parameter :: ZEROl = 0._CUSTOM_REAL
    
    MPI_iglob = 1

    ulimit = ZEROl
    buffer_v = ZEROl

    do ispec=1,nspec
    
    uh = ZEROl
    ! Compute cell averages
    ! Compute => uh = invV*u
    do k=1,NGLL
        prod_ind = 0
        do j=1,NGLLX
          do h=1,NGLLZ
                iglob2 = ibool_DG(j,h,ispec)
                prod_ind = prod_ind + 1
                uh(k)  = uh(k) + &
                        invVandermonde(k,prod_ind)*Q(iglob2)
               ! val_in_elmnt(j,h) = Q(iglob2)
          enddo
      enddo
    enddo
    
    ! Extract linear polynomial
    uhl = uh
    uhl(3:NGLLX) = ZEROl
    uhl(NGLLX+2:NGLL) = ZEROl
    
    ul(ispec,:) = ZEROl
    ! => ul = V*uhl;
    ! Compute => uavg = V*uh
    do k=1,NGLL
        do j=1,NGLL
                ul(ispec,k)  = ul(ispec,k) + Vandermonde(k,j)*uhl(j)
        enddo
    enddo
    
    ! Remove high order coef (> 0)
    uh(2:NGLL) = ZEROl
    
    ! Init
    uavg = ZEROl
    ! Compute => uavg = V*uh
    do k=1,NGLL
        do j=1,NGLL
                uavg(k)  = uavg(k) + Vandermonde(k,j)*uh(j)
        enddo
    enddo
    ! Store cell average
   v(ispec) = uavg(1)
   
   do k=1,NGLLX
          do l=1,NGLLZ
          v_MPI(ibool_DG(k,l,ispec)) = v(ispec)
          enddo
   enddo
   
   enddo ! do ispec=1,NSPEC
   
   ! find cell averages
   !vk = v
   
   call assemble_MPI_vector_DG(v_MPI, buffer_v)
   
   !nb_mod = 0
    do ispec=1,nspec
    
     ! Test with boundary conditions
     call boundary_condition_DG(1, NGLLZ/2, ispec, timelocal, rho_DG_P, rhovx_DG_P, rhovz_DG_P, E_DG_P, &
                veloc_x_DG_P, veloc_z_DG_P, p_DG_P, e1_DG_P)
    
     if(type_var == 1) vkm1_x = rho_DG_P
     if(type_var == 2) vkm1_x = rhovx_DG_P
     if(type_var == 3) vkm1_x = rhovz_DG_P
     if(type_var == 4) vkm1_x = E_DG_P
     
     call boundary_condition_DG(NGLLX, NGLLZ/2, ispec, timelocal, rho_DG_P, rhovx_DG_P, rhovz_DG_P, E_DG_P, &
                veloc_x_DG_P, veloc_z_DG_P, p_DG_P, e1_DG_P)
    
     if(type_var == 1) vkp1_x = rho_DG_P
     if(type_var == 2) vkp1_x = rhovx_DG_P
     if(type_var == 3) vkp1_x = rhovz_DG_P
     if(type_var == 4) vkp1_x = E_DG_P
     
     call boundary_condition_DG(NGLLX/2, 1, ispec, timelocal, rho_DG_P, rhovx_DG_P, rhovz_DG_P, E_DG_P, &
                veloc_x_DG_P, veloc_z_DG_P, p_DG_P, e1_DG_P)
    
     if(type_var == 1) vkm1_z = rho_DG_P
     if(type_var == 2) vkm1_z = rhovx_DG_P
     if(type_var == 3) vkm1_z = rhovz_DG_P
     if(type_var == 4) vkm1_z = E_DG_P
     
     call boundary_condition_DG(NGLLX/2, NGLLZ, ispec, timelocal, rho_DG_P, rhovx_DG_P, rhovz_DG_P, E_DG_P, &
                veloc_x_DG_P, veloc_z_DG_P, p_DG_P, e1_DG_P)
    
     if(type_var == 1) vkp1_z = rho_DG_P
     if(type_var == 2) vkp1_z = rhovx_DG_P
     if(type_var == 3) vkp1_z = rhovz_DG_P
     if(type_var == 4) vkp1_z = E_DG_P
    
     if(neighbor_DG(1,NGLLZ/2,ispec,3) > - 1) vkm1_x = v(neighbor_DG(1,NGLLZ/2,ispec,3))
     
     if(neighbor_DG(NGLLX,NGLLZ/2,ispec,3) > - 1) vkp1_x = v(neighbor_DG(NGLLX,NGLLZ/2,ispec,3))
     
     if(neighbor_DG(NGLLX/2,1,ispec,3) > - 1) vkm1_z = v(neighbor_DG(NGLLX/2,1,ispec,3))
     
     if(neighbor_DG(NGLLX/2,NGLLZ,ispec,3) > - 1) vkp1_z = v(neighbor_DG(NGLLX/2,NGLLZ,ispec,3)) !v(neighbor_DG(NGLLX/2,NGLLZ,ispec,3))
     
     ispec_bound_MPI = .false.
     do k=1,NGLLX
     do l=1,NGLLZ
     
     if( (k == NGLLX/2 .AND. l == NGLLZ) .OR. &
        (k == 1 .AND. l == NGLLZ/2) .OR. &
        (k == NGLLX .AND. l == NGLLZ/2) .OR. &
        (k == NGLLX/2 .AND. l == 1) ) then
     
     iglob = ibool_DG(k,l,ispec)
     !!!!!!!!!!!!!!!!!!!!!!!!
     ! MPI neighbor
     ipoin         = -1
     num_interface = -1
     if(NPROC > 1) then
        ipoin         = MPI_transfer(iglob,1,1)
        num_interface = MPI_transfer(iglob,1,2)
     endif                  
     ! If MPI neighbor, but not diagonal corner element and not outside element
     if(ipoin > -1) then
               
               vkm1_x_save = vkm1_x
               if(k == 1 .AND. l == NGLLZ/2) vkm1_x = buffer_v(ipoin,num_interface)
               if(k == NGLLX .AND. l == NGLLZ/2) vkp1_x = buffer_v(ipoin,num_interface)
               if(k == NGLLX/2 .AND. l == 1) vkm1_z = buffer_v(ipoin,num_interface)
               if(k == NGLLX/2 .AND. l == NGLLZ) vkp1_z = buffer_v(ipoin,num_interface)
               
               ispec_bound_MPI = .true.
               
     endif
     
     endif ! if( (k == NGLLX/2 .AND. l = NGLLZ)
     
     enddo
     enddo
     
     ! Compute Limited flux
     call SlopeLimitLin(ulimit_temp, ispec, ul(ispec,:),v(ispec),&
        vkm1_x,vkp1_x,vkm1_z,vkp1_z,activate,type_var)
     
     ! Trick to avoid BC issues for slope limiter
     ispec_bound = .false.
     do k=1,0!NGLLX
       !if(ispec_bound) exit
       do l=1,NGLLZ
       iglob = ibool_DG(k,l,ispec)
       !!!!!!!!!!!!!!!!!!!!!!!!
       ! MPI neighbor
       ipoin         = -1
       num_interface = -1
       if(NPROC > 1) then
         ipoin         = MPI_transfer(iglob,MPI_iglob(iglob),1)
         num_interface = MPI_transfer(iglob,MPI_iglob(iglob),2)
         MPI_iglob(iglob) = MPI_iglob(iglob) + 1
       endif                  
       
        if((neighbor_DG(k,l,ispec,3) == -1 .OR. (is_corner(k,l) .AND. &
                neighbor_DG_corner(k,l,ispec,3) == -1)) .AND. &
                ipoin == -1) then
        endif
        
       enddo
     enddo
     
     prod_ind = 0
     do k=1,NGLLX
       do l=1,NGLLZ
        iglob = ibool_DG(k,l,ispec)
        prod_ind = prod_ind + 1
        ulimit(iglob) = ulimit_temp(prod_ind)
        
        if(.not. activate) ulimit(iglob) = Q(iglob)
       enddo
     enddo
   enddo
   
   ! Update solution
   Q(1:nglob_DG) = ulimit(1:nglob_DG)
   
  end subroutine SlopeLimit1
  
  subroutine SlopeLimitLin(ulimit, ispec, ul, v0, vkm1_x,vkp1_x,vkm1_z,vkp1_z,activate,type_unknown)

    use constants,only: CUSTOM_REAL,NGLLX,NGLLZ

    use specfem_par, only: ibool_before_perio, coord, &
        xix, xiz, gammax, gammaz, hprime_xx, hprime_zz!, Drx, Drz

    implicit none 

        !double precision, dimension(0:NGLOB_DG+1) :: ulimit
        integer, parameter :: NGLL = NGLLX*NGLLZ
        
        real(kind=CUSTOM_REAL), dimension(1:NGLL) :: x0, ul, xl, ux, uz, ulimit
        real(kind=CUSTOM_REAL), dimension(1:NGLL) :: z0, zl
        real(kind=CUSTOM_REAL), dimension(1:3,1)  :: minmod
        real(kind=CUSTOM_REAL) :: hx, hz, v0, vkm1_x,vkp1_x,vkm1_z,vkp1_z, &
                xixl, xizl, gammaxl, gammazl, du_dxi, du_dgamma
        real(kind=CUSTOM_REAL), dimension(1)  :: ulimit_temp
        integer :: ispec, iglob, j, k, m, n, i, type_unknown
        real(kind=CUSTOM_REAL), dimension(1:NGLLX,1:NGLLZ) :: ul_loc
        
        ! Parameters
        real(kind=CUSTOM_REAL), parameter :: ZEROl = 0._CUSTOM_REAL
        real(kind=CUSTOM_REAL), parameter :: ONEl  = 1._CUSTOM_REAL
        real(kind=CUSTOM_REAL), parameter :: TWOl  = 2._CUSTOM_REAL
        real(kind=CUSTOM_REAL), parameter :: HALFl = ONEl/TWOl
        
        real(kind=CUSTOM_REAL), parameter :: gradient_factor = ONEl!/TWOl

        logical :: activate_temp, activate

        ! function ulimit = SlopeLimitLin(ul,xl,vm1,v0,vp1);
        ! Purpose: Apply slopelimited on linear function ul(Np,1) on x(Np,1)
        !          (vm1,v0,vp1) are cell averages left, center, and right
        k = 0
        do m=1,NGLLX
          do n=1,NGLLZ
                k = k + 1
                iglob = ibool_before_perio(m,n,ispec)
                xl(k) = real(coord(1,iglob), kind=CUSTOM_REAL)
                zl(k) = real(coord(2,iglob), kind=CUSTOM_REAL)
          enddo
        enddo
        
        ! Assume that x, z > 0
        hx  = maxval(xl) - minval(xl)
        !WRITE(*,*) "hx",hx
        hz  = maxval(zl) - minval(zl)
        
        ! Coordinates of the element centeroid
        x0 = minval(xl) + hx*HALFl
        z0 = minval(zl) + hz*HALFl
        
        ! Local version of the limited flux
        m = 0
        do i=1,NGLLX
          do j=1,NGLLZ
                m = m + 1
                ul_loc(i,j) = ul(m)
          enddo
        enddo
              
        ! Compute 1st order derivatives at center        
        ux = ZEROl
        uz = ZEROl
        m = 0
        do i=1,NGLLX
          do j=1,NGLLZ
                m = m + 1
                
                du_dxi    = ZEROl
                du_dgamma = ZEROl
                ! first double loop over GLL points to compute and store gradients
                ! we can merge the two loops because NGLLX == NGLLZ
                do k = 1,NGLLX
                  du_dxi    = du_dxi + ul_loc(k,j) * hprime_xx(i,k)
                  du_dgamma = du_dgamma + ul_loc(i,k) * hprime_zz(j,k)
                enddo

                xixl = xix(i,j,ispec)
                xizl = xiz(i,j,ispec)
                gammaxl = gammax(i,j,ispec)
                gammazl = gammaz(i,j,ispec)

                ! derivatives of potential
                ux(m) = du_dxi * xixl + du_dgamma * gammaxl
                uz(m) = du_dxi * xizl + du_dgamma * gammazl
                
          enddo
        enddo
        
        ulimit = v0
        minmod(1,1) = 0
        do k=1,NGLLX*NGLLZ
                minmod(1,1) = minmod(1,1) + (1/NGLLX*NGLLZ)*ux(k)
        enddo
        minmod(1,1) = ux(NGLLX*NGLLZ/2)!+ux(NGLLZ))!0.5*(ux(1)+ux(NGLLX*(NGLLZ-1)+1))!
        minmod(2,1) = (vkp1_x-v0)/(hx*gradient_factor)
        minmod(3,1) = (v0-vkm1_x)/(hx*gradient_factor)
        
        activate = .false.
        
        call minmod_computeB(ulimit_temp, minmod, 3, 1,hx, activate_temp)
        
        if(activate_temp) activate = .true.
        
        ulimit = ulimit + (xl-x0)*ulimit_temp(1) 
         
        minmod(1,1) = 0
        do k=1,NGLLX*NGLLZ
                minmod(1,1) = minmod(1,1) + (1/NGLLX*NGLLZ)*uz(k)
        enddo
        minmod(1,1) = uz(NGLLX*NGLLZ/2)!*(uz(1)+uz(NGLLX*(NGLLZ-1)+1)) :0.5*(uz(1)+uz(NGLLX*(NGLLZ-1)+1))!
        minmod(2,1) = (vkp1_z-v0)/(hz*gradient_factor)
        minmod(3,1) = (v0-vkm1_z)/(hz*gradient_factor)
        
        call minmod_computeB(ulimit_temp, minmod, 3, 1,hz, activate_temp)
        
        if(activate_temp) activate = .true.
        
        ulimit = ulimit + (zl-z0)*ulimit_temp(1)
        
        ! Positivity preserving ?
        do k=1,NGLLX*NGLLZ
        if((type_unknown == 1 .OR. type_unknown == 4) &
                 .AND. ulimit(k) < 1d-10) ulimit(k) = 1d-10
        enddo
        
  end subroutine SlopeLimitLin
  
  subroutine minmod_computeB(R, v, n, m, h, activate)

     use constants,only: CUSTOM_REAL

     implicit none 

        ! Usually 
        ! m = NSPEC
        ! n = 3
        integer :: n, m
        real(kind=CUSTOM_REAL), dimension(1:n,1:m) :: v
        real(kind=CUSTOM_REAL), dimension(1:m) :: R
        
        real(kind=CUSTOM_REAL) :: h, M_param
        
        real(kind=CUSTOM_REAL), parameter :: ZEROl = 0._CUSTOM_REAL
        real(kind=CUSTOM_REAL), parameter :: ONEl  = 1._CUSTOM_REAL
        
        logical activate
        
        ! function mfunc = minmod(v)
        ! Purpose: Implement the midmod function v is a vector

        M_param = 0.03!1.
        
        ! Find regions where limiter is needed
        if(abs(v(1,1)) > M_param*h**2) then
                call minmod_compute(R, v, n, m)
                activate = .true.
                !WRITE(*,*) "ispec modified", ispec, R, v(1,1)
        else
                R = v(1,1)
                activate = .false.
        endif
        
   end subroutine minmod_computeB
  
  subroutine minmod_compute(R, v, n, m)

     use constants,only: CUSTOM_REAL

     implicit none 

        ! Usually 
        ! m = NSPEC
        ! n = 3
        integer :: n, m, k, j
        real(kind=CUSTOM_REAL), dimension(1:n,1:m) :: v
        real(kind=CUSTOM_REAL), dimension(1:m) :: R, s
        
        real(kind=CUSTOM_REAL), parameter :: ZEROl = 0._CUSTOM_REAL
        real(kind=CUSTOM_REAL), parameter :: ONEl  = 1._CUSTOM_REAL
        
        ! function mfunc = minmod(v)
        ! Purpose: Implement the midmod function v is a vector

        R = ZEROl
        !R = v(1,:)
        
        ! Find regions where limiter is needed
        !if(abs(v(1,1)) > M_param*h_param**2) then
        
        ! Sum of lines
        s = ZEROl
        ! Skip through columns
        do k=1,m
             ! Skip through lines
             do j=1,n
             s(k) = s(k) + sign(ONEl,v(j,k))/REAL(n, kind=CUSTOM_REAL)
             enddo
        enddo
        
        do k=1,m
          if(abs(s(k)) == ONEl) then
                R(k) = s(k)*minval(abs(v(:,k)))
          endif
        enddo
        
        !endif
        
   end subroutine minmod_compute
   
   subroutine compute_viscous_tensors(T_DG, V_DG, rho_DG, rhovx_DG, rhovz_DG, E_DG, timelocal)

! compute forces in the acoustic elements in forward simulation and in adjoint simulation in adjoint inversion
  
  use constants,only: CUSTOM_REAL,NGLLX,NGLLZ,gamma_euler

  use specfem_par, only: nglob_DG,nspec,ispec_is_acoustic, &
                         xix,xiz,gammax,gammaz,jacobian, &
                         wxgll,wzgll, ibool_DG, &
                         hprimewgll_zz, hprimewgll_xx, &
                         hprime_xx, hprime_zz, rmass_inverse_acoustic_DG, &
                         neighbor_DG, neighbor_DG_corner, normal_DG, normal_DG_corner, &
                         weight_DG, weight_DG_corner, dir_normal_DG, dir_normal_DG_corner, is_corner, cnu, &
                         elastic_tensor
                         
  implicit none
  
  ! local parameters
  integer :: ispec,i,j,k,iglob, iglobM, iglobP, it_corner
  integer :: ifirstelem,ilastelem
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
  
  real(kind=CUSTOM_REAL), parameter :: threshold = 0.000000001_CUSTOM_REAL
  
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: temp_Tx_1, temp_Tx_2, &
        temp_Tz_1, temp_Tz_2, temp_Vxx_1, temp_Vxx_2, &
        temp_Vxz_1, temp_Vxz_2, temp_Vzx_1, temp_Vzx_2, temp_Vzz_1, temp_Vzz_2
  
  ifirstelem = 1
  ilastelem = nspec
  
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
  
! loop over spectral elements
  do ispec = ifirstelem,ilastelem

    ! acoustic spectral element
    if (ispec_is_acoustic(ispec)) then

      ! first double loop over GLL points to compute and store gradients
      do j = 1,NGLLZ
        do i = 1,NGLLX
          
          iglob = ibool_DG(i,j,ispec)

          jacobianl = jacobian(i,j,ispec)
        
          xixl = xix(i,j,ispec)
          xizl = xiz(i,j,ispec)
          gammaxl = gammax(i,j,ispec)
          gammazl = gammaz(i,j,ispec)
          
          wzl = wzgll(j)
          wxl = wxgll(i)
          
          if(.true.) then
          
          !!!!!!!!!!!!!!!!!!!!!!!
          ! Viscous stress tensors
          
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
          
          dux_dxi    = ZERO
          dux_dgamma = ZERO
          duz_dxi    = ZERO
          duz_dgamma = ZERO
          dT_dxi    = ZERO
          dT_dgamma = ZERO
          
          ! first double loop over GLL points to compute and store gradients
          ! we can merge the two loops because NGLLX == NGLLZ
          do k = 1,NGLLX
                dux_dxi    = dux_dxi    + veloc_x_DG(ibool_DG(k,j,ispec)) * real(hprime_xx(i,k), kind=CUSTOM_REAL)
                dux_dgamma = dux_dgamma + veloc_x_DG(ibool_DG(i,k,ispec)) * real(hprime_zz(j,k), kind=CUSTOM_REAL)
                duz_dxi    = duz_dxi    + veloc_z_DG(ibool_DG(k,j,ispec)) * real(hprime_xx(i,k), kind=CUSTOM_REAL)
                duz_dgamma = duz_dgamma + veloc_z_DG(ibool_DG(i,k,ispec)) * real(hprime_zz(j,k), kind=CUSTOM_REAL)
                dT_dxi    = dT_dxi    + T(ibool_DG(k,j,ispec)) * real(hprime_xx(i,k), kind=CUSTOM_REAL)
                dT_dgamma = dT_dgamma + T(ibool_DG(i,k,ispec)) * real(hprime_zz(j,k), kind=CUSTOM_REAL)
          enddo

          ! derivatives of velocities
          dux_dx = ( dux_dxi * xixl + dux_dgamma * gammaxl )
          dux_dz = ( dux_dxi * xizl + dux_dgamma * gammazl )
          duz_dx = ( duz_dxi * xixl + duz_dgamma * gammaxl )
          duz_dz = ( duz_dxi * xizl + duz_dgamma * gammazl )
          ! derivative of temperature
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

        enddo ! second loop over the GLL points
      enddo
      
      do j = 1,NGLLZ
        do i = 1,NGLLX
        
          iglob = ibool_DG(i,j,ispec)
          ! along x direction and z direction
          ! and assemble the contributions
          do k = 1,NGLLX
            
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
            
        enddo ! second loop over the GLL points
      enddo
      
    !!!!!!!!!! INTERFACE !!!!!!!!!!!!!!!
      
    it_corner = 0
    MPI_change_cpt = .false.
    
    do  i = 1,NGLLZ
      
      j = 1
      do while (j <= NGLLZ)
      
        ! We skip interior points
        if(i > 1 .AND. i < NGLLX .AND. j > 1 .AND. j < NGLLZ) then
                j = j + 1
                cycle
        endif
        
        ! Recover neighbor location
        neighbor = neighbor_DG(i,j,ispec,:)
        
        chosen_nxnz_forMPI = -1
      
        ! Reinit boolean to know if neighbor exists
        neighbor_exists = .false.
        
        !nx = normal_DG(ispec, ind, 1)
        !nz = normal_DG(ispec, ind, 2)
        nx     = normal_DG(i,j,ispec, 1)
        nz     = normal_DG(i,j,ispec, 2)
        weight = weight_DG(i, j, ispec)
        dir_normal = dir_normal_DG(i,j,ispec)
        chosen_nxnz_forMPI = 0
        
        ! Needs x2 points at corners to correctly map edges
        ! => 2 step on the same point
        if(it_corner == 1) then
                neighbor = neighbor_DG_corner(i,j,ispec,:)
                nx     = normal_DG_corner(i,j,ispec, 1)
                nz     = normal_DG_corner(i,j,ispec, 2)
                weight = weight_DG_corner(i, j, ispec)
                dir_normal = dir_normal_DG_corner(i,j,ispec)
                chosen_nxnz_forMPI = 1
                it_corner = 2
        endif
        
        if(neighbor_DG(i,j,ispec,3) > -1 .OR. &
                neighbor_DG_corner(i,j,ispec,3) > -1) neighbor_exists = .true.
        
        ! If not outer boundary check for corresponding neighbor normal
        if(is_corner(i,j)) then
                
                ! If at least one neighbor exists
                if(neighbor_exists) then
                
                i_ex = neighbor(1)
                j_ex = neighbor(2)
                ispec_ex = neighbor(3)
                ! If corner of an outside edge
                if(it_corner == 2 .AND. neighbor(3) == -1) then
                        i_ex = neighbor_DG(i,j,ispec,1)
                        j_ex = neighbor_DG(i,j,ispec,2)
                        ispec_ex = neighbor_DG(i,j,ispec,3)
                elseif(it_corner < 2 .AND. neighbor(3) == -1) then
                        i_ex = neighbor_DG_corner(i,j,ispec,1)
                        j_ex = neighbor_DG_corner(i,j,ispec,2)
                        ispec_ex = neighbor_DG_corner(i,j,ispec,3)
                endif
                
                ! Cross product to verify if the normal corresponds to the normal
                !normal_DG(i_ex,j_ex,ispec_ex, 1) 
                !normal_DG(i_ex,j_ex,ispec_ex, 2)
                if( dir_normal /= -dir_normal_DG(i_ex,j_ex,ispec_ex) .AND. &
                     dir_normal /= -dir_normal_DG_corner(i_ex,j_ex,ispec_ex) ) then
                        ! Only change normal if inner element
                        if(neighbor(3) > -1 .AND. it_corner < 2) then
                                nx     = normal_DG_corner(i,j,ispec, 1)
                                nz     = normal_DG_corner(i,j,ispec, 2)
                                weight = weight_DG_corner(i, j, ispec)
                                dir_normal = dir_normal_DG_corner(i,j,ispec)
                                ! MODIF for MPI
                                chosen_nxnz_forMPI = 1
                        elseif(neighbor(3) > -1 .AND. it_corner == 2) then
                                nx     = normal_DG(i,j,ispec, 1)
                                nz     = normal_DG(i,j,ispec, 2)
                                weight = weight_DG(i, j, ispec)
                                dir_normal = dir_normal_DG(i,j,ispec)
                                ! MODIF for MPI
                                chosen_nxnz_forMPI = 0
                        endif
                ! If outside element, if the normal corresponds to the one computed here
                ! it means that we should take the other one
                elseif(neighbor(3) == -1) then
                        ! Only change normal if inner element
                        if(it_corner < 2) then
                                nx     = normal_DG_corner(i,j,ispec, 1)
                                nz     = normal_DG_corner(i,j,ispec, 2)
                                weight = weight_DG_corner(i, j, ispec)
                                dir_normal = dir_normal_DG_corner(i,j,ispec)
                                ! MODIF for MPI
                                chosen_nxnz_forMPI = 1
                        elseif(it_corner == 2) then
                                nx     = normal_DG(i,j,ispec, 1)
                                nz     = normal_DG(i,j,ispec, 2)
                                weight = weight_DG(i, j, ispec)
                                dir_normal = dir_normal_DG(i,j,ispec)
                                ! MODIF for MPI
                                chosen_nxnz_forMPI = 0
                        endif
                endif
                
                endif
                
        endif
        
        ! Interior point
        iglobM = ibool_DG(i,j,ispec)
        
        ! If a MPI surface node has been ill referenced and we need to witch between
        ! normal_DG and normal_DG_corner
        if(MPI_change_cpt) then
        
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
                
        endif !if(is_MPI_interface_DG(iglobM) .AND. NPROC > 1)
        
        ! If at corner notify that we will need to go again 
        if( is_corner(i,j) .AND. it_corner == 0) it_corner = 1
        
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

        ! compute dot product
        flux_x = T(iglobM) + T_P
        flux_n = flux_x*nx
        grad_Tx(iglobM) = grad_Tx(iglobM) + weight*flux_n*HALF
        
        ! compute dot product
        flux_z = T(iglobM) + T_P
        flux_n = flux_z*nz
        grad_Tz(iglobM) = grad_Tz(iglobM) + weight*flux_n*HALF
        
        ! compute dot product
        flux_x = veloc_x_DG(iglobM) + veloc_x_DG_P
        flux_n = flux_x*nx
        grad_Vxx(iglobM) = grad_Vxx(iglobM) + weight*flux_n*HALF
        
        ! compute dot product
        flux_z = veloc_x_DG(iglobM) + veloc_x_DG_P
        flux_n = flux_z*nz
        grad_Vxz(iglobM) = grad_Vxz(iglobM) + weight*flux_n*HALF
        
        ! compute dot product
        flux_x = veloc_z_DG(iglobM) + veloc_z_DG_P
        flux_n = flux_x*nx
        grad_Vzx(iglobM) = grad_Vzx(iglobM) + weight*flux_n*HALF
        
        ! compute dot product
        flux_z = veloc_z_DG(iglobM) + veloc_z_DG_P
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
      
    endif ! end of test if acoustic element
    
  enddo
  
  !WRITE(*,*) i_stage,"---------->", minval(grad_Tz), maxval(grad_Tz), timelocal
  if(.true.) then
  grad_Tx(:)  = grad_Tx(:) * rmass_inverse_acoustic_DG(:)
  grad_Tz(:)  = grad_Tz(:) * rmass_inverse_acoustic_DG(:)
  grad_Vxx(:) = grad_Vxx(:) * rmass_inverse_acoustic_DG(:)
  grad_Vxz(:) = grad_Vxz(:) * rmass_inverse_acoustic_DG(:)
  grad_Vzx(:) = grad_Vzx(:) * rmass_inverse_acoustic_DG(:)
  grad_Vzz(:) = grad_Vzz(:) * rmass_inverse_acoustic_DG(:)
  endif
  
  T_DG(1,:) = grad_Tx
  T_DG(2,:) = grad_Tz
  V_DG(1,1,:) = grad_Vxx
  V_DG(1,2,:) = grad_Vxz
  V_DG(2,1,:) = grad_Vzx
  V_DG(2,2,:) = grad_Vzz
  
  end subroutine compute_viscous_tensors
 
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
                         buffer_DG_gamma_P, myrank!, coord, ibool
                         
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
                        
        elseif(ACOUSTIC_FORCING .AND. ispec_is_acoustic_forcing(i,j,ispec)) then
        
                stop 'ACOUSTIC_FORCING obsolete for DG simulations'
                
        ! Elastic coupling        
        elseif(ispec_is_acoustic_coupling_el(i,j,ispec,3) >= 0) then
        
               ! If we already know the "real" flux at boundary
               !exact_interface_flux = .false.
        
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
               normal_v     = veloc_x*nx + veloc_z*nz
               tangential_v = veloc_x_DG_P*tx + veloc_z_DG_P*tz
               
               ! Transformation matrix between mesh coordinates and normal/tangential coordinates
               trans_boundary(1,1) = tz
               trans_boundary(1,2) = -nz
               trans_boundary(2,1) = -tx
               trans_boundary(2,2) = nx
               trans_boundary = trans_boundary/(nx*tz - tx*nz)
               
               ! From free slip and normal velocity continuity
               veloc_x_DG_P = trans_boundary(1,1)*normal_v + trans_boundary(1,2)*tangential_v!veloc_elastic(1,iglob)
               veloc_z_DG_P = trans_boundary(2,1)*normal_v + trans_boundary(2,2)*tangential_v
               
               !tensor_temp = 0
               ! Maybe we should remove the initial velocity field to veloc_x_DG_P/veloc_z_DG_P in the following
               ! Since we should only have the perturbated non-linear stress tensor
               tensor_temp(1,1) = -elastic_tensor(i_el,j_el,ispec_el,1) - rho_DG_P*veloc_x_DG_P**2 
               tensor_temp(1,2) = -elastic_tensor(i_el,j_el,ispec_el,2) - rho_DG_P*veloc_x_DG_P*veloc_z_DG_P
               tensor_temp(2,1) = -elastic_tensor(i_el,j_el,ispec_el,3) - rho_DG_P*veloc_x_DG_P*veloc_z_DG_P
               tensor_temp(2,2) = -elastic_tensor(i_el,j_el,ispec_el,4) - rho_DG_P*veloc_z_DG_P**2 
               
               ! From traction continuity
               p_DG_P = p_DG_init(iglobM) - (&
                  nx*( nx*tensor_temp(1,1) + nz*tensor_temp(1,2) ) &
                + nz*( nx*tensor_temp(2,1) + nz*tensor_temp(2,2) ) )
               
               
               E_DG_P       = p_DG_P/(gammaext_DG(iglobM) - 1.) + HALF*rho_DG_P*( veloc_x_DG_P**2 + veloc_z_DG_P**2 )
               rhovx_DG_P   = rho_DG_P*veloc_x_DG_P
               rhovz_DG_P   = rho_DG_P*veloc_z_DG_P
               
        ! Classical boundary conditions
        else
        
                call boundary_condition_DG(i, j, ispec, timelocal, rho_DG_P, rhovx_DG_P, rhovz_DG_P, E_DG_P, &
                veloc_x_DG_P, veloc_z_DG_P, p_DG_P, e1_DG_P)
                
                tx = -nz
                tz = nx
               ! Normal velocity of the solid perturbation and tangential velocity of the fluid flow
                normal_v     = veloc_x_DG_P*nx + veloc_z_DG_P*nz
                normal_v     = 2*normal_v-(veloc_x_DG_iM*nx + veloc_z_DG_iM*nz)
                tangential_v = -(veloc_x_DG_iM*tx + veloc_z_DG_iM*tz)
                
               ! Transformation matrix between mesh coordinates and normal/tangential coordinates
                trans_boundary(1,1) = tz
                trans_boundary(1,2) = -nz
                trans_boundary(2,1) = -tx
                trans_boundary(2,2) = nx
                trans_boundary = trans_boundary/(nx*tz - tx*nz)
               
               ! From free slip and normal velocity continuity
               veloc_x_DG_P = trans_boundary(1,1)*normal_v + trans_boundary(1,2)*tangential_v!veloc_elastic(1,iglob)
               veloc_z_DG_P = trans_boundary(2,1)*normal_v + trans_boundary(2,2)*tangential_v
                
                E_DG_P       = p_DG_P/(gammaext_DG(iglobM) - ONE) &
                        + rho_DG_P*HALF*( veloc_x_DG_P**2 + veloc_z_DG_P**2 ) 
                
                rhovx_DG_P   = rho_DG_P*veloc_x_DG_P
                rhovz_DG_P   = rho_DG_P*veloc_z_DG_P
                
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

  endif
  
  ! Temperature computation
  T_P = (E_DG_P/rho_DG_P - 0.5*((rhovx_DG_P/rho_DG_P)**2 + (rhovz_DG_P/rho_DG_P)**2))/(cnu)
  
  end subroutine compute_interface_unknowns


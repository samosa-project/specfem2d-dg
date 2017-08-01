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


  !subroutine compute_forces_acoustic_backward_DG(b_rho_DG, b_rhovx_DG, b_rhovz_DG, b_E_DG, &
  !      rho_DG, rhovx_DG, rhovz_DG, E_DG, b_dot_rho, b_dot_rhovx, b_dot_rhovz, b_dot_E)
  subroutine compute_forces_acoustic_backward_DG(it_temp, istage_temp, b_rho_DG, b_rhovx_DG, b_rhovz_DG, b_E_DG, &
        b_dot_rho, b_dot_rhovx, b_dot_rhovz, b_dot_E)

! compute forces in the acoustic elements in forward simulation and in adjoint simulation in adjoint inversion

  use constants,only: CUSTOM_REAL,NGLLX,NGLLZ,PI

  use specfem_par, only: nglob,nspec, &
                         ibool,ibool_DG,ispec_is_acoustic, &
                         xix,xiz,gammax,gammaz,jacobian, &
                         myrank, coord, hprimewgll_xx, hprimewgll_zz, &
                         !vpext,rhoext, &
                         hprime_xx, hprime_zz,wxgll,wzgll!, &
                         !rhoext, windxext, gammaext_DG, gravityext, pext!, assign_external_model, 

  implicit none

  real(kind=CUSTOM_REAL), dimension(nglob) :: b_rho_DG, b_rhovx_DG, b_rhovz_DG, b_E_DG
  !real(kind=CUSTOM_REAL), dimension(nglob_DG) :: rho_DG, rhovx_DG, rhovz_DG, E_DG, &
  real(kind=CUSTOM_REAL), dimension(nglob) :: b_dot_rho, b_dot_rhovx, b_dot_rhovz, b_dot_E

  ! local parameters
  integer :: ispec,i,j,k,iglob,iglob_DG
  integer :: it_temp, istage_temp
  integer :: ifirstelem,ilastelem

  real(kind=CUSTOM_REAL), dimension(NGLLX, NGLLZ) :: temp_b_rho, temp_b_E, temp_b_rhovx, temp_b_rhovz, &
        temp_b_rho_1, temp_b_rho_2, temp_b_p_1, temp_b_p_2, &
        temp_b_vx_1, temp_b_vx_2, temp_b_vz_1, temp_b_vz_2

  ! Jacobian matrix and determinant
  real(kind=CUSTOM_REAL) :: xixl,xizl,gammaxl,gammazl,jacobianl

  ! material properties of the acoustic medium
  double precision :: gammal, gravityl!, c_adiab_l, p0l, v0l,rho0l, 
  
  double precision :: drho_dxi, drho_dgamma, drhovx_dxi, drhovx_dgamma, &
        drhovz_dxi, drhovz_dgamma, dE_dxi, dE_dgamma
        
  double precision :: b_drhovx_dxl, b_drhovx_dzl, b_drhovz_dxl, b_drhovz_dzl, b_dE_dxl, &
        b_dE_dzl, b_drho_dxl, b_drho_dzl, Div_b_rhov, &
        b_rhovx_l, b_rhovz_l, b_E_l, b_rho_l, rhovx_l, &
        rhovz_l, E_l, rho_l, p_l, wzl, wxl, &
        H, veloc_x_l, veloc_z_l, z!, dtrho_l, vx_l, vz_l, 

  double precision :: drho2x_dxi, drho2z_dxi, drho2x_dgamma, drho2z_dgamma, Div_b_rhov_2, &
        dux_dxl, dux_dzl, &
        rho0_l, p0_l, E0_l, &
        mx_l, mx0_l, mz_l, mz0_l, vx_l, vz_l, cpl

  double precision, dimension(2) :: Sigma_drho_x, Sigma_drho_z, Sigma_drhovx_x, Sigma_drhovx_z, &
          Sigma_drhovz_x, Sigma_drhovz_z
  
  !character(len=100) file_name
  
  if(myrank == 5) &
  WRITE(*,*) it_temp, istage_temp, ">>>> 1 >>>>", minval(b_E_DG), maxval(b_E_DG), minval(b_rhovz_DG), maxval(b_rhovz_DG), &
       minval(b_rhovx_DG), maxval(b_rhovx_DG), minval(b_rho_DG), maxval(b_rho_DG)
  
  !if(myrank == 5) &
  !WRITE(*,*) it_temp, istage_temp, ">>>> 2 >>>>", minval(b_dot_E), maxval(b_dot_E), minval(b_dot_rhovz), maxval(b_dot_rhovz), &
  !     minval(b_dot_rhovx), maxval(b_dot_rhovx), minval(b_dot_rho), maxval(b_dot_rho)
  
  ifirstelem = 1
  ilastelem = nspec

! loop over spectral elements
  do ispec = ifirstelem,ilastelem

    ! acoustic spectral element
    if (ispec_is_acoustic(ispec)) then

      ! first double loop over GLL points to compute and store gradients
      do j = 1,NGLLZ
        do i = 1,NGLLX
        
          iglob_DG = ibool_DG(i,j,ispec)
          iglob    = ibool(i,j,ispec)
        
          z = real(coord(2,ibool(i,j,ispec)), kind=CUSTOM_REAL)
          H = 15000
          gammal = 1.4
          gravityl = 9.81
          
          cpl = 340
          rho_l     = exp(-z/H)
          rho_l     = 1.
          !rho_l     = exp(-3750/H)
          veloc_x_l = 0.!50*sin(z*PI/7500)!0.
          veloc_z_l = 0.
          rhovx_l = rho_l*veloc_x_l
          rhovz_l = rho_l*veloc_z_l
          
          p_l     = rho_l*gravityl*H
          !p_l     = (gammal*gravityl*H)*rho_l/gammal
          
          E_l     = p_l/(gammal-1) + (rho_l/2.)*( veloc_x_l**2 + veloc_z_l**2 )
        
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! NO IPP ADJOINT NAVIER STOKES
          if(.false.) then
        
          gravityl = -gravityl
        
          ! For latter
          !muext(i,j,ispec)  = mu_visco
          !etaext(i,j,ispec) = eta_visco
          !kappa_DG(i,j,ispec) = 0.2!10
        
          ! derivative along x and along z
          drho_dxi = 0._CUSTOM_REAL; drho_dgamma = 0._CUSTOM_REAL
          drhovx_dxi = 0._CUSTOM_REAL; drhovx_dgamma = 0._CUSTOM_REAL
          drhovz_dxi = 0._CUSTOM_REAL; drhovz_dgamma = 0._CUSTOM_REAL
          dE_dxi = 0._CUSTOM_REAL; dE_dgamma = 0._CUSTOM_REAL

          drho2x_dxi = 0._CUSTOM_REAL; drho2z_dxi = 0._CUSTOM_REAL
          drho2x_dgamma = 0._CUSTOM_REAL; drho2z_dgamma = 0._CUSTOM_REAL

          ! first double loop over GLL points to compute and store gradients
          ! we can merge the two loops because NGLLX == NGLLZ
          do k = 1,NGLLX
            !drhovx_dxi    = drhovx_dxi    + b_rhovx_DG(ibool_DG(k,j,ispec)) * hprime_xx(i,k)
            !drhovx_dgamma = drhovx_dgamma + b_rhovx_DG(ibool_DG(i,k,ispec)) * hprime_zz(j,k)
            !drhovz_dxi    = drhovz_dxi    + b_rhovz_DG(ibool_DG(k,j,ispec)) * hprime_xx(i,k)
            !drhovz_dgamma = drhovz_dgamma + b_rhovz_DG(ibool_DG(i,k,ispec)) * hprime_zz(j,k)
            !dE_dxi        = dE_dxi        + b_E_DG(ibool_DG(k,j,ispec)) * hprime_xx(i,k)
            !dE_dgamma     = dE_dgamma     + b_E_DG(ibool_DG(i,k,ispec)) * hprime_zz(j,k)
            !drho_dxi      = drho_dxi      + b_rho_DG(ibool_DG(k,j,ispec)) * hprime_xx(i,k)
            !drho_dgamma   = drho_dgamma   + b_rho_DG(ibool_DG(i,k,ispec)) * hprime_zz(j,k)
            
            drhovx_dxi    = drhovx_dxi    + b_rhovx_DG(ibool(k,j,ispec)) * hprime_xx(i,k)
            drhovx_dgamma = drhovx_dgamma + b_rhovx_DG(ibool(i,k,ispec)) * hprime_zz(j,k)
            
            drhovz_dxi    = drhovz_dxi    + b_rhovz_DG(ibool(k,j,ispec)) * hprime_xx(i,k)
            drhovz_dgamma = drhovz_dgamma + b_rhovz_DG(ibool(i,k,ispec)) * hprime_zz(j,k)
            
            dE_dxi        = dE_dxi        + b_E_DG(ibool(k,j,ispec)) * hprime_xx(i,k)
            dE_dgamma     = dE_dgamma     + b_E_DG(ibool(i,k,ispec)) * hprime_zz(j,k)
            
            drho_dxi      = drho_dxi      + b_rho_DG(ibool(k,j,ispec)) * hprime_xx(i,k)
            drho_dgamma   = drho_dgamma   + b_rho_DG(ibool(i,k,ispec)) * hprime_zz(j,k)
            
            z = real(coord(2,ibool(k,j,ispec)), kind=CUSTOM_REAL)
            rho_l     = exp(-z/H)
            drho2x_dxi      = drho2x_dxi      + rho_l*b_rhovx_DG(ibool(k,j,ispec)) * hprime_xx(i,k)
            drho2z_dxi      = drho2z_dxi      + rho_l*b_rhovz_DG(ibool(k,j,ispec)) * hprime_xx(i,k)
            
            z = real(coord(2,ibool(i,k,ispec)), kind=CUSTOM_REAL)
            rho_l     = exp(-z/H)
            drho2x_dgamma   = drho2x_dgamma   + rho_l*b_rhovx_DG(ibool(i,k,ispec)) * hprime_zz(j,k)
            drho2z_dgamma   = drho2z_dgamma   + rho_l*b_rhovz_DG(ibool(i,k,ispec)) * hprime_zz(j,k)
          enddo

          xixl = xix(i,j,ispec)
          xizl = xiz(i,j,ispec)
          gammaxl = gammax(i,j,ispec)
          gammazl = gammaz(i,j,ispec)

          ! derivatives of backward variables
          b_drhovx_dxl = drhovx_dxi * xixl + drhovx_dgamma * gammaxl
          b_drhovx_dzl = drhovx_dxi * xizl + drhovx_dgamma * gammazl
          
          b_drhovz_dxl = drhovz_dxi * xixl + drhovz_dgamma * gammaxl
          b_drhovz_dzl = drhovz_dxi * xizl + drhovz_dgamma * gammazl
          
          b_dE_dxl = dE_dxi * xixl + dE_dgamma * gammaxl
          b_dE_dzl = dE_dxi * xizl + dE_dgamma * gammazl
          
          b_drho_dxl = drho_dxi * xixl + drho_dgamma * gammaxl
          b_drho_dzl = drho_dxi * xizl + drho_dgamma * gammazl
          
          b_drho_dxl = (drho2x_dxi) * xixl + (drho2x_dgamma) * gammaxl
          b_drho_dzl = (drho2z_dxi) * xizl + (drho2z_dgamma) * gammazl

          Div_b_rhov = b_drhovx_dxl + b_drhovz_dzl
          
          ! Temp
          Div_b_rhov_2 = b_drho_dxl + b_drho_dzl

          ! Backward variables
          !b_rhovx_l = b_rhovx_DG(iglob_DG)
          !b_rhovz_l = b_rhovz_DG(iglob_DG)
          !b_E_l     = b_E_DG(iglob_DG)
          !b_rho_l   = b_rho_DG(iglob_DG)
          b_rhovx_l = b_rhovx_DG(iglob)
          b_rhovz_l = b_rhovz_DG(iglob)
          b_E_l     = b_E_DG(iglob)
          b_rho_l   = b_rho_DG(iglob)
          
          ! Forward variables
          !rhovx_l = rhovx_DG(iglob_DG)
          !rhovz_l = rhovz_DG(iglob_DG)
          !E_l     = E_DG(iglob_DG)
          !rho_l   = rho_DG(iglob_DG)
          
          ! Forward stress tensor component
          Sigma_drho_x(1) = -(1/rho_l**2)*rhovx_l**2
          Sigma_drho_x(2) = -(1/rho_l**2)*rhovx_l*rhovz_l
          Sigma_drho_z(1) = -(1/rho_l**2)*rhovx_l*rhovz_l 
          Sigma_drho_z(2) = -(1/rho_l**2)*rhovz_l**2
          
          Sigma_drhovx_x(1) = (1/rho_l)*2*rhovx_l
          Sigma_drhovx_x(2) = (1/rho_l)*rhovz_l
          Sigma_drhovx_z(1) = (1/rho_l)*rhovz_l 
          Sigma_drhovx_z(2) = 0
          
          Sigma_drhovz_x(1) = 0
          Sigma_drhovz_x(2) = (1/rho_l)*rhovx_l
          Sigma_drhovz_z(1) = (1/rho_l)*rhovx_l 
          Sigma_drhovz_z(2) = (1/rho_l)*2*rhovz_l
          
          jacobianl = jacobian(i,j,ispec)
           
           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           !  MOMENT FORMULATION
           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           
           temp_b_rho(i,j) = jacobianl * ( &
                - gravityl*rhovz_l - (Sigma_drho_x(1)*b_drhovx_dxl + Sigma_drho_x(2)*b_drhovx_dzl) &
                - (Sigma_drho_z(1)*b_drhovz_dxl + Sigma_drho_z(2)*b_drhovz_dzl) &
                - ((gammal-1)/(2*rho_l**2))*(rhovx_l**2 + rhovz_l**2)*Div_b_rhov &
                - ((gammal-1)/(2*rho_l**2))*(rhovx_l**2 + rhovz_l**2)*(1/rho_l)*(rhovx_l*b_dE_dxl &
                        + rhovz_l*b_dE_dzl) &
                + (E_l + p_l)*(1/rho_l**2)*(rhovx_l*b_dE_dxl + rhovz_l*b_dE_dzl) &
                )
           
           temp_b_E(i,j) = jacobianl * ( &
                -(gammal - 1.)*Div_b_rhov - (gammal/rho_l)*(rhovx_l*b_dE_dxl + rhovz_l*b_dE_dzl) &
           )
           
           temp_b_rhovx(i,j) = jacobianl * ( &
                - b_drho_dxl - (Sigma_drhovx_x(1)*b_drhovx_dxl + Sigma_drhovx_x(2)*b_drhovx_dzl) &
                - (Sigma_drhovx_z(1)*b_drhovz_dxl + Sigma_drhovx_z(2)*b_drhovz_dzl) &
                + ((gammal-1)*rhovx_l/rho_l)*Div_b_rhov &
                + ((gammal-1)*(rhovx_l)/(rho_l**2))*(rhovx_l*b_dE_dxl + rhovz_l*b_dE_dzl) &
                - (E_l + p_l)*(1/rho_l)*b_dE_dxl &
           )
           
           temp_b_rhovz(i,j) = jacobianl * ( &
                - b_drho_dzl - (Sigma_drhovz_x(1)*b_drhovx_dxl + Sigma_drhovz_x(2)*b_drhovx_dzl) &
                - (Sigma_drhovz_z(1)*b_drhovz_dxl + Sigma_drhovz_z(2)*b_drhovz_dzl) &
                + ((gammal-1)*rhovz_l/rho_l)*Div_b_rhov &
                + ((gammal-1)*(rhovz_l/(rho_l**2)))*(rhovx_l*b_dE_dxl + rhovz_l*b_dE_dzl) &
                - (E_l + p_l)*(1/rho_l)*b_dE_dzl &
                - gravityl*b_E_l &
           )
           
           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           
           z = real(coord(2,ibool(i,j,ispec)), kind=CUSTOM_REAL)
          H = 15000
          gammal = 1.4
          gravityl = 9.81
          
          rho_l     = exp(-z/H)
          !rho_l     = exp(-3750/H)
          veloc_x_l = 0.
          veloc_z_l = 0.
          rhovx_l = rho_l*veloc_x_l
          rhovz_l = rho_l*veloc_z_l
          
          p_l     = rho_l*gravityl*H
          !p_l     = (gammal*gravityl*H)*rho_l/gammal
          
          E_l     = p_l/(gammal-1) + (rho_l/2.)*( veloc_x_l**2 + veloc_z_l**2 )
        
          gravityl = -gravityl
           
           temp_b_rho(i,j) = jacobianl * ( &
                - Div_b_rhov_2 &
                )
           
           temp_b_E(i,j) = jacobianl * ( &
                - p_l*gammal*Div_b_rhov - rho_l*b_rhovz_l*gravityl &
           )
           
           temp_b_rhovx(i,j) = jacobianl * (1/rho_l)*( &
                - b_dE_dxl &
           )
           
           temp_b_rhovz(i,j) = jacobianl * (1/rho_l)*( &
                - b_dE_dzl + b_rho_l*gravityl &
           )
           
           endif
           
           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           ! IPP ADJOINT NAVIER STOKES
           if(.true.) then
           
           
           b_rhovx_l = b_rhovx_DG(iglob)
           b_rhovz_l = b_rhovz_DG(iglob)
           b_E_l     = b_E_DG(iglob)
           b_rho_l   = b_rho_DG(iglob)
          
           jacobianl = jacobian(i,j,ispec)
           
           H = 15000
           z = real(coord(2,ibool(i,j,ispec)), kind=CUSTOM_REAL)
           rho_l     = exp(-z/H)
           rho_l     = 1.
           gravityl = 0!9.81
           p_l   = rho_l*gammal*cpl**2
           !p_l   = rho_l*gravityl*H
           E_l   = p_l/(gammal-1) + (rho_l/2.)*( veloc_x_l**2 + veloc_z_l**2 )
           
           xixl = xix(i,j,ispec)
           xizl = xiz(i,j,ispec)
           gammaxl = gammax(i,j,ispec)
           gammazl = gammaz(i,j,ispec)
           
           !!!!!!!!!!!!!!!!!
           !!!! RHOP
           dux_dxl = -rho_l * (veloc_x_l**2) * b_rhovx_l &
                + ((gammal - 1.)*rho_l/2)*( veloc_x_l**2 ) * b_rhovx_l &
                + ((gammal - 1.)*rho_l/2)*( veloc_x_l**3 ) * b_E_l &
                - ((E_l + p_l)/rho_l) * veloc_x_l * b_E_l
           dux_dzl = 0. &
                + ((gammal - 1.)*rho_l/2)*( veloc_x_l**2 ) * b_rhovx_l &
                + 0. &!((gammal - 1.)*rho_l/2)*( veloc_x_l**3 ) * b_E_l &
                - 0. !((E_l + p_l)/rho_l) * veloc_x_l * b_E_l
           temp_b_rho_1(i,j) = wzgll(j) * jacobianl * ( &
                xixl * dux_dxl + xizl * dux_dzl)
           temp_b_rho_2(i,j) = wxgll(i) * jacobianl * ( &
                gammaxl * dux_dxl + gammazl * dux_dzl)
           
           ! Advection
           !drho_dxi = 0._CUSTOM_REAL
           !drho_dgamma = 0._CUSTOM_REAL
           !! first double loop over GLL points to compute and store gradients
           !! we can merge the two loops because NGLLX == NGLLZ
           !do k = 1,NGLLX
           ! drho_dxi      = drho_dxi      + b_rhovx_DG(ibool(k,j,ispec)) * hprime_xx(i,k)
           ! drho_dgamma   = drho_dgamma   + b_rhovx_DG(ibool(i,k,ispec)) * hprime_zz(j,k)
           !enddo
           !b_drho_dxl = (drho_dxi) * xixl + (drho_dgamma) * gammaxl
           !
           temp_b_rho(i,j) = jacobianl * ( &
                gravityl * b_rhovz_l &
           )
           
           !!!!!!!!!!!!!!!!!
           !!!! P
           drho_dxi = 0._CUSTOM_REAL
           drho_dgamma = 0._CUSTOM_REAL
           dE_dxi = 0.
           dE_dgamma = 0.
           ! first double loop over GLL points to compute and store gradients
           ! we can merge the two loops because NGLLX == NGLLZ
           do k = 1,NGLLX
            
            ! Gradxi rho*v^2
            z = real(coord(2,ibool(k,j,ispec)), kind=CUSTOM_REAL)
            rho_l = exp(-z/H)
            rho_l     = 1.
            p_l   = rho_l*gravityl*H
            gravityl = 0!9.81
            p_l   = rho_l*gammal*cpl**2
            veloc_x_l = 0.!50*sin(z*PI/7500)
            E_l   = p_l/(gammal-1) + (rho_l/2.)*( veloc_x_l**2 + veloc_z_l**2 )
            !p_l   = (340.**2)*rho_l/gammal
            !drho_dxi      = drho_dxi      + p_l*gammal * hprime_xx(i,k)
            dE_dxi      = dE_dxi + (gammal - 1.) * hprime_xx(i,k)
            
            ! Gradxi rho
            !dE_dxi        = dE_dxi        + (1/rho_l) * hprime_xx(i,k)
            drho_dxi    = drho_dxi   + ((E_l + p_l)/rho_l) * hprime_xx(i,k)
            
            ! Gradgamma rho*v^2
            z = real(coord(2,ibool(i,k,ispec)), kind=CUSTOM_REAL)
            rho_l = exp(-z/H)
            rho_l     = 1.
            veloc_x_l = 0.!50*sin(z*PI/7500)
            p_l   = rho_l*gravityl*H
            gravityl = 0!9.81
           p_l   = rho_l*gammal*cpl**2
            E_l   = p_l/(gammal-1) + (rho_l/2.)*( veloc_x_l**2 + veloc_z_l**2 )
            !p_l   = (340.**2)*rho_l/gammal
            !drho_dgamma   = drho_dgamma   + p_l*gammal * hprime_zz(j,k)
            dE_dgamma   = dE_dgamma   + (gammal - 1.) * hprime_zz(j,k)
            ! Gradgamma rho
            !dE_dgamma     = dE_dgamma     + (1/rho_l) * hprime_zz(j,k)
            drho_dgamma     = drho_dgamma     + ((E_l + p_l)/rho_l) * hprime_zz(j,k)
            
           enddo
           
           !b_drho_dxl = (drho_dxi) * xixl + (drho_dgamma) * gammaxl
           b_drho_dzl = (drho_dxi) * xizl + (drho_dgamma) * gammazl
           !b_dE_dxl = (dE_dxi) * xixl + (dE_dgamma) * gammaxl
           b_dE_dzl = (dE_dxi) * xizl + (dE_dgamma) * gammazl
           
           !!!!!!!!!!!!!!!!!!!!!!!!
           ! Background parameters
           H = 15000
           z = real(coord(2,ibool(i,j,ispec)), kind=CUSTOM_REAL)
           rho_l     = exp(-z/H)
           rho_l     = 1.
           veloc_x_l = 0.!50*sin(z*PI/7500)
           gravityl = 9.81
           p_l   = rho_l*gravityl*H
           gravityl = 0!9.81
           p_l   = rho_l*gammal*cpl**2
           E_l   = p_l/(gammal-1) + (rho_l/2.)*( veloc_x_l**2 + veloc_z_l**2 )
           
           temp_b_E(i,j) = jacobianl * ( &
                !(b_drho_dxl*b_rhovx_l + b_drho_dzl*b_rhovz_l) &
                !- rho_l*b_rhovz_l*(-gravityl) &
                b_rhovz_l * b_dE_dzl &
           )
           
           !WRITE(*,*) "b_dE_dzl", b_dE_dzl
           !stop 'b_dE_dzl'
           
           !dux_dxl = p_l*gammal * b_rhovx_l + b_E_l * veloc_x_l
           !dux_dzl = p_l*gammal * b_rhovz_l
           dux_dxl = (gammal - 1.) * b_rhovx_l &
                + gammal * veloc_x_l * b_E_l
           dux_dzl = (gammal - 1.) * b_rhovz_l &
                - 0.
           
           temp_b_p_1(i,j) = wzgll(j) * jacobianl * ( &
                xixl * dux_dxl + xizl * dux_dzl)
           temp_b_p_2(i,j) = wxgll(i) * jacobianl * ( &
                gammaxl * dux_dxl + gammazl * dux_dzl)
           
           !!!!!!!!!!!!!!!!!
           !!!! VX
           dux_dxl = ( b_rho_l + ((p_l + E_l)/rho_l) * b_E_l ) &
                + 2. * veloc_x_l * b_rhovx_l &
                - (gammal - 1.) * veloc_x_l * b_rhovx_l &
                - (gammal - 1.) * (veloc_x_l**2) * b_E_l
           !dux_dxl = b_E_l + rho_l*b_rhovx_l * veloc_x_l
           dux_dzl = 0. &!- b_rho_l + ((p_l + E_l)/rho_l) * b_E_l
                - 0. &
                - (gammal - 1.) * veloc_x_l * b_rhovz_l &
                - 0.
           temp_b_vx_1(i,j) = wzgll(j) * jacobianl * ( &
                xixl * dux_dxl + xizl * dux_dzl)
           temp_b_vx_2(i,j) = wxgll(i) * jacobianl * ( &
                gammaxl * dux_dxl + gammazl * dux_dzl)
           
           temp_b_rhovx(i,j) = jacobianl * ( &
                0 &!b_E_l*b_dE_dxl &
           )
           
           ! Advection
           !drho_dxi = 0._CUSTOM_REAL
           !drho_dgamma = 0._CUSTOM_REAL
           !! first double loop over GLL points to compute and store gradients
           !! we can merge the two loops because NGLLX == NGLLZ
           !do k = 1,NGLLX
           ! drho_dxi      = drho_dxi      + b_rhovx_DG(ibool(k,j,ispec)) * hprime_xx(i,k)
           ! drho_dgamma   = drho_dgamma   + b_rhovx_DG(ibool(i,k,ispec)) * hprime_zz(j,k)
           !enddo
           !b_drho_dxl = (drho_dxi) * xixl + (drho_dgamma) * gammaxl
          ! 
          ! temp_b_rhovx(i,j) = temp_b_rhovx(i,j) + jacobianl * ( &
          !      -rho_l*veloc_x_l*b_drho_dxl &
          ! )
           
           !!!!!!!!!!!!!!!!!
           !!!! VZ
           !dux_dxl = 0. + rho_l*b_rhovz_l * veloc_x_l
           dux_dxl = 0. &
                + veloc_x_l * b_rhovz_l!(1./rho_l)*b_E_l
           dux_dzl = (b_rho_l + ((p_l + E_l)/rho_l) * b_E_l ) &
                + veloc_x_l * b_rhovx_l!b_E_l
           temp_b_vz_1(i,j) = wzgll(j) * jacobianl * ( &
                xixl * dux_dxl + xizl * dux_dzl)
           temp_b_vz_2(i,j) = wxgll(i) * jacobianl * ( &
                gammaxl * dux_dxl + gammazl * dux_dzl)
           
           temp_b_rhovz(i,j) = jacobianl * ( &
                + gravityl * b_E_l &
                + b_drho_dzl * b_E_l &
                !b_E_l*b_dE_dzl &
               ! + (b_rho_l/rho_l)*(-gravityl) &
               !+ (b_rho_l)*(-gravityl) &
           )
           
           ! Advection
           !drho_dxi = 0._CUSTOM_REAL
           !drho_dgamma = 0._CUSTOM_REAL
           !! first double loop over GLL points to compute and store gradients
           !! we can merge the two loops because NGLLX == NGLLZ
           !do k = 1,NGLLX
           ! drho_dxi      = drho_dxi      + b_rhovz_DG(ibool(k,j,ispec)) * hprime_xx(i,k)
           ! drho_dgamma   = drho_dgamma   + b_rhovz_DG(ibool(i,k,ispec)) * hprime_zz(j,k)
           !enddo
           !b_drho_dxl = (drho_dxi) * xixl + (drho_dgamma) * gammaxl
           !
           !temp_b_rhovz(i,j) = temp_b_rhovz(i,j) + jacobianl * ( &
           !     -rho_l*veloc_x_l*b_drho_dxl &
           !)
           
           endif ! if(.false.) then
           
           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           ! COMMON DIRECT NAVIER STOKES
           if(.false.) then
           z = real(coord(2,ibool(i,j,ispec)), kind=CUSTOM_REAL)
           gravityl = 9.81
           H = 15000
           veloc_x_l = 0!10.
           gammal = 1.4
           ! RHO0
           rho0_l = exp(-z/H)
           ! P0
           p0_l   = rho0_l*gravityl*H
           
           rho0_l = 1.
           p0_l   = (340.**2)*rho0_l/gammal
           gravityl =0
           mx0_l = rho0_l*veloc_x_l
           mz0_l = 0.
           
           ! E0
           E0_l   = p0_l/(gammal - 1) + (rho0_l/2.)*( veloc_x_l**2 )
           
           if(it_temp == 1 .AND. istage_temp == 1) then
                b_rhovx_DG(iglob) = rho0_l*veloc_x_l
           endif
           
           
           b_rhovx_l = b_rhovx_DG(iglob)
           b_rhovz_l = b_rhovz_DG(iglob)
           b_E_l     = b_E_DG(iglob)
           b_rho_l   = b_rho_DG(iglob)
          
           rho_l = rho0_l + b_rho_l
           mx_l  = mx0_l + b_rhovx_l
           mz_l  = mz0_l + b_rhovz_l
           
           jacobianl = jacobian(i,j,ispec)
           
           xixl = xix(i,j,ispec)
           xizl = xiz(i,j,ispec)
           gammaxl = gammax(i,j,ispec)
           gammazl = gammaz(i,j,ispec)
           
           !!!!!!!!!!!!!!!!!
           !!!! RHOP
           dux_dxl = mx_l
           dux_dzl = mz_l
           temp_b_rho_1(i,j) = wzgll(j) * jacobianl * ( &
                xixl * dux_dxl + xizl * dux_dzl)
           temp_b_rho_2(i,j) = wxgll(i) * jacobianl * ( &
                gammaxl * dux_dxl + gammazl * dux_dzl)
           
           temp_b_rho(i,j) = 0.
           
           !!!!!!!!!!!!!!!!!
           !!!! E
           vx_l = (mx_l/(rho_l+b_rho_l))
           vz_l = (mz_l/(rho_l+b_rho_l))
           
           !E_l   = gammal * ( b_E_l + E_l  ) - (gammal - 1)*( ((rho_l + b_rho_l)/2.)*(veloc_x_l + ) )
           ! (E + P)
           E_l   = gammal * ( b_E_l + E0_l  ) &
                - (gammal - 1)*( (0.5/(rho_l+b_rho_l))*(mx_l**2 + mz_l**2) )
           p_l = (gammal - 1)*( b_E_l + E0_l - (0.5/(rho_l+b_rho_l))*(mx_l**2 + mz_l**2) )
           p_l = p_l - p0_l
           
           E_l = (p_l + E_l - E0_l)*veloc_x_l + E0_l*vx_l
           dux_dxl = E_l/rho_l
           
           E_l   = gammal * ( b_E_l + E0_l  ) &
                - (gammal - 1)*( (0.5/(rho_l+b_rho_l))*(mx_l**2 + mz_l**2) )
           E_l =  E0_l*vz_l
           dux_dzl = E_l/rho_l
           
           temp_b_p_1(i,j) = wzgll(j) * jacobianl * ( &
                xixl * dux_dxl + xizl * dux_dzl)
           temp_b_p_2(i,j) = wxgll(i) * jacobianl * ( &
                gammaxl * dux_dxl + gammazl * dux_dzl)
           
           temp_b_E(i,j) = jacobianl * ( &
                mz_l*(-gravityl) &
           )
           
           !!!!!!!!!!!!!!!!!
           !!!! VX
           E_l = b_E_l + E0_l
           p_l = (gammal - 1)*( E_l - (0.5/(rho_l+b_rho_l))*(mx_l**2 + mz_l**2) )
           p_l = p_l - p0_l
           
           dux_dxl = p_l + 0*(1./rho_l)*mx_l**2
           dux_dzl = 0*(1./rho_l)*mx_l*mz_l
           temp_b_vx_1(i,j) = wzgll(j) * jacobianl * ( &
                xixl * dux_dxl + xizl * dux_dzl)
           temp_b_vx_2(i,j) = wxgll(i) * jacobianl * ( &
                gammaxl * dux_dxl + gammazl * dux_dzl)
           
           temp_b_rhovx(i,j) = jacobianl * ( &
                0 &!b_E_l*b_dE_dxl &
           )
           
           
           !!!!!!!!!!!!!!!!!
           !!!! VZ
           dux_dxl = 0*(1./rho_l)*mx_l*mz_l
           dux_dzl = p_l + 0*(1./rho_l)*mz_l**2
           temp_b_vz_1(i,j) = wzgll(j) * jacobianl * ( &
                xixl * dux_dxl + xizl * dux_dzl)
           temp_b_vz_2(i,j) = wxgll(i) * jacobianl * ( &
                gammaxl * dux_dxl + gammazl * dux_dzl)
           
           temp_b_rhovz(i,j) = jacobianl * ( &
               ! b_E_l*b_dE_dzl &
               ! + (b_rho_l/rho_l)*(-gravityl) &
               (b_rho_l)*(-gravityl) &
           )
           endif
           
           !WRITE(*,*) it_temp,istage_temp,">>>>>>>>", maxval(temp_b_p_1), maxval(temp_b_p_2), &
           !     maxval(temp_b_rho_1), maxval(temp_b_rho_2), &
           !     rho0_l, p0_l, E0_l, E_l
           
        enddo
      enddo
!
! second double-loop over GLL to compute all the terms
!
      do j = 1,NGLLZ
        do i = 1,NGLLX
        
          iglob = ibool(i,j,ispec)
          
          ! along x direction and z direction
          ! and assemble the contributions
            do k = 1,NGLLX
              b_dot_rho(iglob) = b_dot_rho(iglob) + &
                        (temp_b_rho_1(k,j) * hprimewgll_xx(k,i) + &
                        temp_b_rho_2(i,k) * hprimewgll_zz(k,j))
                        
              b_dot_E(iglob) = b_dot_E(iglob) + &
                        (temp_b_p_1(k,j) * hprimewgll_xx(k,i) + &
                        temp_b_p_2(i,k) * hprimewgll_zz(k,j))
                        
              b_dot_rhovx(iglob) = b_dot_rhovx(iglob) + &
                        (temp_b_vx_1(k,j) * hprimewgll_xx(k,i) + &
                        temp_b_vx_2(i,k) * hprimewgll_zz(k,j))
                        
              b_dot_rhovz(iglob) = b_dot_rhovz(iglob) + &
                        (temp_b_vz_1(k,j) * hprimewgll_xx(k,i) + &
                        temp_b_vz_2(i,k) * hprimewgll_zz(k,j))
            enddo
            
            wzl = real(wzgll(j), kind=CUSTOM_REAL)
            wxl = real(wxgll(i), kind=CUSTOM_REAL)
            
            !b_dot_rho(iglob_DG)   = b_dot_rho(iglob_DG)   + temp_b_rho * wxl * wzl
            !b_dot_rhovx(iglob_DG) = b_dot_rhovx(iglob_DG) + temp_b_rhovx * wxl * wzl
            !b_dot_rhovz(iglob_DG) = b_dot_rhovz(iglob_DG) + temp_b_rhovz * wxl * wzl
            !b_dot_E(iglob_DG)     = b_dot_E(iglob_DG)     + temp_b_E * wxl * wzl
            b_dot_rho(iglob)   = b_dot_rho(iglob)   + temp_b_rho(i,j) * wxl * wzl
            b_dot_rhovx(iglob) = b_dot_rhovx(iglob) + temp_b_rhovx(i,j) * wxl * wzl
            b_dot_rhovz(iglob) = b_dot_rhovz(iglob) + temp_b_rhovz(i,j) * wxl * wzl
            b_dot_E(iglob)     = b_dot_E(iglob)     + temp_b_E(i,j) * wxl * wzl
            
            !if((abs(b_dot_rho(iglob))) > 0) &
            !WRITE(*,*) myrank,">>>",iglob, b_dot_rho(iglob), b_dot_rhovx(iglob), b_dot_rhovz(iglob), b_dot_E(iglob)
            
        enddo ! second loop over the GLL points
      enddo

    endif ! end of test if acoustic element
  enddo ! end of loop over all spectral elements
  
  !if(myrank == 5) &
  !WRITE(*,*) it_temp, istage_temp, ">>>> 2 bis >>>>", minval(b_dot_E), maxval(b_dot_E), minval(b_dot_rhovz), maxval(b_dot_rhovz), &
  !     minval(b_dot_rhovx), maxval(b_dot_rhovx), minval(b_dot_rho), maxval(b_dot_rho)

  end subroutine compute_forces_acoustic_backward_DG

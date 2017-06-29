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

  ! -----------------------------------------------------------------------------------------------
  
  ! Compute initial condition
  subroutine boundary_forcing_function_DG(timelocal, x, z, rho_DG_P, rhovx_DG_P, rhovz_DG_P, E_DG_P)

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
  
  real(kind=CUSTOM_REAL) :: x, z
  ! Forcing
  real(kind=CUSTOM_REAL) :: to, perio, lambdo, xo
  
  real(kind=CUSTOM_REAL) :: Tl, Tu, rho0, p0, RR
  
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
   
  end subroutine boundary_forcing_function_DG
  

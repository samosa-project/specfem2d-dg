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

  subroutine compute_forces_NS_main()

  use specfem_par
  use specfem_par_noise

  implicit none
  
  ! MODIF DG
  double precision, dimension(5) :: rk4a_d, rk4b_d, rk4c_d
  double precision :: timelocal

  if(it == 1 .AND. i_stage == 1) then
          resu_rhoveloc_x  = 0.
          resu_rhoveloc_z  = 0.
          resu_density    = 0.
          resu_rhoenergy  = 0.
          
          rhoveloc_acoustic(1,:) = 0.
          rhoveloc_acoustic(2,:) = 0.
          rhoenergy(:) = 0.
          density_p(:) = 0.
  endif

  rk4a_d(1) = 0d0
          rk4a_d(2) = -567301805773.0/1357537059087.0
          rk4a_d(3) = -2404267990393.0/2016746695238.0
          rk4a_d(4) = -3550918686646.0/2091501179385.0
          rk4a_d(5) = -1275806237668.0/842570457699.0
            
          rk4b_d(1) = 1432997174477.0/9575080441755.0 
          rk4b_d(2) = 5161836677717.0/13612068292357.0 
          rk4b_d(3) = 1720146321549.0/2090206949498.0 
          rk4b_d(4) = 3134564353537.0/4481467310338.0 
          rk4b_d(5) = 2277821191437.0/14882151754819.0
            
          rk4c_d(1) = 0d0
          rk4c_d(2) = 1432997174477.0/9575080441755.0 
          rk4c_d(3) = 2526269341429.0/6820363962896.0 
          rk4c_d(4) = 2006345519317.0/3224310063776.0 
          rk4c_d(5) = 2802321613138.0/2924317926251.0

  ! main solver for the acoustic elements
  timelocal = (it-1)*deltat + rk4c_d(i_stage)*deltat

  WRITE(*,*) "ITERATION", timelocal, minval(rhoveloc_acoustic(2,:)),&
         maxval(rhoveloc_acoustic(2,:))

  ! checks if anything to do in this slice
  if (.not. any_acoustic) return

  ! visco-elastic term
  call compute_forces_NS(dt_rhoveloc_acoustic, dt_density, dt_rhoenergy, &
        rhoveloc_acoustic, rhoenergy, density_p)

  call compute_add_sources_NS(dt_rhoveloc_acoustic,density_p,it,i_stage)

  ! add coupling with the acoustic side
  !if (coupled_acoustic_elastic) then
  !  call compute_coupling_viscoelastic_ac()
  !endif

  ! dt_rhoveloc_acoustic(2,100) = dt_rhoveloc_acoustic(2,100) + 0.00001*exp(-((timelocal - 0.1 )/0.1)**3)

  ! multiply by the inverse of the mass matrix and update velocity
  !! DK DK this should be vectorized
  dt_rhoveloc_acoustic(1,:) = dt_rhoveloc_acoustic(1,:) * rmass_inverse_NS(:)
  dt_rhoveloc_acoustic(2,:) = dt_rhoveloc_acoustic(2,:) * rmass_inverse_NS(:)
  dt_density(:)            = dt_density(:) * rmass_inverse_NS(:)
  dt_rhoenergy(:)          = dt_rhoenergy(:) * rmass_inverse_NS(:)

  if(time_stepping_scheme == 3) then
  
          ! RK5-low dissipation Update
          resu_rhoveloc_x = rk4a_d(i_stage)*resu_rhoveloc_x + deltat*dt_rhoveloc_acoustic(1,:)
          resu_rhoveloc_z = rk4a_d(i_stage)*resu_rhoveloc_z + deltat*dt_rhoveloc_acoustic(2,:)
          resu_density    = rk4a_d(i_stage)*resu_density    + deltat*dt_density(:)
          resu_rhoenergy  = rk4a_d(i_stage)*resu_rhoenergy  + deltat*dt_rhoenergy(:)
    
          rhoveloc_acoustic(1,:) = rhoveloc_acoustic(1,:) + rk4b_d(i_stage)*resu_rhoveloc_x
          rhoveloc_acoustic(2,:) = rhoveloc_acoustic(2,:) + rk4b_d(i_stage)*resu_rhoveloc_z
          rhoenergy(:) = rhoenergy(:) + rk4b_d(i_stage)*resu_density
          density_p(:) = density_p(:) + rk4b_d(i_stage)*resu_rhoenergy  
  
  endif

  end subroutine compute_forces_NS_main


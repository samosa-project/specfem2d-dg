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


  subroutine prepare_timerun_wavefields()

  use specfem_par

  implicit none

  ! local parameters
  integer :: nglob_acoustic_b,nglob_elastic_b,nglob_poroelastic_b
  integer :: nspec_acoustic_b,nspec_elastic_b,nspec_poroelastic_b
  integer :: ier

  ! displacement, velocity, acceleration and inverse of the mass matrix for elastic elements
  if (myrank == 0) then
    write(IMAIN,*) 'Preparing array allocations'
    call flush_IMAIN()
  endif

  !
  ! elastic domains
  !
  ! user info
  if (ELASTIC_SIMULATION) then
    if (myrank == 0) then
      write(IMAIN,*) '  arrays for elastic domains'
      call flush_IMAIN()
    endif
  endif

  ! sets global points in this slice
  if (any_elastic) then
    nglob_elastic = nglob
  else
    ! dummy allocate unused arrays with fictitious size
    nglob_elastic = 1
  endif

  allocate(displ_elastic(NDIM,nglob_elastic), &
           veloc_elastic(NDIM,nglob_elastic), &
           accel_elastic(NDIM,nglob_elastic),stat=ier)
  if (ier /= 0) stop 'Error allocating elastic wavefield arrays'
  
  ! MODIF DG
  allocate(resu_accel_x(nglob_elastic), resu_accel_z(nglob_elastic), &
        resu_veloc_x(nglob_elastic), resu_veloc_z(nglob_elastic))

  ! PML
  allocate(displ_elastic_old(NDIM,nglob_elastic),stat=ier)
  if (ier /= 0) stop 'Error allocating old elastic wavefield arrays'

  if (SIMULATION_TYPE == 3) then
    allocate(accel_elastic_adj_coupling(NDIM,nglob_elastic))
  endif

  allocate(rmass_inverse_elastic_one(nglob_elastic))
  allocate(rmass_inverse_elastic_three(nglob_elastic))

  if (time_stepping_scheme==2) then
    allocate(displ_elastic_LDDRK(NDIM,nglob_elastic), &
             veloc_elastic_LDDRK(NDIM,nglob_elastic), &
             veloc_elastic_LDDRK_temp(NDIM,nglob_elastic),stat=ier)
    if (ier /= 0) stop 'Error allocating elastic LDDRK wavefield arrays'
  endif

  if (time_stepping_scheme == 3) then
    allocate(accel_elastic_rk(NDIM,nglob_elastic,stage_time_scheme), &
             veloc_elastic_rk(NDIM,nglob_elastic,stage_time_scheme), &
             veloc_elastic_initial_rk(NDIM,nglob_elastic), &
             displ_elastic_initial_rk(NDIM,nglob_elastic),stat=ier)
    if (ier /= 0) stop 'Error allocating elastic RK wavefield arrays'
  endif

  ! extra array if adjoint and kernels calculation
  if (SIMULATION_TYPE == 3 .and. any_elastic) then
    nglob_elastic_b = nglob
    nspec_elastic_b = nspec
  else
    ! dummy allocations
    nglob_elastic_b = 1
    nspec_elastic_b = 1
  endif

  allocate(b_displ_elastic(NDIM,nglob_elastic_b), &
           b_veloc_elastic(NDIM,nglob_elastic_b), &
           b_accel_elastic(NDIM,nglob_elastic_b),stat=ier)
  if (ier /= 0) stop 'Error allocating elastic backward wavefield arrays'

  allocate(b_displ_elastic_old(NDIM,nglob_elastic_b))

  ! kernels
  ! on global nodes
  allocate(rho_k(nglob_elastic_b))
  allocate(mu_k(nglob_elastic_b))
  allocate(kappa_k(nglob_elastic_b))
  allocate(c11_k(nglob_elastic_b))
  allocate(c13_k(nglob_elastic_b))
  allocate(c15_k(nglob_elastic_b))
  allocate(c33_k(nglob_elastic_b))
  allocate(c35_k(nglob_elastic_b))
  allocate(c55_k(nglob_elastic_b))
  ! on local nodes
  allocate(rho_kl(NGLLX,NGLLZ,nspec_elastic_b))
  allocate(mu_kl(NGLLX,NGLLZ,nspec_elastic_b))
  allocate(kappa_kl(NGLLX,NGLLZ,nspec_elastic_b))
  allocate(rhop_kl(NGLLX,NGLLZ,nspec_elastic_b))
  allocate(alpha_kl(NGLLX,NGLLZ,nspec_elastic_b))
  allocate(beta_kl(NGLLX,NGLLZ,nspec_elastic_b))
  allocate(bulk_c_kl(NGLLX,NGLLZ,nspec_elastic_b))
  allocate(bulk_beta_kl(NGLLX,NGLLZ,nspec_elastic_b))
  allocate(c11_kl(NGLLX,NGLLZ,nspec_elastic_b))
  allocate(c13_kl(NGLLX,NGLLZ,nspec_elastic_b))
  allocate(c15_kl(NGLLX,NGLLZ,nspec_elastic_b))
  allocate(c33_kl(NGLLX,NGLLZ,nspec_elastic_b))
  allocate(c35_kl(NGLLX,NGLLZ,nspec_elastic_b))
  allocate(c55_kl(NGLLX,NGLLZ,nspec_elastic_b))
  if (APPROXIMATE_HESS_KL) then
    allocate(rhorho_el_hessian_final2(NGLLX,NGLLZ,nspec_elastic_b))
    allocate(rhorho_el_hessian_final1(NGLLX,NGLLZ,nspec_elastic_b))
  endif

  !
  ! poro-elastic domains
  !
  if (POROELASTIC_SIMULATION) then
    if (myrank == 0) then
      write(IMAIN,*) '  arrays for poroelastic domains'
      call flush_IMAIN()
    endif
  endif

  ! sets number of points for this slice
  if (any_poroelastic) then
    nglob_poroelastic = nglob
  else
    ! dummy allocate unused arrays with fictitious size
    nglob_poroelastic = 1
  endif
  allocate(displs_poroelastic(NDIM,nglob_poroelastic))
  allocate(displs_poroelastic_old(NDIM,nglob_poroelastic))
  allocate(velocs_poroelastic(NDIM,nglob_poroelastic))
  allocate(accels_poroelastic(NDIM,nglob_poroelastic))
  if (SIMULATION_TYPE == 3) then
    allocate(accels_poroelastic_adj_coupling(NDIM,nglob_poroelastic))
  endif
  allocate(rmass_s_inverse_poroelastic(nglob_poroelastic))

  allocate(displw_poroelastic(NDIM,nglob_poroelastic))
  allocate(velocw_poroelastic(NDIM,nglob_poroelastic))
  allocate(accelw_poroelastic(NDIM,nglob_poroelastic))
  if (SIMULATION_TYPE == 3) then
    allocate(accelw_poroelastic_adj_coupling(NDIM,nglob_poroelastic))
  endif
  allocate(rmass_w_inverse_poroelastic(nglob_poroelastic))

  if (time_stepping_scheme == 2) then
    allocate(displs_poroelastic_LDDRK(NDIM,nglob_poroelastic))
    allocate(velocs_poroelastic_LDDRK(NDIM,nglob_poroelastic))
    allocate(displw_poroelastic_LDDRK(NDIM,nglob_poroelastic))
    allocate(velocw_poroelastic_LDDRK(NDIM,nglob_poroelastic))
  endif

  if (time_stepping_scheme == 3) then
    allocate(accels_poroelastic_rk(NDIM,nglob_poroelastic,stage_time_scheme))
    allocate(velocs_poroelastic_rk(NDIM,nglob_poroelastic,stage_time_scheme))
    allocate(accelw_poroelastic_rk(NDIM,nglob_poroelastic,stage_time_scheme))
    allocate(velocw_poroelastic_rk(NDIM,nglob_poroelastic,stage_time_scheme))
    allocate(displs_poroelastic_initial_rk(NDIM,nglob_poroelastic))
    allocate(velocs_poroelastic_initial_rk(NDIM,nglob_poroelastic))
    allocate(displw_poroelastic_initial_rk(NDIM,nglob_poroelastic))
    allocate(velocw_poroelastic_initial_rk(NDIM,nglob_poroelastic))
  endif

  ! extra array if adjoint and kernels calculation
  if (SIMULATION_TYPE == 3 .and. any_poroelastic) then
    nglob_poroelastic_b = nglob
    nspec_poroelastic_b = nspec
  else
    ! dummy allocations
    nglob_poroelastic_b = 1
    nspec_poroelastic_b = 1
  endif
  allocate(b_displs_poroelastic(NDIM,nglob_poroelastic_b))
  allocate(b_velocs_poroelastic(NDIM,nglob_poroelastic_b))
  allocate(b_accels_poroelastic(NDIM,nglob_poroelastic_b))
  allocate(b_displw_poroelastic(NDIM,nglob_poroelastic_b))
  allocate(b_velocw_poroelastic(NDIM,nglob_poroelastic_b))
  allocate(b_accelw_poroelastic(NDIM,nglob_poroelastic_b))
  ! kernels
  ! on global nodes
  allocate(rhot_k(nglob_poroelastic_b))
  allocate(rhof_k(nglob_poroelastic_b))
  allocate(sm_k(nglob_poroelastic_b))
  allocate(eta_k(nglob_poroelastic_b))
  allocate(mufr_k(nglob_poroelastic_b))
  allocate(B_k(nglob_poroelastic_b))
  allocate(C_k(nglob_poroelastic_b))
  allocate(M_k(nglob_poroelastic_b))
  ! on local nodes
  allocate(rhot_kl(NGLLX,NGLLZ,nspec_poroelastic_b))
  allocate(rhof_kl(NGLLX,NGLLZ,nspec_poroelastic_b))
  allocate(sm_kl(NGLLX,NGLLZ,nspec_poroelastic_b))
  allocate(eta_kl(NGLLX,NGLLZ,nspec_poroelastic_b))
  allocate(mufr_kl(NGLLX,NGLLZ,nspec_poroelastic_b))
  allocate(B_kl(NGLLX,NGLLZ,nspec_poroelastic_b))
  allocate(C_kl(NGLLX,NGLLZ,nspec_poroelastic_b))
  allocate(M_kl(NGLLX,NGLLZ,nspec_poroelastic_b))

  allocate(phi_kl(NGLLX,NGLLZ,nspec_poroelastic_b))
  allocate(phib_kl(NGLLX,NGLLZ,nspec_poroelastic_b))
  allocate(mufrb_kl(NGLLX,NGLLZ,nspec_poroelastic_b))

  allocate(rhob_kl(NGLLX,NGLLZ,nspec_poroelastic_b))
  allocate(rhobb_kl(NGLLX,NGLLZ,nspec_poroelastic_b))
  allocate(rhofb_kl(NGLLX,NGLLZ,nspec_poroelastic_b))
  allocate(rhofbb_kl(NGLLX,NGLLZ,nspec_poroelastic_b))

  allocate(cpI_kl(NGLLX,NGLLZ,nspec_poroelastic_b))
  allocate(cpII_kl(NGLLX,NGLLZ,nspec_poroelastic_b))
  allocate(cs_kl(NGLLX,NGLLZ,nspec_poroelastic_b))

  allocate(ratio_kl(NGLLX,NGLLZ,nspec_poroelastic_b))

  if (any_poroelastic .and. any_elastic) then
    allocate(icount(nglob))
  else
    allocate(icount(1))
  endif

  if (COMPUTE_INTEGRATED_ENERGY_FIELD) then ! = int_0^t v^2 dt
    allocate(integrated_cinetic_energy_field(nspec),stat=ier)
    if (ier /= 0) stop 'Error allocating integrated_cinetic_energy_field array'
    integrated_cinetic_energy_field(:) = 0._CUSTOM_REAL
    allocate(max_cinetic_energy_field(nspec),stat=ier)
    if (ier /= 0) stop 'Error allocating max_cinetic_energy_field array'
    max_cinetic_energy_field(:) = 0._CUSTOM_REAL
    allocate(integrated_potential_energy_field(nspec),stat=ier)
    if (ier /= 0) stop 'Error allocating integrated_potential_energy_field array'
    integrated_potential_energy_field(:) = 0._CUSTOM_REAL
    allocate(max_potential_energy_field(nspec),stat=ier)
    if (ier /= 0) stop 'Error allocating max_potential_energy_field array'
    max_potential_energy_field(:) = 0._CUSTOM_REAL
    allocate(cinetic_effective_duration_field(nspec),stat=ier)
    if (ier /= 0) stop 'Error allocating cinetic_effective_duration_field array'
    cinetic_effective_duration_field(:) = 0._CUSTOM_REAL
    allocate(potential_effective_duration_field(nspec),stat=ier)
    if (ier /= 0) stop 'Error allocating potential_effective_duration_field array'
    potential_effective_duration_field(:) = 0._CUSTOM_REAL
  endif
  !
  ! acoustic domains
  !
  if (ACOUSTIC_SIMULATION) then
    if (myrank == 0) then
      write(IMAIN,*) '  arrays for acoustic domains'
      call flush_IMAIN()
    endif
  endif

  ! potential, its first and second derivative, and inverse of the mass matrix for acoustic elements
  if (any_acoustic) then
    nglob_acoustic = nglob
  else
    ! dummy allocate unused arrays with fictitious size
    nglob_acoustic = 1
  endif

  ! MODIF NS
  allocate(resu_rhoveloc_x(nglob_acoustic), resu_rhoveloc_z(nglob_acoustic), &
        resu_density(nglob_acoustic), resu_rhoenergy(nglob_acoustic))
        
  allocate(dt_rhoveloc_acoustic(NDIM,nglob_acoustic), &
           rhoveloc_acoustic(NDIM,nglob_acoustic), &
           density_p(nglob_acoustic), rhoenergy(nglob_acoustic), stat=ier)  
   
  allocate(dt_density(nglob_acoustic), &
           dt_rhoenergy(nglob_acoustic),stat=ier)      
           
  allocate(rmass_inverse_NS(nglob_acoustic)) 

  ! Modif DG
  allocate(dot_rho(nglob_DG), dot_rhovx(nglob_DG), dot_rhovz(nglob_DG), dot_E(nglob_DG))
  allocate(resu_rho(nglob_DG), resu_rhovx(nglob_DG), resu_rhovz(nglob_DG), resu_E(nglob_DG))
  !allocate(rho_init_rk(nglob_DG), rhovx_init_rk(nglob_DG), rhovz_init_rk(nglob_DG), E_init_rk(nglob_DG))
  !allocate(dot_rho_rk(4,nglob_DG), dot_rhovx_rk(4,nglob_DG), dot_rhovz_rk(4,nglob_DG), dot_E_rk(4,nglob_DG))
  allocate(rho_DG(nglob_DG), p_DG(nglob_DG), rhovx_DG(nglob_DG), rhovz_DG(nglob_DG), &
        veloc_x_DG(nglob_DG), veloc_z_DG(nglob_DG), E_DG(nglob_DG))
  allocate(p_DG_init(nglob_DG))
  allocate(rmass_inverse_acoustic_DG(nglob_DG))
  allocate(potential_dphi_dx_DG(nglob), potential_dphi_dz_DG(nglob))
  
  allocate(potential_acoustic(nglob_acoustic))
  allocate(potential_acoustic_old(nglob_acoustic))
  if (SIMULATION_TYPE == 3) then
    allocate(potential_acoustic_adj_coupling(nglob_acoustic))
  endif
  allocate(potential_dot_acoustic(nglob_acoustic))
  allocate(potential_dot_dot_acoustic(nglob_acoustic))
  allocate(rmass_inverse_acoustic(nglob_acoustic))
  if (time_stepping_scheme == 2) then
    allocate(potential_acoustic_LDDRK(nglob_acoustic))
    allocate(potential_dot_acoustic_LDDRK(nglob_acoustic))
    allocate(potential_dot_acoustic_temp(nglob_acoustic))
  endif

  if (time_stepping_scheme == 3) then
    allocate(potential_acoustic_init_rk(nglob_acoustic))
    allocate(potential_dot_acoustic_init_rk(nglob_acoustic))
    allocate(potential_dot_dot_acoustic_rk(nglob_acoustic,stage_time_scheme))
    allocate(potential_dot_acoustic_rk(nglob_acoustic,stage_time_scheme))
  endif

  if (SIMULATION_TYPE == 3 .and. any_acoustic) then
    nglob_acoustic_b = nglob
    nspec_acoustic_b = nspec
  else
    ! dummy array allocations
    ! allocates unused arrays with fictitious size
    nglob_acoustic_b = 1
    nspec_acoustic_b = 1
  endif
  allocate(b_potential_acoustic(nglob_acoustic_b))
  allocate(b_potential_acoustic_old(nglob_acoustic_b))
  allocate(b_potential_dot_acoustic(nglob_acoustic_b))
  allocate(b_potential_dot_dot_acoustic(nglob_acoustic_b))
  allocate(b_displ_ac(2,nglob_acoustic_b))
  allocate(b_accel_ac(2,nglob_acoustic_b))
  allocate(accel_ac(2,nglob_acoustic_b))
  ! kernels
  ! on global points
  allocate(rhol_ac_global(nglob_acoustic_b))
  allocate(kappal_ac_global(nglob_acoustic_b))
  ! on local points
  allocate(rho_ac_kl(NGLLX,NGLLZ,nspec_acoustic_b))
  allocate(kappa_ac_kl(NGLLX,NGLLZ,nspec_acoustic_b))
  allocate(rhop_ac_kl(NGLLX,NGLLZ,nspec_acoustic_b))
  allocate(alpha_ac_kl(NGLLX,NGLLZ,nspec_acoustic_b))
  if (APPROXIMATE_HESS_KL) then
    allocate(rhorho_ac_hessian_final2(NGLLX,NGLLZ,nspec_acoustic_b))
    allocate(rhorho_ac_hessian_final1(NGLLX,NGLLZ,nspec_acoustic_b))
  endif

  !
  ! gravito-acoustic domains
  !
  if (GRAVITOACOUSTIC_SIMULATION) then
    if (myrank == 0) then
      write(IMAIN,*) '  arrays for elastic domains'
      call flush_IMAIN()
    endif
  endif

  ! potential, its first and second derivative, and inverse of the mass matrix for gravitoacoustic elements
  if (any_gravitoacoustic) then
    nglob_gravitoacoustic = nglob
  else
    ! allocate unused arrays with fictitious size
    nglob_gravitoacoustic = 1
  endif

  allocate(potential_gravitoacoustic(nglob_gravitoacoustic))
  allocate(potential_dot_gravitoacoustic(nglob_gravitoacoustic))
  allocate(potential_dot_dot_gravitoacoustic(nglob_gravitoacoustic))
  allocate(rmass_inverse_gravitoacoustic(nglob_gravitoacoustic))
  allocate(potential_gravito(nglob_gravitoacoustic))
  allocate(potential_dot_gravito(nglob_gravitoacoustic))
  allocate(potential_dot_dot_gravito(nglob_gravitoacoustic))
  allocate(rmass_inverse_gravito(nglob_gravitoacoustic))

  ! synchronizes all processes
  call synchronize_all()

  ! initializes wavefields
  if (myrank == 0) then
    write(IMAIN,*) '  wavefield initialization'
    call flush_IMAIN()
  endif
  
  resu_rhoveloc_x      = 0._CUSTOM_REAL
  resu_rhoveloc_z      = 0._CUSTOM_REAL
  resu_density         = 0._CUSTOM_REAL
  resu_rhoenergy       = 0._CUSTOM_REAL
  dt_rhoveloc_acoustic = 0._CUSTOM_REAL
  rhoveloc_acoustic    = 0._CUSTOM_REAL
  dt_density           = 0._CUSTOM_REAL
  dt_rhoenergy         = 0._CUSTOM_REAL
  density_p            = 0._CUSTOM_REAL
  rhoenergy            = 0._CUSTOM_REAL

  ! initialize arrays to zero
  displ_elastic = 0._CUSTOM_REAL
  displ_elastic_old = 0._CUSTOM_REAL
  veloc_elastic = 0._CUSTOM_REAL
  accel_elastic = 0._CUSTOM_REAL
  
  ! MODIF DG
  resu_accel_x = 0._CUSTOM_REAL
  resu_accel_z = 0._CUSTOM_REAL
  resu_veloc_x = 0._CUSTOM_REAL
  resu_veloc_z = 0._CUSTOM_REAL

  if (SIMULATION_TYPE == 3 .and. any_elastic) then
    b_displ_elastic_old = 0._CUSTOM_REAL
    b_displ_elastic = 0._CUSTOM_REAL
    b_veloc_elastic = 0._CUSTOM_REAL
    b_accel_elastic = 0._CUSTOM_REAL
  endif

  if (time_stepping_scheme == 2) then
    displ_elastic_LDDRK = 0._CUSTOM_REAL
    veloc_elastic_LDDRK = 0._CUSTOM_REAL
    veloc_elastic_LDDRK_temp = 0._CUSTOM_REAL
  endif

  if (time_stepping_scheme == 3) then
    accel_elastic_rk = 0._CUSTOM_REAL
    veloc_elastic_rk = 0._CUSTOM_REAL
    veloc_elastic_initial_rk = 0._CUSTOM_REAL
    displ_elastic_initial_rk = 0._CUSTOM_REAL
  endif

  displs_poroelastic = 0._CUSTOM_REAL
  displs_poroelastic_old = 0._CUSTOM_REAL
  velocs_poroelastic = 0._CUSTOM_REAL
  accels_poroelastic = 0._CUSTOM_REAL
  displw_poroelastic = 0._CUSTOM_REAL
  velocw_poroelastic = 0._CUSTOM_REAL
  accelw_poroelastic = 0._CUSTOM_REAL

  if (time_stepping_scheme == 2) then
    displs_poroelastic_LDDRK = 0._CUSTOM_REAL
    velocs_poroelastic_LDDRK = 0._CUSTOM_REAL
    displw_poroelastic_LDDRK = 0._CUSTOM_REAL
    velocw_poroelastic_LDDRK = 0._CUSTOM_REAL
  endif

  if (time_stepping_scheme == 3) then
    accels_poroelastic_rk = 0._CUSTOM_REAL
    velocs_poroelastic_rk = 0._CUSTOM_REAL

    accelw_poroelastic_rk = 0._CUSTOM_REAL
    velocw_poroelastic_rk = 0._CUSTOM_REAL

    velocs_poroelastic_initial_rk = 0._CUSTOM_REAL
    displs_poroelastic_initial_rk = 0._CUSTOM_REAL

    velocw_poroelastic_initial_rk = 0._CUSTOM_REAL
    displw_poroelastic_initial_rk = 0._CUSTOM_REAL
  endif

  potential_acoustic = 0._CUSTOM_REAL
  potential_acoustic_old = 0._CUSTOM_REAL
  potential_dot_acoustic = 0._CUSTOM_REAL
  potential_dot_dot_acoustic = 0._CUSTOM_REAL

  if (time_stepping_scheme == 2) then
    potential_acoustic_LDDRK = 0._CUSTOM_REAL
    potential_dot_acoustic_LDDRK = 0._CUSTOM_REAL
    potential_dot_acoustic_temp = 0._CUSTOM_REAL
  endif

  if (time_stepping_scheme == 3) then
    potential_acoustic_init_rk = 0._CUSTOM_REAL
    potential_dot_acoustic_init_rk = 0._CUSTOM_REAL
    potential_dot_dot_acoustic_rk = 0._CUSTOM_REAL
    potential_dot_acoustic_rk = 0._CUSTOM_REAL
  endif

  potential_gravitoacoustic = 0._CUSTOM_REAL
  potential_dot_gravitoacoustic = 0._CUSTOM_REAL
  potential_dot_dot_gravitoacoustic = 0._CUSTOM_REAL
  potential_gravito = 0._CUSTOM_REAL
  potential_dot_gravito = 0._CUSTOM_REAL
  potential_dot_dot_gravito = 0._CUSTOM_REAL

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '  done initialization'
    call flush_IMAIN()
  endif

  ! synchronizes all processes
  call synchronize_all()

  end subroutine prepare_timerun_wavefields


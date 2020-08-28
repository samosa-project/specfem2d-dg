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


  subroutine prepare_timerun_PML()

  use specfem_par
  use specfem_par_lns ! TODO: select variables to use.

  implicit none

  ! local parameters
  integer :: i,ier
  character(len=MAX_STRING_LEN) :: outputname

  ! safety check
  if (GPU_MODE .and. PML_BOUNDARY_CONDITIONS ) stop 'error : PML not implemented on GPU mode. Please use Stacey instead'

  ! PML absorbing conditions
  ! sets global flag for all slices
  call any_all_l(anyabs, anyabs_glob)

  if (PML_BOUNDARY_CONDITIONS .and. anyabs_glob) then
    ! user output
    if (myrank == 0) then
      write(IMAIN,*) 'Preparing PML'
      call flush_IMAIN()
    endif

    ! allocates PML arrays
    allocate(spec_to_PML(nspec),stat=ier)
    if (ier /= 0) stop 'error: not enough memory to allocate array spec_to_PML'

    allocate(which_PML_elem(4,nspec),stat=ier)
    if (ier /= 0) stop 'error: not enough memory to allocate array which_PML_elem'
    which_PML_elem(:,:) = .false.

    if (SIMULATION_TYPE == 3 .or. (SIMULATION_TYPE == 1 .and. SAVE_FORWARD)) then
      allocate(PML_interior_interface(4,nspec),stat=ier)
      if (ier /= 0) stop 'error: not enough memory to allocate array PML_interior_interface'
      PML_interior_interface = .false.
    else
      allocate(PML_interior_interface(4,1))
    endif

    ! add support for using PML in MPI mode with external mesh
    if (read_external_mesh) then
      allocate(mask_ibool_pml(nglob))
    else
      allocate(mask_ibool_pml(1))
    endif

    call pml_init()

    if ((SIMULATION_TYPE == 3 .or. (SIMULATION_TYPE == 1 .and. SAVE_FORWARD)) .and. PML_BOUNDARY_CONDITIONS) then

      if (nglob_interface > 0) then
        allocate(point_interface(nglob_interface),stat=ier)
        if (ier /= 0) stop 'error: not enough memory to allocate array point_interface'
      endif

      if (any_elastic .and. nglob_interface > 0) then
        allocate(pml_interface_history_displ(NDIM,nglob_interface,NSTEP),stat=ier)
        if (ier /= 0) stop 'error: not enough memory to allocate array pml_interface_history_displ'
        allocate(pml_interface_history_veloc(NDIM,nglob_interface,NSTEP),stat=ier)
        if (ier /= 0) stop 'error: not enough memory to allocate array pml_interface_history_veloc'
        allocate(pml_interface_history_accel(NDIM,nglob_interface,NSTEP),stat=ier)
        if (ier /= 0) stop 'error: not enough memory to allocate array pml_interface_history_accel'
      endif

      if (any_acoustic .and. nglob_interface > 0) then
        allocate(pml_interface_history_potential(nglob_interface,NSTEP),stat=ier)
        if (ier /= 0) stop 'error: not enough memory to allocate array pml_interface_history_potential'
        allocate(pml_interface_history_potential_dot(nglob_interface,NSTEP),stat=ier)
        if (ier /= 0) stop 'error: not enough memory to allocate array pml_interface_history_potential_dot'
        allocate(pml_interface_history_potential_dot_dot(nglob_interface,NSTEP),stat=ier)
        if (ier /= 0) stop 'error: not enough memory to allocate array pml_interface_history_potential_dot_dot'
      endif

      if (nglob_interface > 0) then
        call determin_interface_pml_interior()
        deallocate(PML_interior_interface)
        deallocate(mask_ibool_pml)
      endif

      if (any_elastic .and. nglob_interface > 0) then
        write(outputname,'(a,i6.6,a)') 'pml_interface_elastic',myrank,'.bin'
        open(unit=71,file='OUTPUT_FILES/'//outputname,status='unknown',form='unformatted')
      endif

      if (any_acoustic .and. nglob_interface > 0) then
        write(outputname,'(a,i6.6,a)') 'pml_interface_acoustic',myrank,'.bin'
        open(unit=72,file='OUTPUT_FILES/'//outputname,status='unknown',form='unformatted')
      endif
    else
      allocate(point_interface(1))
      allocate(pml_interface_history_displ(NDIM,1,1))
      allocate(pml_interface_history_veloc(NDIM,1,1))
      allocate(pml_interface_history_accel(NDIM,1,1))
      allocate(pml_interface_history_potential(1,1))
      allocate(pml_interface_history_potential_dot(1,1))
      allocate(pml_interface_history_potential_dot_dot(1,1))
    endif

    if (SIMULATION_TYPE == 3 .and. PML_BOUNDARY_CONDITIONS) then

      if (any_elastic .and. nglob_interface > 0) then
        do it = 1,NSTEP
          do i = 1, nglob_interface
            read(71) pml_interface_history_accel(1,i,it),pml_interface_history_accel(2,i,it), &
                     pml_interface_history_veloc(1,i,it),pml_interface_history_veloc(2,i,it), &
                     pml_interface_history_displ(1,i,it),pml_interface_history_displ(2,i,it)
          enddo
        enddo
        close(71)
      endif

      if (any_acoustic .and. nglob_interface > 0) then
        do it = 1,NSTEP
          do i = 1, nglob_interface
            read(72) pml_interface_history_potential_dot_dot(i,it), &
                     pml_interface_history_potential_dot(i,it), &
                     pml_interface_history_potential(i,it)
          enddo
        enddo
        close(72)
      endif
    endif

    deallocate(which_PML_elem)

    if (nspec_PML==0) nspec_PML=1 ! DK DK added this

    if (nspec_PML > 0) then
      allocate(K_x_store(NGLLX,NGLLZ,nspec_PML),stat=ier)
      if (ier /= 0) stop 'error: not enough memory to allocate array K_x_store'
      allocate(K_z_store(NGLLX,NGLLZ,nspec_PML),stat=ier)
      if (ier /= 0) stop 'error: not enough memory to allocate array K_z_store'
      allocate(d_x_store(NGLLX,NGLLZ,nspec_PML),stat=ier)
      if (ier /= 0) stop 'error: not enough memory to allocate array d_x_store'
      allocate(d_z_store(NGLLX,NGLLZ,nspec_PML),stat=ier)
      if (ier /= 0) stop 'error: not enough memory to allocate array d_z_store'
      allocate(alpha_x_store(NGLLX,NGLLZ,nspec_PML),stat=ier)
      if (ier /= 0) stop 'error: not enough memory to allocate array alpha_x_store'
      allocate(alpha_z_store(NGLLX,NGLLZ,nspec_PML),stat=ier)
      if (ier /= 0) stop 'error: not enough memory to allocate array alpha_z_store'
      K_x_store(:,:,:) = ZERO
      K_z_store(:,:,:) = ZERO
      d_x_store(:,:,:) = ZERO
      d_z_store(:,:,:) = ZERO
      alpha_x_store(:,:,:) = ZERO
      alpha_z_store(:,:,:) = ZERO
      call define_PML_coefficients()
    else
      allocate(K_x_store(NGLLX,NGLLZ,1),stat=ier)
      if (ier /= 0) stop 'error: not enough memory to allocate array K_x_store'
      allocate(K_z_store(NGLLX,NGLLZ,1),stat=ier)
      if (ier /= 0) stop 'error: not enough memory to allocate array K_z_store'
      allocate(d_x_store(NGLLX,NGLLZ,1),stat=ier)
      if (ier /= 0) stop 'error: not enough memory to allocate array d_x_store'
      allocate(d_z_store(NGLLX,NGLLZ,1),stat=ier)
      if (ier /= 0) stop 'error: not enough memory to allocate array d_z_store'
      allocate(alpha_x_store(NGLLX,NGLLZ,1),stat=ier)
      if (ier /= 0) stop 'error: not enough memory to allocate array alpha_x_store'
      allocate(alpha_z_store(NGLLX,NGLLZ,1),stat=ier)
      if (ier /= 0) stop 'error: not enough memory to allocate array alpha_z_store'
      K_x_store(:,:,:) = ZERO
      K_z_store(:,:,:) = ZERO
      d_x_store(:,:,:) = ZERO
      d_z_store(:,:,:) = ZERO
      alpha_x_store(:,:,:) = ZERO
      alpha_z_store(:,:,:) = ZERO
    endif

    ! elastic PML memory variables
    if (any_elastic .and. nspec_PML > 0) then
      allocate(rmemory_displ_elastic(2,NDIM,NGLLX,NGLLZ,nspec_PML),stat=ier)
      if (ier /= 0) stop 'error: not enough memory to allocate array rmemory_displ_elastic'
      allocate(rmemory_dux_dx(NGLLX,NGLLZ,nspec_PML,2),stat=ier)
      if (ier /= 0) stop 'error: not enough memory to allocate array rmemory_dux_dx'
      allocate(rmemory_dux_dz(NGLLX,NGLLZ,nspec_PML,2),stat=ier)
      if (ier /= 0) stop 'error: not enough memory to allocate array rmemory_dux_dz'
      allocate(rmemory_duz_dx(NGLLX,NGLLZ,nspec_PML,2),stat=ier)
      if (ier /= 0) stop 'error: not enough memory to allocate array rmemory_duz_dx'
      allocate(rmemory_duz_dz(NGLLX,NGLLZ,nspec_PML,2),stat=ier)
      if (ier /= 0) stop 'error: not enough memory to allocate array rmemory_duz_dz'
      if (any_acoustic .and. num_fluid_solid_edges > 0) then
        allocate(rmemory_fsb_displ_elastic(1,NDIM,NGLLX,NGLLZ,num_fluid_solid_edges),stat=ier)
        if (ier /= 0) stop 'error: not enough memory to allocate array rmemory_fsb_displ_elastic'
        allocate(rmemory_sfb_potential_ddot_acoustic(1,NGLLX,NGLLZ,num_fluid_solid_edges),stat=ier)
        if (ier /= 0) stop 'error: not enough memory to allocate array rmemory_sfb_potential_ddot_acoustic'
      endif

      if (ROTATE_PML_ACTIVATE) then
        allocate(rmemory_dux_dx_prime(NGLLX,NGLLZ,nspec_PML,2),stat=ier)
        if (ier /= 0) stop 'error: not enough memory to allocate array rmemory_dux_dx_prime'
        allocate(rmemory_dux_dz_prime(NGLLX,NGLLZ,nspec_PML,2),stat=ier)
        if (ier /= 0) stop 'error: not enough memory to allocate array rmemory_dux_dz_prime'
        allocate(rmemory_duz_dx_prime(NGLLX,NGLLZ,nspec_PML,2),stat=ier)
        if (ier /= 0) stop 'error: not enough memory to allocate array rmemory_duz_dx_prime'
        allocate(rmemory_duz_dz_prime(NGLLX,NGLLZ,nspec_PML,2),stat=ier)
        if (ier /= 0) stop 'error: not enough memory to allocate array rmemory_duz_dz_prime'
      else
        allocate(rmemory_dux_dx_prime(1,1,1,2),stat=ier)
        if (ier /= 0) stop 'error: not enough memory to allocate array rmemory_dux_dx_prime'
        allocate(rmemory_dux_dz_prime(1,1,1,2),stat=ier)
        if (ier /= 0) stop 'error: not enough memory to allocate array rmemory_dux_dz_prime'
        allocate(rmemory_duz_dx_prime(1,1,1,2),stat=ier)
        if (ier /= 0) stop 'error: not enough memory to allocate array rmemory_duz_dx_prime'
        allocate(rmemory_duz_dz_prime(1,1,1,2),stat=ier)
        if (ier /= 0) stop 'error: not enough memory to allocate array rmemory_duz_dz_prime'
      endif

      if (time_stepping_scheme == 2) then
        allocate(rmemory_displ_elastic_LDDRK(2,NDIM,NGLLX,NGLLZ,nspec_PML),stat=ier)
        if (ier /= 0) stop 'error: not enough memory to allocate array rmemory_displ_elastic'
        allocate(rmemory_dux_dx_LDDRK(NGLLX,NGLLZ,nspec_PML,2),stat=ier)
        if (ier /= 0) stop 'error: not enough memory to allocate array rmemory_dux_dx'
        allocate(rmemory_dux_dz_LDDRK(NGLLX,NGLLZ,nspec_PML,2),stat=ier)
        if (ier /= 0) stop 'error: not enough memory to allocate array rmemory_dux_dz'
        allocate(rmemory_duz_dx_LDDRK(NGLLX,NGLLZ,nspec_PML,2),stat=ier)
        if (ier /= 0) stop 'error: not enough memory to allocate array rmemory_duz_dx'
        allocate(rmemory_duz_dz_LDDRK(NGLLX,NGLLZ,nspec_PML,2),stat=ier)
        if (ier /= 0) stop 'error: not enough memory to allocate array rmemory_duz_dz'
        if (any_acoustic .and. num_fluid_solid_edges > 0) then
          allocate(rmemory_fsb_displ_elastic_LDDRK(1,NDIM,NGLLX,NGLLZ,num_fluid_solid_edges),stat=ier)
          if (ier /= 0) stop 'error: not enough memory to allocate array rmemory_fsb_displ_elastic'
          allocate(rmemory_sfb_potential_ddot_acoustic_LDDRK(1,NGLLX,NGLLZ,num_fluid_solid_edges),stat=ier)
          if (ier /= 0) stop 'error: not enough memory to allocate array rmemory_sfb_potential_ddot_acoustic'
        endif
      else
        allocate(rmemory_displ_elastic_LDDRK(1,1,1,1,1),stat=ier)
        allocate(rmemory_dux_dx_LDDRK(1,1,1,2),stat=ier)
        allocate(rmemory_dux_dz_LDDRK(1,1,1,2),stat=ier)
        allocate(rmemory_duz_dx_LDDRK(1,1,1,2),stat=ier)
        allocate(rmemory_duz_dz_LDDRK(1,1,1,2),stat=ier)
        if (any_acoustic .and. num_fluid_solid_edges > 0) then
          allocate(rmemory_fsb_displ_elastic_LDDRK(1,NDIM,NGLLX,NGLLZ,1),stat=ier)
          if (ier /= 0) stop 'error: not enough memory to allocate array rmemory_fsb_displ_elastic'
          allocate(rmemory_sfb_potential_ddot_acoustic_LDDRK(1,NGLLX,NGLLZ,1),stat=ier)
          if (ier /= 0) stop 'error: not enough memory to allocate array rmemory_sfb_potential_ddot_acoustic'
        endif
      endif

      rmemory_displ_elastic(:,:,:,:,:) = ZERO
      rmemory_dux_dx(:,:,:,:) = ZERO
      rmemory_dux_dz(:,:,:,:) = ZERO
      rmemory_duz_dx(:,:,:,:) = ZERO
      rmemory_duz_dz(:,:,:,:) = ZERO

      if (any_acoustic .and. num_fluid_solid_edges > 0) then
        rmemory_fsb_displ_elastic(:,:,:,:,:) = ZERO
        rmemory_sfb_potential_ddot_acoustic(:,:,:,:) = ZERO
      endif

      if (ROTATE_PML_ACTIVATE) then
        rmemory_dux_dx_prime(:,:,:,:) = ZERO
        rmemory_dux_dz_prime(:,:,:,:) = ZERO
        rmemory_duz_dx_prime(:,:,:,:) = ZERO
        rmemory_duz_dz_prime(:,:,:,:) = ZERO
      endif

      if (time_stepping_scheme == 2) then
        rmemory_displ_elastic_LDDRK(:,:,:,:,:) = ZERO
        rmemory_dux_dx_LDDRK(:,:,:,:) = ZERO
        rmemory_dux_dz_LDDRK(:,:,:,:) = ZERO
        rmemory_duz_dx_LDDRK(:,:,:,:) = ZERO
        rmemory_duz_dz_LDDRK(:,:,:,:) = ZERO
        if (any_acoustic .and. num_fluid_solid_edges > 0) then
          rmemory_fsb_displ_elastic_LDDRK(:,:,:,:,:) = ZERO
          rmemory_sfb_potential_ddot_acoustic_LDDRK(:,:,:,:) = ZERO
        endif
      endif

    else

      allocate(rmemory_displ_elastic(1,1,1,1,1))
      allocate(rmemory_dux_dx(1,1,1,1))
      allocate(rmemory_dux_dz(1,1,1,1))
      allocate(rmemory_duz_dx(1,1,1,1))
      allocate(rmemory_duz_dz(1,1,1,1))
      if (any_acoustic .and. num_fluid_solid_edges > 0) then
        allocate(rmemory_fsb_displ_elastic(1,NDIM,NGLLX,NGLLZ,1))
        allocate(rmemory_sfb_potential_ddot_acoustic(1,NGLLX,NGLLZ,1))
        allocate(rmemory_fsb_displ_elastic_LDDRK(1,NDIM,NGLLX,NGLLZ,1))
        allocate(rmemory_sfb_potential_ddot_acoustic_LDDRK(1,NGLLX,NGLLZ,1))
      endif

      allocate(rmemory_dux_dx_prime(1,1,1,1))
      allocate(rmemory_dux_dz_prime(1,1,1,1))
      allocate(rmemory_duz_dx_prime(1,1,1,1))
      allocate(rmemory_duz_dz_prime(1,1,1,1))

      allocate(rmemory_displ_elastic_LDDRK(1,1,1,1,1))
      allocate(rmemory_dux_dx_LDDRK(1,1,1,1))
      allocate(rmemory_dux_dz_LDDRK(1,1,1,1))
      allocate(rmemory_duz_dx_LDDRK(1,1,1,1))
      allocate(rmemory_duz_dz_LDDRK(1,1,1,1))
    endif

    if (any_acoustic .and. nspec_PML > 0) then
      allocate(rmemory_potential_acoustic(2,NGLLX,NGLLZ,nspec_PML),stat=ier)
      if (ier /= 0) stop 'error: not enough memory to allocate array rmemory_potential_acoustic'
      allocate(rmemory_acoustic_dux_dx(NGLLX,NGLLZ,nspec_PML,2),stat=ier)
      if (ier /= 0) stop 'error: not enough memory to allocate array rmemory_acoustic_dux_dx'
      allocate(rmemory_acoustic_dux_dz(NGLLX,NGLLZ,nspec_PML,2),stat=ier)
      if (ier /= 0) stop 'error: not enough memory to allocate array rmemory_acoustic_dux_dz'

      rmemory_potential_acoustic = ZERO
      rmemory_acoustic_dux_dx = ZERO
      rmemory_acoustic_dux_dz = ZERO

      if (time_stepping_scheme == 2) then
        allocate(rmemory_potential_acoustic_LDDRK(2,NGLLX,NGLLZ,nspec_PML),stat=ier)
        if (ier /= 0) stop 'error: not enough memory to allocate array rmemory_potential_acoustic'
        allocate(rmemory_acoustic_dux_dx_LDDRK(NGLLX,NGLLZ,nspec_PML,2),stat=ier)
        if (ier /= 0) stop 'error: not enough memory to allocate array rmemory_acoustic_dux_dx'
        allocate(rmemory_acoustic_dux_dz_LDDRK(NGLLX,NGLLZ,nspec_PML,2),stat=ier)
        if (ier /= 0) stop 'error: not enough memory to allocate array rmemory_acoustic_dux_dz'
      else
        allocate(rmemory_potential_acoustic_LDDRK(1,1,1,1),stat=ier)
        allocate(rmemory_acoustic_dux_dx_LDDRK(1,1,1,1),stat=ier)
        allocate(rmemory_acoustic_dux_dz_LDDRK(1,1,1,1),stat=ier)
      endif

      rmemory_potential_acoustic_LDDRK = ZERO
      rmemory_acoustic_dux_dx_LDDRK = ZERO
      rmemory_acoustic_dux_dz_LDDRK = ZERO
      
! ************************* !
! TODO: MARKED FOR DELETION !
! ************************* !
      !write(*,*) "we are acoustic and we need pml." ! DEBUG
      if(USE_DISCONTINUOUS_METHOD .and. USE_LNS) then
        !write(*,*) "we are acoustic LNS and we need pml." ! DEBUG
        
        ! See [Xie et al., 2014] Xie, Z., Komatitsch, D., Martin, R., and Matzen, R. (2014). Improved forward wave propagation and adjoint-based sensitivity kernel calculations using a numerically stable finite-element PML. Geophysical Journal International, 198(3):1714–1747.
        
        ! With the formulation of [Xie et al., 2014], we'll need 2 auxiliary ADEs, for each constitutive variable. It is conjectured we'll need 3 for 3D. Thus, allocate the first dimension to NDIM.
        ! If changes are needed for this point, one will also need to change the various loops using those arrays.
        
        !allocate(LNS_PML_drho(NDIM,NGLLX,NGLLZ,nspec_PML)) ! Auxiliary evolution variable for constitutive variable 1 (mass conservation).
        !allocate(LNS_PML_rho0dv(NDIM,NDIM,NGLLX,NGLLZ,nspec_PML)) ! Auxiliary evolution variable for constitutive variable 1 (momenta).
        !allocate(LNS_PML_dE(NDIM,NGLLX,NGLLZ,nspec_PML)) ! Auxiliary evolution variable for constitutive variable 1 (energy).
        ! 
        !allocate(RHS_PML_drho(NDIM,NGLLX,NGLLZ,nspec_PML), &
        !         aux_PML_drho(NDIM,NGLLX,NGLLZ,nspec_PML), &
        !         RHS_PML_dE(NDIM,NGLLX,NGLLZ,nspec_PML), &
        !         aux_PML_dE(NDIM,NGLLX,NGLLZ,nspec_PML))
        !allocate(RHS_PML_rho0dv(NDIM,NDIM,NGLLX,NGLLZ,nspec_PML), &
        !         aux_PML_rho0dv(NDIM,NDIM,NGLLX,NGLLZ,nspec_PML))
        ! 
        ! nglob_PML was computed in pml_init(), routine which was called above
        allocate(LNS_PML(2+2*NDIM, LNS_PML_NAUX, nglob_PML), stat=ier) ! Auxiliary evolution variable for all constitutive variables (dimension 1), all auxiliary variables (dimension 2), and all PML mesh points (dimension 3).
        if (ier /= 0) stop 'Error: could not allocate array LNS_PML (see prepare_timerun_pml.f90).'
        !allocate(LNS_PML_drho(NDIM,NGLLX*NGLLZ*nspec_PML)) ! Auxiliary evolution variable for constitutive variable 1 (mass conservation).
        !allocate(LNS_PML_rho0dv(NDIM,NDIM,NGLLX*NGLLZ*nspec_PML)) ! Auxiliary evolution variable for constitutive variable 1 (momenta).
        !allocate(LNS_PML_dE(NDIM,NGLLX*NGLLZ*nspec_PML)) ! Auxiliary evolution variable for constitutive variable 1 (energy).
        
        ! Same as LNS_PML, but arrays for RK iteration.
        allocate(LNS_PML_RHS(2+2*NDIM, LNS_PML_NAUX, nglob_PML), &
                 LNS_PML_aux(2+2*NDIM, LNS_PML_NAUX, nglob_PML), stat=ier)
        if (ier /= 0) stop 'Error: could not allocate arrays LNS_PML_RHS or LNS_PML_aux (see prepare_timerun_pml.f90).'
        !allocate(RHS_PML_drho(NDIM,NGLLX*NGLLZ*nspec_PML), &
        !         aux_PML_drho(NDIM,NGLLX*NGLLZ*nspec_PML), &
        !         RHS_PML_rho0dv(NDIM,NDIM,NGLLX*NGLLZ*nspec_PML), &
        !         aux_PML_rho0dv(NDIM,NDIM,NGLLX*NGLLZ*nspec_PML), &
        !         RHS_PML_dE(NDIM,NGLLX*NGLLZ*nspec_PML), &
        !         aux_PML_dE(NDIM,NGLLX*NGLLZ*nspec_PML), stat=ier)
        !if (ier /= 0) stop 'Error: could not allocate arrays RHS_PML_* or aux_PML_* (see prepare_timerun_pml.f90).'
        
        allocate(LNS_PML_kapp(NDIM,NGLLX,NGLLZ,nspec_PML), &
                 LNS_PML_alpha(NDIM,NGLLX,NGLLZ,nspec_PML), &
                 LNS_PML_a0(NGLLX,NGLLZ,nspec_PML), &
                 LNS_PML_b(NDIM,NGLLX,NGLLZ,nspec_PML), &
                 LNS_PML_d(NDIM,NGLLX,NGLLZ,nspec_PML),stat=ier) ! With the formulation of [Xie et al., 2014], we'll need 2 auxiliary ADEs, for each constitutive variable. It is conjectured we'll need 3 for 3D. Thus, allocate the first dimension to NDIM.
        if (ier /= 0) stop 'Error: could not allocate arrays LNS_PML_{alpha, a0, b, or d} (see prepare_timerun_pml.f90).'
        
        !allocate(rmass_inverse_acoustic_LNS_PML(nglob_PML),stat=ier)
        !if (ier /= 0) stop 'Error: could not allocate rmass_inverse_acoustic_LNS_PML (see prepare_timerun_pml.f90).'
        
        LNS_PML = ZERO
        LNS_PML_RHS = ZERO
        LNS_PML_aux = ZERO
        !LNS_PML_drho = ZERO
        !LNS_PML_rho0dv = ZERO
        !LNS_PML_dE = ZERO
        !RHS_PML_drho = ZERO
        !aux_PML_drho = ZERO
        !RHS_PML_rho0dv = ZERO
        !aux_PML_rho0dv = ZERO
        !RHS_PML_dE = ZERO
        !aux_PML_dE = ZERO
        
        LNS_PML_kapp  = ZERO
        LNS_PML_alpha = ZERO
        LNS_PML_a0    = ZERO
        LNS_PML_b     = ZERO
        LNS_PML_d     = ZERO
        
        !rmass_inverse_acoustic_LNS_PML = ZERO
        
        ! Savage memory free. Not really optimal, but less invasive.
        deallocate(rmemory_potential_acoustic,rmemory_acoustic_dux_dx,rmemory_acoustic_dux_dz, &
                   rmemory_potential_acoustic_LDDRK,rmemory_acoustic_dux_dx_LDDRK,rmemory_acoustic_dux_dz_LDDRK) 
      endif
! ************************* !
! TODO: MARKED FOR DELETION !
! ************************* !
      
    else
      allocate(rmemory_potential_acoustic(1,1,1,1))
      allocate(rmemory_acoustic_dux_dx(1,1,1,1))
      allocate(rmemory_acoustic_dux_dz(1,1,1,1))
    endif

  else
    allocate(rmemory_dux_dx(1,1,1,1))
    allocate(rmemory_dux_dz(1,1,1,1))
    allocate(rmemory_duz_dx(1,1,1,1))
    allocate(rmemory_duz_dz(1,1,1,1))
    allocate(rmemory_fsb_displ_elastic(1,NDIM,NGLLX,NGLLZ,1))
    allocate(rmemory_sfb_potential_ddot_acoustic(1,NGLLX,NGLLZ,1))
    allocate(rmemory_fsb_displ_elastic_LDDRK(1,NDIM,NGLLX,NGLLZ,1))
    allocate(rmemory_sfb_potential_ddot_acoustic_LDDRK(1,NGLLX,NGLLZ,1))

    allocate(rmemory_dux_dx_prime(1,1,1,1))
    allocate(rmemory_dux_dz_prime(1,1,1,1))
    allocate(rmemory_duz_dx_prime(1,1,1,1))
    allocate(rmemory_duz_dz_prime(1,1,1,1))

    allocate(rmemory_displ_elastic(1,1,1,1,1))

    allocate(rmemory_displ_elastic_LDDRK(1,1,1,1,1))
    allocate(rmemory_dux_dx_LDDRK(1,1,1,1))
    allocate(rmemory_dux_dz_LDDRK(1,1,1,1))
    allocate(rmemory_duz_dx_LDDRK(1,1,1,1))
    allocate(rmemory_duz_dz_LDDRK(1,1,1,1))

    allocate(rmemory_potential_acoustic(1,1,1,1))
    allocate(rmemory_acoustic_dux_dx(1,1,1,1))
    allocate(rmemory_acoustic_dux_dz(1,1,1,1))

    allocate(rmemory_potential_acoustic_LDDRK(1,1,1,1))
    allocate(rmemory_acoustic_dux_dx_LDDRK(1,1,1,1))
    allocate(rmemory_acoustic_dux_dz_LDDRK(1,1,1,1))

    allocate(spec_to_PML(1))

    allocate(K_x_store(1,1,1))
    allocate(K_z_store(1,1,1))
    allocate(d_x_store(1,1,1))
    allocate(d_z_store(1,1,1))
    allocate(alpha_x_store(1,1,1))
    allocate(alpha_z_store(1,1,1))
  endif ! PML_BOUNDARY_CONDITIONS

  ! avoid a potential side effect owing to the "if" statements above: this array may be unallocated,
  ! if so we need to allocate a dummy version in order to be able to use that array as an argument
  ! in some subroutine calls below
  if (.not. allocated(rmemory_fsb_displ_elastic)) allocate(rmemory_fsb_displ_elastic(1,NDIM,NGLLX,NGLLZ,1))
  if (.not. allocated(rmemory_sfb_potential_ddot_acoustic)) allocate(rmemory_sfb_potential_ddot_acoustic(1,NGLLX,NGLLZ,1))
  if (.not. allocated(rmemory_fsb_displ_elastic_LDDRK)) then
    allocate(rmemory_fsb_displ_elastic_LDDRK(1,NDIM,NGLLX,NGLLZ,1))
  endif
  if (.not. allocated(rmemory_sfb_potential_ddot_acoustic_LDDRK)) then
    allocate(rmemory_sfb_potential_ddot_acoustic_LDDRK(1,NGLLX,NGLLZ,1))
  endif
  
  !stop "kek" ! DEBUG
  
  end subroutine prepare_timerun_PML



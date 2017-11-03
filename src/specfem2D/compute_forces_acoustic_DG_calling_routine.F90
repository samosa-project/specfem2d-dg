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

  subroutine compute_forces_acoustic_DG_main()

  use specfem_par

  implicit none

  ! local parameters
  ! for rk44
  !double precision :: weight_rk
  
  real(kind=CUSTOM_REAL) :: timelocal
  
  ! Parameters
  real(kind=CUSTOM_REAL), parameter :: ZEROl = 0._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: ONEl  = 1._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: TWOl  = 2._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: SIXl  = 6._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: HALFl = 0.5_CUSTOM_REAL
  
  real(kind=CUSTOM_REAL), dimension(5) :: rk4a, rk4b, rk4c
  double precision, dimension(5) :: rk4a_d, rk4b_d, rk4c_d
  
  real(kind=CUSTOM_REAL), dimension(nglob_DG) :: veloc_x
  
  !integer :: i, j, ispec, numelem
  
  real(kind=CUSTOM_REAL), parameter :: threshold = 0.0000001_CUSTOM_REAL
  
  character(len=100) file_name
  integer :: i,j,ispec
  
  !integer :: i ,j, ispec
  
  ! Checks if anything has to be done.
  if (.not. any_acoustic_DG) then
    return
  endif
  
  ! Intialisation.
  if(it == 1 .AND. i_stage == 1) then
    ! Slope limiter initialisation.
    ! TODO: put this in an "if(USE_SLOPE_LIMITER)"? In other words, check wether or not the Vandermonde matrix is used when the slope limiter is not.
    call setUpVandermonde()
    
    resu_rhovx   = ZEROl
    resu_rhovz   = ZEROl
    resu_rho     = ZEROl
    resu_E       = ZEROl
    resu_e1      = ZEROl
    dot_rho(:)   = ZEROl
    dot_rhovx(:) = ZEROl
    dot_rhovz(:) = ZEROl
    dot_E(:)     = ZEROl
    dot_e1(:)    = ZEROl
    
    allocate(save_pressure(nglob_DG))
    save_pressure = ZEROl
    
    if(.not. trim(MODEL) == 'external') then
      deallocate(gravityext, muext, etaext, kappa_DG, &
                 tau_epsilon, tau_sigma)
      allocate(gravityext(NGLLX, NGLLZ, nspec), &
               etaext(NGLLX, NGLLZ, nspec), &
               muext(NGLLX, NGLLZ, nspec), &
               kappa_DG(NGLLX, NGLLZ, nspec), &
               tau_epsilon(NGLLX, NGLLZ, nspec), &
               tau_sigma(NGLLX, NGLLZ, nspec)) 
    endif
    
    !rhovx_DG   = ZEROl
    !rhovz_DG   = ZEROl
    !rho_DG     = ONEl
    !gravityext = ZEROl
    !E_DG = ZEROl
    !kappa_DG = ZEROl
    !etaext = ZEROl
    !muext = ZEROl
    allocate(rhovx_init(nglob_DG), &
        rhovz_init(nglob_DG), &
        E_init(nglob_DG), &
        rho_init(nglob_DG))
    
    call prepare_MPI_DG()

    allocate(ispec_is_acoustic_coupling_ac(nglob_DG))
    ispec_is_acoustic_coupling_ac = -1
    if(.not. only_DG_acoustic) call find_DG_acoustic_coupling()
    
    !stop
    T_DG = 0
    V_DG = 0

    call initial_condition_DG()

    !if(timelocal == 0) then
    rhovx_init = rhovx_DG
    rhovz_init = rhovz_DG
    E_init     = E_DG
    rho_init   = rho_DG
    p_DG_init  = (gammaext_DG - ONE)*( E_DG &
        - (HALF/rho_DG)*( rhovx_DG**2 + rhovz_DG**2 ) )
    T_init = (E_DG/rho_DG - 0.5*((rhovx_DG/rho_DG)**2 + (rhovz_DG/rho_DG)**2))/(cnu)
    !endif
    
#ifdef USE_MPI
    if (NPROC > 1 .and. ninterface_acoustic > 0) then
      call assemble_MPI_vector_DG(gammaext_DG, buffer_DG_gamma_P)
    endif
#endif
    
  endif !if(it == 1 .AND. i_stage == 1)
  
  if(it == 1 .AND. i_stage == 1 .AND. .false.) then
    write(file_name, "('./boundaries_MPI_',i3.3)") myrank
    open(10, file=file_name, form='formatted')
    do ispec = 1,nspec
      ! acoustic spectral element
      !if (ispec_is_acoustic(ispec)) then
      !if (ispec_is_acoustic_DG(ispec)) then
      !if(ispec_is_acoustic_DG(ispec)) then
      ! first double loop over GLL points to compute and store gradients
      do j = 1,5
        do i = 1,5
          WRITE(10, *) coord(:, ibool_before_perio(i, j, ispec))
        enddo
      enddo
      !endif
    enddo
    close(10)
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
  
  ! RK3
  !rk4a_d(1) = 0d0
  !rk4a_d(2) = -567301805773.0/1357537059087.0
  !rk4a_d(3) = -2404267990393.0/2016746695238.0
  !rk4a_d(4) = -3550918686646.0/2091501179385.0
  !rk4a_d(5) = -1275806237668.0/842570457699.0
    
  !rk4b_d(1) = 1432997174477.0/9575080441755.0 
  !rk4b_d(2) = 5161836677717.0/13612068292357.0 
  !rk4b_d(3) = 1720146321549.0/2090206949498.0 
  !rk4b_d(4) = 3134564353537.0/4481467310338.0 
  !rk4b_d(5) = 2277821191437.0/14882151754819.0
    
  !rk4c_d(1) = 0d0
  !rk4c_d(2) = 1d0
  !rk4c_d(3) = 1d0/2d0
  !rk4c_d(4) = 0d0
  !rk4c_d(5) = 0d0
  
  rk4a = real(rk4a_d, kind=CUSTOM_REAL)
  rk4b = real(rk4b_d, kind=CUSTOM_REAL)
  rk4c = real(rk4c_d, kind=CUSTOM_REAL)
  !rk4a = rk4a_d
  !rk4b = rk4b_d
  !rk4c = rk4c_d
  
  timelocal = (it-1)*deltat + rk4c(i_stage)*deltat
  
  ! TODO: introduce a verbosity parameter in order to prevent unwanted flooding of the terminal.
  if(myrank == 0 .AND. mod(it, 50)==0) then
    WRITE(*,*) "****************************************************************"
    WRITE(*,"(a,i5,a,i1,a,e23.16,a,i3,a)") "Iteration", it, ", stage ", i_stage, ", local time", timelocal, &
    ". Informations for process number ", myrank, "."
  endif
  
#ifdef USE_MPI
  if (NPROC > 1 .and. ninterface_acoustic > 0) then
    call assemble_MPI_vector_DG(rho_DG, buffer_DG_rho_P)
    call assemble_MPI_vector_DG(rhovx_DG, buffer_DG_rhovx_P)
    call assemble_MPI_vector_DG(rhovz_DG, buffer_DG_rhovz_P)
    call assemble_MPI_vector_DG(E_DG, buffer_DG_E_P)
  endif
#endif

  ! Local Discontinuous Galerkin for viscous fluxes.
  if((maxval(muext) > 0 .OR. maxval(etaext) > 0 .OR. maxval(kappa_DG) > 0)) then
#ifdef USE_MPI
    if (NPROC > 1 .and. ninterface_acoustic > 0) then
      call assemble_MPI_vector_DG(T_DG(1,:), buffer_DG_Tx_P)
      call assemble_MPI_vector_DG(T_DG(2,:), buffer_DG_Tz_P)
      call assemble_MPI_vector_DG(V_DG(1,1,:), buffer_DG_Vxx_P)
      call assemble_MPI_vector_DG(V_DG(2,2,:), buffer_DG_Vzz_P)
      call assemble_MPI_vector_DG(V_DG(1,2,:), buffer_DG_Vxz_P)
      call assemble_MPI_vector_DG(V_DG(2,1,:), buffer_DG_Vzx_P)
    endif
#endif
    call compute_viscous_tensors(T_DG, V_DG, rho_DG, rhovx_DG, rhovz_DG, E_DG, timelocal)
  endif
  
  call compute_forces_acoustic_DG(rho_DG, rhovx_DG, rhovz_DG, E_DG, &
        T_DG, V_DG, e1_DG, &
        dot_rho, dot_rhovx, dot_rhovz, dot_E, dot_e1, &
        timelocal)
  
  if (time_stepping_scheme == 3) then
    
    save_pressure = (gammaext_DG - ONEl)*( E_DG &
        - (HALFl)*(ONEl/rho_DG)*( rhovx_DG**2 + rhovz_DG**2 ) )
    
    ! Inverse mass matrix
    dot_rho(:)   = dot_rho(:)   * rmass_inverse_acoustic_DG(:)
    dot_rhovx(:) = dot_rhovx(:) * rmass_inverse_acoustic_DG(:)
    dot_rhovz(:) = dot_rhovz(:) * rmass_inverse_acoustic_DG(:)
    dot_E(:)     = dot_E(:)     * rmass_inverse_acoustic_DG(:)
    !dot_e1(:)     = dot_e1(:) * rmass_inverse_acoustic_DG(:)
  
    ! RK5-low dissipation Update
    resu_rho = rk4a(i_stage)*resu_rho + deltat*dot_rho
    resu_rhovx = rk4a(i_stage)*resu_rhovx + deltat*dot_rhovx
    resu_rhovz = rk4a(i_stage)*resu_rhovz + deltat*dot_rhovz
    resu_E = rk4a(i_stage)*resu_E + deltat*dot_E
    resu_e1 = rk4a(i_stage)*resu_e1 + deltat*dot_e1
    
    rho_DG   = rho_DG + rk4b(i_stage)*resu_rho
    rhovx_DG = rhovx_DG + rk4b(i_stage)*resu_rhovx
    rhovz_DG = rhovz_DG + rk4b(i_stage)*resu_rhovz
    E_DG     = E_DG + rk4b(i_stage)*resu_E 
    e1_DG    = e1_DG + rk4b(i_stage)*resu_e1 
    
    ! If we want to compute kernels we save regularly.
    if(.false.) then
      call save_forward_solution()
    endif
    
    ! Check non-positivity.
    if(minval(rho_DG) < 10d-14) then
      WRITE(*,*) "*********************************************"
      WRITE(*,*) "*********************************************"
      WRITE(*,*) "*********************************************"
      WRITE(*,*) "CAREFUL, VERY SMALL DENSITY: ", minval(rho_DG)
      WRITE(*,*) "*********************************************"
      WRITE(*,*) "*********************************************"
      WRITE(*,*) "*********************************************"
    endif
  endif
  
  ! --------------------------- !
  ! Remove high-order           !
  ! coefficients of the         !
  ! solution.                   !
  ! --------------------------- !
  if(USE_SLOPE_LIMITER) then
    ! rho.
    if(CONSTRAIN_HYDROSTATIC) then
      veloc_x = rho_DG - rho_init
      call SlopeLimit1(veloc_x, timelocal, 1)
      rho_DG = veloc_x + rho_init
    else
      call SlopeLimit1(rho_DG, timelocal, 1)
    endif
    ! rhovx.
    if(CONSTRAIN_HYDROSTATIC) then
      veloc_x = rhovx_DG - rhovx_init
      call SlopeLimit1(veloc_x, timelocal, 1)
      rhovx_DG = veloc_x + rhovx_init
    else
      call SlopeLimit1(rhovx_DG, timelocal, 2)
    endif
    ! rhovz.
    if(CONSTRAIN_HYDROSTATIC) then
      veloc_x = rhovz_DG - rhovz_init
      call SlopeLimit1(veloc_x, timelocal, 1)
      rhovz_DG = veloc_x + rhovz_init
    else
      call SlopeLimit1(rhovz_DG, timelocal, 3)
    endif
    ! E.
    if(CONSTRAIN_HYDROSTATIC) then
      veloc_x = E_DG - E_init
      call SlopeLimit1(veloc_x, timelocal, 1)
      E_DG = veloc_x + E_init
    else
      call SlopeLimit1(E_DG, timelocal, 4)
    endif
  endif
  
  ! --------------------------- !
  ! Damp solution in a buffer   !
  ! zone.                       !
  ! --------------------------- !
  if(.true.) then ! TODO: add a parameter for this option in parfile.
    call absorb_condition_DG(rho_DG, rhovx_DG, rhovz_DG, E_DG, timelocal) ! See "boundary_terms_DG.f90".
  endif
  
  end subroutine compute_forces_acoustic_DG_main

!
!-------------------------------------------------------------------------------------
!

  subroutine compute_forces_acoustic_main_backward_DG()

  use specfem_par

  implicit none

  ! local parameters
  ! for rk44
  !double precision :: weight_rk
  
  real(kind=CUSTOM_REAL) :: timelocal
  
  ! Parameters
  real(kind=CUSTOM_REAL), parameter :: ZEROl = 0._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: ONEl  = 1._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: TWOl  = 2._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: SIXl  = 6._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: HALFl = 0.5_CUSTOM_REAL
  
  real(kind=CUSTOM_REAL), dimension(5) :: rk4a, rk4b, rk4c
  double precision, dimension(5) :: rk4a_d, rk4b_d, rk4c_d
  
  !integer :: i, j, ispec, numelem
  
  real(kind=CUSTOM_REAL), parameter :: threshold = 0.0000001_CUSTOM_REAL
  !integer :: i ,j, ispec
  ! checks if anything to do in this slice
  if (.not. any_acoustic) return
  ! Intialization
  if(it == 1 .AND. i_stage == 1) then
  
    ! Sloper limiter initialization
    call setUpVandermonde()
  
    resu_rhovx = ZEROl
    resu_rhovz = ZEROl
    resu_rho   = ZEROl
    resu_E     = ZEROl
    dot_rho(:)   = ZEROl
    dot_rhovx(:) = ZEROl
    dot_rhovz(:) = ZEROl
    dot_E(:)    = ZEROl
    
    if(.not. trim(MODEL) == 'external') then
        deallocate(muext, etaext,kappa_DG)
        allocate(etaext(NGLLX, NGLLZ, nspec), &
              muext(NGLLX, NGLLZ, nspec), &
              kappa_DG(NGLLX, NGLLZ, nspec) ) 
        deallocate(gravityext, windxext, &
              rhoext, vpext, pext_DG)
        allocate(gravityext(NGLLX, NGLLZ, nspec), &
              windxext(NGLLX, NGLLZ, nspec), &
              rhoext(NGLLX, NGLLZ, nspec), &
              vpext(NGLLX, NGLLZ, nspec), &
              pext_DG(NGLLX, NGLLZ, nspec) &
              )
    endif
    
    !rhovx_DG   = ZEROl
    !rhovz_DG   = ZEROl
    !rho_DG     = ONEl
    !gravityext = ZEROl
    !E_DG = ZEROl
    !kappa_DG = ZEROl
    !etaext = ZEROl
    !muext = ZEROl
    
    call initial_condition_DG_backward()
    
    call prepare_MPI_DG()
    !call prepare_MPI_DG(my_neighbours)
    
#ifdef USE_MPI
  if (NPROC > 1 .and. ninterface_acoustic > 0) then
    call assemble_MPI_vector_DG(gammaext_DG, buffer_DG_gamma_P)
  endif
#endif
    
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
  
  rk4a = real(rk4a_d, kind=CUSTOM_REAL)
  rk4b = real(rk4b_d, kind=CUSTOM_REAL)
  rk4c = real(rk4c_d, kind=CUSTOM_REAL)
  
  timelocal = (it-1)*deltat + rk4c(i_stage)*deltat
  
  if(myrank == 0) &
  WRITE(*,*) "iter", it, i_stage, timelocal
#ifdef USE_MPI
  if (NPROC > 1 .and. ninterface_acoustic > 0) then
    call assemble_MPI_vector_DG(rho_DG, buffer_DG_rho_P)
    call assemble_MPI_vector_DG(rhovx_DG, buffer_DG_rhovx_P)
    call assemble_MPI_vector_DG(rhovz_DG, buffer_DG_rhovz_P)
    call assemble_MPI_vector_DG(E_DG, buffer_DG_E_P)
  endif
#endif
  
  call compute_forces_acoustic_DG_backward_real(rho_DG, rhovx_DG, rhovz_DG, E_DG, &
        dot_rho, dot_rhovx, dot_rhovz, dot_E, &
        timelocal)
  
  !WRITE(*,*) timelocal,"MINMAX >>", myrank,minval(sqrt(rho_DG)), maxval(sqrt(rho_DG))
        
  if (time_stepping_scheme == 3) then
    
    ! Inverse mass matrix
    dot_rho(:)   = dot_rho(:) * rmass_inverse_acoustic_DG(:)
    dot_rhovx(:) = dot_rhovx(:) * rmass_inverse_acoustic_DG(:)
    dot_rhovz(:) = dot_rhovz(:) * rmass_inverse_acoustic_DG(:)
    dot_E(:)     = dot_E(:) * rmass_inverse_acoustic_DG(:)
  
    ! RK5-low dissipation Update
    resu_rho = rk4a(i_stage)*resu_rho + deltat*dot_rho
    resu_rhovx = rk4a(i_stage)*resu_rhovx + deltat*dot_rhovx
    resu_rhovz = rk4a(i_stage)*resu_rhovz + deltat*dot_rhovz
    resu_E = rk4a(i_stage)*resu_E + deltat*dot_E
    
    rho_DG   = rho_DG + rk4b(i_stage)*resu_rho
    rhovx_DG = rhovx_DG + rk4b(i_stage)*resu_rhovx
    rhovz_DG = rhovz_DG + rk4b(i_stage)*resu_rhovz
    E_DG     = E_DG + rk4b(i_stage)*resu_E 
    
    ! If we want to compute kernels we save regurlaly
    if(.false.) then
    call save_forward_solution()
    endif
    
  endif

  end subroutine compute_forces_acoustic_main_backward_DG

!
!-------------------------------------------------------------------------------------
!

  subroutine compute_forces_acoustic_main_backward_DG_other()

  use specfem_par

  implicit none

  ! local parameters
  integer :: it_temp, istage_temp
  
  ! Parameters
  real(kind=CUSTOM_REAL), parameter :: ZEROl = 0._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: ONEl  = 1._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: TWOl  = 2._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: SIXl  = 6._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: HALFl = 0.5_CUSTOM_REAL

  real(kind=CUSTOM_REAL), dimension(5) :: rk4a, rk4b, rk4c
  double precision, dimension(5) :: rk4a_d, rk4b_d, rk4c_d
  
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
  
  ! RK3
  !rk4a_d(1) = 0d0
  !rk4a_d(2) = -567301805773.0/1357537059087.0
  !rk4a_d(3) = -2404267990393.0/2016746695238.0
  !rk4a_d(4) = -3550918686646.0/2091501179385.0
  !rk4a_d(5) = -1275806237668.0/842570457699.0
    
  !rk4b_d(1) = 1432997174477.0/9575080441755.0 
  !rk4b_d(2) = 5161836677717.0/13612068292357.0 
  !rk4b_d(3) = 1720146321549.0/2090206949498.0 
  !rk4b_d(4) = 3134564353537.0/4481467310338.0 
  !rk4b_d(5) = 2277821191437.0/14882151754819.0
    
  !rk4c_d(1) = 0d0
  !rk4c_d(2) = 1d0
  !rk4c_d(3) = 1d0/2d0
  !rk4c_d(4) = 0d0
  !rk4c_d(5) = 0d0
  
  rk4a = real(rk4a_d, kind=CUSTOM_REAL)
  rk4b = real(rk4b_d, kind=CUSTOM_REAL)
  rk4c = real(rk4c_d, kind=CUSTOM_REAL)

  ! checks
  if (SIMULATION_TYPE /= 3 ) return

  ! checks if anything to do in this slice
  if (.not. any_acoustic) return

  ! timing
  !if (UNDO_ATTENUATION) then
  !  ! time increment
  !  ! example: NSTEP = 800, NT_DUMP_ATTENUATION = 500 -> 1. subset: it_temp = (2-1)*500 + 1 = 501,502,..,800
  !  !                                                 -> 2. subset: it_temp = (2-2)*500 + 1 = 1,2,..,500
  !  it_temp = (NSUBSET_ITERATIONS - iteration_on_subset)*NT_DUMP_ATTENUATION + it_of_this_subset
  !  ! time scheme
  !  istage_temp = i_stage
  !else
    ! time increment
    !it_temp = NSTEP - it + 1
    ! time scheme
    !istage_temp = stage_time_scheme - i_stage + 1
    
    ! time increment
    it_temp = it
    ! time scheme
    istage_temp = i_stage
    
  !WRITE(*,*) "it/istage", it_temp, istage_temp
    
  !endif
  
   if(.false.) then
   
   !if(it_temp == NSTEP .AND. istage_temp == stage_time_scheme) then
  if(it_temp == 1 .AND. istage_temp == 1) then
  
    resu_b_rhovx = 0.
    resu_b_rhovz = 0.
    resu_b_rho   = 0.
    resu_b_E     = 0.
    
    b_rho_DG      = 0.
    b_rhovx_DG    = 0.
    b_rhovz_DG    = 0.
    b_E_DG        = 0.
    
    ! Read forward solution at stations and data (u_obs - u_forward) for adjoint source
    !call read_source_adj_DG()
    
    endif

   b_dot_rho(:)   = 0.
   b_dot_rhovx(:) = 0.
   b_dot_rhovz(:) = 0.
   b_dot_E(:)     = 0.
   
   ! Source function defined by misfit
   call compute_add_sources_acoustic_DG_backward(it_temp, &
        b_dot_rho, b_dot_rhovx, b_dot_rhovz, b_dot_E)
   
   call compute_forces_acoustic_backward_DG(b_rho_DG, b_rhovx_DG, b_rhovz_DG, b_E_DG, &
        rho_DG, rhovx_DG, rhovz_DG, E_DG, b_dot_rho, b_dot_rhovx, b_dot_rhovz, b_dot_E)
   
     ! assembling potential_dot_dot or b_potential_dot_dot for acoustic elements
#ifdef USE_MPI
  if (NPROC > 1 .and. ninterface_acoustic > 0) then
    call assemble_MPI_vector_ac(b_dot_rho)
    call assemble_MPI_vector_ac(b_dot_rhovx)
    call assemble_MPI_vector_ac(b_dot_rhovz)
    call assemble_MPI_vector_ac(b_dot_E)
  endif
#endif

  !WRITE(*,*) "1>>>>>>>>>>>",istage_temp, maxval(rmass_inverse_acoustic_DG_b), minval(rmass_inverse_acoustic_DG_b), &
  !      maxval(b_dot_rho), minval(b_dot_rho), maxval(b_dot_rhovz), minval(b_dot_rhovz)

  ! free surface for an acoustic medium
  !call enforce_acoustic_free_surface(b_dot_rhovx,b_dot_rhovz,b_dot_E)
  !call enforce_acoustic_free_surface(b_dot_rhovx,b_dot_rhovz,b_dot_rho)

  !b_dot_rho(:)   = b_dot_rho(:) * rmass_inverse_acoustic_DG_b(:)
  !b_dot_rhovx(:) = b_dot_rhovx(:) * rmass_inverse_acoustic_DG_b(:)! rmass_inverse_acoustic_DG_b(:)
  !b_dot_rhovz(:) = b_dot_rhovz(:) * rmass_inverse_acoustic_DG_b(:)!rmass_inverse_acoustic_DG_b(:)
  !b_dot_E(:)     = b_dot_E(:) * rmass_inverse_acoustic_DG_b(:)
  
  b_dot_rho(:)   = b_dot_rho(:) * rmass_inverse_acoustic(:)
  b_dot_rhovx(:) = b_dot_rhovx(:) * rmass_inverse_acoustic(:)! rmass_inverse_acoustic_DG_b(:)
  b_dot_rhovz(:) = b_dot_rhovz(:) * rmass_inverse_acoustic(:)!rmass_inverse_acoustic_DG_b(:)
  b_dot_E(:)     = b_dot_E(:) * rmass_inverse_acoustic(:)

  !WRITE(*,*) "2>>>>>>>>>>>",istage_temp, maxval(rmass_inverse_acoustic_DG_b), minval(rmass_inverse_acoustic_DG_b), &
    !    maxval(b_dot_rho), minval(b_dot_rho), maxval(b_dot_rhovz), minval(b_dot_rhovz)
  !if(istage_temp == 2) stop 'KKKK'
  ! multiply by the inverse of the mass matrix and update velocity
  if (time_stepping_scheme == 1) then
    !! DK DK this should be vectorized
    !b_potential_dot_dot_acoustic(:) = b_potential_dot_dot_acoustic(:) * rmass_inverse_acoustic(:)
    b_rho_DG(:)   = b_rho_DG(:)   + b_deltatover2 * b_dot_rho(:)
    b_rhovx_DG(:) = b_rhovx_DG(:) + b_deltatover2 * b_dot_rhovx(:)
    b_rhovz_DG(:) = b_rhovz_DG(:) + b_deltatover2 * b_dot_rhovz(:)
    b_E_DG(:)     = b_E_DG(:)     + b_deltatover2 * b_dot_E(:)
  
  elseif(time_stepping_scheme == 3) then
  
    ! RK5-low dissipation Update
    resu_b_rho   = rk4a(istage_temp)*resu_b_rho   + deltat*b_dot_rho
    resu_b_rhovx = rk4a(istage_temp)*resu_b_rhovx + deltat*b_dot_rhovx
    resu_b_rhovz = rk4a(istage_temp)*resu_b_rhovz + deltat*b_dot_rhovz
    resu_b_E     = rk4a(istage_temp)*resu_b_E     + deltat*b_dot_E
    
    b_rho_DG   = b_rho_DG   + rk4b(istage_temp)*resu_b_rho
    b_rhovx_DG = b_rhovx_DG + rk4b(istage_temp)*resu_b_rhovx
    b_rhovz_DG = b_rhovz_DG + rk4b(istage_temp)*resu_b_rhovz
    b_E_DG     = b_E_DG     + rk4b(istage_temp)*resu_b_E 
  
  endif
   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!
   else
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
   if(it == 1 .AND. i_stage == 1) then
   
    resu_rhovx = ZEROl
    resu_rhovz = ZEROl
    resu_rho   = ZEROl
    resu_E     = ZEROl

    rho_DG = ZEROl
    rhovx_DG = ZEROl
    rhovz_DG   = ZEROl
    E_DG     = ZEROl

    if(.not. trim(MODEL) == 'external') then
        deallocate(gravityext, gammaext_DG, windxext, &
              rhoext, vpext, pext_DG)
        allocate(gravityext(NGLLX, NGLLZ, nspec), &
              gammaext_DG(nglob_DG), &
              windxext(NGLLX, NGLLZ, nspec), &
              rhoext(NGLLX, NGLLZ, nspec), &
              vpext(NGLLX, NGLLZ, nspec), &
              pext_DG(NGLLX, NGLLZ, nspec) &
              )
    endif

    call initial_condition_DG_backward()

    call prepare_MPI_DG()

    endif

    dot_rho(:)   = ZEROl
    dot_rhovx(:) = ZEROl
    dot_rhovz(:) = ZEROl
    dot_E(:)    = ZEROl
   
#ifdef USE_MPI
  if (NPROC > 1 .and. ninterface_acoustic > 0) then
    call assemble_MPI_vector_DG(rho_DG, buffer_DG_rho_P)
    call assemble_MPI_vector_DG(rhovx_DG, buffer_DG_rhovx_P)
    call assemble_MPI_vector_DG(rhovz_DG, buffer_DG_rhovz_P)
    call assemble_MPI_vector_DG(E_DG, buffer_DG_E_P)
  endif
#endif
   
   !call compute_forces_acoustic_DG_backward_real(it_temp, istage_temp, rho_DG, rhovx_DG, rhovz_DG, E_DG, &
   !     dot_rho, dot_rhovx, dot_rhovz, dot_E)
  
  if (time_stepping_scheme == 3) then
    
    ! Inverse mass matrix
    dot_rho(:)   = dot_rho(:) * rmass_inverse_acoustic_DG(:)
    dot_rhovx(:) = dot_rhovx(:) * rmass_inverse_acoustic_DG(:)
    dot_rhovz(:) = dot_rhovz(:) * rmass_inverse_acoustic_DG(:)
    dot_E(:)     = dot_E(:) * rmass_inverse_acoustic_DG(:)
  
    ! RK5-low dissipation Update
    resu_rho = rk4a(i_stage)*resu_rho + deltat*dot_rho
    resu_rhovx = rk4a(i_stage)*resu_rhovx + deltat*dot_rhovx
    resu_rhovz = rk4a(i_stage)*resu_rhovz + deltat*dot_rhovz
    resu_E = rk4a(i_stage)*resu_E + deltat*dot_E
    
    rho_DG   = rho_DG + rk4b(i_stage)*resu_rho
    rhovx_DG = rhovx_DG + rk4b(i_stage)*resu_rhovx
    rhovz_DG = rhovz_DG + rk4b(i_stage)*resu_rhovz
    E_DG     = E_DG + rk4b(i_stage)*resu_E 
   
   
   endif
   
   endif
   
  ! Before computing kernel we read the forward solution saved frames
  !call read_forward_solution(it_temp)

  end subroutine compute_forces_acoustic_main_backward_DG_other

! ------------------------------------------------------------------------------------

  subroutine read_source_adj_DG()
    
    use constants,only: CUSTOM_REAL
    
    use specfem_par, only: NSTEP,myrank, nrec,&!NSOURCES
                         source_time_function_rho_DG, source_time_function_E_DG, &
                         source_time_function_rhovx_DG, source_time_function_rhovz_DG, &
                         deltat, which_proc_receiver,ispec_is_acoustic, ispec_selected_rec!,is_proc_source, ispec_is_acoustic, ispec_selected_source

    use mpi

    implicit none 
    
    include "precision.h"
    
    character(len=100) file_name
    integer :: ok_forward, ok_data, it_tmp, irec, irec_local
    real(kind=CUSTOM_REAL) :: time, time_data, &
        !rho_temp, rhovx_temp, rhovz_temp, E_temp, &
        !rho_temp_data, rhovx_temp_data, rhovz_temp_data, E_temp_data, &
        veloc_z, veloc_z_data
    
   irec_local = 0
   do irec = 1,nrec
   
    ! add the source (only if this proc carries the source)
    if (myrank == which_proc_receiver(irec)) then
    
      irec_local = irec_local + 1
      
      if (ispec_is_acoustic(ispec_selected_rec(irec))) then
    
    !do i_source = 1,NSOURCES
    
    !if(ispec_is_acoustic(ispec_selected_source(i_source))) then
    
    ! if this processor core carries the source and the source element is acoustic
    !if (is_proc_source(i_source) == 1) then
    
        !write(file_name,"('solution_forward_',i3.3,'.txt')") irec
        write(file_name,"('../solution_forward_H15000/AA.S',i3.3,'.BXZ.semd')") irec
        ! Open forward solution source file
        open(10,file=file_name, form='formatted',iostat = ok_forward)
        
        !write(file_name,"('solution_data_',i3.3,'.txt')") irec
        write(file_name,"('../solution_forward_H8417/AA.S',i4.4,'.BXZ.semd')") irec
        ! Open exact solution source file
        open(11,file=file_name, form='formatted',iostat = ok_data)

        it_tmp = -1
        ! COUNT NUMBER OF PARAMS (74)
        do while(ok_forward == 0 .AND. ok_data == 0 .AND. it_tmp <= NSTEP)
          
                read(10,*,iostat = ok_forward)  time, veloc_z
                !        time, rho_temp, rhovx_temp, rhovz_temp, E_temp
                read(11,*,iostat = ok_data) time_data, veloc_z_data
                !        time_data, rho_temp_data, rhovx_temp_data, rhovz_temp_data, E_temp_data
                        
                it_tmp = NSTEP - INT( time/deltat ) + 1
                !it_tmp = NSTEP - it + 1
                
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                !
                !  DESIGN THE MISFIT HERE
                !
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                if(it_tmp <= NSTEP) then
                !adj_sourcearrays_rho_DG(irec_local,it_tmp,1,i,j)   = 0
                !adj_sourcearrays_rhovx_DG(irec_local,it_tmp,1,i,j) = 0
                !adj_sourcearrays_rhovz_DG(irec_local,it_tmp,1,i,j) = veloc_z - veloc_z_data
                !adj_sourcearrays_E_DG(irec_local,it_tmp,1,i,j)     = 0
                source_time_function_rho_DG(irec, it_tmp, :)   = 0!rho_temp - rho_temp_data
                source_time_function_rhovx_DG(irec, it_tmp, :) = 0!rhovx_temp - rhovx_temp_data
                source_time_function_rhovz_DG(irec, it_tmp, :) = veloc_z - veloc_z_data!rhovz_temp - rhovz_temp_data
                source_time_function_E_DG(irec, it_tmp, :)     = 0!E_temp - E_temp_data
                endif
                
        end do 
        
        close(11)
        
        close(10)
    
    endif !ispec_is_acoustic(ispec_selected_rec(irec))
    
    endif !myrank == which_proc_receiver(irec)
    
    enddo !irec = 1,nrec
    
  end subroutine read_source_adj_DG

! ------------------------------------------------------------------------------------

  subroutine save_forward_solution()
    
    use constants,only: CUSTOM_REAL
    
    use specfem_par

    use mpi

    implicit none 
    
    include "precision.h"
    
    character(len=MAX_STRING_LEN) :: outputname
    integer :: ier
    
    integer :: IT_SAVE
    
    ! --------------- TO BE CHANGED ----------------!
    IT_SAVE = 100
    
    !WRITE(*,*) ">>>>>>>>>>>>>>> ", it, IT_SAVE, SIMULATION_TYPE, any_acoustic, SAVE_FORWARD
    
! save last frame
  if (SAVE_FORWARD .and. SIMULATION_TYPE == 1 .and. any_acoustic &
        .and. mod(it, IT_SAVE) == 0 .AND. i_stage == 1) then
  
    if (myrank == 0) then
      write(IMAIN,*)
      write(IMAIN,*) 'Saving DG acoustic frame...'
      write(IMAIN,*)
      call flush_IMAIN()
    endif

    write(outputname,'(a,i6.6,a,i6.6,a)') 'frame_DG_',myrank,'_',it,'.bin'
    open(unit=55,file='OUTPUT_FILES/output_frames/'//trim(outputname),status='unknown',form='unformatted',iostat=ier)
    if (ier /= 0) call exit_MPI(myrank,'Error opening file frame_DG...')

    write(55) rho_DG
    write(55) rhovx_DG
    write(55) rhovz_DG
    write(55) E_DG

    close(55)
    
  endif

  end subroutine save_forward_solution

! ------------------------------------------------------------------------------------

  subroutine read_forward_solution(it_temp)
    
    use constants,only: CUSTOM_REAL
    
    use specfem_par

    use mpi

    implicit none 
    
    include "precision.h"
    
    character(len=MAX_STRING_LEN) :: outputname
    integer :: ier
    
    integer :: it_temp, IT_SAVE
    
    ! --------------- TO BE CHANGED ----------------!
    IT_SAVE = 500
    
! save last frame
  if (SAVE_FORWARD .and. SIMULATION_TYPE == 1 .and. any_acoustic &
        .and. mod(it, IT_SAVE) == 0 .AND. i_stage == 1) then
  
    if (myrank == 0) then
      write(IMAIN,*)
      write(IMAIN,*) 'Saving DG acoustic frame...'
      write(IMAIN,*)
      call flush_IMAIN()
    endif

    write(outputname,'(a,i6.6,a,i6.6,a)') 'frame_DG_',myrank,'_',it_temp,'.bin'
    open(unit=55,file='OUTPUT_FILES/'//trim(outputname),status='unknown',form='unformatted',iostat=ier)
    if (ier /= 0) call exit_MPI(myrank,'Error opening file frame_DG...')

    read(55) rho_DG
    read(55) rhovx_DG
    read(55) rhovz_DG
    read(55) E_DG

    close(55)
    
  endif

  end subroutine read_forward_solution

! ------------------------------------------------------------------------------------

  !subroutine prepare_MPI_DG(my_neighbours_loc)
  subroutine prepare_MPI_DG()
    
    use constants,only: CUSTOM_REAL,NGLLX,NGLLZ
    
    use specfem_par

    use mpi

    implicit none 
    
    include "precision.h"
    
    integer :: iinterface, ipoin, num_interface, iglob, cpt, &
        i, j, ispec, iglob_DG, nb_values, ier
    integer, dimension(nglob_DG, 3) :: link_ij_iglob
    integer, dimension(nglob_DG) :: MPI_iglob
    
    !integer, dimension(ninterface) :: my_neighbours_loc
    
    !integer, dimension(max_ibool_interfaces_size_ac,ninterface_acoustic) :: &
    !    buffer_recv_faces_vector_DG_i, &
    !    buffer_send_faces_vector_DG_i, &
    !    buffer_recv_faces_vector_DG_j, &
    !    buffer_send_faces_vector_DG_j
    
    integer, dimension(NGLLX*max_interface_size,ninterface) :: &
        buffer_recv_faces_vector_DG_i, &
        buffer_send_faces_vector_DG_i, &
        buffer_recv_faces_vector_DG_j, &
        buffer_send_faces_vector_DG_j    
    
    integer :: i_2, j_2
    
   ! character(len=100) file_name
    !write(file_name,"('./boundaries_elastic2_MPI_',i3.3)") myrank
    ! Open output forcing file
    !open(100,file=file_name,form='formatted')
  
    allocate(MPI_transfer(nglob_DG, 2, 4))
    !allocate(diag_MPI(max_ibool_interfaces_size_ac,ninterface_acoustic))
    
    do ispec = 1, nspec
    
    do i = 1, NGLLX
    do j = 1, NGLLZ
    
        iglob_DG = ibool_DG(i, j, ispec)
    
        link_ij_iglob(iglob_DG,1) = i
        link_ij_iglob(iglob_DG,2) = j
        link_ij_iglob(iglob_DG,3) = ispec
    
    enddo
    enddo
    
    enddo
    
    buffer_recv_faces_vector_DG_i = -1
    buffer_send_faces_vector_DG_i = -1  
    buffer_recv_faces_vector_DG_j = -1
    buffer_send_faces_vector_DG_j = -1   
    
    ! MPI SEND INFO ABOUT DIAG ELEMENT OR NOT
    do iinterface = 1, ninterface_acoustic_DG

    ! gets interface index in the range of all interfaces [1,ninterface]
    num_interface = inum_interfaces_acoustic_DG(iinterface)
    
    ! loops over all interface points
    do ipoin = 1, nibool_interfaces_acoustic_DG(num_interface)
    
        iglob = ibool_interfaces_acoustic_DG(ipoin,num_interface)
        
        i = link_ij_iglob(iglob,1)
        j = link_ij_iglob(iglob,2)
        
        buffer_send_faces_vector_DG_i(ipoin,iinterface) = i
        buffer_send_faces_vector_DG_j(ipoin,iinterface) = j
        
    enddo
    
    enddo
    
    do iinterface = 1, ninterface_acoustic_DG

    ! gets interface index in the range of all interfaces [1,ninterface]
    num_interface = inum_interfaces_acoustic_DG(iinterface)
        
    nb_values = nibool_interfaces_acoustic_DG(num_interface)
    
    call MPI_ISEND( buffer_send_faces_vector_DG_i(1,iinterface), &
             nb_values, MPI_INTEGER, &
             my_neighbours(num_interface), 14, MPI_COMM_WORLD, &
             tab_requests_send_recv_DG(iinterface), ier)
        
    if (ier /= MPI_SUCCESS) then
      call exit_MPI(myrank,'MPI_ISEND unsuccessful in assemble_MPI_vector_start')
    endif

    ! starts a non-blocking receive
    call MPI_IRECV ( buffer_recv_faces_vector_DG_i(1,iinterface), &
             nb_values, MPI_INTEGER, &
             my_neighbours(num_interface), 14, MPI_COMM_WORLD, &
             tab_requests_send_recv_DG(ninterface_acoustic_DG+iinterface), ier)
     
    if (ier /= MPI_SUCCESS) then
      call exit_MPI(myrank,'MPI_IRECV unsuccessful in assemble_MPI_vector')
    endif
    
    call MPI_ISEND( buffer_send_faces_vector_DG_j(1,iinterface), &
             nb_values, MPI_INTEGER, &
             my_neighbours(num_interface), 14, MPI_COMM_WORLD, &
             tab_requests_send_recv_DG(ninterface_acoustic_DG*2+iinterface), ier)
        
    if (ier /= MPI_SUCCESS) then
      call exit_MPI(myrank,'MPI_ISEND unsuccessful in assemble_MPI_vector_start')
    endif

    ! starts a non-blocking receive
    call MPI_IRECV ( buffer_recv_faces_vector_DG_j(1,iinterface), &
             nb_values, MPI_INTEGER, &
             my_neighbours(num_interface), 14, MPI_COMM_WORLD, &
             tab_requests_send_recv_DG(ninterface_acoustic_DG*3+iinterface), ier)
     
    if (ier /= MPI_SUCCESS) then
      call exit_MPI(myrank,'MPI_IRECV unsuccessful in assemble_MPI_vector')
    endif
    
    WRITE(*,*) "-------------->", myrank, iinterface, my_neighbours(num_interface), num_interface, nb_values
    
    enddo
    
    ! waits for MPI requests to complete (recv)
    ! each wait returns once the specified MPI request completed
    do iinterface = 1, 4*ninterface_acoustic_DG
    !do iinterface = 1, ninterface_acoustic_DG
      call MPI_Wait(tab_requests_send_recv_DG(iinterface), &
                  MPI_STATUS_IGNORE, ier)
    
      !call MPI_Wait(tab_requests_send_recv_DG(iinterface + 3*ninterface_acoustic_DG), &
      !            MPI_STATUS_IGNORE, ier)
    enddo
    
    MPI_transfer = -1
    MPI_iglob    = 0
    
    do iinterface = 1, ninterface_acoustic_DG

    ! gets interface index in the range of all interfaces [1,ninterface]
    num_interface = inum_interfaces_acoustic_DG(iinterface)

    ! loops over all interface points
    do ipoin = 1, nibool_interfaces_acoustic_DG(num_interface)
    
       ! WRITE(*,*) "-------------->", myrank, ipoin, iinterface
    
        iglob = ibool_interfaces_acoustic_DG(ipoin,num_interface)
        
        i = link_ij_iglob(iglob,1)
        j = link_ij_iglob(iglob,2)
        ispec = link_ij_iglob(iglob,3)
        
        !if(ispec_is_acoustic_coupling_el(i,j,ispec,3) >= 0) then
        !        WRITE(100,*) coord(:,ibool_before_perio(i, j, ispec))
        !endif
        
        i_2     = buffer_recv_faces_vector_DG_i(ipoin, iinterface)
        j_2     = buffer_recv_faces_vector_DG_j(ipoin, iinterface)
        
        !WRITE(100,*) myrank,i,j,ispec,coord(:,ibool(i, j, ispec))
        
        if( (neighbor_DG(i,j,ispec,3) == -1 .OR. &
                neighbor_DG_corner(i,j,ispec,3) == -1) .AND. (i == i_2 .OR. j == j_2) ) then
        
        !WRITE(*,*) myrank,"------------>", i,j,ispec,i_2,j_2,coord(:,ibool(i, j, ispec))
        !WRITE(*,*) "-----", neighbor_DG(i,j,ispec,3), neighbor_DG_corner(i,j,ispec,3)
        !WRITE(*,*) "***********"
        
        MPI_iglob(iglob) = MPI_iglob(iglob) + 1
        cpt = MPI_iglob(iglob)
        MPI_transfer(iglob,cpt,1) = ipoin 
        MPI_transfer(iglob,cpt,2) = num_interface 
        
        MPI_transfer(iglob,cpt,3) = i_2
        MPI_transfer(iglob,cpt,4) = j_2
       
       endif
        
    enddo
    
    enddo
    
    !do iinterface = 1, ninterface_acoustic_DG
    !  call MPI_Wait(tab_requests_send_recv_DG(iinterface + 0), &
    !              MPI_STATUS_IGNORE, ier)
    !
    !  call MPI_Wait(tab_requests_send_recv_DG(iinterface + 2*ninterface_acoustic_DG), &
    !              MPI_STATUS_IGNORE, ier)
    !enddo
    
    !close(100)
    
  end subroutine prepare_MPI_DG

  ! Setup Vandermonde matrices and derivation matrix for slope limiter
  ! purposes
  subroutine setUpVandermonde()
    
    use constants,only: CUSTOM_REAL,NGLLX,NGLLZ
    
    use specfem_par,only: Vandermonde, invVandermonde, Drx, Drz

    implicit none 
    
    call compute_Vander_matrices(NGLLX*NGLLZ,&
       Vandermonde, invVandermonde, Drx, Drz )
       
       !WRITE(*,*) "maxval Vandermonde", maxval(Vandermonde), minval(Vandermonde)
       
  end subroutine setUpVandermonde
  
  subroutine compute_Vander_matrices(np,&
       V1D, V1D_inv, Drx, Drz )

    use constants,only: CUSTOM_REAL,NGLLX,NGLLZ,GAUSSALPHA,GAUSSBETA
    use specfem_par,only: xigll, zigll!, hprime_xx, hprime_zz

    implicit none 

  ! Parameters
  real(kind=CUSTOM_REAL), parameter :: ZEROl = 0._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: ONEl  = 1._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: TWOl  = 2._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: SIXl  = 6._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: HALFl = 0.5_CUSTOM_REAL

  !double precision, parameter :: zero=0.d0,one=1.d0,two=2.d0

  integer np
  !real(kind=CUSTOM_REAL), parameter :: threshold = 0.000000001_CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: threshold = 0.00001_CUSTOM_REAL
  real(kind=CUSTOM_REAL) V1D(np,np), V1D_inv(np, np),&
         Drx(np, np), Drz(np, np)!,  diagMass(np,np)!, V1D_der_x(np, np),&
          !V1D_der_z(np, np), V1D_t(np, np)!, diagMass(np,np)!, R(np, np)
  double precision, dimension(np,np) :: V1D_d, V1D_inv_d
  !double precision, dimension(NGLLX) :: hx_d, hprimex_d
  !double precision, dimension(NGLLZ) :: hz_d, hprimez_d   
  !real(kind=CUSTOM_REAL), dimension(NGLLX) :: hx, hprimex
  !real(kind=CUSTOM_REAL), dimension(NGLLZ) :: hz, hprimez  

  double precision :: p,pd,pm1,pdm1,pm2,pdm2,p2,pd2

  integer i,j,k,errorflag,l,m,n
  
  !double precision, external :: hgll

  !!!!!!!!!!!!!!!!!!!!!!!
  !!!! => NEEDS TO BE FINISHED
  !!!!!!!!!!!!!!!!!!!!!!!
  ! Init
  V1D     = ZEROl
  !V1D_t   = ZEROl
  V1D_inv = ZEROl
  Drx     = ZEROl
  Drz     = ZEROl
  V1D_d   = 0d0
  
  ! NEW MATRIcES
 ! do k=1,np
 k = 0
  do m = 1,NGLLX
    do n = 1,NGLLZ
    k = k+1
  l = 0
  do i = 1,NGLLX
    do j = 1,NGLLZ
      l = l + 1
      !WRITE(*,*) "1*********"
      !call lagrange_any(xigll(j),NGLLX,xigll,hx_d,hprimex_d)
      !hx_d(i) = hgll(i-1,xigll(j-1),xigll,NGLLX) !lagrange_deriv_GLL(i1-1,i2-1,xigll,NGLLX)
      call jacobf(p,pd,pm1,pdm1,pm2,pdm2,i-1,GAUSSALPHA,GAUSSBETA,xigll(m))
      call jacobf(p2,pd2,pm1,pdm1,pm2,pdm2,j-1,GAUSSALPHA,GAUSSBETA,zigll(n))
      !call jacobf(p,pd,pm1,pdm1,pm2,pdm2,m-1,GAUSSALPHA,GAUSSBETA,xigll(i))
      !call jacobf(p2,pd2,pm1,pdm1,pm2,pdm2,n-1,GAUSSALPHA,GAUSSBETA,zigll(j))
  
      !WRITE(*,*) "TEST LAGRANGE h", h(i2)
      !WRITE(*,*) "TEST LAGRANGE hprime", hprimex_d(i), xigll(j),NGLLX
      !WRITE(*,*) "TEST LAGRANGE hprime_xx(i2,i1)", hprime_xx(i2,i1)
      !WRITE(*,*) "*********"
      !hprimewgll_xx(i2,i1) = wxgll(i2) * hprime_xx(i2,i1)
      !Drx(k,j) = REAL(hprimex_d(i), kind=CUSTOM_REAL)
      !V1D(k,j) = REAL(hx_d(i), kind=CUSTOM_REAL)
      !V1D_der_x(k,l) = REAL(pd, kind=CUSTOM_REAL)
      !V1D_der_z(k,l) = REAL(pd2, kind=CUSTOM_REAL)
      V1D(k,l) = REAL(p*p2, kind=CUSTOM_REAL)
      if(abs(V1D(k,l)) < threshold) V1D(k,l) = ZEROl
      !WRITE(*,*) "V1D(",k,",",l,") = ", V1D(k,l)
      V1D_d(k,l) = p*p2
      !V1D_t(l,k) = V1D(k,l)
    enddo
  enddo
  !enddo
    enddo
  enddo
  
  call FINDInv(V1D_d, V1D_inv_d, np, errorflag)
  do j=1,np
    do i=1,np
        V1D_inv(i,j) = REAL(V1D_inv_d(i,j), kind=CUSTOM_REAL)
        if(abs(V1D_inv(i,j)) < threshold) V1D_inv(i,j) = ZEROl
    enddo
  enddo
  
  !do j=1,np
  !  do i=1,np
  !      do k=1,np
  !      Drx(i,j) = Drx(i,j) + V1D_der_x(i,k)*V1D_inv(k,j)
  !      !WRITE(*,*) "A(",k,",",l,") = ", V1D_der_x(i,k), REAL(V1D_inv(k,j), kind=CUSTOM_REAL)
  !      Drz(i,j) = Drz(i,j) + V1D_der_z(i,k)*V1D_inv(k,j)
  !      enddo
  !     ! WRITE(*,*) "Dr", Dr(i,j) 
  !  enddo
  !enddo
  
  !TEST INVERSION
  !R = ZEROl
  !do k=1,np
  !      do l=1,np
  !              do m=1,np
  !              R(k,l)  = R(k,l) + V1D(k,m)*V1D_inv(m,l)
  !              enddo
  !              WRITE(*,*) "R(k,l)",k,l,R(k,l)
  !    enddo
  !    
  !enddo

  end subroutine compute_Vander_matrices
  
  !Subroutine to find the inverse of a square matrix
!Author : Louisda16th a.k.a Ashwith J. Rego
!Reference : Algorithm has been well explained in:
!http://math.uww.edu/~mcfarlat/inverse.htm 
!http://www.tutor.ms.unimelb.edu.au/matrix/matrix_inverse.html

SUBROUTINE FINDInv(matrix, inverse, n, errorflag)

use constants,only: CUSTOM_REAL

IMPLICIT NONE
!Declarations
INTEGER, INTENT(IN) :: n
INTEGER, INTENT(OUT) :: errorflag !Return error status. -1 for error, 0 for normal
double precision, INTENT(IN), DIMENSION(n,n) :: matrix !Input matrix
double precision, INTENT(OUT), DIMENSION(n,n) :: inverse !Inverted matrix

LOGICAL :: FLAG = .TRUE.
INTEGER :: i, j, k
double precision :: m
double precision, DIMENSION(n,2*n) :: augmatrix !augmented matrix

!Augment input matrix with an identity matrix
DO i = 1, n
DO j = 1, 2*n
IF (j <= n ) THEN
augmatrix(i,j) = matrix(i,j)
ELSE IF ((i+n) == j) THEN
augmatrix(i,j) = 1
Else
augmatrix(i,j) = 0
ENDIF
END DO
END DO

!Reduce augmented matrix to upper traingular form
DO k =1, n-1
IF (augmatrix(k,k) == 0) THEN
        FLAG = .FALSE.
        DO i = k+1, n
                        IF (augmatrix(i,k) /= 0) THEN
                        DO j = 1,2*n
                                augmatrix(k,j) = augmatrix(k,j)+augmatrix(i,j)
                        END DO
                        FLAG = .TRUE.
                        EXIT
                ENDIF
                IF (FLAG .EQV. .FALSE.) THEN
                        PRINT*, "Matrix is non - invertible"
                        inverse = 0
                        errorflag = -1
                        return
                ENDIF
        END DO
ENDIF

DO j = k+1, n 
        m = augmatrix(j,k)/augmatrix(k,k)
        DO i = k, 2*n
                augmatrix(j,i) = augmatrix(j,i) - m*augmatrix(k,i)
        END DO
END DO

END DO

!Test for invertibility
DO i = 1, n
IF (augmatrix(i,i) == 0) THEN
PRINT*, "Matrix is non - invertible"
inverse = 0
errorflag = -1
return
ENDIF
END DO

!Make diagonal elements as 1
DO i = 1 , n
m = augmatrix(i,i)
DO j = i , (2 * n) 
augmatrix(i,j) = (augmatrix(i,j) / m)
END DO
END DO

!Reduced right side half of augmented matrix to identity matrix
DO k = n-1, 1, -1
DO i =1, k
m = augmatrix(i,k+1)
DO j = k, (2*n)
augmatrix(i,j) = augmatrix(i,j) -augmatrix(k+1,j) * m
END DO
END DO
END DO 

!store answer
DO i =1, n
DO j = 1, n
inverse(i,j) = augmatrix(i,j+n)
END DO
END DO
errorflag = 0
RETURN
END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine build_veloc_boundary_DG(veloc_vector_acoustic_DG_coupling)    
    
    use constants,only: CUSTOM_REAL,NGLLX,NGLLZ
    
    use specfem_par,only: ispec_is_acoustic, ispec_is_acoustic_DG, nspec, nglob, ibool, potential_dot_acoustic, &
        hprime_xx, hprime_zz, xix, xiz, gammax, gammaz

    implicit none 
    
    real(kind=CUSTOM_REAL), dimension(nglob,2) :: veloc_vector_acoustic_DG_coupling
    integer :: i, j, k, ispec
    real(kind=CUSTOM_REAL) :: xixl, xizl, gammaxl, gammazl, dux_dxi, dux_dgamma

    veloc_vector_acoustic_DG_coupling = 0.

    do ispec = 1,nspec

    ! acoustic spectral element
    if (ispec_is_acoustic(ispec) .AND. .not. ispec_is_acoustic_DG(ispec)) then

      ! first double loop over GLL points to compute and store gradients
      do j = 1,NGLLZ
        do i = 1,NGLLX

            xixl = xix(i, j, ispec)
            xizl = xiz(i, j, ispec)
            gammaxl = gammax(i, j, ispec)
            gammazl = gammaz(i, j, ispec)

            do k = 1,NGLLX
                dux_dxi    = dux_dxi    + potential_dot_acoustic(ibool(k,j,ispec)) * hprime_xx(i,k)
                dux_dgamma = dux_dgamma + potential_dot_acoustic(ibool(i,k,ispec)) * hprime_zz(j,k)
            enddo

            ! derivatives of potential
            veloc_vector_acoustic_DG_coupling(ibool(i, j, ispec), 1) = dux_dxi * xixl + dux_dgamma * gammaxl
            veloc_vector_acoustic_DG_coupling(ibool(i, j, ispec), 2) = dux_dxi * xizl + dux_dgamma * gammazl
        enddo
     enddo
     
    endif
    
    enddo

  end subroutine 


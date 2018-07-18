module specfem_par_LNS

  use constants

  implicit none
  
  ! Switch enabling the use of LNS instead of FNS. To be added in parfile later.
  logical, parameter :: USE_LNS = .true.
  
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: LNS_rho0, LNS_v0x,     LNS_v0z,     LNS_E0 ! Initial state.
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: LNS_drho, LNS_rho0dvx, LNS_rho0dvz, LNS_dE ! State.
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: RHS_drho, RHS_rho0dvx, RHS_rho0dvz, RHS_dE ! RHS.
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: aux_drho, aux_rho0dvx, aux_rho0dvz, aux_dE ! Auxiliary registers.
  
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: nabla_v0
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: LNS_p0, LNS_T0
  
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: sigma_v_0 ! Initial viscous stress tensor.
  
  integer(4), parameter :: LNS_VERBOSE = 31 ! Verbosity parameter, ranging from 0 to 31.
  
  ! MPI: Transfers' buffers.
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable  :: buffer_LNS_drho_P, buffer_LNS_rho0dvx_P, &
                                                          buffer_LNS_rho0dvz_P, buffer_LNS_dE_P

end module specfem_par_LNS

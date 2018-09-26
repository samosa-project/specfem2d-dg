module specfem_par_LNS

  use constants

  implicit none
  
  ! Switch enabling the use of LNS instead of FNS. To be added in parfile later.
  logical, parameter :: USE_LNS = .false.
  
  ! 2D/3D generalisation pre-work.
  integer(kind=selected_int_kind(1)), parameter :: SPACEDIM = 2 ! Spatial dimension (for later generalisation) Min/Maximum values: [-10^1+1=-9, 10^1-1=9] (ok since we only need 2 or 3, and maybe 1).
  
  ! Physical parameters.
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: LNS_g, LNS_mu, LNS_eta, LNS_kappa
  
  ! State registers.
  real(kind=CUSTOM_REAL), dimension(:),   allocatable :: LNS_rho0, LNS_E0 ! Initial state.
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: LNS_v0 ! Initial state.
  real(kind=CUSTOM_REAL), dimension(:),   allocatable :: LNS_drho, LNS_dE ! State.
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: LNS_rho0dv ! State.
  real(kind=CUSTOM_REAL), dimension(:),   allocatable :: LNS_dp,LNS_dT ! Pressure and temperature perturbation.
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: LNS_dm ! Momentum (1st order) perturbation.
  real(kind=CUSTOM_REAL), dimension(:),   allocatable :: RHS_drho, RHS_dE ! RHS.
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: RHS_rho0dv ! RHS.
  real(kind=CUSTOM_REAL), dimension(:),   allocatable :: aux_drho, aux_dE ! Auxiliary registers.
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: aux_rho0dv ! Auxiliary registers.
  
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: nabla_v0 ! Gradient of initial velocity.
  real(kind=CUSTOM_REAL), dimension(:),     allocatable :: LNS_p0, LNS_T0 ! Initial pressure and temperature.
  
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: sigma_v_0 ! Initial viscous stress tensor. Symmetric, thus only need to save 3 entries: (1,:) index should correspond to the index (1,1) of the actual tensor, (2,:) <-> (1,2) and (2,1) of the actual tensor, and (3,:) <-> (2,2) of the actual tensor.
  
  integer(kind=selected_int_kind(2)), parameter :: LNS_VERBOSE = 99 ! Verbosity parameter. Min/Maximum values: [-10^2+1=-99, 10^2-1=99].
  ! LNS_VERBOSE>= 1: printing iteration, stages, and local times, every 100 iterations
  ! LNS_VERBOSE>=51: printing min/max values for each constitutive variable on CPU 0 every 100 iterations
  
  ! MPI: Transfers' buffers.
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable  :: buffer_LNS_drho_P, buffer_LNS_dE_P
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable  :: buffer_LNS_rho0dv_P

  ! Dummy variables.
  real(kind=CUSTOM_REAL), dimension(:),     allocatable :: LNS_dummy_1d ! Dummy variables are not optimal, but prevent from duplicating subroutines.
  real(kind=CUSTOM_REAL), dimension(:,:),   allocatable :: LNS_dummy_2d ! Dummy variables are not optimal, but prevent from duplicating subroutines.
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: LNS_dummy_3d ! Dummy variables are not optimal, but prevent from duplicating subroutines.
end module specfem_par_LNS

module specfem_par_LNS

  use constants

  implicit none
  
  ! Switch enabling the use of LNS instead of FNS.
  logical :: USE_LNS
  
  ! 2D/3D generalisation pre-work.
  integer(kind=selected_int_kind(1)), parameter :: SPACEDIM = 2 ! Spatial dimension (for later generalisation) Min/Maximum values: [-10^1+1=-9, 10^1-1=9] (ok since we only need 2 or 3, and maybe 1).
  integer, parameter :: NVALSIGMA=int(0.5*SPACEDIM*(SPACEDIM+1)) ! Number of distinct values in the symmetric viscous tensor.
  
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
  
  ! Initial quantities.
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: nabla_v0 ! Gradient of initial velocity.
  real(kind=CUSTOM_REAL), dimension(:),     allocatable :: LNS_p0, LNS_T0 ! Initial pressure and temperature.
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: sigma_v_0 ! Initial viscous stress tensor. Symmetric, thus only need to save few (NVALSIGMA) entries. In 2D : (1,:) index must correspond to the index (1,1) of the actual tensor, (2,:) <-> (1,2) and (2,1) of the actual tensor, and (3,:) <-> (2,2) of the actual tensor.
  
  ! Auxiliary quantities.
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: nabla_dT ! Gradient of temperature perturbation.
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: sigma_dv ! Viscous stress tensor perturbation. Symmetric, thus only need to save few entries (see sigma_v_0).
  
  integer(kind=selected_int_kind(2)), parameter :: LNS_VERBOSE = 99 ! Verbosity parameter. Min/Maximum values: [-10^2+1=-99, 10^2-1=99].
  ! LNS_VERBOSE>= 1: printing iteration, stages, and local times, every LNS_MODPRINT iterations
  ! LNS_VERBOSE>=51: printing min/max values for each constitutive variable on CPU 0 every LNS_MODPRINT iterations
  
  ! Integer used in modulo for condition on printing min/max values for each constitutive variable on CPU 0 every 100 iterations. To be added in parfile later.
  integer, parameter :: LNS_MODPRINT = 500
  
  ! MPI: Transfers' buffers.
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable  :: buffer_LNS_drho_P, buffer_LNS_dE_P
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable  :: buffer_LNS_rho0dv_P
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable  :: buffer_LNS_nabla_dT, buffer_LNS_sigma_dv

  ! Dummy variables.
  real(kind=CUSTOM_REAL), dimension(:),     allocatable :: LNS_dummy_1d ! Dummy variables are not optimal, but prevent from duplicating subroutines.
  real(kind=CUSTOM_REAL), dimension(:,:),   allocatable :: LNS_dummy_2d ! Dummy variables are not optimal, but prevent from duplicating subroutines.
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: LNS_dummy_3d ! Dummy variables are not optimal, but prevent from duplicating subroutines.
  
  contains
  
  ! Redefine the Fortran 2008 routine norm2, alongside one dedicated to rank 1 vectors.
  function norm2(r2arr)
    use constants, only: CUSTOM_REAL
    implicit none
    ! Input/Output.
    real(kind=CUSTOM_REAL), dimension(:,:), intent(in) :: r2arr ! Input: rank 2 array (tensor field).
    real(kind=CUSTOM_REAL), dimension(size(r2arr,2)) :: norm2
    ! Local.
    !N./A.
    norm2 = sum(r2arr**2,1)
  end function norm2
  function norm2r1(r1arr)
    use constants, only: CUSTOM_REAL
    implicit none
    ! Input/Output.
    real(kind=CUSTOM_REAL), dimension(:), intent(in) :: r1arr ! Input: rank 1 array (vector).
    real(kind=CUSTOM_REAL) :: norm2r1
    ! Local.
    !N./A.
    norm2r1 = sum(r1arr**2)
  end function norm2r1
end module specfem_par_LNS

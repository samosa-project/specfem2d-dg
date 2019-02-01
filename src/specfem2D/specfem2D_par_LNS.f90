module specfem_par_LNS

  use constants

  implicit none
  
  ! Constants (maybe move to constants.h.in).
  real(kind=CUSTOM_REAL), parameter :: R_ADIAB = 8.3144598
  
  ! LSRK coefficients.
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: LNS_scheme_A, LNS_scheme_B, LNS_scheme_C ! Size should be stage_time_scheme.
  
  ! Switch enabling the use of LNS instead of FNS.
  logical :: USE_LNS
  
  logical, parameter :: LNS_switch_gradient=.false. ! Switch to activate the use of the "desintegration method" for gradient computation methods, and to desactivate to use the SEM definition of the gradient. ! Warning: LNS_switch_gradient = .true. is not yet fully implemented.
  
  ! 2D/3D generalisation pre-work.
  !integer(kind=selected_int_kind(1)), parameter :: SPACEDIM = 2 ! Spatial dimension (for later generalisation) Min/Maximum values: [-10^1+1=-9, 10^1-1=9] (ok since we only need 2 or 3, and maybe 1).
  integer, parameter :: NVALSIGMA=int(0.5*NDIM*(NDIM+1)) ! Number of distinct values in the symmetric viscous tensor.
  
  ! Physical parameters.
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: LNS_g, LNS_mu, LNS_eta, LNS_kappa, LNS_c0
  logical :: LNS_viscous ! General switch being true if (maxval(LNS_mu) > 0. .OR. maxval(LNS_eta) > 0. .OR. maxval(LNS_kappa) > 0.) and false if not. See compute_forces_acoustic_LNS_calling_routine. In the latter case, enables faster verification and thus faster skipping of some parts of the code.
  
  ! State registers.
  ! Pretty much all these arrays are allocated in 'prepare_timerun_wavefields.f90'.
  real(kind=CUSTOM_REAL), dimension(:),   allocatable :: LNS_rho0, LNS_E0 ! Initial state.
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: LNS_v0 ! Initial state.
  real(kind=CUSTOM_REAL), dimension(:),   allocatable :: LNS_drho, LNS_dE ! State.
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: LNS_rho0dv ! State.
  real(kind=CUSTOM_REAL), dimension(:),   allocatable :: LNS_dp,LNS_dT ! Pressure and temperature perturbation.
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: LNS_dv, LNS_dm ! Velocity perturbation, and momentum (1st order) perturbation.
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
  
  ! PMLs.
  ! Pretty much all these arrays are allocated in 'prepare_timerun_pml.f90'.
  integer :: nglob_PML ! Number of PML points (spatial duplicates included).
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: rmass_inverse_acoustic_LNS_PML ! Inverse mass matrix. Size allocated should be (nglob_PML).
  integer, dimension(:,:,:), allocatable :: ibool_LNS_PML ! Same as ibool_DG (see 'specfem2D_par'), but for PML only.
  !real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: LNS_PML_drho, LNS_PML_dE ! Size allocated should be (NADE, NGLLX, NGLLZ, nspec_PML).
  !real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: LNS_PML_rho0dv ! Size allocated should be (NADE, NDIM, NGLLX, NGLLZ, nspec_PML).
  !real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: RHS_PML_drho,aux_PML_drho, RHS_PML_dE,aux_PML_dE
  !real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: RHS_PML_rho0dv, aux_PML_rho0dv
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: LNS_PML_drho, LNS_PML_dE ! Size allocated should be (NADE, NGLLX*NGLLZ*nspec_PML).
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: LNS_PML_rho0dv ! Size allocated should be (NADE, NDIM, NGLLX*NGLLZ*nspec_PML).
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: RHS_PML_drho,aux_PML_drho, RHS_PML_dE,aux_PML_dE
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: RHS_PML_rho0dv, aux_PML_rho0dv
  

  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: LNS_PML_kapp ! Constant coefficient in the stretching s. Size allocated should be (NDIM, NGLLX, NGLLZ, nspec_PML).
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: LNS_PML_alpha ! Coefficient in front of the auxiliary variables. For classical formulation, only 2=NDIM ADE are to be solved for each variable, hence the first dimension.
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: LNS_PML_a0 ! Coefficient in front of the \delta, that is in front of the q in the updated strong form. Size allocated should be (NGLLX,NGLLZ,nspec_PML) in order to save memory.
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: LNS_PML_b ! Coefficient in front of each auxiliary variable (ADEs). For classical formulation, only 2=NDIM ADE are to be solved for each variable, hence the first dimension.
  
  integer(kind=selected_int_kind(2)), parameter :: LNS_VERBOSE = 99 ! Verbosity parameter. Min/Maximum values: [-10^2+1=-99, 10^2-1=99].
  ! LNS_VERBOSE>= 1: printing iteration, stages, and local times, every LNS_MODPRINT iterations
  ! LNS_VERBOSE>=51: printing min/max values for each constitutive variable on CPU 0 every LNS_MODPRINT iterations
  
  ! Integer used in modulo for condition on printing min/max values for each constitutive variable on CPU 0 every 100 iterations. To be added in parfile later.
  integer, parameter :: LNS_MODPRINT = 500
  
  ! MPI: Transfers' buffers.
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: buffer_LNS_drho_P, buffer_LNS_dE_P
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: buffer_LNS_rho0dv_P
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: buffer_LNS_nabla_dT, buffer_LNS_sigma_dv

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

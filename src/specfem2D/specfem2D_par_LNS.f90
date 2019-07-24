module specfem_par_LNS

  use constants

  implicit none
  
  ! ---------------------------- !
  ! Switches.
  ! ---------------------------- !
  ! Switch enabling the use of LNS instead of FNS.
  logical :: USE_LNS
  
  ! Switch to activate the use of the "desintegration method" for gradient computation methods, and to desactivate to use the SEM definition of the gradient. ! Warning: LNS_switch_gradient = .true. is not yet fully implemented.
  logical, parameter :: LNS_switch_gradient = .false.
  
  ! General switch being true if (maxval(LNS_mu) > 0. .OR. maxval(LNS_eta) > 0. .OR. maxval(LNS_kappa) > 0.) and false if not. See compute_forces_acoustic_LNS_calling_routine. In the latter case, enables faster verification and thus faster skipping of some parts of the code.
  logical :: LNS_viscous
  ! ------------------------------------------------------------ !
  
  ! ---------------------------- !
  ! Constants and coefficients.
  ! ---------------------------- !
  ! Constants (TODO maybe move to constants.h.in).
  real(kind=CUSTOM_REAL), parameter :: R_ADIAB = 8.3144598
  
  ! LSRK coefficients.
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: LNS_scheme_A, LNS_scheme_B, LNS_scheme_C ! Size should be stage_time_scheme.
  ! ------------------------------------------------------------ !
  
  ! ---------------------------- !
  ! 2D/3D generalisation
  ! pre-work.
  ! ---------------------------- !
  integer, parameter :: NVALSIGMA = int(0.5*NDIM*(NDIM+1)) ! Number of distinct values in the symmetric viscous tensor.
  ! ------------------------------------------------------------ !
  
  ! ---------------------------- !
  ! Printing parameters.
  ! ---------------------------- !
  integer(kind=selected_int_kind(2)), parameter :: LNS_VERBOSE = 99 ! Verbosity parameter. Min/Maximum values: [-10^2+1=-99, 10^2-1=99].
  ! LNS_VERBOSE>= 1: printing iteration, stages, and local times, every LNS_MODPRINT iterations
  ! LNS_VERBOSE>=51: printing min/max values for each constitutive variable on CPU 0 every LNS_MODPRINT iterations
  
  ! Integer used in modulo for condition on printing min/max values for each constitutive variable on CPU 0 every 100 iterations. To be added in parfile later.
  integer, parameter :: LNS_MODPRINT = 500
  ! ------------------------------------------------------------ !
  
  ! ---------------------------- !
  ! Physical quantities.
  ! ---------------------------- !
  ! Pretty much all these arrays are allocated in 'prepare_timerun_wavefields.f90'.
  ! Physical parameters.
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: LNS_g, LNS_mu, LNS_eta, LNS_kappa, LNS_c0
  
  ! Initial state.
  real(kind=CUSTOM_REAL), dimension(:),   allocatable :: LNS_rho0, LNS_E0
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: LNS_v0
  
  ! Auxiliary initial state quantities.
  real(kind=CUSTOM_REAL), dimension(:),     allocatable :: LNS_p0, LNS_T0
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: nabla_v0 ! Gradient of initial velocity.
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: sigma_v_0 ! Initial viscous stress tensor. Symmetric, thus only need to save few (NVALSIGMA) entries. In 2D : (1,:) index must correspond to the index (1, 1) of the actual tensor, (2, :) <-> (1, 2) and (2, 1) of the actual tensor, and (3, :) <-> (2, 2) of the actual tensor.
  
  ! State.
  real(kind=CUSTOM_REAL), dimension(:),   allocatable :: LNS_drho, LNS_dE
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: LNS_rho0dv
  
  ! Auxiliary state quantities.
  real(kind=CUSTOM_REAL), dimension(:),   allocatable :: LNS_dp, LNS_dT
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: LNS_dv, LNS_dm ! Velocity perturbation, and momentum (1st order) perturbation.
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: nabla_dT ! Gradient of temperature perturbation.
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: sigma_dv ! Viscous stress tensor perturbation. Symmetric, thus only need to save few entries (see sigma_v_0).
  ! ------------------------------------------------------------ !
  
  ! ---------------------------- !
  ! LSRK time iteration.
  ! ---------------------------- !
  ! RHS registers.
  real(kind=CUSTOM_REAL), dimension(:),   allocatable :: RHS_drho, RHS_dE
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: RHS_rho0dv
  
  ! Auxiliary registers.
  real(kind=CUSTOM_REAL), dimension(:),   allocatable :: aux_drho, aux_dE
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: aux_rho0dv
  ! ------------------------------------------------------------ !
  
  ! ---------------------------- !
  ! Absorbing boundary
  ! conditions.
  ! ---------------------------- !
  ! Real stretching.
  ! See 'specfem2D_par.f90', variables:
  ! ABC_STRETCH_*,
  ! iy_image_color_*_buffer, ix_image_color_*_buffer,
  ! stretching_ya,
  ! stretching_buffer
  
  ! PMLs.
  ! Pretty much all these arrays are allocated in 'prepare_timerun_pml.f90'.
  ! /!\ PMLs are still under developement.
  integer :: nglob_PML ! Number of PML points (spatial duplicates included). Computed in 'pml_init.f90'.
  !real(kind=CUSTOM_REAL), dimension(:), allocatable :: rmass_inverse_acoustic_LNS_PML ! Inverse mass matrix. Size allocated should be (nglob_PML).
  integer, dimension(:,:,:), allocatable :: ibool_LNS_PML ! Same as ibool_DG (see 'specfem2D_par'), but for PML only.
  !real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: LNS_PML_drho, LNS_PML_dE ! Size allocated should be (NADE, NGLLX, NGLLZ, nspec_PML).
  !real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: LNS_PML_rho0dv ! Size allocated should be (NADE, NDIM, NGLLX, NGLLZ, nspec_PML).
  !real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: RHS_PML_drho,aux_PML_drho, RHS_PML_dE,aux_PML_dE
  !real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: RHS_PML_rho0dv, aux_PML_rho0dv
  !real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: LNS_PML_drho, LNS_PML_dE ! Size allocated should be (NADE, NGLLX*NGLLZ*nspec_PML).
  !real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: LNS_PML_rho0dv ! Size allocated should be (NADE, NDIM, NGLLX*NGLLZ*nspec_PML).
  !real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: RHS_PML_drho,aux_PML_drho, RHS_PML_dE,aux_PML_dE
  !real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: RHS_PML_rho0dv, aux_PML_rho0dv
  integer, parameter :: LNS_PML_NAUX=6 ! Number of auxiliary variables per constitutive variable for PMLs.
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: LNS_PML, LNS_PML_RHS, LNS_PML_aux ! Dimensions: 1<-> # constitutive variables, 2<-> # auxiliary variables, 3<-> #PML mesh points.
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: LNS_PML_kapp ! Constant coefficient in the stretching s. Size allocated should be (NDIM, NGLLX, NGLLZ, nspec_PML).
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: LNS_PML_alpha ! Coefficient in front of the auxiliary variables. For classical formulation, only 2=NDIM ADE are to be solved for each variable, hence the first dimension.
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: LNS_PML_a0 ! Coefficient in front of the \delta, that is in front of the q in the updated strong form. Size allocated should be (NGLLX,NGLLZ,nspec_PML) in order to save memory.
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: LNS_PML_b, LNS_PML_d ! Coefficient in front of each auxiliary variable (ADEs). For classical formulation, only 2=NDIM ADE are to be solved for each variable, hence the first dimension.
  ! ------------------------------------------------------------ !
  
  ! ---------------------------- !
  ! Parallel computing.
  ! ---------------------------- !
  ! MPI: Transfers' buffers.
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: buffer_LNS_drho_P, buffer_LNS_dE_P
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: buffer_LNS_rho0dv_P
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: buffer_LNS_nabla_dT, buffer_LNS_sigma_dv
  ! ------------------------------------------------------------ !
  
  ! ---------------------------- !
  ! Dummy variables.
  ! ---------------------------- !
  ! Those are useful for using the 'DG' extension routines without re-coding them.
  real(kind=CUSTOM_REAL), dimension(:),     allocatable :: LNS_dummy_1d ! Dummy variables are not optimal, but prevent from duplicating subroutines.
  real(kind=CUSTOM_REAL), dimension(:,:),   allocatable :: LNS_dummy_2d ! Dummy variables are not optimal, but prevent from duplicating subroutines.
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: LNS_dummy_3d ! Dummy variables are not optimal, but prevent from duplicating subroutines.
  ! ------------------------------------------------------------ !
  
  ! ---------------------------- !
  ! Functions.
  ! ---------------------------- !
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
  
  ! Implement explicitely the inverse of a 2x2 matrix.
  function inverse_2x2(A) result(Am1)
    use constants, only: CUSTOM_REAL
    implicit none
    ! Input/Output.
    real(kind=CUSTOM_REAL), dimension(2,2), intent(in) :: A
    real(kind=CUSTOM_REAL), dimension(2,2) :: Am1
    ! Local.
    !N./A.
    Am1(1, 1) =  A(2, 2)
    Am1(1, 2) = -A(1, 2)
    Am1(2, 1) = -A(2, 1)
    Am1(2, 2) =  A(1, 1)
    Am1 = Am1/(A(1, 1)*A(2, 2)-A(1, 2)*A(2, 1))
  end function inverse_2x2

!  ! ------------------------------------------------------------ !
!  ! RusanovFlux                                                  !
!  ! ------------------------------------------------------------ !
!  ! Implements the Rusanov approximation for the computation of the flux across an interface where the following Riemann problem is posed:
!  !   \partial_t q + \nabla\cdot\Sigma = S
!  function RusanovFlux(Sigma_M, Sigma_P, q_M, q_P, normal_vec, exact_interface_flux, lambda) result(flux)
!    use constants, only: CUSTOM_REAL, NDIM
!    implicit none
!    ! Input/Output.
!    real(kind=CUSTOM_REAL), dimension(NDIM), intent(in) :: Sigma_M, Sigma_P ! tensor field on either side
!    real(kind=CUSTOM_REAL), intent(in) :: q_M, q_P ! constitutive variable on either side
!    real(kind=CUSTOM_REAL), dimension(NDIM), intent(in) :: normal_vec ! normal vector at the interface
!    logical, intent(in) :: exact_interface_flux ! whether to implement a jump or not
!    real(kind=CUSTOM_REAL), intent(in) :: lambda ! a pre-computed approximate maximum linearised acoustic wave speed
!    real(kind=CUSTOM_REAL) :: flux
!    ! Local.
!    !N./A.
!    if(exact_interface_flux) then
!      flux = 0.5_CUSTOM_REAL * DOT_PRODUCT(normal_vec, Sigma_M+Sigma_P)
!    else
!      flux = 0.5_CUSTOM_REAL * (DOT_PRODUCT(normal_vec, Sigma_M+Sigma_P) - lambda*(q_P-q_M))
!    endif
!  end function RusanovFlux
  ! ------------------------------------------------------------ !
  
end module specfem_par_LNS

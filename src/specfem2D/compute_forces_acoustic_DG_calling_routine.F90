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
! but WITHout ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
! The full text of the license is available in file "LICENSE".
!
!========================================================================

! ------------------------------------------------------------ !
! compute_forces_acoustic_DG_main                              !
! ------------------------------------------------------------ !
! Main routine taking care of the DG elements to be solved with the FNS (Full Navier-Stokes) module.
! This routine is called by iterate_time().

subroutine compute_forces_acoustic_DG_main()

  use specfem_par
  use constants, only: rk4a_d, rk4b_d, rk4c_d

  implicit none

  ! Input/Output.
  ! N. A.
  
  ! Local variables.
  real(kind=CUSTOM_REAL), parameter :: ZEROl = 0._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: ONEl = 1._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: threshold = 0.0000001_CUSTOM_REAL
  real(kind=CUSTOM_REAL) :: timelocal
  real(kind=CUSTOM_REAL), dimension(stage_time_scheme) :: scheme_A, scheme_B, scheme_C
  real(kind=CUSTOM_REAL), dimension(nglob_DG) :: veloc_x
  integer :: i, j, ispec
  
  logical CHECK_NONPOSITIVITY, CHECK_NONPOSITIVITY_ON_ALL_PROCS, CHECK_NONPOSITIVITY_FinD_POinT
  
  ! Checks if anything has to be done.
  if (.not. any_acoustic_DG) then
    return
  endif
  
  ! Those are debug switches. They are computationaly very heavy and should not be used on every simulation.
  ! CHECK_NONPOSITIVITY_ON_ALL_PROCS=.true. and CHECK_NONPOSITIVITY_FinD_POinT=.true. are particularly heavy.
  CHECK_NONPOSITIVITY              = .false. ! Set to .true. to enable nonpositivity checking.
  CHECK_NONPOSITIVITY_ON_ALL_PROCS = .false. ! Only used if CHECK_NONPOSITIVITY==.true.. Set to .false. for checking only on proc ZEROl Set to .true. for checking on all procs.
  CHECK_NONPOSITIVITY_FinD_POinT   = .false. ! Only used if CHECK_NONPOSITIVITY==.true.. Set to .true. to find where nonpositivity was encountered.
  
  ! intialisation.
  if(it == 1 .and. i_stage == 1) then
    if(USE_SLOPE_LIMITER) then
      ! The Vandermonde matrices are only used when the slope limiter is.
      call setUpVandermonde()
    endif
    
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
    if(IONOSPHERIC_COUPLING) then
      dot_Ni(:)    = ZEROl
    endif
    
    ! Allocate model arrays.
    if(.not. assign_external_model) then
      deallocate(gravityext, muext, etaext, kappa_DG, tau_epsilon, tau_sigma)
      if(IONOSPHERIC_COUPLING) then
        ! Safeguard: re-allocate.
        deallocate(N0ext)
        allocate(N0ext(NGLLX, NGLLZ, nspec))
      endif
      allocate(gravityext(NGLLX, NGLLZ, nspec), &
               etaext(NGLLX, NGLLZ, nspec), &
               muext(NGLLX, NGLLZ, nspec), &
               kappa_DG(NGLLX, NGLLZ, nspec), &
               tau_epsilon(NGLLX, NGLLZ, nspec), &
               tau_sigma(NGLLX, NGLLZ, nspec))
    endif
    
    if(NPROC > 1) then
      ! Prepare MPI buffers.
      call prepare_MPI_DG()
    endif
    
    ! Allocate acoustic coupling array.
    allocate(ispec_is_acoustic_coupling_ac(nglob_DG))
    ispec_is_acoustic_coupling_ac = -1
    if(.not. only_DG_acoustic) then
      call find_DG_acoustic_coupling()
    endif
    
    ! initialise auxiliary tensors.
    T_DG = ZEROl
    V_DG = ZEROl
    drho_DG = ZEROl
    dEp_DG = ZEROl
    if(IONOSPHERIC_COUPLING) then
      Ni_DG = ZEROl
    endif
    
    call initial_condition_DG()

    ! Allocate and set initial fields.
    allocate(rho_init(nglob_DG), rhovx_init(nglob_DG), rhovz_init(nglob_DG), E_init(nglob_DG))
    rhovx_init = rhovx_DG
    rhovz_init = rhovz_DG
    E_init     = E_DG
    rho_init   = rho_DG
    p_DG_init  = (gammaext_DG - ONE)*( E_DG - (HALF/rho_DG)*( rhovx_DG**2 + rhovz_DG**2 ) )
    T_init = (E_DG/rho_DG - HALF*((rhovx_DG/rho_DG)**2 + (rhovz_DG/rho_DG)**2))/c_V
    
    ! When elastic-DG simulations, p_DG = 0 in elastic elements and 1/p_DG will not be properly defined. Use a hack.
    where(p_DG_init <= ZEROl) p_DG_init = ONEl
    where(rho_init <= ZEROl) rho_init = ONEl
    
#ifdef USE_MPI
    if (NPROC > 1 .and. ninterface_acoustic > 0) then
      call assemble_MPI_vector_DG(gammaext_DG, buffer_DG_gamma_P)
    endif
#endif
    
    if(myrank==0) then
      write(*,*) "min, max RHO0 on proc ", myrank, ": ", minval(rho_init), maxval(rho_init)
      write(*,*) "min, max P0   on proc ", myrank, ": ", minval(p_DG_init), maxval(p_DG_init)
      write(*,*) "min, max E0   on proc ", myrank, ": ", minval(E_init), maxval(E_init)
    endif
    
  endif ! Endif on (it == 1) and (i_stage == 1).
  
  ! Call time scheme coefficients and cast them.
  if (time_stepping_scheme == 3) then
    ! 5 stages.
    scheme_A = real(rk4a_d, kind=CUSTOM_REAL)
    scheme_B = real(rk4b_d, kind=CUSTOM_REAL)
    scheme_C = real(rk4c_d, kind=CUSTOM_REAL)
  else if (time_stepping_scheme == 4) then
    ! 3 stages.
    scheme_A = real(ls33rk_a, kind=CUSTOM_REAL)
    scheme_B = real(ls33rk_b, kind=CUSTOM_REAL)
    scheme_C = real(ls33rk_c, kind=CUSTOM_REAL)
  else
    if(myrank==0) then
      write(*,*) "********************************"
      write(*,*) "*            ERROR             *"
      write(*,*) "********************************"
      write(*,*) "* This time_stepping_scheme is *"
      write(*,*) "* not implemented for DG yet.  *"
      write(*,*) "* time_stepping_scheme = ", time_stepping_scheme
      write(*,*) "********************************"
      stop
    endif
  endif
  
  ! Compute current time.
  timelocal = (it-1)*deltat + scheme_C(i_stage)*deltat
  
  if(myrank == 0 .and. mod(it, 100)==0) then
    WRITE(*,*) "****************************************************************"
    WRITE(*,"(a,i9,a,i1,a,e23.16,a,i3,a)") "Iteration", it, ", stage ", i_stage, ", local time", timelocal, &
    ". informations for process number ", myrank, "."
  endif
  
#ifdef USE_MPI
  if (NPROC > 1 .and. ninterface_acoustic > 0) then
    call assemble_MPI_vector_DG(rho_DG, buffer_DG_rho_P)
    call assemble_MPI_vector_DG(rhovx_DG, buffer_DG_rhovx_P)
    call assemble_MPI_vector_DG(rhovz_DG, buffer_DG_rhovz_P)
    call assemble_MPI_vector_DG(E_DG, buffer_DG_E_P)
    if(IONOSPHERIC_COUPLING) call assemble_MPI_vector_DG(Ni_DG, buffer_DG_Ni_P)
  endif
#endif

  ! Local Discontinuous Galerkin for viscous fluxes.
  if((maxval(muext) > ZEROl .or. maxval(etaext) > ZEROl .or. maxval(kappa_DG) > ZEROl) .or. CONSTRAin_HYDROSTATIC) then
#ifdef USE_MPI
    if (NPROC > 1 .and. ninterface_acoustic > 0) then
      call assemble_MPI_vector_DG(T_DG(1, :), buffer_DG_Tx_P)
      call assemble_MPI_vector_DG(T_DG(2, :), buffer_DG_Tz_P)
      call assemble_MPI_vector_DG(V_DG(1, 1, :), buffer_DG_Vxx_P)
      call assemble_MPI_vector_DG(V_DG(2, 2, :), buffer_DG_Vzz_P)
      call assemble_MPI_vector_DG(V_DG(1, 2, :), buffer_DG_Vxz_P)
      call assemble_MPI_vector_DG(V_DG(2, 1, :), buffer_DG_Vzx_P)
      call assemble_MPI_vector_DG(drho_DG(:), buffer_DG_drho_P)
      call assemble_MPI_vector_DG(dEp_DG(:), buffer_DG_dEp_P)
    endif
#endif
    call compute_viscous_tensors(T_DG, V_DG, drho_DG, dEp_DG, rho_DG, rhovx_DG, rhovz_DG, E_DG, timelocal)
    ! TODO: If CONSTRAin_HYDROSTATIC=.true. and muext=etaext=kappa_DG=0, only T_DG is needed. Consider computing only T_DG in that case.
  endif
  
  call compute_forces_acoustic_DG(rho_DG, rhovx_DG, rhovz_DG, E_DG, &
                                  T_DG, V_DG, drho_DG, dEp_DG, e1_DG, Ni_DG, &
                                  dot_rho, dot_rhovx, dot_rhovz, dot_E, dot_e1, dot_Ni, &
                                  timelocal)
  
  if (time_stepping_scheme == 3 .or. time_stepping_scheme == 4) then
    ! inverse mass matrix multiplication, in order to obtain actual RHS.
    dot_rho(:)   = dot_rho(:)   * rmass_inverse_acoustic_DG(:)
    dot_rhovx(:) = dot_rhovx(:) * rmass_inverse_acoustic_DG(:)
    dot_rhovz(:) = dot_rhovz(:) * rmass_inverse_acoustic_DG(:)
    dot_E(:)     = dot_E(:)     * rmass_inverse_acoustic_DG(:)
    dot_e1(:)    = dot_e1(:)    * rmass_inverse_acoustic_DG(:)
    if(IONOSPHERIC_COUPLING) then
      dot_Ni(:)    = dot_Ni(:)    * rmass_inverse_acoustic_DG(:)
    endif
    
    ! Time scheme low-storage update.
    ! See, e.g., Eq. (2) of J. Berland, C. Bogey, and C. Bailly, “Low-dissipation and low-dispersion fourth-order Runge-Kutta algorithm,” Comput. Fluids, vol. 35, no. 10, pp. 1459–1463, 2006.
    resu_rho   = scheme_A(i_stage)*resu_rho   + deltat*dot_rho
    resu_rhovx = scheme_A(i_stage)*resu_rhovx + deltat*dot_rhovx
    resu_rhovz = scheme_A(i_stage)*resu_rhovz + deltat*dot_rhovz
    resu_E     = scheme_A(i_stage)*resu_E     + deltat*dot_E
    resu_e1    = scheme_A(i_stage)*resu_e1    + deltat*dot_e1
    rho_DG     = rho_DG   + scheme_B(i_stage)*resu_rho
    rhovx_DG   = rhovx_DG + scheme_B(i_stage)*resu_rhovx
    rhovz_DG   = rhovz_DG + scheme_B(i_stage)*resu_rhovz
    E_DG       = E_DG     + scheme_B(i_stage)*resu_E
    e1_DG      = e1_DG    + scheme_B(i_stage)*resu_e1
    if(IONOSPHERIC_COUPLING) then
      Ni_DG    = Ni_DG    + scheme_B(i_stage)*resu_Ni
    endif

    ! If we want to compute kernels, we save regularly.
    if(.false.) then
      call save_forward_solution()
    endif
    
    ! Eventually check non-positivity.
    if(CHECK_NONPOSITIVITY) then
      if(     CHECK_NONPOSITIVITY_ON_ALL_PROCS &
         .or. ((.not. CHECK_NONPOSITIVITY_ON_ALL_PROCS) .and. myrank==0) &
        ) then
        ! If:    we check on all procs,
        !     or we check only on proc 0 and we are on proc 0.
        if(minval(rho_DG) < 1d-14) then
          WRITE(*,*) "***************************************************************"
          WRITE(*,*) "* CAREFUL, VERY SMALL DENSITY: ", minval(rho_DG), "    *"
          if(CHECK_NONPOSITIVITY_FinD_POinT) then
            ! Find where density is low.
            do ispec = 1,nspec
              do j = 1,NGLLZ
                do i = 1,NGLLX
                  if(rho_DG(ibool_DG(i, j, ispec))==minval(rho_DG)) then
                    WRITE(*, *) "* Element", ispec, ", GLL", i, j, ".         *"
                    write(*, *) "* Coords", coord(1, ibool_before_perio(i, j, ispec)), coord(2, ibool_before_perio(i, j, ispec)), &
                                ".*"
                  endif
                enddo
              enddo
            enddo
          endif ! Endif on CHECK_NONPOSITIVITY_FinD_POinT.
          WRITE(*,*) "***************************************************************"
        endif ! Endif on minval(rho_DG).
      endif ! Endif on CHECK_NONPOSITIVITY_ON_ALL_PROCS.
    endif ! Endif on CHECK_NONPOSITIVITY.
  else
    if(myrank==0) then
      ! Safeguard only, as normally the error prompted before (see above).
      write(*,*) "********************************"
      write(*,*) "*            ERROR             *"
      write(*,*) "********************************"
      write(*,*) "* This time_stepping_scheme is *"
      write(*,*) "* not implemented for DG yet.  *"
      write(*,*) "* time_stepping_scheme ", time_stepping_scheme
      write(*,*) "********************************"
      stop
    endif
  endif
  
  ! --------------------------- !
  ! Remove high-order           !
  ! coefficients of the         !
  ! solution.                   !
  ! --------------------------- !
  if(USE_SLOPE_LIMITER) then
    ! rho.
    if(CONSTRAin_HYDROSTATIC) then
      veloc_x = rho_DG - rho_init
      call SlopeLimit1(veloc_x, timelocal, 1)
      rho_DG = veloc_x + rho_init
    else
      call SlopeLimit1(rho_DG, timelocal, 1)
    endif
    ! rhovx.
    if(CONSTRAin_HYDROSTATIC) then
      veloc_x = rhovx_DG - rhovx_init
      call SlopeLimit1(veloc_x, timelocal, 1)
      rhovx_DG = veloc_x + rhovx_init
    else
      call SlopeLimit1(rhovx_DG, timelocal, 2)
    endif
    ! rhovz.
    if(CONSTRAin_HYDROSTATIC) then
      veloc_x = rhovz_DG - rhovz_init
      call SlopeLimit1(veloc_x, timelocal, 1)
      rhovz_DG = veloc_x + rhovz_init
    else
      call SlopeLimit1(rhovz_DG, timelocal, 3)
    endif
    ! E.
    if(CONSTRAin_HYDROSTATIC) then
      veloc_x = E_DG - E_init
      call SlopeLimit1(veloc_x, timelocal, 1)
      E_DG = veloc_x + E_init
    else
      call SlopeLimit1(E_DG, timelocal, 4)
    endif
  endif
end subroutine compute_forces_acoustic_DG_main


#if 0
! ------------------------------------------------------------ !
! save_forward_solution                                        !
! ------------------------------------------------------------ !
! Save some DG fields at some iterations to file, to be used in conjunction with the 'read_forward_solution' routine.
! Prefer using wavefield dumps.
! We keep this code here in case it is needed.

subroutine save_forward_solution()
  use constants, only: CUSTOM_REAL
  use specfem_par
  use mpi
  implicit none
  include "precision.h"
  ! Input/Output.
  ! N. A.
  ! Local Variables.
  character(len=MAX_STRinG_LEN) :: outputname
  integer :: ier
  integer :: IT_SAVE
  
  ! Save last frame.
  IT_SAVE = 100
  if (SAVE_FORWARD .and. SIMULATION_TYPE == 1 .and. any_acoustic &
        .and. mod(it, IT_SAVE) == 0 .and. i_stage == 1) then
    if (myrank == 0) then
      write(IMAin,*)
      write(IMAin,*) 'Saving DG acoustic frame...'
      write(IMAin,*)
      call flush_IMAin()
    endif
    write(outputname,'(a,i6.6,a,i6.6,a)') 'frame_DG_',myrank,'_',it,'.bin'
    open(unit=55,file='outPUT_FILES/output_frames/'//trim(outputname),status='unknown',form='unformatted',iostat=ier)
    if (ier /= 0) call exit_MPI(myrank,'Error opening file frame_DG...')
    write(55) rho_DG
    write(55) rhovx_DG
    write(55) rhovz_DG
    write(55) E_DG
    close(55)
  endif
end subroutine save_forward_solution


! ------------------------------------------------------------ !
! read_forward_solution                                        !
! ------------------------------------------------------------ !
! Read some DG fields from file, to be used in conjunction with the 'save_forward_solution' routine.

subroutine read_forward_solution(it_temp)
  use constants, only: CUSTOM_REAL
  use specfem_par
  use mpi
  implicit none
  include "precision.h"
  ! Input/Output.
  ! N. A.
  ! Local Variables.
  character(len=MAX_STRinG_LEN) :: outputname
  integer :: ier
  integer :: it_temp, IT_SAVE
  
  IT_SAVE = 500
    
  ! Save last frame.
  if (SAVE_FORWARD .and. SIMULATION_TYPE == 1 .and. any_acoustic &
        .and. mod(it, IT_SAVE) == 0 .and. i_stage == 1) then
    if (myrank == 0) then
      write(IMAin,*)
      write(IMAin,*) 'Saving DG acoustic frame...'
      write(IMAin,*)
      call flush_IMAin()
    endif
    write(outputname,'(a,i6.6,a,i6.6,a)') 'frame_DG_',myrank,'_',it_temp,'.bin'
    open(unit=55,file='outPUT_FILES/'//trim(outputname),status='unknown',form='unformatted',iostat=ier)
    if (ier /= 0) call exit_MPI(myrank,'Error opening file frame_DG...')
    read(55) rho_DG
    read(55) rhovx_DG
    read(55) rhovz_DG
    read(55) E_DG
    close(55)
  endif
end subroutine read_forward_solution
#endif


! ------------------------------------------------------------ !
! prepare_MPI_DG                                               !
! ------------------------------------------------------------ !
! Prepare MPI-related exchange buffers and links between neighbouring elements across CPUs.

subroutine prepare_MPI_DG()
  use constants, only: CUSTOM_REAL,NGLLX,NGLLZ,TinYVAL
  use specfem_par, only: MPI_transfer_iface, is_corner, link_ijispec_iface, &
                         max_interface_size, ninterface, ninterface_acoustic_dg, inum_interfaces_acoustic_dg, &
                         link_iface_ijispec, my_neighbours, nibool_interfaces_acoustic_dg, ibool_interfaces_acoustic_dg, &
                         tab_requests_send_recv_dg, neighbor_dg_iface, &
                         nspec, ibool, ibool_dg, nglob_DG, coord, &
                         myrank
  use mpi

  implicit none
  include "precision.h"
  
  ! Input/Output.
  ! N. A.
  
  ! Local Variables.
  real(kind=CUSTOM_REAL), parameter :: ONEl = 1._CUSTOM_REAL
  integer :: iinterface, ipoin, num_interface, iglob, i, j, ispec, iglob_DG, nb_values, ier
  integer, dimension(nglob_DG, 3) :: link_ij_iglob
  double precision, dimension(NGLLX*max_interface_size,ninterface) :: &
      buffer_recv_faces_vector_DG_i, &
      buffer_send_faces_vector_DG_i, &
      buffer_recv_faces_vector_DG_j, &
      buffer_send_faces_vector_DG_j, &
      buffer_send_faces_vector_DG_i1try, &
      buffer_recv_faces_vector_DG_i1try, &
      buffer_send_faces_vector_DG_j1try, &
      buffer_recv_faces_vector_DG_j1try, &
      buffer_send_faces_vector_DG_i2try, &
      buffer_recv_faces_vector_DG_i2try, &
      buffer_send_faces_vector_DG_j2try, &
      buffer_recv_faces_vector_DG_j2try
  integer :: iface1, iface, neighbor, neighbor_corner, iface_corner, iface1_corner
  double precision :: coord_i_1, coord_j_1, coord_i_2, coord_j_2, &
      coord_i_21_try, coord_j_21_try, coord_i_22_try, coord_j_22_try
  integer :: i_try, j_try
  logical :: one_other_node_is_found, one_other_node_is_found_corner
  
  do ispec = 1, nspec
    do i = 1, NGLLX
      do j = 1, NGLLZ
        iglob_DG = ibool_DG(i, j, ispec)
        link_ij_iglob(iglob_DG,1) = i
        link_ij_iglob(iglob_DG,2) = j
        link_ij_iglob(iglob_DG,3) = ispec
        if(     (i == 1 .and. (j == 1 .or. j == NGLLZ)) &
           .or. (i == NGLLX .and. (j == 1 .or. j == NGLLZ))) then
          is_corner(i, j) = .true.
        endif
      enddo
    enddo
  enddo
  
  buffer_send_faces_vector_DG_i = -ONEl
  buffer_send_faces_vector_DG_j = -ONEl
  buffer_send_faces_vector_DG_i2try = -ONEl
  buffer_send_faces_vector_DG_j2try = -ONEl
  buffer_send_faces_vector_DG_i1try = -ONEl
  buffer_send_faces_vector_DG_j1try = -ONEl
  buffer_recv_faces_vector_DG_i = -ONEl
  buffer_recv_faces_vector_DG_j = -ONEl
  buffer_recv_faces_vector_DG_i2try = -ONEl
  buffer_recv_faces_vector_DG_j2try = -ONEl
  buffer_recv_faces_vector_DG_i1try = -ONEl
  buffer_recv_faces_vector_DG_j1try = -ONEl
  
  ! MPI SEND inFO ABout DIAG ELEMENT OR NOT
  do iinterface = 1, ninterface_acoustic_DG
    ! gets interface index in the range of all interfaces [1,ninterface]
    num_interface = inum_interfaces_acoustic_DG(iinterface)
    
    ! loops over all interface points
    do ipoin = 1, nibool_interfaces_acoustic_DG(num_interface)
      iglob = ibool_interfaces_acoustic_DG(ipoin,num_interface)
      i = link_ij_iglob(iglob,1)
      j = link_ij_iglob(iglob,2)
      ispec = link_ij_iglob(iglob,3)
      buffer_send_faces_vector_DG_i(ipoin,iinterface) = coord(1,ibool(i,j,ispec))
      buffer_send_faces_vector_DG_j(ipoin,iinterface) = coord(2,ibool(i,j,ispec))
      
      ! RECOVER CURRENT ELEMENT PROPERTIES (IFACE1, IFACE)
      ! QUENTIN BRISSAUD: Here, we store one by one the coordinates of the interface nodes in the current MPI partition in table: buffer_send_faces_vector_DG_i, buffer_send_faces_vector_DG_j and to be sent to the neighbors
      ! QUENTIN BRISSAUD: in addition, we store the coordinates another point (iface1,iface) belonging to the same element face than the current node (i,j,ispec) to be able to discriminate later on and after the I_SEND on the neighbor's parition to which face the current point belongs to.
      ! QUENTIN BRISSAUD: Note that if it's a corner, there are two element faces that come into play. We thus store an additional point (iface1_corner,iface_corner) corresponding to the other face
      iface1 = link_ijispec_iface(i,j,ispec,1,1)
      iface1 = iface1 + 1
      if(iface1 > 5) then
        iface1 = link_ijispec_iface(i,j,ispec,1,1) - 1
      endif
      iface = link_ijispec_iface(i,j,ispec,2,1)
      iface1_corner = -1
      iface_corner = -1
      neighbor_corner = -1
      if(is_corner(i,j) .and. link_ijispec_iface(i,j,ispec,1,2) > -1) then
        iface1_corner = link_ijispec_iface(i,j,ispec,1,2)
        iface1_corner = iface1_corner + 1
        if(iface1_corner > 5) then
          iface1_corner = link_ijispec_iface(i,j,ispec,1,2) - 1
        endif
        iface_corner = link_ijispec_iface(i,j,ispec,2,2)
        i = link_iface_ijispec(iface1_corner, iface_corner, ispec,1)
        j = link_iface_ijispec(iface1_corner, iface_corner, ispec,2)
        buffer_send_faces_vector_DG_i2try(ipoin,iinterface) = coord(1,ibool(i,j,ispec))
        buffer_send_faces_vector_DG_j2try(ipoin,iinterface) = coord(2,ibool(i,j,ispec))
      else
        buffer_send_faces_vector_DG_i2try(ipoin,iinterface) = -ONEl
        buffer_send_faces_vector_DG_j2try(ipoin,iinterface) = -ONEl
      endif
      i = link_iface_ijispec(iface1, iface, ispec,1)
      j = link_iface_ijispec(iface1, iface, ispec,2)
      buffer_send_faces_vector_DG_i1try(ipoin,iinterface) = coord(1,ibool(i,j,ispec))
      buffer_send_faces_vector_DG_j1try(ipoin,iinterface) = coord(2,ibool(i,j,ispec))
    enddo ! Enddo on ipoin.
  enddo ! Enddo on iinterface.
  
  do iinterface = 1, ninterface_acoustic_DG
    ! gets interface index in the range of all interfaces [1,ninterface]
    num_interface = inum_interfaces_acoustic_DG(iinterface)
    nb_values = nibool_interfaces_acoustic_DG(num_interface)
    
    call MPI_ISEND( buffer_send_faces_vector_DG_i(1,iinterface), &
             nb_values, MPI_DOUBLE, &
             my_neighbours(num_interface), 14, MPI_COMM_WORLD, &
             tab_requests_send_recv_DG(iinterface), ier)
    if (ier /= MPI_SUCCESS) then
      call exit_MPI(myrank,'MPI_ISEND unsuccessful in prepare_MPI_DG - i')
    endif

    ! starts a non-blocking receive
    call MPI_IRECV ( buffer_recv_faces_vector_DG_i(1,iinterface), &
             nb_values, MPI_DOUBLE, &
             my_neighbours(num_interface), 14, MPI_COMM_WORLD, &
             tab_requests_send_recv_DG(ninterface_acoustic_DG+iinterface), ier)
    if (ier /= MPI_SUCCESS) then
      call exit_MPI(myrank,'MPI_IRECV unsuccessful in prepare_MPI_DG - i')
    endif
    
    call MPI_ISEND( buffer_send_faces_vector_DG_j(1,iinterface), &
             nb_values, MPI_DOUBLE, &
             my_neighbours(num_interface), 14, MPI_COMM_WORLD, &
             tab_requests_send_recv_DG(ninterface_acoustic_DG*2+iinterface), ier)
    if (ier /= MPI_SUCCESS) then
      call exit_MPI(myrank,'MPI_ISEND unsuccessful in prepare_MPI_DG - j')
    endif

    ! starts a non-blocking receive
    call MPI_IRECV ( buffer_recv_faces_vector_DG_j(1,iinterface), &
             nb_values, MPI_DOUBLE, &
             my_neighbours(num_interface), 14, MPI_COMM_WORLD, &
             tab_requests_send_recv_DG(ninterface_acoustic_DG*3+iinterface), ier)
    if (ier /= MPI_SUCCESS) then
      call exit_MPI(myrank,'MPI_IRECV unsuccessful in prepare_MPI_DG - j')
    endif
             
    call MPI_ISEND( buffer_send_faces_vector_DG_i1try(1,iinterface), &
             nb_values, MPI_DOUBLE, &
             my_neighbours(num_interface), 14, MPI_COMM_WORLD, &
             tab_requests_send_recv_DG(ninterface_acoustic_DG*4+iinterface), ier)
    if (ier /= MPI_SUCCESS) then
      call exit_MPI(myrank,'MPI_ISEND unsuccessful in prepare_MPI_DG - i1try')
    endif

    ! starts a non-blocking receive
    call MPI_IRECV ( buffer_recv_faces_vector_DG_i1try(1,iinterface), &
             nb_values, MPI_DOUBLE, &
             my_neighbours(num_interface), 14, MPI_COMM_WORLD, &
             tab_requests_send_recv_DG(ninterface_acoustic_DG*5+iinterface), ier)
    if (ier /= MPI_SUCCESS) then
      call exit_MPI(myrank,'MPI_IRECV unsuccessful in prepare_MPI_DG - i1try')
    endif
    
    call MPI_ISEND( buffer_send_faces_vector_DG_j1try(1,iinterface), &
             nb_values, MPI_DOUBLE, &
             my_neighbours(num_interface), 14, MPI_COMM_WORLD, &
             tab_requests_send_recv_DG(ninterface_acoustic_DG*6+iinterface), ier)
    if (ier /= MPI_SUCCESS) then
      call exit_MPI(myrank,'MPI_ISEND unsuccessful in prepare_MPI_DG - j1try')
    endif

    ! starts a non-blocking receive
    call MPI_IRECV ( buffer_recv_faces_vector_DG_j1try(1,iinterface), &
             nb_values, MPI_DOUBLE, &
             my_neighbours(num_interface), 14, MPI_COMM_WORLD, &
             tab_requests_send_recv_DG(ninterface_acoustic_DG*7+iinterface), ier)
    if (ier /= MPI_SUCCESS) then
      call exit_MPI(myrank,'MPI_IRECV unsuccessful in prepare_MPI_DG - j1try')
    endif
    
    call MPI_ISEND( buffer_send_faces_vector_DG_i2try(1,iinterface), &
             nb_values, MPI_DOUBLE, &
             my_neighbours(num_interface), 14, MPI_COMM_WORLD, &
             tab_requests_send_recv_DG(ninterface_acoustic_DG*8+iinterface), ier)
    if (ier /= MPI_SUCCESS) then
      call exit_MPI(myrank,'MPI_ISEND unsuccessful in prepare_MPI_DG - i2try')
    endif

    ! starts a non-blocking receive
    call MPI_IRECV ( buffer_recv_faces_vector_DG_i2try(1,iinterface), &
             nb_values, MPI_DOUBLE, &
             my_neighbours(num_interface), 14, MPI_COMM_WORLD, &
             tab_requests_send_recv_DG(ninterface_acoustic_DG*9+iinterface), ier)
    if (ier /= MPI_SUCCESS) then
      call exit_MPI(myrank,'MPI_IRECV unsuccessful in prepare_MPI_DG - i2try')
    endif
    
    call MPI_ISEND( buffer_send_faces_vector_DG_j2try(1,iinterface), &
             nb_values, MPI_DOUBLE, &
             my_neighbours(num_interface), 14, MPI_COMM_WORLD, &
             tab_requests_send_recv_DG(ninterface_acoustic_DG*10+iinterface), ier)
    if (ier /= MPI_SUCCESS) then
      call exit_MPI(myrank,'MPI_ISEND unsuccessful in prepare_MPI_DG - j2try')
    endif

    ! starts a non-blocking receive
    call MPI_IRECV ( buffer_recv_faces_vector_DG_j2try(1,iinterface), &
             nb_values, MPI_DOUBLE, &
             my_neighbours(num_interface), 14, MPI_COMM_WORLD, &
             tab_requests_send_recv_DG(ninterface_acoustic_DG*11+iinterface), ier)
    if (ier /= MPI_SUCCESS) then
      call exit_MPI(myrank,'MPI_IRECV unsuccessful in prepare_MPI_DG - j2try')
    endif
  enddo
  
  ! waits for MPI requests to complete (recv)
  ! each wait returns once the specified MPI request completed
  do iinterface = 1, 12*ninterface_acoustic_DG
    call MPI_Wait(tab_requests_send_recv_DG(iinterface), MPI_STATUS_IGNORE, ier)
  enddo
  
  MPI_transfer_iface = -1 ! By default, specify that for the triplet (iface1,iface,ispec), the values should not be sought in another partition.
  
  ! Loop on interfaces between partitions to find correspondances of faces between partitions.
  do iinterface = 1, ninterface_acoustic_DG
    ! gets interface index in the range of all interfaces [1,ninterface]
    num_interface = inum_interfaces_acoustic_DG(iinterface)
    ! loops over all interface points
    do ipoin = 1, nibool_interfaces_acoustic_DG(num_interface)
      iglob = ibool_interfaces_acoustic_DG(ipoin,num_interface)
      
      i = link_ij_iglob(iglob,1)
      j = link_ij_iglob(iglob,2)
      ispec = link_ij_iglob(iglob,3)
      
      coord_i_2      = buffer_recv_faces_vector_DG_i(ipoin, iinterface)
      coord_j_2      = buffer_recv_faces_vector_DG_j(ipoin, iinterface)
      coord_i_21_try = buffer_recv_faces_vector_DG_i1try(ipoin, iinterface)
      coord_j_21_try = buffer_recv_faces_vector_DG_j1try(ipoin, iinterface)
      coord_i_22_try = buffer_recv_faces_vector_DG_i2try(ipoin, iinterface)
      coord_j_22_try = buffer_recv_faces_vector_DG_j2try(ipoin, iinterface)
      
      ! Define variables "neighbor" and "neighbor_corner".
      iface  = link_ijispec_iface(i,j,ispec,2,1) ! Face number on which the point (i,j,ispec) is.
      iface1 = link_ijispec_iface(i,j,ispec,1,1) ! GLL point on current face (iface), corresponding to the (i,j,ispec) point is.
      neighbor = neighbor_DG_iface(iface1,iface,ispec,3) ! To check whether or not the considered point has a neighbour in current CPU or in another CPU.
      iface1_corner = -1
      iface_corner  = -1
      neighbor_corner = -1
      if(is_corner(i, j)) then
        iface1_corner   = link_ijispec_iface(i,j,ispec,1,2)
        iface_corner    = link_ijispec_iface(i,j,ispec,2,2)
        neighbor_corner = neighbor_DG_iface(iface1_corner,iface_corner,ispec,3)
      endif
              
      ! QUENTIN BRISSAUD: Here we try to check to which face (on the neighbor's parition) the current point (i,j,ispec) belongs to
      ! QUENTIN BRISSAUD: If one of the point on the current face matches the coordinates of the coordinates received in buffer_recv_faces_vector_DG_i1try, buffer_recv_faces_vector_DG_j1try, we found it => one_other_node_is_found = .true.
      ! QUENTIN BRISSAUD: If it's a corner we need to check an extra face which coordinates are given by coord_i_22_try, coord_j_22_try. If there is a match => one_other_node_is_found_corner = .true.
      one_other_node_is_found = .false.
      one_other_node_is_found_corner = .false.
      do iface1 = 1,NGLLX
        i_try = link_iface_ijispec(iface1, iface, ispec, 1)
        j_try = link_iface_ijispec(iface1, iface, ispec, 2)
        
        coord_i_1 = coord(1,ibool(i_try,j_try,ispec))
        coord_j_1 = coord(2,ibool(i_try,j_try,ispec))
        
        if(     (abs(coord_i_1-coord_i_21_try)<TinYVAL .and. abs(coord_j_1-coord_j_21_try)<TinYVAL) &
           .or. (abs(coord_i_1-coord_i_22_try)<TinYVAL .and. abs(coord_j_1-coord_j_22_try)<TinYVAL)) then
          one_other_node_is_found = .true.
        endif
                
        if(is_corner(i, j)) then
          i_try = link_iface_ijispec(iface1, iface_corner, ispec, 1)
          j_try = link_iface_ijispec(iface1, iface_corner, ispec, 2)
          coord_i_1 = coord(1,ibool(i_try,j_try,ispec))
          coord_j_1 = coord(2,ibool(i_try,j_try,ispec))
          if(     (abs(coord_i_1-coord_i_21_try)<TinYVAL .and. abs(coord_j_1-coord_j_21_try)<TinYVAL) &
             .or. (abs(coord_i_1-coord_i_22_try)<TinYVAL .and. abs(coord_j_1-coord_j_22_try)<TinYVAL)) then
            one_other_node_is_found_corner = .true.
          endif
        endif
      enddo ! Enddo on iface1.
      
      iface1 = link_ijispec_iface(i,j,ispec,1,1)
      coord_i_1 = coord(1,ibool(i,j,ispec))
      coord_j_1 = coord(2,ibool(i,j,ispec))
      
      ! [QUENTIN BRISSAUD]
      !   If:
      !           the current point has no DG neighbors (i.e. no element neighbor in the same parititon)
      !     .and. we are on the right element face (i.e. one_other_node_is_found .or. one_other_node_is_found_corner)
      !     .and. the current point's coordinates are the same than the ones received in the I_RECV from the MPI neighbor
      !   then:
      !     we store this information (link between (iface1, iface, ispec) abd (ipoin, num_interface)) in table MPI_transfer_iface.
      if(      (neighbor == -1 .or. neighbor_corner == -1) &
         .and. (abs(coord_i_1-coord_i_2)<TinYVAL .and. abs(coord_j_1-coord_j_2)<TinYVAL) &
         .and. (one_other_node_is_found .or. one_other_node_is_found_corner)) then
        if(one_other_node_is_found) then
          MPI_transfer_iface(iface1,iface,ispec,1) = ipoin ! Specifies that point ispec on face (iface1,iface) is point ipoin (on interface num_interface).
          MPI_transfer_iface(iface1,iface,ispec,2) = num_interface ! Specifies that point ispec on face (iface1,iface) is (point ipoin) on interface num_interface.
        endif
        if(one_other_node_is_found_corner) then
          ! Same as above, but for corner point.
          MPI_transfer_iface(iface1_corner,iface_corner,ispec,1) = ipoin 
          MPI_transfer_iface(iface1_corner,iface_corner,ispec,2) = num_interface
        endif
      endif
    enddo ! Enddo on ipoin.
  enddo ! Enddo on iinterface.
end subroutine prepare_MPI_DG


! ------------------------------------------------------------ !
! setUpVandermonde                                             !
! ------------------------------------------------------------ !
! Setup Vandermonde matrices. Only used when the slope limiter is.

subroutine setUpVandermonde()
  use constants, only: CUSTOM_REAL,NGLLX,NGLLZ
  use specfem_par, only: Vandermonde, invVandermonde, Drx, Drz
  implicit none
  ! Input/Output.
  ! N. A.
  ! Local Variables.
  ! N. A.
  call compute_Vander_matrices(NGLLX*NGLLZ, Vandermonde, invVandermonde, Drx, Drz )
end subroutine setUpVandermonde


! ------------------------------------------------------------ !
! compute_Vander_matrices                                      !
! ------------------------------------------------------------ !
! Compute Vandermonde matrices. Only used when the slope limiter is.

subroutine compute_Vander_matrices(np, V1D, V1D_inv, Drx, Drz)
  use constants, only: CUSTOM_REAL,NGLLX,NGLLZ,GAUSSALPHA,GAUSSBETA
  use specfem_par, only: xigll, zigll
  
  implicit none
  
  ! Input/Output.
  integer np
  real(kind=CUSTOM_REAL) V1D(np,np), V1D_inv(np, np), Drx(np, np), Drz(np, np)
  
  ! Local Variables.
  real(kind=CUSTOM_REAL), parameter :: ZEROl = 0._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: threshold = 0.00001_CUSTOM_REAL
  double precision, dimension(np, np) :: V1D_d, V1D_inv_d
  double precision :: p, pd, pm1, pdm1, pm2, pdm2, p2, pd2
  integer i, j, k, errorflag, l, m, n
  
  ! Initialise.
  V1D     = ZEROl
  V1D_inv = ZEROl
  Drx     = ZEROl
  Drz     = ZEROl
  V1D_d   = ZEROl
  
  ! New matrices.
  k = 0
  do m = 1,NGLLX
    do n = 1,NGLLZ
      k = k+1
      l = 0
      do i = 1,NGLLX
        do j = 1,NGLLZ
          l = l + 1
          call jacobf(p,pd,pm1,pdm1,pm2,pdm2,i-1,GAUSSALPHA,GAUSSBETA,xigll(m))
          call jacobf(p2,pd2,pm1,pdm1,pm2,pdm2,j-1,GAUSSALPHA,GAUSSBETA,zigll(n))
          V1D(k,l) = REAL(p*p2, kind=CUSTOM_REAL)
          if(abs(V1D(k,l)) < threshold) then
            V1D(k,l) = ZEROl
          endif
          V1D_d(k,l) = p*p2
        enddo
      enddo
    enddo
  enddo
  call findInv(V1D_d, V1D_inv_d, np, errorflag)
  do j=1,np
    do i=1,np
        V1D_inv(i,j) = REAL(V1D_inv_d(i,j), kind=CUSTOM_REAL)
        if(abs(V1D_inv(i,j)) < threshold) then
          V1D_inv(i,j) = ZEROl
        endif
    enddo
  enddo
end subroutine compute_Vander_matrices


! ------------------------------------------------------------ !
! findInv                                                      !
! ------------------------------------------------------------ !
! Subroutine to find the inverse of a square matrix.
! Author : Louisda16th a.k.a Ashwith J. Rego
! Reference : Algorithm has been well explained at http://math.uww.edu/~mcfarlat/inverse.htm or http://www.tutor.ms.unimelb.edu.au/matrix/matrix_inverse.html.
! TODO: Replace the inner code of this routine (or better, directly the calls to it) by a LAPACK call.

subroutine findInv(matrix, inverse, n, errorflag)
  use constants, only: CUSTOM_REAL

  implicit none
  
  ! Input/Output.
  integer, intent(in) :: n
  integer, intent(out) :: errorflag !Return error status. -1 for error, 0 for normal
  double precision, intent(in), dimension(n,n) :: matrix !input matrix
  double precision, intent(out), dimension(n,n) :: inverse !inverted matrix

  ! Local Variables.
  logical :: FLAG = .TRUE.
  integer :: i, j, k
  double precision :: m
  double precision, dimension(n, 2*n) :: augmatrix

  ! Augment input matrix with an identity matrix.
  DO i = 1, n
    DO j = 1, 2*n
      IF (j <= n ) THEN
        augmatrix(i,j) = matrix(i,j)
      ELSE IF ((i+n) == j) THEN
        augmatrix(i,j) = 1
      Else
        augmatrix(i,j) = 0
      endif
    enddo
  enddo

  ! Reduce augmented matrix to upper traingular form.
  DO k =1, n-1
    IF (augmatrix(k,k) == 0) THEN
      FLAG = .FALSE.
      DO i = k+1, n
        IF (augmatrix(i,k) /= 0) THEN
          DO j = 1,2*n
            augmatrix(k,j) = augmatrix(k,j)+augmatrix(i,j)
          enddo
          FLAG = .TRUE.
          EXIT
        endif
        IF (FLAG .EQV. .FALSE.) THEN
          print*, "Matrix is non - invertible"
          inverse = 0
          errorflag = -1
          return
        endif
      enddo
    endif

    DO j = k+1, n 
      m = augmatrix(j,k)/augmatrix(k,k)
      DO i = k, 2*n
        augmatrix(j,i) = augmatrix(j,i) - m*augmatrix(k,i)
      enddo
    enddo
  enddo

  ! Test for invertibility.
  DO i = 1, n
    IF (augmatrix(i,i) == 0) THEN
      print*, "Matrix is non - invertible"
      inverse = 0
      errorflag = -1
      return
    endif
  enddo

  ! Make diagonal elements as 1.
  DO i = 1 , n
    m = augmatrix(i,i)
    DO j = i , (2 * n) 
      augmatrix(i,j) = (augmatrix(i,j) / m)
    enddo
  enddo

  ! Reduced right side half of augmented matrix to identity matrix.
  DO k = n-1, 1, -1
    DO i =1, k
      m = augmatrix(i,k+1)
      DO j = k, (2*n)
        augmatrix(i,j) = augmatrix(i,j) -augmatrix(k+1,j) * m
      enddo
    enddo
  enddo 

  ! Store answer.
  DO i =1, n
    DO j = 1, n
      inverse(i,j) = augmatrix(i,j+n)
    enddo
  enddo
  errorflag = 0
  return
end subroutine findInv

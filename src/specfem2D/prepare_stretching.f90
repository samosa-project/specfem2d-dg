! ------------------------------------------------------------ !
! prepare_stretching                                           !
! ------------------------------------------------------------ !
! Prepares the stretching factors (in stretching_ya) for the buffer-based real stretching absorbing boundary conditions.
! In particular:
! - Initialises the flag variable encoding where buffers are (stretching_buffer).
! - Sets the values of the stretching in all directions (via the 'stretching_function' routine).

subroutine prepare_stretching()
  use specfem_par, only: coord, ibool_before_perio, &
                         myrank, nglob,nspec, &
                         ispec_is_acoustic_DG, &
                         ABC_STRETCH_LEFT, ABC_STRETCH_RIGHT, ABC_STRETCH_TOP, ABC_STRETCH_BOTTOM, &
                         ABC_STRETCH_TOP_LBUF, ABC_STRETCH_LEFT_LBUF, ABC_STRETCH_BOTTOM_LBUF, ABC_STRETCH_RIGHT_LBUF, &
                         stretching_ya, any_elastic,&
                         mesh_xmin, mesh_xmax, mesh_zmin, mesh_zmax, stretching_buffer
  
  use constants, only: CUSTOM_REAL, NDIM, NGLLX, NGLLZ

  implicit none
  
  ! Input/Output.
  ! N. A.
  
  ! Local Variables.
  real(kind=CUSTOM_REAL), parameter :: ZERO = 0._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: ONE  = 1._CUSTOM_REAL
  integer iglob_unique, ispec, i, j
  real(kind=CUSTOM_REAL) :: x,z
  real(kind=CUSTOM_REAL) :: r_l ! Relative buffer coordinate (which is equal to 0 at beginning of buffer, and to 1 at the end).
  
  ! Error checking.
  if(myrank==0) then
    if(ABC_STRETCH_LEFT_LBUF+ABC_STRETCH_RIGHT_LBUF>=mesh_xmax-mesh_xmin) then
      write(*,*) "********************************"
      write(*,*) "*            ERROR             *"
      write(*,*) "********************************"
      write(*,*) "* Sum of left and right        *"
      write(*,*) "* stretching ABC buffer        *"
      write(*,*) "* lengths must be lesser than  *"
      write(*,*) "* the mesh's horizontal size:  *"
      write(*,*) "*  left buffer length :  ", ABC_STRETCH_LEFT_LBUF
      write(*,*) "*  right buffer length : ", ABC_STRETCH_RIGHT_LBUF
      write(*,*) "*  horizontal size:      ", mesh_xmax-mesh_xmin
      write(*,*) "********************************"
      call exit_MPI(myrank, " ")
    endif
    if(ABC_STRETCH_TOP_LBUF+ABC_STRETCH_BOTTOM_LBUF>=mesh_zmax-mesh_zmin) then
      write(*,*) "********************************"
      write(*,*) "*            ERROR             *"
      write(*,*) "********************************"
      write(*,*) "* Sum of left and right        *"
      write(*,*) "* stretching ABC buffer        *"
      write(*,*) "* lengths must be lesser than  *"
      write(*,*) "* the mesh's horizontal size:  *"
      write(*,*) "*  top buffer length :    ", ABC_STRETCH_TOP_LBUF
      write(*,*) "*  bottom buffer length : ", ABC_STRETCH_BOTTOM_LBUF
      write(*,*) "*  vertical size:         ", mesh_zmax-mesh_zmin
      write(*,*) "********************************"
      call exit_MPI(myrank, " ")
    endif
    ! TODO: Check that the user did not ask for a bottom buffer if they have an elastic medium under the acoustic medium.
    ! Quick hack, assuming the elastic media is always under the acoustic media.
    if(any_elastic .and. ABC_STRETCH_BOTTOM) then
      write(*,*) "********************************"
      write(*,*) "*            ERROR             *"
      write(*,*) "********************************"
      write(*,*) "* Cannot use stretching on the *"
      write(*,*) "* bottom boundary if there is  *"
      write(*,*) "* an elastic media.            *"
      write(*,*) "********************************"
      call exit_MPI(myrank, " ")
    endif
  endif
  
  ! Warnings and general user outputs.
  if(myrank==0) then
    if(any_elastic) then
      write(*,*) "********************************"
      write(*,*) "*           WARNING            *"
      write(*,*) "********************************"
      write(*,*) "* Stretching absorbing         *"
      write(*,*) "* boundary conditions are      *"
      write(*,*) "* activated and an elastic     *"
      write(*,*) "* media is present.            *"
      write(*,*) "* Implementation is            *"
      write(*,*) "* approximate, be careful.     *"
      write(*,*) "********************************"
    endif
    call flush_IMAIN()
  endif
  
  allocate(stretching_ya(NDIM, nglob)) ! Two-dimensionnal stretching. Change it to 3 for 3D.
  allocate(stretching_buffer(nglob)) ! Stretching buffers code: see specfem2D_par.
  stretching_ya(:, :) = ONE ! By default, mesh is not stretched.
  stretching_buffer(:) = 0 ! By default, mesh is not stretched.
  
  do ispec = 1, nspec
    if(ispec_is_acoustic_DG(ispec)) then
      do j = 1, NGLLZ
        do i = 1, NGLLX
          iglob_unique = ibool_before_perio(i, j, ispec)
          x=coord(1, iglob_unique)
          z=coord(2, iglob_unique)
          if(     (ABC_STRETCH_LEFT   .and. x < mesh_xmin + ABC_STRETCH_LEFT_LBUF) & ! left stretching and in left buffer zone
             .or. (ABC_STRETCH_RIGHT  .and. x > mesh_xmax - ABC_STRETCH_RIGHT_LBUF) & ! right stretching and in right buffer zone
             .or. (ABC_STRETCH_BOTTOM .and. z < mesh_zmin + ABC_STRETCH_BOTTOM_LBUF) & ! bottom stretching and in bottom buffer zone
             .or. (ABC_STRETCH_TOP    .and. z > mesh_zmax - ABC_STRETCH_TOP_LBUF)) then ! top stretching and in top buffer zone
            
            if(ABC_STRETCH_TOP) then
              r_l = (z - mesh_zmax)/ABC_STRETCH_TOP_LBUF + ONE
              if(r_l>ZERO .and. r_l<=ONE) then
                call stretching_function(r_l, stretching_ya(2, iglob_unique))
                stretching_buffer(iglob_unique) = ibset(stretching_buffer(iglob_unique), 0) ! Set 1st LSB to 1.
              endif
            endif
            if(ABC_STRETCH_LEFT) then
              r_l = ONE - (x - mesh_xmin)/ABC_STRETCH_LEFT_LBUF
              if(r_l>ZERO .and. r_l<=ONE) then
                call stretching_function(r_l, stretching_ya(1, iglob_unique))
                stretching_buffer(iglob_unique) = ibset(stretching_buffer(iglob_unique), 1) ! Set 2nd LSB to 1.
              endif
            endif
            if(ABC_STRETCH_BOTTOM) then
              r_l = ONE - (z - mesh_zmin)/ABC_STRETCH_BOTTOM_LBUF
              if(r_l>ZERO .and. r_l<=ONE) then
                call stretching_function(r_l, stretching_ya(2, iglob_unique))
                stretching_buffer(iglob_unique) = ibset(stretching_buffer(iglob_unique), 2) ! Set 3rd LSB to 1.
                !write(*,*) 'OMEGAKEK', stretching_buffer(ibool_before_perio(i,j,ispec))
              endif
            endif
            if(ABC_STRETCH_RIGHT) then
              r_l = (x - mesh_xmax)/ABC_STRETCH_RIGHT_LBUF + ONE
              if(r_l>ZERO .and. r_l<=ONE) then
                call stretching_function(r_l, stretching_ya(1, iglob_unique))
                stretching_buffer(iglob_unique) = ibset(stretching_buffer(iglob_unique), 3) ! Set 4th LSB to 1.
              endif
            endif
          endif ! Endif.
        enddo ! Enddo on i.
      enddo ! Enddo on j.
    endif ! Endif on ispec_is_acoustic_DG(ispec).
  enddo ! Enddo on ispec.
  call synchronize_all()
end subroutine prepare_stretching


! ------------------------------------------------------------ !
! stretching_function                                          !
! ------------------------------------------------------------ !
! Implementation of ya(x). For now, no distinction has to be made in each direction.
! New expressions can be implemented here. Some compatibility conditions should be respected. The function has be 1 when r_l==0. The derivative of the function should be 0 when r_l==0 and when r_l==1.
! r_l is the relative coordinate in the buffer, going from 0 at its beginning to 1 at its end.
! ya contains the computed value.

subroutine stretching_function(r_l, ya)
  use constants, only: CUSTOM_REAL
  
  implicit none
  
  ! Input/Output.
  real(kind=CUSTOM_REAL), intent(in) :: r_l
  real(kind=CUSTOM_REAL), intent(out) :: ya
  
  ! Local Variables.
  real(kind=CUSTOM_REAL), parameter :: ONE  = 1._CUSTOM_REAL
  real(kind=CUSTOM_REAL) :: eps_l, p, q ! Arina's stretching.
  
  ! Coefficients for the stretching function.
  ! Arina's stretching
  !eps_l = 1.d-4 ! Ending value of stretching. 1.d-4 in Arina's paper
  eps_l = 1.d-3 ! Ending value of stretching. 1.d-4 in Arina's paper
  !eps_l = 0.16 ! Ending value of stretching. 1.d-4 in Arina's paper..
  !eps_l = 0.25 ! Ending value of stretching. 1.d-4 in Arina's paper.
  !eps_l = 0.5 ! Ending value of stretching. 1.d-4 in Arina's paper.
  p = 3.25
  q = 6.
  
  ! Arina's stretching. See 10.3390/aerospace3010007.
  ya = ONE - (ONE - eps_l) * (ONE - (ONE - r_l)**p)**q ! Stretching function.
end subroutine stretching_function


! ------------------------------------------------------------ !
! damp_function                                                !
! ------------------------------------------------------------ !
! Damping function that can be used for priori and posteriori dampings (see Chapter 2 Appendices in Martire's thesis).
! For an example on poseriori damping, see the commented-out portion at the end of the 'compute_forces_acoustic_LNS_main' routine in 'compute_forces_acoustic_LNS_calling_routine.F90'.

subroutine damp_function(r_l, sigma)
  use constants,only: CUSTOM_REAL
  implicit none  
  ! Input/Output.
  real(kind=CUSTOM_REAL), intent(in):: r_l
  real(kind=CUSTOM_REAL), intent(out):: sigma
  ! Local
  real(kind=CUSTOM_REAL), parameter :: ONE  = 1._CUSTOM_REAL
  real(kind=CUSTOM_REAL) :: C_1, C_2 ! Arina's damping.
  real(kind=CUSTOM_REAL) :: beta, sigma_max ! Richards' damping.
  
  ! Coefficients for the stretching function.
  ! Arina's damping coefficients.
  C_1 = 0.!0.0d0 ! 0 in Arina's paper, 0<=C_1<=0.1 in Wasistho's.
  C_2 = 13.!6.!10.!13.0d0 ! 13 in Arina's paper, 10<=C_2<=20 in Wasistho's.
  ! Richards' damping coefficients.
  beta = 2.!3.0d0 ! 4 in Richards' paper.
  sigma_max = -1.0d0
  
  ! Arina's damping. See 10.3390/aerospace3010007.
  sigma = (1.0d0 - C_1*r_l**2.0d0)*(1.0d0 - ( 1.0d0 - exp(C_2*(r_l)**2.0d0) )/( 1.0d0 - exp(C_2) ))
  ! Richards' damping. See 10.1016/j.jsv.2003.09.042.
  !sigma = 1.0d0 + sigma_max * r_l ** beta
end subroutine damp_function


! ------------------------------------------------------------ !
! change_visco                                                 !
! ------------------------------------------------------------ !
! Another tentative function for damping waves in absorbing boundary condition buffers.
! Increase viscosity to "physically" damp waves.
! This routine remains unfinished.

subroutine change_visco(i, j, ispec, x, z)
  use constants, only: CUSTOM_REAL
  use specfem_par, only: etaext, muext, mesh_xmin, ABC_STRETCH_LEFT, ABC_STRETCH_LEFT_LBUF
  implicit none
  ! Input/Output.
  integer, intent(in) :: i, j, ispec
  real(kind=CUSTOM_REAL), intent(in) :: x, z
  ! Local Variables.
  real(kind=CUSTOM_REAL) :: r_l
  real(kind=CUSTOM_REAL), parameter :: ZERO  = 0._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: ONE  = 1._CUSTOM_REAL
  if(ABC_STRETCH_LEFT) then
    r_l = -(x-mesh_xmin)/ABC_STRETCH_LEFT_LBUF + ONE ! 0 at beginning of buffer, 1 at end.
    if(r_l>ZERO .and. r_l<=ONE)then
      muext(i, j, ispec) = 1.25d-5*(1.+1.d3*r_l**2.)
      etaext(i, j, ispec) = (4./3.)*muext(i, j, ispec)
    endif
  endif
  ! TODO: implement for other sides.
  if(.false.) write(*, *) x, z ! Horrible hack, I'm so sorry.
  
end subroutine change_visco


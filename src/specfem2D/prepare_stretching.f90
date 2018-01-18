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

! ------------------------------------------------------------ !
! prepare_stretching                                           !
! ------------------------------------------------------------ !
! TODO: Description.

subroutine prepare_stretching()
  use specfem_par, only: coord, ibool_before_perio,&!ibool_DG,&
                         myrank, nglob,nspec,&
                         ispec_is_acoustic_DG,&
                         ADD_PERIODIC_CONDITIONS, &
                         ABC_STRETCH_LEFT, ABC_STRETCH_RIGHT, ABC_STRETCH_TOP, ABC_STRETCH_BOTTOM, &
                         ABC_STRETCH_LBUF, &
                         stretching_ya,any_elastic,&
                         mesh_xmin, mesh_xmax, mesh_zmin, mesh_zmax!,stretching_buffer
  
  use constants,only: CUSTOM_REAL,NGLLX,NGLLZ

  implicit none
  
  ! Local
  integer iglob_unique,ispec,i,j!,iglob
  real(kind=CUSTOM_REAL), parameter :: ZERO = 0._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: ONE  = 1._CUSTOM_REAL
  !real(kind=CUSTOM_REAL) :: mesh_xmin_local, mesh_xmax_local, mesh_zmin_local, mesh_zmax_local, &
  !                          mesh_xmin, mesh_xmax, mesh_zmin, mesh_zmax
  real(kind=CUSTOM_REAL) :: x,z
  real(kind=CUSTOM_REAL) :: r_l ! Relative buffer coordinate (which is equal to 0 at beginning of buffer, and to 1 at the end).
  
  ! Error checking.
  if(myrank==0) then
    if(ADD_PERIODIC_CONDITIONS .and. (ABC_STRETCH_LEFT .or. ABC_STRETCH_RIGHT)) then
      write(*,*) "********************************"
      write(*,*) "*            ERROR             *"
      write(*,*) "********************************"
      write(*,*) "* Cannot use (horizontal) both *"
      write(*,*) "* periodic boundary conditions *"
      write(*,*) "* and stretching BC on one or  *"
      write(*,*) "* both lateral boundaries.     *"
      write(*,*) "********************************"
      stop
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
      stop
    endif
  endif
  
  allocate(stretching_ya(2,nglob)) ! Two-dimensionnal stretching. Change it to 3 for 3D.
  !allocate(stretching_buffer(nglob)) ! Stretching buffers code: -1 outside buffers, 1 in top buffer, 2 in left buffer, 3 in bottom buffer, 4 in right buffer.
  stretching_ya(:, :) = ONE ! By default, mesh is not stretched.
  !stretching_buffer(:) = 0 ! By default, mesh is not stretched.
  
  do ispec = 1, nspec
    if(ispec_is_acoustic_DG(ispec)) then
      do j = 1, NGLLZ
        do i = 1, NGLLX
          !iglob = ibool_DG(i, j, ispec)
          iglob_unique = ibool_before_perio(i, j, ispec)
          x=coord(1, iglob_unique)
          z=coord(2, iglob_unique)
          if(     (ABC_STRETCH_LEFT   .and. x < mesh_xmin + ABC_STRETCH_LBUF) & ! left stretching and in left buffer zone
             .or. (ABC_STRETCH_RIGHT  .and. x > mesh_xmax - ABC_STRETCH_LBUF) & ! right stretching and in right buffer zone
             .or. (ABC_STRETCH_BOTTOM .and. z < mesh_zmin + ABC_STRETCH_LBUF) & ! bottom stretching and in bottom buffer zone
             .or. (ABC_STRETCH_TOP    .and. z > mesh_zmax - ABC_STRETCH_LBUF)) then ! top stretching and in top buffer zone
            
            if(ABC_STRETCH_TOP) then
              r_l = (z - mesh_zmax)/ABC_STRETCH_LBUF + ONE
              if(r_l>ZERO .and. r_l<=ONE) call stretching_function(r_l, stretching_ya(2, iglob_unique))
              !stretching_buffer(iglob) = ibset(stretching_buffer(iglob), 0) ! Set 1st LSB to 1.
            endif
            if(ABC_STRETCH_LEFT) then
              r_l = ONE - (x - mesh_xmin)/ABC_STRETCH_LBUF
              if(r_l>ZERO .and. r_l<=ONE) call stretching_function(r_l, stretching_ya(1, iglob_unique))
              !stretching_buffer(iglob) = ibset(stretching_buffer(iglob), 1) ! Set 2nd LSB to 1.
            endif
            if(ABC_STRETCH_BOTTOM) then
              r_l = ONE - (z - mesh_zmin)/ABC_STRETCH_LBUF
              if(r_l>ZERO .and. r_l<=ONE) call stretching_function(r_l, stretching_ya(2, iglob_unique))
              !stretching_buffer(iglob) = ibset(stretching_buffer(iglob), 2) ! Set 3rd LSB to 1.
            endif
            if(ABC_STRETCH_RIGHT) then
              r_l = (x - mesh_xmax)/ABC_STRETCH_LBUF + ONE
              if(r_l>ZERO .and. r_l<=ONE) call stretching_function(r_l, stretching_ya(1, iglob_unique))
              !stretching_buffer(iglob) = ibset(stretching_buffer(iglob), 3) ! Set 4th LSB to 1.
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
  use constants,only: CUSTOM_REAL
  
  implicit none
  
  ! Input/output.
  real(kind=CUSTOM_REAL), intent(in) :: r_l
  real(kind=CUSTOM_REAL), intent(out) :: ya
  
  ! Local variables.
  real(kind=CUSTOM_REAL), parameter :: ONE  = 1._CUSTOM_REAL
  real(kind=CUSTOM_REAL) :: eps_l, p, q ! Arina's stretching.
  
  ! Coefficients for the stretching function.
  ! Arina's stretching
  eps_l = 0.2!1.0d-4!0.2!0.3!1.0d-4 ! 1.d-4 in Arina's paper.
  p = 3.25d0
  q = 6.!1.75!1.75!8.!5.!1.75d0
  
  ! Arina's stretching.
  ya = ONE - (ONE - eps_l) * (ONE - (ONE - r_l)**p)**q ! Stretching function.
  !ya=ONE
end subroutine stretching_function

subroutine damp_function(r_l, sigma)
  use constants,only: CUSTOM_REAL
  implicit none  
  ! Input/output.
  real(kind=CUSTOM_REAL), intent(in):: r_l
  real(kind=CUSTOM_REAL), intent(out):: sigma
  ! Local
  real(kind=CUSTOM_REAL), parameter :: ONE  = 1._CUSTOM_REAL
  real(kind=CUSTOM_REAL) :: C_1, C_2 ! Arina's damping.
  real(kind=CUSTOM_REAL) :: beta, sigma_max ! Richards' damping.
  
  ! Coefficients for the stretching function.
  ! Arina's damping coefficients.
  C_1 = 0.08!0.0d0 ! 0 in Arina's paper, 0<=C_1<=0.1 in Wasistho's.
  C_2 = 13.!6.!10.!13.0d0 ! 13 in Arina's paper, 10<=C_2<=20 in Wasistho's.
  ! Richards' damping coefficients.
  beta = 4.!3.0d0 ! 4 in Richards' paper.
  sigma_max = -1.0d0
  
  ! Arina's damping.
  sigma = (1.0d0 - C_1*r_l**2.0d0)*(1.0d0 - ( 1.0d0 - exp(C_2*(r_l)**2.0d0) )/( 1.0d0 - exp(C_2) ))
  ! Richards' damping.
  !sigma = 1.0d0 + sigma_max * r_l ** beta
  !write(*, *) "z", z, "coef_stretch_z", coef_stretch_z ! DEBUG
end subroutine damp_function

! ------------------------------------------------------------ !
! change_visco                                                 !
! ------------------------------------------------------------ !
! TODO: Description.

subroutine change_visco(i, j, ispec, x, z)
  use constants,only: CUSTOM_REAL
  
  use specfem_par, only: etaext, muext,&
                         !ABC_STRETCH,&
                         ABC_STRETCH_TOP, ABC_STRETCH_LBUF,&
                         !mesh_xmin, mesh_xmax, mesh_zmin,&
                         mesh_zmax

  implicit none
  
  ! Input/output.
  integer, intent(in) :: i, j, ispec
  real(kind=CUSTOM_REAL), intent(in) :: x, z
  
  ! Local variables.
  real(kind=CUSTOM_REAL) :: r_l
  real(kind=CUSTOM_REAL), parameter :: ZERO  = 0._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: ONE  = 1._CUSTOM_REAL
            
  ! TEST VISCO
  if(ABC_STRETCH_TOP) then
    r_l = (z - mesh_zmax)/ABC_STRETCH_LBUF + ONE
    if(r_l>ZERO .and. r_l<=ONE)then
      muext(i, j, ispec) = 1.25d-5*(1.+1000.*r_l**2.)
      etaext(i, j, ispec) = (4./3.)*muext(i, j, ispec)
      !write(*,*) z, muext(i,j,ispec)
    endif
    !if(r_l>-0.5 .and. r_l<0) 
  endif
end subroutine change_visco

! ------------------------------------------------------------ !
! damp_solution_DG                                             !
! ------------------------------------------------------------ !
! TODO: This is a test routine.

subroutine damp_solution_DG(rho_DG, rhovx_DG, rhovz_DG, E_DG, timelocal)

  use specfem_par, only: nspec, coord, ibool_DG, nglob_DG, &
        ibool_before_perio,ABC_STRETCH_LBUF,&
        ABC_STRETCH_LEFT, ABC_STRETCH_RIGHT, ABC_STRETCH_TOP, ABC_STRETCH_BOTTOM,&
        ispec_is_acoustic_DG,&
        mesh_xmin, mesh_xmax, mesh_zmin, mesh_zmax
  
  use constants,only: CUSTOM_REAL,NGLLX,NGLLZ

  implicit none
  
  ! Input/output.
  real(kind=CUSTOM_REAL), dimension(nglob_DG) :: rho_DG, rhovx_DG, rhovz_DG, E_DG
  
  ! Local
  integer :: iglob
  real(kind=CUSTOM_REAL) :: sigma
  real(kind=CUSTOM_REAL) :: timelocal
  real(kind=CUSTOM_REAL) :: rho_DG_P, rhovx_DG_P, rhovz_DG_P, E_DG_P, &
      veloc_x_DG_P, veloc_z_DG_P, p_DG_P, e1_DG_P
  integer iglob_unique,ispec,i,j
  real(kind=CUSTOM_REAL), parameter :: ZERO = 0._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: ONE  = 1._CUSTOM_REAL
  real(kind=CUSTOM_REAL) :: x,z
  real(kind=CUSTOM_REAL) :: r_l
  
  ! Determine mesh's min/max coordinates and collect them.
  do ispec = 1, nspec
    if(ispec_is_acoustic_DG(ispec)) then
      do j = 1, NGLLZ
        do i = 1, NGLLX
          iglob_unique = ibool_before_perio(i, j, ispec)
          iglob = ibool_DG(i, j, ispec)
          x=coord(1, iglob_unique)
          z=coord(2, iglob_unique)
          if(     (ABC_STRETCH_LEFT   .and. x < mesh_xmin + ABC_STRETCH_LBUF) & ! left stretching and in left buffer zone
             .or. (ABC_STRETCH_RIGHT  .and. x > mesh_xmax - ABC_STRETCH_LBUF) & ! right stretching and in right buffer zone
             .or. (ABC_STRETCH_BOTTOM .and. z < mesh_zmin + ABC_STRETCH_LBUF) & ! bottom stretching and in bottom buffer zone
             .or. (ABC_STRETCH_TOP    .and. z > mesh_zmax - ABC_STRETCH_LBUF)) then ! top stretching and in top buffer zone
            
            sigma = ONE ! In case something bad happens.
            ! Load damping value.
            if(ABC_STRETCH_LEFT) then
              r_l = ONE - (x - mesh_xmin)/ABC_STRETCH_LBUF
              if(r_l>ZERO .and. r_l<=ONE) call damp_function(r_l, sigma)
            endif
            if(ABC_STRETCH_RIGHT) then
              r_l = (x - mesh_xmax)/ABC_STRETCH_LBUF + ONE
              if(r_l>ZERO .and. r_l<=ONE) call damp_function(r_l, sigma)
            endif
            if(ABC_STRETCH_BOTTOM) then
              r_l = ONE - (z - mesh_zmin)/ABC_STRETCH_LBUF
              if(r_l>ZERO .and. r_l<=ONE) call damp_function(r_l, sigma)
            endif
            if(ABC_STRETCH_TOP) then
              r_l = (z - mesh_zmax)/ABC_STRETCH_LBUF + ONE
              if(r_l>ZERO .and. r_l<=ONE) call damp_function(r_l, sigma)
            endif
            ! Load quiet state.
            call boundary_condition_DG(i, j, ispec, timelocal, rho_DG_P, rhovx_DG_P, rhovz_DG_P, E_DG_P, &
                    veloc_x_DG_P, veloc_z_DG_P, p_DG_P, e1_DG_P)
            ! Damp perturbation.
            rho_DG(iglob)   = rho_DG_P + sigma*( rho_DG(iglob) - rho_DG_P)
            rhovx_DG(iglob) = rhovx_DG_P + sigma*( rhovx_DG(iglob) - rhovx_DG_P)
            rhovz_DG(iglob) = rhovz_DG_P + sigma*( rhovz_DG(iglob) - rhovz_DG_P)
            E_DG(iglob)     = E_DG_P + sigma*( E_DG(iglob) - E_DG_P)
          endif
        enddo
      enddo
    endif
  enddo
end subroutine damp_solution_DG

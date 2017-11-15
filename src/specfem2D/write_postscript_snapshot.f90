!========================================================================
!
!                   S P E C F E M 2 D  Version 7 . 0
!                   --------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!                and CNRS / University of Marseille, France
!                 (there are currently maNZ_IMAGE_color more authors!)
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


  subroutine write_postscript_snapshot()

  use constants,only: IMAIN

  use specfem_par, only: myrank,P_SV,it, &
                         potential_acoustic,potential_gravitoacoustic, &
                         potential_gravito,displ_elastic,displs_poroelastic, &
                         potential_dot_acoustic,potential_dot_gravitoacoustic, &
                         potential_dot_gravito,veloc_elastic,velocs_poroelastic, &
                         potential_dot_dot_acoustic,potential_dot_dot_gravitoacoustic, &
                         potential_dot_dot_gravito,accel_elastic,accels_poroelastic, &
                         potential_dphi_dx_DG!, potential_dphi_dz_DG ! Modification for DG.

  use specfem_par_movie,only: imagetype_postscript

  implicit none

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'Writing PostScript vector plot for time step ',it
    call flush_IMAIN()
  endif

  ! determines postscript output type
  if (imagetype_postscript == 1 .and. P_SV) then

    if (myrank == 0) write(IMAIN,*) 'drawing displacement vector as small arrows...'
    !call compute_vector_whole_medium(potential_acoustic,potential_gravitoacoustic, &
    !                                 potential_gravito,displ_elastic,displs_poroelastic)
    ! Previous call prevents ifort compilation (but, strangely, does not bother gfortran compilation). Thus, we make the following call instead.
    ! TODO: Do something here instead of this poor patch.
    call compute_vector_whole_medium(potential_acoustic,potential_gravitoacoustic, &
                                     potential_gravito,displ_elastic,displs_poroelastic, &
                                     potential_dphi_dx_DG)

    call plot_post()

  else if (imagetype_postscript == 2 .and. P_SV) then

    if (myrank == 0) write(IMAIN,*) 'drawing velocity vector as small arrows...'
    !call compute_vector_whole_medium(potential_dot_acoustic,potential_dot_gravitoacoustic, &
    !                                 potential_dot_gravito,veloc_elastic,velocs_poroelastic)
    ! Previous call prevents ifort compilation (but, strangely, does not bother gfortran compilation). Thus, we make the following call instead.
    ! TODO: Do something here instead of this poor patch.
    call compute_vector_whole_medium(potential_dot_acoustic,potential_dot_gravitoacoustic, &
                                     potential_dot_gravito,veloc_elastic,velocs_poroelastic, &
                                     potential_dphi_dx_DG)

    call plot_post()

  else if (imagetype_postscript == 3 .and. P_SV) then

    if (myrank == 0) write(IMAIN,*) 'drawing acceleration vector as small arrows...'
    !call compute_vector_whole_medium(potential_dot_dot_acoustic,potential_dot_dot_gravitoacoustic, &
    !                                 potential_dot_dot_gravito,accel_elastic,accels_poroelastic)
    ! Previous call prevents ifort compilation (but, strangely, does not bother gfortran compilation). Thus, we make the following call instead.
    ! TODO: Do something here instead of this poor patch.
    call compute_vector_whole_medium(potential_dot_dot_acoustic,potential_dot_dot_gravitoacoustic, &
                                     potential_dot_dot_gravito,accel_elastic,accels_poroelastic, &
                                     potential_dphi_dx_DG)

    call plot_post()

  else if (.not. P_SV) then
    call exit_MPI(myrank,'cannot draw a SH scalar field as a vector plot, turn PostScript plots off')
  else
    call exit_MPI(myrank,'wrong type for PostScript snapshots')
  endif

  ! user output
  if (myrank == 0 .and. imagetype_postscript /= 4 .and. P_SV ) then
    write(IMAIN,*) 'PostScript file written'
    call flush_IMAIN()
  endif

  end subroutine write_postscript_snapshot


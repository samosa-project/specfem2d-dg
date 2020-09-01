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


! Partitioning using METIS

  subroutine metis_partitioning()


#ifdef USE_METIS
  use part_unstruct_par,only: nelmnts,part,nb_edges,vwgt,adjwgt, &
    xadj => xadj_g,adjncy => adjncy_g

  use parameter_file_par,only: nparts => NPROC
#endif

  implicit none

  include "constants.h"

!  integer, intent(in)  :: nelmnts, nparts, nb_edges
!  integer, dimension(0:nelmnts), intent(in)  :: xadj
!  integer, dimension(0:MAX_NEIGHBORS*nelmnts-1), intent(in)  :: adjncy
!  integer, dimension(0:nelmnts-1), intent(in)  :: vwgt
!  integer, dimension(0:nb_edges-1), intent(in)  :: adjwgt
!  integer, dimension(:), pointer  :: part

#ifdef USE_METIS
  integer, dimension(0:4)  :: metis_options
#endif

  integer :: wgtflag
  integer :: remove_min_to_start_at_zero
  integer :: edgecut

!! DK DK support for METIS now removed, we use SCOTCH instead
  stop 'support for the METIS graph partitioner has been discontinued, please use SCOTCH (option 3) instead'

  ! initializes
  remove_min_to_start_at_zero = 0
  wgtflag = 0
  edgecut = 0

#ifdef USE_METIS
  call METIS_PartGraphRecursive(nelmnts, xadj(0), adjncy(0), vwgt(0), adjwgt(0), wgtflag, remove_min_to_start_at_zero, nparts, &
      metis_options, edgecut, part(0))

  !call METIS_PartGraphVKway(nelmnts, xadj(0), adjncy(0), vwgt(0), adjwgt(0), wgtflag, remove_min_to_start_at_zero, nparts, &
  !     options, edgecut, part(0))
#else
  ! safety stop
  print *, 'This version of SPECFEM was not compiled with support of METIS.'
  print *, 'Please recompile with -DUSE_METIS in order to enable use of METIS.'
  stop 'Metis partitioning not compiled'
#endif

 end subroutine metis_partitioning



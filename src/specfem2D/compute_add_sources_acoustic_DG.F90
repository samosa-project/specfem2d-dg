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
!=====================================================================

! ------------------------------------------------------------ !
! compute_add_sources_acoustic_DG_backward                     !
! ------------------------------------------------------------ !
! Add sources in the backward method.

subroutine compute_add_sources_acoustic_DG_backward(it_tmp, b_dot_rho, b_dot_rhovx, b_dot_rhovz, b_dot_E)

  use constants, only: NGLLX, NGLLZ, CUSTOM_REAL

  use specfem_par, only: myrank, ispec_is_acoustic, &
                         nrec, which_proc_receiver, ispec_selected_rec, adj_sourcearrays, &
                         ibool, nglob
  
  implicit none
  
  ! Input/Output.
  integer :: it_tmp
  real(kind=CUSTOM_REAL), dimension(nglob) :: b_dot_rho, b_dot_rhovx, b_dot_rhovz, b_dot_E
  
  ! Local Variables.
  integer :: irec_local, irec, i, j, iglob
  character(len=100) file_name
  
  write(file_name, "('source_', i3.3, '.txt')") myrank
  
  open(11, file=file_name, form='formatted', position='append')
  
  irec_local = 0
  do irec = 1, nrec
    ! add the source (only if this proc carries the source)
    if (myrank == which_proc_receiver(irec)) then
      irec_local = irec_local + 1
      if (ispec_is_acoustic(ispec_selected_rec(irec))) then
        ! add source array
        do j = 1, NGLLZ
          do i = 1, NGLLX
            iglob = ibool(i, j, ispec_selected_rec(irec))
            
            b_dot_rho(iglob)   = b_dot_rho(iglob) &
                                + 0!adj_sourcearrays(irec_local, it_tmp, 1, i, j) &
                                !!ZN becareful the following line is new added, thus when do comparison
                                !!ZN of the new code with the old code, you will have big difference if you
                                !!ZN do not tune the source
                                !/ kappastore(i, j, ispec_selected_rec(irec))
            b_dot_rhovx(iglob) = b_dot_rhovx(iglob) &
                                + 0!adj_sourcearrays(irec_local, it_tmp, 1, i, j) &
                                !!ZN becareful the following line is new added, thus when do comparison
                                !!ZN of the new code with the old code, you will have big difference if you
                                !!ZN do not tune the source
                                !/ kappastore(i, j, ispec_selected_rec(irec))
                                
            if(i == 1 .AND. j == 1) &
                WRITE(11, *)  adj_sourcearrays(irec_local, it_tmp, 2, 3, 5)  
                                
            b_dot_rhovz(iglob) = b_dot_rhovz(iglob) + adj_sourcearrays(irec_local, it_tmp, 2, i, j)!&
                                !+ adj_sourcearrays(irec_local, it_tmp, 2, i, j)! &
                                !ZN becareful the following line is new added, thus when do comparison
                                !ZN of the new code with the old code, you will have big difference if you
                                !ZN do not tune the source
                                !/ kappastore(i, j, ispec_selected_rec(irec))
            !if(i == 2 .AND. j == 2) &
            b_dot_E(iglob)     = b_dot_E(iglob) !&
                                !adj_sourcearrays(irec_local, it_tmp, 1, i, j) &
                                !!ZN becareful the following line is new added, thus when do comparison
                                !!ZN of the new code with the old code, you will have big difference if you
                                !!ZN do not tune the source
                                !/ kappastore(i, j, ispec_selected_rec(irec))

          enddo
        enddo
      endif ! if element acoustic
    endif ! if this processor core carries the adjoint source
  enddo ! irec = 1, nrec

  close(11)
 
end subroutine compute_add_sources_acoustic_DG_backward
  
  
! ------------------------------------------------------------ !
! compute_add_sources_acoustic_DG_backward_real                !
! ------------------------------------------------------------ !
! Same as 'compute_add_sources_acoustic_DG_backward', but over the whole DG mesh (nglob_DG instead of nglob).
  
subroutine compute_add_sources_acoustic_DG_backward_real(it_tmp, b_dot_rho, b_dot_rhovx, b_dot_rhovz, b_dot_E)

  use constants, only: NGLLX, NGLLZ, CUSTOM_REAL

  use specfem_par, only: myrank, ispec_is_acoustic, &
                         nrec, which_proc_receiver, ispec_selected_rec, adj_sourcearrays, &
                         ibool_DG, nglob_DG
  
  implicit none


  ! Input/Output.
  integer :: it_tmp
  real(kind=CUSTOM_REAL), dimension(nglob_DG) :: b_dot_rho, b_dot_rhovx, b_dot_rhovz, b_dot_E
  
  ! Local Variables.
  integer :: irec_local, irec, i, j, iglob
  character(len=100) file_name
  
  write(file_name, "('source_', i3.3, '.txt')") myrank
  
  open(11, file=file_name, form='formatted', position='append')

  irec_local = 0
  do irec = 1, nrec
    ! add the source (only if this proc carries the source)
    if (myrank == which_proc_receiver(irec)) then
      irec_local = irec_local + 1
      if (ispec_is_acoustic(ispec_selected_rec(irec))) then
        ! add source array
        do j = 1, NGLLZ
          do i = 1, NGLLX
            iglob = ibool_DG(i, j, ispec_selected_rec(irec))
            
            b_dot_rho(iglob)   = b_dot_rho(iglob) &
                                + 0!adj_sourcearrays(irec_local, it_tmp, 1, i, j) &
                                !!ZN becareful the following line is new added, thus when do comparison
                                !!ZN of the new code with the old code, you will have big difference if you
                                !!ZN do not tune the source
                                !/ kappastore(i, j, ispec_selected_rec(irec))
            b_dot_rhovx(iglob) = b_dot_rhovx(iglob) &
                                + 0!adj_sourcearrays(irec_local, it_tmp, 1, i, j) &
                                !!ZN becareful the following line is new added, thus when do comparison
                                !!ZN of the new code with the old code, you will have big difference if you
                                !!ZN do not tune the source
                                !/ kappastore(i, j, ispec_selected_rec(irec))
                                
            if(i == 1 .AND. j == 1) &
                WRITE(11, *)  adj_sourcearrays(irec_local, it_tmp, 2, 3, 5)  
                                
            b_dot_rhovz(iglob) = b_dot_rhovz(iglob) + &
                                !exp(-((deltat*(it-1) - 5 )/4)**2)
                                adj_sourcearrays(irec_local, it_tmp, 2, i, j)!&
                                !+ adj_sourcearrays(irec_local, it_tmp, 2, i, j)! &
                                !ZN becareful the following line is new added, thus when do comparison
                                !ZN of the new code with the old code, you will have big difference if you
                                !ZN do not tune the source
                                !/ kappastore(i, j, ispec_selected_rec(irec))
            b_dot_E(iglob)     = b_dot_E(iglob) !&
                                !adj_sourcearrays(irec_local, it_tmp, 1, i, j) &
                                !!ZN becareful the following line is new added, thus when do comparison
                                !!ZN of the new code with the old code, you will have big difference if you
                                !!ZN do not tune the source
                                !/ kappastore(i, j, ispec_selected_rec(irec))

          enddo
        enddo
      endif ! if element acoustic
    endif ! if this processor core carries the adjoint source
  enddo ! irec = 1, nrec

  close(11)
 
end subroutine compute_add_sources_acoustic_DG_backward_real


! ------------------------------------------------------------ !
! compute_add_sources_acoustic_DG_spread                       !
! ------------------------------------------------------------ !
! Adds sources contributions to the constitutive variable of interest.
! Variable "variable_DG" is input/output-intended: it corresponds to the variable to which add source terms.
  
subroutine compute_add_sources_acoustic_DG_spread(variable_DG, it, i_stage)

  use constants, only: CUSTOM_REAL, NGLLX, NGLLZ, PI, HUGEVAL

  use specfem_par, only: nglob_DG, ispec_is_acoustic_DG, &
                         NSOURCES, source_type, source_time_function, &
                         is_proc_source, ispec_selected_source, &
                         ibool_DG, &
                         coord, &
                         jacobian, wxgll, wzgll, ibool_before_perio, &
                         USE_SPREAD_SSF, nspec, source_spatial_function_DG, &
                         ABC_STRETCH, stretching_ya
  
  implicit none

  ! Input/Output.
  real(kind=CUSTOM_REAL), dimension(nglob_DG), intent(inout) :: variable_DG
  integer, intent(in) :: it, i_stage
  
  ! Local Variables.
  real(kind=CUSTOM_REAL) :: x, y, r, accuracy, sigma, dist_min, dist, temp_source
  real(kind=CUSTOM_REAL), dimension(2) :: X1, X2, X3, X4, X0
  real(kind=CUSTOM_REAL), dimension(4, 2) :: Xc
  real(kind=CUSTOM_REAL) :: jacobianl, wxl, wzl
  integer :: i, j, i_source, ispec, iglob, iglob_unique
  real(kind=CUSTOM_REAL) :: stf ! In order to store the source time function at current timestep outside the many loops.
  
  do i_source = 1, NSOURCES ! Loop on sources.
    stf = source_time_function(i_source, it, i_stage) ! Store the source time function outside the many loops.
    ! See "prepare_source_time_function.f90" for the subroutine initialising the vector "source_time_function".
    
    if(USE_SPREAD_SSF) then
      ! Case in which a source spatially distributed over more than one element was initialised.
      do ispec = 1, nspec
        if(ispec_is_acoustic_DG(ispec)) then
          do j = 1, NGLLZ
            do i = 1, NGLLX
              iglob_unique = ibool_before_perio(i, j, ispec)
              iglob = ibool_DG(i, j, ispec)
              temp_source = stf * source_spatial_function_DG(i_source, iglob_unique)
              jacobianl = jacobian(i, j, ispec)
              
              if(ABC_STRETCH) then
                !call virtual_stretch(i, j, ispec, coef_stretch_x_ij, coef_stretch_z_ij)
                ! Add jacobian of stretching into the integrand (artifically).
                jacobianl = stretching_ya(1, iglob_unique)*stretching_ya(2, iglob_unique)*jacobianl
                ! Though, this should NOT have any importance, since NOONE should ever have sources in or near the buffers.
              endif
              
              wzl = real(wzgll(j), kind=CUSTOM_REAL)
              wxl = real(wxgll(i), kind=CUSTOM_REAL)
              variable_DG(iglob) = variable_DG(iglob) + temp_source * wxl * wzl * jacobianl
            enddo
          enddo
        endif
      enddo
    else
      ! Base case in which the source spatial function is distributed on one element only.
      if(ispec_is_acoustic_DG(ispec_selected_source(i_source))) then
        ! If the source element is acoustic.
        if(is_proc_source(i_source) == 1) then
          ! If this processor core carries the source.
          if(source_type(i_source) == 1) then
            ! If the source is an elastic force or an acoustic pressure.
            ispec = ispec_selected_source(i_source)
            
            ! Spatial source function is under the form exp(-(r/sigma)^2) where r is the distance to the source center point.
            ! Here, ispec is the element which is chosen to carry the source. Hence, the source center point must be chosen to
            ! be at the center of this element. Hence, find it and save its coordinates in the variable X0.
            X1(:) = coord(:, ibool_before_perio(1, 1, ispec)) ! Source element's bottom left corner.
            X2(:) = coord(:, ibool_before_perio(NGLLX, 1, ispec)) ! Source element's bottom right corner.
            X3(:) = coord(:, ibool_before_perio(1, NGLLZ, ispec)) ! Source element's top left corner.
            X4(:) = coord(:, ibool_before_perio(NGLLX, NGLLZ, ispec)) ! Source element's top right corner.
            Xc(1, :) = (X1(:) + X2(:))/2.
            Xc(2, :) = (X1(:) + X3(:))/2.
            Xc(3, :) = (X3(:) + X4(:))/2.
            Xc(4, :) = (X4(:) + X2(:))/2.
            X0(:) = (Xc(1, :) + Xc(3, :))/2.
            
            ! Choose sigma such that the value at the edge of the element is very small, in particular roughly equal to
            ! 10^(-accuracy), where accuracy is chosen below.
            accuracy = 7.
            dist_min = HUGEVAL
            do j = 1, 4
              dist = sqrt( (Xc(j, 1) - X0(1))**2 + (Xc(j, 2) - X0(2))**2 )
              if(dist < dist_min) then
                dist_min = dist
              endif
            enddo
            sigma = dist_min/sqrt(accuracy*log(10.))
            
            ! At each GLL point of the source element, add to the variable (variable_DG) the value of the source function.
            do j = 1, NGLLZ
              do i = 1, NGLLX
                iglob = ibool_DG(i, j, ispec)
                x = coord(1, ibool_before_perio(i, j, ispec))
                y = coord(2, ibool_before_perio(i, j, ispec))
                r = sqrt( (x - X0(1))**2 + (y - X0(2))**2 )
                temp_source = stf * exp(-(r/sigma)**2)
                jacobianl = jacobian(i, j, ispec)
                wzl = real(wzgll(j), kind=CUSTOM_REAL)
                wxl = real(wxgll(i), kind=CUSTOM_REAL)
                variable_DG(iglob) = variable_DG(iglob) + temp_source * wxl * wzl * jacobianl
              enddo
            enddo
            
          endif ! Endif on source_type.
        endif ! Endif on is_proc_source
      endif ! Endif on ispec_is_acoustic_DG.
    endif ! Endif on SIGMA_SSF.
  enddo ! Enddo on i_source.
end subroutine compute_add_sources_acoustic_DG_spread


! ------------------------------------------------------------ !
! compute_add_sources_acoustic_DG_spread                       !
! ------------------------------------------------------------ !
! Implements a source on the mass conservation equation:
! 1) Add said source on the mass equation.
! 2) Ensure compatibility by adding it also to the other equations (see Chapter 2 of Martire's thesis).

subroutine compute_add_sources_acoustic_DG_mass(d_rho, d_rhovx, d_rhovz, d_E, rho, vx, vz, E, it, i_stage)

  use constants, only: CUSTOM_REAL, NGLLX, NGLLZ, PI, HUGEVAL
  
  use specfem_par, only: nglob_DG, ispec_is_acoustic_DG, &
                         NSOURCES, myrank, &
                         source_time_function, &
                         ibool_DG, gammaext_DG, &
                         jacobian, wxgll, wzgll, ibool_before_perio, &
                         USE_SPREAD_SSF, nspec, source_spatial_function_DG
  
  implicit none

  ! Input/Output.
  real(kind=CUSTOM_REAL), dimension(nglob_DG), intent(inout) :: d_rho
  real(kind=CUSTOM_REAL), dimension(nglob_DG), intent(inout) :: d_rhovx
  real(kind=CUSTOM_REAL), dimension(nglob_DG), intent(inout) :: d_rhovz
  real(kind=CUSTOM_REAL), dimension(nglob_DG), intent(inout) :: d_E
  real(kind=CUSTOM_REAL), dimension(nglob_DG), intent(in) :: rho
  real(kind=CUSTOM_REAL), dimension(nglob_DG), intent(in) :: vx
  real(kind=CUSTOM_REAL), dimension(nglob_DG), intent(in) :: vz
  real(kind=CUSTOM_REAL), dimension(nglob_DG), intent(in) :: E
  integer, intent(in) :: it, i_stage
  
  ! Local Variables.
  real(kind=CUSTOM_REAL) :: temp_sourcewxlwzljacobianl
  integer :: i, j, i_source, ispec, iglob, iglob_unique
  real(kind=CUSTOM_REAL) :: stf ! In order to store the source time function at current timestep outside the many loops.
  
  do i_source = 1, NSOURCES ! Loop on sources.
    stf = source_time_function(i_source, it, i_stage) ! Store the source time function outside the many loops.
    if(USE_SPREAD_SSF) then
      ! Case in which a source spatially distributed over more than one element was initialised.
      do ispec = 1, nspec
        if(ispec_is_acoustic_DG(ispec)) then
          do j = 1, NGLLZ
            do i = 1, NGLLX
              iglob_unique = ibool_before_perio(i, j, ispec)
              iglob = ibool_DG(i, j, ispec)
              ! See "prepare_source_spatial_function.f90" for the subroutine initialising the vector "source_spatial_function_DG".
              temp_sourcewxlwzljacobianl =   stf &
                                           * source_spatial_function_DG(i_source, iglob_unique) &
                                           * real(wxgll(i), kind=CUSTOM_REAL) &
                                           * real(wzgll(j), kind=CUSTOM_REAL) &
                                           * jacobian(i, j, ispec)
              
              d_rho(iglob) = d_rho(iglob) + temp_sourcewxlwzljacobianl
              
              d_rhovx(iglob) = d_rhovx(iglob) + vx(iglob) * temp_sourcewxlwzljacobianl
              d_rhovz(iglob) = d_rhovz(iglob) + vz(iglob) * temp_sourcewxlwzljacobianl
              
              d_E(iglob) =   d_E(iglob) &
                           + (   gammaext_DG(iglob)*E(iglob)/rho(iglob) &
                               - 0.5*(gammaext_DG(iglob)-1.)*(vx(iglob)**2+vz(iglob)**2) &
                             ) * temp_sourcewxlwzljacobianl
            enddo
          enddo
        endif
      enddo
      
    else
      if(myrank==0) then
        write(*, *) "********************************"
        write(*, *) "*            ERROR             *"
        write(*, *) "********************************"
        write(*, *) "* Mass source is not yet       *"
        write(*, *) "* implemented with             *"
        write(*, *) "* USE_SPREAD_SSF=.false.. See  *"
        write(*, *) "* compute_add_sources_acoustic_DG.f90."
        write(*, *) "********************************"
        stop
      endif
      
    endif ! Endif on SIGMA_SSF.
  enddo ! Enddo on i_source.
end subroutine compute_add_sources_acoustic_DG_mass


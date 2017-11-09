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

! for acoustic solver

 

! ---------------------------------------------------------------

  !subroutine compute_add_sources_acoustic_DG_backward_2(it,i_stage)
!
!  use constants,only: CUSTOM_REAL,NGLLX,NGLLZ!!
!
!  use specfem_par, only: ispec_is_acoustic, nrec,&
!                         !is_proc_source,ispec_selected_source,&
!                         hxis_store,hgammas_store,ibool_DG,myrank, &
!                         !neighbor_DG, neighbor_DG_corner, is_corner, &
!                         ispec_is_acoustic, myrank, &
!                         which_proc_receiver, ispec_selected_rec, &
!                         b_rho_DG, b_rhovx_DG, b_rhovz_DG, b_E_DG, &
!                         source_time_function_rho_DG, source_time_function_E_DG, &
!                         source_time_function_rhovx_DG, source_time_function_rhovz_DG!, nglob_DG, source_time_function, coord, ibool, link_DG_CG
!  
!  implicit none!!
!
!  integer,intent(in) :: it,i_stage!
!
!  !local variables
!  integer :: irec_local,irec,i, j, ispec!,j_neighbor,ispec_neighbor,k!, iglob, cpt_loc, i2, i_neighbor
!  !integer, dimension(4) :: iglob_surface
!  double precision :: hlagrange
!
!  !do i_source= 1,NSOURCES
!  !  if(ispec_is_acoustic(ispec_selected_source(i_source))) then
!  !  
!  !  ! if this processor core carries the source and the source element is acoustic
!  !  if (is_proc_source(i_source) == 1) then
!  !     if (source_type(i_source) == 1) then
! 
!   irec_local = 0
!   do irec = 1,nrec
!   
!    ! add the source (only if this proc carries the source)
!    if (myrank == which_proc_receiver(irec)) then
!    
!      irec_local = irec_local + 1
!      
!      if (ispec_is_acoustic(ispec_selected_rec(irec))) then
! 
! 
!      ! collocated force
!      ! beware, for acoustic medium, source is: pressure divided by Kappa of the fluid
!      ! the sign is negative because pressure p = - Chi_dot_dot therefore we need
!      ! to add minus the source to Chi_dot_dot to get plus the source in pressure
!        ! forward wavefield
!        do j = 1,NGLLZ
!          do i = 1,NGLLX
!            !hlagrange = hxis_store(i_source,i) * hgammas_store(i_source,j)
!            
!            ispec = ispec_selected_rec(irec_local)
!            
!            b_rho_DG(ibool_DG(i, j, ispec))   = b_rho_DG(ibool_DG(i, j, ispec)) &
!                - source_time_function_rho_DG(irec_local,it,i_stage)*hlagrange
!            b_rhovx_DG(ibool_DG(i, j, ispec)) = b_rhovx_DG(ibool_DG(i, j, ispec)) &
!                - source_time_function_rhovx_DG(irec_local,it,i_stage)*hlagrange
!            b_rhovz_DG(ibool_DG(i, j, ispec)) = b_rhovz_DG(ibool_DG(i, j, ispec)) &
!                - source_time_function_rhovz_DG(irec_local,it,i_stage)*hlagrange
!            b_E_DG(ibool_DG(i, j, ispec))     = b_E_DG(ibool_DG(i, j, ispec)) &
!                - source_time_function_E_DG(irec_local,it,i_stage)*hlagrange
!            
!          enddo
!        enddo
!        
!      endif
!      
!      endif ! if(ispec_is_acoustic(ispec_selected_source(i_source)))
!      
!  enddo ! do i_source= 1,NSOURCES
!
!  end subroutine compute_add_sources_acoustic_DG_backward_2

! ---------------------------------------------------------------

  subroutine compute_add_sources_acoustic_DG_backward(it_tmp, &
        b_dot_rho, b_dot_rhovx, b_dot_rhovz, b_dot_E)

  use constants,only: NGLLX,NGLLZ,CUSTOM_REAL

  use specfem_par, only: myrank,ispec_is_acoustic,&
                         nrec,which_proc_receiver,ispec_selected_rec,adj_sourcearrays,&
                         ibool, nglob!ibool_DG
  implicit none

  real(kind=CUSTOM_REAL), dimension(nglob) :: b_dot_rho, b_dot_rhovx, b_dot_rhovz, b_dot_E

  !local variables
  integer :: irec_local,irec,i,j,iglob
  integer :: it_tmp

  character(len=100) file_name
  
  write(file_name,"('source_',i3.3,'.txt')") myrank
  
  open(11,file=file_name,form='formatted',position='append')
  
  ! time step index
  !it_tmp = NSTEP - it + 1

  irec_local = 0
  do irec = 1,nrec
    ! add the source (only if this proc carries the source)
    if (myrank == which_proc_receiver(irec)) then
      irec_local = irec_local + 1
      if (ispec_is_acoustic(ispec_selected_rec(irec))) then
        ! add source array
        do j = 1,NGLLZ
          do i = 1,NGLLX
            !iglob = ibool_DG(i, j, ispec_selected_rec(irec))
            iglob = ibool(i, j, ispec_selected_rec(irec))
            
            b_dot_rho(iglob)   = b_dot_rho(iglob) &
                                + 0!adj_sourcearrays(irec_local,it_tmp,1,i,j) &
                                !!ZN becareful the following line is new added, thus when do comparison
                                !!ZN of the new code with the old code, you will have big difference if you
                                !!ZN do not tune the source
                                !/ kappastore(i, j, ispec_selected_rec(irec))
            b_dot_rhovx(iglob) = b_dot_rhovx(iglob) &
                                + 0!adj_sourcearrays(irec_local,it_tmp,1,i,j) &
                                !!ZN becareful the following line is new added, thus when do comparison
                                !!ZN of the new code with the old code, you will have big difference if you
                                !!ZN do not tune the source
                                !/ kappastore(i, j, ispec_selected_rec(irec))
                                
            if(i == 1 .AND. j == 1) &
                WRITE(11,*)  adj_sourcearrays(irec_local,it_tmp,2,3,5)  
                                
            b_dot_rhovz(iglob) = b_dot_rhovz(iglob) + adj_sourcearrays(irec_local,it_tmp,2,i,j)!&
                                !+ adj_sourcearrays(irec_local,it_tmp,2,i,j)! &
                                !ZN becareful the following line is new added, thus when do comparison
                                !ZN of the new code with the old code, you will have big difference if you
                                !ZN do not tune the source
                                !/ kappastore(i, j, ispec_selected_rec(irec))
            !if(i == 2 .AND. j == 2) &
            b_dot_E(iglob)     = b_dot_E(iglob) !&
                                !adj_sourcearrays(irec_local,it_tmp,1,i,j) &
                                !!ZN becareful the following line is new added, thus when do comparison
                                !!ZN of the new code with the old code, you will have big difference if you
                                !!ZN do not tune the source
                                !/ kappastore(i, j, ispec_selected_rec(irec))

          enddo
        enddo
        !WRITE(*,*) i, j, ispec_selected_rec(irec),irec,myrank,"-- SOURCE -> ", &
        !        minval(adj_sourcearrays(irec_local,it_tmp,2,:,:)), maxval(adj_sourcearrays(irec_local,it_tmp,2,:,:)), &
        !        which_proc_receiver
      endif ! if element acoustic
    endif ! if this processor core carries the adjoint source
  enddo ! irec = 1,nrec

  close(11)
 
  end subroutine compute_add_sources_acoustic_DG_backward
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine compute_add_sources_acoustic_DG_backward_real(it_tmp, &
        b_dot_rho, b_dot_rhovx, b_dot_rhovz, b_dot_E)

  use constants,only: NGLLX,NGLLZ,CUSTOM_REAL

  use specfem_par, only: myrank,ispec_is_acoustic,&
                         nrec,which_proc_receiver,ispec_selected_rec,adj_sourcearrays,&
                         ibool_DG, nglob_DG!ibool_DG
  implicit none

  real(kind=CUSTOM_REAL), dimension(nglob_DG) :: b_dot_rho, b_dot_rhovx, b_dot_rhovz, b_dot_E

  !local variables
  integer :: irec_local,irec,i,j,iglob
  integer :: it_tmp

  character(len=100) file_name
  
  write(file_name,"('source_',i3.3,'.txt')") myrank
  
  open(11,file=file_name,form='formatted',position='append')
  
  ! time step index
  !it_tmp = NSTEP - it + 1

  irec_local = 0
  do irec = 1,nrec
    ! add the source (only if this proc carries the source)
    if (myrank == which_proc_receiver(irec)) then
      irec_local = irec_local + 1
      if (ispec_is_acoustic(ispec_selected_rec(irec))) then
        ! add source array
        do j = 1,NGLLZ
          do i = 1,NGLLX
            !iglob = ibool_DG(i, j, ispec_selected_rec(irec))
            iglob = ibool_DG(i, j, ispec_selected_rec(irec))
            
            b_dot_rho(iglob)   = b_dot_rho(iglob) &
                                + 0!adj_sourcearrays(irec_local,it_tmp,1,i,j) &
                                !!ZN becareful the following line is new added, thus when do comparison
                                !!ZN of the new code with the old code, you will have big difference if you
                                !!ZN do not tune the source
                                !/ kappastore(i, j, ispec_selected_rec(irec))
            b_dot_rhovx(iglob) = b_dot_rhovx(iglob) &
                                + 0!adj_sourcearrays(irec_local,it_tmp,1,i,j) &
                                !!ZN becareful the following line is new added, thus when do comparison
                                !!ZN of the new code with the old code, you will have big difference if you
                                !!ZN do not tune the source
                                !/ kappastore(i, j, ispec_selected_rec(irec))
                                
            if(i == 1 .AND. j == 1) &
                WRITE(11,*)  adj_sourcearrays(irec_local,it_tmp,2,3,5)  
                                
            b_dot_rhovz(iglob) = b_dot_rhovz(iglob) + &
                                !exp(-((deltat*(it-1) - 5 )/4)**2)
                                adj_sourcearrays(irec_local,it_tmp,2,i,j)!&
                                !+ adj_sourcearrays(irec_local,it_tmp,2,i,j)! &
                                !ZN becareful the following line is new added, thus when do comparison
                                !ZN of the new code with the old code, you will have big difference if you
                                !ZN do not tune the source
                                !/ kappastore(i, j, ispec_selected_rec(irec))
            !if(i == 2 .AND. j == 2) &
            b_dot_E(iglob)     = b_dot_E(iglob) !&
                                !adj_sourcearrays(irec_local,it_tmp,1,i,j) &
                                !!ZN becareful the following line is new added, thus when do comparison
                                !!ZN of the new code with the old code, you will have big difference if you
                                !!ZN do not tune the source
                                !/ kappastore(i, j, ispec_selected_rec(irec))

          enddo
        enddo
        !WRITE(*,*) i, j, ispec_selected_rec(irec),irec,myrank,"-- SOURCE -> ", &
        !        minval(adj_sourcearrays(irec_local,it_tmp,2,:,:)), maxval(adj_sourcearrays(irec_local,it_tmp,2,:,:)), &
        !        which_proc_receiver
      endif ! if element acoustic
    endif ! if this processor core carries the adjoint source
  enddo ! irec = 1,nrec

  close(11)
 
  end subroutine compute_add_sources_acoustic_DG_backward_real
  
! ------------------------------------------------------------ !
! compute_add_sources_acoustic_DG_spread                       !
! ------------------------------------------------------------ !
  
  subroutine compute_add_sources_acoustic_DG_spread(variable_DG, it, i_stage)

  use constants,only: CUSTOM_REAL,NGLLX,NGLLZ,PI

  use specfem_par, only: nglob_DG,ispec_is_acoustic_DG,&!ispec_is_acoustic
                         NSOURCES,source_type,source_time_function,&
                         is_proc_source,ispec_selected_source,&
                         ibool_DG, &
                         coord,  &
                         jacobian, wxgll, wzgll, ibool_before_perio
  implicit none

  real(kind=CUSTOM_REAL), dimension(nglob_DG),intent(inout) :: variable_DG
  integer,intent(in) :: it,i_stage
  
  real(kind=CUSTOM_REAL) :: x, y, r, accuracy, sigma, dist_min, dist, temp_source
  real(kind=CUSTOM_REAL), dimension(2) :: X1, X2, X3, X4, X0
  real(kind=CUSTOM_REAL), dimension(4,2) :: Xc
  
  real(kind=CUSTOM_REAL) :: jacobianl, wxl, wzl
  integer :: i, j, i_source, ispec, iglob
  
  !do ispec = ifirstelem,ilastelem
  do i_source= 1, NSOURCES ! Loop on sources.
    if(ispec_is_acoustic_DG(ispec_selected_source(i_source))) then
      ! If the source element is acoustic.
      if(is_proc_source(i_source) == 1) then
        ! If this processor core carries the source.
        if(source_type(i_source) == 1) then
          ! If the source is an elastic force or an acoustic pressure.
          ispec = ispec_selected_source(i_source)

          ! Find coordinates of the center
          X1(:) = coord(:, ibool_before_perio(1, 1, ispec))
          X2(:) = coord(:, ibool_before_perio(NGLLX, 1, ispec))
          X3(:) = coord(:, ibool_before_perio(1, NGLLZ, ispec))
          X4(:) = coord(:, ibool_before_perio(NGLLX, NGLLZ, ispec))
          Xc(1, :) = (X1(:) + X2(:))/2.
          Xc(2, :) = (X1(:) + X3(:))/2.
          Xc(3, :) = (X3(:) + X4(:))/2.
          Xc(4, :) = (X4(:) + X2(:))/2.
          X0(:) = (Xc(1, :) + Xc(3, :))/2.
          dist_min = 1d10
          do j = 1,4
            dist = sqrt( (Xc(j, 1) - X0(1))**2 + (Xc(j,2) - X0(2))**2 )
            if(dist < dist_min) then
              dist_min = dist
            endif
          enddo
          
          accuracy = 7.
          sigma    = dist_min/sqrt(accuracy*log(10.))
          
          do j = 1, NGLLZ
            do i = 1, NGLLX
              iglob = ibool_DG(i, j, ispec)
              x = coord(1, ibool_before_perio(i, j, ispec))
              y = coord(2, ibool_before_perio(i, j, ispec))
              r = sqrt( (x - X0(1))**2 + (y - X0(2))**2 )
              temp_source = source_time_function(i_source, it, i_stage) * exp(-(r/sigma)**2) ! See "prepare_source_time_function.f90" for the subroutine initialising the vector "source_time_function".
              
              jacobianl = jacobian(i, j, ispec)
              wzl = real(wzgll(j), kind=CUSTOM_REAL)
              wxl = real(wxgll(i), kind=CUSTOM_REAL)
              variable_DG(iglob) = variable_DG(iglob) + temp_source * wxl * wzl * jacobianl
              !WRITE(*,*) ">>>>", temp_source, wxl * wzl * jacobianl ! DEBUG
            enddo
          enddo
        endif
        ! TODO: Implement the case source_type = 2.
      endif
    endif
  enddo

  end subroutine compute_add_sources_acoustic_DG_spread


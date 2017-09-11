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


  subroutine write_color_image_snaphot()

#ifdef USE_MPI
  use mpi
#endif

  use constants,only: IMAIN,NGLLX,NGLLZ,REMOVE_PMLS_FROM_JPEG_IMAGES

  use specfem_par,only: myrank,nspec,it,NPROC, &
                        assign_external_model,ibool,kmato,density,rhoext,P_SV, &
                        potential_acoustic,potential_dot_acoustic,potential_dot_dot_acoustic, &
                        potential_gravito,potential_dot_gravito,potential_dot_dot_gravito, &
                        potential_gravitoacoustic,potential_dot_gravitoacoustic,potential_dot_dot_gravitoacoustic, &
                        displ_elastic,veloc_elastic,accel_elastic, &
                        displs_poroelastic,velocs_poroelastic,accels_poroelastic,&! density_p, &
                        E_DG,rhovz_DG,rhovx_DG,rho_DG, nglob_DG, gamma_euler, &
                        rho_DG, rhovz_DG, E_DG, rhovx_DG, cnu, T_init,p_DG_init,gammaext_DG, T_init, E_init, &
                        ispec_is_acoustic_DG, nglob, any_acoustic_DG!, this_iglob_is_acous, ispec_is_acoustic,b_rhovz_DG

  ! PML arrays
  use specfem_par,only: PML_BOUNDARY_CONDITIONS,ispec_is_PML,CONSTRAIN_HYDROSTATIC

  use specfem_par_movie

  implicit none

  !local variables
  integer :: i,j,k,ispec,iglob,iproc
  integer :: ier
  double precision :: rhol, coef
  
  
  logical, dimension(nglob) :: this_iglob_is_acous
  real(kind=CUSTOM_REAL), dimension(nglob_DG) :: vector_DG_temp 
  
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'Creating color image of size ',NX_IMAGE_color,' x ',NZ_IMAGE_color,' for time step ',it
    call flush_IMAIN()
  endif

  if (imagetype_JPEG >= 1 .and. imagetype_JPEG <= 3) then
    if (myrank == 0) write(IMAIN,*) 'drawing scalar image of part of the displacement vector...'!, maxval(rhoveloc_acoustic(2,:))
    
    if(any_acoustic_DG) then
    if(imagetype_JPEG == 1) vector_DG_temp = rhovx_DG/rho_DG
    if(imagetype_JPEG == 2) vector_DG_temp = rhovz_DG/rho_DG
    if(imagetype_JPEG == 3) vector_DG_temp = sqrt((rhovx_DG/rho_DG)**2 + (rhovz_DG/rho_DG)**2)
    !vector_DG_temp = (((gammaext_DG - 1.)*( E_DG &
    !                                 - (0.5)*rho_DG*( (rhovz_DG/rho_DG)**2 + (rhovx_DG/rho_DG)**2 ) )) - p_DG_init)
    endif
    
    call compute_vector_whole_medium(potential_acoustic,potential_gravitoacoustic, &
                                     potential_gravito,veloc_elastic,displs_poroelastic, &
                                     vector_DG_temp)
                                     
  else if (imagetype_JPEG >= 4 .and. imagetype_JPEG <= 6) then
  
    coef = 0.
    if(CONSTRAIN_HYDROSTATIC) coef = 1.  
  
    if(any_acoustic_DG) then
    if(imagetype_JPEG == 4) vector_DG_temp = (((gammaext_DG - 1.)*( E_DG &
                                     - (0.5)*rho_DG*( (rhovz_DG/rho_DG)**2 + (rhovx_DG/rho_DG)**2 ) )) - coef*p_DG_init)
    if(imagetype_JPEG == 5) vector_DG_temp = E_DG - coef*E_init
    if(imagetype_JPEG == 6) vector_DG_temp = ((E_DG/rho_DG - 0.5*((rhovx_DG/rho_DG)**2 + (rhovz_DG/rho_DG)**2))/(cnu) - coef*T_init)
    endif
    
    !WRITE(*,*) "TEST", imagetype_JPEG, coef, CONSTRAIN_HYDROSTATIC
  
    if (myrank == 0) write(IMAIN,*) 'drawing scalar image of part of the velocity vector...'
    call compute_vector_whole_medium(potential_dot_acoustic,potential_dot_gravitoacoustic, &
                                     potential_dot_gravito,veloc_elastic,velocs_poroelastic, &
                                     vector_DG_temp )

  else if (imagetype_JPEG >= 7 .and. imagetype_JPEG <= 9) then
  
    if(any_acoustic_DG) then
    if(imagetype_JPEG == 7) vector_DG_temp = rho_DG 
    if(imagetype_JPEG == 8) vector_DG_temp = rhovx_DG/sqrt(rho_DG)
    if(imagetype_JPEG == 9) vector_DG_temp = rhovz_DG/sqrt(rho_DG)
    endif
  
    if (myrank == 0) write(IMAIN,*) 'drawing scalar image of part of the acceleration vector...'
    call compute_vector_whole_medium(potential_dot_dot_acoustic,potential_dot_dot_gravitoacoustic, &
                                     potential_dot_dot_gravito,accel_elastic,accels_poroelastic, &
                                     vector_DG_temp)

  else if (imagetype_JPEG >= 11 .and. imagetype_JPEG <= 13) then
    ! allocation for normalized representation in JPEG image
    ! for an atmosphere model
    if (myrank == 0) write(IMAIN,*) 'drawing scalar image of part of normalized displacement vector...'
    call compute_vector_whole_medium(potential_acoustic,potential_gravitoacoustic, &
                                     potential_gravito,displ_elastic,displs_poroelastic, &
                                     rho_DG)

    do ispec = 1,nspec
      do j = 1,NGLLZ
        do i = 1,NGLLX
          if (assign_external_model) then
            rhol = rhoext(i,j,ispec)
          else
            rhol = density(1,kmato(ispec))
          endif
          iglob = ibool(i,j,ispec)
          vector_field_display(1,iglob) = sqrt(rhol) * vector_field_display(1,iglob)
          vector_field_display(2,iglob) = sqrt(rhol) * vector_field_display(2,iglob)
        enddo
      enddo
    enddo

  else if (imagetype_JPEG >= 14 .and. imagetype_JPEG <= 16) then
    ! allocation for normalized representation in JPEG image
    ! for an atmosphere model
    call compute_vector_whole_medium(potential_dot_acoustic,potential_dot_gravitoacoustic, &
                                     potential_dot_gravito,veloc_elastic,velocs_poroelastic, &
                                     rho_DG)

    do ispec = 1,nspec
      do j = 1,NGLLZ
        do i = 1,NGLLX
          if (assign_external_model) then
            rhol = rhoext(i,j,ispec)
          else
            rhol = density(1,kmato(ispec))
          endif
          iglob = ibool(i,j,ispec)
          vector_field_display(1,iglob) = sqrt(rhol) * vector_field_display(1,iglob)
          vector_field_display(2,iglob) = sqrt(rhol) * vector_field_display(2,iglob)
        enddo
      enddo
    enddo

  else if (imagetype_JPEG == 10 .and. P_SV) then
    if (myrank == 0) write(IMAIN,*) 'drawing image of pressure field...'
    call compute_pressure_whole_medium()

  else if (imagetype_JPEG == 10 .and. .not. P_SV) then
    call exit_MPI(myrank,'cannot draw pressure field for SH (membrane) waves')

  else
    call exit_MPI(myrank,'wrong type for JPEG snapshots')
  endif

!! DK DK quick hack to remove the PMLs from JPEG images if needed: set the vector field to zero there
  if (PML_BOUNDARY_CONDITIONS .and. REMOVE_PMLS_FROM_JPEG_IMAGES) then
    do ispec = 1,nspec
      if (ispec_is_PML(ispec)) then
        do j = 1,NGLLZ
          do i = 1,NGLLX
            iglob = ibool(i,j,ispec)
            vector_field_display(1,iglob) = 0.d0
            vector_field_display(2,iglob) = 0.d0
          enddo
        enddo
      endif
    enddo
  endif
  
  ! CREATE FLAG FOR NORMALIZATION
  !if(it == 1) then
  this_iglob_is_acous(:) = .false.
  if(any_acoustic_DG) then
  do ispec = 1,nspec
        do j = 1,NGLLZ
          do i = 1,NGLLX
            iglob = ibool(i,j,ispec)
            !if (ispec_is_acoustic_DG(ispec)) then
                this_iglob_is_acous(iglob) = ispec_is_acoustic_DG(ispec)!.true.
           ! endif
          enddo
        enddo
    enddo
  endif
  
!! DK DK quick hack to remove the PMLs from JPEG images if needed

  image_color_data(:,:) = 0.d0

  do k = 1, nb_pixel_loc
    j = ceiling(real(num_pixel_loc(k)) / real(NX_IMAGE_color))
    i = num_pixel_loc(k) - (j-1)*NX_IMAGE_color

    ! avoid edge effects
    if (i < 1 ) i = 1
    if (j < 1 ) j = 1

    if (i > NX_IMAGE_color ) i = NX_IMAGE_color
    if (j > NZ_IMAGE_color ) j = NZ_IMAGE_color

    if (P_SV) then ! P-SH waves, plot a component of vector, its norm, or else pressure
      if (iglob_image_color(i,j) /= -1) then
        ! To normalize independently
        this_ij_image_acous(i,j) = this_iglob_is_acous(iglob_image_color(i,j))
        if (imagetype_JPEG == 1  .or. imagetype_JPEG == 4 .or. imagetype_JPEG == 7 .or. &
            imagetype_JPEG == 11 .or. imagetype_JPEG == 14) then
          ! draw the X component of the vector
          image_color_data(i,j) = vector_field_display(1,iglob_image_color(i,j))

        else if (imagetype_JPEG == 2 .or. imagetype_JPEG == 5 .or. imagetype_JPEG == 8 .or. &
                imagetype_JPEG == 12 .or. imagetype_JPEG == 15) then
          ! draw the Z component of the vector
          image_color_data(i,j) = vector_field_display(2,iglob_image_color(i,j))

        else if (imagetype_JPEG == 3 .or. imagetype_JPEG == 6 .or. imagetype_JPEG == 9 .or. &
                imagetype_JPEG == 13 .or. imagetype_JPEG == 16) then
          ! draw the norm of the vector
          image_color_data(i,j) = sqrt(vector_field_display(1,iglob_image_color(i,j))**2  &
                                     + vector_field_display(2,iglob_image_color(i,j))**2)

        else if (imagetype_JPEG == 10) then
          ! by convention we have stored pressure in the 2. component of the array
          image_color_data(i,j) = vector_field_display(2,iglob_image_color(i,j))

        else
          call exit_MPI(myrank,'wrong type for JPEG snapshots')
        endif
      endif

    else
      ! SH (membrane) waves, plot y-component
      if (iglob_image_color(i,j) /= -1) image_color_data(i,j) = vector_field_display(1,iglob_image_color(i,j))
    endif
  enddo
  !stop
  ! assembling array image_color_data on process zero for color output
#ifdef USE_MPI
  if (NPROC > 1) then
    if (myrank == 0) then
      do iproc = 1, NPROC-1
        call MPI_RECV(data_pixel_recv(1),nb_pixel_per_proc(iproc+1), MPI_DOUBLE_PRECISION, &
                      iproc, 43, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier)

        call MPI_RECV(data_pixel_recv_ij(1),nb_pixel_per_proc(iproc+1), MPI_LOGICAL, &
                      iproc, 44, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier)

        do k = 1, nb_pixel_per_proc(iproc+1)
          j = ceiling(real(num_pixel_recv(k,iproc+1)) / real(NX_IMAGE_color))
          i = num_pixel_recv(k,iproc+1) - (j-1)*NX_IMAGE_color

          ! avoid edge effects
          if (i < 1) i = 1
          if (j < 1) j = 1

          if (i > NX_IMAGE_color) i = NX_IMAGE_color
          if (j > NZ_IMAGE_color) j = NZ_IMAGE_color

          image_color_data(i,j) = data_pixel_recv(k)
          this_ij_image_acous(i,j) = data_pixel_recv_ij(k)
        enddo
      enddo
    else
      do k = 1, nb_pixel_loc
        j = ceiling(real(num_pixel_loc(k)) / real(NX_IMAGE_color))
        i = num_pixel_loc(k) - (j-1)*NX_IMAGE_color

        ! avoid edge effects
        if (i < 1) i = 1
        if (j < 1) j = 1

        if (i > NX_IMAGE_color) i = NX_IMAGE_color
        if (j > NZ_IMAGE_color) j = NZ_IMAGE_color

        if (P_SV) then ! P-SH waves, plot a component of vector, its norm, or else pressure
          data_pixel_send_ij(k) = this_iglob_is_acous(iglob_image_color(i,j))
          if (imagetype_JPEG == 1 .or. imagetype_JPEG == 4 .or. imagetype_JPEG == 7 .or. &
              imagetype_JPEG == 11 .or. imagetype_JPEG == 14) then
             data_pixel_send(k) = vector_field_display(1,iglob_image_color(i,j))  ! draw the X component of the vector

          else if (imagetype_JPEG == 2 .or. imagetype_JPEG == 5 .or. imagetype_JPEG == 8 .or. &
                  imagetype_JPEG == 12 .or. imagetype_JPEG == 15) then
             data_pixel_send(k) = vector_field_display(2,iglob_image_color(i,j))  ! draw the Z component of the vector

          else if (imagetype_JPEG == 3 .or. imagetype_JPEG == 6 .or. imagetype_JPEG == 9 .or. &
                  imagetype_JPEG == 13 .or. imagetype_JPEG == 16) then
            data_pixel_send(k) = sqrt(vector_field_display(1,iglob_image_color(i,j))**2 + &
                                      vector_field_display(2,iglob_image_color(i,j))**2)  ! draw the norm of the vector

          else if (imagetype_JPEG == 10) then
            ! by convention we have stored pressure in the 2. component of the array
            data_pixel_send(k) = vector_field_display(2,iglob_image_color(i,j))

          else
            call exit_MPI(myrank,'wrong type for JPEG snapshots')
          endif

        else ! SH (membrane) waves, plot y-component
          if (iglob_image_color(i,j) /= -1) data_pixel_send(k) = vector_field_display(1,iglob_image_color(i,j))
        endif
      enddo
      
      call MPI_SEND(data_pixel_send(1),nb_pixel_loc,MPI_DOUBLE_PRECISION, 0, 43, MPI_COMM_WORLD, ier)
      call MPI_SEND(data_pixel_send_ij(1),nb_pixel_loc,MPI_LOGICAL, 0, 44, MPI_COMM_WORLD, ier)
      
    endif
  endif
#else
  ! dummy to avoid compiler warning
  ier = 0
  iproc = NPROC
#endif

  ! creates image
  if (myrank == 0) then
    call create_color_image()

    ! user output
    write(IMAIN,*) 'Color image created'
    call flush_IMAIN()
  endif

  end subroutine write_color_image_snaphot


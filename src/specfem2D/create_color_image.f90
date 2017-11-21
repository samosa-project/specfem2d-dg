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

  subroutine create_color_image()

! display a given field as a red and blue color JPEG image

! to display the snapshots : display image*.jpg

  use constants,only: TINYVAL,TINYVAL_DG,HUGEVAL,STABILITY_THRESHOLD

  use specfem_par, only: myrank,it,NSOURCES,P_SV,nrec

  use specfem_par_movie,only: image_color_data,iglob_image_color,NX_IMAGE_color,NZ_IMAGE_color, &
    isnapshot_number,cutsnaps,image_color_vp_display, &
    USE_SNAPSHOT_NUMBER_IN_FILENAME,POWER_DISPLAY_COLOR, &
    DRAW_SOURCES_AND_RECEIVERS, &
    ix_image_color_source,iy_image_color_source,ix_image_color_receiver,iy_image_color_receiver, &
    USE_CONSTANT_MAX_AMPLITUDE,CONSTANT_MAX_AMPLITUDE_TO_USE, this_ij_image_acous

  implicit none

  ! Local variables.
  integer :: i

  ! for the JPEG library
  character(len=1), dimension(3,NX_IMAGE_color,NZ_IMAGE_color) :: JPEG_raw_image
  integer :: ix, iy, R, G, B
  double precision :: amplitude_max, normalized_value, vpmin, vpmax, x1
  character(len=100) :: filename
  logical :: do_warning
  
  double precision, dimension(2) :: amplitude_max_var
  integer, dimension(2,2) :: save_pos

  ! size of cross and square in pixels drawn to represent the source and the receivers in JPEG pictures
  integer :: half_width_cross, thickness_cross, half_size_square

  ! Make the size of the source and receiver symbols depend on the size of the picture using a rule of thumb.
  thickness_cross = 1
  if (NX_IMAGE_color > 2000 .or. NZ_IMAGE_color > 2000) then
    half_width_cross = 6
    half_size_square = 4
  else if (NX_IMAGE_color <= 100 .or. NZ_IMAGE_color <= 100) then
    half_width_cross = 2
    half_size_square = 1
  else if (NX_IMAGE_color <= 250 .or. NZ_IMAGE_color <= 250) then
    half_width_cross = 3
    half_size_square = 2
  else
    half_width_cross = 5
    half_size_square = 3
  endif

! open the image file
! slightly change the beginning of the file name depending if we use the time step of the image number, to avoid confusion
  if (USE_SNAPSHOT_NUMBER_IN_FILENAME) then
    isnapshot_number = isnapshot_number + 1
    write(filename,"('OUTPUT_FILES/img',i7.7,'.jpg')") isnapshot_number
  else
    write(filename,"('OUTPUT_FILES/image',i7.7,'.jpg')") it
  endif

  ! Compute maximum amplitude.
  if (.not. USE_CONSTANT_MAX_AMPLITUDE) then
    !amplitude_max = maxval(abs(image_color_data))
    amplitude_max_var(:) = 0
    ! Different max amplitudes in function of medium*
    do iy = 1, NZ_IMAGE_color
      do ix = 1, NX_IMAGE_color
        if(this_ij_image_acous(ix, iy)) then
          if(abs(image_color_data(ix, iy)) > amplitude_max_var(1)) then
            amplitude_max_var(1) = abs(image_color_data(ix, iy))
            save_pos(1, 1) = ix
            save_pos(1, 2) = iy
          endif
        else
          if(abs(image_color_data(ix, iy)) > amplitude_max_var(2)) then
            amplitude_max_var(2) = abs(image_color_data(ix, iy))
            save_pos(2, 1) = ix
            save_pos(2, 2) = iy
          endif
        endif
        !WRITE(*,*) "-------->", ix, iy,amplitude_max_var, save_pos(1,:), save_pos(2,:) ! DEBUG
      enddo
    enddo
    !stop
  else
    amplitude_max = CONSTANT_MAX_AMPLITUDE_TO_USE
    ! In case of a pre-defined and constant maximum, truncate all values that are outside that constant range.
    where(image_color_data > +CONSTANT_MAX_AMPLITUDE_TO_USE) image_color_data = +CONSTANT_MAX_AMPLITUDE_TO_USE
    where(image_color_data < -CONSTANT_MAX_AMPLITUDE_TO_USE) image_color_data = -CONSTANT_MAX_AMPLITUDE_TO_USE
  endif

  ! checks value (isNaN)
  if (amplitude_max > STABILITY_THRESHOLD .or. amplitude_max < 0 .or. amplitude_max /= amplitude_max) then
    print *,'Warning: failed creating color image, maximum value of amplitude in image color is invalid'
    print *,'amplitude max = ',amplitude_max,' with threshold at ', STABILITY_THRESHOLD
    print *,'Please check your simulation setup...'
    call exit_MPI(myrank,'code became unstable and blew up (image_color_data)')
  endif
  do_warning = .false.

  vpmin = HUGEVAL
  vpmax = TINYVAL
  do iy = 1, NZ_IMAGE_color
    do ix = 1, NX_IMAGE_color
! negative values in image_color_vp_display are a flag indicating a water layer to color in light blue later
      if (iglob_image_color(ix, iy) > -1 .and. image_color_vp_display(ix, iy) >= 0) then
        vpmin = min(vpmin, image_color_vp_display(ix, iy))
        vpmax = max(vpmax, image_color_vp_display(ix, iy))
      endif
    enddo
  enddo

! in the image format, the image starts in the upper-left corner
  do iy = NZ_IMAGE_color, 1, -1
    do ix = 1, NX_IMAGE_color
      if(this_ij_image_acous(ix, iy)) then
        amplitude_max = amplitude_max_var(1)
        ! Trick to avoid very small amplitudes
        if(abs(image_color_data(ix, iy)) < TINYVAL_DG) image_color_data(ix, iy) = 0.
        !WRITE(*,*) "acous ix/iy", ix, iy, amplitude_max
      else
        amplitude_max = amplitude_max_var(2)
       ! WRITE(*,*) "elas ix/iy", ix, iy, amplitude_max
       ! Trick to avoid very small amplitudes
       if(abs(image_color_data(ix, iy)) < TINYVAL) image_color_data(ix, iy) = 0.
      endif

! check if pixel is defined or not (can be above topography for instance)
      if (iglob_image_color(ix, iy) == -1) then

        ! White to display undefined region above topography to avoid visual confusion with a water layer.
        R = 255
        G = 255
        B = 255

      else if (abs(image_color_data(ix, iy)) < amplitude_max * cutsnaps) then
      ! Suppress small amplitudes considered as noise and display the background velocity model instead.
! use P velocity model as background where amplitude is negligible
        if ((P_SV) .and. ((vpmax-vpmin)/max(vpmin, TINYVAL) > 0.02d0)) then
          x1 = (image_color_vp_display(ix, iy)-vpmin)/(vpmax-vpmin)
        else
          x1 = 0.5d0
        endif
! rescale to avoid very dark gray levels
        x1 = x1*0.7 + 0.2
        if (x1 > 1.d0) x1=1.d0
! invert scale: white = vpmin, dark gray = vpmax
        x1 = 1.d0 - x1
! map to [0,255]
        x1 = x1 * 255.d0
        R = nint(x1)
        if (R < 0) R = 0
        if (R > 255) R = 255
        G = R
        B = R
! negative values in image_color_vp_display are a flag indicating a water layer to color in light blue
        if (image_color_vp_display(ix, iy) < 0) then
! use light blue to display water
!!!!!!    R = 204
!!!!!!    G = 255
!!!!!!    B = 255
          R = 135 !!! LightSkyBlue
          G = 206
          B = 250
        endif

      else

! define normalized image data in [-1:1] and convert to nearest integer
! keeping in mind that data values can be negative
        if (amplitude_max >= TINYVAL) then
          normalized_value = image_color_data(ix, iy) / amplitude_max
        else
          normalized_value = image_color_data(ix, iy) / TINYVAL
        endif
        ! check value (isNaN)
        if (normalized_value /= normalized_value) then
          ! will be set to zero
          normalized_value = 0.d0
          do_warning = .true.
        endif

        ! suppress values outside of [-1:+1]
        if (normalized_value < -1.d0) normalized_value = -1.d0
        if (normalized_value > 1.d0) normalized_value = 1.d0
        
        ! TODO: Why was this block here? amplitude_max is already defined above.
        !if(this_ij_image_acous(ix, iy)) then
        !  amplitude_max = amplitude_max_var(1)
        !  !WRITE(*,*) "acous ix/iy", ix, iy, amplitude_max ! DEBUG
        !else
        !  amplitude_max = amplitude_max_var(2)
        !  !WRITE(*,*) "elas ix/iy", ix, iy, amplitude_max ! DEBUG
        !endif

        ! Red if positive value, blue if negative value, no green whatsoever.
        if (normalized_value >= 0.d0) then
          R = nint(255.d0*normalized_value**POWER_DISPLAY_COLOR)
          G = 0
          B = 0
        else
          R = 0
          G = 0
          B = nint(255.d0*abs(normalized_value)**POWER_DISPLAY_COLOR)
        endif

      endif
      ! For JPEG.
      JPEG_raw_image(1, ix, NZ_IMAGE_color-iy+1) = char(R)
      JPEG_raw_image(2, ix, NZ_IMAGE_color-iy+1) = char(G)
      JPEG_raw_image(3, ix, NZ_IMAGE_color-iy+1) = char(B)
    enddo
  enddo

  if (do_warning) then
    print *,'Warning: normalized value is invalid, likely due to unstable simulation! process rank is ',myrank
    print *,'Warning: normalized_value = ',normalized_value,' will be set to zero'
  endif
  
  if (DRAW_SOURCES_AND_RECEIVERS) then
  ! Draw position of the sources and receivers.
! draw position of the sources with orange crosses
    do i = 1, NSOURCES
! avoid edge effects for source or receiver symbols that can be partly outside of the image
      do iy = max(iy_image_color_source(i) - half_width_cross,1), min(iy_image_color_source(i) + half_width_cross, NZ_IMAGE_color)
        do ix = max(ix_image_color_source(i) - thickness_cross,1), min(ix_image_color_source(i) + thickness_cross, NX_IMAGE_color)
          ! Orange.
          R = 255
          G = 157
          B = 0
          ! For JPEG.
          JPEG_raw_image(1, ix, NZ_IMAGE_color-iy+1) = char(R)
          JPEG_raw_image(2, ix, NZ_IMAGE_color-iy+1) = char(G)
          JPEG_raw_image(3, ix, NZ_IMAGE_color-iy+1) = char(B)
        enddo
      enddo
! avoid edge effects for source or receiver symbols that can be partly outside of the image
      do iy = max(iy_image_color_source(i) - thickness_cross,1), min(iy_image_color_source(i) + thickness_cross, NZ_IMAGE_color)
        do ix = max(ix_image_color_source(i) - half_width_cross,1), min(ix_image_color_source(i) + half_width_cross, NX_IMAGE_color)
          ! Orange.
          R = 255
          G = 157
          B = 0
          ! For JPEG.
          JPEG_raw_image(1, ix, NZ_IMAGE_color-iy+1) = char(R)
          JPEG_raw_image(2, ix, NZ_IMAGE_color-iy+1) = char(G)
          JPEG_raw_image(3, ix, NZ_IMAGE_color-iy+1) = char(B)
        enddo
      enddo
    enddo

! draw position of the receivers with green squares
    do i = 1, nrec
! avoid edge effects for source or receiver symbols that can be partly outside of the image
      do iy = max(iy_image_color_receiver(i) - half_size_square,1), &
               min(iy_image_color_receiver(i) + half_size_square,NZ_IMAGE_color)
        do ix = max(ix_image_color_receiver(i) - half_size_square,1), &
                 min(ix_image_color_receiver(i) + half_size_square,NX_IMAGE_color)
          ! Dark green.
          R = 30
          G = 180
          B = 60
          ! For JPEG.
          JPEG_raw_image(1, ix, NZ_IMAGE_color-iy+1) = char(R)
          JPEG_raw_image(2, ix, NZ_IMAGE_color-iy+1) = char(G)
          JPEG_raw_image(3, ix, NZ_IMAGE_color-iy+1) = char(B)
        enddo
      enddo
    enddo
  endif

  ! For JPEG.
  call write_jpeg_image(JPEG_raw_image,NX_IMAGE_color,NZ_IMAGE_color,filename)

  end subroutine create_color_image


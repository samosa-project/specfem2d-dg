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


  subroutine write_wavefield_dumps()

  use constants,only: IMAIN,SIZE_REAL,NGLLX,NGLLZ

  use specfem_par, only: myrank,nglob,nspec, &
                         ibool,coord,P_SV,it,SIMULATION_TYPE, &
                         potential_acoustic,potential_gravitoacoustic, &
                         potential_gravito,displ_elastic,displs_poroelastic, &
                         potential_dot_acoustic,veloc_elastic,velocs_poroelastic,  &
                         potential_dot_dot_acoustic,accel_elastic,accels_poroelastic, &
                         rho_DG, rhovx_DG, rhovz_DG, E_DG, c_V,T_init, &
                         potential_dphi_dx_DG, USE_DISCONTINUOUS_METHOD!, potential_dphi_dz_DG ! Modification for DG.
  
  use specfem_par_lns, only: USE_LNS, LNS_dv, LNS_drho
  
  use specfem_par_movie,only: this_is_the_first_time_we_dump,mask_ibool,imagetype_wavefield_dumps, &
    use_binary_for_wavefield_dumps,vector_field_display

  implicit none

  !local variables
  integer :: i,j,ispec,iglob,icounter,nb_of_values_to_save
  integer :: ier
  ! name of wavefield snapshot file
  character(len=150) :: wavefield_file

  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'Dumping the wave field to a file for time step ',it
    call flush_IMAIN()
  endif
  !WRITE(*,*) ">>>>>>>>>>> this_is_the_first_time_we_dump", this_is_the_first_time_we_dump
  if (this_is_the_first_time_we_dump) then

    if (.not. allocated(mask_ibool)) allocate(mask_ibool(nglob))
        
! save the grid separately once and for all
    if (use_binary_for_wavefield_dumps) then
      write(wavefield_file,"('OUTPUT_FILES/wavefield_grid_for_dumps_',i3.3,'.bin')") myrank
      open(unit=27,file=wavefield_file,form='unformatted',access='direct',status='unknown', &
           action='write',recl=2*SIZE_REAL,iostat=ier)
      if (ier /= 0) call exit_MPI(myrank,'Error opening file wavefield_grid_for_dumps_**.bin')
    else
      write(wavefield_file,"('OUTPUT_FILES/wavefield_grid_for_dumps_',i3.3,'.txt')") myrank
      open(unit=27,file=wavefield_file,status='unknown',action='write',iostat=ier)
      if (ier /= 0) call exit_MPI(myrank,'Error opening file wavefield_grid_for_dumps_**.txt')
    endif

    icounter = 0
    mask_ibool(:) = .false.
    do ispec = 1,nspec
      do j = 1,NGLLZ
        do i = 1,NGLLX
           iglob = ibool(i,j,ispec)
           if (.not. mask_ibool(iglob)) then
             icounter = icounter + 1
             mask_ibool(iglob) = .true.
             if (use_binary_for_wavefield_dumps) then
               write(27,rec=icounter) sngl(coord(1,iglob)),sngl(coord(2,iglob))
             else
               write(27,'(2e16.6)') coord(1,iglob),coord(2,iglob)
             endif
           endif
        enddo
      enddo
    enddo

    close(27)

    ! save nglob to a file once and for all
    write(wavefield_file,"('OUTPUT_FILES/wavefield_grid_value_of_nglob_',i3.3,'.txt')") myrank
    open(unit=27,file=wavefield_file,status='unknown',action='write',iostat=ier)
    if (ier /= 0) call exit_MPI(myrank,'Error opening file wavefield_grid_value_of_nglob_**.txt')
    write(27,*) icounter
    close(27)
    !if (icounter /= nglob) stop 'error: should have icounter == nglob in wavefield dumps'

    this_is_the_first_time_we_dump = .false.

  endif

  if (imagetype_wavefield_dumps == 1) then
    if(USE_DISCONTINUOUS_METHOD) then
      if(USE_LNS) then
        if (myrank == 0) write(IMAIN,*) 'Dumping the displacement vector in elastic elements, and density in LNS elements...'
        call compute_vector_whole_medium(potential_acoustic,potential_gravitoacoustic, &
                                         potential_gravito,displ_elastic,displs_poroelastic,&
                                         LNS_drho)
      else
        if (myrank == 0) write(IMAIN,*) 'Dumping the displacement vector in elastic elements, and density in DG elements...'
        call compute_vector_whole_medium(potential_acoustic,potential_gravitoacoustic, &
                                         potential_gravito,displ_elastic,displs_poroelastic,&
                                         rho_DG)
      endif
    else
      ! Full classical SPECFEM.
      if (myrank == 0) write(IMAIN,*) 'dumping the displacement vector...'
      ! The DG field in compute_vector_whole_medium serves no purpose.
      ! TODO: Do something here instead of this poor patch. Maybe introduce a dummy variable, or a call to a dedicated method (without the DG argument).
      call compute_vector_whole_medium(potential_acoustic,potential_gravitoacoustic, &
                                       potential_gravito,displ_elastic,displs_poroelastic,&
                                       potential_dphi_dx_DG)
    endif
  else if (imagetype_wavefield_dumps == 2) then
    if (myrank == 0) write(IMAIN,*) 'dumping the velocity vector...'
    !call compute_vector_whole_medium(potential_dot_acoustic,potential_gravitoacoustic, &
    !                                 potential_gravito,veloc_elastic,velocs_poroelastic)
    ! Previous call prevents ifort compilation (but, strangely, does not bother gfortran compilation). Thus, we make the following call instead.
    if(USE_DISCONTINUOUS_METHOD) then
      if(USE_LNS) then
        ! Send vx at least.
        ! LM: This is not a very complete method, since compute_vector_whole_medium has the capability of sending a vector to vector_field_display. However for DG, compute_vector_whole_medium is implemented in such a way that this becomes impossible, and duplicates a scalar value over the whole vector attributed variable. This is convoluted, and should be modified. I do not have the time to do it right now.
        call compute_vector_whole_medium(potential_dot_acoustic,potential_gravitoacoustic, &
                                         potential_gravito,veloc_elastic,velocs_poroelastic, &
                                         LNS_dv(1, :))
      else
        ! Send vx at least.
        ! LM: This is not a very complete method, since compute_vector_whole_medium has the capability of sending a vector to vector_field_display. However for DG, compute_vector_whole_medium is implemented in such a way that this becomes impossible, and duplicates a scalar value over the whole vector attributed variable. This is convoluted, and should be modified. I do not have the time to do it right now.
        call compute_vector_whole_medium(potential_dot_acoustic,potential_gravitoacoustic, &
                                         potential_gravito,veloc_elastic,velocs_poroelastic, &
                                         rhovx_DG/rho_DG)
      endif
    else
      ! The DG field in compute_vector_whole_medium serves no purpose.
      ! TODO: Do something here instead of this poor patch. Maybe introduce a dummy variable, or a call to a dedicated method (without the DG argument).
      call compute_vector_whole_medium(potential_dot_acoustic,potential_gravitoacoustic, &
                                       potential_gravito,veloc_elastic,velocs_poroelastic, &
                                       potential_dphi_dx_DG)
    endif

  else if (imagetype_wavefield_dumps == 3) then
    if (myrank == 0) write(IMAIN,*) 'dumping the acceleration vector...'
    !call compute_vector_whole_medium(potential_dot_dot_acoustic,potential_gravitoacoustic, &
    !                                 potential_gravito,accel_elastic,accels_poroelastic)
    ! Previous call prevents ifort compilation (but, strangely, does not bother gfortran compilation). Thus, we make the following call instead.
    ! TODO: Switches as above (imagetype_wavefield_dumps == 2) for DG.
    ! TODO: Do something here instead of this poor patch. Maybe introduce a dummy variable, or a call to a dedicated method (without the DG argument).
    call compute_vector_whole_medium(potential_dot_dot_acoustic,potential_gravitoacoustic, &
                                     potential_gravito,accel_elastic,accels_poroelastic, &
                                     potential_dphi_dx_DG)

  else if (imagetype_wavefield_dumps == 4 .and. P_SV) then
    if (myrank == 0) write(IMAIN,*) 'dumping the pressure field...'
    call compute_pressure_whole_medium()

  else if (imagetype_wavefield_dumps == 4 .and. .not. P_SV) then
    call exit_MPI(myrank,'cannot dump the pressure field for SH (membrane) waves')
  
  else if (imagetype_wavefield_dumps == 5) then
    if(USE_DISCONTINUOUS_METHOD) then
      if (myrank == 0) write(IMAIN,*) 'Dumping the displacement vector in elastic elements, and temperature in DG elements...'
      call compute_vector_whole_medium(potential_acoustic,potential_gravitoacoustic, &
                                       potential_gravito,displ_elastic,displs_poroelastic,&
                                       T_init - (E_DG/rho_DG - 0.5*((rhovx_DG/rho_DG)**2 + (rhovz_DG/rho_DG)**2))/c_V)
    else
      call exit_MPI(myrank,"Can't use imagetype_wavefield_dumps when USE_DISCONTINUOUS_METHOD isn't used.")
    endif

  else
    call exit_MPI(myrank,'wrong type of flag for wavefield dumping')
  endif

  if (use_binary_for_wavefield_dumps) then
    if (P_SV .and. .not. imagetype_wavefield_dumps == 4) then
      nb_of_values_to_save = 2
    else
      nb_of_values_to_save = 1
    endif
    write(wavefield_file,"('OUTPUT_FILES/wavefield',i7.7,'_',i2.2,'_',i3.3,'.bin')") it,SIMULATION_TYPE,myrank
    open(unit=27,file=wavefield_file,form='unformatted',access='direct',status='unknown', &
             action='write',recl=nb_of_values_to_save*SIZE_REAL,iostat=ier)
    if (ier /= 0) call exit_MPI(myrank,'Error opening file wavefield**.bin')
  else
    write(wavefield_file,"('OUTPUT_FILES/wavefield',i7.7,'_',i2.2,'_',i3.3,'.txt')") it,SIMULATION_TYPE,myrank
    open(unit=27,file=wavefield_file,status='unknown',action='write',iostat=ier)
    if (ier /= 0) call exit_MPI(myrank,'Error opening file wavefield**.txt')
  endif

  icounter = 0
  mask_ibool(:) = .false.
  do ispec = 1,nspec
    do j = 1,NGLLZ
      do i = 1,NGLLX
        iglob = ibool(i,j,ispec)
        if (.not. mask_ibool(iglob)) then
          icounter = icounter + 1
          mask_ibool(iglob) = .true.
          if (use_binary_for_wavefield_dumps) then
            if (P_SV) then
              ! P-SV waves
              if (imagetype_wavefield_dumps == 4) then
                ! by convention we use the 2. component of the array to store the pressure above
                write(27,rec=icounter) sngl(vector_field_display(2,iglob))
              else
                write(27,rec=icounter) sngl(vector_field_display(1,iglob)),sngl(vector_field_display(2,iglob))
              endif
            else
              ! SH case
              write(27,rec=icounter) sngl(vector_field_display(1,iglob))
            endif
          else
            if (P_SV) then
              ! P-SV waves
              if (imagetype_wavefield_dumps == 4) then
                ! by convention we use the 2. component of the array to store the pressure above
                write(27,*) sngl(vector_field_display(2,iglob))
              else
                write(27,*) sngl(vector_field_display(1,iglob)),sngl(vector_field_display(2,iglob))
              endif
            else
              ! SH case
              write(27,*) sngl(vector_field_display(1,iglob))
            endif
          endif
        endif
      enddo
    enddo
  enddo

  close(27)

  if (myrank == 0) then
    write(IMAIN,*) 'Wave field dumped'
    call flush_IMAIN()
  endif

  end subroutine write_wavefield_dumps


! ------------------------------------------------------------ !
! DG_WholeDump                                                 !
! ------------------------------------------------------------ !
! Routine to dump all interesting DG wavefields at once.
! Does not account for elastic wavefields.

subroutine DG_WholeDump()
  use constants,only: IMAIN,SIZE_REAL,NGLLX,NGLLZ,CUSTOM_REAL,NDIM
  use specfem_par, only: myrank,nglob, nglob_DG,nspec, &
                         ibool,coord,it,USE_DISCONTINUOUS_METHOD, &
                         rho_DG, rhovx_DG, rhovz_DG, E_DG, gammaext_DG,ispec_is_acoustic_DG!, c_V,T_init, &
!                         potential_dphi_dx_DG, USE_DISCONTINUOUS_METHOD!, potential_dphi_dz_DG ! Modification for DG.
  
  use specfem_par_movie,only: this_is_the_first_time_we_dump,mask_ibool,imagetype_wavefield_dumps, &
    use_binary_for_wavefield_dumps
  
  implicit none
  ! Input/output.
  ! Local variables.
  integer :: i,j,ispec,iglob,icounter,ier
  character(len=150) :: wavefield_file
  character(len=10), parameter :: wvflddmp_tag_rho = 'rho'
  character(len=10), parameter :: wvflddmp_tag_vel = 'vel'
  character(len=10), parameter :: wvflddmp_tag_pre = 'pre'
  integer, parameter :: THEUNIT = 301
  real(kind=CUSTOM_REAL), dimension(NDIM, nglob_DG) :: DG_VEL
  
  if(.not. USE_DISCONTINUOUS_METHOD) then
    if(myrank==0) then
      write(*,*) "********************************"
      write(*,*) "*            ERROR             *"
      write(*,*) "********************************"
      write(*,*) "* The subroutine               *"
      write(*,*) "* DG_WholeDump, "
      write(*,*) "* should only be called when   *"
      write(*,*) "* USE_DISCONTINUOUS_METHOD is  *"
      write(*,*) "* set to .true..               *"
      write(*,*) "********************************"
    endif
    call exit_MPI(myrank,'ERROR.')
  endif
  if(.not.(imagetype_wavefield_dumps==10)) then
    if(myrank==0) then
      write(*,*) "********************************"
      write(*,*) "*            ERROR             *"
      write(*,*) "********************************"
      write(*,*) "* The subroutine               *"
      write(*,*) "* DG_WholeDump, "
      write(*,*) "* should only be called when   *"
      write(*,*) "* imagetype_wavefield_dumps is *"
      write(*,*) "* set to 10.                   *"
      write(*,*) "********************************"
    endif
    call exit_MPI(myrank,'ERROR.')
  endif
  
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'Dumping all wavefields to files, for time step ',it,'.'
    call flush_IMAIN()
  endif
  
  ! This following 'if' is basically the same as in 'write_wavefield_dumps' above, and should not change. Did not factorise code to prevent modifying 'write_wavefield_dumps' too much.
  if(this_is_the_first_time_we_dump) then
    if(.not. allocated(mask_ibool)) allocate(mask_ibool(nglob))
    ! Save the grid, separately, once and for all.
    if(use_binary_for_wavefield_dumps) then
      write(wavefield_file,"('OUTPUT_FILES/wavefield_grid_for_dumps_',i3.3,'.bin')") myrank
      open(unit=THEUNIT,file=wavefield_file,form='unformatted',access='direct',status='unknown', &
           action='write',recl=2*SIZE_REAL,iostat=ier)
      if(ier /= 0) call exit_MPI(myrank,'Error opening file wavefield_grid_for_dumps_**.bin')
    else
      write(wavefield_file,"('OUTPUT_FILES/wavefield_grid_for_dumps_',i3.3,'.txt')") myrank
      open(unit=THEUNIT,file=wavefield_file,status='unknown',action='write',iostat=ier)
      if(ier /= 0) call exit_MPI(myrank,'Error opening file wavefield_grid_for_dumps_**.txt')
    endif
    icounter = 0
    mask_ibool(:) = .false.
    do ispec=1,nspec
      if(ispec_is_acoustic_DG(ispec)) then
        do j=1,NGLLZ; do i=1,NGLLX
          iglob = ibool(i,j,ispec)
          if (.not. mask_ibool(iglob)) then
            icounter = icounter + 1
            mask_ibool(iglob) = .true.
            if(use_binary_for_wavefield_dumps) then
              !write(THEUNIT,rec=icounter) sngl(coord(1,iglob)), sngl(coord(2,iglob))
              write(THEUNIT,rec=icounter) sngl(coord(1:NDIM,iglob))
            else
              write(THEUNIT,'(2e16.6)') coord(1, iglob), coord(2, iglob)
            endif
          endif
        enddo; enddo
      endif ! Endif on ispec_is_acoustic_DG
    enddo
    close(THEUNIT)
    ! Save nglob to a separate file, once and for all.
    write(wavefield_file,"('OUTPUT_FILES/wavefield_grid_value_of_nglob_',i3.3,'.txt')") myrank
    open(unit=THEUNIT,file=wavefield_file,status='unknown',action='write',iostat=ier)
    if (ier /= 0) call exit_MPI(myrank,'Error opening file wavefield_grid_value_of_nglob_**.txt')
    write(THEUNIT,*) icounter
    close(THEUNIT)
    !if (icounter /= nglob) stop 'error: should have icounter == nglob in wavefield dumps'
    this_is_the_first_time_we_dump = .false.
  endif ! Endif on this_is_the_first_time_we_dump.
  
  ! Now starts the real DG_WholeDump implementation.
  if(use_binary_for_wavefield_dumps) then
      ! RHO
      call DG_WholeDump_OneWaveField(1, wvflddmp_tag_rho, rho_DG)
      ! VEL
      DG_VEL(1,:) = rhovx_DG/rho_DG
      DG_VEL(2,:) = rhovz_DG/rho_DG
      call DG_WholeDump_OneWaveField(2, wvflddmp_tag_vel, DG_VEL)
      ! PRE.
      call DG_WholeDump_OneWaveField(1, wvflddmp_tag_pre, &
             (   (gammaext_DG-1.)*E_DG                &
               - 0.5*rho_DG*(   (rhovz_DG/rho_DG)**2   &
                              + (rhovx_DG/rho_DG)**2 )))
    
    if(myrank == 0) then
      write(IMAIN,*) 'Dumping ended.'
      call flush_IMAIN()
    endif
    
  else
    write(*,*) "********************************"
    write(*,*) "*            ERROR             *"
    write(*,*) "********************************"
    write(*,*) "* When using                   *"
    write(*,*) "* DG_WholeDump, "
    write(*,*) "* only binary outputs are      *"
    write(*,*) "* authorised (for now).        *"
    write(*,*) "********************************"
    stop
  endif ! Endif on use_binary_for_wavefield_dumps.
end subroutine DG_WholeDump


! ------------------------------------------------------------ !
! DG_WholeDump_OneWaveField                                    !
! ------------------------------------------------------------ !
! Code facotrisation used by DG_WholeDump: open
! an output file, call DG_WholeDump_WriteToUnit,
! and close the file.

subroutine DG_WholeDump_OneWaveField(nbValuesPerPoint, tag, THEFIELD_DG_1D)
  use constants,only: IMAIN,SIZE_REAL,CUSTOM_REAL
  use specfem_par, only: myrank, it, nglob_DG, SIMULATION_TYPE
  ! Input/output.
  character(len=10), intent(in) :: tag
  integer, intent(in) :: nbValuesPerPoint
  real(kind=CUSTOM_REAL), dimension(nbValuesPerPoint, nglob_DG), intent(in) :: THEFIELD_DG_1D
  ! Local variables.
  integer, parameter :: THEUNIT = 301
  character(len=150) :: wavefield_file
  integer :: ier
  ! Name file.
  write(wavefield_file,"('OUTPUT_FILES/wavefield_',A,'_',i8.8,'_',i2.2,'_',i5.5,'.bin')") &
        trim(tag), it, SIMULATION_TYPE, myrank
  ! Open file.
  open(unit=THEUNIT, file=wavefield_file, form='unformatted', access='direct', status='unknown', &
           action='write', recl=nbValuesPerPoint*SIZE_REAL, iostat=ier)
  if (ier /= 0) call exit_MPI(myrank,'Error opening file for PRE dumping.')
  ! Print to file.
  call DG_WholeDump_WriteToUnit(THEUNIT, nbValuesPerPoint, THEFIELD_DG_1D)
  close(THEUNIT)
  if(myrank == 0) then
    write(IMAIN,*) '> Dumping: ', trim(tag), ' wavefield dumped.'
    call flush_IMAIN()
  endif
end subroutine DG_WholeDump_OneWaveField


! ------------------------------------------------------------ !
! DG_WholeDump_WriteToUnit                                     !
! ------------------------------------------------------------ !
! Code factorisation used by DG_WholeDump_OneWaveField. Dump
! some wavefield with some file already open.

subroutine DG_WholeDump_WriteToUnit(THEUNIT, nbValuesPerPoint, THEFIELD_DG_1D)
  use constants,only: IMAIN,SIZE_REAL,NGLLX,NGLLZ,CUSTOM_REAL,NDIM
  use specfem_par, only: myrank, ibool_DG, nspec, nglob_DG, ispec_is_acoustic_DG, ibool
  use specfem_par_movie,only: mask_ibool
  ! Input/output.
  integer, intent(in) :: THEUNIT
  integer, intent(in) :: nbValuesPerPoint
  real(kind=CUSTOM_REAL), dimension(nbValuesPerPoint, nglob_DG), intent(in) :: THEFIELD_DG_1D
  ! Local variables.
  integer :: i, j, ispec, iglob, iglob_DG, icounter
  icounter = 0
  mask_ibool(:) = .false.
  do ispec=1,nspec
    if(ispec_is_acoustic_DG(ispec)) then
      do j=1,NGLLZ; do i=1,NGLLX
        iglob = ibool(i, j, ispec)
        iglob_DG = ibool_DG(i, j, ispec)
        if (.not. mask_ibool(iglob)) then
          icounter = icounter + 1
          mask_ibool(iglob) = .true.
          if(nbValuesPerPoint==1) then
            write(THEUNIT, rec=icounter) sngl(THEFIELD_DG_1D(1, iglob_DG))
          else if(nbValuesPerPoint==NDIM) then
            write(THEUNIT, rec=icounter) sngl(THEFIELD_DG_1D(1:NDIM, iglob_DG))
          else
            call exit_MPI(myrank,"Error in 'DG_WholeDump_WriteToUnit': wrong nbValuesPerPoint.")
          endif
        endif
      enddo; enddo
    endif ! Endif on ispec_is_acoustic_DG
  enddo
end subroutine DG_WholeDump_WriteToUnit


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
!               Pieyre Le Loher, pieyre DOT le-loher aT inria.fr
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


  subroutine read_external_model()

  use constants,only: CUSTOM_REAL,NGLLX,NGLLZ,TINYVAL

  use specfem_par, only: any_acoustic,any_gravitoacoustic,any_elastic,any_poroelastic, &
    ispec_is_acoustic,ispec_is_gravitoacoustic,ispec_is_elastic,ispec_is_poroelastic,ispec_is_anisotropic, &
    nspec,nglob,ibool, &
    READ_VELOCITIES_AT_f0,inv_tau_sigma_nu1_sent,&
    phi_nu1_sent,inv_tau_sigma_nu2_sent,phi_nu2_sent,Mu_nu1_sent,Mu_nu2_sent, &
    inv_tau_sigma_nu1,inv_tau_sigma_nu2,phi_nu1,phi_nu2,Mu_nu1,Mu_nu2,&
    coord,kmato,rhoext,vpext,vsext,gravityext,Nsqext, &
    QKappa_attenuationext,Qmu_attenuationext, &
    c11ext,c13ext,c15ext,c33ext,c35ext,c55ext,c12ext,c23ext,c25ext, &
    MODEL,ATTENUATION_VISCOELASTIC_SOLID,P_SV,&
    tomo_material,myrank, &
    ! MODIF DG
    windxext, windzext, pext_DG, gammaext_DG, etaext, muext, kappa_DG, USE_DISCONTINUOUS_METHOD,&
    Bxext, Bzext, N0ext, &
    EXTERNAL_DG_ONLY_MODEL_FILENAME, &!ispec_is_acoustic_DG,&
    !TEST
    ibool_before_perio
  use specfem_par_lns, only: BCKGRD_MDL_LNS_is_binary, BCKGRD_MDL_LNS_FILENAME

  implicit none

  ! Local variables
  integer :: i,j,ispec
  integer :: ier
  real(kind=CUSTOM_REAL) :: previous_vsext
  real(kind=CUSTOM_REAL) :: tmp1, tmp2,tmp3
  double precision :: rho_dummy,vp_dummy,vs_dummy,mu_dummy,lambda_dummy,vs_val,vp_val,rho_val
  character(len=150) :: inputname
  integer :: nlines_header, nblines_model
  
  nlines_header = -1
  nblines_model = -1
  
  ! Initialise to prevent uninitialised conditionnal jumps.
  c11ext = 0._CUSTOM_REAL
  c13ext = 0._CUSTOM_REAL
  c15ext = 0._CUSTOM_REAL
  c33ext = 0._CUSTOM_REAL
  c35ext = 0._CUSTOM_REAL
  c55ext = 0._CUSTOM_REAL
  c12ext = 0._CUSTOM_REAL
  c23ext = 0._CUSTOM_REAL
  c25ext = 0._CUSTOM_REAL
  vsext = 0._CUSTOM_REAL
  gravityext = 0._CUSTOM_REAL

  if (tomo_material > 0) MODEL = 'tomo'

  if (trim(MODEL) == 'legacy') then

    write(inputname,'(a,i6.6,a)') 'DATA/proc',myrank,'_model_velocity.dat_input'
    open(unit=1001,file=inputname,status='unknown')
    do ispec = 1,nspec
      do j = 1,NGLLZ
        do i = 1,NGLLX
          read(1001,*) tmp1,tmp2,tmp3,rho_val,vp_val,vs_val
          rhoext(i,j,ispec) = rho_val
          vpext(i,j,ispec) = vp_val
          vsext(i,j,ispec) = vs_val
          QKappa_attenuationext(i,j,ispec) = 9999.d0
          Qmu_attenuationext(i,j,ispec) = 9999.d0
        enddo
      enddo
    enddo
    close(1001)

  else if (trim(MODEL)=='ascii') then
    write(inputname,'(a,i6.6,a)') 'DATA/proc',myrank,'_rho_vp_vs.dat'
    open(unit=1001,file= inputname,status='unknown')
    do ispec = 1,nspec
      do j = 1,NGLLZ
        do i = 1,NGLLX
          read(1001,*) tmp1,tmp2,rho_val,vp_val,vs_val
          rhoext(i,j,ispec) = rho_val
          vpext(i,j,ispec) = vp_val
          vsext(i,j,ispec) = vs_val
          QKappa_attenuationext(i,j,ispec) = 9999.d0
          Qmu_attenuationext(i,j,ispec) = 9999.d0
        enddo
      enddo
    enddo
    close(1001)

  else if ((trim(MODEL) == 'binary') .or. (trim(MODEL) == 'gll')) then
      write(inputname,'(a,i6.6,a)') 'DATA/proc',myrank,'_rho.bin'
      open(unit = 1001, file = inputname, status='old',action='read',form='unformatted', iostat=ier)
      if (ier /= 0) stop 'Error opening rho.bin file.'

      read(1001) rhoext
      close(1001)
      print *, 'rho', minval(rhoext), maxval(rhoext)

      write(inputname,'(a,i6.6,a)') 'DATA/proc',myrank,'_vp.bin'
      open(unit = 1001, file = inputname, status='old',action='read',form='unformatted', iostat=ier)
      if (ier /= 0) stop 'Error opening vp.bin file.'

      read(1001) vpext
      close(1001)

      write(inputname,'(a,i6.6,a)') 'DATA/proc',myrank,'_vs.bin'
      open(unit = 1001, file = inputname, status='old',action='read',form='unformatted', iostat=ier)
      if (ier /= 0) stop 'Error opening vs.bin file.'

      read(1001) vsext
      close(1001)

      QKappa_attenuationext(:,:,:) = 9999.d0
      Qmu_attenuationext(:,:,:) = 9999.d0

  else if (trim(MODEL)=='binary_voigt') then
      write(inputname,'(a,i6.6,a)') 'DATA/proc',myrank,'_rho.bin'
      open(unit = 1001, file = inputname,status='old',action='read',form='unformatted', iostat=ier)
      if (ier /= 0) stop 'Error opening rho.bin file.'

      read(1001) rhoext
      close(1001)
      print *, 'rho', minval(rhoext), maxval(rhoext)

      write(inputname,'(a,i6.6,a)') 'DATA/proc',myrank,'_c11.bin'
      open(unit = 1001, file = inputname,status='old',action='read',form='unformatted', iostat=ier)
      if (ier /= 0) stop 'Error opening c11.bin file.'

      read(1001) c11ext
      close(1001)
      print *, 'c11ext', minval(c11ext), maxval(c11ext)

      write(inputname,'(a,i6.6,a)') 'DATA/proc',myrank,'_c13.bin'
      open(unit = 1001, file = inputname,status='old',action='read',form='unformatted', iostat=ier)
      if (ier /= 0) stop 'Error opening c13.bin file.'

      read(1001) c13ext
      close(1001)
      print *, 'c13ext', minval(c13ext), maxval(c13ext)

      write(inputname,'(a,i6.6,a)') 'DATA/proc',myrank,'_c15.bin'
      open(unit = 1001, file = inputname,status='old',action='read',form='unformatted', iostat=ier)
      if (ier /= 0) stop 'Error opening c15.bin file.'

      read(1001) c15ext
      close(1001)
      print *, 'c15ext', minval(c15ext), maxval(c15ext)

      write(inputname,'(a,i6.6,a)') 'DATA/proc',myrank,'_c33.bin'
      open(unit = 1001, file = inputname,status='old',action='read',form='unformatted', iostat=ier)
      if (ier /= 0) stop 'Error opening c33.bin file.'

      read(1001) c33ext
      close(1001)
      print *, 'c33ext', minval(c33ext), maxval(c33ext)

      write(inputname,'(a,i6.6,a)') 'DATA/proc',myrank,'_c35.bin'
      open(unit = 1001, file = inputname,status='old',action='read',form='unformatted', iostat=ier)
      if (ier /= 0) stop 'Error opening c35.bin file.'

      read(1001) c35ext
      close(1001)
      print *, 'c35ext', minval(c35ext), maxval(c35ext)

      write(inputname,'(a,i6.6,a)') 'DATA/proc',myrank,'_c55.bin'
      open(unit = 1001, file = inputname,status='old',action='read',form='unformatted', iostat=ier)
      if (ier /= 0) stop 'Error opening c55.bin file.'

      read(1001) c55ext
      close(1001)
      print *, 'c55ext', minval(c55ext), maxval(c55ext)

      do ispec = 1,nspec
        do j = 1,NGLLZ
          do i = 1,NGLLX
            if(c55ext(i,j,ispec) == 0.0)then
               c33ext(i,j,ispec) = 0.0
               c12ext(i,j,ispec) = 0.0
               c23ext(i,j,ispec) = 0.0
               vpext(i,j,ispec) = 1500.0
               vsext(i,j,ispec) = 0.0
            else
               c12ext(i,j,ispec) = 1.d-6
               c23ext(i,j,ispec) = 1.d-6
               vpext(i,j,ispec) = sqrt(c33ext(i,j,ispec)/rhoext(i,j,ispec))
               vsext(i,j,ispec) = sqrt(c55ext(i,j,ispec)/rhoext(i,j,ispec))
            endif
            c25ext(i,j,ispec) = 0.0
          enddo
        enddo
      enddo

  else if (trim(MODEL)=='external') then
    !call define_external_model(coord,kmato,ibool,rhoext,vpext,vsext, &
    call define_external_model_atmos_tabular_gravitoacoustic(coord,kmato,ibool,rhoext,vpext,vsext, &
                               QKappa_attenuationext,Qmu_attenuationext,gravityext,Nsqext, &
                               c11ext,c13ext,c15ext,c33ext,c35ext,c55ext,c12ext,c23ext,c25ext,nspec,nglob, &
                               ! MODIF DG
                               windxext, windzext, pext_DG, gammaext_DG, etaext, muext, kappa_DG, Bxext, Bzext, N0ext)
  
  
  else if (trim(MODEL)=='external_DG') then
    call external_dg_check_and_get_nblines(EXTERNAL_DG_ONLY_MODEL_FILENAME, nlines_header, nblines_model)
    call define_external_model_DG_only(nlines_header, nblines_model)
    ! DEBUG PRINT MODEL
    if(.false. .and. myrank==1) then ! DEBUG
      open(unit=504,file='OUTPUT_FILES/TESTMODEL',status='unknown',action='write', position="append")
      do ispec = 1,nspec
        do j = 1,NGLLZ
          do i = 1,NGLLX
            write(504,*) coord(1, ibool_before_perio(i, j, ispec)), coord(2, ibool_before_perio(i, j, ispec)),&
                         pext_DG(i, j, ispec)!,merge(1.,0.,ispec_is_acoustic_DG(ispec))!
          enddo
        enddo
      enddo
      close(504)
      stop
      ! Matlab one-liner plot:
      !a=importdata("/home/l.martire/Documents/SPECFEM/specfem-dg-master/EXAMPLES/test_lns_load_external/OUTPUT_FILES/TESTMODEL");X=a(:,1);Y=a(:,2);V=a(:,3);scatter(X,Y,20,V,'filled'); colorbar
    endif
    
  else if (trim(MODEL)=='LNS_generalised') then
    if(BCKGRD_MDL_LNS_is_binary) then
      call LNS_generalised_binary_grab_nlines(nblines_model)
    else
      call external_dg_check_and_get_nblines(BCKGRD_MDL_LNS_FILENAME, nlines_header, nblines_model)
    endif
    call lns_load_background_model(nlines_header, nblines_model)
!    stop 'kek'

  else if (trim(MODEL)=='tomo') then
    call define_external_model_from_tomo_file()
  
  else
    write(*,*) "********************************"
    write(*,*) "*            ERROR             *"
    write(*,*) "********************************"
    write(*,*) "* Parfile MODEL variable has   *"
    write(*,*) "* unknown value:               *"
    write(*,*) "* ", trim(MODEL)
    write(*,*) "********************************"
    stop
  endif ! Endif on trim(MODEL).

  if (trim(MODEL)=='external' .or. trim(MODEL)=='tomo' .or. trim(MODEL)=='external_DG') then
    ! check that the external model that has just been defined makes sense
    do ispec = 1,nspec
      do j = 1,NGLLZ
        do i = 1,NGLLX

          if (c11ext(i,j,ispec) > TINYVAL .or. c13ext(i,j,ispec) > TINYVAL .or. c15ext(i,j,ispec) > TINYVAL .or. &
              c33ext(i,j,ispec) > TINYVAL .or. c35ext(i,j,ispec) > TINYVAL .or. c55ext(i,j,ispec) > TINYVAL) then
            ! vp, vs : assign dummy values, trick to avoid floating point errors in the case of an anisotropic medium
            vpext(i,j,ispec) = 20.d0
            vsext(i,j,ispec) = 10.d0
          endif

          ! check that the element type is not redefined compared to what is defined initially in DATA/Par_file
          if ((c11ext(i,j,ispec) > TINYVAL .or. c13ext(i,j,ispec) > TINYVAL .or. c15ext(i,j,ispec) > TINYVAL .or. &
               c33ext(i,j,ispec) > TINYVAL .or. c35ext(i,j,ispec) > TINYVAL .or. c55ext(i,j,ispec) > TINYVAL) &
              .and. .not. ispec_is_anisotropic(ispec)) &
            stop 'Error: non anisotropic material in DATA/Par_file or &
                 &external mesh redefined as anisotropic in define_external_model()'

          if (vsext(i,j,ispec) < TINYVAL .and. (ispec_is_elastic(ispec) .or. ispec_is_anisotropic(ispec))) &
            stop 'Error: non acoustic material in DATA/Par_file or &
                 &external mesh redefined as acoustic in define_external_model()'

          if (vsext(i,j,ispec) > TINYVAL .and. .not. ispec_is_elastic(ispec)) &
            stop 'Error: acoustic material in DATA/Par_file or &
                 &external mesh redefined as non acoustic in define_external_model()'

        enddo
      enddo
    enddo
  endif

  if(trim(MODEL)=='external_DG') then
    ! In that case, flag re-assigning is not needed since the DG and elastic parts stay unchanged.
    return
  endif

  ! re-assigns flags
  ! initializes
  any_acoustic = .false.
  any_gravitoacoustic = .false.
  any_elastic = .false.
  any_poroelastic = .false.

  ispec_is_acoustic(:) = .false.
  ispec_is_gravitoacoustic(:) = .false.
  ispec_is_anisotropic(:) = .false.
  ispec_is_elastic(:) = .false.
  ispec_is_poroelastic(:) = .false.

  ! initialize to dummy values
  ! convention to indicate that Q = 9999 in that element i.e. that there is no viscoelasticity in that element
  inv_tau_sigma_nu1(:,:,:,:) = -1._CUSTOM_REAL
  phi_nu1(:,:,:,:) = -1._CUSTOM_REAL
  inv_tau_sigma_nu2(:,:,:,:) = -1._CUSTOM_REAL
  phi_nu2(:,:,:,:) = -1._CUSTOM_REAL
  Mu_nu1(:,:,:) = -1._CUSTOM_REAL
  Mu_nu2(:,:,:) = -1._CUSTOM_REAL

  do ispec = 1,nspec
    previous_vsext = vsext(1,1,ispec)
    do j = 1,NGLLZ
      do i = 1,NGLLX
        !print *,"vsext(i,j,ispec)",vsext(i,j,ispec)
        !print *,"gravityext(i,j,ispec)",gravityext(i,j,ispec)
        if (P_SV .and. (.not. (i == 1 .and. j == 1)) .and. &
          ((vsext(i,j,ispec) >= TINYVAL .and. previous_vsext < TINYVAL) .or. &
           (vsext(i,j,ispec) < TINYVAL  .and. previous_vsext >= TINYVAL)))  &
          call exit_MPI(myrank,'external velocity model cannot be both fluid and solid inside the same spectral element')

        if (c11ext(i,j,ispec) > TINYVAL .or. c13ext(i,j,ispec) > TINYVAL .or. c15ext(i,j,ispec) > TINYVAL .or. &
            c33ext(i,j,ispec) > TINYVAL .or. c35ext(i,j,ispec) > TINYVAL .or. c55ext(i,j,ispec) > TINYVAL) then
          ispec_is_anisotropic(ispec) = .true.
          ispec_is_poroelastic(ispec) = .false.
          ispec_is_elastic(ispec) = .true.
          any_elastic = .true.
          QKappa_attenuationext(i,j,ispec) = 9999.d0
          Qmu_attenuationext(i,j,ispec) = 9999.d0
        else if ((vsext(i,j,ispec) < TINYVAL) .and. (gravityext(i,j,ispec) < TINYVAL)) then
          ispec_is_elastic(ispec) = .false.
          ispec_is_poroelastic(ispec) = .false.
          ispec_is_gravitoacoustic(ispec) = .false.
          ispec_is_acoustic(ispec) = .true.
          any_acoustic = .true.
        else if ((vsext(i,j,ispec) < TINYVAL) .and. (gravityext(i,j,ispec) >= TINYVAL)) then
          ispec_is_elastic(ispec) = .false.
          ispec_is_poroelastic(ispec) = .false.
          ispec_is_acoustic(ispec)=.false.
          ispec_is_gravitoacoustic(ispec) = .true.
          any_gravitoacoustic = .true.
        else
          ispec_is_poroelastic(ispec) = .false.
          ispec_is_elastic(ispec) = .true.
          any_elastic = .true.
        endif
        
        ! DG modif => overwrite an acoustic medium
         if (vsext(i,j,ispec) < TINYVAL .AND. USE_DISCONTINUOUS_METHOD) then
                  ispec_is_elastic(ispec) = .false.
                  ispec_is_poroelastic(ispec) = .false.
                  ispec_is_gravitoacoustic(ispec) = .false.
                  ispec_is_acoustic(ispec) = .true.
                  any_acoustic = .true.
                  any_gravitoacoustic = .false.
         endif

        ! attenuation is not implemented in acoustic (i.e. fluid) media for now, only in viscoelastic (i.e. solid) media
        if (ispec_is_acoustic(ispec)) cycle

        ! check that attenuation values entered by the user make sense
        if ((QKappa_attenuationext(i,j,ispec) <= 9998.999d0 .and. Qmu_attenuationext(i,j,ispec) >  9998.999d0) .or. &
            (QKappa_attenuationext(i,j,ispec) >  9998.999d0 .and. Qmu_attenuationext(i,j,ispec) <= 9998.999d0)) &
          stop 'need to have Qkappa and Qmu both above or both below 9999 for a given material; &
               &trick: use 9998 if you want to turn off one'

        ! if no attenuation in that elastic element
        if (QKappa_attenuationext(i,j,ispec) > 9998.999d0) cycle

        ! attenuation
        call attenuation_model(dble(QKappa_attenuationext(i,j,ispec)),dble(Qmu_attenuationext(i,j,ispec)))

        inv_tau_sigma_nu1(i,j,ispec,:) = inv_tau_sigma_nu1_sent(:)
        phi_nu1(i,j,ispec,:) = phi_nu1_sent(:)
        inv_tau_sigma_nu2(i,j,ispec,:) = inv_tau_sigma_nu2_sent(:)
        phi_nu2(i,j,ispec,:) = phi_nu2_sent(:)
        Mu_nu1(i,j,ispec) = Mu_nu1_sent
        Mu_nu2(i,j,ispec) = Mu_nu2_sent

        if (ATTENUATION_VISCOELASTIC_SOLID .and. READ_VELOCITIES_AT_f0) then
          if (ispec_is_anisotropic(ispec) .or. ispec_is_poroelastic(ispec) .or. ispec_is_gravitoacoustic(ispec)) &
            stop 'READ_VELOCITIES_AT_f0 only implemented for non anisotropic, non poroelastic, &
                  &non gravitoacoustic materials for now'

          vp_dummy = dble(vpext(i,j,ispec))
          vs_dummy = dble(vpext(i,j,ispec))
          rho_dummy = dble(rhoext(i,j,ispec))

          call shift_velocities_from_f0(vp_dummy,vs_dummy,rho_dummy,mu_dummy,lambda_dummy)

        endif

        previous_vsext = vsext(i,j,ispec)

      enddo
    enddo
  enddo

  end subroutine read_external_model


subroutine external_dg_check_and_get_nblines(FILENAME, nlines_header, nblines_model)
  !use specfem_par, only: EXTERNAL_DG_ONLY_MODEL_FILENAME
  implicit none
  ! Input/output.
  character(len=100), intent(in) :: FILENAME
  integer, intent(out) :: nlines_header, nblines_model
  ! Local variables.
  logical fileexists
  ! Check existence of model file.
  fileexists=.false.
  INQUIRE(File=FILENAME, Exist=fileexists)
  fileexists=.true.
  if(.not. fileexists) then
    write(*,*) "********************************"
    write(*,*) "*            ERROR             *"
    write(*,*) "********************************"
    write(*,*) "* External model file          *"
    write(*,*) "* ", trim(FILENAME)
    write(*,*) "* does not exist in folder.    *"
    write(*,*) "* Get the model file in the    *"
    write(*,*) "* simulation folder, or set    *"
    write(*,*) "* MODEL to a different value   *"
    write(*,*) "* in the parameter file.       *"
    write(*,*) "********************************"
    stop
  endif
  call external_model_DG_only_find_nblines(FILENAME, nlines_header, nblines_model)
end subroutine external_dg_check_and_get_nblines

subroutine LNS_generalised_binary_grab_nlines(nblines_model)
  use specfem_par_lns, only: BCKGRD_MDL_LNS_FILENAME, BCKGRD_MDL_LNS_FILENAME_HEADER
  implicit none
  ! Input/output.
  integer(kind=4), intent(out) :: nblines_model
  ! Local variables.
  integer io
  OPEN(100, file=BCKGRD_MDL_LNS_FILENAME_HEADER, iostat=io)
  if (io/=0) stop "Error opening background model header file."
  read(100, *, iostat=io) nblines_model
!  write(*,*) nblines_model
  close(100)
!  stop 'kek'
end subroutine LNS_generalised_binary_grab_nlines

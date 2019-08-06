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

  subroutine check_grid()

! check the mesh, stability and number of points per wavelength

#ifdef USE_MPI
  use mpi
#endif

  use specfem_par
  use specfem_par_movie

  implicit none

! option to display only part of the mesh and not the whole mesh,
! for instance to analyze Cuthill-McKee mesh partitioning etc.

  ! local parameters
  double precision :: vpIImax_local,vpIImin_local
  double precision :: vsmin,vsmax,densmin,densmax,vpImax_local,vpImin_local,vsmin_local

  double precision :: kappa_s,kappa_f,kappa_fr,mu_s,mu_fr,rho_s,rho_f,eta_f,phi,tort,rho_bar
  double precision :: D_biot,H_biot,C_biot,M_biot

  double precision :: cpIloc,cpIIloc,csloc
  double precision :: cpIsquare,cpIIsquare,cssquare
  double precision :: f0max,w_c,perm_xx
  double precision :: denst
  double precision :: lambdaplus2mu,mu
  double precision :: distance_min,distance_max,distance_min_local,distance_max_local
  double precision :: courant_stability_number_max,lambdaPImin,lambdaPImax,lambdaPIImin,lambdaPIImax,lambdaSmin,lambdaSmax
  double precision :: distance_1,distance_2,distance_3,distance_4

! for the stability condition
! maximum polynomial degree for which we can compute the stability condition
  integer, parameter :: NGLLX_MAX_STABILITY = 15
  double precision :: percent_GLL(NGLLX_MAX_STABILITY)

! for slice totals
  double precision  :: vpImin_glob,vpImax_glob,vsmin_glob,vsmax_glob,densmin_glob,densmax_glob
  double precision  :: vpIImin_glob,vpIImax_glob
  double precision  :: distance_min_glob,distance_max_glob
  double precision  :: courant_stability_max_glob,lambdaPImin_glob,lambdaPImax_glob,&
                       lambdaPIImin_glob,lambdaPIImax_glob,lambdaSmin_glob,lambdaSmax_glob, &
                       lambdaPmin_in_fluid_histo_glob,lambdaPmax_in_fluid_histo_glob

  integer :: ier
  integer :: i,j,ispec,material

! for histogram of number of points per wavelength
  logical :: any_fluid_histo,any_fluid_histo_glob
  logical :: create_wavelength_histogram
  double precision :: lambdaPmin_in_fluid_histo,lambdaPmax_in_fluid_histo
  double precision :: lambdaSmin_histo,lambdaSmax_histo

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'Checking mesh and stability'
    call flush_IMAIN()
  endif

! define percentage of smallest distance between GLL points for NGLLX points
! percentages were computed by calling the GLL points routine for each degree
  call check_grid_setup_GLLper(percent_GLL,NGLLX_MAX_STABILITY)


!---- compute parameters for the spectral elements

  vpImin = HUGEVAL
  vpImax = -HUGEVAL

  if (any_elastic .or. any_poroelastic) then
    vsmin = HUGEVAL
    vsmax = -HUGEVAL
  else
    vsmin = 0
    vsmax = 0
  endif

  if (any_poroelastic) then
    vpIImin = HUGEVAL
    vpIImax = -HUGEVAL
  else
    vpIImin = 0
    vpIImax = 0
  endif

  densmin = HUGEVAL
  densmax = -HUGEVAL

  distance_min = HUGEVAL
  distance_max = -HUGEVAL

  courant_stability_number_max = -HUGEVAL

  lambdaPImin = HUGEVAL
  lambdaPImax = -HUGEVAL

  if (any_elastic .or. any_poroelastic) then
    lambdaSmin = HUGEVAL
    lambdaSmax = -HUGEVAL
  else
    lambdaSmin = 0
    lambdaSmax = 0
  endif

  if (any_poroelastic) then
    lambdaPIImin = HUGEVAL
    lambdaPIImax = -HUGEVAL
  else
    lambdaPIImin = 0
    lambdaPIImax = 0
  endif

  lambdaPmin_in_fluid_histo = HUGEVAL
  lambdaPmax_in_fluid_histo = -HUGEVAL

  any_fluid_histo = .false.

  do ispec= 1,nspec

    if (ispec_is_poroelastic(ispec)) then

      ! gets poroelastic material
      call get_poroelastic_material(ispec,phi,tort,mu_s,kappa_s,rho_s,kappa_f,rho_f,eta_f,mu_fr,kappa_fr,rho_bar)
      denst = rho_s

      ! Biot coefficients for the input phi
      call get_poroelastic_Biot_coeff(phi,kappa_s,kappa_f,kappa_fr,mu_fr,D_biot,H_biot,C_biot,M_biot)

      ! permeability xx
      perm_xx = permeability(1,kmato(ispec))

      ! computes velocities
      call get_poroelastic_velocities(cpIsquare,cpIIsquare,cssquare,H_biot,C_biot,M_biot,mu_fr,phi, &
             tort,rho_s,rho_f,eta_f,perm_xx,f0_source(1),freq0,Q0,w_c,ATTENUATION_PORO_FLUID_PART)

      cpIloc = sqrt(cpIsquare)
      cpIIloc = sqrt(cpIIsquare)
      csloc = sqrt(cssquare)
    else
      material = kmato(ispec)
      mu = poroelastcoef(2,1,material)
      lambdaplus2mu  = poroelastcoef(3,1,material)
      denst = density(1,material)

      cpIloc = sqrt(lambdaplus2mu/denst)
      cpIIloc = 0.d0
      csloc = sqrt(mu/denst)
    endif

    vpImax_local = -HUGEVAL
    vpImin_local = HUGEVAL
    vpIImax_local = -HUGEVAL
    vpIImin_local = HUGEVAL
    vsmin_local = HUGEVAL

    distance_min_local = HUGEVAL
    distance_max_local = -HUGEVAL

    do j = 1,NGLLZ
      do i = 1,NGLLX

        !--- if heterogeneous formulation with external velocity model
        if (assign_external_model) then
          cpIloc = vpext(i,j,ispec)
          csloc = vsext(i,j,ispec)
          denst = rhoext(i,j,ispec)
        endif

        !--- compute min and max of velocity and density models
        vpImin = min(vpImin,cpIloc)
        vpImax = max(vpImax,cpIloc)

        ! ignore acoustic and elastic regions with cpII = 0
        if (cpIIloc > 1.d-20) vpIImin = min(vpIImin,cpIIloc)
        vpIImax = max(vpIImax,cpIIloc)

        ! ignore fluid regions with Vs = 0
        if (csloc > 1.d-20) vsmin = min(vsmin,csloc)
        vsmax = max(vsmax,csloc)

        densmin = min(densmin,denst)
        densmax = max(densmax,denst)

        vpImax_local = max(vpImax_local,cpIloc)
        vpImin_local = min(vpImin_local,cpIloc)
        vpIImax_local = max(vpIImax_local,cpIIloc)
        vpIImin_local = min(vpIImin_local,cpIIloc)
        vsmin_local = min(vsmin_local,csloc)

      enddo
    enddo

    ! compute minimum and maximum size of edges of this grid cell
    distance_1 = sqrt((coord(1,ibool(1,1,ispec)) - coord(1,ibool(NGLLX,1,ispec)))**2 + &
               (coord(2,ibool(1,1,ispec)) - coord(2,ibool(NGLLX,1,ispec)))**2)

    distance_2 = sqrt((coord(1,ibool(NGLLX,1,ispec)) - coord(1,ibool(NGLLX,NGLLZ,ispec)))**2 + &
               (coord(2,ibool(NGLLX,1,ispec)) - coord(2,ibool(NGLLX,NGLLZ,ispec)))**2)

    distance_3 = sqrt((coord(1,ibool(NGLLX,NGLLZ,ispec)) - coord(1,ibool(1,NGLLZ,ispec)))**2 + &
               (coord(2,ibool(NGLLX,NGLLZ,ispec)) - coord(2,ibool(1,NGLLZ,ispec)))**2)

    distance_4 = sqrt((coord(1,ibool(1,NGLLZ,ispec)) - coord(1,ibool(1,1,ispec)))**2 + &
               (coord(2,ibool(1,NGLLZ,ispec)) - coord(2,ibool(1,1,ispec)))**2)

    distance_min_local = min(distance_1,distance_2,distance_3,distance_4)
    distance_max_local = max(distance_1,distance_2,distance_3,distance_4)

    distance_min = min(distance_min,distance_min_local)
    distance_max = max(distance_max,distance_max_local)

    courant_stability_number_max = max(courant_stability_number_max, &
                                       vpImax_local * deltat / (distance_min_local * percent_GLL(NGLLX)))
    
    ! Safeguard: check if CFL is > 1. If so, warn user, so that they can know where the issue is and do something about it.
    ! Do not stop program.
    if(vpImax_local*deltat/(distance_min_local*percent_GLL(NGLLX))>1.) then
        write(*,*) "********************************"
        write(*,*) "*           WARNING            *"
        write(*,*) "********************************"
        write(*,*) "* CFL > 1 somewhere in the     *"
        write(*,*) "* mesh:                        *"
        write(*,*) "* (x, z) = (", coord(1,ibool(3,3,ispec)), ", ", coord(2,ibool(3,3,ispec)), ")"
        write(*,*) "* vpImax_local = ", vpImax_local
        write(*,*) "* deltat = ", deltat
        write(*,*) "* distance_min_local = ", distance_min_local
        write(*,*) "* Reduce vpImax_local, reduce  *"
        write(*,*) "* deltat, or increase          *"
        write(*,*) "* distance_min_local.          *"
        write(*,*) "********************************"
    endif

    ! check if fluid region with Vs = 0
    if (vsmin_local > 1.d-20) then
      lambdaSmin = min(lambdaSmin,vsmin_local / (distance_max_local / (NGLLX - 1)))
      lambdaSmax = max(lambdaSmax,vsmin_local / (distance_max_local / (NGLLX - 1)))
    else
      any_fluid_histo = .true.
      lambdaPmin_in_fluid_histo = min(lambdaPmin_in_fluid_histo,vpImin_local / (distance_max_local / (NGLLX - 1)))
      lambdaPmax_in_fluid_histo = max(lambdaPmax_in_fluid_histo,vpImin_local / (distance_max_local / (NGLLX - 1)))
    endif

    lambdaPImin = min(lambdaPImin,vpImin_local / (distance_max_local / (NGLLX - 1)))
    lambdaPImax = max(lambdaPImax,vpImin_local / (distance_max_local / (NGLLX - 1)))

    if (cpIIloc > 1.d-20) then
      lambdaPIImin = min(lambdaPIImin,vpIImin_local / (distance_max_local / (NGLLX - 1)))
      lambdaPIImax = max(lambdaPIImax,vpIImin_local / (distance_max_local / (NGLLX - 1)))
    endif

  enddo ! ispec

  ! global statistics
  call min_all_all_dp(vpImin, vpImin_glob)
  call max_all_all_dp(vpImax, vpImax_glob)
  call min_all_all_dp(vpIImin, vpIImin_glob)
  call max_all_all_dp(vpIImax, vpIImax_glob)
  call min_all_all_dp(vsmin, vsmin_glob)
  call max_all_all_dp(vsmax, vsmax_glob)
  call min_all_all_dp(densmin, densmin_glob)
  call max_all_all_dp(densmax, densmax_glob)
  call min_all_all_dp(distance_min, distance_min_glob)
  call max_all_all_dp(distance_max, distance_max_glob)
  call min_all_all_dp(lambdaPImin, lambdaPImin_glob)
  call max_all_all_dp(lambdaPImax, lambdaPImax_glob)
  call min_all_all_dp(lambdaPIImin, lambdaPIImin_glob)
  call max_all_all_dp(lambdaPIImax, lambdaPIImax_glob)
  call min_all_all_dp(lambdaSmin, lambdaSmin_glob)
  call max_all_all_dp(lambdaSmax, lambdaSmax_glob)
  call min_all_all_dp(lambdaPmin_in_fluid_histo, lambdaPmin_in_fluid_histo_glob)
  call max_all_all_dp(lambdaPmax_in_fluid_histo, lambdaPmax_in_fluid_histo_glob)
  call max_all_all_dp(courant_stability_number_max, courant_stability_max_glob)

  vpImin = vpImin_glob
  vpImax = vpImax_glob
  vpIImin = vpIImin_glob
  vpIImax = vpIImax_glob
  vsmin = vsmin_glob
  vsmax = vsmax_glob
  densmin = densmin_glob
  densmax = densmax_glob
  distance_min = distance_min_glob
  distance_max = distance_max_glob
  lambdaPImin = lambdaPImin_glob
  lambdaPImax = lambdaPImax_glob
  lambdaPIImin = lambdaPIImin_glob
  lambdaPIImax = lambdaPIImax_glob
  lambdaSmin = lambdaSmin_glob
  lambdaSmax = lambdaSmax_glob
  lambdaPmin_in_fluid_histo = lambdaPmin_in_fluid_histo_glob
  lambdaPmax_in_fluid_histo = lambdaPmax_in_fluid_histo_glob
  courant_stability_number_max = courant_stability_max_glob


  ! sets if any slice has fluid histogram
  call any_all_l(any_fluid_histo,any_fluid_histo_glob)

  if (myrank == 0) then
    if (.not. all_anisotropic) then
      write(IMAIN,*)
      write(IMAIN,*) '********'
      write(IMAIN,*) 'Model: P (or PI) velocity min,max = ',vpImin,vpImax
      write(IMAIN,*) 'Model: PII velocity min,max = ',vpIImin,vpIImax
      write(IMAIN,*) 'Model: S velocity min,max = ',vsmin,vsmax
      write(IMAIN,*) 'Model: density min,max = ',densmin,densmax
      write(IMAIN,*) '********'
      write(IMAIN,*)

      write(IMAIN,*)
      write(IMAIN,*) '*********************************************'
      write(IMAIN,*) '*** Verification of simulation parameters ***'
      write(IMAIN,*) '*********************************************'
      write(IMAIN,*)
      write(IMAIN,*) '*** Max grid size = ',distance_max
      write(IMAIN,*) '*** Min grid size = ',distance_min
      write(IMAIN,*) '*** Max/min ratio = ',distance_max / distance_min
      write(IMAIN,*)
      write(IMAIN,*) '*** Max CFL stability condition of the time scheme &
                         &based on P wave velocity (must be below about 0.50 or so) = ',courant_stability_number_max
      if(courant_stability_number_max>0.5 .or. courant_stability_number_max<0.4) then
        ! If CFL isn't optimal, advise user a better choice.
        write(IMAIN,*) "*** (Advised DT to get to 0.49: ", (0.49*deltat/courant_stability_number_max), ". ",&
                       " If done, multiply NSTEP by ",courant_stability_number_max/0.49," to get same end time.)"
      endif
      write(IMAIN,*)
      call flush_IMAIN()
    endif

    create_wavelength_histogram = .false.

    ! only if time source is not a Dirac or Heaviside (otherwise maximum frequency of spectrum undefined)
    ! and if source is not an initial field, for the same reason
    if (.not. initialfield) then
      f0max = -HUGEVAL

      do i = 1,NSOURCES

        ! excludes Dirac and Heaviside sources
        if (time_function_type(i) /= 4 .and. time_function_type(i) /= 5) then

          ! sets min/max frequency
          if (f0_source(i) > f0max) f0max = f0_source(i)

          if (i == NSOURCES) then
            write(IMAIN,*) '----'
            write(IMAIN,*) ' Nb pts / lambdaPI_fmax min = ',lambdaPImin/(2.5d0*f0max)
            write(IMAIN,*) ' Nb pts / lambdaPI_fmax max = ',lambdaPImax/(2.5d0*f0max)
            write(IMAIN,*) '----'
            write(IMAIN,*) ' Nb pts / lambdaPII_fmax min = ',lambdaPIImin/(2.5d0*f0max)
            write(IMAIN,*) ' Nb pts / lambdaPII_fmax max = ',lambdaPIImax/(2.5d0*f0max)
            write(IMAIN,*) '----'
            write(IMAIN,*) ' Nb pts / lambdaS_fmax min = ',lambdaSmin/(2.5d0*f0max)
            write(IMAIN,*) ' Nb pts / lambdaS_fmax max = ',lambdaSmax/(2.5d0*f0max)
            write(IMAIN,*) '----'

            ! for histogram
            lambdaPmin_in_fluid_histo = lambdaPmin_in_fluid_histo/(2.5d0*f0max)
            lambdaPmax_in_fluid_histo = lambdaPmax_in_fluid_histo/(2.5d0*f0max)

            lambdaSmin_histo = lambdaSmin/(2.5d0*f0max)
            lambdaSmax_histo = lambdaSmax/(2.5d0*f0max)

            create_wavelength_histogram = .true.

          endif

        endif
      enddo
    endif
  endif

  ! master sends to all others
  ier = 0
#ifdef USE_MPI
  call MPI_BCAST(lambdaPmin_in_fluid_histo,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(lambdaPmax_in_fluid_histo,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(lambdaSmin_histo,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(lambdaSmax_histo,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(f0max,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(create_wavelength_histogram,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ier)
#endif

!!!!!!!!!!!!!!!! DK DK: added histogram of minimum number of points per wavelength

!! DK DK take into account the fact that there is no S velocity in the fluid
!! DK DK in this case, for the fluid, use the P wave data

  if (create_wavelength_histogram) then
    ! create statistics about mesh sampling (number of points per wavelength)
    call check_grid_create_histogram(any_fluid_histo_glob,lambdaPmin_in_fluid_histo,lambdaPmax_in_fluid_histo, &
                                         lambdaSmin_histo,lambdaSmax_histo,f0max)
  endif


  ! creates a PostScript file with stability condition
  if (output_postscript_snapshot) then
    call check_grid_create_postscript(courant_stability_number_max,lambdaPImin,lambdaPImax,lambdaSmin,lambdaSmax)
  endif

  end subroutine check_grid


!
!-------------------------------------------------------------------------------------------------
!


  subroutine check_grid_create_histogram(any_fluid_histo_glob,lambdaPmin_in_fluid_histo,lambdaPmax_in_fluid_histo, &
                                         lambdaSmin_histo,lambdaSmax_histo,f0max)

! create statistics about mesh sampling (number of points per wavelength)

#ifdef USE_MPI
  use mpi
#endif

  use specfem_par
  use specfem_par_movie

  implicit none

  logical,intent(in) :: any_fluid_histo_glob
  double precision,intent(in) :: lambdaPmin_in_fluid_histo,lambdaPmax_in_fluid_histo
  double precision,intent(in) :: lambdaSmin_histo,lambdaSmax_histo
  double precision,intent(in) :: f0max

! option to display only part of the mesh and not the whole mesh,
! for instance to analyze Cuthill-McKee mesh partitioning etc.

  ! local parameters
  double precision :: vpImin_local,vpIImin_local,vsmin_local

  double precision :: kappa_s,kappa_f,kappa_fr,mu_s,mu_fr,rho_s,rho_f,eta_f,phi,tort,rho_bar
  double precision :: D_biot,H_biot,C_biot,M_biot

  double precision :: cpIloc,cpIIloc,csloc
  double precision :: cpIsquare,cpIIsquare,cssquare
  double precision :: w_c,perm_xx
  double precision :: denst
  double precision :: lambdaplus2mu,mu
  double precision :: distance_min_local,distance_max_local
  double precision :: distance_1,distance_2,distance_3,distance_4

  integer :: ier
  integer :: i,j,ispec,material

! for histogram of number of points per wavelength
  double precision :: min_nb_of_points_per_wavelength,max_nb_of_points_per_wavelength,nb_of_points_per_wavelength, &
                   scaling_factor,scaling_factor_S,scaling_factor_P

  integer, parameter :: NCLASSES = 20
  integer, dimension(0:NCLASSES-1) :: classes_wavelength,classes_wavelength_all

  integer :: iclass,nspec_all,ipass,nspec_counted,nspec_counted_all,nspec_counted_all_solid,nspec_counted_all_fluid
  double precision :: current_percent,total_percent


  nspec_counted_all_solid = 0
  nspec_counted_all_fluid = 0

  ! first pass is for S wave sampling in solid, second pass is for P wave sampling in fluid
  do ipass = 1,2

    nspec_counted = 0

    if (ipass == 1) then
      min_nb_of_points_per_wavelength = lambdaSmin_histo
      max_nb_of_points_per_wavelength = lambdaSmax_histo
      ! do not create this histogram if the model is entirely fluid
      if (.not. ELASTIC_SIMULATION .and. .not. POROELASTIC_SIMULATION) cycle
    else
      ! do not create this histogram if the model is entirely solid
      if (.not. any_fluid_histo_glob) cycle
      min_nb_of_points_per_wavelength = lambdaPmin_in_fluid_histo
      max_nb_of_points_per_wavelength = lambdaPmax_in_fluid_histo
    endif

    ! when the grid is regular and the medium is homogeneous, the minimum and the maximum are equal
    ! and thus we cannot create an histogram; in such a case, let us artificially create a non-empty range
    if (abs(max_nb_of_points_per_wavelength - min_nb_of_points_per_wavelength) < 1.d-10) then
      min_nb_of_points_per_wavelength = min_nb_of_points_per_wavelength * 0.99d0
      max_nb_of_points_per_wavelength = max_nb_of_points_per_wavelength * 1.01d0
    endif

    ! erase histogram of wavelength
    classes_wavelength(:) = 0

    ! loop on all the elements
    do ispec = 1,nspec

      material = kmato(ispec)

      if (ispec_is_poroelastic(ispec)) then

        ! gets poroelastic material
        call get_poroelastic_material(ispec,phi,tort,mu_s,kappa_s,rho_s,kappa_f,rho_f,eta_f,mu_fr,kappa_fr,rho_bar)
        denst = rho_s

        ! Biot coefficients for the input phi
        call get_poroelastic_Biot_coeff(phi,kappa_s,kappa_f,kappa_fr,mu_fr,D_biot,H_biot,C_biot,M_biot)

        ! permeability xx
        perm_xx = permeability(1,kmato(ispec))

        ! computes velocities
        call get_poroelastic_velocities(cpIsquare,cpIIsquare,cssquare,H_biot,C_biot,M_biot,mu_fr,phi, &
               tort,rho_s,rho_f,eta_f,perm_xx,f0_source(1),freq0,Q0,w_c,ATTENUATION_PORO_FLUID_PART)

        cpIloc = sqrt(cpIsquare)
        cpIIloc = sqrt(cpIIsquare)
        csloc = sqrt(cssquare)
      else
        mu = poroelastcoef(2,1,material)
        lambdaplus2mu  = poroelastcoef(3,1,material)
        denst = density(1,material)

        cpIloc = sqrt(lambdaplus2mu/denst)
        cpIIloc = 0.d0
        csloc = sqrt(mu/denst)
      endif

      vpImin_local = HUGEVAL
      vpIImin_local = HUGEVAL
      vsmin_local = HUGEVAL

      do j = 1,NGLLZ
        do i = 1,NGLLX

          !--- if heterogeneous formulation with external velocity model
          if (assign_external_model) then
            cpIloc = vpext(i,j,ispec)
            csloc = vsext(i,j,ispec)
          endif

          vpImin_local = min(vpImin_local,cpIloc)
          vpIImin_local = min(vpIImin_local,cpIIloc)
          vsmin_local = min(vsmin_local,csloc)

        enddo
      enddo

      ! compute minimum and maximum size of edges of this grid cell
      distance_1 = sqrt((coord(1,ibool(1,1,ispec)) - coord(1,ibool(NGLLX,1,ispec)))**2 + &
                 (coord(2,ibool(1,1,ispec)) - coord(2,ibool(NGLLX,1,ispec)))**2)

      distance_2 = sqrt((coord(1,ibool(NGLLX,1,ispec)) - coord(1,ibool(NGLLX,NGLLZ,ispec)))**2 + &
                 (coord(2,ibool(NGLLX,1,ispec)) - coord(2,ibool(NGLLX,NGLLZ,ispec)))**2)

      distance_3 = sqrt((coord(1,ibool(NGLLX,NGLLZ,ispec)) - coord(1,ibool(1,NGLLZ,ispec)))**2 + &
                 (coord(2,ibool(NGLLX,NGLLZ,ispec)) - coord(2,ibool(1,NGLLZ,ispec)))**2)

      distance_4 = sqrt((coord(1,ibool(1,NGLLZ,ispec)) - coord(1,ibool(1,1,ispec)))**2 + &
                 (coord(2,ibool(1,NGLLZ,ispec)) - coord(2,ibool(1,1,ispec)))**2)

      distance_min_local = min(distance_1,distance_2,distance_3,distance_4)
      distance_max_local = max(distance_1,distance_2,distance_3,distance_4)

      if (ipass == 1) then
        ! in first pass, only solid regions, thus ignore fluid regions with Vs = 0
        if (vsmin_local > 1.d-20) then
          nb_of_points_per_wavelength = vsmin_local / (distance_max_local / (NGLLX - 1))

          nspec_counted = nspec_counted + 1

          nb_of_points_per_wavelength = nb_of_points_per_wavelength/(2.5d0*f0max)

          ! store number of points per wavelength in histogram
          iclass = int((nb_of_points_per_wavelength - min_nb_of_points_per_wavelength) / &
                       (max_nb_of_points_per_wavelength - min_nb_of_points_per_wavelength) * dble(NCLASSES))
          if (iclass < 0) iclass = 0
          if (iclass > NCLASSES-1) iclass = NCLASSES-1
          classes_wavelength(iclass) = classes_wavelength(iclass) + 1

        endif

      else
        ! in second pass, only fluid regions, thus ignore solid regions with Vs > 0
        if (abs(vsmin_local) < 1.d-20) then
          if (vpIImin_local <= ZERO) then
            nb_of_points_per_wavelength = vpImin_local / (distance_max_local / (NGLLX - 1))
          else
            nb_of_points_per_wavelength = min(vpImin_local,vpIImin_local) / (distance_max_local / (NGLLX - 1))
          endif

          nspec_counted = nspec_counted + 1

          nb_of_points_per_wavelength = nb_of_points_per_wavelength/(2.5d0*f0max)

          ! store number of points per wavelength in histogram
          iclass = int((nb_of_points_per_wavelength - min_nb_of_points_per_wavelength) / &
                       (max_nb_of_points_per_wavelength - min_nb_of_points_per_wavelength) * dble(NCLASSES))
          if (iclass < 0) iclass = 0
          if (iclass > NCLASSES-1) iclass = NCLASSES-1
          classes_wavelength(iclass) = classes_wavelength(iclass) + 1
        endif

      endif

    enddo

#ifdef USE_MPI
    call MPI_REDUCE(classes_wavelength, classes_wavelength_all, NCLASSES, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ier)
#else
    ier = 0
    classes_wavelength_all(:) = classes_wavelength(:)
#endif

    ! gets total for all slices
    call sum_all_i(nspec,nspec_all)
    call sum_all_i(nspec_counted,nspec_counted_all)

    if (ipass == 1) then
      nspec_counted_all_solid = nspec_counted_all
    else
      nspec_counted_all_fluid = nspec_counted_all
    endif

    ! create histogram of wavelength and save in Gnuplot file
    if (myrank == 0) then
      ! user output
      write(IMAIN,*)
      write(IMAIN,*) '-----------------------------------------'
      write(IMAIN,*)
      if (ipass == 1) then
        write(IMAIN,*) 'histogram of min number of points per S wavelength in solid regions:'
        write(IMAIN,*)
        write(IMAIN,*) 'there are ',nspec_counted_all,' elements out of ',nspec_all,' in solid regions'
      else
        write(IMAIN,*) 'histogram of min number of points per P wavelength in fluid regions:'
        write(IMAIN,*)
        write(IMAIN,*) 'there are ',nspec_counted_all,' elements out of ',nspec_all,' in fluid regions'
      endif
      write(IMAIN,*) '  (i.e., ',sngl(100.d0*nspec_counted_all/dble(nspec_all)),'% of the total)'
      write(IMAIN,*)
      write(IMAIN,*) '(too small = poor resolution of calculations -'
      write(IMAIN,*) ' too big = wasting memory and CPU time)'
      write(IMAIN,*) '(threshold value is around 4.5 points per S wavelength'
      write(IMAIN,*) ' in elastic regions and 5.5 per P wavelength in fluid regions):'
      write(IMAIN,*)

      total_percent = 0.
      scaling_factor = max_nb_of_points_per_wavelength - min_nb_of_points_per_wavelength

      if (ipass == 1) then
        open(unit=14,file='OUTPUT_FILES/points_per_wavelength_histogram_S_in_solid.txt',status='unknown')
        scaling_factor_S = scaling_factor
      else
        open(unit=14,file='OUTPUT_FILES/points_per_wavelength_histogram_P_in_fluid.txt',status='unknown')
        scaling_factor_P = scaling_factor
      endif
      do iclass = 0,NCLASSES-1
        current_percent = 100.*dble(classes_wavelength_all(iclass))/dble(nspec_counted_all)
        total_percent = total_percent + current_percent
        write(IMAIN,*) sngl(min_nb_of_points_per_wavelength + scaling_factor*iclass/dble(NCLASSES)),' - ', &
            sngl(min_nb_of_points_per_wavelength + scaling_factor*(iclass+1)/dble(NCLASSES)),classes_wavelength_all(iclass), &
            ' ',sngl(current_percent),' %'
        write(14,*) 0.5*(sngl(min_nb_of_points_per_wavelength + scaling_factor*iclass/dble(NCLASSES)) + &
            sngl(min_nb_of_points_per_wavelength + scaling_factor*(iclass+1)/dble(NCLASSES))),' ',sngl(current_percent)
      enddo
      close(14)

      if (total_percent < 99.9d0 .or. total_percent > 100.1d0) then
        write(IMAIN,*) 'total percentage = ',total_percent,' %'
        stop 'total percentage should be 100%'
      else
        write(IMAIN,*)
        write(IMAIN,*) 'total percentage = ',total_percent,' %'
      endif

    endif ! of if myrank == 0

  enddo ! end of the two passes on S wavelength data and P wavelength data

  ! create script for Gnuplot histogram file
  if (myrank == 0) then

    open(unit=14,file='OUTPUT_FILES/plot_points_per_wavelength_histogram.gnu',status='unknown')
    write(14,*) 'set term wxt'

    if (nspec_counted_all_solid > 0) then
      write(14,*) '#set term gif'
      write(14,*) '#set output "points_per_wavelength_histogram_S_in_solid.gif"'
      write(14,*)
      write(14,*) 'set boxwidth ',real(scaling_factor_S/NCLASSES)
      write(14,*) 'set xlabel "Range of min number of points per S wavelength in solid"'
      write(14,*) 'set ylabel "Percentage of elements (%)"'
      write(14,*) 'set loadpath "./OUTPUT_FILES"'
      write(14,*) 'plot "points_per_wavelength_histogram_S_in_solid.txt" with boxes'
      write(14,*) 'pause -1 "hit any key..."'
    endif

    if (nspec_counted_all_fluid > 0) then
      write(14,*) '#set term gif'
      write(14,*) '#set output "points_per_wavelength_histogram_P_in_fluid.gif"'
      write(14,*)
      write(14,*) 'set boxwidth ',real(scaling_factor_P/NCLASSES)
      write(14,*) 'set xlabel "Range of min number of points per P wavelength in fluid"'
      write(14,*) 'set ylabel "Percentage of elements (%)"'
      write(14,*) 'set loadpath "./OUTPUT_FILES"'
      write(14,*) 'plot "points_per_wavelength_histogram_P_in_fluid.txt" with boxes'
      write(14,*) 'pause -1 "hit any key..."'
    endif

    close(14)

    write(IMAIN,*)
    write(IMAIN,*)
    write(IMAIN,*) 'total number of elements in fluid and solid regions = ',nspec_all
    write(IMAIN,*)
    call flush_IMAIN()

  endif ! of if myrank == 0


  end subroutine check_grid_create_histogram


!
!-------------------------------------------------------------------------------------------------
!


  subroutine check_grid_setup_GLLper(percent_GLL,NGLLX_MAX_STABILITY)

  use constants,only: NGLLX

  implicit none

  integer :: NGLLX_MAX_STABILITY
  double precision :: percent_GLL(NGLLX_MAX_STABILITY)

  if (NGLLX_MAX_STABILITY /= 15 ) stop 'check NGLLX_MAX_STABILITY is equal to 15 in check_grid.f90'

! define percentage of smallest distance between GLL points for NGLLX points
! percentages were computed by calling the GLL points routine for each degree
  percent_GLL(:) = 100.d0

  percent_GLL(2) = 100.d0
  percent_GLL(3) = 50.d0
  percent_GLL(4) = 27.639320225002102d0
  percent_GLL(5) = 17.267316464601141d0
  percent_GLL(6) = 11.747233803526763d0
  percent_GLL(7) = 8.4888051860716516d0
  percent_GLL(8) = 6.4129925745196719d0
  percent_GLL(9) = 5.0121002294269914d0
  percent_GLL(10) = 4.0233045916770571d0
  percent_GLL(11) = 3.2999284795970416d0
  percent_GLL(12) = 2.7550363888558858d0
  percent_GLL(13) = 2.3345076678918053d0
  percent_GLL(14) = 2.0032477366369594d0
  percent_GLL(15) = 1.7377036748080721d0

! convert to real percentage
  percent_GLL(:) = percent_GLL(:) / 100.d0

  if (NGLLX > NGLLX_MAX_STABILITY) then
    stop 'cannot estimate the stability condition for degree NGLLX > NGLLX_MAX_STABILITY'
  endif

  end subroutine check_grid_setup_GLLper


!
!-------------------------------------------------------------------------------------------------
!


  subroutine check_grid_create_postscript(courant_stability_number_max,lambdaPImin,lambdaPImax,lambdaSmin,lambdaSmax)

#ifdef USE_MPI
  use mpi
#endif

  use specfem_par
  use specfem_par_movie

  implicit none

! option to display only part of the mesh and not the whole mesh,
! for instance to analyze Cuthill-McKee mesh partitioning etc.

  double precision,intent(in) :: courant_stability_number_max,lambdaPImin,lambdaPImax,lambdaSmin,lambdaSmax

  ! local parameters
  double precision :: vpImax_local,vpImin_local,vsmin_local
  double precision :: kappa_s,kappa_f,kappa_fr,mu_s,mu_fr,rho_s,rho_f,eta_f,phi,tort,rho_bar
  double precision :: D_biot,H_biot,C_biot,M_biot

  double precision :: cpIloc,csloc
  double precision :: cpIsquare,cpIIsquare,cssquare
  double precision :: w_c,perm_xx
  double precision :: denst
  double precision :: lambdaplus2mu,mu
  double precision :: lambdaS_local,lambdaPI_local

  double precision :: distance_min_local,distance_max_local
  double precision :: distance_1,distance_2,distance_3,distance_4

! for the stability condition
! maximum polynomial degree for which we can compute the stability condition
  integer, parameter :: NGLLX_MAX_STABILITY = 15
  double precision :: percent_GLL(NGLLX_MAX_STABILITY)

! color palette
  integer, parameter :: NUM_COLORS = 236
  double precision, dimension(NUM_COLORS) :: red,green,blue

  double precision :: xmax,zmax,height,usoffset,sizex,sizez,courant_stability_number
  double precision :: x1,z1,x2,z2,ratio_page,xmin,zmin

  double precision  :: xmin_glob, xmax_glob, zmin_glob, zmax_glob
#ifdef USE_MPI
  integer  :: icol
#endif

  double precision, dimension(2,nspec*5)  :: coorg_send
  double precision, dimension(:,:), allocatable  :: coorg_recv
  integer, dimension(nspec)  :: RGB_send
  integer, dimension(:), allocatable  :: RGB_recv
  real, dimension(nspec)  :: greyscale_send
  real, dimension(:), allocatable  :: greyscale_recv
  integer :: nspec_recv
  integer :: num_ispec
  integer :: iproc
  integer :: ier
  integer :: i,j,ispec,material
  integer :: is,ir,in,nnum
  integer :: UPPER_LIMIT_DISPLAY

  ! sets display limit
  if (DISPLAY_SUBSET_OPTION == 1) then
    UPPER_LIMIT_DISPLAY = nspec
  else if (DISPLAY_SUBSET_OPTION == 2) then
    UPPER_LIMIT_DISPLAY = nspec_inner
  else if (DISPLAY_SUBSET_OPTION == 3) then
    UPPER_LIMIT_DISPLAY = nspec_outer
  else if (DISPLAY_SUBSET_OPTION == 4) then
    UPPER_LIMIT_DISPLAY = NSPEC_DISPLAY_SUBSET
  else
    stop 'incorrect value of DISPLAY_SUBSET_OPTION'
  endif

  ! checks limit
  if (UPPER_LIMIT_DISPLAY > nspec) &
    call exit_MPI(myrank,'cannot have UPPER_LIMIT_DISPLAY > nspec in postscript creation of check_grid.F90')

  ! define color palette in random order
  call set_color_palette(red,green,blue,NUM_COLORS)

  ! define percentage of smallest distance between GLL points for NGLLX points
  ! percentages were computed by calling the GLL points routine for each degree
  call check_grid_setup_GLLper(percent_GLL,NGLLX_MAX_STABILITY)

! A4 or US letter paper
  if (US_LETTER) then
    usoffset = 1.75d0
    sizex = 27.94d0
    sizez = 21.59d0
  else
    usoffset = 0.d0
    sizex = 29.7d0
    sizez = 21.d0
  endif

! height of domain numbers in centimeters
  height = 0.25d0

! get minimum and maximum values of mesh coordinates
  xmin = minval(coord(1,:))
  zmin = minval(coord(2,:))
  xmax = maxval(coord(1,:))
  zmax = maxval(coord(2,:))

  call min_all_all_dp(xmin, xmin_glob)
  call max_all_all_dp(xmax, xmax_glob)
  call min_all_all_dp(zmin, zmin_glob)
  call max_all_all_dp(zmax, zmax_glob)
  xmin = xmin_glob
  xmax = xmax_glob
  zmin = zmin_glob
  zmax = zmax_glob

! ratio of physical page size/size of the domain meshed
  ratio_page = min(RPERCENTZ*sizez/(zmax-zmin),RPERCENTX*sizex/(xmax-xmin)) / 100.d0


  if (myrank == 0) then

    write(IMAIN,*)
    write(IMAIN,*) 'Creating PostScript file with stability condition'

!
!---- open PostScript file
!
    open(unit=24,file='OUTPUT_FILES/mesh_stability.ps',status='unknown')

!
!---- write PostScript header
!
    write(24,10) simulation_title
    write(24,*) '/CM {28.5 mul} def'
    write(24,*) '/LR {rlineto} def'
    write(24,*) '/LT {lineto} def'
    write(24,*) '/L {lineto} def'
    write(24,*) '/MR {rmoveto} def'
    write(24,*) '/MV {moveto} def'
    write(24,*) '/M {moveto} def'
    write(24,*) '/ST {stroke} def'
    write(24,*) '/CP {closepath} def'
    write(24,*) '/RG {setrgbcolor} def'
    write(24,*) '/GF {gsave fill grestore} def'
    write(24,*) '% different useful symbols'
    write(24,*) '/Point {2 0 360 arc CP 0 setgray fill} def'
    write(24,*) '/VDot {-0.75 -1.5 MR 1.5 0 LR 0 3. LR -1.5 0 LR'
    write(24,*) 'CP fill} def'
    write(24,*) '/HDot {-1.5 -0.75 MR 3. 0 LR 0 1.5 LR -3. 0 LR'
    write(24,*) 'CP fill} def'
    write(24,*) '/Cross {gsave 0.05 CM setlinewidth'
    write(24,*) 'gsave 3 3 MR -6. -6. LR ST grestore'
    write(24,*) 'gsave 3 -3 MR -6. 6. LR ST grestore'
    write(24,*) '0.01 CM setlinewidth} def'
    write(24,*) '/SmallLine {MV 0.07 CM 0 rlineto} def'
    write(24,*) '/Diamond {gsave 0.05 CM setlinewidth 0 4.2 MR'
    write(24,*) '-3 -4.2 LR 3 -4.2 LR 3 4.2 LR CP ST'
    write(24,*) 'grestore 0.01 CM setlinewidth} def'
    write(24,*) '%'
    write(24,*) '% macro to draw the contour of the elements'
    write(24,*) '/CO {M counttomark 2 idiv {L} repeat cleartomark CP} def'
    write(24,*) '%'
    write(24,*) '.01 CM setlinewidth'
    write(24,*) '/Times-Roman findfont'
    write(24,*) '.35 CM scalefont setfont'
    write(24,*) '%'
    write(24,*) '/vshift ',-height/2,' CM def'
    write(24,*) '/Rshow { currentpoint stroke MV'
    write(24,*) 'dup stringwidth pop neg vshift MR show } def'
    write(24,*) '/Cshow { currentpoint stroke MV'
    write(24,*) 'dup stringwidth pop -2 div vshift MR show } def'
    write(24,*) '/fN {/Helvetica-Bold findfont ',height,' CM scalefont setfont} def'
    write(24,*) '%'
    write(24,*) 'gsave newpath 90 rotate'
    write(24,*) '0 ',-sizez,' CM translate 1. 1. scale'
    write(24,*) '%'

    !
    !--- write captions of PostScript figure
    !
    write(24,*) '0 setgray'
    write(24,*) '/Times-Roman findfont'
    write(24,*) '.5 CM scalefont setfont'

    write(24,*) '%'
    write(24,*) '/Times-Roman findfont'
    write(24,*) '.6 CM scalefont setfont'
    write(24,*) '.4 .9 .9 setrgbcolor'
    write(24,*) '11 CM 1.1 CM MV'
    write(24,*) '(X axis) show'
    write(24,*) '%'
    write(24,*) '1.4 CM 9.5 CM MV'
    write(24,*) 'currentpoint gsave translate 90 rotate 0 0 moveto'
    write(24,*) '(Z axis) show'
    write(24,*) 'grestore'
    write(24,*) '%'
    write(24,*) '/Times-Roman findfont'
    write(24,*) '.7 CM scalefont setfont'
    write(24,*) '.8 0 .8 setrgbcolor'
    write(24,*) '24.35 CM 18.9 CM MV'
    write(24,*) usoffset,' CM 2 div neg 0 MR'
    write(24,*) 'currentpoint gsave translate -90 rotate 0 0 moveto'
    write(24,*) '(Mesh stability condition \(red = bad\)) show'
    write(24,*) 'grestore'
    write(24,*) '25.35 CM 18.9 CM MV'
    write(24,*) usoffset,' CM 2 div neg 0 MR'
    write(24,*) 'currentpoint gsave translate -90 rotate 0 0 moveto'
    write(24,*) '(',simulation_title,') show'
    write(24,*) 'grestore'
    write(24,*) '26.45 CM 18.9 CM MV'
    write(24,*) usoffset,' CM 2 div neg 0 MR'
    write(24,*) 'currentpoint gsave translate -90 rotate 0 0 moveto'
    write(24,*) '(2D Spectral Element Method) show'
    write(24,*) 'grestore'

    write(24,*) '%'
    write(24,*) '1 1 scale'
    write(24,*) '%'

    !
    !---- draw the spectral element mesh
    !
    write(24,*) '%'
    write(24,*) '% spectral element mesh'
    write(24,*) '%'
    write(24,*) '0 setgray'

    num_ispec = 0
  endif

  do ispec = 1, nspec
    if (myrank == 0) then
      num_ispec = num_ispec + 1
      write(24,*) '% elem ',num_ispec
    endif

    do i = 1,pointsdisp
      do j = 1,pointsdisp
        xinterp(i,j) = 0.d0
        zinterp(i,j) = 0.d0
        do in = 1,ngnod
          nnum = knods(in,ispec)
          xinterp(i,j) = xinterp(i,j) + shape2D_display(in,i,j)*coorg(1,nnum)
          zinterp(i,j) = zinterp(i,j) + shape2D_display(in,i,j)*coorg(2,nnum)
        enddo
      enddo
    enddo

    is = 1
    ir = 1
    x1 = (xinterp(ir,is)-xmin)*ratio_page + ORIG_X
    z1 = (zinterp(ir,is)-zmin)*ratio_page + ORIG_Z
    x1 = x1 * CENTIM
    z1 = z1 * CENTIM
    if (myrank == 0) then
      write(24,*) 'mark'
      write(24,681) x1,z1
    else
      coorg_send(1,(ispec-1)*5+1) = x1
      coorg_send(2,(ispec-1)*5+1) = z1
    endif

    ! draw straight lines if elements have 4 nodes

    ir=pointsdisp
    x2 = (xinterp(ir,is)-xmin)*ratio_page + ORIG_X
    z2 = (zinterp(ir,is)-zmin)*ratio_page + ORIG_Z
    x2 = x2 * CENTIM
    z2 = z2 * CENTIM
    if (myrank == 0) then
      write(24,681) x2,z2
    else
      coorg_send(1,(ispec-1)*5+2) = x2
      coorg_send(2,(ispec-1)*5+2) = z2
    endif

    ir=pointsdisp
    is=pointsdisp
    x2 = (xinterp(ir,is)-xmin)*ratio_page + ORIG_X
    z2 = (zinterp(ir,is)-zmin)*ratio_page + ORIG_Z
    x2 = x2 * CENTIM
    z2 = z2 * CENTIM
    if (myrank == 0) then
      write(24,681) x2,z2
    else
      coorg_send(1,(ispec-1)*5+3) = x2
      coorg_send(2,(ispec-1)*5+3) = z2
    endif

    is=pointsdisp
    ir=1
    x2 = (xinterp(ir,is)-xmin)*ratio_page + ORIG_X
    z2 = (zinterp(ir,is)-zmin)*ratio_page + ORIG_Z
    x2 = x2 * CENTIM
    z2 = z2 * CENTIM
    if (myrank == 0) then
      write(24,681) x2,z2
    else
      coorg_send(1,(ispec-1)*5+4) = x2
      coorg_send(2,(ispec-1)*5+4) = z2
    endif

    ir=1
    is=2
    x2 = (xinterp(ir,is)-xmin)*ratio_page + ORIG_X
    z2 = (zinterp(ir,is)-zmin)*ratio_page + ORIG_Z
    x2 = x2 * CENTIM
    z2 = z2 * CENTIM
    if (myrank == 0) then
      write(24,681) x2,z2
      write(24,*) 'CO'
    else
      coorg_send(1,(ispec-1)*5+5) = x2
      coorg_send(2,(ispec-1)*5+5) = z2
    endif


    if (ispec_is_poroelastic(ispec)) then
      ! poroelastic material
      call get_poroelastic_material(ispec,phi,tort,mu_s,kappa_s,rho_s,kappa_f,rho_f,eta_f,mu_fr,kappa_fr,rho_bar)
      denst = rho_s

      ! Biot coefficients for the input phi
      call get_poroelastic_Biot_coeff(phi,kappa_s,kappa_f,kappa_fr,mu_fr,D_biot,H_biot,C_biot,M_biot)

      ! permeability xx
      perm_xx = permeability(1,kmato(ispec))

      ! computes velocities
      call get_poroelastic_velocities(cpIsquare,cpIIsquare,cssquare,H_biot,C_biot,M_biot,mu_fr,phi, &
           tort,rho_s,rho_f,eta_f,perm_xx,f0_source(1),freq0,Q0,w_c,ATTENUATION_PORO_FLUID_PART)

      cpIloc = sqrt(cpIsquare)
    else
      material = kmato(ispec)

      lambdaplus2mu  = poroelastcoef(3,1,material)
      denst = density(1,material)

      cpIloc = sqrt(lambdaplus2mu/denst)
    endif

    vpImax_local = -HUGEVAL

    do j = 1,NGLLZ
      do i = 1,NGLLX

        !--- if heterogeneous formulation with external velocity model
        if (assign_external_model) then
          cpIloc = vpext(i,j,ispec)
          denst = rhoext(i,j,ispec)
        endif

        vpImax_local = max(vpImax_local,cpIloc)

      enddo
    enddo

! compute minimum and maximum size of edges of this grid cell
    distance_1 = sqrt((coord(1,ibool(1,1,ispec)) - coord(1,ibool(NGLLX,1,ispec)))**2 + &
             (coord(2,ibool(1,1,ispec)) - coord(2,ibool(NGLLX,1,ispec)))**2)

    distance_2 = sqrt((coord(1,ibool(NGLLX,1,ispec)) - coord(1,ibool(NGLLX,NGLLZ,ispec)))**2 + &
             (coord(2,ibool(NGLLX,1,ispec)) - coord(2,ibool(NGLLX,NGLLZ,ispec)))**2)

    distance_3 = sqrt((coord(1,ibool(NGLLX,NGLLZ,ispec)) - coord(1,ibool(1,NGLLZ,ispec)))**2 + &
             (coord(2,ibool(NGLLX,NGLLZ,ispec)) - coord(2,ibool(1,NGLLZ,ispec)))**2)

    distance_4 = sqrt((coord(1,ibool(1,NGLLZ,ispec)) - coord(1,ibool(1,1,ispec)))**2 + &
             (coord(2,ibool(1,NGLLZ,ispec)) - coord(2,ibool(1,1,ispec)))**2)

    distance_min_local = min(distance_1,distance_2,distance_3,distance_4)
    distance_max_local = max(distance_1,distance_2,distance_3,distance_4)

    courant_stability_number = vpImax_local * deltat / (distance_min_local * percent_GLL(NGLLX))

! display bad elements that are above the threshold
    if (courant_stability_number >= THRESHOLD_POSTSCRIPT * courant_stability_number_max) then
      if (myrank == 0) then
        write(24,*) '1 0 0 RG GF 0 setgray ST'
      else
        RGB_send(ispec) = 1
      endif
    else
! do not color the elements if below the threshold
      if (myrank == 0) then
        write(24,*) 'ST'
      else
        RGB_send(ispec) = 0
      endif
    endif

  enddo ! end of loop on all the spectral elements

#ifdef USE_MPI
  if (myrank == 0) then

    do iproc = 1, NPROC-1
      call MPI_RECV (nspec_recv, 1, MPI_INTEGER, &
              iproc, 42, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier)
      allocate(coorg_recv(2,nspec_recv*5))
      allocate(RGB_recv(nspec_recv))
      call MPI_RECV (coorg_recv(1,1), nspec_recv*5*2, MPI_DOUBLE_PRECISION, &
              iproc, 42, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier)
      call MPI_RECV (RGB_recv(1), nspec_recv, MPI_INTEGER, &
              iproc, 42, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier)

      do ispec = 1, nspec_recv
        num_ispec = num_ispec + 1
        write(24,*) '% elem ',num_ispec
        write(24,*) 'mark'
        write(24,681) coorg_recv(1,(ispec-1)*5+1), coorg_recv(2,(ispec-1)*5+1)
        write(24,681) coorg_recv(1,(ispec-1)*5+2), coorg_recv(2,(ispec-1)*5+2)
        write(24,681) coorg_recv(1,(ispec-1)*5+3), coorg_recv(2,(ispec-1)*5+3)
        write(24,681) coorg_recv(1,(ispec-1)*5+4), coorg_recv(2,(ispec-1)*5+4)
        write(24,681) coorg_recv(1,(ispec-1)*5+5), coorg_recv(2,(ispec-1)*5+5)
        write(24,*) 'CO'
        if (RGB_recv(ispec)  == 1) then
          write(24,*) '1 0 0 RG GF 0 setgray ST'
        else
          write(24,*) 'ST'
        endif
      enddo
      deallocate(coorg_recv)
      deallocate(RGB_recv)

    enddo

  else
    call MPI_SEND (nspec, 1, MPI_INTEGER, 0, 42, MPI_COMM_WORLD, ier)
    call MPI_SEND (coorg_send, nspec*5*2, MPI_DOUBLE_PRECISION, 0, 42, MPI_COMM_WORLD, ier)
    call MPI_SEND (RGB_send, nspec, MPI_INTEGER, 0, 42, MPI_COMM_WORLD, ier)
  endif
#else
  ! dummy statements to avoid compiler warnings
  allocate(coorg_recv(1,1))
  allocate(RGB_recv(1))
  nspec_recv = 0
  ier = 0
  iproc = NPROC
  deallocate(coorg_recv)
  deallocate(RGB_recv)
#endif

  if (myrank == 0) then
    write(24,*) '%'
    write(24,*) 'grestore'
    write(24,*) 'showpage'

    close(24)

    write(IMAIN,*) 'End of creation of PostScript file with stability condition'
  endif

!
!--------------------------------------------------------------------------------
!

  if (myrank == 0) then

    write(IMAIN,*)
    write(IMAIN,*) 'Creating PostScript file with mesh dispersion'

!
!---- open PostScript file
!
    if (ELASTIC_SIMULATION .or. POROELASTIC_SIMULATION) then
      open(unit=24,file='OUTPUT_FILES/mesh_S_wave_dispersion.ps',status='unknown')
    else
      open(unit=24,file='OUTPUT_FILES/mesh_P_wave_dispersion.ps',status='unknown')
    endif

!
!---- write PostScript header
!
    write(24,10) simulation_title
    write(24,*) '/CM {28.5 mul} def'
    write(24,*) '/LR {rlineto} def'
    write(24,*) '/LT {lineto} def'
    write(24,*) '/L {lineto} def'
    write(24,*) '/MR {rmoveto} def'
    write(24,*) '/MV {moveto} def'
    write(24,*) '/M {moveto} def'
    write(24,*) '/ST {stroke} def'
    write(24,*) '/CP {closepath} def'
    write(24,*) '/RG {setrgbcolor} def'
    write(24,*) '/GF {gsave fill grestore} def'
    write(24,*) '% different useful symbols'
    write(24,*) '/Point {2 0 360 arc CP 0 setgray fill} def'
    write(24,*) '/VDot {-0.75 -1.5 MR 1.5 0 LR 0 3. LR -1.5 0 LR'
    write(24,*) 'CP fill} def'
    write(24,*) '/HDot {-1.5 -0.75 MR 3. 0 LR 0 1.5 LR -3. 0 LR'
    write(24,*) 'CP fill} def'
    write(24,*) '/Cross {gsave 0.05 CM setlinewidth'
    write(24,*) 'gsave 3 3 MR -6. -6. LR ST grestore'
    write(24,*) 'gsave 3 -3 MR -6. 6. LR ST grestore'
    write(24,*) '0.01 CM setlinewidth} def'
    write(24,*) '/SmallLine {MV 0.07 CM 0 rlineto} def'
    write(24,*) '/Diamond {gsave 0.05 CM setlinewidth 0 4.2 MR'
    write(24,*) '-3 -4.2 LR 3 -4.2 LR 3 4.2 LR CP ST'
    write(24,*) 'grestore 0.01 CM setlinewidth} def'
    write(24,*) '%'
    write(24,*) '% macro to draw the contour of the elements'
    write(24,*) '/CO {M counttomark 2 idiv {L} repeat cleartomark CP} def'
    write(24,*) '%'
    write(24,*) '.01 CM setlinewidth'
    write(24,*) '/Times-Roman findfont'
    write(24,*) '.35 CM scalefont setfont'
    write(24,*) '%'
    write(24,*) '/vshift ',-height/2,' CM def'
    write(24,*) '/Rshow { currentpoint stroke MV'
    write(24,*) 'dup stringwidth pop neg vshift MR show } def'
    write(24,*) '/Cshow { currentpoint stroke MV'
    write(24,*) 'dup stringwidth pop -2 div vshift MR show } def'
    write(24,*) '/fN {/Helvetica-Bold findfont ',height,' CM scalefont setfont} def'
    write(24,*) '%'
    write(24,*) 'gsave newpath 90 rotate'
    write(24,*) '0 ',-sizez,' CM translate 1. 1. scale'
    write(24,*) '%'

!
!--- write captions of PostScript figure
!
    write(24,*) '0 setgray'
    write(24,*) '/Times-Roman findfont'
    write(24,*) '.5 CM scalefont setfont'

    write(24,*) '%'
    write(24,*) '/Times-Roman findfont'
    write(24,*) '.6 CM scalefont setfont'
    write(24,*) '.4 .9 .9 setrgbcolor'
    write(24,*) '11 CM 1.1 CM MV'
    write(24,*) '(X axis) show'
    write(24,*) '%'
    write(24,*) '1.4 CM 9.5 CM MV'
    write(24,*) 'currentpoint gsave translate 90 rotate 0 0 moveto'
    write(24,*) '(Z axis) show'
    write(24,*) 'grestore'
    write(24,*) '%'
    write(24,*) '/Times-Roman findfont'
    write(24,*) '.7 CM scalefont setfont'
    write(24,*) '.8 0 .8 setrgbcolor'
    write(24,*) '24.35 CM 18.9 CM MV'
    write(24,*) usoffset,' CM 2 div neg 0 MR'
    write(24,*) 'currentpoint gsave translate -90 rotate 0 0 moveto'
    if (ELASTIC_SIMULATION) then
      write(24,*) '(Mesh elastic S-wave dispersion \(red = good, blue = bad\)) show'
    else
      write(24,*) '(Mesh acoustic P-wave dispersion \(red = good, blue = bad\)) show'
    endif
    write(24,*) 'grestore'
    write(24,*) '25.35 CM 18.9 CM MV'
    write(24,*) usoffset,' CM 2 div neg 0 MR'
    write(24,*) 'currentpoint gsave translate -90 rotate 0 0 moveto'
    write(24,*) '(',simulation_title,') show'
    write(24,*) 'grestore'
    write(24,*) '26.45 CM 18.9 CM MV'
    write(24,*) usoffset,' CM 2 div neg 0 MR'
    write(24,*) 'currentpoint gsave translate -90 rotate 0 0 moveto'
    write(24,*) '(2D Spectral Element Method) show'
    write(24,*) 'grestore'

    write(24,*) '%'
    write(24,*) '1 1 scale'
    write(24,*) '%'

!
!---- draw the spectral element mesh
!
    write(24,*) '%'
    write(24,*) '% spectral element mesh'
    write(24,*) '%'
    write(24,*) '0 setgray'

    num_ispec = 0
  endif

  do ispec = 1, nspec
    if (myrank == 0) then
      num_ispec = num_ispec + 1
      write(24,*) '% elem ',num_ispec
    endif

    do i = 1,pointsdisp
      do j = 1,pointsdisp
        xinterp(i,j) = 0.d0
        zinterp(i,j) = 0.d0
        do in = 1,ngnod
          nnum = knods(in,ispec)
          xinterp(i,j) = xinterp(i,j) + shape2D_display(in,i,j)*coorg(1,nnum)
          zinterp(i,j) = zinterp(i,j) + shape2D_display(in,i,j)*coorg(2,nnum)
        enddo
      enddo
    enddo

    is = 1
    ir = 1
    x1 = (xinterp(ir,is)-xmin)*ratio_page + ORIG_X
    z1 = (zinterp(ir,is)-zmin)*ratio_page + ORIG_Z
    x1 = x1 * CENTIM
    z1 = z1 * CENTIM
    if (myrank == 0) then
      write(24,*) 'mark'
      write(24,681) x1,z1
    else
      coorg_send(1,(ispec-1)*5+1) = x1
      coorg_send(2,(ispec-1)*5+1) = z1
    endif

! draw straight lines if elements have 4 nodes

    ir=pointsdisp
    x2 = (xinterp(ir,is)-xmin)*ratio_page + ORIG_X
    z2 = (zinterp(ir,is)-zmin)*ratio_page + ORIG_Z
    x2 = x2 * CENTIM
    z2 = z2 * CENTIM
    if (myrank == 0) then
       write(24,681) x2,z2
    else
       coorg_send(1,(ispec-1)*5+2) = x2
       coorg_send(2,(ispec-1)*5+2) = z2
    endif

    ir=pointsdisp
    is=pointsdisp
    x2 = (xinterp(ir,is)-xmin)*ratio_page + ORIG_X
    z2 = (zinterp(ir,is)-zmin)*ratio_page + ORIG_Z
    x2 = x2 * CENTIM
    z2 = z2 * CENTIM
    if (myrank == 0) then
       write(24,681) x2,z2
    else
       coorg_send(1,(ispec-1)*5+3) = x2
       coorg_send(2,(ispec-1)*5+3) = z2
    endif

    is=pointsdisp
    ir=1
    x2 = (xinterp(ir,is)-xmin)*ratio_page + ORIG_X
    z2 = (zinterp(ir,is)-zmin)*ratio_page + ORIG_Z
    x2 = x2 * CENTIM
    z2 = z2 * CENTIM
    if (myrank == 0) then
       write(24,681) x2,z2
    else
       coorg_send(1,(ispec-1)*5+4) = x2
       coorg_send(2,(ispec-1)*5+4) = z2
    endif

    ir=1
    is=2
    x2 = (xinterp(ir,is)-xmin)*ratio_page + ORIG_X
    z2 = (zinterp(ir,is)-zmin)*ratio_page + ORIG_Z
    x2 = x2 * CENTIM
    z2 = z2 * CENTIM
    if (myrank == 0) then
       write(24,681) x2,z2
       write(24,*) 'CO'
    else
       coorg_send(1,(ispec-1)*5+5) = x2
       coorg_send(2,(ispec-1)*5+5) = z2
    endif

    if (ispec_is_poroelastic(ispec)) then
      ! gets poroelastic material
      call get_poroelastic_material(ispec,phi,tort,mu_s,kappa_s,rho_s,kappa_f,rho_f,eta_f,mu_fr,kappa_fr,rho_bar)
      denst = rho_s

      ! Biot coefficients for the input phi
      call get_poroelastic_Biot_coeff(phi,kappa_s,kappa_f,kappa_fr,mu_fr,D_biot,H_biot,C_biot,M_biot)

      ! permeability xx
      perm_xx = permeability(1,kmato(ispec))

      ! computes velocities
      call get_poroelastic_velocities(cpIsquare,cpIIsquare,cssquare,H_biot,C_biot,M_biot,mu_fr,phi, &
                                      tort,rho_s,rho_f,eta_f,perm_xx,f0_source(1),freq0,Q0,w_c,ATTENUATION_PORO_FLUID_PART)

      cpIloc = sqrt(cpIsquare)
      csloc = sqrt(cssquare)
    else
      material = kmato(ispec)
      mu = poroelastcoef(2,1,material)
      lambdaplus2mu  = poroelastcoef(3,1,material)
      denst = density(1,material)

      cpIloc = sqrt(lambdaplus2mu/denst)
      csloc = sqrt(mu/denst)
    endif

    vpImax_local = -HUGEVAL
    vpImin_local = HUGEVAL
    vsmin_local = HUGEVAL

    do j = 1,NGLLZ
      do i = 1,NGLLX
!--- if heterogeneous formulation with external velocity model
        if (assign_external_model) then
          cpIloc = vpext(i,j,ispec)
          csloc = vsext(i,j,ispec)
          denst = rhoext(i,j,ispec)
        endif

        vpImax_local = max(vpImax_local,cpIloc)
        vpImin_local = min(vpImin_local,cpIloc)
        vsmin_local = min(vsmin_local,csloc)
      enddo
    enddo

! compute minimum and maximum size of edges of this grid cell
    distance_1 = sqrt((coord(1,ibool(1,1,ispec)) - coord(1,ibool(NGLLX,1,ispec)))**2 + &
                 (coord(2,ibool(1,1,ispec)) - coord(2,ibool(NGLLX,1,ispec)))**2)

    distance_2 = sqrt((coord(1,ibool(NGLLX,1,ispec)) - coord(1,ibool(NGLLX,NGLLZ,ispec)))**2 + &
                 (coord(2,ibool(NGLLX,1,ispec)) - coord(2,ibool(NGLLX,NGLLZ,ispec)))**2)

    distance_3 = sqrt((coord(1,ibool(NGLLX,NGLLZ,ispec)) - coord(1,ibool(1,NGLLZ,ispec)))**2 + &
                 (coord(2,ibool(NGLLX,NGLLZ,ispec)) - coord(2,ibool(1,NGLLZ,ispec)))**2)

    distance_4 = sqrt((coord(1,ibool(1,NGLLZ,ispec)) - coord(1,ibool(1,1,ispec)))**2 + &
                 (coord(2,ibool(1,NGLLZ,ispec)) - coord(2,ibool(1,1,ispec)))**2)

    distance_min_local = min(distance_1,distance_2,distance_3,distance_4)
    distance_max_local = max(distance_1,distance_2,distance_3,distance_4)

! display mesh dispersion for S waves if there is at least one elastic element in the mesh
    if (ELASTIC_SIMULATION .or. POROELASTIC_SIMULATION) then

! ignore fluid regions with Vs = 0
      if (vsmin_local > 1.d-20) then

        lambdaS_local = vsmin_local / (distance_max_local / (NGLLX - 1))

! display very good elements that are above the threshold in red
        if (lambdaS_local >= THRESHOLD_POSTSCRIPT * lambdaSmax) then
          if (myrank == 0) then
            write(24,*) '1 0 0 RG GF 0 setgray ST'
          else
            RGB_send(ispec) = 1
          endif

! display bad elements that are below the threshold in blue
        else if (lambdaS_local <= (1. + (1. - THRESHOLD_POSTSCRIPT)) * lambdaSmin) then
          if (myrank == 0) then
            write(24,*) '0 0 1 RG GF 0 setgray ST'
          else
            RGB_send(ispec) = 3
          endif

        else
! do not color the elements if not close to the threshold
          if (myrank == 0) then
            write(24,*) 'ST'
          else
            RGB_send(ispec) = 0
          endif
        endif

      else
! do not color the elements if S-wave velocity undefined
        if (myrank == 0) then
          write(24,*) 'ST'
        else
          RGB_send(ispec) = 0
        endif
      endif

! display mesh dispersion for P waves if there is no elastic element in the mesh
    else

      lambdaPI_local = vpImin_local / (distance_max_local / (NGLLX - 1))

! display very good elements that are above the threshold in red
      if (lambdaPI_local >= THRESHOLD_POSTSCRIPT * lambdaPImax) then
        if (myrank == 0) then
          write(24,*) '1 0 0 RG GF 0 setgray ST'
        else
          RGB_send(ispec) = 1
        endif

! display bad elements that are below the threshold in blue
      else if (lambdaPI_local <= (1. + (1. - THRESHOLD_POSTSCRIPT)) * lambdaPImin) then
        if (myrank == 0) then
          write(24,*) '0 0 1 RG GF 0 setgray ST'
        else
          RGB_send(ispec) = 3
        endif

      else
! do not color the elements if not close to the threshold
        if (myrank == 0) then
          write(24,*) 'ST'
        else
          RGB_send(ispec) = 0
        endif
      endif

    endif

  enddo ! end of loop on all the spectral elements

#ifdef USE_MPI
  if (myrank == 0) then
    do iproc = 1, NPROC-1
      call MPI_RECV (nspec_recv, 1, MPI_INTEGER,iproc, 42, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier)
      allocate(coorg_recv(2,nspec_recv*5))
      allocate(RGB_recv(nspec_recv))
      call MPI_RECV (coorg_recv(1,1), nspec_recv*5*2, MPI_DOUBLE_PRECISION, &
            iproc, 42, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier)
      call MPI_RECV (RGB_recv(1), nspec_recv, MPI_INTEGER, &
            iproc, 42, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier)

      do ispec = 1, nspec_recv
        num_ispec = num_ispec + 1
        write(24,*) '% elem ',num_ispec
        write(24,*) 'mark'
        write(24,681) coorg_recv(1,(ispec-1)*5+1), coorg_recv(2,(ispec-1)*5+1)
        write(24,681) coorg_recv(1,(ispec-1)*5+2), coorg_recv(2,(ispec-1)*5+2)
        write(24,681) coorg_recv(1,(ispec-1)*5+3), coorg_recv(2,(ispec-1)*5+3)
        write(24,681) coorg_recv(1,(ispec-1)*5+4), coorg_recv(2,(ispec-1)*5+4)
        write(24,681) coorg_recv(1,(ispec-1)*5+5), coorg_recv(2,(ispec-1)*5+5)
        write(24,*) 'CO'
        if (RGB_recv(ispec)  == 1) then
          write(24,*) '1 0 0 RG GF 0 setgray ST'
        endif
        if (RGB_recv(ispec)  == 3) then
          write(24,*) '0 0 1 RG GF 0 setgray ST'
        endif
        if (RGB_recv(ispec)  == 0) then
          write(24,*) 'ST'
        endif

      enddo
      deallocate(coorg_recv)
      deallocate(RGB_recv)
    enddo
  else
    call MPI_SEND (nspec, 1, MPI_INTEGER, 0, 42, MPI_COMM_WORLD, ier)
    call MPI_SEND (coorg_send, nspec*5*2, MPI_DOUBLE_PRECISION, 0, 42, MPI_COMM_WORLD, ier)
    call MPI_SEND (RGB_send, nspec, MPI_INTEGER, 0, 42, MPI_COMM_WORLD, ier)
  endif
#endif

  if (myrank == 0) then
    write(24,*) '%'
    write(24,*) 'grestore'
    write(24,*) 'showpage'

    close(24)

    write(IMAIN,*) 'End of creation of PostScript file with mesh dispersion'

  endif

!
!--------------------------------------------------------------------------------
!

  if (myrank == 0) then

    write(IMAIN,*)
    write(IMAIN,*) 'Creating PostScript file with velocity model'

!
!---- open PostScript file
!
    open(unit=24,file='OUTPUT_FILES/P_velocity_model.ps',status='unknown')

!
!---- write PostScript header
!
    write(24,10) simulation_title
    write(24,*) '/CM {28.5 mul} def'
    write(24,*) '/LR {rlineto} def'
    write(24,*) '/LT {lineto} def'
    write(24,*) '/L {lineto} def'
    write(24,*) '/MR {rmoveto} def'
    write(24,*) '/MV {moveto} def'
    write(24,*) '/M {moveto} def'
    write(24,*) '/ST {stroke} def'
    write(24,*) '/CP {closepath} def'
    write(24,*) '/RG {setrgbcolor} def'
    write(24,*) '/GF {gsave fill grestore} def'
    write(24,*) '% different useful symbols'
    write(24,*) '/Point {2 0 360 arc CP 0 setgray fill} def'
    write(24,*) '/VDot {-0.75 -1.5 MR 1.5 0 LR 0 3. LR -1.5 0 LR'
    write(24,*) 'CP fill} def'
    write(24,*) '/HDot {-1.5 -0.75 MR 3. 0 LR 0 1.5 LR -3. 0 LR'
    write(24,*) 'CP fill} def'
    write(24,*) '/Cross {gsave 0.05 CM setlinewidth'
    write(24,*) 'gsave 3 3 MR -6. -6. LR ST grestore'
    write(24,*) 'gsave 3 -3 MR -6. 6. LR ST grestore'
    write(24,*) '0.01 CM setlinewidth} def'
    write(24,*) '/SmallLine {MV 0.07 CM 0 rlineto} def'
    write(24,*) '/Diamond {gsave 0.05 CM setlinewidth 0 4.2 MR'
    write(24,*) '-3 -4.2 LR 3 -4.2 LR 3 4.2 LR CP ST'
    write(24,*) 'grestore 0.01 CM setlinewidth} def'
    write(24,*) '%'
    write(24,*) '% macro to draw the contour of the elements'
    write(24,*) '/CO {M counttomark 2 idiv {L} repeat cleartomark CP} def'
    write(24,*) '%'
    write(24,*) '.01 CM setlinewidth'
    write(24,*) '/Times-Roman findfont'
    write(24,*) '.35 CM scalefont setfont'
    write(24,*) '%'
    write(24,*) '/vshift ',-height/2,' CM def'
    write(24,*) '/Rshow { currentpoint stroke MV'
    write(24,*) 'dup stringwidth pop neg vshift MR show } def'
    write(24,*) '/Cshow { currentpoint stroke MV'
    write(24,*) 'dup stringwidth pop -2 div vshift MR show } def'
    write(24,*) '/fN {/Helvetica-Bold findfont ',height,' CM scalefont setfont} def'
    write(24,*) '%'
    write(24,*) 'gsave newpath 90 rotate'
    write(24,*) '0 ',-sizez,' CM translate 1. 1. scale'
    write(24,*) '%'

!
!--- write captions of PostScript figure
!
    write(24,*) '0 setgray'
    write(24,*) '/Times-Roman findfont'
    write(24,*) '.5 CM scalefont setfont'

    write(24,*) '%'
    write(24,*) '/Times-Roman findfont'
    write(24,*) '.6 CM scalefont setfont'
    write(24,*) '.4 .9 .9 setrgbcolor'
    write(24,*) '11 CM 1.1 CM MV'
    write(24,*) '(X axis) show'
    write(24,*) '%'
    write(24,*) '1.4 CM 9.5 CM MV'
    write(24,*) 'currentpoint gsave translate 90 rotate 0 0 moveto'
    write(24,*) '(Z axis) show'
    write(24,*) 'grestore'
    write(24,*) '%'
    write(24,*) '/Times-Roman findfont'
    write(24,*) '.7 CM scalefont setfont'
    write(24,*) '.8 0 .8 setrgbcolor'
    write(24,*) '24.35 CM 18.9 CM MV'
    write(24,*) usoffset,' CM 2 div neg 0 MR'
    write(24,*) 'currentpoint gsave translate -90 rotate 0 0 moveto'
    write(24,*) '(P-velocity model \(dark = fast, light = slow\)) show'
    write(24,*) 'grestore'
    write(24,*) '25.35 CM 18.9 CM MV'
    write(24,*) usoffset,' CM 2 div neg 0 MR'
    write(24,*) 'currentpoint gsave translate -90 rotate 0 0 moveto'
    write(24,*) '(',simulation_title,') show'
    write(24,*) 'grestore'
    write(24,*) '26.45 CM 18.9 CM MV'
    write(24,*) usoffset,' CM 2 div neg 0 MR'
    write(24,*) 'currentpoint gsave translate -90 rotate 0 0 moveto'
    write(24,*) '(2D Spectral Element Method) show'
    write(24,*) 'grestore'

    write(24,*) '%'
    write(24,*) '1 1 scale'
    write(24,*) '%'

!
!---- draw the spectral element mesh
!
    write(24,*) '%'
    write(24,*) '% spectral element mesh'
    write(24,*) '%'
    write(24,*) '0 setgray'

    num_ispec = 0
  endif

  do ispec = 1, UPPER_LIMIT_DISPLAY
    if (myrank == 0) then
      num_ispec = num_ispec + 1
      write(24,*) '% elem ',num_ispec
    endif

    do i = 1,pointsdisp
      do j = 1,pointsdisp
        xinterp(i,j) = 0.d0
        zinterp(i,j) = 0.d0
        do in = 1,ngnod
          nnum = knods(in,ispec)
          xinterp(i,j) = xinterp(i,j) + shape2D_display(in,i,j)*coorg(1,nnum)
          zinterp(i,j) = zinterp(i,j) + shape2D_display(in,i,j)*coorg(2,nnum)
        enddo
      enddo
    enddo

    is = 1
    ir = 1
    x1 = (xinterp(ir,is)-xmin)*ratio_page + ORIG_X
    z1 = (zinterp(ir,is)-zmin)*ratio_page + ORIG_Z
    x1 = x1 * CENTIM
    z1 = z1 * CENTIM
    if (myrank == 0) then
       write(24,*) 'mark'
       write(24,681) x1,z1
    else
       coorg_send(1,(ispec-1)*5+1) = x1
       coorg_send(2,(ispec-1)*5+1) = z1
    endif

! draw straight lines if elements have 4 nodes

    ir=pointsdisp
    x2 = (xinterp(ir,is)-xmin)*ratio_page + ORIG_X
    z2 = (zinterp(ir,is)-zmin)*ratio_page + ORIG_Z
    x2 = x2 * CENTIM
    z2 = z2 * CENTIM
    if (myrank == 0) then
       write(24,681) x2,z2
    else
       coorg_send(1,(ispec-1)*5+2) = x2
       coorg_send(2,(ispec-1)*5+2) = z2
    endif

    ir=pointsdisp
    is=pointsdisp
    x2 = (xinterp(ir,is)-xmin)*ratio_page + ORIG_X
    z2 = (zinterp(ir,is)-zmin)*ratio_page + ORIG_Z
    x2 = x2 * CENTIM
    z2 = z2 * CENTIM
    if (myrank == 0) then
       write(24,681) x2,z2
    else
       coorg_send(1,(ispec-1)*5+3) = x2
       coorg_send(2,(ispec-1)*5+3) = z2
    endif

    is=pointsdisp
    ir=1
    x2 = (xinterp(ir,is)-xmin)*ratio_page + ORIG_X
    z2 = (zinterp(ir,is)-zmin)*ratio_page + ORIG_Z
    x2 = x2 * CENTIM
    z2 = z2 * CENTIM
    if (myrank == 0) then
       write(24,681) x2,z2
    else
       coorg_send(1,(ispec-1)*5+4) = x2
       coorg_send(2,(ispec-1)*5+4) = z2
    endif

    ir=1
    is=2
    x2 = (xinterp(ir,is)-xmin)*ratio_page + ORIG_X
    z2 = (zinterp(ir,is)-zmin)*ratio_page + ORIG_Z
    x2 = x2 * CENTIM
    z2 = z2 * CENTIM
    if (myrank == 0) then
       write(24,681) x2,z2
       write(24,*) 'CO'
    else
       coorg_send(1,(ispec-1)*5+5) = x2
       coorg_send(2,(ispec-1)*5+5) = z2
    endif

    if ((vpImax-vpImin)/vpImin > 0.02d0) then
      if (assign_external_model) then
! use lower-left corner
        x1 = (vpext(1,1,ispec)-vpImin) / (vpImax-vpImin)
      else
        if (ispec_is_poroelastic(ispec)) then
          ! gets poroelastic material
          call get_poroelastic_material(ispec,phi,tort,mu_s,kappa_s,rho_s,kappa_f,rho_f,eta_f,mu_fr,kappa_fr,rho_bar)
          denst = rho_s

          ! Biot coefficients for the input phi
          call get_poroelastic_Biot_coeff(phi,kappa_s,kappa_f,kappa_fr,mu_fr,D_biot,H_biot,C_biot,M_biot)

          ! permeability xx
          perm_xx = permeability(1,kmato(ispec))

          ! computes velocities
          call get_poroelastic_velocities(cpIsquare,cpIIsquare,cssquare,H_biot,C_biot,M_biot,mu_fr,phi, &
                                          tort,rho_s,rho_f,eta_f,perm_xx,f0_source(1),freq0,Q0,w_c,ATTENUATION_PORO_FLUID_PART)

          cpIloc = sqrt(cpIsquare)
        else
          material = kmato(ispec)

          lambdaplus2mu  = poroelastcoef(3,1,material)
          denst = density(1,material)
          cpIloc = sqrt(lambdaplus2mu/denst)
        endif
        x1 = (cpIloc-vpImin)/(vpImax-vpImin)
      endif
    else
      x1 = 0.5d0
    endif

! rescale to avoid very dark gray levels
    x1 = x1*0.7 + 0.2
    if (x1 > 1.d0) x1=1.d0

! invert scale: white = vpmin, dark gray = vpmax
    x1 = 1.d0 - x1

! display P-velocity model using gray levels
    if (myrank == 0) then
      write(24,*) sngl(x1),' setgray GF 0 setgray ST'
    else
      greyscale_send(ispec) = sngl(x1)
    endif

  enddo ! end of loop on all the spectral elements

#ifdef USE_MPI
  if (myrank == 0) then
    do iproc = 1, NPROC-1
      call MPI_RECV (nspec_recv, 1, MPI_INTEGER,iproc, 42, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier)
      allocate(coorg_recv(2,nspec_recv*5))
      allocate(greyscale_recv(nspec_recv))
      call MPI_RECV (coorg_recv(1,1), nspec_recv*5*2, MPI_DOUBLE_PRECISION, &
            iproc, 42, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier)
      call MPI_RECV (greyscale_recv(1), nspec_recv, MPI_REAL, &
            iproc, 42, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier)

      do ispec = 1, nspec_recv
        num_ispec = num_ispec + 1
        write(24,*) '% elem ',num_ispec
        write(24,*) 'mark'
        write(24,681) coorg_recv(1,(ispec-1)*5+1), coorg_recv(2,(ispec-1)*5+1)
        write(24,681) coorg_recv(1,(ispec-1)*5+2), coorg_recv(2,(ispec-1)*5+2)
        write(24,681) coorg_recv(1,(ispec-1)*5+3), coorg_recv(2,(ispec-1)*5+3)
        write(24,681) coorg_recv(1,(ispec-1)*5+4), coorg_recv(2,(ispec-1)*5+4)
        write(24,681) coorg_recv(1,(ispec-1)*5+5), coorg_recv(2,(ispec-1)*5+5)
        write(24,*) 'CO'
        write(24,*) greyscale_recv(ispec), ' setgray GF 0 setgray ST'
      enddo
      deallocate(coorg_recv)
      deallocate(greyscale_recv)
    enddo
  else
    call MPI_SEND (UPPER_LIMIT_DISPLAY, 1, MPI_INTEGER, 0, 42, MPI_COMM_WORLD, ier)
    call MPI_SEND (coorg_send, UPPER_LIMIT_DISPLAY*5*2, MPI_DOUBLE_PRECISION, 0, 42, MPI_COMM_WORLD, ier)
    call MPI_SEND (greyscale_send, UPPER_LIMIT_DISPLAY, MPI_INTEGER, 0, 42, MPI_COMM_WORLD, ier)
  endif
#else
  ! dummy statements to avoid compiler warnings
  allocate(greyscale_recv(1))
  ier = 0
  deallocate(greyscale_recv)
#endif

  if (myrank == 0) then

    write(24,*) '%'
    write(24,*) 'grestore'
    write(24,*) 'showpage'

    close(24)

    write(IMAIN,*) 'End of creation of PostScript file with velocity model'

  endif

  if (myrank == 0) then

    write(IMAIN,*)
    write(IMAIN,*) 'Creating PostScript file with mesh partitioning'

!
!---- open PostScript file
!
    open(unit=24,file='OUTPUT_FILES/mesh_partitioning.ps',status='unknown')

!
!---- write PostScript header
!
    write(24,10) simulation_title
    write(24,*) '/CM {28.5 mul} def'
    write(24,*) '/LR {rlineto} def'
    write(24,*) '/LT {lineto} def'
    write(24,*) '/L {lineto} def'
    write(24,*) '/MR {rmoveto} def'
    write(24,*) '/MV {moveto} def'
    write(24,*) '/M {moveto} def'
    write(24,*) '/ST {stroke} def'
    write(24,*) '/CP {closepath} def'
    write(24,*) '/RG {setrgbcolor} def'
    write(24,*) '/GF {gsave fill grestore} def'
    write(24,*) '% different useful symbols'
    write(24,*) '/Point {2 0 360 arc CP 0 setgray fill} def'
    write(24,*) '/VDot {-0.75 -1.5 MR 1.5 0 LR 0 3. LR -1.5 0 LR'
    write(24,*) 'CP fill} def'
    write(24,*) '/HDot {-1.5 -0.75 MR 3. 0 LR 0 1.5 LR -3. 0 LR'
    write(24,*) 'CP fill} def'
    write(24,*) '/Cross {gsave 0.05 CM setlinewidth'
    write(24,*) 'gsave 3 3 MR -6. -6. LR ST grestore'
    write(24,*) 'gsave 3 -3 MR -6. 6. LR ST grestore'
    write(24,*) '0.01 CM setlinewidth} def'
    write(24,*) '/SmallLine {MV 0.07 CM 0 rlineto} def'
    write(24,*) '/Diamond {gsave 0.05 CM setlinewidth 0 4.2 MR'
    write(24,*) '-3 -4.2 LR 3 -4.2 LR 3 4.2 LR CP ST'
    write(24,*) 'grestore 0.01 CM setlinewidth} def'
    write(24,*) '%'
    write(24,*) '% macro to draw the contour of the elements'
    write(24,*) '/CO {M counttomark 2 idiv {L} repeat cleartomark CP} def'
    write(24,*) '%'
    write(24,*) '.01 CM setlinewidth'
    write(24,*) '/Times-Roman findfont'
    write(24,*) '.35 CM scalefont setfont'
    write(24,*) '%'
    write(24,*) '/vshift ',-height/2,' CM def'
    write(24,*) '/Rshow { currentpoint stroke MV'
    write(24,*) 'dup stringwidth pop neg vshift MR show } def'
    write(24,*) '/Cshow { currentpoint stroke MV'
    write(24,*) 'dup stringwidth pop -2 div vshift MR show } def'
    write(24,*) '/fN {/Helvetica-Bold findfont ',height,' CM scalefont setfont} def'
    write(24,*) '%'
    write(24,*) 'gsave newpath 90 rotate'
    write(24,*) '0 ',-sizez,' CM translate 1. 1. scale'
    write(24,*) '%'

!
!--- write captions of PostScript figure
!
    write(24,*) '0 setgray'
    write(24,*) '/Times-Roman findfont'
    write(24,*) '.5 CM scalefont setfont'

    write(24,*) '%'
    write(24,*) '/Times-Roman findfont'
    write(24,*) '.6 CM scalefont setfont'
    write(24,*) '.4 .9 .9 setrgbcolor'
    write(24,*) '11 CM 1.1 CM MV'
    write(24,*) '(X axis) show'
    write(24,*) '%'
    write(24,*) '1.4 CM 9.5 CM MV'
    write(24,*) 'currentpoint gsave translate 90 rotate 0 0 moveto'
    write(24,*) '(Z axis) show'
    write(24,*) 'grestore'
    write(24,*) '%'
    write(24,*) '/Times-Roman findfont'
    write(24,*) '.7 CM scalefont setfont'
    write(24,*) '.8 0 .8 setrgbcolor'
    write(24,*) '24.35 CM 18.9 CM MV'
    write(24,*) usoffset,' CM 2 div neg 0 MR'
    write(24,*) 'currentpoint gsave translate -90 rotate 0 0 moveto'
    write(24,*) '(Mesh partitioning) show'
    write(24,*) 'grestore'
    write(24,*) '25.35 CM 18.9 CM MV'
    write(24,*) usoffset,' CM 2 div neg 0 MR'
    write(24,*) 'currentpoint gsave translate -90 rotate 0 0 moveto'
    write(24,*) '(',simulation_title,') show'
    write(24,*) 'grestore'
    write(24,*) '26.45 CM 18.9 CM MV'
    write(24,*) usoffset,' CM 2 div neg 0 MR'
    write(24,*) 'currentpoint gsave translate -90 rotate 0 0 moveto'
    write(24,*) '(2D Spectral Element Method) show'
    write(24,*) 'grestore'

    write(24,*) '%'
    write(24,*) '1 1 scale'
    write(24,*) '%'

!
!---- draw the spectral element mesh
!
    write(24,*) '%'
    write(24,*) '% spectral element mesh'
    write(24,*) '%'
    write(24,*) '0 setgray'

    num_ispec = 0
  endif

  do ispec = 1, UPPER_LIMIT_DISPLAY

    if (myrank == 0) then
      num_ispec = num_ispec + 1
      write(24,*) '% elem ',num_ispec
    endif

    do i = 1,pointsdisp
      do j = 1,pointsdisp
        xinterp(i,j) = 0.d0
        zinterp(i,j) = 0.d0
        do in = 1,ngnod
          nnum = knods(in,ispec)
          xinterp(i,j) = xinterp(i,j) + shape2D_display(in,i,j)*coorg(1,nnum)
          zinterp(i,j) = zinterp(i,j) + shape2D_display(in,i,j)*coorg(2,nnum)
        enddo
      enddo
    enddo

    is = 1
    ir = 1
    x1 = (xinterp(ir,is)-xmin)*ratio_page + ORIG_X
    z1 = (zinterp(ir,is)-zmin)*ratio_page + ORIG_Z
    x1 = x1 * CENTIM
    z1 = z1 * CENTIM
    if (myrank == 0) then
       write(24,*) 'mark'
       write(24,681) x1,z1
    else
       coorg_send(1,(ispec-1)*5+1) = x1
       coorg_send(2,(ispec-1)*5+1) = z1
    endif

! draw straight lines if elements have 4 nodes

    ir=pointsdisp
    x2 = (xinterp(ir,is)-xmin)*ratio_page + ORIG_X
    z2 = (zinterp(ir,is)-zmin)*ratio_page + ORIG_Z
    x2 = x2 * CENTIM
    z2 = z2 * CENTIM
    if (myrank == 0) then
       write(24,681) x2,z2
    else
       coorg_send(1,(ispec-1)*5+2) = x2
       coorg_send(2,(ispec-1)*5+2) = z2
    endif

    ir=pointsdisp
    is=pointsdisp
    x2 = (xinterp(ir,is)-xmin)*ratio_page + ORIG_X
    z2 = (zinterp(ir,is)-zmin)*ratio_page + ORIG_Z
    x2 = x2 * CENTIM
    z2 = z2 * CENTIM
    if (myrank == 0) then
       write(24,681) x2,z2
    else
       coorg_send(1,(ispec-1)*5+3) = x2
       coorg_send(2,(ispec-1)*5+3) = z2
    endif

    is=pointsdisp
    ir=1
    x2 = (xinterp(ir,is)-xmin)*ratio_page + ORIG_X
    z2 = (zinterp(ir,is)-zmin)*ratio_page + ORIG_Z
    x2 = x2 * CENTIM
    z2 = z2 * CENTIM
    if (myrank == 0) then
       write(24,681) x2,z2
    else
       coorg_send(1,(ispec-1)*5+4) = x2
       coorg_send(2,(ispec-1)*5+4) = z2
    endif

    ir=1
    is=2
    x2 = (xinterp(ir,is)-xmin)*ratio_page + ORIG_X
    z2 = (zinterp(ir,is)-zmin)*ratio_page + ORIG_Z
    x2 = x2 * CENTIM
    z2 = z2 * CENTIM
    if (myrank == 0) then
       write(24,681) x2,z2
       write(24,*) 'CO'
    else
       coorg_send(1,(ispec-1)*5+5) = x2
       coorg_send(2,(ispec-1)*5+5) = z2
    endif

    if (myrank == 0) then
      write(24,*) red(1), green(1), blue(1), 'RG GF 0 setgray ST'
    endif

  enddo ! end of loop on all the spectral elements

#ifdef USE_MPI
  if (myrank == 0) then
    do iproc = 1, NPROC-1
      ! use a different color for each material set
      icol = mod(iproc, NUM_COLORS) + 1

      call MPI_RECV (nspec_recv, 1, MPI_INTEGER, &
          iproc, 42, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier)
      allocate(coorg_recv(2,nspec_recv*5))
      call MPI_RECV (coorg_recv(1,1), nspec_recv*5*2, MPI_DOUBLE_PRECISION, &
          iproc, 42, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier)

      do ispec = 1, nspec_recv
         num_ispec = num_ispec + 1
         write(24,*) '% elem ',num_ispec
         write(24,*) 'mark'
         write(24,681) coorg_recv(1,(ispec-1)*5+1), coorg_recv(2,(ispec-1)*5+1)
         write(24,681) coorg_recv(1,(ispec-1)*5+2), coorg_recv(2,(ispec-1)*5+2)
         write(24,681) coorg_recv(1,(ispec-1)*5+3), coorg_recv(2,(ispec-1)*5+3)
         write(24,681) coorg_recv(1,(ispec-1)*5+4), coorg_recv(2,(ispec-1)*5+4)
         write(24,681) coorg_recv(1,(ispec-1)*5+5), coorg_recv(2,(ispec-1)*5+5)
         write(24,*) 'CO'

         write(24,*) red(icol), green(icol), blue(icol), ' RG GF 0 setgray ST'
      enddo
      deallocate(coorg_recv)
    enddo

  else
    call MPI_SEND (UPPER_LIMIT_DISPLAY, 1, MPI_INTEGER,0, 42, MPI_COMM_WORLD, ier)
    call MPI_SEND (coorg_send, UPPER_LIMIT_DISPLAY*5*2, MPI_DOUBLE_PRECISION,0, 42, MPI_COMM_WORLD, ier)
  endif
#endif

  if (myrank == 0) then
    write(24,*) '%'
    write(24,*) 'grestore'
    write(24,*) 'showpage'

    close(24)

    write(IMAIN,*) 'End of creation of PostScript file with partitioning'
    write(IMAIN,*)
  endif

 10  format('%!PS-Adobe-2.0',/,'%%',/,'%% Title: ',a100,/,'%% Created by: Specfem2D',/,'%% Author: Dimitri Komatitsch',/,'%%')

 681 format(f6.2,1x,f6.2)

  end subroutine check_grid_create_postscript

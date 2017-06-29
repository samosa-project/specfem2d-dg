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

  subroutine compute_stacey_elastic(accel_elastic,veloc_elastic)

! absorbing boundaries
!
! Clayton-Engquist condition if elastic

  use constants,only: CUSTOM_REAL,NGLLX,NGLLZ,NDIM, &
    ZERO,ONE,TWO,FOUR_THIRDS,IEDGE1,IEDGE2,IEDGE3,IEDGE4

  use specfem_par, only: AXISYM,nglob,nelemabs,it,any_elastic, &
                         assign_external_model,ibool,kmato,numabs,ispec_is_elastic, &
                         codeabs,codeabs_corner,density,poroelastcoef,xix,xiz,gammax,gammaz,jacobian, &
                         vpext,vsext,rhoext,wxgll,wzgll,P_SV, &
                         SIMULATION_TYPE,SAVE_FORWARD, &
                         b_absorb_elastic_left,b_absorb_elastic_right, &
                         b_absorb_elastic_bottom,b_absorb_elastic_top,&
                         ib_left,ib_right,ib_bottom,ib_top, &
                         STACEY_BOUNDARY_CONDITIONS,deltat

  ! initialfield
  use specfem_par,only: v0x_left,v0z_left,v0x_right,v0z_right,v0x_bot,v0z_bot, &
                        t0x_left,t0z_left,t0x_right,t0z_right,t0x_bot,t0z_bot, &
                        add_Bielak_conditions_bottom,add_Bielak_conditions_right, &
                        add_Bielak_conditions_top,add_Bielak_conditions_left, &
                        initialfield,over_critical_angle, &
                        anglesource,anglesource_refl,A_plane,B_plane,C_plane,c_inc,c_refl,time_offset

  ! for Bielak
  use specfem_par, only: x_source,z_source,f0_source,coord

  implicit none

  real(kind=CUSTOM_REAL), dimension(NDIM,nglob),intent(inout) :: accel_elastic
  real(kind=CUSTOM_REAL), dimension(NDIM,nglob),intent(in) :: veloc_elastic

  ! local parameters
  integer :: ispecabs,ispec,i,j,iglob
  integer :: ibegin,iend
  real(kind=CUSTOM_REAL) :: weight,xxi,zxi,xgamma,zgamma,jacobian1D
  real(kind=CUSTOM_REAL) :: nx,nz,vx,vy,vz,vn,rho_vp,rho_vs,tx,ty,tz

  ! material properties of the elastic medium
  real(kind=CUSTOM_REAL) :: mul_unrelaxed_elastic,lambdal_unrelaxed_elastic, &
    lambdaplus2mu_unrelaxed_elastic,kappal,cpl,csl,rhol

  ! for analytical initial plane wave for Bielak's conditions
  double precision :: veloc_horiz,veloc_vert,dxUx,dzUx,dxUz,dzUz,traction_x_t0,traction_z_t0
  integer :: count_left,count_right,count_bottom

  ! checks if anything to do
  if (.not. STACEY_BOUNDARY_CONDITIONS) return
  if (.not. any_elastic) return

  ! Clayton-Engquist condition if elastic
  count_left=1
  count_right=1
  count_bottom=1
  do ispecabs = 1,nelemabs

    ispec = numabs(ispecabs)

    ! only for elastic elements, skip others
    if (.not. ispec_is_elastic(ispec) ) cycle

    ! get elastic parameters of current spectral element
    lambdal_unrelaxed_elastic = poroelastcoef(1,1,kmato(ispec))
    mul_unrelaxed_elastic = poroelastcoef(2,1,kmato(ispec))
    lambdaplus2mu_unrelaxed_elastic = lambdal_unrelaxed_elastic + 2._CUSTOM_REAL * mul_unrelaxed_elastic
    rhol  = density(1,kmato(ispec))
    if (AXISYM) then ! CHECK kappa
      kappal  = lambdal_unrelaxed_elastic + TWO*mul_unrelaxed_elastic/3._CUSTOM_REAL
    else
      kappal  = lambdal_unrelaxed_elastic + mul_unrelaxed_elastic
    endif
    cpl = sqrt((kappal + FOUR_THIRDS * mul_unrelaxed_elastic)/rhol) ! CHECK kappa
    csl = sqrt(mul_unrelaxed_elastic/rhol)

    !--- left absorbing boundary
    if (codeabs(IEDGE4,ispecabs)) then
      i = 1
      do j = 1,NGLLZ
        ! Clayton-Engquist condition if elastic
        iglob = ibool(i,j,ispec)

        veloc_horiz = 0._CUSTOM_REAL
        veloc_vert = 0._CUSTOM_REAL
        traction_x_t0 = 0._CUSTOM_REAL
        traction_z_t0 = 0._CUSTOM_REAL

        ! for analytical initial plane wave for Bielak's conditions
        ! left or right edge, horizontal normal vector
        if (add_Bielak_conditions_left .and. initialfield) then
          if (.not. over_critical_angle) then
            call compute_Bielak_conditions(coord,iglob,nglob,it,deltat,dxUx,dxUz,dzUx,dzUz,veloc_horiz,veloc_vert, &
                        x_source(1), z_source(1), A_plane, B_plane, C_plane, anglesource(1), anglesource_refl, &
                        c_inc, c_refl, time_offset,f0_source(1))

            traction_x_t0 = lambdaplus2mu_unrelaxed_elastic * dxUx + lambdal_unrelaxed_elastic * dzUz
            traction_z_t0 = mul_unrelaxed_elastic * (dxUz + dzUx)
          else
            veloc_horiz = v0x_left(count_left,it)
            veloc_vert = v0z_left(count_left,it)
            traction_x_t0 = t0x_left(count_left,it)
            traction_z_t0 = t0z_left(count_left,it)
            count_left = count_left+1
          endif
        endif

        ! external velocity model
        if (assign_external_model) then
          cpl = vpext(i,j,ispec)
          csl = vsext(i,j,ispec)
          rhol = rhoext(i,j,ispec)
        endif

        rho_vp = rhol*cpl
        rho_vs = rhol*csl

        xgamma = - xiz(i,j,ispec) * jacobian(i,j,ispec)
        zgamma = + xix(i,j,ispec) * jacobian(i,j,ispec)
        jacobian1D = sqrt(xgamma**2 + zgamma**2)
        nx = - zgamma / jacobian1D
        nz = + xgamma / jacobian1D

        weight = jacobian1D * wzgll(j)

        if (P_SV) then
          ! P_SV case
          vx = veloc_elastic(1,iglob) - veloc_horiz
          vz = veloc_elastic(2,iglob) - veloc_vert
          vn = nx*vx+nz*vz
          tx = rho_vp*vn*nx+rho_vs*(vx-vn*nx)
          tz = rho_vp*vn*nz+rho_vs*(vz-vn*nz)
        else
          ! SH case
          vy = veloc_elastic(1,iglob)
          ty = rho_vs*vy
        endif

        if (P_SV) then
          ! P_SV case
          accel_elastic(1,iglob) = accel_elastic(1,iglob) - (tx + traction_x_t0)*weight
          accel_elastic(2,iglob) = accel_elastic(2,iglob) - (tz + traction_z_t0)*weight
        else
          ! SH case
          accel_elastic(1,iglob) = accel_elastic(1,iglob) - ty*weight
        endif

        if (SAVE_FORWARD .and. SIMULATION_TYPE == 1) then
          if (P_SV) then
            ! P-SV waves
            b_absorb_elastic_left(1,j,ib_left(ispecabs),it) = (tx + traction_x_t0)*weight
            b_absorb_elastic_left(2,j,ib_left(ispecabs),it) = (tz + traction_z_t0)*weight
          else
            ! SH (membrane) waves
            b_absorb_elastic_left(1,j,ib_left(ispecabs),it) = ty*weight
          endif
        endif
      enddo
    endif  !  end of left absorbing boundary

    !--- right absorbing boundary
    if (codeabs(IEDGE2,ispecabs)) then
      i = NGLLX
      do j = 1,NGLLZ
        ! Clayton-Engquist condition if elastic
        iglob = ibool(i,j,ispec)

        veloc_horiz = 0._CUSTOM_REAL
        veloc_vert = 0._CUSTOM_REAL
        traction_x_t0 = 0._CUSTOM_REAL
        traction_z_t0 = 0._CUSTOM_REAL

        ! for analytical initial plane wave for Bielak's conditions
        ! left or right edge, horizontal normal vector
        if (add_Bielak_conditions_right .and. initialfield) then
          if (.not. over_critical_angle) then
            call compute_Bielak_conditions(coord,iglob,nglob,it,deltat,dxUx,dxUz,dzUx,dzUz,veloc_horiz,veloc_vert, &
                        x_source(1), z_source(1), A_plane, B_plane, C_plane, anglesource(1), anglesource_refl, &
                        c_inc, c_refl, time_offset,f0_source(1))

            traction_x_t0 = lambdaplus2mu_unrelaxed_elastic * dxUx + lambdal_unrelaxed_elastic * dzUz
            traction_z_t0 = mul_unrelaxed_elastic * (dxUz + dzUx)
          else
            veloc_horiz = v0x_right(count_right,it)
            veloc_vert = v0z_right(count_right,it)
            traction_x_t0 = t0x_right(count_right,it)
            traction_z_t0 = t0z_right(count_right,it)
            count_right = count_right+1
          endif
        endif

        ! external velocity model
        if (assign_external_model) then
          cpl = vpext(i,j,ispec)
          csl = vsext(i,j,ispec)
          rhol = rhoext(i,j,ispec)
        endif

        rho_vp = rhol*cpl
        rho_vs = rhol*csl

        xgamma = - xiz(i,j,ispec) * jacobian(i,j,ispec)
        zgamma = + xix(i,j,ispec) * jacobian(i,j,ispec)
        jacobian1D = sqrt(xgamma**2 + zgamma**2)
        nx = + zgamma / jacobian1D
        nz = - xgamma / jacobian1D

        weight = jacobian1D * wzgll(j)

        if (P_SV) then
          ! P_SV case
          vx = veloc_elastic(1,iglob) - veloc_horiz
          vz = veloc_elastic(2,iglob) - veloc_vert
          vn = nx*vx+nz*vz
          tx = rho_vp*vn*nx+rho_vs*(vx-vn*nx)
          tz = rho_vp*vn*nz+rho_vs*(vz-vn*nz)
        else
          ! SH case
          vy = veloc_elastic(1,iglob)
          ty = rho_vs*vy
        endif

        if (P_SV) then
          ! P_SV case
          accel_elastic(1,iglob) = accel_elastic(1,iglob) - (tx - traction_x_t0)*weight
          accel_elastic(2,iglob) = accel_elastic(2,iglob) - (tz - traction_z_t0)*weight
        else
          ! SH case
          accel_elastic(1,iglob) = accel_elastic(1,iglob) - ty*weight
        endif

        if (SAVE_FORWARD .and. SIMULATION_TYPE == 1) then
          if (P_SV) then
            ! P-SV waves
            b_absorb_elastic_right(1,j,ib_right(ispecabs),it) = (tx - traction_x_t0)*weight
            b_absorb_elastic_right(2,j,ib_right(ispecabs),it) = (tz - traction_z_t0)*weight
          else
            ! SH (membrane) waves
            b_absorb_elastic_right(1,j,ib_right(ispecabs),it) = ty*weight
          endif
        endif
      enddo
    endif  !  end of right absorbing boundary

    !--- bottom absorbing boundary
    if (codeabs(IEDGE1,ispecabs)) then
      j = 1
!! DK DK not needed           ! exclude corners to make sure there is no contradiction on the normal
      ibegin = 1
      iend = NGLLX
!! DK DK not needed           if (codeabs(IEDGE4,ispecabs)) ibegin = 2
!! DK DK not needed           if (codeabs(IEDGE2,ispecabs)) iend = NGLLX-1
      do i = ibegin,iend
        ! Clayton-Engquist condition if elastic
        iglob = ibool(i,j,ispec)

        veloc_horiz = 0._CUSTOM_REAL
        veloc_vert = 0._CUSTOM_REAL
        traction_x_t0 = 0._CUSTOM_REAL
        traction_z_t0 = 0._CUSTOM_REAL

        ! for analytical initial plane wave for Bielak's conditions
        ! top or bottom edge, vertical normal vector
        if (add_Bielak_conditions_bottom .and. initialfield) then
          if (.not. over_critical_angle) then
            call compute_Bielak_conditions(coord,iglob,nglob,it,deltat,dxUx,dxUz,dzUx,dzUz,veloc_horiz,veloc_vert, &
                        x_source(1), z_source(1), A_plane, B_plane, C_plane, anglesource(1), anglesource_refl, &
                        c_inc, c_refl, time_offset,f0_source(1))

            traction_x_t0 = mul_unrelaxed_elastic * (dxUz + dzUx)
            traction_z_t0 = lambdal_unrelaxed_elastic * dxUx + lambdaplus2mu_unrelaxed_elastic * dzUz
          else
            veloc_horiz = v0x_bot(count_bottom,it)
            veloc_vert = v0z_bot(count_bottom,it)
            traction_x_t0 = t0x_bot(count_bottom,it)
            traction_z_t0 = t0z_bot(count_bottom,it)
            count_bottom = count_bottom+1
          endif
        endif

        ! external velocity model
        if (assign_external_model) then
          cpl = vpext(i,j,ispec)
          csl = vsext(i,j,ispec)
          rhol = rhoext(i,j,ispec)
        endif

        rho_vp = rhol*cpl
        rho_vs = rhol*csl

        xxi = + gammaz(i,j,ispec) * jacobian(i,j,ispec)
        zxi = - gammax(i,j,ispec) * jacobian(i,j,ispec)
        jacobian1D = sqrt(xxi**2 + zxi**2)
        nx = + zxi / jacobian1D
        nz = - xxi / jacobian1D

        weight = jacobian1D * wxgll(i)

        if (P_SV) then
          ! P_SV case
          vx = veloc_elastic(1,iglob) - veloc_horiz
          vz = veloc_elastic(2,iglob) - veloc_vert
          vn = nx*vx+nz*vz
          tx = rho_vp*vn*nx+rho_vs*(vx-vn*nx)
          tz = rho_vp*vn*nz+rho_vs*(vz-vn*nz)
        else
          ! SH case
          vy = veloc_elastic(1,iglob)
          ty = rho_vs*vy
        endif

! exclude corners to make sure there is no contradiction on the normal
! for Stacey absorbing conditions but not for incident plane waves;
! thus subtract nothing i.e. zero in that case
        if ((codeabs_corner(1,ispecabs) .and. i == 1) .or. (codeabs_corner(2,ispecabs) .and. i == NGLLX)) then
          tx = 0._CUSTOM_REAL
          ty = 0._CUSTOM_REAL
          tz = 0._CUSTOM_REAL
        endif

        if (P_SV) then
          ! P_SV case
          accel_elastic(1,iglob) = accel_elastic(1,iglob) - (tx + traction_x_t0)*weight
          accel_elastic(2,iglob) = accel_elastic(2,iglob) - (tz + traction_z_t0)*weight
        else
          ! SH case
          accel_elastic(1,iglob) = accel_elastic(1,iglob) - ty*weight
        endif

        if (SAVE_FORWARD .and. SIMULATION_TYPE == 1) then
          if (P_SV) then
            ! P-SV waves
            b_absorb_elastic_bottom(1,i,ib_bottom(ispecabs),it) = (tx + traction_x_t0)*weight
            b_absorb_elastic_bottom(2,i,ib_bottom(ispecabs),it) = (tz + traction_z_t0)*weight
          else
            ! SH (membrane) waves
            b_absorb_elastic_bottom(1,i,ib_bottom(ispecabs),it) = ty*weight
          endif
        endif
      enddo
    endif  !  end of bottom absorbing boundary

    !--- top absorbing boundary
    if (codeabs(IEDGE3,ispecabs)) then
      j = NGLLZ
!! DK DK not needed           ! exclude corners to make sure there is no contradiction on the normal
      ibegin = 1
      iend = NGLLX
!! DK DK not needed           if (codeabs(IEDGE4,ispecabs)) ibegin = 2
!! DK DK not needed           if (codeabs(IEDGE2,ispecabs)) iend = NGLLX-1
      do i = ibegin,iend
        ! Clayton-Engquist condition if elastic
        iglob = ibool(i,j,ispec)

        veloc_horiz = 0._CUSTOM_REAL
        veloc_vert = 0._CUSTOM_REAL
        traction_x_t0 = 0._CUSTOM_REAL
        traction_z_t0 = 0._CUSTOM_REAL

        ! for analytical initial plane wave for Bielak's conditions
        ! top or bottom edge, vertical normal vector
        if (add_Bielak_conditions_top .and. initialfield) then
        ! at the top there is no test for whether we are above the critical angle
        ! because a critical angle can only exist when the top edge is a free surface, not in an infinite medium
          call compute_Bielak_conditions(coord,iglob,nglob,it,deltat,dxUx,dxUz,dzUx,dzUz,veloc_horiz,veloc_vert, &
                      x_source(1), z_source(1), A_plane, B_plane, C_plane, anglesource(1), anglesource_refl, &
                      c_inc, c_refl, time_offset,f0_source(1))

          traction_x_t0 = mul_unrelaxed_elastic * (dxUz + dzUx)
          traction_z_t0 = lambdal_unrelaxed_elastic * dxUx + lambdaplus2mu_unrelaxed_elastic * dzUz
        endif

        ! external velocity model
        if (assign_external_model) then
          cpl = vpext(i,j,ispec)
          csl = vsext(i,j,ispec)
          rhol = rhoext(i,j,ispec)
        endif

        rho_vp = rhol*cpl
        rho_vs = rhol*csl

        xxi = + gammaz(i,j,ispec) * jacobian(i,j,ispec)
        zxi = - gammax(i,j,ispec) * jacobian(i,j,ispec)
        jacobian1D = sqrt(xxi**2 + zxi**2)
        nx = - zxi / jacobian1D
        nz = + xxi / jacobian1D

        weight = jacobian1D * wxgll(i)

        if (P_SV) then
          ! P_SV case
          vx = veloc_elastic(1,iglob) - veloc_horiz
          vz = veloc_elastic(2,iglob) - veloc_vert
          vn = nx*vx+nz*vz
          tx = rho_vp*vn*nx+rho_vs*(vx-vn*nx)
          tz = rho_vp*vn*nz+rho_vs*(vz-vn*nz)
        else
          ! SH case
          vy = veloc_elastic(1,iglob)
          ty = rho_vs*vy
        endif

! exclude corners to make sure there is no contradiction on the normal
! for Stacey absorbing conditions but not for incident plane waves;
! thus subtract nothing i.e. zero in that case
        if ((codeabs_corner(3,ispecabs) .and. i == 1) .or. (codeabs_corner(4,ispecabs) .and. i == NGLLX)) then
          tx = 0._CUSTOM_REAL
          ty = 0._CUSTOM_REAL
          tz = 0._CUSTOM_REAL
        endif

        if (P_SV) then
          ! P_SV case
          accel_elastic(1,iglob) = accel_elastic(1,iglob) - (tx - traction_x_t0)*weight
          accel_elastic(2,iglob) = accel_elastic(2,iglob) - (tz - traction_z_t0)*weight
        else
          ! SH case
          accel_elastic(1,iglob) = accel_elastic(1,iglob) - ty*weight
        endif

        if (SAVE_FORWARD .and. SIMULATION_TYPE == 1) then
          if (P_SV) then
            ! P-SV waves
            b_absorb_elastic_top(1,i,ib_top(ispecabs),it) = (tx - traction_x_t0)*weight
            b_absorb_elastic_top(2,i,ib_top(ispecabs),it) = (tz - traction_z_t0)*weight
          else
            ! SH (membrane) waves
            b_absorb_elastic_top(1,i,ib_top(ispecabs),it) = ty*weight
          endif
        endif
      enddo
    endif  !  end of top absorbing boundary

  enddo

  end subroutine compute_stacey_elastic

!
!------------------------------------------------------------------------------------------
!

  subroutine compute_stacey_elastic_backward(b_accel_elastic)

! absorbing boundaries
!
! Clayton-Engquist condition if elastic

  use constants,only: CUSTOM_REAL,NGLLX,NGLLZ,NDIM, &
    ZERO,ONE,TWO,FOUR_THIRDS,IEDGE1,IEDGE2,IEDGE3,IEDGE4

  use specfem_par, only: nglob,any_elastic,ibool,ispec_is_elastic, &
                         NSTEP,it, &
                         nelemabs,numabs,codeabs, &
                         b_absorb_elastic_left,b_absorb_elastic_right, &
                         b_absorb_elastic_bottom,b_absorb_elastic_top, &
                         ib_left,ib_right,ib_bottom,ib_top, &
                         STACEY_BOUNDARY_CONDITIONS,P_SV

  implicit none

  real(kind=CUSTOM_REAL), dimension(NDIM,nglob),intent(inout) :: b_accel_elastic

  ! local parameters
  integer :: ispecabs,ispec,i,j,iglob,it_tmp
  integer :: ibegin,iend

  ! checks if anything to do
  if (.not. STACEY_BOUNDARY_CONDITIONS) return
  if (.not. any_elastic) return

  ! time increment index
  it_tmp = NSTEP - it + 1

  ! Clayton-Engquist condition if elastic
  do ispecabs = 1,nelemabs

    ispec = numabs(ispecabs)
    if (.not. ispec_is_elastic(ispec) ) cycle

    !--- left absorbing boundary
    if (codeabs(IEDGE4,ispecabs)) then
      i = 1
      do j = 1,NGLLZ
        ! Clayton-Engquist condition if elastic
        iglob = ibool(i,j,ispec)
        if (P_SV) then
          ! P-SV waves
          b_accel_elastic(1,iglob) = b_accel_elastic(1,iglob) - b_absorb_elastic_left(1,j,ib_left(ispecabs),it_tmp)
          b_accel_elastic(2,iglob) = b_accel_elastic(2,iglob) - b_absorb_elastic_left(2,j,ib_left(ispecabs),it_tmp)
        else
          ! SH (membrane) waves
          b_accel_elastic(1,iglob) = b_accel_elastic(1,iglob) - b_absorb_elastic_left(1,j,ib_left(ispecabs),it_tmp)
        endif
      enddo
    endif  !  end of left absorbing boundary

    !--- right absorbing boundary
    if (codeabs(IEDGE2,ispecabs)) then
      i = NGLLX
      do j = 1,NGLLZ
        iglob = ibool(i,j,ispec)
        if (P_SV) then
          ! P-SV waves
          b_accel_elastic(1,iglob) = b_accel_elastic(1,iglob) - b_absorb_elastic_right(1,j,ib_right(ispecabs),it_tmp)
          b_accel_elastic(2,iglob) = b_accel_elastic(2,iglob) - b_absorb_elastic_right(2,j,ib_right(ispecabs),it_tmp)
        else
          ! SH (membrane) waves
          b_accel_elastic(1,iglob) = b_accel_elastic(1,iglob) - b_absorb_elastic_right(1,j,ib_right(ispecabs),it_tmp)
        endif
      enddo
    endif  !  end of right absorbing boundary

    !--- bottom absorbing boundary
    if (codeabs(IEDGE1,ispecabs)) then
      j = 1
!! DK DK not needed           ! exclude corners to make sure there is no contradiction on the normal
      ibegin = 1
      iend = NGLLX
!! DK DK not needed           if (codeabs(IEDGE4,ispecabs)) ibegin = 2
!! DK DK not needed           if (codeabs(IEDGE2,ispecabs)) iend = NGLLX-1
      do i = ibegin,iend
        iglob = ibool(i,j,ispec)
        if (P_SV) then
          ! P-SV waves
          b_accel_elastic(1,iglob) = b_accel_elastic(1,iglob) - b_absorb_elastic_bottom(1,i,ib_bottom(ispecabs),it_tmp)
          b_accel_elastic(2,iglob) = b_accel_elastic(2,iglob) - b_absorb_elastic_bottom(2,i,ib_bottom(ispecabs),it_tmp)
        else
          ! SH (membrane) waves
          b_accel_elastic(1,iglob) = b_accel_elastic(1,iglob) - b_absorb_elastic_bottom(1,i,ib_bottom(ispecabs),it_tmp)
        endif
      enddo
    endif  !  end of bottom absorbing boundary

    !--- top absorbing boundary
    if (codeabs(IEDGE3,ispecabs)) then
      j = NGLLZ
!! DK DK not needed           ! exclude corners to make sure there is no contradiction on the normal
      ibegin = 1
      iend = NGLLX
!! DK DK not needed           if (codeabs(IEDGE4,ispecabs)) ibegin = 2
!! DK DK not needed           if (codeabs(IEDGE2,ispecabs)) iend = NGLLX-1
      do i = ibegin,iend
        iglob = ibool(i,j,ispec)
        if (P_SV) then
          ! P-SV waves
          b_accel_elastic(1,iglob) = b_accel_elastic(1,iglob) - b_absorb_elastic_top(1,i,ib_top(ispecabs),it_tmp)
          b_accel_elastic(2,iglob) = b_accel_elastic(2,iglob) - b_absorb_elastic_top(2,i,ib_top(ispecabs),it_tmp)
        else
          ! SH (membrane) waves
          b_accel_elastic(1,iglob) = b_accel_elastic(1,iglob) - b_absorb_elastic_top(1,i,ib_top(ispecabs),it_tmp)
        endif
      enddo
    endif  !  end of top absorbing boundary

  enddo

  end subroutine compute_stacey_elastic_backward


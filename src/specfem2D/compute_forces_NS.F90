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

! TODO: Unused, decide what to do with this.

  subroutine compute_forces_NS(dt_rhoveloc_acoustic, dt_density, dt_rhoenergy, &
        rhoveloc_acoustic, rhoenergy, density_p)

  ! compute forces for the acoustic elements
  use constants,only: CUSTOM_REAL,NGLLX,NGLLZ,NDIM

  use specfem_par, only: nglob,nspec, &
                         ibool, jacobian,pext,rhoext, &
                         hprimewgll_xx,hprimewgll_zz,wxgll,wzgll, &
                         ispec_is_acoustic, gamma_euler, &
                         gammax, gammaz, xix, xiz

  implicit none

  real(kind=CUSTOM_REAL), dimension(NDIM,nglob) :: dt_rhoveloc_acoustic, rhoveloc_acoustic
  real(kind=CUSTOM_REAL), dimension(nglob) :: dt_density, dt_rhoenergy, rhoenergy, density_p

  !---
  !--- local variables
  !---

  integer :: ispec,i,j,k,iglob

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: &
        tempx1_accel,tempx2_accel,tempz1_accel,tempz2_accel, &
        temp1_density, temp2_density, temp1_energy, temp2_energy

  ! Jacobian matrix and determinant
  real(kind=CUSTOM_REAL) :: xixl,xizl,gammaxl,gammazl,jacobianl

  ! material properties of the acoustic medium
  real(kind=CUSTOM_REAL) :: rhol, pl
  
  real(kind=CUSTOM_REAL) :: sigma_xx, sigma_xz, sigma_zx, sigma_zz, &
        rho_total, veloc_x, veloc_z, e0, rhoe_total, pressure_total, pressure

  integer :: ifirstelem,ilastelem
  
  dt_density           = 0._CUSTOM_REAL
  dt_rhoveloc_acoustic = 0._CUSTOM_REAL
  dt_rhoenergy         = 0._CUSTOM_REAL
  
  ifirstelem = 1
  ilastelem = nspec
  
  ! loop over spectral elements
  do ispec = ifirstelem,ilastelem
  
    tempx1_accel(:,:) = 0._CUSTOM_REAL
    tempx2_accel(:,:) = 0._CUSTOM_REAL
    tempz1_accel(:,:) = 0._CUSTOM_REAL
    tempz2_accel(:,:) = 0._CUSTOM_REAL
    temp1_density(:,:) = 0._CUSTOM_REAL
    temp2_density(:,:) = 0._CUSTOM_REAL
    temp1_energy(:,:) = 0._CUSTOM_REAL
    temp2_energy(:,:) = 0._CUSTOM_REAL

    !--- elastic spectral element
    if (ispec_is_acoustic(ispec)) then

      ! first double loop over GLL points to compute and store gradients
      do j = 1,NGLLZ
        do i = 1,NGLLX
        
          !--- if external medium, get elastic parameters of current grid point
          rhol = 1.!rhoext(i,j,ispec)
          pl   = 1.!pext(i,j,ispec)
          if(.false.) WRITE(*,*) rhoext, pext
          
          jacobianl = jacobian(i,j,ispec)
          xixl = xix(i,j,ispec)
          xizl = xiz(i,j,ispec)
          gammaxl = gammax(i,j,ispec)
          gammazl = gammaz(i,j,ispec)
          
          ! Total density contributions
          rho_total      = rhol + density_p(ibool(i,j,ispec))
          
          ! Non linear velocity Perturbation
          veloc_x = rhoveloc_acoustic(1,ibool(i,j,ispec))/rho_total
          veloc_z = rhoveloc_acoustic(2,ibool(i,j,ispec))/rho_total
          
          ! Energy Background
          e0 = pl/(gamma_euler - 1.)! &
                !+ (rhol/2.)*( veloc_x**2 + veloc_z**2 )
          
          ! Total rho*energy contribution
          rhoe_total = rhoenergy(ibool(i,j,ispec)) + rho_total*e0
          
          ! Total pressure contribution
          pressure_total = (gamma_euler - 1.)*(rhoe_total &
                - (rho_total/2.)*( veloc_x**2 + veloc_z**2 ) )

          ! Perturbation pressure contribution
          !pressure       = pressure_total - pl
          pressure = (gamma_euler - 1.)*(rhoenergy(ibool(i,j,ispec)) &
                - (density_p(ibool(i,j,ispec))/2.)*( veloc_x**2 + veloc_z**2 ) )
          
          ! Mass conservation
          sigma_xx = rhoveloc_acoustic(1,ibool(i,j,ispec))
          sigma_zx = rhoveloc_acoustic(2,ibool(i,j,ispec))
          temp1_density(i,j) = wzgll(j)*jacobianl*(sigma_xx*xixl+sigma_zx*xizl)
          temp2_density(i,j) = wxgll(i)*jacobianl*(sigma_xx*gammaxl+sigma_zx*gammazl)

          ! Momentum
          sigma_xx = rho_total*veloc_x**2 + pressure
          sigma_xz = rho_total*veloc_x*veloc_z
          sigma_zz = rho_total*veloc_z**2 + pressure
          sigma_zx = rho_total*veloc_x*veloc_z
          
          tempx1_accel(i,j) = wzgll(j)*jacobianl*(sigma_xx*xixl+sigma_zx*xizl) ! this goes to accel_x
          tempz1_accel(i,j) = wzgll(j)*jacobianl*(sigma_xz*xixl+sigma_zz*xizl) ! this goes to accel_z

          tempx2_accel(i,j) = wxgll(i)*jacobianl*(sigma_xx*gammaxl+sigma_zx*gammazl) ! this goes to accel_x
          tempz2_accel(i,j) = wxgll(i)*jacobianl*(sigma_xz*gammaxl+sigma_zz*gammazl) ! this goes to accel_z
          
          ! Energy conservation
          sigma_xx = (rhoe_total + pressure_total)*veloc_x
          sigma_zx = (rhoe_total + pressure_total)*veloc_z
          temp1_energy(i,j) = wzgll(j)*jacobianl*(sigma_xx*xixl+sigma_zx*xizl)
          temp2_energy(i,j) = wxgll(i)*jacobianl*(sigma_xx*gammaxl+sigma_zx*gammazl)

        enddo
      enddo  ! end of the loops on the collocation points i,j

      !
      ! second double-loop over GLL to compute all the terms
      !
      do j = 1,NGLLZ
        do i = 1,NGLLX
          iglob = ibool(i,j,ispec)
          ! along x direction and z direction
          ! and assemble the contributions
          ! we can merge the two loops because NGLLX == NGLLZ

            !if AXISYM == false
            do k = 1,NGLLX
              dt_density(iglob) = dt_density(iglob) + &
                (temp1_density(k,j)*hprimewgll_xx(k,i) + temp2_density(i,k)*hprimewgll_zz(k,j))
                
              dt_rhoveloc_acoustic(1,iglob) = dt_rhoveloc_acoustic(1,iglob) + &
                (tempx1_accel(k,j)*hprimewgll_xx(k,i) + tempx2_accel(i,k)*hprimewgll_zz(k,j))
                
              dt_rhoveloc_acoustic(2,iglob) = dt_rhoveloc_acoustic(2,iglob) + &
                (tempz1_accel(k,j)*hprimewgll_xx(k,i) + tempz2_accel(i,k)*hprimewgll_zz(k,j))
                
              dt_rhoenergy(iglob) = dt_rhoenergy(iglob) + &
                (temp1_energy(k,j)*hprimewgll_xx(k,i) + temp2_energy(i,k)*hprimewgll_zz(k,j))
            enddo

        enddo
      enddo ! second loop over the GLL points
   
      endif ! end of test if acoustic element
   
   enddo ! end of loop over all elements
   
   WRITE(*,*) "MAXVAL", maxval(dt_rhoveloc_acoustic(2,:)), minval(dt_rhoveloc_acoustic(2,:))

  end subroutine compute_forces_NS


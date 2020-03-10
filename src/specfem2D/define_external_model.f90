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

  subroutine define_external_model_dummy(coord,material_element,ibool, &
                                         rho,vp,vs,QKappa_attenuation,Qmu_attenuation,gravity,Nsq, &
                                         c11,c13,c15,c33,c35,c55,c12,c23,c25,nspec,nglob)

  use constants,only: CUSTOM_REAL,NGLLX,NGLLZ,NDIM,IMAIN

  use specfem_par, only: poroelastcoef,density,kmato

  implicit none

! -------------------------------------------------------------------------------------
! Dummy example of this routine, to be used as a template that users can modify.
! To use it you will need to rename it as define_external_model() (i.e., get rid of "_dummy")
! and suppress the existing define_external_model() routine provided below for the AK135F global Earth model.
! -------------------------------------------------------------------------------------

! users can modify this routine to assign any different external model (rho, vp, vs)
! based on the x and y coordinates of that grid point and the material number of the region it belongs to

  integer, intent(in) :: nspec,nglob

  double precision, dimension(NDIM,nglob), intent(in) :: coord

  integer, dimension(nspec), intent(in) :: material_element

  integer, dimension(NGLLX,NGLLZ,nspec), intent(in) :: ibool

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ,nspec), intent(out) :: rho,vp,vs,QKappa_attenuation,Qmu_attenuation,gravity,Nsq, &
                                                                 c11,c15,c13,c33,c35,c55,c12,c23,c25

  integer :: i,j,ispec,iglob

  double precision :: x,z

! dummy routine here, just to demonstrate how the model can be assigned
! and how such a routine can be written

! remove gravity
! leave these arrays here even if you do not assign them to use them because they need to be cleared
  gravity(:,:,:) = 0.d0
  Nsq(:,:,:) = 0.d0

! loop on all the elements of the mesh, and inside each element loop on all the GLL points
  do ispec = 1,nspec
    do j = 1,NGLLZ
      do i = 1,NGLLX

        iglob = ibool(i,j,ispec)

        x = coord(1,iglob)
        z = coord(2,iglob)

        if (material_element(ispec) == 1 .or. x < 1700.d0 .or. z >= 2300.d0) then
          rho(i,j,ispec) = 2000.d0 + 10.d0*x - 32.17d0*z**3
          vp(i,j,ispec) = 3000.d0 - 7.12d0*cos(x)
          vs(i,j,ispec) = vp(i,j,ispec) / sqrt(3.d0)
          QKappa_attenuation(i,j,ispec) = 9999. ! this means no attenuation
          Qmu_attenuation(i,j,ispec)    = 9999. ! this means no attenuation
          c11(i,j,ispec) = 169.d9
          c13(i,j,ispec) = 122.d9
          c15(i,j,ispec) = 0.d0
          c33(i,j,ispec) = c11(i,j,ispec)
          c35(i,j,ispec) = 0.d0
          c55(i,j,ispec) = 75.3d9
          c12(i,j,ispec) = 0.d0
          c23(i,j,ispec) = 0.d0
          c25(i,j,ispec) = 0.d0

        else if (material_element(ispec) == 2) then
          rho(i,j,ispec) = 2500.d0
          vp(i,j,ispec) = 3600.d0
          vs(i,j,ispec) = vp(i,j,ispec) / 2.d0
          QKappa_attenuation(i,j,ispec) = 120.
          Qmu_attenuation(i,j,ispec) = 120.
          c11(i,j,ispec) = 0.d0   ! this means no anisotropy
          c13(i,j,ispec) = 0.d0
          c15(i,j,ispec) = 0.d0
          c33(i,j,ispec) = 0.d0
          c35(i,j,ispec) = 0.d0
          c55(i,j,ispec) = 0.d0
          c12(i,j,ispec) = 0.d0
          c23(i,j,ispec) = 0.d0
          c25(i,j,ispec) = 0.d0

        else
          ! this case should not occur
          write(IMAIN,*) 'Error: invalid flag number in external model is equal to ',material_element(ispec)
          stop 'wrong flag number in external model; exiting...'
        endif

        !! AB AB Do not forget these 3 lines otherwise PML may not work !!
        density(1,kmato(ispec)) = rho(i,j,ispec)
        poroelastcoef(3,1,kmato(ispec)) = rho(i,j,ispec) * vp(i,j,ispec) * vp(i,j,ispec)
        poroelastcoef(2,1,kmato(ispec)) =  rho(i,j,ispec) * vs(i,j,ispec) * vs(i,j,ispec)
        !! AB AB Do not forget these 3 lines otherwise PML may not work !!

      enddo
    enddo
  enddo

  end subroutine define_external_model_dummy


!========================================================================
!
! another example below, to define the AK135F global Earth model
!
!========================================================================
  subroutine define_external_model(coord,material_element,ibool, &
                                   rho,vp,vs,QKappa_attenuation,Qmu_attenuation,gravity,Nsq, &
                                   c11,c13,c15,c33,c35,c55,c12,c23,c25,nspec,nglob)

  use constants,only: CUSTOM_REAL,NGLLX,NGLLZ,NDIM,IMAIN

  implicit none

!--------------------------------------------------------------------------------------------------
!
!          taken from S p e c f e m 3 D  G l o b e  V e r s i o n  5 . 1
!          -------------------------------------------------------------

! Modified AK135 model:
!
! Spherically symmetric isotropic AK135 model (Kennett et al., 1995).
! modified to use the density and Q attenuation models of Montagner and Kennett (1995).
! That modified model is traditionally called AK135-F,
! see http://rses.anu.edu.au/seismology/ak135/ak135f.html for more details.
! As we do not want to use the 300 m-thick mud layer from that model nor the ocean layer,
! above the d120 discontinuity we switch back to the classical AK135 model of Kennett et al. (1995),
! i.e., we use AK135-F below and AK135 above.

! B. L. N. Kennett, E. R. Engdahl and R. Buland,
! Constraints on seismic velocities in the Earth from traveltimes,
! Geophysical Journal International, volume 122, issue 1, pages 108-124 (1995),
! DOI: 10.1111/j.1365-246X.1995.tb03540.x

! J. P. Montagner and B. L. N. Kennett,
! How to reconcile body-wave and normal-mode reference Earth models?,
! Geophysical Journal International, volume 122, issue 1, pages 229-248 (1995)

!! DK DK values below entirely checked and fixed by Dimitri Komatitsch in December 2012.

!--------------------------------------------------------------------------------------------------

  integer, intent(in) :: nspec,nglob

  double precision, dimension(NDIM,nglob), intent(in) :: coord

  integer, dimension(nspec), intent(in) :: material_element

  integer, dimension(NGLLX,NGLLZ,nspec), intent(in) :: ibool

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ,nspec), intent(out) :: rho,vp,vs,QKappa_attenuation,Qmu_attenuation,gravity,Nsq, &
                                                                 c11,c15,c13,c33,c35,c55,c12,c23,c25

! number of layers in ak135-f
  integer, parameter :: NR_AK135F_NO_MUD = 136

  double precision, dimension(NR_AK135F_NO_MUD) :: radius_ak135
  double precision, dimension(NR_AK135F_NO_MUD) :: density_ak135
  double precision, dimension(NR_AK135F_NO_MUD) :: vp_ak135
  double precision, dimension(NR_AK135F_NO_MUD) :: vs_ak135
  double precision, dimension(NR_AK135F_NO_MUD) :: Qkappa_ak135
  double precision, dimension(NR_AK135F_NO_MUD) :: Qmu_ak135

! region flags to assign the AK135F_NO_MUD Earth model
  integer, parameter :: IREGION_MANTLE_CRUST_ABOVE_d670 = 1
  integer, parameter :: IREGION_MANTLE_BELOW_d670 = 2
  integer, parameter :: IREGION_OUTER_CORE = 3
  integer, parameter :: IREGION_INNER_CORE = 4

  integer :: i,j,ispec,iglob,ii

  double precision :: x,z,r,frac

! remove gravity
! leave these arrays here even if you do not assign them to use them because they need to be cleared
  gravity(:,:,:) = 0.d0
  Nsq(:,:,:)     = 0.d0

! define all the values in the model once and for all

  radius_ak135(  1) =  0.000000000000000E+000
  radius_ak135(  2) =   50710.0000000000
  radius_ak135(  3) =   101430.000000000
  radius_ak135(  4) =   152140.000000000
  radius_ak135(  5) =   202850.000000000
  radius_ak135(  6) =   253560.000000000
  radius_ak135(  7) =   304280.000000000
  radius_ak135(  8) =   354990.000000000
  radius_ak135(  9) =   405700.000000000
  radius_ak135( 10) =   456410.000000000
  radius_ak135( 11) =   507130.000000000
  radius_ak135( 12) =   557840.000000000
  radius_ak135( 13) =   659260.000000000
  radius_ak135( 14) =   710000.000000000
  radius_ak135( 15) =   760690.000000000
  radius_ak135( 16) =   811400.000000000
  radius_ak135( 17) =   862110.000000000
  radius_ak135( 18) =   912830.000000000
  radius_ak135( 19) =   963540.000000000
  radius_ak135( 20) =   1014250.00000000
  radius_ak135( 21) =   1064960.00000000
  radius_ak135( 22) =   1115680.00000000
  radius_ak135( 23) =   1166390.00000000
  radius_ak135( 24) =   1217500.00000000
  radius_ak135( 25) =   1217500.00000000
  radius_ak135( 26) =   1267430.00000000
  radius_ak135( 27) =   1317760.00000000
  radius_ak135( 28) =   1368090.00000000
  radius_ak135( 29) =   1418420.00000000
  radius_ak135( 30) =   1468760.00000000
  radius_ak135( 31) =   1519090.00000000
  radius_ak135( 32) =   1569420.00000000
  radius_ak135( 33) =   1670080.00000000
  radius_ak135( 34) =   1720410.00000000
  radius_ak135( 35) =   1770740.00000000
  radius_ak135( 36) =   1821070.00000000
  radius_ak135( 37) =   1871400.00000000
  radius_ak135( 38) =   1921740.00000000
  radius_ak135( 39) =   1972070.00000000
  radius_ak135( 40) =   2022400.00000000
  radius_ak135( 41) =   2072730.00000000
  radius_ak135( 42) =   2123060.00000000
  radius_ak135( 43) =   2173390.00000000
  radius_ak135( 44) =   2223720.00000000
  radius_ak135( 45) =   2274050.00000000
  radius_ak135( 46) =   2324380.00000000
  radius_ak135( 47) =   2374720.00000000
  radius_ak135( 48) =   2425050.00000000
  radius_ak135( 49) =   2475380.00000000
  radius_ak135( 50) =   2525710.00000000
  radius_ak135( 51) =   2576040.00000000
  radius_ak135( 52) =   2626370.00000000
  radius_ak135( 53) =   2676700.00000000
  radius_ak135( 54) =   2727030.00000000
  radius_ak135( 55) =   2777360.00000000
  radius_ak135( 56) =   2827700.00000000
  radius_ak135( 57) =   2878030.00000000
  radius_ak135( 58) =   2928360.00000000
  radius_ak135( 59) =   2978690.00000000
  radius_ak135( 60) =   3029020.00000000
  radius_ak135( 61) =   3079350.00000000
  radius_ak135( 62) =   3129680.00000000
  radius_ak135( 63) =   3180010.00000000
  radius_ak135( 64) =   3230340.00000000
  radius_ak135( 65) =   3280680.00000000
  radius_ak135( 66) =   3331010.00000000
  radius_ak135( 67) =   3381340.00000000
  radius_ak135( 68) =   3431670.00000000
  radius_ak135( 69) =   3479500.00000000
  radius_ak135( 70) =   3479500.00000000
  radius_ak135( 71) =   3531670.00000000
  radius_ak135( 72) =   3581330.00000000
  radius_ak135( 73) =   3631000.00000000
  radius_ak135( 74) =   3631000.00000000
  radius_ak135( 75) =   3681000.00000000
  radius_ak135( 76) =   3731000.00000000
  radius_ak135( 77) =   3779500.00000000
  radius_ak135( 78) =   3829000.00000000
  radius_ak135( 79) =   3878500.00000000
  radius_ak135( 80) =   3928000.00000000
  radius_ak135( 81) =   3977500.00000000
  radius_ak135( 82) =   4027000.00000000
  radius_ak135( 83) =   4076500.00000000
  radius_ak135( 84) =   4126000.00000000
  radius_ak135( 85) =   4175500.00000000
  radius_ak135( 86) =   4225000.00000000
  radius_ak135( 87) =   4274500.00000000
  radius_ak135( 88) =   4324000.00000000
  radius_ak135( 89) =   4373500.00000000
  radius_ak135( 90) =   4423000.00000000
  radius_ak135( 91) =   4472500.00000000
  radius_ak135( 92) =   4522000.00000000
  radius_ak135( 93) =   4571500.00000000
  radius_ak135( 94) =   4621000.00000000
  radius_ak135( 95) =   4670500.00000000
  radius_ak135( 96) =   4720000.00000000
  radius_ak135( 97) =   4769500.00000000
  radius_ak135( 98) =   4819000.00000000
  radius_ak135( 99) =   4868500.00000000
  radius_ak135(100) =   4918000.00000000
  radius_ak135(101) =   4967500.00000000
  radius_ak135(102) =   5017000.00000000
  radius_ak135(103) =   5066500.00000000
  radius_ak135(104) =   5116000.00000000
  radius_ak135(105) =   5165500.00000000
  radius_ak135(106) =   5215000.00000000
  radius_ak135(107) =   5264500.00000000
  radius_ak135(108) =   5314000.00000000
  radius_ak135(109) =   5363500.00000000
  radius_ak135(110) =   5413000.00000000
  radius_ak135(111) =   5462500.00000000
  radius_ak135(112) =   5512000.00000000
  radius_ak135(113) =   5561500.00000000
  radius_ak135(114) =   5611000.00000000
  radius_ak135(115) =   5661000.00000000
  radius_ak135(116) =   5711000.00000000
  radius_ak135(117) =   5711000.00000000
  radius_ak135(118) =   5761000.00000000
  radius_ak135(119) =   5811000.00000000
  radius_ak135(120) =   5861000.00000000
  radius_ak135(121) =   5911000.00000000
  radius_ak135(122) =   5961000.00000000
  radius_ak135(123) =   5961000.00000000
  radius_ak135(124) =   6011000.00000000
  radius_ak135(125) =   6061000.00000000
  radius_ak135(126) =   6111000.00000000
  radius_ak135(127) =   6161000.00000000
  radius_ak135(128) =   6161000.00000000
  radius_ak135(129) =   6206000.00000000
  radius_ak135(130) =   6251000.00000000
  radius_ak135(131) =   6293500.00000000
  radius_ak135(132) =   6336000.00000000
  radius_ak135(133) =   6336000.00000000
  radius_ak135(134) =   6351000.00000000
  radius_ak135(135) =   6351000.00000000
  radius_ak135(136) =   6371000.00000000

  density_ak135(  1) =   13.0122000000000
  density_ak135(  2) =   13.0117000000000
  density_ak135(  3) =   13.0100000000000
  density_ak135(  4) =   13.0074000000000
  density_ak135(  5) =   13.0036000000000
  density_ak135(  6) =   12.9988000000000
  density_ak135(  7) =   12.9929000000000
  density_ak135(  8) =   12.9859000000000
  density_ak135(  9) =   12.9779000000000
  density_ak135( 10) =   12.9688000000000
  density_ak135( 11) =   12.9586000000000
  density_ak135( 12) =   12.9474000000000
  density_ak135( 13) =   12.9217000000000
  density_ak135( 14) =   12.9070000000000
  density_ak135( 15) =   12.8917000000000
  density_ak135( 16) =   12.8751000000000
  density_ak135( 17) =   12.8574000000000
  density_ak135( 18) =   12.8387000000000
  density_ak135( 19) =   12.8188000000000
  density_ak135( 20) =   12.7980000000000
  density_ak135( 21) =   12.7760000000000
  density_ak135( 22) =   12.7530000000000
  density_ak135( 23) =   12.7289000000000
  density_ak135( 24) =   12.7037000000000
  density_ak135( 25) =   12.1391000000000
  density_ak135( 26) =   12.1133000000000
  density_ak135( 27) =   12.0867000000000
  density_ak135( 28) =   12.0593000000000
  density_ak135( 29) =   12.0311000000000
  density_ak135( 30) =   12.0001000000000
  density_ak135( 31) =   11.9722000000000
  density_ak135( 32) =   11.9414000000000
  density_ak135( 33) =   11.8772000000000
  density_ak135( 34) =   11.8437000000000
  density_ak135( 35) =   11.8092000000000
  density_ak135( 36) =   11.7737000000000
  density_ak135( 37) =   11.7373000000000
  density_ak135( 38) =   11.6998000000000
  density_ak135( 39) =   11.6612000000000
  density_ak135( 40) =   11.6216000000000
  density_ak135( 41) =   11.5809000000000
  density_ak135( 42) =   11.5391000000000
  density_ak135( 43) =   11.4962000000000
  density_ak135( 44) =   11.4521000000000
  density_ak135( 45) =   11.4069000000000
  density_ak135( 46) =   11.3604000000000
  density_ak135( 47) =   11.3127000000000
  density_ak135( 48) =   11.2639000000000
  density_ak135( 49) =   11.2137000000000
  density_ak135( 50) =   11.1623000000000
  density_ak135( 51) =   11.1095000000000
  density_ak135( 52) =   11.0555000000000
  density_ak135( 53) =   11.0001000000000
  density_ak135( 54) =   10.9434000000000
  density_ak135( 55) =   10.8852000000000
  density_ak135( 56) =   10.8257000000000
  density_ak135( 57) =   10.7647000000000
  density_ak135( 58) =   10.7023000000000
  density_ak135( 59) =   10.6385000000000
  density_ak135( 60) =   10.5731000000000
  density_ak135( 61) =   10.5062000000000
  density_ak135( 62) =   10.4378000000000
  density_ak135( 63) =   10.3679000000000
  density_ak135( 64) =   10.2964000000000
  density_ak135( 65) =   10.2233000000000
  density_ak135( 66) =   10.1485000000000
  density_ak135( 67) =   10.0722000000000
  density_ak135( 68) =   9.99420000000000
  density_ak135( 69) =   9.91450000000000
  density_ak135( 70) =   5.77210000000000
  density_ak135( 71) =   5.74580000000000
  density_ak135( 72) =   5.71960000000000
  density_ak135( 73) =   5.69340000000000
  density_ak135( 74) =   5.43870000000000
  density_ak135( 75) =   5.41760000000000
  density_ak135( 76) =   5.39620000000000
  density_ak135( 77) =   5.37480000000000
  density_ak135( 78) =   5.35310000000000
  density_ak135( 79) =   5.33130000000000
  density_ak135( 80) =   5.30920000000000
  density_ak135( 81) =   5.28700000000000
  density_ak135( 82) =   5.26460000000000
  density_ak135( 83) =   5.24200000000000
  density_ak135( 84) =   5.21920000000000
  density_ak135( 85) =   5.19630000000000
  density_ak135( 86) =   5.17320000000000
  density_ak135( 87) =   5.14990000000000
  density_ak135( 88) =   5.12640000000000
  density_ak135( 89) =   5.10270000000000
  density_ak135( 90) =   5.07890000000000
  density_ak135( 91) =   5.05480000000000
  density_ak135( 92) =   5.03060000000000
  density_ak135( 93) =   5.00620000000000
  density_ak135( 94) =   4.98170000000000
  density_ak135( 95) =   4.95700000000000
  density_ak135( 96) =   4.93210000000000
  density_ak135( 97) =   4.90690000000000
  density_ak135( 98) =   4.88170000000000
  density_ak135( 99) =   4.85620000000000
  density_ak135(100) =   4.83070000000000
  density_ak135(101) =   4.80500000000000
  density_ak135(102) =   4.77900000000000
  density_ak135(103) =   4.75280000000000
  density_ak135(104) =   4.72660000000000
  density_ak135(105) =   4.70010000000000
  density_ak135(106) =   4.67350000000000
  density_ak135(107) =   4.64670000000000
  density_ak135(108) =   4.61980000000000
  density_ak135(109) =   4.59260000000000
  density_ak135(110) =   4.56540000000000
  density_ak135(111) =   4.51620000000000
  density_ak135(112) =   4.46500000000000
  density_ak135(113) =   4.41180000000000
  density_ak135(114) =   4.35650000000000
  density_ak135(115) =   4.29860000000000
  density_ak135(116) =   4.23870000000000
  density_ak135(117) =   3.92010000000000
  density_ak135(118) =   3.92060000000000
  density_ak135(119) =   3.92180000000000
  density_ak135(120) =   3.92330000000000
  density_ak135(121) =   3.92730000000000
  density_ak135(122) =   3.93170000000000
  density_ak135(123) =   3.50680000000000
  density_ak135(124) =   3.45770000000000
  density_ak135(125) =   3.41100000000000
  density_ak135(126) =   3.36630000000000
  density_ak135(127) =   3.32430000000000
  density_ak135(128) =   3.32430000000000
  density_ak135(129) =   3.37110000000000
  density_ak135(130) =   3.42680000000000
  density_ak135(131) =   3.34500000000000
  density_ak135(132) =   3.32000000000000
  density_ak135(133) =   2.92000000000000
  density_ak135(134) =   2.92000000000000
  density_ak135(135) =   2.72000000000000
  density_ak135(136) =   2.72000000000000

  vp_ak135(  1) =   11.2622000000000
  vp_ak135(  2) =   11.2618000000000
  vp_ak135(  3) =   11.2606000000000
  vp_ak135(  4) =   11.2586000000000
  vp_ak135(  5) =   11.2557000000000
  vp_ak135(  6) =   11.2521000000000
  vp_ak135(  7) =   11.2477000000000
  vp_ak135(  8) =   11.2424000000000
  vp_ak135(  9) =   11.2364000000000
  vp_ak135( 10) =   11.2295000000000
  vp_ak135( 11) =   11.2219000000000
  vp_ak135( 12) =   11.2134000000000
  vp_ak135( 13) =   11.1941000000000
  vp_ak135( 14) =   11.1830000000000
  vp_ak135( 15) =   11.1715000000000
  vp_ak135( 16) =   11.1590000000000
  vp_ak135( 17) =   11.1457000000000
  vp_ak135( 18) =   11.1316000000000
  vp_ak135( 19) =   11.1166000000000
  vp_ak135( 20) =   11.0983000000000
  vp_ak135( 21) =   11.0850000000000
  vp_ak135( 22) =   11.0718000000000
  vp_ak135( 23) =   11.0585000000000
  vp_ak135( 24) =   11.0427000000000
  vp_ak135( 25) =   10.2890000000000
  vp_ak135( 26) =   10.2854000000000
  vp_ak135( 27) =   10.2745000000000
  vp_ak135( 28) =   10.2565000000000
  vp_ak135( 29) =   10.2329000000000
  vp_ak135( 30) =   10.2049000000000
  vp_ak135( 31) =   10.1739000000000
  vp_ak135( 32) =   10.1415000000000
  vp_ak135( 33) =   10.0768000000000
  vp_ak135( 34) =   10.0439000000000
  vp_ak135( 35) =   10.0103000000000
  vp_ak135( 36) =   9.97610000000000
  vp_ak135( 37) =   9.94100000000000
  vp_ak135( 38) =   9.90510000000000
  vp_ak135( 39) =   9.86820000000000
  vp_ak135( 40) =   9.83040000000000
  vp_ak135( 41) =   9.79140000000000
  vp_ak135( 42) =   9.75130000000000
  vp_ak135( 43) =   9.71000000000000
  vp_ak135( 44) =   9.66730000000000
  vp_ak135( 45) =   9.62320000000000
  vp_ak135( 46) =   9.57770000000000
  vp_ak135( 47) =   9.53060000000000
  vp_ak135( 48) =   9.48140000000000
  vp_ak135( 49) =   9.42970000000000
  vp_ak135( 50) =   9.37600000000000
  vp_ak135( 51) =   9.32050000000000
  vp_ak135( 52) =   9.26340000000000
  vp_ak135( 53) =   9.20420000000000
  vp_ak135( 54) =   9.14260000000000
  vp_ak135( 55) =   9.07920000000000
  vp_ak135( 56) =   9.01380000000000
  vp_ak135( 57) =   8.94610000000000
  vp_ak135( 58) =   8.87610000000000
  vp_ak135( 59) =   8.80360000000000
  vp_ak135( 60) =   8.72830000000000
  vp_ak135( 61) =   8.64960000000000
  vp_ak135( 62) =   8.56920000000000
  vp_ak135( 63) =   8.48610000000000
  vp_ak135( 64) =   8.40010000000000
  vp_ak135( 65) =   8.31220000000000
  vp_ak135( 66) =   8.22130000000000
  vp_ak135( 67) =   8.12830000000000
  vp_ak135( 68) =   8.03820000000000
  vp_ak135( 69) =   8.00000000000000
  vp_ak135( 70) =   13.6601000000000
  vp_ak135( 71) =   13.6570000000000
  vp_ak135( 72) =   13.6533000000000
  vp_ak135( 73) =   13.6498000000000
  vp_ak135( 74) =   13.6498000000000
  vp_ak135( 75) =   13.5899000000000
  vp_ak135( 76) =   13.5311000000000
  vp_ak135( 77) =   13.4741000000000
  vp_ak135( 78) =   13.4156000000000
  vp_ak135( 79) =   13.3584000000000
  vp_ak135( 80) =   13.3017000000000
  vp_ak135( 81) =   13.2465000000000
  vp_ak135( 82) =   13.1895000000000
  vp_ak135( 83) =   13.1337000000000
  vp_ak135( 84) =   13.0786000000000
  vp_ak135( 85) =   13.0226000000000
  vp_ak135( 86) =   12.9663000000000
  vp_ak135( 87) =   12.9093000000000
  vp_ak135( 88) =   12.8524000000000
  vp_ak135( 89) =   12.7956000000000
  vp_ak135( 90) =   12.7384000000000
  vp_ak135( 91) =   12.6807000000000
  vp_ak135( 92) =   12.6226000000000
  vp_ak135( 93) =   12.5638000000000
  vp_ak135( 94) =   12.5030000000000
  vp_ak135( 95) =   12.4427000000000
  vp_ak135( 96) =   12.3813000000000
  vp_ak135( 97) =   12.3181000000000
  vp_ak135( 98) =   12.2558000000000
  vp_ak135( 99) =   12.1912000000000
  vp_ak135(100) =   12.1247000000000
  vp_ak135(101) =   12.0571000000000
  vp_ak135(102) =   11.9891000000000
  vp_ak135(103) =   11.9208000000000
  vp_ak135(104) =   11.8491000000000
  vp_ak135(105) =   11.7768000000000
  vp_ak135(106) =   11.7020000000000
  vp_ak135(107) =   11.6265000000000
  vp_ak135(108) =   11.5493000000000
  vp_ak135(109) =   11.4704000000000
  vp_ak135(110) =   11.3897000000000
  vp_ak135(111) =   11.3068000000000
  vp_ak135(112) =   11.2228000000000
  vp_ak135(113) =   11.1355000000000
  vp_ak135(114) =   11.0553000000000
  vp_ak135(115) =   10.9222000000000
  vp_ak135(116) =   10.7909000000000
  vp_ak135(117) =   10.2000000000000
  vp_ak135(118) =   10.0320000000000
  vp_ak135(119) =   9.86400000000000
  vp_ak135(120) =   9.69620000000000
  vp_ak135(121) =   9.52800000000000
  vp_ak135(122) =   9.36010000000000
  vp_ak135(123) =   9.03020000000000
  vp_ak135(124) =   8.84760000000000
  vp_ak135(125) =   8.66500000000000
  vp_ak135(126) =   8.48220000000000
  vp_ak135(127) =   8.30070000000000
  vp_ak135(128) =   8.30070000000000
  vp_ak135(129) =   8.17500000000000
  vp_ak135(130) =   8.05050000000000
  vp_ak135(131) =   8.04500000000000
  vp_ak135(132) =   8.04000000000000
  vp_ak135(133) =   6.50000000000000
  vp_ak135(134) =   6.50000000000000
  vp_ak135(135) =   5.80000000000000
  vp_ak135(136) =   5.80000000000000

  vs_ak135(  1) =   3.66780000000000
  vs_ak135(  2) =   3.66750000000000
  vs_ak135(  3) =   3.66670000000000
  vs_ak135(  4) =   3.66530000000000
  vs_ak135(  5) =   3.66330000000000
  vs_ak135(  6) =   3.66080000000000
  vs_ak135(  7) =   3.65770000000000
  vs_ak135(  8) =   3.65400000000000
  vs_ak135(  9) =   3.64980000000000
  vs_ak135( 10) =   3.64500000000000
  vs_ak135( 11) =   3.63960000000000
  vs_ak135( 12) =   3.63370000000000
  vs_ak135( 13) =   3.62020000000000
  vs_ak135( 14) =   3.61300000000000
  vs_ak135( 15) =   3.60440000000000
  vs_ak135( 16) =   3.59570000000000
  vs_ak135( 17) =   3.58640000000000
  vs_ak135( 18) =   3.57650000000000
  vs_ak135( 19) =   3.56610000000000
  vs_ak135( 20) =   3.55510000000000
  vs_ak135( 21) =   3.54350000000000
  vs_ak135( 22) =   3.53140000000000
  vs_ak135( 23) =   3.51870000000000
  vs_ak135( 24) =   3.50430000000000
  vs_ak135( 25) =  0.000000000000000E+000
  vs_ak135( 26) =  0.000000000000000E+000
  vs_ak135( 27) =  0.000000000000000E+000
  vs_ak135( 28) =  0.000000000000000E+000
  vs_ak135( 29) =  0.000000000000000E+000
  vs_ak135( 30) =  0.000000000000000E+000
  vs_ak135( 31) =  0.000000000000000E+000
  vs_ak135( 32) =  0.000000000000000E+000
  vs_ak135( 33) =  0.000000000000000E+000
  vs_ak135( 34) =  0.000000000000000E+000
  vs_ak135( 35) =  0.000000000000000E+000
  vs_ak135( 36) =  0.000000000000000E+000
  vs_ak135( 37) =  0.000000000000000E+000
  vs_ak135( 38) =  0.000000000000000E+000
  vs_ak135( 39) =  0.000000000000000E+000
  vs_ak135( 40) =  0.000000000000000E+000
  vs_ak135( 41) =  0.000000000000000E+000
  vs_ak135( 42) =  0.000000000000000E+000
  vs_ak135( 43) =  0.000000000000000E+000
  vs_ak135( 44) =  0.000000000000000E+000
  vs_ak135( 45) =  0.000000000000000E+000
  vs_ak135( 46) =  0.000000000000000E+000
  vs_ak135( 47) =  0.000000000000000E+000
  vs_ak135( 48) =  0.000000000000000E+000
  vs_ak135( 49) =  0.000000000000000E+000
  vs_ak135( 50) =  0.000000000000000E+000
  vs_ak135( 51) =  0.000000000000000E+000
  vs_ak135( 52) =  0.000000000000000E+000
  vs_ak135( 53) =  0.000000000000000E+000
  vs_ak135( 54) =  0.000000000000000E+000
  vs_ak135( 55) =  0.000000000000000E+000
  vs_ak135( 56) =  0.000000000000000E+000
  vs_ak135( 57) =  0.000000000000000E+000
  vs_ak135( 58) =  0.000000000000000E+000
  vs_ak135( 59) =  0.000000000000000E+000
  vs_ak135( 60) =  0.000000000000000E+000
  vs_ak135( 61) =  0.000000000000000E+000
  vs_ak135( 62) =  0.000000000000000E+000
  vs_ak135( 63) =  0.000000000000000E+000
  vs_ak135( 64) =  0.000000000000000E+000
  vs_ak135( 65) =  0.000000000000000E+000
  vs_ak135( 66) =  0.000000000000000E+000
  vs_ak135( 67) =  0.000000000000000E+000
  vs_ak135( 68) =  0.000000000000000E+000
  vs_ak135( 69) =  0.000000000000000E+000
  vs_ak135( 70) =   7.28170000000000
  vs_ak135( 71) =   7.27000000000000
  vs_ak135( 72) =   7.25930000000000
  vs_ak135( 73) =   7.24850000000000
  vs_ak135( 74) =   7.24850000000000
  vs_ak135( 75) =   7.22530000000000
  vs_ak135( 76) =   7.20310000000000
  vs_ak135( 77) =   7.18040000000000
  vs_ak135( 78) =   7.15840000000000
  vs_ak135( 79) =   7.13680000000000
  vs_ak135( 80) =   7.11440000000000
  vs_ak135( 81) =   7.09320000000000
  vs_ak135( 82) =   7.07220000000000
  vs_ak135( 83) =   7.05040000000000
  vs_ak135( 84) =   7.02860000000000
  vs_ak135( 85) =   7.00690000000000
  vs_ak135( 86) =   6.98520000000000
  vs_ak135( 87) =   6.96250000000000
  vs_ak135( 88) =   6.94160000000000
  vs_ak135( 89) =   6.91940000000000
  vs_ak135( 90) =   6.89720000000000
  vs_ak135( 91) =   6.87430000000000
  vs_ak135( 92) =   6.85170000000000
  vs_ak135( 93) =   6.82890000000000
  vs_ak135( 94) =   6.80560000000000
  vs_ak135( 95) =   6.78200000000000
  vs_ak135( 96) =   6.75790000000000
  vs_ak135( 97) =   6.73230000000000
  vs_ak135( 98) =   6.70700000000000
  vs_ak135( 99) =   6.68130000000000
  vs_ak135(100) =   6.65540000000000
  vs_ak135(101) =   6.62850000000000
  vs_ak135(102) =   6.60090000000000
  vs_ak135(103) =   6.57280000000000
  vs_ak135(104) =   6.54310000000000
  vs_ak135(105) =   6.51310000000000
  vs_ak135(106) =   6.48220000000000
  vs_ak135(107) =   6.45140000000000
  vs_ak135(108) =   6.41820000000000
  vs_ak135(109) =   6.38600000000000
  vs_ak135(110) =   6.35190000000000
  vs_ak135(111) =   6.31640000000000
  vs_ak135(112) =   6.27990000000000
  vs_ak135(113) =   6.24240000000000
  vs_ak135(114) =   6.21000000000000
  vs_ak135(115) =   6.08980000000000
  vs_ak135(116) =   5.96070000000000
  vs_ak135(117) =   5.61040000000000
  vs_ak135(118) =   5.50470000000000
  vs_ak135(119) =   5.39890000000000
  vs_ak135(120) =   5.29220000000000
  vs_ak135(121) =   5.18640000000000
  vs_ak135(122) =   5.08060000000000
  vs_ak135(123) =   4.87020000000000
  vs_ak135(124) =   4.78320000000000
  vs_ak135(125) =   4.69640000000000
  vs_ak135(126) =   4.60940000000000
  vs_ak135(127) =   4.51840000000000
  vs_ak135(128) =   4.51840000000000
  vs_ak135(129) =   4.50900000000000
  vs_ak135(130) =   4.50000000000000
  vs_ak135(131) =   4.49000000000000
  vs_ak135(132) =   4.48000000000000
  vs_ak135(133) =   3.85000000000000
  vs_ak135(134) =   3.85000000000000
  vs_ak135(135) =   3.46000000000000
  vs_ak135(136) =   3.46000000000000

  Qkappa_ak135(  1) =   601.270000000000
  Qkappa_ak135(  2) =   601.320000000000
  Qkappa_ak135(  3) =   601.460000000000
  Qkappa_ak135(  4) =   601.700000000000
  Qkappa_ak135(  5) =   602.050000000000
  Qkappa_ak135(  6) =   602.490000000000
  Qkappa_ak135(  7) =   603.040000000000
  Qkappa_ak135(  8) =   603.690000000000
  Qkappa_ak135(  9) =   604.440000000000
  Qkappa_ak135( 10) =   605.280000000000
  Qkappa_ak135( 11) =   606.260000000000
  Qkappa_ak135( 12) =   607.310000000000
  Qkappa_ak135( 13) =   609.740000000000
  Qkappa_ak135( 14) =   611.180000000000
  Qkappa_ak135( 15) =   612.620000000000
  Qkappa_ak135( 16) =   614.210000000000
  Qkappa_ak135( 17) =   615.930000000000
  Qkappa_ak135( 18) =   617.780000000000
  Qkappa_ak135( 19) =   619.710000000000
  Qkappa_ak135( 20) =   621.500000000000
  Qkappa_ak135( 21) =   624.080000000000
  Qkappa_ak135( 22) =   626.870000000000
  Qkappa_ak135( 23) =   629.890000000000
  Qkappa_ak135( 24) =   633.260000000000
  Qkappa_ak135( 25) =   57822.0000000000
  Qkappa_ak135( 26) =   57822.0000000000
  Qkappa_ak135( 27) =   57822.0000000000
  Qkappa_ak135( 28) =   57822.0000000000
  Qkappa_ak135( 29) =   57822.0000000000
  Qkappa_ak135( 30) =   57822.0000000000
  Qkappa_ak135( 31) =   57822.0000000000
  Qkappa_ak135( 32) =   57822.0000000000
  Qkappa_ak135( 33) =   57822.0000000000
  Qkappa_ak135( 34) =   57822.0000000000
  Qkappa_ak135( 35) =   57822.0000000000
  Qkappa_ak135( 36) =   57822.0000000000
  Qkappa_ak135( 37) =   57822.0000000000
  Qkappa_ak135( 38) =   57822.0000000000
  Qkappa_ak135( 39) =   57822.0000000000
  Qkappa_ak135( 40) =   57822.0000000000
  Qkappa_ak135( 41) =   57822.0000000000
  Qkappa_ak135( 42) =   57822.0000000000
  Qkappa_ak135( 43) =   57822.0000000000
  Qkappa_ak135( 44) =   57822.0000000000
  Qkappa_ak135( 45) =   57822.0000000000
  Qkappa_ak135( 46) =   57822.0000000000
  Qkappa_ak135( 47) =   57822.0000000000
  Qkappa_ak135( 48) =   57822.0000000000
  Qkappa_ak135( 49) =   57822.0000000000
  Qkappa_ak135( 50) =   57822.0000000000
  Qkappa_ak135( 51) =   57822.0000000000
  Qkappa_ak135( 52) =   57822.0000000000
  Qkappa_ak135( 53) =   57822.0000000000
  Qkappa_ak135( 54) =   57822.0000000000
  Qkappa_ak135( 55) =   57822.0000000000
  Qkappa_ak135( 56) =   57822.0000000000
  Qkappa_ak135( 57) =   57822.0000000000
  Qkappa_ak135( 58) =   57822.0000000000
  Qkappa_ak135( 59) =   57822.0000000000
  Qkappa_ak135( 60) =   57822.0000000000
  Qkappa_ak135( 61) =   57822.0000000000
  Qkappa_ak135( 62) =   57822.0000000000
  Qkappa_ak135( 63) =   57822.0000000000
  Qkappa_ak135( 64) =   57822.0000000000
  Qkappa_ak135( 65) =   57822.0000000000
  Qkappa_ak135( 66) =   57822.0000000000
  Qkappa_ak135( 67) =   57822.0000000000
  Qkappa_ak135( 68) =   57822.0000000000
  Qkappa_ak135( 69) =   57822.0000000000
  Qkappa_ak135( 70) =   723.120000000000
  Qkappa_ak135( 71) =   725.110000000000
  Qkappa_ak135( 72) =   726.870000000000
  Qkappa_ak135( 73) =   722.730000000000
  Qkappa_ak135( 74) =   933.210000000000
  Qkappa_ak135( 75) =   940.880000000000
  Qkappa_ak135( 76) =   952.000000000000
  Qkappa_ak135( 77) =   960.360000000000
  Qkappa_ak135( 78) =   968.460000000000
  Qkappa_ak135( 79) =   976.810000000000
  Qkappa_ak135( 80) =   985.630000000000
  Qkappa_ak135( 81) =   990.770000000000
  Qkappa_ak135( 82) =   999.440000000000
  Qkappa_ak135( 83) =   1008.79000000000
  Qkappa_ak135( 84) =   1018.38000000000
  Qkappa_ak135( 85) =   1032.14000000000
  Qkappa_ak135( 86) =   1042.07000000000
  Qkappa_ak135( 87) =   1048.09000000000
  Qkappa_ak135( 88) =   1058.03000000000
  Qkappa_ak135( 89) =   1064.23000000000
  Qkappa_ak135( 90) =   1070.38000000000
  Qkappa_ak135( 91) =   1085.97000000000
  Qkappa_ak135( 92) =   1097.16000000000
  Qkappa_ak135( 93) =   1108.58000000000
  Qkappa_ak135( 94) =   1120.09000000000
  Qkappa_ak135( 95) =   1127.02000000000
  Qkappa_ak135( 96) =   1134.01000000000
  Qkappa_ak135( 97) =   1141.32000000000
  Qkappa_ak135( 98) =   1148.76000000000
  Qkappa_ak135( 99) =   1156.04000000000
  Qkappa_ak135(100) =   1163.16000000000
  Qkappa_ak135(101) =   1170.53000000000
  Qkappa_ak135(102) =   1178.19000000000
  Qkappa_ak135(103) =   1186.06000000000
  Qkappa_ak135(104) =   1193.99000000000
  Qkappa_ak135(105) =   1202.04000000000
  Qkappa_ak135(106) =   1210.02000000000
  Qkappa_ak135(107) =   1217.91000000000
  Qkappa_ak135(108) =   1226.52000000000
  Qkappa_ak135(109) =   1234.54000000000
  Qkappa_ak135(110) =   1243.02000000000
  Qkappa_ak135(111) =   1251.69000000000
  Qkappa_ak135(112) =   1260.68000000000
  Qkappa_ak135(113) =   1269.44000000000
  Qkappa_ak135(114) =   1277.93000000000
  Qkappa_ak135(115) =   1311.17000000000
  Qkappa_ak135(116) =   1350.54000000000
  Qkappa_ak135(117) =   428.690000000000
  Qkappa_ak135(118) =   425.510000000000
  Qkappa_ak135(119) =   422.550000000000
  Qkappa_ak135(120) =   419.940000000000
  Qkappa_ak135(121) =   417.320000000000
  Qkappa_ak135(122) =   413.660000000000
  Qkappa_ak135(123) =   377.930000000000
  Qkappa_ak135(124) =   366.340000000000
  Qkappa_ak135(125) =   355.850000000000
  Qkappa_ak135(126) =   346.370000000000
  Qkappa_ak135(127) =   338.470000000000
  Qkappa_ak135(128) =   200.970000000000
  Qkappa_ak135(129) =   188.720000000000
  Qkappa_ak135(130) =   182.570000000000
  Qkappa_ak135(131) =   182.030000000000
  Qkappa_ak135(132) =   182.030000000000
  Qkappa_ak135(133) =   972.770000000000
  Qkappa_ak135(134) =   972.770000000000
  Qkappa_ak135(135) =   1368.02000000000
  Qkappa_ak135(136) =   1368.02000000000

  Qmu_ak135(  1) =   85.0300000000000
  Qmu_ak135(  2) =   85.0300000000000
  Qmu_ak135(  3) =   85.0300000000000
  Qmu_ak135(  4) =   85.0300000000000
  Qmu_ak135(  5) =   85.0300000000000
  Qmu_ak135(  6) =   85.0300000000000
  Qmu_ak135(  7) =   85.0300000000000
  Qmu_ak135(  8) =   85.0300000000000
  Qmu_ak135(  9) =   85.0300000000000
  Qmu_ak135( 10) =   85.0300000000000
  Qmu_ak135( 11) =   85.0300000000000
  Qmu_ak135( 12) =   85.0300000000000
  Qmu_ak135( 13) =   85.0300000000000
  Qmu_ak135( 14) =   85.0300000000000
  Qmu_ak135( 15) =   85.0300000000000
  Qmu_ak135( 16) =   85.0300000000000
  Qmu_ak135( 17) =   85.0300000000000
  Qmu_ak135( 18) =   85.0300000000000
  Qmu_ak135( 19) =   85.0300000000000
  Qmu_ak135( 20) =   85.0300000000000
  Qmu_ak135( 21) =   85.0300000000000
  Qmu_ak135( 22) =   85.0300000000000
  Qmu_ak135( 23) =   85.0300000000000
  Qmu_ak135( 24) =   85.0300000000000
  Qmu_ak135( 25) =  0.000000000000000E+000
  Qmu_ak135( 26) =  0.000000000000000E+000
  Qmu_ak135( 27) =  0.000000000000000E+000
  Qmu_ak135( 28) =  0.000000000000000E+000
  Qmu_ak135( 29) =  0.000000000000000E+000
  Qmu_ak135( 30) =  0.000000000000000E+000
  Qmu_ak135( 31) =  0.000000000000000E+000
  Qmu_ak135( 32) =  0.000000000000000E+000
  Qmu_ak135( 33) =  0.000000000000000E+000
  Qmu_ak135( 34) =  0.000000000000000E+000
  Qmu_ak135( 35) =  0.000000000000000E+000
  Qmu_ak135( 36) =  0.000000000000000E+000
  Qmu_ak135( 37) =  0.000000000000000E+000
  Qmu_ak135( 38) =  0.000000000000000E+000
  Qmu_ak135( 39) =  0.000000000000000E+000
  Qmu_ak135( 40) =  0.000000000000000E+000
  Qmu_ak135( 41) =  0.000000000000000E+000
  Qmu_ak135( 42) =  0.000000000000000E+000
  Qmu_ak135( 43) =  0.000000000000000E+000
  Qmu_ak135( 44) =  0.000000000000000E+000
  Qmu_ak135( 45) =  0.000000000000000E+000
  Qmu_ak135( 46) =  0.000000000000000E+000
  Qmu_ak135( 47) =  0.000000000000000E+000
  Qmu_ak135( 48) =  0.000000000000000E+000
  Qmu_ak135( 49) =  0.000000000000000E+000
  Qmu_ak135( 50) =  0.000000000000000E+000
  Qmu_ak135( 51) =  0.000000000000000E+000
  Qmu_ak135( 52) =  0.000000000000000E+000
  Qmu_ak135( 53) =  0.000000000000000E+000
  Qmu_ak135( 54) =  0.000000000000000E+000
  Qmu_ak135( 55) =  0.000000000000000E+000
  Qmu_ak135( 56) =  0.000000000000000E+000
  Qmu_ak135( 57) =  0.000000000000000E+000
  Qmu_ak135( 58) =  0.000000000000000E+000
  Qmu_ak135( 59) =  0.000000000000000E+000
  Qmu_ak135( 60) =  0.000000000000000E+000
  Qmu_ak135( 61) =  0.000000000000000E+000
  Qmu_ak135( 62) =  0.000000000000000E+000
  Qmu_ak135( 63) =  0.000000000000000E+000
  Qmu_ak135( 64) =  0.000000000000000E+000
  Qmu_ak135( 65) =  0.000000000000000E+000
  Qmu_ak135( 66) =  0.000000000000000E+000
  Qmu_ak135( 67) =  0.000000000000000E+000
  Qmu_ak135( 68) =  0.000000000000000E+000
  Qmu_ak135( 69) =  0.000000000000000E+000
  Qmu_ak135( 70) =   273.970000000000
  Qmu_ak135( 71) =   273.970000000000
  Qmu_ak135( 72) =   273.970000000000
  Qmu_ak135( 73) =   271.740000000000
  Qmu_ak135( 74) =   350.880000000000
  Qmu_ak135( 75) =   354.610000000000
  Qmu_ak135( 76) =   359.710000000000
  Qmu_ak135( 77) =   363.640000000000
  Qmu_ak135( 78) =   367.650000000000
  Qmu_ak135( 79) =   371.750000000000
  Qmu_ak135( 80) =   375.940000000000
  Qmu_ak135( 81) =   378.790000000000
  Qmu_ak135( 82) =   383.140000000000
  Qmu_ak135( 83) =   387.600000000000
  Qmu_ak135( 84) =   392.160000000000
  Qmu_ak135( 85) =   398.410000000000
  Qmu_ak135( 86) =   403.230000000000
  Qmu_ak135( 87) =   406.500000000000
  Qmu_ak135( 88) =   411.520000000000
  Qmu_ak135( 89) =   414.940000000000
  Qmu_ak135( 90) =   418.410000000000
  Qmu_ak135( 91) =   425.530000000000
  Qmu_ak135( 92) =   431.030000000000
  Qmu_ak135( 93) =   436.680000000000
  Qmu_ak135( 94) =   442.480000000000
  Qmu_ak135( 95) =   446.430000000000
  Qmu_ak135( 96) =   450.450000000000
  Qmu_ak135( 97) =   454.550000000000
  Qmu_ak135( 98) =   458.720000000000
  Qmu_ak135( 99) =   462.960000000000
  Qmu_ak135(100) =   467.290000000000
  Qmu_ak135(101) =   471.700000000000
  Qmu_ak135(102) =   476.190000000000
  Qmu_ak135(103) =   480.770000000000
  Qmu_ak135(104) =   485.440000000000
  Qmu_ak135(105) =   490.200000000000
  Qmu_ak135(106) =   495.050000000000
  Qmu_ak135(107) =   500.000000000000
  Qmu_ak135(108) =   505.050000000000
  Qmu_ak135(109) =   510.200000000000
  Qmu_ak135(110) =   515.460000000000
  Qmu_ak135(111) =   520.830000000000
  Qmu_ak135(112) =   526.320000000000
  Qmu_ak135(113) =   531.910000000000
  Qmu_ak135(114) =   537.630000000000
  Qmu_ak135(115) =   543.480000000000
  Qmu_ak135(116) =   549.450000000000
  Qmu_ak135(117) =   172.930000000000
  Qmu_ak135(118) =   170.820000000000
  Qmu_ak135(119) =   168.780000000000
  Qmu_ak135(120) =   166.800000000000
  Qmu_ak135(121) =   164.870000000000
  Qmu_ak135(122) =   162.500000000000
  Qmu_ak135(123) =   146.570000000000
  Qmu_ak135(124) =   142.760000000000
  Qmu_ak135(125) =   139.380000000000
  Qmu_ak135(126) =   136.380000000000
  Qmu_ak135(127) =   133.720000000000
  Qmu_ak135(128) =   79.4000000000000
  Qmu_ak135(129) =   76.5500000000000
  Qmu_ak135(130) =   76.0600000000000
  Qmu_ak135(131) =   75.6000000000000
  Qmu_ak135(132) =   75.6000000000000
  Qmu_ak135(133) =   403.930000000000
  Qmu_ak135(134) =   403.930000000000
  Qmu_ak135(135) =   599.990000000000
  Qmu_ak135(136) =   599.990000000000

! strip the crust and replace it with mantle if needed
! if (SUPPRESS_CRUSTAL_MESH .or. USE_EXTERNAL_CRUSTAL_MODEL) then
!   vp_ak135(133:136) = vp_ak135(132)
!   vs_ak135(133:136) = vs_ak135(132)
!   density_ak135(133:136) = density_ak135(132)
!   Qkappa_ak135(133:136) = Qkappa_ak135(132)
!   Qmu_ak135(133:136) = Qmu_ak135(132)
! endif

! loop on all the elements of the mesh, and inside each element loop on all the GLL points
  do ispec = 1,nspec

  if (material_element(ispec) /= IREGION_MANTLE_CRUST_ABOVE_d670 .and. &
     material_element(ispec) /= IREGION_MANTLE_BELOW_d670 .and. &
     material_element(ispec) /= IREGION_OUTER_CORE .and. &
     material_element(ispec) /= IREGION_INNER_CORE) stop 'wrong flag number in external model'

    do j = 1,NGLLZ
      do i = 1,NGLLX

   iglob = ibool(i,j,ispec)

   x = coord(1,iglob)
   z = coord(2,iglob)

! compute the radius
  r = sqrt(x**2 + z**2)

  ii = 1
  do while(r >= radius_ak135(ii) .and. ii /= NR_AK135F_NO_MUD)
    ii = ii + 1
  enddo

! make sure we stay in the right region and never take a point above
! and a point below the ICB or the CMB and interpolate between them,
! which would lead to a wrong value (keeping in mind that we interpolate
! between points i-1 and i below)
  if (material_element(ispec) == IREGION_INNER_CORE .and. ii > 24) ii = 24

  if (material_element(ispec) == IREGION_OUTER_CORE .and. ii < 26) ii = 26
  if (material_element(ispec) == IREGION_OUTER_CORE .and. ii > 69) ii = 69

  if ((material_element(ispec) == IREGION_MANTLE_CRUST_ABOVE_d670 .or. &
      material_element(ispec) == IREGION_MANTLE_BELOW_d670) .and. ii < 71) ii = 71

  if (ii == 1) then
    rho(i,j,ispec) = density_ak135(1)
    vp(i,j,ispec) = vp_ak135(1)
    vs(i,j,ispec) = vs_ak135(1)
    Qmu_attenuation(i,j,ispec) = Qmu_ak135(1)
    Qkappa_attenuation(i,j,ispec) = Qkappa_ak135(1)
  else

! interpolate from radius_ak135(ii-1) to r using the values at ii-1 and ii
    frac = (r-radius_ak135(ii-1))/(radius_ak135(ii)-radius_ak135(ii-1))

    rho(i,j,ispec) = density_ak135(ii-1) + frac*(density_ak135(ii)-density_ak135(ii-1))
    vp(i,j,ispec) = vp_ak135(ii-1) + frac*(vp_ak135(ii)-vp_ak135(ii-1))
    vs(i,j,ispec) = vs_ak135(ii-1) + frac*(vs_ak135(ii)-vs_ak135(ii-1))
    Qmu_attenuation(i,j,ispec) = Qmu_ak135(ii-1) + frac*(Qmu_ak135(ii)-Qmu_ak135(ii-1))
    Qkappa_attenuation(i,j,ispec) = Qkappa_ak135(ii-1) + frac*(Qkappa_ak135(ii)-Qkappa_ak135(ii-1))

  endif

! make sure Vs is zero in the outer core even if roundoff errors on depth
! also set fictitious attenuation to a very high value (attenuation is not used in the fluid)
  if (material_element(ispec) == IREGION_OUTER_CORE) then
    vs(i,j,ispec) = 0.d0
    Qkappa_attenuation(i,j,ispec) = 9999.d0
    Qmu_attenuation(i,j,ispec) = 9999.d0
  endif

      enddo
    enddo
  enddo

! convert to m/s
  vp(:,:,:)=vp(:,:,:)*1000.0d0
  vs(:,:,:)=vs(:,:,:)*1000.0d0

! no anisotropy
  c11(:,:,:) = 0.d0
  c13(:,:,:) = 0.d0
  c15(:,:,:) = 0.d0
  c33(:,:,:) = 0.d0
  c35(:,:,:) = 0.d0
  c55(:,:,:) = 0.d0
  c12(:,:,:) = 0.d0
  c23(:,:,:) = 0.d0
  c25(:,:,:) = 0.d0

  end subroutine define_external_model


!========================================================================
! another example below, to read data from a 1D atmosphere model
! including gravity
!========================================================================

  subroutine define_external_model_atmos_tabular_gravitoacoustic(coord,material_element,ibool, &
                                                                 rho,vp,vs,QKappa_attenuation,Qmu_attenuation,gravity,Nsq, &
                                                                 c11,c13,c15,c33,c35,c55,c12,c23,c25,nspec,nglob, &
                                                                 windx, windz, pext, gammaext, etaext, muext, &
                                                                 kappa_DG, Bxext, Bzext, N0ext)

  use constants,only: CUSTOM_REAL,NGLLX,NGLLZ,NDIM,IMAIN
  use specfem_par,only: nglob_DG, ibool_DG, Htabext_DG, cp, c_V, coord_interface, &
        poroelastcoef,density,kmato, ispec_is_acoustic_DG, ispec_is_elastic, IONOSPHERIC_COUPLING

  implicit none

!--------------------------------------------------------------------------------------------------
!
!          Model is either test one or being based on the Nrl_MSISE2000 atmosphere model
!          ----------------------------------------------------------

! 1D model:
!
! The 1D model will be used to create a tabular model of atmosphere

!--------------------------------------------------------------------------------------------------

  integer, intent(in) :: nspec,nglob

  double precision, dimension(NDIM,nglob), intent(in) :: coord

  integer, dimension(nspec), intent(in) :: material_element

  integer, dimension(NGLLX,NGLLZ,nspec), intent(in) :: ibool

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ,nspec), intent(out) :: rho,vp,vs,QKappa_attenuation,Qmu_attenuation,gravity,Nsq, &
                                              c11,c15,c13,c33,c35,c55,c12,c23,c25, &
                                              ! MODIF DG
                                              windx, windz, pext, etaext, muext, kappa_DG, N0ext
  real(kind=CUSTOM_REAL), dimension(nglob_DG), intent(out) :: gammaext, Bxext, Bzext

! number of layers in the model
  integer, parameter :: NR_LAYER = 99000!4900!

  double precision, dimension(NR_LAYER) :: z_atmos
  double precision, dimension(NR_LAYER) :: density_atmos
  double precision, dimension(NR_LAYER) :: vp_atmos
  double precision, dimension(NR_LAYER) :: gravity_atmos
  double precision, dimension(NR_LAYER) :: Nsq_atmos
  ! MODIF DG
  double precision, dimension(NR_LAYER) :: wx_atmos
  double precision, dimension(NR_LAYER) :: wz_atmos
  double precision, dimension(NR_LAYER) :: eta_atmos
  double precision, dimension(NR_LAYER) :: mu_atmos
  double precision, dimension(NR_LAYER) :: p_atmos
  double precision, dimension(NR_LAYER) :: Htab
  double precision, dimension(NR_LAYER) :: kappa_atmos
  double precision, dimension(NR_LAYER) :: gamma_atmos
  ! Electron density
  double precision, dimension(NR_LAYER) :: Bx_atmos, Bz_atmos, N0_atmos
  
! region flag to assign the Atmosphere model
  integer, parameter :: IREGION_AIR = 1

  integer :: i,j,ispec,iglob,ii

  double precision :: x,z,frac,tmp1, gamma_temp_prev, gamma_temp, pii, piim1, piim2,piip1, max_z

! read all the values in the 1D model once and for all
  !open(10,file='./msise_enhanced_2.txt', &
  !open(10,file='./msise_enhanced_nonlinear_2.txt', &
  !open(10,file='./msise_winds_Sumatra_cross2.txt', &
  !open(10,file='./msise_hydrostatic.txt', &
  open(10,file='./msise_Sumatra_model_DG.dat', &
  !open(10,file='./msise_enhanced_lat16_modrho3.txt', &
  form='formatted')
  
  !open(11,file='./extract_profile.txt', &
  !form='formatted')
  !if(z < 100000 .AND. z > 0) &
  !  WRITE(11,*) z, rho(i,j,ispec), pext(i,j,ispec), p_atmos(ii_aux), (z-z_atmos(ii_aux)), density_atmos(ii_aux), &
  !      z_save, ii, ii_aux, ii_prev

  if(.false.) WRITE(*,*) material_element

  do i = 1,NR_LAYER

  !read(10,*) z_atmos(i),density_atmos(i),tmp1,vp_atmos(i),tmp2,gravity_atmos(i),Nsq_atmos(i)
  !!!!!!!!!!!!!!!!
  ! msise sumatra
  ! tmp1 = Temperature
  ! tmp2 = Pressure
  ! tmp3 = Scale height
  ! tmp4 = kappa
  gamma_atmos(i) = 0.
  if(IONOSPHERIC_COUPLING) then
  read(10,*) z_atmos(i),density_atmos(i),tmp1,vp_atmos(i),p_atmos(i),Htab(i),gravity_atmos(i),Nsq_atmos(i),kappa_atmos(i),&
                mu_atmos(i),eta_atmos(i), wx_atmos(i), wz_atmos(i), cp, c_V, gamma_atmos(i), &
                Bx_atmos(i), Bz_atmos(i), N0_atmos(i)
  else
  read(10,*) z_atmos(i),density_atmos(i),tmp1,vp_atmos(i),p_atmos(i),Htab(i),gravity_atmos(i),Nsq_atmos(i),kappa_atmos(i),&
                mu_atmos(i),eta_atmos(i), wx_atmos(i), wz_atmos(i), cp, c_V, gamma_atmos(i)
  endif
  enddo
  allocate(Htabext_DG(nglob_DG))

  close(10)

  windx = 0
  windz = 0

! loop on all the elements of the mesh, and inside each element loop on all the GLL points
  do ispec = 1,nspec

    ! Recover current element's max height along z to discriminate between regions
    max_z = 0.
    do j = 1,NGLLZ
        iglob = ibool(1,j,ispec)
        z = coord(2,iglob) - coord_interface
        if(j == 1 .OR. z > max_z) max_z = z
    enddo

    do j = 1,NGLLZ
      do i = 1,NGLLX
      
       iglob = ibool(i,j,ispec)

       x = coord(1,iglob)
       z = coord(2,iglob) - coord_interface

       ii = 1
       do while(z >= z_atmos(ii) .and. ii /= NR_LAYER)
         ii = ii + 1
       enddo
      
   !if (material_element(ispec) /= IREGION_AIR ) then
   !if(z <= 0 .AND. max_z <= 0.) then
   if(ispec_is_elastic(ispec)) then
   
          rho(i,j,ispec) = 2500.d0
          vp(i,j,ispec) = 3600.d0
          vs(i,j,ispec) = vp(i,j,ispec) / 2.d0
          QKappa_attenuation(i,j,ispec) = 120.
          Qmu_attenuation(i,j,ispec) = 120.
          c11(i,j,ispec) = 0.d0   ! this means no anisotropy
          c13(i,j,ispec) = 0.d0
          c15(i,j,ispec) = 0.d0
          c33(i,j,ispec) = 0.d0
          c35(i,j,ispec) = 0.d0
          c55(i,j,ispec) = 0.d0
          c12(i,j,ispec) = 0.d0
          c23(i,j,ispec) = 0.d0
          c25(i,j,ispec) = 0.d0
        !! AB AB Do not forget these 3 lines otherwise PML may not work !!
        density(1,kmato(ispec)) = rho(i,j,ispec)
        poroelastcoef(3,1,kmato(ispec)) = rho(i,j,ispec) * vp(i,j,ispec) * vp(i,j,ispec)
        poroelastcoef(2,1,kmato(ispec)) =  rho(i,j,ispec) * vs(i,j,ispec) * vs(i,j,ispec)
  
  else
      
  ! To be changed if not DG simulation
  ispec_is_acoustic_DG(ispec) = .true.

  if (ii == 1) then
    rho(i,j,ispec) = density_atmos(1)
    vp(i,j,ispec) = vp_atmos(1)
    gravity(i,j,ispec) = gravity_atmos(1)
    vs(i,j,ispec) = 0.d0
    Qmu_attenuation(i,j,ispec) = 9999.d0
    Qkappa_attenuation(i,j,ispec) = 9999.d0
    
    ! MODIF DG
    windx(i,j,ispec) = wx_atmos(1)
    windz(i,j,ispec) = wz_atmos(1)
    pext(i,j,ispec) = p_atmos(1)
    gamma_temp = 1.4!gamma_atmos(1)!(vp_atmos(1)**2)*density_atmos(1)/p_atmos(1)!(vp_atmos(1)**2)/(gravity_atmos(1)*Htab(1))
    gammaext(ibool_DG(i,j,ispec)) = gamma_temp
    Htabext_DG(ibool_DG(i,j,ispec)) = Htab(1)
    etaext(i,j,ispec) = eta_atmos(1)
    muext(i,j,ispec)  = mu_atmos(1)
    kappa_DG(i,j,ispec)  = kappa_atmos(1)
    
    if(IONOSPHERIC_COUPLING) then
    Bxext(ibool_DG(i,j,ispec)) = Bx_atmos(1)
    Bzext(ibool_DG(i,j,ispec))  = Bz_atmos(1)
    N0ext(i,j,ispec)  = N0_atmos(1)
    endif
    
  elseif (ii == 2) then

! interpolate from radius_ak135(ii-1) to r using the values at ii-1 and ii
    frac = (z-z_atmos(ii-1))/(z_atmos(ii)-z_atmos(ii-1))

    rho(i,j,ispec) = density_atmos(ii-1) + frac*(density_atmos(ii)-density_atmos(ii-1))
    !rho(i,j,ispec) = exp(log(density_atmos(ii-1)) + frac*(log(density_atmos(ii))-log(density_atmos(ii-1))))
    
    vp(i,j,ispec) = vp_atmos(ii-1) + frac*(vp_atmos(ii)-vp_atmos(ii-1))
    gravity(i,j,ispec) = gravity_atmos(ii-1) + frac*(gravity_atmos(ii)-gravity_atmos(ii-1))
    Nsq(i,j,ispec) = Nsq_atmos(ii-1) + frac*(Nsq_atmos(ii)-Nsq_atmos(ii-1))
    if (Nsq(i,j,ispec) <= 0.0) then
       write(*,*) 'STOP Negative Nsquare !!! :', &
       i,j,coord(1,iglob),coord(2,iglob),Nsq(i,j,ispec),gravity(i,j,ispec),vp(i,j,ispec)
       stop
    endif
    vs(i,j,ispec) = 0.d0
    Qmu_attenuation(i,j,ispec) = 9999.d0
    Qkappa_attenuation(i,j,ispec) = 9999.d0
    
    ! MODIF DG
    windx(i,j,ispec) = wx_atmos(ii-1) + frac*(wx_atmos(ii)-wx_atmos(ii-1))
    windz(i,j,ispec) = wz_atmos(ii-1) + frac*(wz_atmos(ii)-wz_atmos(ii-1))

    pext(i,j,ispec) = p_atmos(ii-1) + frac*(p_atmos(ii)-p_atmos(ii-1))
    !pext(i,j,ispec) = exp(log(p_atmos(ii-1)) + frac*(log(p_atmos(ii))-log(p_atmos(ii-1))))
    
    !if(z > z_save) ii_aux = ii_prev
    !pext(i,j,ispec) = p_atmos(ii_aux-1) - (z-z_atmos(ii_aux-1))*density_atmos(ii_aux-1)*9.81
    !p_atmos(ii)     = pext(i,j,ispec)
    !ii_prev         = ii
    !z_save = z
    !if(z < 100000 .AND. z > 0) &
    !WRITE(*,*) "---------->",z, rho(i,j,ispec), pext(i,j,ispec), p_atmos(ii-1), (z-z_atmos(ii-1)), density_atmos(ii-1), &
    !    z_save, ii, ii_aux, ii_prev
    !pext(i,j,ispec) = p_atmos(ii-1) + frac*(p_atmos(ii)-p_atmos(ii-1))
    gamma_temp_prev = gamma_temp
    gamma_temp = (vp_atmos(ii)**2)*density_atmos(ii)/p_atmos(ii)!(vp_atmos(ii)**2)/(gravity_atmos(ii)*Htab(ii))
    gammaext(ibool_DG(i,j,ispec)) = 1.4!        gamma_atmos(ii-1) + frac*(gamma_atmos(ii)-gamma_atmos(ii-1))! (vp(i,j,ispec)**2)*rho(i,j,ispec)/pext(i,j,ispec)!gamma_temp_prev + frac*(gamma_temp-gamma_temp_prev)

    etaext(i,j,ispec) = eta_atmos(ii-1) + frac*(eta_atmos(ii)-eta_atmos(ii-1))
    muext(i,j,ispec)  = mu_atmos(ii-1) + frac*(mu_atmos(ii)-mu_atmos(ii-1))

    Htabext_DG(ibool_DG(i,j,ispec)) = Htab(ii-1) + frac*(Htab(ii)-Htab(ii-1))

    kappa_DG(i,j,ispec)  = kappa_atmos(ii-1) + frac*(kappa_atmos(ii)-kappa_atmos(ii-1))
    
    if(IONOSPHERIC_COUPLING) then
    Bxext(ibool_DG(i,j,ispec)) = Bx_atmos(ii-1) + frac*(Bx_atmos(ii)-Bx_atmos(ii-1))
    Bzext(ibool_DG(i,j,ispec)) = Bz_atmos(ii-1) + frac*(Bz_atmos(ii)-Bz_atmos(ii-1))
    N0ext(i,j,ispec) = N0_atmos(ii-1) + frac*(N0_atmos(ii)-N0_atmos(ii-1))
    endif

  else
  
  if(.false.) then
  pii   = (z-z_atmos(ii+1))*(z-z_atmos(ii-1))*(z-z_atmos(ii-2))&
        /((z_atmos(ii)-z_atmos(ii-1))*(z_atmos(ii)-z_atmos(ii-2))*(z_atmos(ii)-z_atmos(ii+1)))
  piim1 = (z-z_atmos(ii+1))*(z-z_atmos(ii))*(z-z_atmos(ii-2))&
        /((z_atmos(ii-1)-z_atmos(ii))*(z_atmos(ii-1)-z_atmos(ii-2))*(z_atmos(ii-1)-z_atmos(ii+1)))
  piim2 = (z-z_atmos(ii+1))*(z-z_atmos(ii-1))*(z-z_atmos(ii))&
        /((z_atmos(ii-2)-z_atmos(ii))*(z_atmos(ii-2)-z_atmos(ii-1))*(z_atmos(ii-2)-z_atmos(ii+1)))
  piip1 = (z-z_atmos(ii-2))*(z-z_atmos(ii-1))*(z-z_atmos(ii))&
        /((z_atmos(ii+1)-z_atmos(ii))*(z_atmos(ii+1)-z_atmos(ii-1))*(z_atmos(ii+1)-z_atmos(ii-2)))
  endif
        
  pii   = (z-z_atmos(ii-1))*(z-z_atmos(ii-2))/((z_atmos(ii)-z_atmos(ii-1))*(z_atmos(ii)-z_atmos(ii-2)))
  piim1 = (z-z_atmos(ii))*(z-z_atmos(ii-2))/((z_atmos(ii-1)-z_atmos(ii))*(z_atmos(ii-1)-z_atmos(ii-2)))
  piim2 = (z-z_atmos(ii-1))*(z-z_atmos(ii))/((z_atmos(ii-2)-z_atmos(ii))*(z_atmos(ii-2)-z_atmos(ii-1)))
  
  rho(i,j,ispec) = density_atmos(ii)*pii + density_atmos(ii-1)*piim1 + density_atmos(ii-2)*piim2 !+ density_atmos(ii+1)*piip1
  !rho(i,j,ispec) = log(density_atmos(ii))*pii + log(density_atmos(ii-1))*piim1 + log(density_atmos(ii-2))*piim2
   ! rho(i,j,ispec) = exp(rho(i,j,ispec))
    
    vp(i,j,ispec)      = vp_atmos(ii)*pii + vp_atmos(ii-1)*piim1 + vp_atmos(ii-2)*piim2 !+ vp_atmos(ii+1)*piip1
    gravity(i,j,ispec) = gravity_atmos(ii)*pii + gravity_atmos(ii-1)*piim1 + gravity_atmos(ii-2)*piim2 !+ gravity_atmos(ii+1)*piip1
    Nsq(i,j,ispec)     = Nsq_atmos(ii)*pii + Nsq_atmos(ii-1)*piim1 + Nsq_atmos(ii-2)*piim2 !+ Nsq_atmos(ii+1)*piip1
    vs(i,j,ispec) = 0.d0
    Qmu_attenuation(i,j,ispec) = 9999.d0
    Qkappa_attenuation(i,j,ispec) = 9999.d0
    
    ! MODIF DG
    windx(i,j,ispec) = wx_atmos(ii)*pii + wx_atmos(ii-1)*piim1 + wx_atmos(ii-2)*piim2 !+ wx_atmos(ii+1)*piip1
    windz(i,j,ispec) = wz_atmos(ii)*pii + wz_atmos(ii-1)*piim1 + wz_atmos(ii-2)*piim2 !+ wz_atmos(ii+1)*piip1

    pext(i,j,ispec) = p_atmos(ii)*pii + p_atmos(ii-1)*piim1 + p_atmos(ii-2)*piim2! + p_atmos(ii+1)*piip1
    !pext(i,j,ispec) = log(p_atmos(ii))*pii + log(p_atmos(ii-1))*piim1 + log(p_atmos(ii-2))*piim2
    !pext(i,j,ispec) = exp(pext(i,j,ispec))
    
    gammaext(ibool_DG(i,j,ispec)) = 1.4

    etaext(i,j,ispec) = eta_atmos(ii)*pii + eta_atmos(ii-1)*piim1 + eta_atmos(ii-2)*piim2! + eta_atmos(ii+1)*piip1
    muext(i,j,ispec)  = mu_atmos(ii)*pii + mu_atmos(ii-1)*piim1 + mu_atmos(ii-2)*piim2! + mu_atmos(ii+1)*piip1

    Htabext_DG(ibool_DG(i,j,ispec)) = Htab(ii)*pii + Htab(ii-1)*piim1 + Htab(ii-2)*piim2 !+ Htab(ii+1)*piip1

    kappa_DG(i,j,ispec)  = kappa_atmos(ii)*pii + kappa_atmos(ii-1)*piim1 + kappa_atmos(ii-2)*piim2! + kappa_atmos(ii+1)*piip1
    
    if(IONOSPHERIC_COUPLING) then
    Bxext(ibool_DG(i,j,ispec)) = Bx_atmos(ii)*pii + Bx_atmos(ii-1)*piim1 + Bx_atmos(ii-2)*piim2
    Bzext(ibool_DG(i,j,ispec)) = Bz_atmos(ii)*pii + Bz_atmos(ii-1)*piim1 + Bz_atmos(ii-2)*piim2
    N0ext(i,j,ispec) = N0_atmos(ii)*pii + N0_atmos(ii-1)*piim1 + N0_atmos(ii-2)*piim2
    endif

  endif ! IREGION_AIR
  
    endif
  
      enddo
    enddo
    
  enddo

  !WRITE(*,*) "MAXVAL/MINVAL", maxval(rho), minval(rho)
  !stop

! remove gravity for acoustic-only simulations
  !gravity(:,:,:) = 0.d0

  close(11)

! no anisotropy
  c11(:,:,:) = 0.d0
  c13(:,:,:) = 0.d0
  c15(:,:,:) = 0.d0
  c33(:,:,:) = 0.d0
  c35(:,:,:) = 0.d0
  c55(:,:,:) = 0.d0
  c12(:,:,:) = 0.d0
  c23(:,:,:) = 0.d0
  c25(:,:,:) = 0.d0

  end subroutine define_external_model_atmos_tabular_gravitoacoustic
  
! ------------------------------------------------------------ !
! external_model_DG_testfile                                   !
! ------------------------------------------------------------ !
! Finds the total number of lines contained in the 'atmospheric_model.dat' file.
subroutine external_model_DG_testfile(FILENAME)
  implicit none
  ! Input/output.
  character(len=100), intent(in) :: FILENAME
  ! Local variables.
  integer io
  OPEN(100, file=FILENAME)
  READ(100,*,iostat=io)
  if(io/=0) then
    write(*,*) "********************************"
    write(*,*) "*            ERROR             *"
    write(*,*) "********************************"
    write(*,*) "* Cannot read external         *"
    write(*,*) "* model file                   *"
    write(*,*) "* ", trim(FILENAME)
    write(*,*) "********************************"
    stop
  endif
  CLOSE(100)
end subroutine external_model_DG_testfile
! ------------------------------------------------------------ !
! external_model_DG_only_find_nblines                          !
! ------------------------------------------------------------ !
! Finds the total number of lines contained in the 'atmospheric_model.dat' file.
subroutine external_model_DG_only_find_nblines(FILENAME, nlines_header, nlines_model)
  !use specfem_par, only: EXTERNAL_DG_ONLY_MODEL_FILENAME
  implicit none
  ! Input/output.
  integer, intent(out) :: nlines_header, nlines_model
  character(len=100), intent(in) :: FILENAME
  ! Local variables.
  integer nlines_file, io
  nlines_header = 3 ! Number of header lines in model file. ! TODO: Find a way of discriminating this automatically.
  nlines_file = 0
  OPEN(100, file=FILENAME)
  DO
    READ(100,*,iostat=io)
    IF (io/=0) EXIT
    nlines_file = nlines_file+1
  ENDDO 
  CLOSE(100)
  nlines_model=nlines_file-nlines_header
end subroutine external_model_DG_only_find_nblines
! ------------------------------------------------------------ !
! external_model_DG_only_find_nbcols                           !
! ------------------------------------------------------------ !
! Finds the number of columns contained in the 'atmospheric_model.dat' file.
subroutine external_model_DG_only_find_nbcols(FILENAME, nlines_header, ncols)
  !use specfem_par, only: EXTERNAL_DG_ONLY_MODEL_FILENAME
  implicit none
  ! Input/output.
  character(len=100), intent(in) :: FILENAME
  integer, intent(in) :: nlines_header
  integer, intent(out) :: ncols
  ! Local variables.
  integer, parameter :: MAX_NUM_OF_COLS=30
  integer, parameter :: MAX_LINE_LENGTH=1000
  character(len=MAX_LINE_LENGTH) line
  integer i, io, reading_unit
  double precision, dimension(MAX_NUM_OF_COLS) :: kek
  reading_unit=100
  OPEN(reading_unit, file=trim(FILENAME))
  ! Skip header.
  DO i=1,nlines_header
    READ(reading_unit,*,iostat=io)
    IF (io/=0) then
      write(*,*) "Error reading line in file '",trim(FILENAME),"' (iostat=",io,")."
      stop "FATAL ERROR."
    endif
  enddo
  ! Get first line of model.
  DO
    READ(reading_unit,'(A)',iostat=io) line
    IF (io/=0) then
      write(*,*) "Error reading line in file '",trim(FILENAME),"' (iostat=",io,")."
      stop "FATAL ERROR."
    endif
    exit
  ENDDO
  CLOSE(reading_unit)
  !write(*,*) line
  do i=1,MAX_NUM_OF_COLS!,1,-1
    READ(line,*,iostat=io) kek(1:i)
    !write(*,*) 'io', io
    if(io==-1) exit
  enddo
  !write(*,*) 'i',(i-1)
  !stop 'kek'
  ncols=i-1
end subroutine external_model_DG_only_find_nbcols

! ------------------------------------------------------------ !
! define_external_model_DG_only                                !
! ------------------------------------------------------------ !
! Defines model from a vertical model file, but only for DG elements.

subroutine define_external_model_DG_only(nlines_header, nlines_model)

  use constants,only: CUSTOM_REAL,NGLLX,NGLLZ,NDIM,IMAIN,FOUR_THIRDS,TINYVAL,PI
  use specfem_par,only: ibool_DG, cp, c_V, coord_interface, &
        ispec_is_elastic,ispec_is_acoustic_DG,&!Htabext_DG&
        mesh_zmax,tau_sigma, tau_epsilon,&
        rhoext,vpext,vsext,&
        QKappa_attenuationext,Qmu_attenuationext,gravityext,&!Nsqext,&
        nspec,coord,ibool,myrank,&
        windxext, windzext, pext_DG, gammaext_DG, etaext, muext, kappa_DG,&
        c11ext,c13ext,c15ext,c33ext,c35ext,c55ext,c12ext,c23ext,c25ext,&
        EXTERNAL_DG_ONLY_MODEL_FILENAME,&
        USE_DISCONTINUOUS_METHOD,deltat,&
        ADD_PERIODIC_CONDITIONS,& ! LNS: keep ADD_PERIODIC_CONDITIONS to throw an error if ground wind is nonzero
        IONOSPHERIC_COUPLING,Bxext, Bzext, N0ext
  
  implicit none
  
  ! Input/output.
  integer, intent(in) :: nlines_header, nlines_model
  
  ! Local variables.
  real(kind=CUSTOM_REAL), parameter :: ZERO = 0._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: ONE = 1._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: HUGEVAL = 9999._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: THRSHLD_CO2 = 1000000._CUSTOM_REAL ! Threshold for CO2 relaxation parameters.
  real(kind=CUSTOM_REAL), dimension(nlines_model) :: z_model
  real(kind=CUSTOM_REAL), dimension(nlines_model) :: density_model
  real(kind=CUSTOM_REAL), dimension(nlines_model) :: vp_model
  real(kind=CUSTOM_REAL), dimension(nlines_model) :: gravity_model
  !real(kind=CUSTOM_REAL), dimension(nlines_model) :: Nsq_model
  real(kind=CUSTOM_REAL), dimension(nlines_model) :: wx_model
  real(kind=CUSTOM_REAL), dimension(nlines_model) :: wz_model
  real(kind=CUSTOM_REAL), dimension(nlines_model) :: eta_model
  real(kind=CUSTOM_REAL), dimension(nlines_model) :: mu_model
  real(kind=CUSTOM_REAL), dimension(nlines_model) :: p_model
  !real(kind=CUSTOM_REAL), dimension(nlines_model) :: Htab_model
  real(kind=CUSTOM_REAL), dimension(nlines_model) :: kappa_model
  real(kind=CUSTOM_REAL), dimension(nlines_model) :: cp_model
  real(kind=CUSTOM_REAL), dimension(nlines_model) :: cv_model
  real(kind=CUSTOM_REAL), dimension(nlines_model) :: gamma_model
  real(kind=CUSTOM_REAL), dimension(nlines_model) :: fr_model
  real(kind=CUSTOM_REAL), dimension(nlines_model) :: svib_model
  real(kind=CUSTOM_REAL), dimension(nlines_model) :: tau_sigma_model
  real(kind=CUSTOM_REAL), dimension(nlines_model) :: tau_epsilon_model
  ! Electron density
  real(kind=CUSTOM_REAL), dimension(nlines_model) :: Bx_atmos
  real(kind=CUSTOM_REAL), dimension(nlines_model) :: Bz_atmos
  real(kind=CUSTOM_REAL), dimension(nlines_model) :: N0_atmos
  
  integer :: i, j, ispec, ii, io, indglob_DG
  integer ncolz
  real(kind=CUSTOM_REAL) dummy1, dummy2, dummy3, dummy4, dummy5, dummy6 ! Dummy reals for reading parameters which we do not care about.
  real(kind=CUSTOM_REAL) ONE_over_twopifr ! For CO2 relaxation.
  double precision :: z, frac, pii, piim1, piim2!, piip1!,tmp1, gamma_temp,gamma_temp_prev,x,max_z
  
  z_model=ZERO
  density_model=ZERO
  vp_model=ZERO
  gravity_model=ZERO
  !Nsq_model=ZERO
  wx_model=ZERO
  wz_model=ZERO
  eta_model=ZERO
  mu_model=ZERO
  p_model=ZERO
  !Htab_model=ZERO
  kappa_model=ZERO
  cp_model=ZERO
  cv_model=ZERO
  gamma_model=ZERO
  fr_model=ZERO
  svib_model=ZERO
  tau_sigma_model=ONE
  tau_epsilon_model=ONE
  
  rhoext = ZERO
  vpext = ZERO
  gravityext = ZERO
  vsext = ZERO
  Qmu_attenuationext = ZERO
  QKappa_attenuationext = ZERO
  windxext = ZERO
  pext_DG = ZERO
  gammaext_DG(indglob_DG) = ZERO
  etaext = ZERO
  muext = ZERO
  kappa_DG = ZERO
  tau_sigma = ONE ! Relaxation deactivated by default.
  tau_epsilon = ONE ! Relaxation deactivated by default.
  !if(myrank==0) then
  !  write(*,*) "********************************"
  !  write(*,*) "*           WARNING            *"
  !  write(*,*) "********************************"
  !  write(*,*) "* The support for atmospheric  *"
  !  write(*,*) "* relaxation when using        *"
  !  write(*,*) "* external DG models is not    *"
  !  write(*,*) "* implemented yet. Relaxation  *"
  !  write(*,*) "* is deactivated by setting    *"
  !  write(*,*) "* both relaxation times to 1.  *"
  !  write(*,*) "* See                          *"
  !  write(*,*) "* 'define_external_model.f90'. *"
  !  write(*,*) "********************************"
  !endif
  
  ! Safeguard.
  if(.not. USE_DISCONTINUOUS_METHOD) then
    write(*,*) "********************************"
    write(*,*) "*            ERROR             *"
    write(*,*) "********************************"
    write(*,*) "* Currently cannot use the     *"
    write(*,*) "* 'external_DG' model if not   *"
    write(*,*) "* using the DG method.         *"
    write(*,*) "********************************"
    stop
  endif
  call external_model_DG_testfile(EXTERNAL_DG_ONLY_MODEL_FILENAME) ! Test-reading the file.
  
  if(myrank==0) then
    write(*,*) "> Reading atmospheric model file '", trim(EXTERNAL_DG_ONLY_MODEL_FILENAME),&
               "' and setting the external model accordingly by performing a Lagrange interpolation."
  endif
  
  call external_model_DG_only_find_nbcols(EXTERNAL_DG_ONLY_MODEL_FILENAME, nlines_header, ncolz)
  
  if(ncolz/=17 .and. ncolz/=19 .and. ncolz/=20) then
    write(*,*) "********************************"
    write(*,*) "*            ERROR             *"
    write(*,*) "********************************"
    write(*,*) "* Atmospheric models can only  *"
    write(*,*) "* take either 17 or 19         *"
    write(*,*) "* columns. See                 *"
    write(*,*) "* 'define_external_model.f90'. *"
    write(*,*) "********************************"
    stop
  endif
  if(myrank==0) then
    write(*,*) "> Atmospheric model file contains ",ncolz," columns."
  endif
  
  ! Read and store values of model.
  OPEN(100, file=EXTERNAL_DG_ONLY_MODEL_FILENAME)
  i=1
  DO i=1,nlines_header
    ! Read and skip in header.
    READ(100,*,iostat=io)
    IF (io/=0) stop "Error reading line in atmospheric model file."
  enddo
  do i=1,nlines_model
    ! Read values.
    if(ncolz==17) then
      ! Default model, do not try to read fr and svib.
      read(100,*,iostat=io) z_model(i),density_model(i),dummy1,&
                            vp_model(i),p_model(i),dummy2,&
                            gravity_model(i),dummy6,&!Nsq_model(i),&
                            kappa_model(i),mu_model(i),dummy3,&
                            dummy4,dummy5,wx_model(i),&
                            cp_model(i),cv_model(i),gamma_model(i)
    elseif(ncolz==19) then
      ! Read also fr and svib.
      read(100,*,iostat=io) z_model(i),density_model(i),dummy1,&
                            vp_model(i),p_model(i),dummy2,&
                            gravity_model(i),dummy6,&!Nsq_model(i),&
                            kappa_model(i),mu_model(i),dummy3,&
                            dummy4,dummy5,wx_model(i),&
                            cp_model(i),cv_model(i),gamma_model(i),&
                            fr_model(i), svib_model(i)
    
      ! If columns were found, compute tau_*.
      ! If not, do nothing. It will leave tau_sigma_model=tau_epsilon_model=1.
      ONE_over_twopifr = ONE/(2.*PI*fr_model(i))
      tau_sigma_model(i) = 0.5*ONE_over_twopifr*(-svib_model(i) + sqrt(svib_model(i)**2.+4.)) ! See 10.1007/s11214-016-0324-6, equation (11).
      tau_epsilon_model(i) = tau_sigma_model(i) + svib_model(i)*ONE_over_twopifr ! See 10.1007/s11214-016-0324-6, equation (11).
      
      if(fr_model(i) >= THRSHLD_CO2*0.5/deltat) then
        write(*,*) "********************************"
        write(*,*) "*           WARNING            *"
        write(*,*) "********************************"
        write(*,*) "* In the model, there exists a *"
        write(*,*) "* f_r > THRSHLD_CO2*f_max      *"
        write(*,*) "*       = THRSHLD_CO2*0.5*/dt  *"
        write(*,*) "* (this run's highest possible *"
        write(*,*) "* frequency, with an empirical *"
        write(*,*) "* threshold). Because of this, *"
        write(*,*) "* we disable attenuation here: *"
        write(*,*) "* \tau_{\epsilon,\sigma}=1.    *"
        write(*,*) "********************************"
        write(*,*) "* (z,THRSHLD_CO2) = ", z_model(i), THRSHLD_CO2
        write(*,*) "********************************"
        tau_sigma_model(i)   = ONE
        tau_epsilon_model(i) = ONE
      endif
      
      !if(.false.) then
      !  !if(z_model(i)<=(mesh_zmax-coord_interface)) then ! Only check relaxation inside domain.
      !    if(min(tau_sigma_model(i), tau_epsilon_model(i)) < deltat*THRSHLD_CO2) then
      !      ! THRSHLD_CO2=21.4 was empirically found to work, on a case where (...)=2.
      !      ! Trying to explain it:
      !      ! if     f_r                       > 67*f_max
      !      !    <=> f_r/f_max                 > 67
      !      !    <=> (...)/(4*pi/tau) * (2*dt) > 67                 (formula above)
      !      !    <=> tau                       < 67*(...)/(2*pi)*dt
      !      !    <=> tau                       < 21.4 * dt          (if (...)\simeq2)
      !      ! then there is a problem.
      !      ! Maybe consider the substeps and not only the steps...
      !      ! TODO: clarify all that.
      !      write(*,*) "********************************"
      !      write(*,*) "*           WARNING            *"
      !      write(*,*) "********************************"
      !      write(*,*) "* In the model, there exists a *"
      !      write(*,*) "* tau < Dt*THRSHLD_CO2.      *"
      !      !write(*,*) "* z =            ",z_model(i)
      !      !write(*,*) "* \tau_\epsilon =",tau_epsilon_model(i)
      !      !write(*,*) "* \tau_\sigma =  ",tau_sigma_model(i)
      !      !write(*,*) "* Dt =           ",deltat
      !      !write(*,*) "* THRSHLD_CO2 =",THRSHLD_CO2
      !      write(*,*) "* z, \tau_{\epsilon,\sigma}, Dt, THRSHLD_CO2 = ",&
      !                 z_model(i),tau_epsilon_model(i),tau_sigma_model(i),deltat,THRSHLD_CO2
      !      !write(*,*) "* svib_factor(i) =",(-svib_model(i) + sqrt(svib_model(i)**2.+4.))
      !      write(*,*) "* This has proven not to work  *"
      !      write(*,*) "* well. We hereby force        *"
      !      write(*,*) "* \tau_{\epsilon,\sigma}=1.    *"
      !      write(*,*) "********************************"
      !      write(*,*) "* To prevent this, either make *"
      !      write(*,*) "* sure that                    *"
      !      write(*,*) "* \tau_{\epsilon,\sigma} is    *"
      !      write(*,*) "* high enough (> Dt*THRSHLD_CO2 = ",(deltat*THRSHLD_CO2),")"
      !      write(*,*) "* or make sure Dt is small     *"
      !      write(*,*) "* enough (< min(\tau_{\epsilon,\sigma})/THRSHLD_CO2 = ",&
      !                 (min(tau_sigma_model(i), tau_epsilon_model(i))/THRSHLD_CO2),")."
      !      write(*,*) "********************************"
      !      tau_sigma_model(i)=ONE
      !      tau_epsilon_model(i)=ONE
      !    endif
      !  !endif
      !endif
    elseif(ncolz==20) then
      ! Read also fr and svib.
      read(100,*,iostat=io) z_model(i),density_model(i),dummy1,&
                            vp_model(i),p_model(i),dummy2,&
                            gravity_model(i),dummy6,&!Nsq_model(i),&
                            kappa_model(i),mu_model(i),dummy3,&
                            dummy4,wx_model(i),dummy5,&
                            cp_model(i),cv_model(i),gamma_model(i),&
                            !fr_model(i), svib_model(i),&
                            Bx_atmos(i), Bz_atmos(i), N0_atmos(i)

      tau_epsilon_model(i) = 1.
      tau_sigma_model(i)   = 1.

      !if(density_model(i) < 0.) then
      !WRITE(*,*) "density_model",i, z_model(i), density_model(i)
      !stop 'ToTO'
      !endif

    endif
    
    ! Compute eta, in all cases.
    eta_model(i) = (FOUR_THIRDS) * mu_model(i)
    
    IF(io/=0) EXIT
    !if(.true.) then ! DEBUG
    !  write(*,*) 'z_model(i),density_model(i)',&!,tmp1',&
    !             'vp_model(i),p_model(i),',&!Htab_model(i)',&
    !             'gravity_model(i),',&!Nsq_model(i)',&
    !             'kappa_model(i),mu_model(i),dummy1,',&
    !             'dummy2, dummy3, wx_model(i),',&
    !             'cp_model(i),cv_model(i),gamma_model(i)',&
    !             'fr_model(i), svib_model(i)'
    !  write(*,*) z_model(i),density_model(i),&!,tmp1,&
    !             vp_model(i),p_model(i),&!Htab_model(i),&
    !             gravity_model(i),&!Nsq_model(i),&
    !             kappa_model(i),mu_model(i),dummy1,&
    !             dummy2, dummy3, wx_model(i),&
    !             cp_model(i),cv_model(i),gamma_model(i),&
    !             fr_model(i), svib_model(i)
    !endif
  ENDDO 
  CLOSE(100)

  ! Quickly check values of the model.
  if(myrank==0) then
    if(mesh_zmax-coord_interface>z_model(nlines_model)-z_model(1)) then
      write(*,*) "********************************"
      write(*,*) "*           WARNING            *"
      write(*,*) "********************************"
      write(*,*) "* The mesh vertically spans    *"
      write(*,*) "* more than the model gives.   *"
      write(*,*) "* In other words, the          *"
      write(*,*) "* atmospheric part of the mesh *"
      write(*,*) "* is thicker than what is      *"
      write(*,*) "* given in the model file.     *"
      write(*,*) "*  mesh z_max:             ", mesh_zmax
      write(*,*) "*  mesh altitude 0:        ", coord_interface
      write(*,*) "*  mesh vertical span:     ", mesh_zmax-coord_interface
      write(*,*) "*  model lowest altitude:  ", z_model(1)
      write(*,*) "*  model highest altitude: ", z_model(nlines_model)
      write(*,*) "*  model vertical span:    ", z_model(nlines_model)-z_model(1)
      write(*,*) "* Interpolation could work,    *"
      write(*,*) "* but be careful. See          *"
      write(*,*) "* 'define_external_model.f90'. *"
      write(*,*) "********************************"
    endif
    if(z_model(nlines_model)-z_model(1)>10.*(mesh_zmax-coord_interface)) then
      write(*,*) "********************************"
      write(*,*) "*           WARNING            *"
      write(*,*) "********************************"
      write(*,*) "* The model gives more than 10 *"
      write(*,*) "* times as much as the mesh    *"
      write(*,*) "* vertically spans. Consider   *"
      write(*,*) "* generating a model over a    *"
      write(*,*) "* smaller vertical span.       *"
      write(*,*) "*  mesh z_max:             ", mesh_zmax
      write(*,*) "*  mesh altitude 0:        ", coord_interface
      write(*,*) "*  mesh vertical span:     ", mesh_zmax-coord_interface
      write(*,*) "*  model lowest altitude:  ", z_model(1)
      write(*,*) "*  model highest altitude: ", z_model(nlines_model)
      write(*,*) "*  model vertical span:    ", z_model(nlines_model)-z_model(1)
      write(*,*) "* If the vertical step of the  *"
      write(*,*) "* model is not small enough,   *"
      write(*,*) "* interpolation could lead to  *"
      write(*,*) "* absurd values. This could    *"
      write(*,*) "* work, but be careful. See    *"
      write(*,*) "* 'define_external_model.f90'. *"
      write(*,*) "********************************"
    endif
  endif

  ! THIS CONDITION SHOULD REMAIN ON LNS BRANCH ONLY, BECAUSE QUENTIN WANTS TO REMOVE IT FROM MASTER BRANCH.
  ! ON LNS BRANCH, I LEAVE IT HERE, BECAUSE I PREFER TO DO SO.
  if(.not. ADD_PERIODIC_CONDITIONS) then
    do i=1, nlines_model
      if(myrank==0 .and. abs(z_model(i))==0. .and. abs(wx_model(i))>0.) then
        write(*,*) "********************************"
        write(*,*) "*            ERROR             *"
        write(*,*) "********************************"
        write(*,*) "* In the loaded model, wind is *"
        write(*,*) "* non-zero at altitude 0. This *"
        write(*,*) "* is not physically            *"
        write(*,*) "* acceptable. Activating       *"
        write(*,*) "* horizontal periodic boundary *"
        write(*,*) "* conditions can make this     *"
        write(*,*) "* work, but this will remain   *"
        write(*,*) "* questionnable.               *"
        write(*,*) "********************************"
        stop
      endif
    enddo
  endif
  if(minval(density_model)<=ZERO) then
    write(*,*) "********************************"
    write(*,*) "*            ERROR             *"
    write(*,*) "********************************"
    write(*,*) "* A negative density was       *"
    write(*,*) "* found in the model.          *"
    write(*,*) "********************************"
    stop
  endif
  if(minval(p_model)<=ZERO) then
    write(*,*) "********************************"
    write(*,*) "*            ERROR             *"
    write(*,*) "********************************"
    write(*,*) "* A negative pressure was      *"
    write(*,*) "* found in the model.          *"
    write(*,*) "********************************"
    stop
  endif
  if(minval(cp_model)<=ZERO .or. minval(cv_model)<=ZERO) then
    write(*,*) "********************************"
    write(*,*) "*            ERROR             *"
    write(*,*) "********************************"
    write(*,*) "* A negative c_p or c_v was    *"
    write(*,*) "* found in the model.          *"
    write(*,*) "********************************"
    stop
  endif
  
  ! TODO: Correctly affect those variables.
  cp = cp_model(1)
  c_V = cv_model(1)
  
  rhoext = ZERO
  vpext = ZERO
  vsext = ZERO
  QKappa_attenuationext = ZERO
  Qmu_attenuationext = ZERO
  c11ext = ZERO
  c13ext = ZERO
  c15ext = ZERO
  c33ext = ZERO
  c35ext = ZERO
  c55ext = ZERO
  c12ext = ZERO
  c23ext = ZERO
  c25ext = ZERO
  gravityext = ZERO
  windxext = ZERO
  windzext = ZERO
  pext_DG = ZERO
  gammaext_DG = ZERO
  !Htabext_DG = ZERO
  etaext = ZERO
  muext = ZERO
  kappa_DG = ZERO
  if(IONOSPHERIC_COUPLING) then
  Bxext = 1.
  Bzext = 1.
  N0ext = 0.
  endif
  
  do ispec = 1, nspec
    if(ispec_is_elastic(ispec)) then
      cycle ! Elastic regions are treated separately (see call to 'external_DG_update_elastic_from_parfile' below).
    else if(ispec_is_acoustic_DG(ispec)) then
      ! For DG elements, go through GLL points one by one.
      do j = 1, NGLLZ
        do i = 1, NGLLX
          indglob_DG = ibool_DG(i, j, ispec)
          
          z = coord(2, ibool(i, j, ispec)) ! Get z-coordinate of GLL point.
          z = z - coord_interface ! Get altitude of GLL point.
          
          ! Find first index ii (from bottom) in model such that altitude_model(ii) > z.
          ii = 1
          do while(z >= z_model(ii) .and. ii /= nlines_model)
            ii = ii + 1
          enddo
          
          if (ii == 1) then
            ! Altitude of point is <= the altitude of first line of model.
            ! Use first line of model.
            
            if(.false.) then ! DEBUG
              write(*,*) density_model(1)
              write(*,*) vp_model(1)
              write(*,*) gravity_model(1)
              write(*,*) wx_model(1)
              write(*,*) wz_model(1)
              write(*,*) p_model(1)
              write(*,*) gamma_model(1)
              !write(*,*) Htab_model(1)
              write(*,*) eta_model(1)
              write(*,*) mu_model(1)
              write(*,*) kappa_model(1)
              write(*,*) rhoext(i, j, ispec)
              write(*,*) vpext(i, j, ispec)
              write(*,*) gravityext(i, j, ispec)
              write(*,*) vsext(i, j, ispec)
              write(*,*) Qmu_attenuationext(i, j, ispec)
              write(*,*) QKappa_attenuationext(i, j, ispec)
              write(*,*) windxext(i, j, ispec)
              write(*,*) windzext(i, j, ispec)
              write(*,*) pext_DG(i, j, ispec)
              write(*,*) "gamma", gammaext_DG(indglob_DG)
              !write(*,*) "Htab", Htabext_DG(indglob_DG)
              write(*,*) etaext(i, j, ispec)
              write(*,*) muext(i, j, ispec)
              write(*,*) kappa_DG(i, j, ispec)
            endif
            
            rhoext(i, j, ispec) = density_model(1)
            vpext(i, j, ispec) = vp_model(1)
            gravityext(i, j, ispec) = gravity_model(1)
            !Nsqext(i, j, ispec) = Nsq_model(1)
            vsext(i, j, ispec) = ZERO
            Qmu_attenuationext(i, j, ispec) = HUGEVAL
            QKappa_attenuationext(i, j, ispec) = HUGEVAL
            windxext(i, j, ispec) = wx_model(1)
            !windzext(i, j, ispec) = wz_model(1)
            pext_DG(i, j, ispec) = p_model(1)
            gammaext_DG(indglob_DG) = gamma_model(1)
            !Htabext_DG(indglob_DG) = Htab_model(1)
            etaext(i, j, ispec) = eta_model(1)
            muext(i, j, ispec) = mu_model(1)
            kappa_DG(i, j, ispec) = kappa_model(1)
            tau_sigma(i, j, ispec) = tau_sigma_model(i)
            tau_epsilon(i, j, ispec) = tau_epsilon_model(i)

            if(IONOSPHERIC_COUPLING) then
            Bxext(ibool_DG(i,j,ispec)) = Bx_atmos(1)
            Bzext(ibool_DG(i,j,ispec))  = Bz_atmos(1)
            N0ext(i,j,ispec)  = N0_atmos(1)
            endif
          elseif (ii == 2) then
            ! Altitude of point is > 1st line of model, and <= 2nd line of model.
            ! Interpolate using the values at lines 1 (ii-1) and 2 (ii) of model.
            
            frac = (z-z_model(ii-1))/(z_model(ii)-z_model(ii-1))

            rhoext(i, j, ispec)          = density_model(ii-1)     + frac*(density_model(ii)-density_model(ii-1))
            !rhoext(i, j, ispec)          = exp(log(density_model(ii-1)) + frac*(log(density_model(ii))-log(density_model(ii-1)))) ! DEBUG: exponential interpolation
            vpext(i, j, ispec)           = vp_model(ii-1)          + frac*(vp_model(ii)-vp_model(ii-1))
            gravityext(i, j, ispec)      = gravity_model(ii-1)     + frac*(gravity_model(ii)-gravity_model(ii-1))
            !Nsqext(i, j, ispec)          = Nsq_model(ii-1) + frac*(Nsq_model(ii)-Nsq_model(ii-1))
            vsext(i, j, ispec)           = ZERO
            Qmu_attenuationext(i, j, ispec) = HUGEVAL
            QKappa_attenuationext(i, j, ispec) = HUGEVAL
            windxext(i, j, ispec)        = wx_model(ii-1)          + frac*(wx_model(ii)-wx_model(ii-1))
            !windzext(i, j, ispec)        = wz_model(ii-1)          + frac*(wz_model(ii)-wz_model(ii-1))
            pext_DG(i, j, ispec)         = p_model(ii-1)           + frac*(p_model(ii)-p_model(ii-1))
            gammaext_DG(indglob_DG)      = gamma_model(ii-1)       + frac*(gamma_model(ii)-gamma_model(ii-1))
            etaext(i, j, ispec)          = eta_model(ii-1)         + frac*(eta_model(ii)-eta_model(ii-1))
            muext(i, j, ispec)           = mu_model(ii-1)          + frac*(mu_model(ii)-mu_model(ii-1))
            !Htabext_DG(indglob_DG)      = Htab_model(ii-1) + frac*(Htab_model(ii)-Htab_model(ii-1))
            kappa_DG(i, j, ispec)        = kappa_model(ii-1)       + frac*(kappa_model(ii)-kappa_model(ii-1))
            tau_sigma(i, j, ispec)       = tau_sigma_model(ii-1)   + frac*(tau_sigma_model(ii)-tau_sigma_model(ii-1))
            tau_epsilon(i, j, ispec)     = tau_epsilon_model(ii-1) + frac*(tau_epsilon_model(ii)-tau_epsilon_model(ii-1))
            
            if(IONOSPHERIC_COUPLING) then
              Bxext(ibool_DG(i, j, ispec)) = Bx_atmos(ii-1)          + frac*(Bx_atmos(ii)-Bx_atmos(ii-1))
              Bzext(ibool_DG(i, j, ispec)) = Bz_atmos(ii-1)          + frac*(Bz_atmos(ii)-Bz_atmos(ii-1))
              N0ext(i, j, ispec)           = N0_atmos(ii-1)          + frac*(N0_atmos(ii)-N0_atmos(ii-1))
            endif

          else
            ! Altitude of point is > (ii-1)-th line of model, and <= (ii)-th line of model.
            ! Interpolate using the values at lines (ii-2), (ii-1), and ii of model.
            
            !if(.false.) then ! DEBUG
            !pii   = (z-z_model(ii+1))*(z-z_model(ii-1))*(z-z_model(ii-2))&
            !/((z_model(ii)-z_model(ii-1))*(z_model(ii)-z_model(ii-2))*(z_model(ii)-z_model(ii+1)))
            !piim1 = (z-z_model(ii+1))*(z-z_model(ii))*(z-z_model(ii-2))&
            !/((z_model(ii-1)-z_model(ii))*(z_model(ii-1)-z_model(ii-2))*(z_model(ii-1)-z_model(ii+1)))
            !piim2 = (z-z_model(ii+1))*(z-z_model(ii-1))*(z-z_model(ii))&
            !/((z_model(ii-2)-z_model(ii))*(z_model(ii-2)-z_model(ii-1))*(z_model(ii-2)-z_model(ii+1)))
            !piip1 = (z-z_model(ii-2))*(z-z_model(ii-1))*(z-z_model(ii))&
            !/((z_model(ii+1)-z_model(ii))*(z_model(ii+1)-z_model(ii-1))*(z_model(ii+1)-z_model(ii-2)))
            !endif

            pii   = (z-z_model(ii-1))*(z-z_model(ii-2))/((z_model(ii)-z_model(ii-1))*(z_model(ii)-z_model(ii-2)))
            piim1 = (z-z_model(ii))  *(z-z_model(ii-2))/((z_model(ii-1)-z_model(ii))*(z_model(ii-1)-z_model(ii-2)))
            piim2 = (z-z_model(ii-1))*(z-z_model(ii))  /((z_model(ii-2)-z_model(ii))*(z_model(ii-2)-z_model(ii-1)))

            rhoext(i, j, ispec)          = density_model(ii)*pii     + density_model(ii-1)*piim1 + density_model(ii-2)*piim2
            vpext(i, j, ispec)           = vp_model(ii)*pii          + vp_model(ii-1)*piim1 + vp_model(ii-2)*piim2
            gravityext(i, j, ispec)      = gravity_model(ii)*pii     + gravity_model(ii-1)*piim1 + gravity_model(ii-2)*piim2
            !Nsqext(i, j, ispec)          = Nsq_model(ii)*pii + Nsq_model(ii-1)*piim1 + Nsq_model(ii-2)*piim2
            vsext(i, j, ispec)           = ZERO
            Qmu_attenuationext(i, j, ispec) = HUGEVAL
            QKappa_attenuationext(i, j, ispec) = HUGEVAL
            windxext(i, j, ispec)        = wx_model(ii)*pii          + wx_model(ii-1)*piim1 + wx_model(ii-2)*piim2
            !windzext(i, j, ispec)        = wz_model(ii)*pii          + wz_model(ii-1)*piim1 + wz_model(ii-2)*piim2
            pext_DG(i, j, ispec)         = p_model(ii)*pii           + p_model(ii-1)*piim1 + p_model(ii-2)*piim2
            gammaext_DG(indglob_DG)      = gamma_model(ii)*pii       + gamma_model(ii-1)*piim1 + gamma_model(ii-2)*piim2
            etaext(i, j, ispec)          = eta_model(ii)*pii         + eta_model(ii-1)*piim1 + eta_model(ii-2)*piim2
            muext(i, j, ispec)           = mu_model(ii)*pii          + mu_model(ii-1)*piim1 + mu_model(ii-2)*piim2
            !Htabext_DG(indglob_DG)       = Htab_model(ii)*pii        + Htab_model(ii-1)*piim1 + Htab_model(ii-2)*piim2
            kappa_DG(i, j, ispec)        = kappa_model(ii)*pii       + kappa_model(ii-1)*piim1 + kappa_model(ii-2)*piim2
            tau_sigma(i, j, ispec)       = tau_sigma_model(ii)*pii   + tau_sigma_model(ii-1)*piim1 + tau_sigma_model(ii-2)*piim2
            tau_epsilon(i, j, ispec)     = tau_epsilon_model(ii)*pii + tau_epsilon_model(ii-1)*piim1 + tau_epsilon_model(ii-2)*piim2

            if(IONOSPHERIC_COUPLING) then
              Bxext(ibool_DG(i, j, ispec)) = Bx_atmos(ii)*pii          + Bx_atmos(ii-1)*piim1 + Bx_atmos(ii-2)*piim2
              Bzext(ibool_DG(i, j, ispec)) = Bz_atmos(ii)*pii          + Bz_atmos(ii-1)*piim1 + Bz_atmos(ii-2)*piim2
              N0ext(i, j, ispec)           = N0_atmos(ii)*pii          + N0_atmos(ii-1)*piim1 + N0_atmos(ii-2)*piim2
            endif

          endif ! Endif on ii.
          
          if(rhoext(i, j, ispec) <= TINYVAL) then
            write(*,*) "********************************"
            write(*,*) "*            ERROR             *"
            write(*,*) "********************************"
            write(*,*) "* A negative (<=0) density was *"
            write(*,*) "* found after interpolation.   *"
            write(*,*) "********************************"
            write(*,*) ispec,i,j,ii,&
                       coord(1,indglob_DG),coord(2,indglob_DG),&
                       rhoext(i, j, ispec)
            write(*,*) "********************************"
            stop
          endif
          if(tau_sigma(i, j, ispec) <= TINYVAL .or. tau_epsilon(i, j, ispec) <= TINYVAL) then
            write(*,*) "********************************"
            write(*,*) "*            ERROR             *"
            write(*,*) "********************************"
            write(*,*) "* A negative (<=0) relaxation  *"
            write(*,*) "* time (tau_sigma or           *"
            write(*,*) "* tau_epsilon) was found after *"
            write(*,*) "* interpolation.               *"
            write(*,*) "********************************"
            write(*,*) ispec,i,j,ii,&
                       coord(1,indglob_DG),coord(2,indglob_DG),&
                       tau_sigma(i, j, ispec),&
                       tau_epsilon(i, j, ispec)
            write(*,*) "********************************"
            stop
          endif
          !if(Nsqext(i, j, ispec) <= ZERO) then
          !  write(*,*) "********************************"
          !  write(*,*) "*            ERROR             *"
          !  write(*,*) "********************************"
          !  write(*,*) "* A negative Nsq (Nsquare) was *"
          !  write(*,*) "* found in the model.          *"
          !  write(*,*) "********************************"
          !  write(*,*) ispec,i,j,ii,&
          !             coord(1,indglob_DG),coord(2,indglob_DG),&
          !             Nsqext(i, j, ispec),gravityext(i, j, ispec),vpext(i, j, ispec)
          !  write(*,*) "********************************"
          !  stop
          !endif
        enddo ! Enddo on i.
      enddo ! Enddo on j.
    else
      write(*,*) "********************************"
      write(*,*) "*            ERROR             *"
      write(*,*) "********************************"
      write(*,*) "* The external model loading   *"
      write(*,*) "* is not supported for         *"
      write(*,*) "* elements which are neither   *"
      write(*,*) "* DG acoustic nor elastic.     *"
      write(*,*) "*  ispec: ", ispec
      write(*,*) "*  x: ", coord(1, ibool(i, j, ispec))
      write(*,*) "*  z: ", coord(2, ibool(i, j, ispec))
      write(*,*) "********************************"
      stop
    endif ! Endif on ispec_is_acoustic_DG.
  enddo ! Enddo on ispec.
  
  ! CHECKING
  !if(.true.) then ! DEBUG
  !  do ispec = 1, nspec
  !    !do j = 1, NGLLZ
  !    do j = 1, 1
  !      !do i = 1, NGLLX
  !      do i = 1, 1
  !        if( &
  !            !      abs(coord(2, ibool(i, j, ispec)))<=20. &
  !            !.and. &
  !                 abs(coord(1, ibool(i, j, ispec)))<=1e-2) then
  !          write(*,*) coord(1, ibool(i, j, ispec)), coord(2, ibool(i, j, ispec)),&
  !                     ispec_is_acoustic_DG(ispec), ispec_is_elastic(ispec),&
  !                     "rho", rhoext(i, j, ispec),&
  !                     "vp", vpext(i, j, ispec),&
  !                     "g", gravityext(i, j, ispec),&
  !                     "vs", vsext(i, j, ispec),&
  !                     "Qm", Qmu_attenuationext(i, j, ispec),&
  !                     "Qk", QKappa_attenuationext(i, j, ispec),&
  !                     "wx", windxext(i, j, ispec),&
  !                    "ceff_right", vpext(i, j, ispec)+windxext(i, j, ispec),&
  !                     "wz", windzext(i, j, ispec),&
  !                     "p", pext_DG(i, j, ispec),&
  !                     "gamma", gammaext_DG(ibool_DG(i, j, ispec)),&
  !                     !"Htab", Htabext_DG(ibool_DG(i, j, ispec)),&
  !                     "eta", etaext(i, j, ispec),&
  !                     "mu", muext(i, j, ispec),&
  !                     "kap", kappa_DG(i, j, ispec),&
  !                     "tau_sig", tau_sigma(i, j, ispec),&
  !                     "tau_eps", tau_epsilon(i, j, ispec)
  !                     ""
  !        endif
  !      enddo
  !    enddo
  !  enddo ! Enddo on ispec.
  !  stop
  !endif
  
  if(IONOSPHERIC_COUPLING) then
    where(abs(Bxext) > 0._CUSTOM_REAL) Bxext = Bxext / sqrt(Bxext**2 + Bzext**2)
    where(abs(Bzext) > 0._CUSTOM_REAL) Bzext = Bzext / sqrt(Bxext**2 + Bzext**2)
  endif

  call external_DG_update_elastic_from_parfile() ! Update elastic regions by reading parameters directly from parfile.
  !call define_external_model_PREM()
  
  if(myrank==0) then
    write(*,*) "********************************"
    write(*,*) "*           WARNING            *"
    write(*,*) "********************************"
    write(*,*) "* When using external DG       *"
    write(*,*) "* models without periodic      *"
    write(*,*) "* conditions, non-zero         *"
    write(*,*) "* horizontal wind at           *"
    write(*,*) "* fluid-solid interface        *"
    write(*,*) "* produces instabilities at    *"
    write(*,*) "* the fluid/solid boundary.    *"
    write(*,*) "********************************"
  endif
  
  if(myrank==0) then
    write(*,*) "> Done reading atmospheric model file and setting the external model accordingly."
  endif
  
end subroutine define_external_model_DG_only

! ------------------------------------------------------------ !
! external_DG_update_elastic_from_parfile                      !
! ------------------------------------------------------------ !
! Only used for trim(MODEL)=='external_DG': sets the variables for elastic model according to what is specified in parfile.

subroutine external_DG_update_elastic_from_parfile()

  use constants,only: CUSTOM_REAL,NGLLX,NGLLZ,NDIM,IMAIN
  use specfem_par,only: &
        poroelastcoef,density,kmato, ispec_is_elastic,&
        rhoext,vpext,vsext,&
        QKappa_attenuationext,Qmu_attenuationext,&
        nspec,myrank,&
        QKappa_attenuation,Qmu_attenuation,anisotropy,&
        c11ext,c13ext,c15ext,c33ext,c35ext,c55ext,c12ext,c23ext,c25ext
  
  implicit none
  
  ! Local variables.
  real(kind=CUSTOM_REAL), parameter :: ZERO = 0._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: ONE = 1._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: HUGEVAL = 9999._CUSTOM_REAL
  integer :: i, j, ispec
  
  if(myrank==0) then
    write(*,*) "> > Updating elastic part of external model: read values from ",&
               "materials variables (which initially come from parfile) to fill the 'external' variables."
  endif
  
  do ispec = 1, nspec
    if(ispec_is_elastic(ispec)) then
      ! Use quantities loaded from parfile to fill "external" model for elastic media.
      rhoext(:, :, ispec) = density(1,kmato(ispec))
      vpext(:, :, ispec) = sqrt(poroelastcoef(3,1,kmato(ispec))/density(1,kmato(ispec))) ! Since poroelastcoef(3,1,n) = lambda+2mu = rho cp^2 (see gmat01.f90).
      vsext(:, :, ispec) = sqrt(poroelastcoef(2,1,kmato(ispec))/density(1,kmato(ispec))) ! Since poroelastcoef(2,1,n) = mu = rho vs^2 (see gmat01.f90).
      QKappa_attenuationext(:, :, ispec) = QKappa_attenuation(kmato(ispec))
      Qmu_attenuationext(:, :, ispec) = Qmu_attenuation(kmato(ispec))
      c11ext(:, :, ispec) = anisotropy(1,kmato(ispec))
      c13ext(:, :, ispec) = anisotropy(2,kmato(ispec))
      c15ext(:, :, ispec) = anisotropy(3,kmato(ispec))
      c33ext(:, :, ispec) = anisotropy(4,kmato(ispec))
      c35ext(:, :, ispec) = anisotropy(5,kmato(ispec))
      c55ext(:, :, ispec) = anisotropy(6,kmato(ispec))
      c12ext(:, :, ispec) = anisotropy(7,kmato(ispec))
      c23ext(:, :, ispec) = anisotropy(8,kmato(ispec))
      c25ext(:, :, ispec) = anisotropy(9,kmato(ispec))
      
      !write(*,*) 'lul' ! DEBUG
      if(.false.) then ! DEBUG
        do j = 1, NGLLZ
          do i = 1, NGLLX
            rhoext(i,j,ispec) = 2500.d0
            vpext(i,j,ispec) = 3600.d0
            vsext(i,j,ispec) = vpext(i,j,ispec) / 2.d0
            QKappa_attenuationext(i,j,ispec) = 120.
            Qmu_attenuationext(i,j,ispec) = 120.
            c11ext(i,j,ispec) = 0.d0   ! this means no anisotropy
            c13ext(i,j,ispec) = 0.d0
            c15ext(i,j,ispec) = 0.d0
            c33ext(i,j,ispec) = 0.d0
            c35ext(i,j,ispec) = 0.d0
            c55ext(i,j,ispec) = 0.d0
            c12ext(i,j,ispec) = 0.d0
            c23ext(i,j,ispec) = 0.d0
            c25ext(i,j,ispec) = 0.d0
            !! AB AB Do not forget these 3 lines otherwise PML may not work !!
            density(1,kmato(ispec)) = rhoext(i,j,ispec)
            poroelastcoef(3,1,kmato(ispec)) = rhoext(i,j,ispec) * vpext(i,j,ispec) * vpext(i,j,ispec)
            poroelastcoef(2,1,kmato(ispec)) =  rhoext(i,j,ispec) * vsext(i,j,ispec) * vsext(i,j,ispec)
          enddo
        enddo
      endif
    else
      cycle
    endif  ! Endif on ispec_is_elastic.
  enddo ! Enddo on ispec.
end subroutine external_DG_update_elastic_from_parfile

! ------------

   subroutine define_external_model_PREM()

  use constants,only: CUSTOM_REAL,NGLLX,NGLLZ,NDIM,IMAIN
  use specfem_par,only: poroelastcoef,density,kmato, ispec_is_elastic, &
        coord,ibool, rhoext,vpext,vsext,QKappa_attenuationext,Qmu_attenuationext, &
        c11ext,c13ext,c15ext,c33ext,c35ext,c55ext,c12ext,c23ext,c25ext,nspec

  implicit none

!--------------------------------------------------------------------------------------------------
!
!          Model is either test one or being based on the Nrl_MSISE2000 atmosphere model
!          ----------------------------------------------------------

! 1D model:
!
! The 1D model will be used to create a tabular model of atmosphere

!--------------------------------------------------------------------------------------------------

  integer, parameter :: NR_LAYER = 198
  real(kind=CUSTOM_REAL), dimension(NR_LAYER) :: zz,rho_aux,vp_aux,vs_aux,Qkappa_aux,Qs_aux

! number of layers in the model

  integer :: i,j,ispec,iglob,k

  double precision :: x,z
  
  open(10,file='./PREM.txt', form='formatted')

  do i = 1,NR_LAYER
  read(10,*) zz(i),rho_aux(i),vp_aux(i),vs_aux(i),Qkappa_aux(i),Qs_aux(i)
  !WRITE(*,*) ">>>> z ", zz(i),rho_aux(i),vp_aux(i),vs_aux(i),Qkappa_aux(i),Qs_aux(i)
  enddo
  close(10)
  
! loop on all the elements of the mesh, and inside each element loop on all the GLL points
  do ispec = 1,nspec

    ! Recover current element's max height along z to discriminate between regions
    do j = 1,NGLLZ
      do i = 1,NGLLX
      
       iglob = ibool(i,j,ispec)

       x = coord(1,iglob)
       z = -coord(2,iglob)
       
       k = 1
       do while (zz(k) < z .AND. k < NR_LAYER)
        k = k + 1
       enddo

   if(ispec_is_elastic(ispec)) then
   
        ! From external file
        density(1,kmato(ispec)) = rho_aux(k)
        rhoext(i,j,ispec)          = rho_aux(k)
        vpext(i,j,ispec)      = vp_aux(k)
        vsext(i,j,ispec)      = vs_aux(k)
        QKappa_attenuationext(i,j,ispec) = Qkappa_aux(k)
        Qmu_attenuationext(i,j,ispec)    = Qs_aux(k)
        
        poroelastcoef(3,1,kmato(ispec)) = rhoext(i,j,ispec) * vpext(i,j,ispec) * vpext(i,j,ispec)
        poroelastcoef(2,1,kmato(ispec)) =  rhoext(i,j,ispec) * vsext(i,j,ispec) * vsext(i,j,ispec)

     endif
  
      enddo
    enddo
    
  enddo

! no anisotropy
  c11ext(:,:,:) = 0.d0
  c13ext(:,:,:) = 0.d0
  c15ext(:,:,:) = 0.d0
  c33ext(:,:,:) = 0.d0
  c35ext(:,:,:) = 0.d0
  c55ext(:,:,:) = 0.d0
  c12ext(:,:,:) = 0.d0
  c23ext(:,:,:) = 0.d0
  c25ext(:,:,:) = 0.d0

  end subroutine define_external_model_PREM

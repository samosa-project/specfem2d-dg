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

  subroutine define_derivation_matrices()

  use specfem_par, only: xigll,zigll,wxgll,wzgll,hprime_xx,hprime_zz,hprimewgll_xx,hprimewgll_zz

  use constants,only: GAUSSALPHA,GAUSSBETA,NGLLX,NGLLZ,ZERO

  implicit none

! function for calculating derivatives of Lagrange polynomials
  double precision, external :: lagrange_deriv_GLL
  !double precision, dimension(NGLLX) :: h, hprime

  ! Local
  integer :: i1,i2,k1,k2

  ! set up coordinates of the Gauss-Lobatto-Legendre points
  call zwgljd(xigll,wxgll,NGLLX,GAUSSALPHA,GAUSSBETA)
  call zwgljd(zigll,wzgll,NGLLZ,GAUSSALPHA,GAUSSBETA)

  ! if number of points is odd, the middle abscissa is exactly zero
  if (mod(NGLLX,2) /= 0) xigll((NGLLX-1)/2+1) = ZERO
  if (mod(NGLLZ,2) /= 0) zigll((NGLLZ-1)/2+1) = ZERO
  
  ! TODO: Two loops are needed if NGLLX != NGLLZ. Maybe consider merging them since NGLLX will always equal NGLLZ.
  
  ! Calculate derivatives of the Lagrange polynomials and precalculate some products in double precision:
  ! hprime(i,j) = h'_j(xigll_i) (by definition of the derivation matrix).
  do i1 = 1,NGLLX
    do i2 = 1,NGLLX
      hprime_xx(i2,i1) = lagrange_deriv_GLL(i1-1,i2-1,xigll,NGLLX)
      !call lagrange_any(xigll(i2),NGLLX,xigll,h,hprime) ! DEBUG
      !WRITE(*,*) "TEST LAGRANGE h", h(i2) ! DEBUG
      !WRITE(*,*) "TEST LAGRANGE hprime", hprime(i1) ! DEBUG
      !WRITE(*,*) "TEST LAGRANGE hprime_xx(i2,i1)", hprime_xx(i2,i1) ! DEBUG
      !WRITE(*,*) "*********" ! DEBUG
      hprimewgll_xx(i2,i1) = wxgll(i2) * hprime_xx(i2,i1)
    enddo
  enddo

  do k1 = 1,NGLLZ
    do k2 = 1,NGLLZ
      hprime_zz(k2,k1) = lagrange_deriv_GLL(k1-1,k2-1,zigll,NGLLZ)
      hprimewgll_zz(k2,k1) = wzgll(k2) * hprime_zz(k2,k1)
    enddo
  enddo

  end subroutine define_derivation_matrices

!____________________________________________________________________________________
!
!subroutine define_GLJ_derivation_matrix(xiglj,wxglj,hprimeBar_xx,hprimeBarwglj_xx)
! Calculate all that we need for the GLJ quadrature on axial elements :
! Weights, GLJ points and derivatives of polynomials at GLJ points.
!____________________________________________________________________________________
!

  subroutine define_GLJ_derivation_matrix()

  use constants,only: NGLJ

  use specfem_par, only: xiglj,wxglj,hprimeBar_xx,hprimeBarwglj_xx

  implicit none

  double precision, parameter    :: alphaGLJ = 0.d0,betaGLJ = 1.d0

! function for calculating derivatives of GLJ polynomials
  double precision, external :: poly_deriv_GLJ

  integer i1,i2

! set up coordinates of the Gauss-Lobatto-Jacobi points
  call zwgljd(xiglj,wxglj,NGLJ,alphaGLJ,betaGLJ)

! calculate derivatives of the GLJ quadrature polynomials
! and precalculate some products in double precision
! hprimeBar(i,j) = hBar'_j(xiglj_i) by definition of the derivation matrix
  do i1 = 1,NGLJ
    do i2 = 1,NGLJ
      hprimeBar_xx(i2,i1) = poly_deriv_GLJ(i1-1,i2-1,xiglj,NGLJ)
      hprimeBarwglj_xx(i2,i1) = wxglj(i2) * hprimeBar_xx(i2,i1)
    enddo
  enddo

  end subroutine define_GLJ_derivation_matrix


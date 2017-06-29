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

! for acoustic solver

  subroutine compute_coupling_acoustic_DG_el(i_current, j_current, ispec_current, &
                timelocal, rho_DG_P, rhovx_DG_P, &
                rhovz_DG_P, E_DG_P, &
                veloc_x_DG_P, veloc_z_DG_P, p_DG_P)

  use constants,only: CUSTOM_REAL,NGLLX,NGLLZ,NGLJ,NDIM, &
    CPML_X_ONLY,CPML_Z_ONLY,IRIGHT,ILEFT,IBOTTOM,ITOP,ONE,ALPHA_LDDRK,BETA_LDDRK

  use specfem_par, only: num_fluid_solid_edges,ibool,wxgll,wzgll,xix,xiz,&
                         gammax,gammaz,jacobian,ivalue,jvalue,ivalue_inverse,jvalue_inverse,&
                         fluid_solid_acoustic_ispec,fluid_solid_acoustic_iedge, &
                         fluid_solid_elastic_ispec,fluid_solid_elastic_iedge,&
                         AXISYM,coord,is_on_the_axis,xiglj,wxglj,&
                         rmemory_fsb_displ_elastic,timeval,deltat,&
                         rmemory_fsb_displ_elastic_LDDRK,i_stage,stage_time_scheme, &
                         nglob_acoustic,nglob_elastic

  ! PML arrays
  use specfem_par, only: PML_BOUNDARY_CONDITIONS,ispec_is_PML,nspec_PML,spec_to_PML,region_CPML, &
                         K_x_store,K_z_store,d_x_store,d_z_store,alpha_x_store,alpha_z_store

  implicit none

  real(kind=CUSTOM_REAL),dimension(NDIM,nglob_elastic) :: displ_elastic,displ_elastic_old
  real(kind=CUSTOM_REAL),dimension(nglob_acoustic) :: potential_dot_dot_acoustic

  !local variable
  real(kind=CUSTOM_REAL), dimension(NGLJ,NGLLZ) :: r_xiplus1
  integer :: inum,ispec_acoustic,ispec_elastic,iedge_acoustic,iedge_elastic,ipoin1D,i,j,iglob
  real(kind=CUSTOM_REAL) :: displ_x,displ_z,displ_n,&
                            xxi,zxi,xgamma,zgamma,jacobian1D,nx,nz,weight
  ! PML
  integer :: ispec_PML,CPML_region_local,singularity_type_xz
  double precision :: kappa_x,kappa_z,d_x,d_z,alpha_x,alpha_z,beta_x,beta_z, &
                      A8,A9,A10,bb_xz_1,bb_xz_2,coef0_xz_1,coef1_xz_1,coef2_xz_1,coef0_xz_2,coef1_xz_2,coef2_xz_2

  ! loop on all the coupling edges
  do inum = 1,num_fluid_solid_edges

    ! get the edge of the acoustic element
    ispec_acoustic = fluid_solid_acoustic_ispec(inum)
    iedge_acoustic = fluid_solid_acoustic_iedge(inum)

    ! get the corresponding edge of the elastic element
    ispec_elastic = fluid_solid_elastic_ispec(inum)
    iedge_elastic = fluid_solid_elastic_iedge(inum)

    ! Check if current element matches inum 
    if(ispec_acoustic /= ispec_current) cycle

    ! implement 1D coupling along the edge
    do ipoin1D = 1,NGLLX
    
      ! get point values for the elastic side, which matches our side in the inverse direction
      i = ivalue_inverse(ipoin1D,iedge_elastic)
      j = jvalue_inverse(ipoin1D,iedge_elastic)
      iglob = ibool(i,j,ispec_elastic)
      
      displ_x = displ_elastic(1,iglob)
      displ_z = displ_elastic(2,iglob)

      ! get point values for the acoustic side
      i = ivalue(ipoin1D,iedge_acoustic)
      j = jvalue(ipoin1D,iedge_acoustic)
      iglob = ibool(i,j,ispec_acoustic)
      
      ! Check if current node matches ipoin1D 
      if(i /= i_current .AND. j /= j_current) cycle

      ! compute the 1D Jacobian and the normal to the edge: for their expression see for instance
      ! O. C. Zienkiewicz and R. L. Taylor, The Finite Element Method for Solid and Structural Mechanics,
      ! Sixth Edition, electronic version, www.amazon.com, p. 204 and Figure 7.7(a),
      ! or Y. K. Cheung, S. H. Lo and A. Y. T. Leung, Finite Element Implementation,
      ! Blackwell Science, page 110, equation (4.60).

      if (iedge_acoustic == ITOP) then

        xxi = + gammaz(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
        zxi = - gammax(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
        jacobian1D = sqrt(xxi**2 + zxi**2)
        nx = - zxi / jacobian1D
        nz = + xxi / jacobian1D
        weight = jacobian1D * wxgll(i)

      else if (iedge_acoustic == IBOTTOM) then

        xxi = + gammaz(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
        zxi = - gammax(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
        jacobian1D = sqrt(xxi**2 + zxi**2)
        nx = + zxi / jacobian1D
        nz = - xxi / jacobian1D
        weight = jacobian1D * wxgll(i)

      else if (iedge_acoustic ==ILEFT) then

        xgamma = - xiz(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
        zgamma = + xix(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
        jacobian1D = sqrt(xgamma**2 + zgamma**2)
        nx = - zgamma / jacobian1D
        nz = + xgamma / jacobian1D
        weight = jacobian1D * wzgll(j)

      else if (iedge_acoustic ==IRIGHT) then

        xgamma = - xiz(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
        zgamma = + xix(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
        jacobian1D = sqrt(xgamma**2 + zgamma**2)
        nx = + zgamma / jacobian1D
        nz = - xgamma / jacobian1D
        weight = jacobian1D * wzgll(j)

      endif

      ! compute dot product
      displ_n = displ_x*nx + displ_z*nz
      potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) + weight*displ_n

    enddo
    
  enddo

  end subroutine compute_coupling_acoustic_DG_el

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

! ------------------------------------------------------------ !
! find_normals                                                 !
! ------------------------------------------------------------ !
! Find normal vectors for every edge.

subroutine find_normals()

  use constants,only: CUSTOM_REAL,NGLLX,NGLLZ

  use specfem_par, only : nspec,ibool_before_perio, coord, wzgll, &
        jacobian, gammax, gammaz, xix, xiz, wxgll, wzgll, &
        ispec_is_acoustic, &
        ! MODIF MAJEURE
        nx_iface, nz_iface, weight_iface, link_iface_ijispec, neighbor_DG_iface, &
        AXISYM, is_on_the_axis, wxglj, xiglj

  implicit none
  
  ! Local variables.
  integer :: i, j, numelem
  double precision :: sign_cond1, sign_cond2, sign_cond3, sign_cond4
  double precision :: test_node_x, test_node_z
  double precision :: xxi, zxi, xgamma, zgamma, &
                      jacobian1D_1, jacobian1D_2, jacobian1D, &
                      nx1, nz1, nx2, nz2, nx, nz, weight, coef
  double precision :: cond_1, cond_2, prod1, prod2
  double precision, dimension(4,2) :: pos, vd, mid
  integer, dimension(NGLLX, NGLLZ, nspec) :: ibool
  real(kind=CUSTOM_REAL) :: nx_neighbor, nz_neighbor
  integer :: iface, iglob_first, iglob_last, iface1
  double precision, dimension(2) :: mid_iface, vd_iface
  
  real(kind=CUSTOM_REAL), parameter :: ONE  = 1._CUSTOM_REAL
  real(kind=CUSTOM_REAL), dimension(NGLLX, NGLLZ) :: r_xiplus1
  
  allocate(nx_iface(4,nspec), nz_iface(4,nspec), weight_iface(NGLLX,4,nspec))
  
  ibool = ibool_before_perio ! For normal computation => use old ibool before periodic cond was applied
  
  do numelem = 1,nspec
    if (.not. ispec_is_acoustic(numelem)) cycle ! Skip non acoustic (thus non fluid) elements.
    
    ! Find corner location.
    pos(1,:) = coord(:,ibool(1,1,numelem))
    pos(2,:) = coord(:,ibool(NGLLX,1,numelem))
    pos(3,:) = coord(:,ibool(1,NGLLZ,numelem))
    pos(4,:) = coord(:,ibool(NGLLX,NGLLZ,numelem)) 
    
    ! Direction vector of each faces.
    vd(1,:)  = pos(1,:) - pos(2,:)
    vd(1,:)  = vd(1,:)/sqrt(vd(1,1)**2 + vd(1,2)**2 )
    mid(1,:) = 0.5d0*DBLE(pos(1,:) + pos(2,:))
    vd(2,:) = pos(3,:) - pos(1,:)
    vd(2,:)  = vd(2,:)/sqrt(vd(2,1)**2 + vd(2,2)**2 )
    mid(2,:) = 0.5d0*DBLE(pos(1,:) + pos(3,:))
    vd(3,:) = pos(4,:) - pos(3,:)
    vd(3,:)  = vd(3,:)/sqrt(vd(3,1)**2 + vd(3,2)**2 )
    mid(3,:) = 0.5d0*DBLE(pos(4,:) + pos(3,:))
    vd(4,:) = pos(2,:) - pos(4,:)
    vd(4,:)  = vd(4,:)/sqrt(vd(4,1)**2 + vd(4,2)**2 )
    mid(4,:) = 0.5d0*DBLE(pos(4,:) + pos(2,:))
  
    ! Edge orientation finding.
    do iface = 1, 4
      iglob_first = ibool(link_iface_ijispec(1,iface,numelem,1),link_iface_ijispec(1,iface,numelem,2),numelem)
      iglob_last  = ibool(link_iface_ijispec(NGLLX,iface,numelem,1),link_iface_ijispec(NGLLX,iface,numelem,2),numelem)
      
      mid_iface(:) = 0.5d0*DBLE(coord(:,iglob_first) + coord(:,iglob_last))
      vd_iface(:)  = coord(:,iglob_first) - coord(:,iglob_last)
      vd_iface(:)  = vd_iface(:)/sqrt(vd_iface(1)**2 + vd_iface(2)**2)
      do iface1 = 1, NGLLX
        i = link_iface_ijispec(iface1,iface,numelem,1)
        j = link_iface_ijispec(iface1,iface,numelem,2)
        
        ! Top/Bottom normal directions
        xxi = + gammaz(i,j,numelem) * jacobian(i,j,numelem)
        zxi = - gammax(i,j,numelem) * jacobian(i,j,numelem)
        jacobian1D_1 = sqrt(xxi**2 + zxi**2)
        nx1 = - zxi / jacobian1D_1
        nz1 = + xxi / jacobian1D_1
        
        !prod1 = nx1*vd(ind,1) + nz1*vd(ind,2)
        prod1 = nx1*vd_iface(1) + nz1*vd_iface(2)
        
        ! Right/left normal directions
        xgamma = - xiz(i,j,numelem) * jacobian(i,j,numelem)
        zgamma = + xix(i,j,numelem) * jacobian(i,j,numelem)
        jacobian1D_2 = sqrt(xgamma**2 + zgamma**2)
        nx2 = + zgamma / jacobian1D_2
        nz2 = - xgamma / jacobian1D_2
       
        !prod2 = nx2*vd(ind,1) + nz2*vd(ind,2)
        prod2 = nx2*vd_iface(1) + nz2*vd_iface(2)
        
        if (AXISYM) then
! axial elements are always rotated by the mesher to make sure it is their i == 1 edge that is on the axis
! and thus the case of an edge located along j == constant does not need to be considered here
        if (is_on_the_axis(numelem) .and. i == 1) then
          xxi = + gammaz(i,j,numelem) * jacobian(i,j,numelem)
          r_xiplus1(i,j) = xxi
        else if (is_on_the_axis(numelem)) then
          r_xiplus1(i,j) = coord(1,ibool(i,j,numelem))/(xiglj(i)+ONE)
        endif
        endif
        
        if(abs(prod1) < abs(prod2)) then
          nx = nx1
          nz = nz1
          jacobian1D = jacobian1D_1
          !weight     = jacobian1D*wxgll(i)
          
          ! Axysimmetric implementation
          if (AXISYM) then
          if (is_on_the_axis(numelem)) then
             weight = jacobian1D * wxglj(i) * r_xiplus1(i,j)
             !WRITE(*,*) "weight", weight, nx, nz
          else
             weight = jacobian1D * wxgll(i) * coord(1,ibool(i,j,numelem))
          endif
          else
            weight  = jacobian1D * wxgll(i)
          endif
          ! BIG MODIF INSTEAD OF  
          !nx = -vd_iface(2)
          !nz = vd_iface(1)
        else
          nx = nx2
          nz = nz2
          jacobian1D = jacobian1D_2
          
          if (AXISYM) then
          if (is_on_the_axis(numelem)) then
             weight = jacobian1D * wzgll(j) * r_xiplus1(i,j)
             !if(nx == -1) then
             !   WRITE(*,*) "weight", weight
             !   stop 'TOTO'
             !endif
             !WRITE(*,*) "weight", weight, nx, nz
          else
             weight = jacobian1D * wzgll(j) * coord(1,ibool(i,j,numelem))
          endif
          else
            weight  = jacobian1D * wzgll(j)
          endif
          
          !nx = -vd_iface(2)
          !nz = vd_iface(1)
        endif
        
        ! Remove small values
        coef = 1d-18
        if(abs(nx) < coef) nx = 0d0
        if(abs(nz) < coef) nz = 0d0
        
        !coef = sqrt(vd_iface(1)**2 + vd_iface(2)**2)/10.
        coef = sqrt((coord(1,iglob_first) - coord(1,iglob_last))**2 +&
          (coord(2,iglob_first) - coord(2,iglob_last))**2 )/100.
        test_node_x = mid_iface(1) + coef*nx
        test_node_z = mid_iface(2) + coef*nz
        
        ! Verify if the test node is inside or outside the convec shape
        ! French technique : https://fr.wikipedia.org/wiki/Quadrilat%C3%A8re
        cond_1 = (pos(2,2) - pos(1,2))*test_node_x - (pos(2,1) - pos(1,1))*test_node_z &
                 - pos(1,1)*pos(2,2) + pos(2,1)*pos(1,2)
        cond_2 = (pos(2,2) - pos(1,2))*pos(4,1) - (pos(2,1) - pos(1,1))*pos(4,2) &
                 - pos(1,1)*pos(2,2) + pos(2,1)*pos(1,2)
        sign_cond1 = sign(1.d0, cond_1*cond_2)
           
        cond_1 = (pos(4,2) - pos(2,2))*test_node_x - (pos(4,1) - pos(2,1))*test_node_z &
                 - pos(2,1)*pos(4,2) + pos(4,1)*pos(2,2)
        cond_2 = (pos(4,2) - pos(2,2))*pos(3,1) - (pos(4,1) - pos(2,1))*pos(3,2) &
                 - pos(2,1)*pos(4,2) + pos(4,1)*pos(2,2)
        sign_cond2 = sign(1.d0, cond_1*cond_2)
          
        cond_1 = (pos(3,2) - pos(4,2))*test_node_x - (pos(3,1) - pos(4,1))*test_node_z &
                 - pos(4,1)*pos(3,2) + pos(3,1)*pos(4,2)
        cond_2 = (pos(3,2) - pos(4,2))*pos(1,1) - (pos(3,1) - pos(4,1))*pos(1,2) &
                 - pos(4,1)*pos(3,2) + pos(3,1)*pos(4,2)
        sign_cond3 = sign(1.d0, cond_1*cond_2)
        
        cond_1 = (pos(1,2) - pos(3,2))*test_node_x - (pos(1,1) - pos(3,1))*test_node_z &
                 - pos(3,1)*pos(1,2) + pos(1,1)*pos(3,2)
        cond_2 = (pos(1,2) - pos(3,2))*pos(2,1) - (pos(1,1) - pos(3,1))*pos(2,2) &
                 - pos(3,1)*pos(1,2) + pos(1,1)*pos(3,2)
        sign_cond4 = sign(1.d0, cond_1*cond_2)
           
        ! If point inside the convex shape => take the normal in the opposite direction
        ! to get a normal pointing outside
        if(sign_cond1 + sign_cond2 + sign_cond3 + sign_cond4 >= 4d0) then
          nx = -nx
          nz = -nz
        endif
        
        if (AXISYM) then
          if (is_on_the_axis(numelem) .AND. coord(1,ibool(i,j,numelem)) == 0.) then
                WRITE(*,*) "coord:", weight,coord(:,ibool(i,j,numelem))
                !weight = 0.
          endif
        endif
        
        nx_iface(iface,numelem) = real(nx, kind=CUSTOM_REAL)
        nz_iface(iface,numelem) = real(nz, kind=CUSTOM_REAL)
        weight_iface(iface1,iface,numelem) = real(weight, kind=CUSTOM_REAL)
        
      enddo ! Enddo on iface1.
    enddo ! Enddo on iface.
  enddo ! Enddo on numelem.
  
  do numelem = 1,nspec
    if (.not. ispec_is_acoustic(numelem)) cycle ! Skip non acoustic (thus non fluid) elements.
    do iface = 1,4
      do iface1 = 1,NGLLX
        i = link_iface_ijispec(iface1,iface,numelem,1)
        j = link_iface_ijispec(iface1,iface,numelem,2)
        nx = nx_iface(iface,numelem)
        nz = nz_iface(iface,numelem)
        if(neighbor_DG_iface(iface1,iface,numelem,1) > -1) then
          nx_neighbor = nx_iface( neighbor_DG_iface(iface1,iface,numelem,2), neighbor_DG_iface(iface1,iface,numelem,3) )
          nz_neighbor = nz_iface( neighbor_DG_iface(iface1,iface,numelem,2), neighbor_DG_iface(iface1,iface,numelem,3) )
        endif
      enddo ! Enddo on iface1.
    enddo ! Enddo on iface.
  enddo ! Enddo on numelem.
end subroutine find_normals

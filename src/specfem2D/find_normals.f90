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

  subroutine find_normals()

! find normal at each edge

  use constants,only: CUSTOM_REAL,NGLLX,NGLLZ

  use specfem_par, only : nspec,ibool_before_perio, coord, normal_DG, wzgll, &
        neighbor_DG, weight_DG, jacobian, gammax, gammaz, xix, xiz, wxgll, wzgll, &
        is_corner, ispec_is_acoustic, normal_DG_corner, weight_DG_corner, &
        dir_normal_DG, dir_normal_DG_corner, &
        DIR_RIGHT, DIR_LEFT, DIR_UP, DIR_DOWN, neighbor_DG, neighbor_DG_corner

  implicit none
  
  integer :: i, j, ind, ind_corner, ind_classic, numelem, corner_loop
  !real(kind=CUSTOM_REAL) :: sign_cond1, sign_cond2, sign_cond3, sign_cond4
  !real(kind=CUSTOM_REAL) :: test_node_x, test_node_z
  double precision :: sign_cond1, sign_cond2, sign_cond3, sign_cond4
  double precision :: test_node_x, test_node_z
  !real(kind=CUSTOM_REAL) :: xxi, zxi, xgamma, zgamma, &
  !      jacobian1D_1, jacobian1D_2, jacobian1D, &
  !      nx1, nz1, nx2, nz2, nx, nz, weight, coef
  double precision :: xxi, zxi, xgamma, zgamma, &
        jacobian1D_1, jacobian1D_2, jacobian1D, &
        nx1, nz1, nx2, nz2, nx, nz, weight, coef
  double precision :: cond_1, cond_2, prod1, prod2
  !real(kind=CUSTOM_REAL), dimension(4,2) :: pos, vd
  double precision, dimension(4,2) :: pos, vd, mid
  integer, dimension(NGLLX,NGLLZ,nspec) :: ibool
 
  ! Check
  integer :: inc, inc_tot, ispec, inc_tot_surf
  integer, dimension(3) :: neighbor
  
  logical, parameter :: CHECK_NORMALS = .false.
  
  integer :: dir_normal
  
  !integer :: ii
  !double precision :: prod3, prod4
 
  ! For normal computation => use old ibool before periodic cond was applied
  ibool = ibool_before_perio
 
  WRITE(*,*) "-- compute normals for DG simulations"
 
  dir_normal_DG_corner = 0
  dir_normal_DG        = 0
 
! Temporary way to find orientation
  do numelem = 1,nspec
  
    ! Skip non acoustic (thus non fluid) elements
    if (.not. ispec_is_acoustic(numelem)) cycle
    
    ! Find corner location
    pos(1,:) = coord(:,ibool(1,1,numelem))
    pos(2,:) = coord(:,ibool(NGLLX,1,numelem))
    pos(3,:) = coord(:,ibool(1,NGLLZ,numelem))
    pos(4,:) = coord(:,ibool(NGLLX,NGLLZ,numelem)) 
    
    !do i = 1,NGLLX
    !do j = 1,NGLLZ
    !    WRITE(*,*) "POS >> (",i,j,") : ", coord(:,ibool(i,j,numelem))
    !enddo
    !enddo
    !if(pos(4,1) == 10.) WRITE(*,*) "TOTO", numelem, pos(4,:)
    !stop
    ! Direction vector of each faces
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
  
    !WRITE(*,*) "POS", pos(1,:), pos(2,:), pos(3,:), pos(4,:)
  
    do i = 1,NGLLX
    do j = 1,NGLLZ
    
  ! Temporary EDGE ORIENTATION FINDING
    if(i == 1 .OR. j == 1 .OR. j == NGLLZ .OR. i == NGLLX) then

        ! Find indices to recover element sides and corner orientations
        ind_corner = -1
        if(j == 1) then
                ind_classic = 1
                if(i == 1) ind_corner = 1
        endif
        if(i == 1) then
                ind_classic = 2
                if(j == NGLLZ) ind_corner = 2
        endif
        if(j == NGLLZ) then
                ind_classic = 3
                if(i == NGLLX) ind_corner = 3
        endif
        if(i == NGLLX) then
                ind_classic = 4
                if(j == 1) ind_corner = 1
        endif
        
        ! Loop over corner various normals
        do corner_loop = 1,1+INT(max(ind_corner, 0)/ind_corner)
        
        ind = ind_classic
        if(corner_loop > 1) ind = ind_corner
        
        ! Top/Bottom normal directions
        xxi = + gammaz(i,j,numelem) * jacobian(i,j,numelem)
        zxi = - gammax(i,j,numelem) * jacobian(i,j,numelem)
        jacobian1D_1 = sqrt(xxi**2 + zxi**2)
        nx1 = - zxi / jacobian1D_1
        nz1 = + xxi / jacobian1D_1
        ! Test
        !nx1 = vd(ind,2)
        !nz1 = -vd(ind,1)
        ! Trick to avoid non normal findings
        !!if(abs(vd(ind,1)) > 1d-8) then
        !        nz1 = -nx1*vd(ind,1)/vd(ind,2)
        !elseif(abs(vd(ind,2)) > 1d-8) then
        !        nx1 = -nz1*vd(ind,2)/vd(ind,1)
        !endif
        !nx1 = nx1/sqrt(nx1**2 + nz1**2)
        !nz1 = nz1/sqrt(nx1**2 + nz1**2)
        prod1 = nx1*vd(ind,1) + nz1*vd(ind,2)
        
        ! Right/left normal directions
        xgamma = - xiz(i,j,numelem) * jacobian(i,j,numelem)
        zgamma = + xix(i,j,numelem) * jacobian(i,j,numelem)
        jacobian1D_2 = sqrt(xgamma**2 + zgamma**2)
        nx2 = + zgamma / jacobian1D_2
        nz2 = - xgamma / jacobian1D_2
        ! Trick to avoid non normal findings
        !if(abs(vd(ind,1)) > 1d-8) then
        !        nz2 = -nx2*vd(ind,2)/vd(ind,1)
        !elseif(abs(vd(ind,2)) > 1d-8) then
        !        nx2 = -nz2*vd(ind,1)/vd(ind,2)
        !endif
        !nx2 = nx1/sqrt(nx1**2 + nz1**2)
        !nz2 = nz1/sqrt(nx1**2 + nz1**2)
        prod2 = nx2*vd(ind,1) + nz2*vd(ind,2)
        
        !if(numelem == 2013) &
        !WRITE(*,*) "NORMAL COMPUTATION :",i,j,numelem,nx1, nz1,nx2, nz2, prod1, prod2!, coord(:,ibool(i,j,numelem))
        
        !if(abs(prod1) > 1d-8 .AND. abs(prod2) > 1d-8) then
        if(.false. .AND. numelem == 2015) then
        
        WRITE(*,*) "POS >> (",numelem,i,j,") : ", coord(:,ibool(i,j,numelem)), ibool(i,j,numelem)
        WRITE(*,*) "nx/nz ------ >", nx1, nz1, nx2, nz2
        WRITE(*,*) "vd ------ >", vd(1,:), vd(2,:), vd(3,:), vd(4,:)
        WRITE(*,*) "prod ------ >", prod1, prod2
        WRITE(*,*) "jacobian1D_1 ------ >", jacobian1D_1*wxgll(i), jacobian1D_1*wzgll(j), &
                jacobian1D_2*wxgll(i), jacobian1D_2*wzgll(j)
        WRITE(*,*) "**********"
        endif
        !!do ii = 1,4
        !!        prod3 = nx1*vd(ii,1) + nz1*vd(ii,2)
        !!        prod4 = nx2*vd(ii,1) + nz2*vd(ii,2)
        !!        WRITE(*,*) "prod ------ >", ii, prod3, prod4
        !!enddo
        !WRITE(*,*) "*************"
        !endif
        ! Choose the right normal
        ! => Dot product between direction vector and normal should be zero
        !WRITE(*,*) "TEST PROD", i,j,numelem, prod1, prod2, &
        !        nx1, nz1, nx2, nz2, coord(:,ibool(i,j,numelem))
        !WRITE(*,*) "---------------"
        if(abs(prod1) < abs(prod2)) then
                nx = nx1
                nz = nz1
                jacobian1D = jacobian1D_1
                weight     = jacobian1D*wxgll(i)
                dir_normal = DIR_UP
        else
                nx = nx2
                nz = nz2
                jacobian1D = jacobian1D_2
                weight     = jacobian1D*wzgll(j)
                dir_normal = DIR_RIGHT
        endif
        
        ! Remove small values
        !coef = 0.0000001_CUSTOM_REAL
        !if(abs(nx) < coef) nx = 0._CUSTOM_REAL
        !if(abs(nz) < coef) nz = 0._CUSTOM_REAL
        coef = 1d-18
        if(abs(nx) < coef) nx = 0d0
        if(abs(nz) < coef) nz = 0d0
        
        !if(numelem == 1) then
        !WRITE(*,*) "DEBUG NORMALS 2 :",i, j, numelem, coord(:,ibool(i,j,numelem))
        !WRITE(*,*) "DEBUG NORMALS 2 :", nx1, nz1, nx2, nz2, coord(:,ibool(i,j,numelem)), i, j, vd(ind,:)
        !WRITE(*,*) "DEBUG NORMALS 23 :",nx, nz
        !WRITE(*,*) "DEBUG NORMALS 232 :",mid(ind,1), mid(ind,2)
        !endif
        
        ! 10d-4 is a small positive random value that should be several times smaller 
        ! than the smallest edge of the element in order to give a point inside the element
        ! in the direction of the normal if you start at mid dsitance of the considered edge
        !coef = 0.01_CUSTOM_REAL
        coef = min( sqrt( (pos(1,1) - pos(2,1))**2 + (pos(1,2) - pos(2,2))**2 )/10,&
                 sqrt( (pos(1,2) - pos(3,2))**2 + (pos(1,1) - pos(3,1))**2 )/10 )
        test_node_x = mid(ind,1) + coef*nx
        test_node_z = mid(ind,2) + coef*nz
        !if(numelem == 1) then
        !WRITE(*,*) "DEBUG NORMALS 24 :",test_node_x, test_node_z
        !WRITE(*,*) "DEBUG NORMALS 22 :", pos(1,:), pos(2,:), pos(3,:), pos(4,:)
        !endif
        
        ! Verify if the test node is inside or outside the convec shape
        ! French technique : https://fr.wikipedia.org/wiki/Quadrilat%C3%A8re
        cond_1 = (pos(2,2) - pos(1,2))*test_node_x - (pos(2,1) - pos(1,1))*test_node_z &
                 - pos(1,1)*pos(2,2) + pos(2,1)*pos(1,2)
        cond_2 = (pos(2,2) - pos(1,2))*pos(4,1) - (pos(2,1) - pos(1,1))*pos(4,2) &
                 - pos(1,1)*pos(2,2) + pos(2,1)*pos(1,2)
        sign_cond1 = sign(1.d0, cond_1*cond_2)
           ! if(numelem == 6400) WRITE(*,*) "SIGNCOND1", sign_cond
           
        cond_1 = (pos(4,2) - pos(2,2))*test_node_x - (pos(4,1) - pos(2,1))*test_node_z &
                 - pos(2,1)*pos(4,2) + pos(4,1)*pos(2,2)
        cond_2 = (pos(4,2) - pos(2,2))*pos(3,1) - (pos(4,1) - pos(2,1))*pos(3,2) &
                 - pos(2,1)*pos(4,2) + pos(4,1)*pos(2,2)
        sign_cond2 = sign(1.d0, cond_1*cond_2)
          ! if(numelem == 6400) WRITE(*,*) "SIGNCOND2", sign_cond      
          
        cond_1 = (pos(3,2) - pos(4,2))*test_node_x - (pos(3,1) - pos(4,1))*test_node_z &
                 - pos(4,1)*pos(3,2) + pos(3,1)*pos(4,2)
        cond_2 = (pos(3,2) - pos(4,2))*pos(1,1) - (pos(3,1) - pos(4,1))*pos(1,2) &
                 - pos(4,1)*pos(3,2) + pos(3,1)*pos(4,2)
        sign_cond3 = sign(1.d0, cond_1*cond_2)
        !     if(numelem == 1) WRITE(*,*) "SIGNCOND3", cond_1, cond_2    , &
        !     (pos(3,2) - pos(4,2))*test_node_x,-(pos(3,1) - pos(4,1))*test_node_z, &
        !     -pos(4,1)*pos(3,2), pos(3,1)*pos(4,2)
        
        cond_1 = (pos(1,2) - pos(3,2))*test_node_x - (pos(1,1) - pos(3,1))*test_node_z &
                 - pos(3,1)*pos(1,2) + pos(1,1)*pos(3,2)
        cond_2 = (pos(1,2) - pos(3,2))*pos(2,1) - (pos(1,1) - pos(3,1))*pos(2,2) &
                 - pos(3,1)*pos(1,2) + pos(1,1)*pos(3,2)
        sign_cond4 = sign(1.d0, cond_1*cond_2)
        !   if(numelem == 6400) WRITE(*,*) "SIGNCOND4", cond_1, cond_2 
           
        ! If point inside the convex shape => take the normal in the opposite direction
        ! to get a normal pointing outside
        if(sign_cond1 + sign_cond2 + sign_cond3 + sign_cond4 >= 4d0) then
        !if(cond_1 >= 0._CUSTOM_REAL .AND. cond_2 >= 0._CUSTOM_REAL) then
                nx = -nx
                nz = -nz
                ! Change orientation
                if(dir_normal == DIR_RIGHT) dir_normal = DIR_LEFT
                if(dir_normal == DIR_UP) dir_normal = DIR_DOWN
        endif
        
        !if(numelem == 1) then
        !WRITE(*,*) "DEBUG NORMALS 25 :",i,j,numelem,nx, nz,weight, cond_1, cond_2, &
        !        pos(1,:),pos(2,:),pos(3,:),pos(4,:)
        !WRITE(*,*) "***********"
        !endif
        
        
        !normal_DG(numelem, ind, 1) = real(nx, kind=CUSTOM_REAL)
        !normal_DG(numelem, ind, 2) = real(nz, kind=CUSTOM_REAL)
        if(corner_loop == 1)  then
        
        normal_DG(i,j,numelem, 1) = real(nx, kind=CUSTOM_REAL)
        normal_DG(i,j,numelem, 2) = real(nz, kind=CUSTOM_REAL)
        
        weight_DG(i, j, numelem) = real(weight, kind=CUSTOM_REAL)
        
        dir_normal_DG(i, j, numelem) = dir_normal
        !if(numelem == 2015 .AND. j == 1) then
        !WRITE(*,*) "POS >> (",numelem,i,j,") : ", ibool(i,j,numelem)
        !WRITE(*,*) "weight_DG ------ >", weight_DG(i, j, numelem)
        !WRITE(*,*) "normal_DG ------ >", normal_DG(i, j, numelem,:)
        !WRITE(*,*) "coord ------ >", coord(:,ibool(i,j,numelem))
        !endif
        
        else
        
        normal_DG_corner(i,j,numelem, 1) = real(nx, kind=CUSTOM_REAL)
        normal_DG_corner(i,j,numelem, 2) = real(nz, kind=CUSTOM_REAL)
        
        weight_DG_corner(i, j, numelem) = real(weight, kind=CUSTOM_REAL)
        
        dir_normal_DG_corner(i, j, numelem) = dir_normal
        !if(numelem == 2015 .AND. j == 1) then
        !WRITE(*,*) "normal_DG correr ------ >", normal_DG_corner(i, j, numelem,:), weight
        !endif
        
        endif
        
        !WRITE(*,*) "NORMAL COMPUTATOPN :",i,j,numelem,nx, nz,weight, coord(:,ibool(i,j,numelem))
        
        enddo ! corner loop
       
    endif
    
    enddo
    enddo
    
    !if(numelem == 2015) then
    !    
    !    WRITE(*,*) "nx ------ >", normal_DG(numelem, :, 1)
    !    WRITE(*,*) "nx ------ >", normal_DG(numelem, :, 2)
    !    
    !    WRITE(*,*) "**********"
    !    endif
    
  enddo
  
  !stop
  
  if(CHECK_NORMALS) then
  ! Check normals
  inc = 0
  inc_tot_surf = 0
  inc_tot = 0
  do ispec = 1, nspec
        
        do i = 1,NGLLX
        do j = 1,NGLLZ
        
            if(i == 1 .OR. j == 1 .OR. j == NGLLZ .OR. i == NGLLX) then
        
            inc_tot = inc_tot + 1
        
            neighbor = neighbor_DG(i,j,ispec,:)
        
            if ( .not. is_corner(i,j) .AND. &
                normal_DG(i,j,ispec,1) == &
                -normal_DG(neighbor(1), neighbor(2), neighbor(3),1) .AND. &
                normal_DG(i,j,ispec,2) == &
                -normal_DG(neighbor(1), neighbor(2), neighbor(3),2) ) inc = inc + 1    
            
            if ( .not. is_corner(i,j) .AND. neighbor(3) == -1) inc_tot_surf = inc_tot_surf + 1
                
            if ( is_corner(i,j) ) then
            
                inc_tot = inc_tot + 1
                
                neighbor = neighbor_DG(i,j,ispec,:)
                
                if(neighbor(3) > -1) then
                
                if ( ( normal_DG(i,j,ispec,1) == &
                -normal_DG(neighbor(1), neighbor(2), neighbor(3),1) .AND. &
                     normal_DG(i,j,ispec,2) == &
                -normal_DG(neighbor(1), neighbor(2), neighbor(3),2) ) .OR. &
                 ( normal_DG(i,j,ispec,1) == &
                -normal_DG_corner(neighbor(1), neighbor(2), neighbor(3),1) .AND. &
                     normal_DG(i,j,ispec,2) == &
                -normal_DG_corner(neighbor(1), neighbor(2), neighbor(3),2) ) .OR. &
                !-----------------------------------------
                ( normal_DG_corner(i,j,ispec,1) == &
                -normal_DG(neighbor(1), neighbor(2), neighbor(3),1) .AND. &
                     normal_DG_corner(i,j,ispec,2) == &
                -normal_DG(neighbor(1), neighbor(2), neighbor(3),2) ) .OR. &
                 ( normal_DG_corner(i,j,ispec,1) == &
                -normal_DG_corner(neighbor(1), neighbor(2), neighbor(3),1) .AND. &
                     normal_DG_corner(i,j,ispec,2) == &
                -normal_DG_corner(neighbor(1), neighbor(2), neighbor(3),2) ) ) inc = inc + 1
                
                else
                
                inc_tot_surf = inc_tot_surf + 1
                
                endif
                
                neighbor = neighbor_DG_corner(i,j,ispec,:)
                
                if(neighbor(3) > -1) then
                
                if ( ( normal_DG(i,j,ispec,1) == &
                -normal_DG(neighbor(1), neighbor(2), neighbor(3),1) .AND. &
                     normal_DG(i,j,ispec,2) == &
                -normal_DG(neighbor(1), neighbor(2), neighbor(3),2) ) .OR. &
                 ( normal_DG(i,j,ispec,1) == &
                -normal_DG_corner(neighbor(1), neighbor(2), neighbor(3),1) .AND. &
                     normal_DG(i,j,ispec,2) == &
                -normal_DG_corner(neighbor(1), neighbor(2), neighbor(3),2) ) .OR. &
                !-----------------------------------------
                ( normal_DG_corner(i,j,ispec,1) == &
                -normal_DG(neighbor(1), neighbor(2), neighbor(3),1) .AND. &
                     normal_DG_corner(i,j,ispec,2) == &
                -normal_DG(neighbor(1), neighbor(2), neighbor(3),2) ) .OR. &
                 ( normal_DG_corner(i,j,ispec,1) == &
                -normal_DG_corner(neighbor(1), neighbor(2), neighbor(3),1) .AND. &
                     normal_DG_corner(i,j,ispec,2) == &
                -normal_DG_corner(neighbor(1), neighbor(2), neighbor(3),2) ) ) inc = inc + 1
                
                else
                
                inc_tot_surf = inc_tot_surf + 1
                
                endif
                
                
                
            endif 
            
            !WRITE(*,*) i,j,ispec," : ",inc,"/",inc_tot
            
            endif
            
        enddo
        enddo

  enddo
  
  WRITE(*,*) "CHECK NORMAL : ", inc, "/", inc_tot, " nodes fine"," and ", inc_tot_surf, "surfaces normals"
  
  !stop
  
  endif
  
  !stop
  
  end subroutine find_normals

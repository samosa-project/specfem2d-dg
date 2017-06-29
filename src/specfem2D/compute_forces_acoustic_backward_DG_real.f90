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


  subroutine compute_forces_acoustic_DG_backward_real(rho_DG_main, rhovx_DG_main, rhovz_DG_main, &
   E_DG_main, dot_rho, dot_rhovx, dot_rhovz, dot_E, timelocal)

! compute forces in the acoustic elements in forward simulation and in adjoint simulation in adjoint inversion

  use constants,only: CUSTOM_REAL,NGLLX,NGLLZ,gamma_euler

  use specfem_par, only: nglob_DG,nspec,ispec_is_acoustic, &
                         xix,xiz,gammax,gammaz,jacobian, &
                         hprimewgll_xx, hprimewgll_zz, wxgll, wzgll, &
                         normal_DG, normal_DG_corner, ibool_DG, weight_DG, weight_DG_corner, &
                         neighbor_DG, neighbor_DG_corner, is_corner, ibool, &
                         elastic_tensor, &!ispec_is_acoustic_coupling_el, veloc_elastic,&
                         dir_normal_DG, dir_normal_DG_corner, &
                         DIR_RIGHT, DIR_LEFT, DIR_UP, DIR_DOWN, &
                         myrank, gammaext_DG, muext, etaext, cnu, coord, &
                         hprime_xx, hprime_zz, &
                         windxext, gammaext_DG, gravityext, rhoext, vpext, pext_DG, &
                         it, i_stage
                         
  implicit none

  real(kind=CUSTOM_REAL), dimension(nglob_DG) :: rho_DG_main, rhovx_DG_main, rhovz_DG_main, E_DG_main
  real(kind=CUSTOM_REAL), dimension(nglob_DG) :: b_rho_DG, b_rhovx_DG, b_rhovz_DG, b_E_DG
  real(kind=CUSTOM_REAL), dimension(nglob_DG) :: dot_rho, dot_rhovx, dot_rhovz, dot_E
  
  ! local parameters
  integer :: ispec,i,j,k,iglob, it_corner!, it, i_stage
  integer :: ifirstelem,ilastelem

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: temp_rho_1, temp_rho_2, &
        temp_rhovx_1, temp_rhovx_2, temp_rhovz_1, temp_rhovz_2, &
        temp_E_1, temp_E_2, &
        temp_rho_gravi, temp_rhovx_gravi, temp_rhovz_gravi, temp_E_gravi

  ! Jacobian matrix and determinant
  real(kind=CUSTOM_REAL) :: xixl,xizl,gammaxl,gammazl,jacobianl

  real(kind=CUSTOM_REAL) :: lambda, nx, nz, weight, &
        temp_unknown_M, temp_unknown_P, temp_unknown2_M, temp_unknown2_P, &
        temp_unknown, temp_unknown2, &
        flux_x, flux_z, flux_n, jump, &
        rho_DG_P, veloc_x_DG_P, veloc_z_DG_P, &
        E_DG_P, p_DG_P, rhovx_DG_P, rhovz_DG_P, timelocal, &
        ! TEST
        gamma_P!, &
        !templ_rhovx_gravi, templ_rhovz_gravi, templ_rho_gravi, templ_E_gravi
    
  !real(kind=CUSTOM_REAL) :: veloc_n_M, veloc_n_P
        
  integer :: iglobM, iglobP
  
  integer, dimension(3) :: neighbor
  
  ! Local
  !real(kind=CUSTOM_REAL), dimension(nglob_DG) :: b_veloc_x_DG, b_veloc_z_DG, b_p_DG
  
  ! Viscosity
  !real(kind=CUSTOM_REAL) :: dux_dxi, dux_dgamma, duz_dxi, duz_dgamma
  !real(kind=CUSTOM_REAL) :: dux_dx, dux_dz, duz_dx, duz_dz
  real(kind=CUSTOM_REAL) :: wxl, wzl
  
  ! Parameters
  real(kind=CUSTOM_REAL), parameter :: ZERO = 0._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: ONE  = 1._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: TWO  = 2._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: HALF = 0.5_CUSTOM_REAL
  
  real(kind=CUSTOM_REAL), parameter :: threshold = 0.000000001_CUSTOM_REAL
  
  ! Temporary
  !real(kind=CUSTOM_REAL), dimension(nglob_DG) :: rho0_DG, E0_DG, p0_DG
  logical, parameter :: ONLY_PERTURBATION = .false.
  
  logical :: neighbor_exists
  integer :: dir_normal, i_ex, j_ex, ispec_ex
  
  ! TEmporary test
  integer :: chosen_nxnz_forMPI
  !character(len=100) file_name
  
  logical :: MPI_change_cpt, exact_interface_flux
  
  integer, dimension(nglob_DG) :: MPI_iglob
  
  real(kind=CUSTOM_REAL) :: &!flux_rho, flux_rhovx, flux_rhovz, flux_E, &
        b_rho_l, b_rhovx_l, b_rhovz_l, b_E_l, gammal, &
        veloc_x_l, rho_l, cpl, p_l, E_l, gravityl, &
        drho_dxi, drho_dgamma, dE_dxi, dE_dgamma, &
        b_drho_dzl, b_dE_dzl
  
  ! for MPI transfer
  MPI_iglob = 1
  
  ifirstelem = 1
  ilastelem = nspec
  
  b_rho_DG   = rho_DG_main
  b_rhovx_DG = rhovx_DG_main
  b_rhovz_DG = rhovz_DG_main
  b_E_DG     = E_DG_main
        
  ! Initialization
  dot_rho   = ZERO
  dot_rhovx = ZERO
  dot_rhovz = ZERO
  dot_E     = ZERO
  
  ! add adjoint force source
  call compute_add_sources_acoustic_DG_backward_real(it, &
        dot_rho, dot_rhovx, dot_rhovz, dot_E)
        
  if(myrank == 0) then
  WRITE(*,*) it,i_stage,"MAXVAL ", maxval(b_rho_DG), minval(b_rho_DG)
  WRITE(*,*) it,i_stage,"MAXVAL ", maxval(b_rhovx_DG), minval(b_rhovx_DG)
  WRITE(*,*) it,i_stage,"MAXVAL ", maxval(b_rhovz_DG), minval(b_rhovz_DG)
  WRITE(*,*) it,i_stage,"MAXVAL ", maxval(b_E_DG), minval(b_E_DG)
  WRITE(*,*) "*****************"
  endif
  
! loop over spectral elements
  do ispec = ifirstelem,ilastelem

    ! acoustic spectral element
    if (ispec_is_acoustic(ispec)) then

      ! first double loop over GLL points to compute and store gradients
      do j = 1,NGLLZ
        do i = 1,NGLLX
        
          iglob = ibool_DG(i,j,ispec)
          
          jacobianl = jacobian(i,j,ispec)
        
          xixl = xix(i,j,ispec)
          xizl = xiz(i,j,ispec)
          gammaxl = gammax(i,j,ispec)
          gammazl = gammaz(i,j,ispec)
          
          wzl = wzgll(j)
          wxl = wxgll(i)
          
          b_rhovx_l = b_rhovx_DG(iglob)
          b_rhovz_l = b_rhovz_DG(iglob)
          b_E_l     = b_E_DG(iglob)
          b_rho_l   = b_rho_DG(iglob)
          
          ! Background parameters
          veloc_x_l = windxext(i,j,ispec)
          gammal = gammaext_DG(ibool_DG(i,j,ispec))
          gravityl = gravityext(i,j,ispec)
          rho_l  = rhoext(i,j,ispec)
          cpl    = vpext(i,j,ispec)
          p_l    = pext_DG(i,j,ispec)
          E_l    = p_l/(gammal - 1.) + (rho_l/2.)*(veloc_x_l**2)
          
          !!!!!!!!!!!!!!!!!!!!!!!!
          ! Inviscid stress tensor
          temp_unknown = -rho_l * (veloc_x_l**2) * b_rhovx_l &
                + ((gammal - 1.)*rho_l/2)*( veloc_x_l**2 ) * b_rhovx_l &
                + ((gammal - 1.)*rho_l/2)*( veloc_x_l**3 ) * b_E_l &
                - ((E_l + p_l)/rho_l) * veloc_x_l * b_E_l
          temp_unknown2 = 0. &
                + ((gammal - 1.)*rho_l/2)*( veloc_x_l**2 ) * b_rhovx_l &
                + 0. &!((gammal - 1.)*rho_l/2)*( veloc_x_l**3 ) * b_E_l &
                - 0. !((E_l + p_l)/rho_l) * veloc_x_l * b_E_l
                
          temp_rho_1(i,j) = wzl * jacobianl * (xixl * temp_unknown + xizl * temp_unknown2) 
          temp_rho_2(i,j) = wxl * jacobianl * (gammaxl * temp_unknown + gammazl * temp_unknown2) 
          
          temp_unknown = ( b_rho_l + ((p_l + E_l)/rho_l) * b_E_l ) &
                + 2. * veloc_x_l * b_rhovx_l &
                - (gammal - 1.) * veloc_x_l * b_rhovx_l &
                - (gammal - 1.) * (veloc_x_l**2) * b_E_l
          temp_unknown2 = 0. &!- b_rho_l + ((p_l + E_l)/rho_l) * b_E_l
                - 0. &
                - (gammal - 1.) * veloc_x_l * b_rhovz_l &
                - 0.
          
          temp_rhovx_1(i,j) = wzl * jacobianl * (xixl * temp_unknown + xizl * temp_unknown2) 
          temp_rhovx_2(i,j) = wxl * jacobianl * (gammaxl * temp_unknown + gammazl * temp_unknown2) 
          
          temp_unknown = 0. &
                + veloc_x_l * b_rhovz_l!(1./rho_l)*b_E_l
          temp_unknown2 = (b_rho_l + ((p_l + E_l)/rho_l) * b_E_l ) &
                + veloc_x_l * b_rhovx_l!b_E_l
          
          temp_rhovz_1(i,j) = wzl * jacobianl * (xixl * temp_unknown + xizl * temp_unknown2) 
          temp_rhovz_2(i,j) = wxl * jacobianl * (gammaxl * temp_unknown + gammazl * temp_unknown2) 
          
          temp_unknown = (gammal - 1.) * b_rhovx_l &
                + gammal * veloc_x_l * b_E_l
          temp_unknown2 = (gammal - 1.) * b_rhovz_l &
                - 0.
                
          temp_E_1(i,j) = wzl * jacobianl * (xixl * temp_unknown + xizl * temp_unknown2) 
          temp_E_2(i,j) = wxl * jacobianl * (gammaxl * temp_unknown + gammazl * temp_unknown2) 
          
          !!!!!!!!!!!!!!!!!!!!!
          ! Gravity potentials
          
           xixl = xix(i,j,ispec)
           xizl = xiz(i,j,ispec)
           gammaxl = gammax(i,j,ispec)
           gammazl = gammaz(i,j,ispec)
           
           !!!!!!!!!!!!!!!!!
           !!!! P
           drho_dxi = 0._CUSTOM_REAL
           drho_dgamma = 0._CUSTOM_REAL
           dE_dxi = 0.
           dE_dgamma = 0.
           ! first double loop over GLL points to compute and store gradients
           ! we can merge the two loops because NGLLX == NGLLZ
           
           do k = 1,NGLLX
            
            ! Gradxi rho*v^2
            veloc_x_l = windxext(k,j,ispec)
            gammal = gammaext_DG(ibool_DG(k,j,ispec))
            gravityl = gravityext(k,j,ispec)
            rho_l  = rhoext(k,j,ispec)
            cpl    = vpext(k,j,ispec)
            p_l    = pext_DG(k,j,ispec)
            E_l    = p_l/(gammal - 1.) + (rho_l/2.)*(veloc_x_l**2)
            
            !WRITE(*,*) ">>>>>>>>>>", veloc_x_l, gammal, gravityl, rho_l, cpl, p_l, E_l
            !stop 'TUTU'
            
            dE_dxi      = dE_dxi + (gammal - 1.) * hprime_xx(i,k)
            drho_dxi    = drho_dxi   + ((E_l + p_l)/rho_l) * hprime_xx(i,k)
            
            ! Gradgamma rho*v^2
            veloc_x_l = windxext(i,k,ispec)
            gammal = gammaext_DG(ibool_DG(i,k,ispec))
            gravityl = gravityext(i,k,ispec)
            rho_l  = rhoext(i,k,ispec)
            cpl    = vpext(i,k,ispec)
            p_l    = pext_DG(i,k,ispec)
            E_l    = p_l/(gammal - 1.) + (rho_l/2.)*(veloc_x_l**2)
            
            dE_dgamma   = dE_dgamma   + (gammal - 1.) * hprime_zz(j,k)
            drho_dgamma     = drho_dgamma     + ((E_l + p_l)/rho_l) * hprime_zz(j,k)
            
           enddo
           
           b_drho_dzl = (drho_dxi) * xizl + (drho_dgamma) * gammazl
           b_dE_dzl   = (dE_dxi) * xizl + (dE_dgamma) * gammazl
          
          temp_rho_gravi(i,j)   = gravityl * b_rhovz_l * jacobianl !&ZERO
          temp_rhovx_gravi(i,j) = 0.
          temp_rhovz_gravi(i,j) = ( gravityl * b_E_l + b_drho_dzl * b_E_l )* jacobianl
          temp_E_gravi(i,j)     = b_rhovz_l * b_dE_dzl * jacobianl
          
        enddo
      enddo
      !WRITE(*,*) "dot_rhovx1", maxval(dot_rhovx), minval(dot_rhovx)
!    
! second double-loop over GLL to compute all the terms
!
      do j = 1,NGLLZ
        do i = 1,NGLLX
          iglob = ibool_DG(i,j,ispec)
          ! along x direction and z direction
          ! and assemble the contributions
            do k = 1,NGLLX
            
              dot_rho(iglob) = dot_rho(iglob) + &
                       (temp_rho_1(k,j) * real(hprimewgll_xx(k,i), kind=CUSTOM_REAL) + &
                       temp_rho_2(i,k) * real(hprimewgll_zz(k,j), kind=CUSTOM_REAL))
                       
              dot_rhovx(iglob) = dot_rhovx(iglob) + &
                       (temp_rhovx_1(k,j) * real(hprimewgll_xx(k,i), kind=CUSTOM_REAL) + &
                       temp_rhovx_2(i,k) * real(hprimewgll_zz(k,j), kind=CUSTOM_REAL))
              
              dot_rhovz(iglob) = dot_rhovz(iglob) + &
                       (temp_rhovz_1(k,j) * real(hprimewgll_xx(k,i), kind=CUSTOM_REAL) + &
                       temp_rhovz_2(i,k) * real(hprimewgll_zz(k,j), kind=CUSTOM_REAL))
              
              dot_E(iglob) = dot_E(iglob) + &
                       (temp_E_1(k,j) * real(hprimewgll_xx(k,i), kind=CUSTOM_REAL) + &
                       temp_E_2(i,k) * real(hprimewgll_zz(k,j), kind=CUSTOM_REAL))
                       
            enddo
            
            wzl = real(wzgll(j), kind=CUSTOM_REAL)
            wxl = real(wxgll(i), kind=CUSTOM_REAL)
            
            dot_rho(iglob)   = dot_rho(iglob)   + temp_rho_gravi(i,j) * wxl * wzl
            dot_rhovx(iglob) = dot_rhovx(iglob) + temp_rhovx_gravi(i,j) * wxl * wzl
            dot_rhovz(iglob) = dot_rhovz(iglob) + temp_rhovz_gravi(i,j) * wxl * wzl
            dot_E(iglob)     = dot_E(iglob)     + temp_E_gravi(i,j) * wxl * wzl

        enddo ! second loop over the GLL points
      enddo
      
      !WRITE(*,*) "1 - dot_rhovx2", maxval(dot_rhovx), minval(dot_rhovx)
      
    it_corner = 0
    MPI_change_cpt = .false.
    
    do  i = 1,NGLLX
      
      j = 1
      do while (j <= NGLLZ)
      
        ! We skip interior points
        if(i > 1 .AND. i < NGLLX .AND. j > 1 .AND. j < NGLLZ) then
                j = j + 1
                cycle
        endif
        
        ! Recover neighbor location
        neighbor = neighbor_DG(i,j,ispec,:)
        
        chosen_nxnz_forMPI = -1
      
        ! Reinit boolean to know if neighbor exists
        neighbor_exists = .false.
        
        !nx = normal_DG(ispec, ind, 1)
        !nz = normal_DG(ispec, ind, 2)
        nx     = normal_DG(i,j,ispec, 1)
        nz     = normal_DG(i,j,ispec, 2)
        weight = weight_DG(i, j, ispec)
        dir_normal = dir_normal_DG(i,j,ispec)
        chosen_nxnz_forMPI = 0
        
        ! Needs x2 points at corners to correctly map edges
        ! => 2 step on the same point
        if(it_corner == 1) then
                neighbor = neighbor_DG_corner(i,j,ispec,:)
                nx     = normal_DG_corner(i,j,ispec, 1)
                nz     = normal_DG_corner(i,j,ispec, 2)
                weight = weight_DG_corner(i, j, ispec)
                dir_normal = dir_normal_DG_corner(i,j,ispec)
                chosen_nxnz_forMPI = 1
                !if(j == 1) ind = 1
                !if(i == 1) ind = 2
                !if(j == NGLLZ) ind = 3
                !if(i == NGLLX) ind = 4
                it_corner = 2
        endif
        
        if(neighbor_DG(i,j,ispec,3) > -1 .OR. &
                neighbor_DG_corner(i,j,ispec,3) > -1) neighbor_exists = .true.
        
        ! If not outer boundary check for corresponding neighbor normal
        if(is_corner(i,j)) then
                
                ! If at least one neighbor exists
                if(neighbor_exists) then
                
                i_ex = neighbor(1)
                j_ex = neighbor(2)
                ispec_ex = neighbor(3)
                ! If corner of an outside edge
                if(it_corner == 2 .AND. neighbor(3) == -1) then
                        i_ex = neighbor_DG(i,j,ispec,1)
                        j_ex = neighbor_DG(i,j,ispec,2)
                        ispec_ex = neighbor_DG(i,j,ispec,3)
                elseif(it_corner < 2 .AND. neighbor(3) == -1) then
                        i_ex = neighbor_DG_corner(i,j,ispec,1)
                        j_ex = neighbor_DG_corner(i,j,ispec,2)
                        ispec_ex = neighbor_DG_corner(i,j,ispec,3)
                endif
                
                ! Cross product to verify if the normal corresponds to the normal
                if( dir_normal /= -dir_normal_DG(i_ex,j_ex,ispec_ex) .AND. &
                     dir_normal /= -dir_normal_DG_corner(i_ex,j_ex,ispec_ex) ) then
                        ! Only change normal if inner element
                        if(neighbor(3) > -1 .AND. it_corner < 2) then
                                nx     = normal_DG_corner(i,j,ispec, 1)
                                nz     = normal_DG_corner(i,j,ispec, 2)
                                weight = weight_DG_corner(i, j, ispec)
                                dir_normal = dir_normal_DG_corner(i,j,ispec)
                                ! MODIF for MPI
                                chosen_nxnz_forMPI = 1
                        elseif(neighbor(3) > -1 .AND. it_corner == 2) then
                                nx     = normal_DG(i,j,ispec, 1)
                                nz     = normal_DG(i,j,ispec, 2)
                                weight = weight_DG(i, j, ispec)
                                dir_normal = dir_normal_DG(i,j,ispec)
                                ! MODIF for MPI
                                chosen_nxnz_forMPI = 0
                        endif
                ! If outside element, if the normal corresponds to the one computed here
                ! it means that we should take the other one
                elseif(neighbor(3) == -1) then
                        ! Only change normal if inner element
                        if(it_corner < 2) then
                                nx     = normal_DG_corner(i,j,ispec, 1)
                                nz     = normal_DG_corner(i,j,ispec, 2)
                                weight = weight_DG_corner(i, j, ispec)
                                dir_normal = dir_normal_DG_corner(i,j,ispec)
                                ! MODIF for MPI
                                chosen_nxnz_forMPI = 1
                        elseif(it_corner == 2) then
                                nx     = normal_DG(i,j,ispec, 1)
                                nz     = normal_DG(i,j,ispec, 2)
                                weight = weight_DG(i, j, ispec)
                                dir_normal = dir_normal_DG(i,j,ispec)
                                ! MODIF for MPI
                                chosen_nxnz_forMPI = 0
                        endif
                endif
                
                endif
                
        endif
        
        ! Interior point
        iglobM = ibool_DG(i,j,ispec)
        
        ! If a MPI surface node has been ill referenced and we need to witch between
        ! normal_DG and normal_DG_corner
        if(MPI_change_cpt) then
        
                if(chosen_nxnz_forMPI == 1) then
                                nx     = normal_DG(i,j,ispec, 1)
                                nz     = normal_DG(i,j,ispec, 2)
                                weight = weight_DG(i, j, ispec)
                                dir_normal = dir_normal_DG(i,j,ispec)
                elseif(chosen_nxnz_forMPI == 0) then
                                nx     = normal_DG_corner(i,j,ispec, 1)
                                nz     = normal_DG_corner(i,j,ispec, 2)
                                weight = weight_DG_corner(i, j, ispec)
                                dir_normal = dir_normal_DG_corner(i,j,ispec)
                endif
                
        endif !if(is_MPI_interface_DG(iglobM) .AND. NPROC > 1)
        
        ! If at corner notify that we will need to go again 
        if( is_corner(i,j) .AND. it_corner == 0) it_corner = 1
        
        rho_DG_P     = ZERO
        rhovx_DG_P   = ZERO
        rhovz_DG_P   = ZERO
        E_DG_P       = ZERO
        veloc_x_DG_P = ZERO
        veloc_z_DG_P = ZERO
        p_DG_P       = ZERO
        
        iglobP = 1
        if(neighbor(1) > -1) &
        iglobP = ibool_DG(neighbor(1),neighbor(2),neighbor(3))
        
        call compute_interface_unknowns_backward(i, j, ispec, rho_DG_P, rhovx_DG_P, &
                rhovz_DG_P, E_DG_P, veloc_x_DG_P, veloc_z_DG_P, p_DG_P, gamma_P,&
                neighbor, MPI_change_cpt, &
                exact_interface_flux, &
                b_rho_DG(iglobM), b_E_DG(iglobM), b_rhovx_DG(iglobM), b_rhovz_DG(iglobM), &
                b_rho_DG(iglobP), b_E_DG(iglobP), b_rhovx_DG(iglobP), b_rhovz_DG(iglobP), &
                MPI_iglob, chosen_nxnz_forMPI, dir_normal, nx, nz, weight, timelocal)
        
        ! Approximate local maximum linearized acoustic wave speed
        ! Lax-Friedrich
        lambda = 0.
        jump   = 0.
        !if(.not. exact_interface_flux) then
        !
        !! TEST
        !!gamma_P = gammaext_DG(iglobM)
        !veloc_n_M = sqrt(b_veloc_x_DG(iglobM)**2+b_veloc_z_DG(iglobM)**2)
        !veloc_n_P = sqrt(veloc_x_DG_P**2+veloc_z_DG_P**2)
        lambda = 0.!windxext(i,j,ispec) + sqrt(abs(gammaext_DG(iglobM)*pext_DG(i,j,ispec)/rhoext(i,j,ispec)))
                
       ! endif
        
        b_rho_l   = b_rho_DG(iglobM)
        b_rhovx_l = b_rhovx_DG(iglobM)
        b_rhovz_l = b_rhovz_DG(iglobM)
        b_E_l     = b_E_DG(iglobM)
        
        veloc_x_l = windxext(i,j,ispec)
        gammal = gammaext_DG(ibool_DG(i,j,ispec))
        gravityl = gravityext(i,j,ispec)
        rho_l  = rhoext(i,j,ispec)
        cpl    = vpext(i,j,ispec)
        p_l    = pext_DG(i,j,ispec)
        E_l    = p_l/(gammal - 1.) + (rho_l/2.)*(veloc_x_l**2)
        
        !WRITE(*,*) ">>>", veloc_x_l, gammal, rho_l, cpl, p_l, E_l, gravityl
        !stop 'RRR'

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Mass conservation equation
        temp_unknown_M = -rho_l * (veloc_x_l**2) * b_rhovx_l &
                + ((gammal - 1.)*rho_l/2)*( veloc_x_l**2 ) * b_rhovx_l &
                + ((gammal - 1.)*rho_l/2)*( veloc_x_l**3 ) * b_E_l &
                - ((E_l + p_l)/rho_l) * veloc_x_l * b_E_l
        temp_unknown_P = -rho_l * (veloc_x_l**2) * rhovx_DG_P &
                + ((gammal - 1.)*rho_l/2)*( veloc_x_l**2 ) * rhovx_DG_P &
                + ((gammal - 1.)*rho_l/2)*( veloc_x_l**3 ) * E_DG_P &
                - ((E_l + p_l)/rho_l) * veloc_x_l * E_DG_P
                
        temp_unknown2_M = ((gammal - 1.)*rho_l/2)*( veloc_x_l**2 ) * b_rhovx_l 
        temp_unknown2_P = ((gammal - 1.)*rho_l/2)*( veloc_x_l**2 ) * rhovx_DG_P 
                
                ! compute dot product
                flux_x = temp_unknown_M + temp_unknown_P
                flux_z = temp_unknown2_M + temp_unknown2_P
                flux_n = flux_x*nx + flux_z*nz
                jump   = b_rho_l - rho_DG_P
        
        dot_rho(iglobM) = dot_rho(iglobM) &!- weight*flux_n*HALF
                -weight*(flux_n + lambda*jump)*HALF
        !!!!!!!!!!!!!!!!!!!!!
        ! x-Momentum equation
        
        temp_unknown_M = ( b_rho_l + ((p_l + E_l)/rho_l) * b_E_l ) &
                + 2. * veloc_x_l * b_rhovx_l &
                - (gammal - 1.) * veloc_x_l * b_rhovx_l &
                - (gammal - 1.) * (veloc_x_l**2) * b_E_l
        temp_unknown_P = ( rho_DG_P + ((p_l + E_l)/rho_l) * E_DG_P ) &
                + 2. * veloc_x_l * rhovx_DG_P &
                - (gammal - 1.) * veloc_x_l * rhovx_DG_P &
                - (gammal - 1.) * (veloc_x_l**2) * E_DG_P
                
        temp_unknown2_M = - (gammal - 1.) * veloc_x_l * b_rhovz_l 
        temp_unknown2_P = - (gammal - 1.) * veloc_x_l * rhovz_DG_P 
        
                ! compute dot product
                flux_x = temp_unknown_M + temp_unknown_P
                flux_z = temp_unknown2_M + temp_unknown2_P
                flux_n = flux_x*nx + flux_z*nz
                jump   = b_rhovx_l - rhovx_DG_P
        
        dot_rhovx(iglobM) = dot_rhovx(iglobM) &!- weight*flux_n*HALF
                -weight*(flux_n + lambda*jump)*HALF
                
        !!!!!!!!!!!!!!!!!!!!!
        ! z-Momentum equation
        temp_unknown_M = veloc_x_l * b_rhovz_l
        temp_unknown_P = veloc_x_l * rhovz_DG_P
                
        temp_unknown2_M = (b_rho_l + ((p_l + E_l)/rho_l) * b_E_l ) &
                + veloc_x_l * b_rhovx_l!b_E_l
        temp_unknown2_P = (rho_DG_P + ((p_l + E_l)/rho_l) * E_DG_P ) &
                + veloc_x_l * rhovx_DG_P!b_E_l
        
                ! compute dot product
                flux_x = temp_unknown_M + temp_unknown_P
                flux_z = temp_unknown2_M + temp_unknown2_P
                flux_n = flux_x*nx + flux_z*nz
                jump   = b_rhovz_l - rhovz_DG_P
                
        dot_rhovz(iglobM) = dot_rhovz(iglobM) &!- weight*flux_n*HALF
                -weight*(flux_n + lambda*jump)*HALF
                
        !!!!!!!!!!!!!!!!!
        ! Energy equation
        temp_unknown_M = (gammal - 1.) * b_rhovx_l &
                + gammal * veloc_x_l * b_E_l
        temp_unknown_P = (gammal - 1.) * rhovx_DG_P &
                + gammal * veloc_x_l * E_DG_P
        
        temp_unknown2_M = (gammal - 1.) * b_rhovz_l 
        temp_unknown2_P = (gammal - 1.) * rhovz_DG_P 
        
                ! compute dot product
                flux_x = temp_unknown_M + temp_unknown_P
                flux_z = temp_unknown2_M + temp_unknown2_P
                flux_n = flux_x*nx + flux_z*nz
                jump   = b_E_l - E_DG_P

        dot_E(iglobM) = dot_E(iglobM) &!- weight*flux_n*HALF
                -weight*(flux_n + lambda*jump)*HALF
        ! Increment NGLLZ counter
        j = j + 1
        
        ! If at corner and first step => go again 
        if(it_corner == 1) j = j - 1
        ! Reset corner notification
        if(it_corner == 2) it_corner = 0

      enddo
    enddo
    
    endif ! end of test if acoustic element
    
  enddo
  
  end subroutine compute_forces_acoustic_DG_backward_real
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  ! Compute initial condition
  subroutine initial_condition_DG_backward()

! compute forces in the acoustic elements in forward simulation and in adjoint simulation in adjoint inversion

  use constants,only: CUSTOM_REAL,NGLLX,NGLLZ,PI

  use specfem_par, only: nglob_DG, ibool_DG, &
        rho_DG, rhovx_DG, rhovz_DG, E_DG, nspec, &
        !ispec_is_acoustic, vpext, &
        hprime_xx, hprime_zz, jacobian, xix, xiz, gammax, gammaz, ibool, &
        gammaext_DG, rhoext, windxext, gravityext, pext_DG, vpext, MODEL, myrank!, &!gravityext, coord &
        !myrank, my_neighbours

  implicit none
  
  !real(kind=CUSTOM_REAL), dimension(nglob_DG) :: veloc_x_DG, veloc_z_DG, p_DG
  
  ! Parameters
  real(kind=CUSTOM_REAL), parameter :: ZEROl = 0._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: ONEl  = 1._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: TWOl  = 2._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: SIXl  = 6._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: HALFl = 0.5_CUSTOM_REAL
  
  integer :: ispec, iglob, i, j
  
  real(kind=CUSTOM_REAL) :: rho_DG_P, rhovx_DG_P, rhovz_DG_P, E_DG_P, &
      veloc_x_DG_P, veloc_z_DG_P, p_DG_P
      
  !character(len=100) file_name
  !  
  !write(file_name,"('gravitest_',i3.3,'.txt')") myrank
  !! Open output forcing file
  !open(10,file=file_name,form='formatted')
  
  do ispec = 1,nspec

        ! first double loop over GLL points to compute and store gradients
        do j = 1,NGLLZ
          do i = 1,NGLLX
          
           iglob = ibool_DG(i,j,ispec)
  
           call boundary_condition_DG_backward(i, j, ispec, ZEROl, rho_DG_P, rhovx_DG_P, rhovz_DG_P, E_DG_P, &
                veloc_x_DG_P, veloc_z_DG_P, p_DG_P)
  
           rho_DG(iglob)     = 0.!rho_DG_P
           rhovx_DG(iglob)   = 0.!rhovx_DG_P
           rhovz_DG(iglob)   = 0.!rhovz_DG_P
           E_DG(iglob)       = 0.!E_DG_P
           
          enddo
       enddo
       
  enddo  
  
  
  end subroutine initial_condition_DG_backward
  
  ! Compute initial condition
  subroutine boundary_condition_DG_backward(i, j, ispec, timelocal, rho_DG_P, rhovx_DG_P, rhovz_DG_P, E_DG_P, &
      veloc_x_DG_P, veloc_z_DG_P, p_DG_P)

! compute forces in the acoustic elements in forward simulation and in adjoint simulation in adjoint inversion

  use constants,only: CUSTOM_REAL,PI

  use specfem_par, only: ibool_before_perio, ibool_DG, coord, MODEL, &
        rhoext, windxext, pext_DG, gravityext, gammaext_DG, &
        etaext, muext, coord_interface, kappa_DG, cp, cnu, vpext, Htabext_DG

  implicit none
  
  integer :: i, j, ispec
  
  real(kind=CUSTOM_REAL) :: rho_DG_P, rhovx_DG_P, rhovz_DG_P, E_DG_P, &
      veloc_x_DG_P, veloc_z_DG_P, p_DG_P, timelocal, gamma_euler
  
  ! Parameters
  real(kind=CUSTOM_REAL), parameter :: ZEROl = 0._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: ONEl  = 1._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: TWOl  = 2._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: SIXl  = 6._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: HALFl = 0.5_CUSTOM_REAL
  
  real(kind=CUSTOM_REAL) :: x, z, H, G, A
  !real(kind=CUSTOM_REAL) :: x0, z0, r, beta
  ! Forcing
  real(kind=CUSTOM_REAL) :: to, perio, lambdo!, xo
  
  real(kind=CUSTOM_REAL) :: mu_visco, eta_visco
  !real(kind=CUSTOM_REAL) :: Tl, Tu, rho0!, p0, RR
  !real(kind=CUSTOM_REAL) :: x0, z0, r, theta, RR, p0, exner, theta0
  !real(kind=CUSTOM_REAL) :: &
  !      lz, vs, vs_x, vs_z, Ms, p_1, rho_1, veloc_x_1, veloc_z_1, p_2, rho_2, veloc_x_2, veloc_z_2 
  
  coord_interface = 0.
  
  x = real(coord(1,ibool_before_perio(i,j,ispec)), kind=CUSTOM_REAL)
  z = real(coord(2,ibool_before_perio(i,j,ispec)), kind=CUSTOM_REAL)
  z = z - coord_interface
  
  ! If an external data file is given for initial conditions
  if(trim(MODEL) == 'external') then
  
  rho_DG_P     = rhoext(i,j,ispec)
  p_DG_P       = pext_DG(i,j,ispec)!rho_DG_P*gravityext(i,j,ispec)*Htabext_DG(ibool_DG(i,j,ispec))!pext_DG(i,j,ispec)
  veloc_x_DG_P = windxext(i,j,ispec)
  veloc_z_DG_P = 0.!windzext(i,j,ispec)
  
    H = 15000!30964!30964.2932
   G = 9.81
   rho_DG_P = ONEl*exp(-z/H)
   p_DG_P   = rho_DG_P*G*H
   
   if(timelocal == 0) then
   gravityext(i,j,ispec) = G!9.831
   gammaext_DG(ibool_DG(i,j,ispec))   = gamma_euler!(rho_DG_P/p_DG_P)*(vpext(i,j,ispec)**2)!gamma_euler!
   muext(i,j,ispec)  = 2e-05
   etaext(i,j,ispec)  = 0!(4/3)*muext(i,j,ispec)
   kappa_DG(i,j,ispec) = 0.2
   cp = 7/2
   cnu = 5/2
   endif
   
   ! Forcing
   lambdo = 100
   perio  = 8*PI*lambdo/600000
   A  = 0.001
   if(z == ZEROl .AND. .true.) then     
        veloc_z_DG_P = A*cos(perio*(timelocal - x/lambdo)) 
        ! Amplitude for smooth starting of forcing
        A = 1e-08
        ! Time at which we have smooth starting of forcing
        lambdo = 10
        ! Time at which we have full forcing
        to = 35
        perio = sqrt(-log(A)/(lambdo-to)**2)
        if(timelocal <= to) &
        veloc_z_DG_P = veloc_z_DG_P*exp(-(perio*(timelocal - to))**2)
   endif
  
  ! If not an external model => means that no initial conditions specified
  else
  
  ! At initial time if no external model we create an uniform gravity fields
  if(timelocal == 0) then
        H = 10000
        rhoext(i,j,ispec) = 1!exp(-z/H)
        gravityext(i,j,ispec) = 0!9.81!81
        gammaext_DG(ibool_DG(i,j,ispec))   = 1.4!gamma_euler
        windxext(i,j,ispec) = 100
        mu_visco  = 0!.0001!2!1e-05!0!0!1e-05!0!
        eta_visco = mu_visco*(4/3)!4!(4/3)*mu_visco!1!4e-05!1!4e-05!
        muext(i,j,ispec)  = 0!mu_visco
        etaext(i,j,ispec) = eta_visco
        kappa_DG(i,j,ispec) = 0
        vpext(i,j,ispec) = 340
        
        pext_DG(i,j,ispec) = (vpext(i,j,ispec)**2)*rhoext(i,j,ispec)/gammaext_DG(ibool_DG(i,j,ispec))
        if(gravityext(i,j,ispec) > 0) pext_DG(i,j,ispec) = &
                rhoext(i,j,ispec)*gravityext(i,j,ispec)*H
  endif
  
  !WRITE(*,*) i,j,ispec,"TEST>>>", rhoext(i,j,ispec), gravityext(i,j,ispec), gammaext_DG(ibool_DG(i,j,ispec)), &
  !      windxext(i,j,ispec), vpext(i,j,ispec), pext_DG(i,j,ispec)
  
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! BASIC FORCING
   H = 10000!30000!30964.2932!15000!30964.2932!15000!30964.2932!12000!
   G = gravityext(i,j,ispec)!9.831!9.81!
   ! http://catalog.conveyorspneumatic.com/Asset/FLS%20Specific%20Heat%20Capacities%20of%20Gases.pdf
   cp = 7!1010!7/2
   cnu = 5!718!5/2
   rho_DG_P     = 1.!exp(-z/H)
   
   veloc_x_DG_P = ZEROl!20
   veloc_z_DG_P = ZEROl
   
   gamma_euler = 1.4
   p_DG_P = (340**2)*rho_DG_P/gamma_euler
   
   to = 3.
   perio = 3.     
   if(z == ZEROl .AND. .false.) &     
   veloc_z_DG_P = 1
  
   endif
        
   E_DG_P       = p_DG_P/(gammaext_DG(ibool_DG(i,j,ispec)) - ONEl) &
                        + rho_DG_P*HALFl*( veloc_x_DG_P**2 + veloc_z_DG_P**2 ) 
   rhovx_DG_P   = rho_DG_P*veloc_x_DG_P
   rhovz_DG_P   = rho_DG_P*veloc_z_DG_P
  
  end subroutine boundary_condition_DG_backward
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
  subroutine compute_interface_unknowns_backward(i, j, ispec, rho_DG_P, rhovx_DG_P, &
                rhovz_DG_P, E_DG_P, veloc_x_DG_P, veloc_z_DG_P, p_DG_P, gamma_P, &
                neighbor, MPI_change_cpt, &
                exact_interface_flux, &
                rho_DG_iM, E_DG_iM, rhovx_DG_iM, rhovz_DG_iM, &
                rho_DG_iP, E_DG_iP, rhovx_DG_iP, rhovz_DG_iP, &
                MPI_iglob, chosen_nxnz_forMPI, dir_normal, nx, nz, weight, timelocal)
  
! compute forces in the acoustic elements in forward simulation and in adjoint simulation in adjoint inversion
  
  use constants,only: CUSTOM_REAL,NGLLX,NGLLZ,gamma_euler

  use specfem_par,only: nglob_DG,nspec, &!ispec_is_acoustic, &
                         normal_DG, normal_DG_corner, ibool_DG, weight_DG, weight_DG_corner, &
                         ispec_is_acoustic_forcing, &
                         ACOUSTIC_FORCING, is_corner, &
                          ispec_is_acoustic_coupling_el, veloc_elastic,&
                         !coord,ibool_before_perio, &
                         dir_normal_DG, dir_normal_DG_corner, &
                         DIR_RIGHT, DIR_LEFT, DIR_UP, DIR_DOWN, &
                         !myrank, &
                         buffer_DG_rho_P, buffer_DG_rhovx_P, buffer_DG_rhovz_P, buffer_DG_E_P, NPROC, &
                         buffer_DG_Vxx_P, buffer_DG_Vzz_P, buffer_DG_Vxz_P, buffer_DG_Vzx_P, buffer_DG_Tz_P, buffer_DG_Tx_P, &
                         MPI_transfer, p_DG_init, gammaext_DG, muext, etaext, kappa_DG, ibool, cnu, &
                         ! TEST
                         buffer_DG_gamma_P
                         
  implicit none
  
  !real(kind=CUSTOM_REAL), dimension(NGLLX, NGLLZ, nspec, 4) :: elastic_tensor
  
  integer, intent(in) :: i, j, ispec, chosen_nxnz_forMPI
  
  integer, dimension(3), intent(in) :: neighbor
  
  real(kind=CUSTOM_REAL), intent(out) :: rho_DG_P, rhovx_DG_P, rhovz_DG_P, &
        E_DG_P, veloc_x_DG_P, veloc_z_DG_P, p_DG_P
        
  real(kind=CUSTOM_REAL), intent(in) :: timelocal
        
  logical, intent(inout) :: exact_interface_flux, MPI_change_cpt
  integer, dimension(nglob_DG), intent(inout) :: MPI_iglob
  integer, intent(inout) :: dir_normal
  
  real(kind=CUSTOM_REAL) :: tx, tz, normal_v, tangential_v, &
        nx, nz, veloc_x, veloc_z, weight, &
        veloc_x_DG_iM, veloc_z_DG_iM, gamma_P
        
  real(kind=CUSTOM_REAL), intent(in) :: rho_DG_iM, E_DG_iM, rhovx_DG_iM, rhovz_DG_iM, &
                rho_DG_iP, E_DG_iP, rhovx_DG_iP, rhovz_DG_iP
        
  ! Local variables     
  integer :: iglobM, i_el, j_el, ispec_el, iglob, iglobP, ipoin, num_interface
  real(kind=CUSTOM_REAL), dimension(2,2) :: trans_boundary
  
  ! Parameters
  real(kind=CUSTOM_REAL), parameter :: ZERO = 0._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: ONE  = 1._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: TWO  = 2._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: HALF = 0.5_CUSTOM_REAL
  
  iglobM = ibool_DG(i,j,ispec)
  
  gamma_P = gammaext_DG(iglobM)
  
  exact_interface_flux = .false.
  MPI_change_cpt = .false.
  if(neighbor(3) == -1) then
        
        ! MPI neighbor
        ipoin         = -1
        num_interface = -1
        if(NPROC > 1) then
        ipoin         = MPI_transfer(iglobM,MPI_iglob(iglobM),1)
        num_interface = MPI_transfer(iglobM,MPI_iglob(iglobM),2)
        endif                  
        ! If MPI neighbor, but not diagonal corner element and not outside element
        if(ipoin > -1) then
                        
                        ! Check for mistakes at corners
                        if(is_corner(i,j) .AND. &
                                ( (MPI_transfer(iglobM,MPI_iglob(iglobM),3) == i .AND. &
                                        (dir_normal == DIR_LEFT .OR. dir_normal == DIR_RIGHT) ) .OR. &
                                  (MPI_transfer(iglobM,MPI_iglob(iglobM),4) == j .AND. &
                                        (dir_normal == DIR_UP .OR. dir_normal == DIR_DOWN) ) ) ) then
                        
                                MPI_change_cpt = .true.
                        
                                if(chosen_nxnz_forMPI == 1) then
                                nx     = normal_DG(i,j,ispec, 1)
                                nz     = normal_DG(i,j,ispec, 2)
                                weight = weight_DG(i, j, ispec)
                                dir_normal = dir_normal_DG(i,j,ispec)
                                elseif(chosen_nxnz_forMPI == 0) then
                                nx     = normal_DG_corner(i,j,ispec, 1)
                                nz     = normal_DG_corner(i,j,ispec, 2)
                                weight = weight_DG_corner(i, j, ispec)
                                dir_normal = dir_normal_DG_corner(i,j,ispec)
                                endif
                                
                        endif
                        
                        MPI_iglob(iglobM) = MPI_iglob(iglobM) + 1
                        
                        rho_DG_P     = buffer_DG_rho_P(ipoin,num_interface)
                        E_DG_P       = buffer_DG_E_P(ipoin,num_interface)
                        rhovx_DG_P   = buffer_DG_rhovx_P(ipoin,num_interface)
                        rhovz_DG_P   = buffer_DG_rhovz_P(ipoin,num_interface)
                        
                        gamma_P   = buffer_DG_gamma_P(ipoin,num_interface)
                              
        elseif(ACOUSTIC_FORCING .AND. ispec_is_acoustic_forcing(i,j,ispec)) then
        
                stop 'ACOUSTIC_FORCING obsolete for DG simulations'
                
        ! Elastic coupling        
        elseif(ispec_is_acoustic_coupling_el(i,j,ispec,3) >= 0) then
        
               ! If we already know the "real" flux at boundary
               !exact_interface_flux = .false.
        
               ! Coordinates of elastic element
               i_el     = ispec_is_acoustic_coupling_el(i,j,ispec,1)
               j_el     = ispec_is_acoustic_coupling_el(i,j,ispec,2)
               ispec_el = ispec_is_acoustic_coupling_el(i,j,ispec,3)
               
               iglob = ibool(i_el,j_el,ispec_el)
        
               ! Only for density
               call boundary_condition_DG_backward(i, j, ispec, timelocal, rho_DG_P, rhovx_DG_P, rhovz_DG_P, E_DG_P, &
                veloc_x_DG_P, veloc_z_DG_P, p_DG_P)
        
               rho_DG_P = rho_DG_iM
               
               ! Elastic velocities
               veloc_x = veloc_elastic(1,iglob)
               veloc_z = veloc_elastic(2,iglob)
               
               ! Tangential vector
               ! Since only bottom topography nz > 0
               tx = -nz
               tz = nx
               
               ! Normal velocity of the solid perturbation and tangential velocity of the fluid flow
               normal_v     = veloc_x*nx + veloc_z*nz
               tangential_v = veloc_x_DG_P*tx + veloc_z_DG_P*tz
               
               E_DG_P       = E_DG_iM
               rhovx_DG_P   = rhovx_DG_iM
               rhovz_DG_P   = rhovz_DG_iM
               
        ! Classical boundary conditions
        else
        
                !exact_interface_flux = .true.
                
                call boundary_condition_DG_backward(i, j, ispec, timelocal, rho_DG_P, rhovx_DG_P, rhovz_DG_P, E_DG_P, &
                veloc_x_DG_P, veloc_z_DG_P, p_DG_P)
                
                rho_DG_P     = -rho_DG_iM
                
                tx = -nz
                tz = nx
                
               ! Normal velocity of the solid perturbation and tangential velocity of the fluid flow
                normal_v     = veloc_x_DG_P*nx + veloc_z_DG_P*nz
                normal_v     = 2*normal_v-(veloc_x_DG_iM*nx + veloc_z_DG_iM*nz)
                tangential_v = veloc_x_DG_P*tx + veloc_z_DG_P*tz
                
               ! Transformation matrix between mesh coordinates and normal/tangential coordinates
                trans_boundary(1,1) = tz
                trans_boundary(1,2) = -nz
                trans_boundary(2,1) = -tx
                trans_boundary(2,2) = nx
                trans_boundary = trans_boundary/(nx*tz - tx*nz)
               
               ! From free slip and normal velocity continuity
               veloc_x_DG_P = trans_boundary(1,1)*normal_v + trans_boundary(1,2)*tangential_v!veloc_elastic(1,iglob)
               veloc_z_DG_P = trans_boundary(2,1)*normal_v + trans_boundary(2,2)*tangential_v
                
                E_DG_P       = 0!E_DG_iM
                
                rhovx_DG_P   = -rhovx_DG_iM
                rhovz_DG_P   = -rhovz_DG_iM
                
        endif
        
   ! Not an outside edge
   else
                iglobP       = ibool_DG(neighbor(1),neighbor(2),neighbor(3))
                
                gamma_P = gammaext_DG(iglobP)
                
                rho_DG_P     = rho_DG_iP
                E_DG_P       = E_DG_iP
                rhovx_DG_P   = rhovx_DG_iP
                rhovz_DG_P   = rhovz_DG_iP
                
  endif
  
  end subroutine compute_interface_unknowns_backward

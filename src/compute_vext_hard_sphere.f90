!> Compute external potential in a basic way
!> Solvation of a hard sphere in a hard sphere fluid
SUBROUTINE compute_Vext_hard_sphere
    use module_grid, only: grid
! use precision_kinds,only : i2b , dp
! use module_input,only : input_line,getinput%dp
! use system,only : solute, grid
! use external_potential,only : Vext_total , Vext_hard
! use quadrature, only: angGrid, molRotGrid
! use hardspheres ,ONLY: hs
!
! IMPLICIT NONE
! integer(i2b):: i , j , k !> @var dummy
! integer(i2b):: species ! dummy between 1 and solvent(1)%nspec
! integer(i2b):: solute_site ! dummy between 1 and
! real(dp):: x_grid , y_grid , z_grid !> @var coordinates of grid mesh nodes
! real(dp):: x_nm2 , y_nm2 , z_nm2 , r_nm2  ! distance between solute and grid point
! real(dp):: hard_sphere_solute_radius ! radius of the hard sphere solute. tag in dft.in
! real(dp):: sum_of_solute_and_solvent_radius ! sum of solute and solvent radius
! INTEGER(i2b) :: nfft1, nfft2, nfft3
! REAL(dp) :: lx, ly, lz, deltax, deltay, deltaz
!
!     lx= grid%length(1)
!     ly= grid%length(2)
!     lz= grid%length(3)
!     deltax= grid%dl(1)
!     deltay= grid%dl(2)
!     deltaz= grid%dl(3)
!     nfft1= grid%n_nodes(1)
!     nfft2= grid%n_nodes(2)
!     nfft3= grid%n_nodes(3)
!
!
! ! tell user
! write (*,*) '*************************'
! write (*,*) '>>> Compute HS Vext for hard sphere in hs fluid'
! write (*,*) 'there are as many hard sphere solute as solutes in solute.in'
! ! be sure Vext_total is allocated
! if ( .not. allocated ( Vext_total ) ) then
!   allocate ( Vext_total ( nfft1 , nfft2 , nfft3 , angGrid%n_angles , molRotGrid%n_angles , solvent(1)%nspec ) )
!   Vext_total = 0.0_dp
! END IF
! ! look for tag 'hard_sphere_solute_radius' which gives the radius of the sphere for the external potential
! hard_sphere_solute_radius=getinput%dp('radius_of_hard_sphere_solute')
! write (*,*) 'the radius of these solute hard spheres is defined as a constant for now and is ', hard_sphere_solute_radius
!
! !check if radius is allocated. if not it may be because you wand an hard sphere solute but no hardsphere solvent, in that case
! !radius=0, if you want an hard sphere solvent there is a bug
! if ( .not. allocated ( hs ) ) then
!     print*,'!!!!!!!!!!!!WARNING!!!!!!!!!!!!!!!!!!WARNING!!!!!!!!!!!!!!!!!!WARNING!!!!!!!!!!!!!!!!WARNING!!!!!!!!!!!!!!WARNING!!'
!     print*,'WARNING solvent radius is not defined in comput_vext_hard_sphere it is now equal to 0, which is ok ONLY if you want '
!     print*,'an HS solute but no HS Solvent'
!     print*,'!!!!!!!!!!!!WARNING!!!!!!!!!!!!!!!!!!WARNING!!!!!!!!!!!!!!!!!!WARNING!!!!!!!!!!!!!!!!WARNING!!!!!!!!!!!!!!WARNING!!'
! END IF
! ! compute the sum of solute and solvent radius
! do solute_site = 1 , solute%nsite
!   do species = 1 , solvent(1)%nspec
!     if ( .not. allocated ( hs ) ) then
!         sum_of_solute_and_solvent_radius = ( hard_sphere_solute_radius)**2
!     else
!         sum_of_solute_and_solvent_radius = ( hard_sphere_solute_radius + HS(species)%R ) ** 2
!     end if
!     do k = 1 , nfft3
!       z_grid = real( k - 1 , dp ) * deltaz
!       z_nm2 = ( z_grid - solute%site(solute_site)%r(3) ) ** 2
!       do j = 1 , nfft2
!         y_grid = real ( j - 1 , dp ) * deltay
!         y_nm2 = ( y_grid - solute%site(solute_site)%r(2) ) ** 2
!         do i = 1 , nfft1
!           x_grid = real ( i - 1 , dp ) * deltax
!           x_nm2 = ( x_grid - solute%site(solute_site)%r(1) ) ** 2
!           r_nm2 = x_nm2 + y_nm2 + z_nm2
!           ! test if distance is smaller than sum of radius of solvent and solute
!           if ( r_nm2 <= sum_of_solute_and_solvent_radius ) then
!             ! we're inside the hard potential. vext is very high
!             Vext_total ( i , j , k , : , : , species ) = huge ( 1.0_dp ) ! no dependancy over Omega here
!             ! if not , Vext_total is unchanged
!           END IF
!         END DO
!       END DO
!     END DO
!   END DO
! END DO
END SUBROUTINE compute_Vext_hard_sphere

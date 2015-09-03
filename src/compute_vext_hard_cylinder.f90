! compute the external potential as created by a hard cylinder solute
! hard cylinder is along axe z=0, x=Lx/2, y=Ly/2
SUBROUTINE compute_Vext_hard_cylinder
    use module_grid, only: grid
! use precision_kinds,only : i2b , dp
! use module_input,only : getinput%dp
! use system,only : solute
! use external_potential,only : Vext_total
! use quadrature, only: angGrid, molRotGrid
! use hardspheres ,ONLY: hs
! IMPLICIT NONE
! integer(i2b):: i , j , n ! dummy
! real(dp):: x_nm2 , y_nm2 , r_nm2 ! distance**2 between solute and grid point
! real(dp):: hard_cylinder_radius ! radius of the hard cylinder solute
! real(dp):: deltax , deltay ! == Lx / nfft1 , Ly / nfft2 , Lz / nfft3
! real(dp):: sum_rad2 ! sum of solute cylinder radius and radius of solvent hard sphere
! integer(i2b):: species ! dummy between 1 and solvent(1)%nspec
!     INTEGER(i2b) :: nfft1, nfft2, nfft3
!     REAL(dp) :: lx, ly, lz
!
!     lx= grid%length(1)
!     ly= grid%length(2)
!     lz= grid%length(3)
!     nfft1= grid%n_nodes(1)
!     nfft2= grid%n_nodes(2)
!     nfft3= grid%n_nodes(3)
! ! init variables
! deltax = Lx / real ( nfft1 , dp )
! deltay = Ly / real ( nfft2 , dp )
! ! be sure Vext_total is allocated
! if ( .not. allocated ( Vext_total ) ) then
!   allocate ( Vext_total ( nfft1 , nfft2 , nfft3 , angGrid%n_angles , molRotGrid%n_angles, solvent(1)%nspec ) )
!   Vext_total = 0.0_dp
! END IF
! ! tell user
! write(*,*)'>>> Compute HS Vext Daniel for cylinder in z dir in hs fluid'
! ! look for tag 'hard_cylinder_radius' in dft.in
! hard_cylinder_radius=getinput%dp('radius_of_hard_cylinder')
! write ( * , * ) 'the hard cylinder solute has a radius of ' , hard_cylinder_radius
! write ( * , * ) 'this cylinder is along Z. position in x and y are extracted from input/solute.in solute%site%r'
! write ( * , * ) 'there are ' , solute%nsite , ' cylinders in the system. They all have same radius'
! ! compute if grid point is or not inside the hard cylinder
! ! hard cylinder is along axe x=Lx/2, y=Ly/2
! do n = 1 , solute%nsite
!   do species = 1 , solvent(1)%nspec
!     sum_rad2 = ( hard_cylinder_radius + hs(species)%R ) ** 2
!     do j = 1 , nfft2
!       y_nm2 = ( real(j-1,dp) * deltay - solute%site(n)%r(2) )**2
!       do i = 1 , nfft1
!         x_nm2 = ( real(i-1,dp) * deltax - solute%site(n)%r(1) )**2
!         ! relative distances between solute site and solvent site
!         r_nm2 = x_nm2 + y_nm2
!         ! if we are inside , vext is huge, ELSE vext is unchanged
!         if ( r_nm2 <= sum_rad2 ) then
!           Vext_total ( i , j , : , : , : , species ) = huge ( 1.0_dp )
!         END IF
!       END DO
!     END DO
!   END DO
! END DO
END SUBROUTINE compute_Vext_hard_cylinder

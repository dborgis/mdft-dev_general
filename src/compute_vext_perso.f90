! this is the routine where you may want to build a personnal vext, something that is not implemented automaticaly
SUBROUTINE compute_vext_perso
    use module_grid, only: grid
!
! use precision_kinds,only : i2b , dp
! use module_input,only : input_line,getinput%dp
! use module_solute, only : solute
! use external_potential,only : Vext_total
! use hardspheres ,ONLY: hs
!
! IMPLICIT NONE
!
!     integer(i2b):: i , j , k , n !> @var dummy
!     real(dp):: x_grid , y_grid , z_grid !> @var coordinates of grid mesh nodes
!     real(dp):: x_nm2 , y_nm2 , z_nm2 , r_nm2  ! distance between solute and grid point
!     real(dp):: solute_radius ! radius of the hard sphere solute. tag in dft.in
!     real(dp):: radius_sum_sq ! sum of solute and solvent radius
!     real(dp):: deltax , deltay , deltaz ! == Lx / nfft1 , Ly / nfft2 , Lz / nfft3
!     integer(i2b):: species ! dummy for loops over species from 1 to solvent(1)%nspec
!     INTEGER(i2b) :: nfft1, nfft2, nfft3
!     REAL(dp) :: lx, ly, lz
!
!     lx= grid%length(1)
!     ly= grid%length(2)
!     lz= grid%length(3)
!     nfft1= grid%n_nodes(1)
!     nfft2= grid%n_nodes(2)
!     nfft3= grid%n_nodes(3)
!
! ! tell user
! write (*,*) '*************************'
! write (*,*) 'Compute personnal_vext'
! ! look for tag 'hard_sphere_solute_radius' which gives the radius of the sphere for the external potential
! solute_radius=getinput%dp('hard_sphere_solute_radius')
! write (*,*) 'the radius of these solute hard spheres is defined as a constant for now and is ', solute_radius
! ! compute vext hard sphere
! deltax = Lx / real ( nfft1 , dp )
! deltay = Ly / real ( nfft2 , dp )
! deltaz = Lz / real ( nfft3 , dp )
! ! check if radius is allocated. if not there is a bug.
! if ( .not. allocated ( hs ) ) then
!   write (*,*) 'hs is not allocated. should be. critial stop in compute_vext_perso.f90'
!   stop
! END IF
! ! compute the sum of solute and solvent radius
! do species = 1 , solvent(1)%nspec
! radius_sum_sq = ( solute_radius + hs(species)%R ) ** 2
! do n = 1 , solute%nsite
!   do i = 1 , nfft1
!     x_grid = real ( i - 1 , dp) * deltax
!     x_nm2 = ( x_grid - solute%site(n)%r(1) )**2
!     do j = 1 , nfft2
!       y_grid = real ( j - 1 , dp ) * deltay
!       y_nm2 = ( y_grid - solute%site(n)%r(2) ) ** 2
!       do k = 1 , nfft3
!         z_grid = real ( k - 1 , dp ) * deltaz
!         z_nm2 = ( z_grid - solute%site(n)%r(3) ) ** 2
!         r_nm2 = x_nm2 + y_nm2 + z_nm2
!         ! test if distance is smaller than sum of radius of solvent and solute
!         if ( r_nm2 <= radius_sum_sq ) then
!           ! we're inside the hard potential. vext is very high
!           Vext_total ( i , j , k , : , : , species ) = huge ( 1.0_dp )
!           ! if not, vext is unchanged
!         END IF
!       END DO
!     END DO
!   END DO
! END DO ! solutes sites
! END DO ! species
END SUBROUTINE compute_vext_perso

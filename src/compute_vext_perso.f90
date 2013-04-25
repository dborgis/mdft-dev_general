! this is the routine where you may want to build a personnal vext, something that is not implemented automaticaly

subroutine compute_vext_perso

use precision_kinds , only : i2b , dp

use input , only : input_line,input_dp

use system , only : nfft1 , nfft2 , nfft3 , x_mol , y_mol , z_mol , nb_solute_sites , Lx , Ly , Lz , radius , nb_species

use external_potential , only : Vext_total


implicit none

integer(i2b):: i , j , k , n !> @var dummy

real(dp):: x_grid , y_grid , z_grid !> @var coordinates of grid mesh nodes

real(dp):: x_nm2 , y_nm2 , z_nm2 , r_nm2  ! distance between solute and grid point

real(dp):: solute_radius ! radius of the hard sphere solute. tag in dft.in

real(dp):: radius_sum_sq ! sum of solute and solvent radius

real(dp):: deltax , deltay , deltaz ! == Lx / nfft1 , Ly / nfft2 , Lz / nfft3

integer(i2b):: species ! dummy for loops over species from 1 to nb_species



! tell user

write (*,*) '*************************'

write (*,*) 'Compute personnal_vext'

! look for tag 'hard_sphere_solute_radius' which gives the radius of the sphere for the external potential
solute_radius=input_dp('hard_sphere_solute_radius')
write (*,*) 'the radius of these solute hard spheres is defined as a constant for now and is ', solute_radius

! compute vext hard sphere

deltax = Lx / real ( nfft1 , kind = dp )

deltay = Ly / real ( nfft2 , kind = dp )

deltaz = Lz / real ( nfft3 , kind = dp )

! check if radius is allocated. if not there is a bug.

if ( .not. allocated ( radius ) ) then

  write (*,*) 'radius is not allocated. should be. critial stop in compute_vext_perso.f90'

  stop

end if

! compute the sum of solute and solvent radius

do species = 1 , nb_species

radius_sum_sq = ( solute_radius + radius ( species ) ) ** 2

do n = 1 , nb_solute_sites

  do i = 1 , nfft1

    x_grid = real ( i - 1 , kind = dp) * deltax

    x_nm2 = ( x_grid - x_mol ( n ) )**2

    do j = 1 , nfft2

      y_grid = real ( j - 1 , kind = dp ) * deltay

      y_nm2 = ( y_grid - y_mol ( n ) ) ** 2

      do k = 1 , nfft3

        z_grid = real ( k - 1 , kind = dp ) * deltaz

        z_nm2 = ( z_grid - z_mol ( n ) ) ** 2

        r_nm2 = x_nm2 + y_nm2 + z_nm2

        ! test if distance is smaller than sum of radius of solvent and solute

        if ( r_nm2 <= radius_sum_sq ) then

          ! we're inside the hard potential. vext is very high

          Vext_total ( i , j , k , : , : , species ) = huge ( 1.0_dp )

          ! if not, vext is unchanged

        end if

      end do

    end do

  end do

end do ! solutes sites

end do ! species

end subroutine compute_vext_perso

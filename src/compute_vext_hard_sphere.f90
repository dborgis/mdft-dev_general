!> Compute external potential in a basic way
!> Solvation of a hard sphere in a hard sphere fluid
SUBROUTINE compute_Vext_hard_sphere
USE precision_kinds , only : i2b , dp
use input , only : input_line,input_dp
use system , only : nfft1 , nfft2 , nfft3 , x_mol , y_mol , z_mol , nb_solute_sites , Lx , Ly , Lz , radius , nb_species &
                    , deltax , deltay , deltaz
use external_potential , only : Vext_total , Vext_hard
use quadrature, only: angGrid, molRotGrid
IMPLICIT NONE
integer(i2b):: i , j , k !> @var dummy
integer(i2b):: species ! dummy between 1 and nb_species
integer(i2b):: solute_site ! dummy between 1 and nb_solute_sites
real(dp):: x_grid , y_grid , z_grid !> @var coordinates of grid mesh nodes
real(dp):: x_nm2 , y_nm2 , z_nm2 , r_nm2  ! distance between solute and grid point
real(dp):: hard_sphere_solute_radius ! radius of the hard sphere solute. tag in dft.in
real(dp):: sum_of_solute_and_solvent_radius ! sum of solute and solvent radius
! tell user
write (*,*) '*************************'
write (*,*) '>>> Compute HS Vext for hard sphere in hs fluid'
write (*,*) 'there are as many hard sphere solute as solutes in solute.in'
! be sure Vext_total is allocated
if ( .not. allocated ( Vext_total ) ) then
  allocate ( Vext_total ( nfft1 , nfft2 , nfft3 , angGrid%n_angles , molRotGrid%n_angles , nb_species ) )
  Vext_total = 0.0_dp
END IF
! look for tag 'hard_sphere_solute_radius' which gives the radius of the sphere for the external potential
hard_sphere_solute_radius=input_dp('radius_of_hard_sphere_solute')
write (*,*) 'the radius of these solute hard spheres is defined as a constant for now and is ', hard_sphere_solute_radius
 
!check if radius is allocated. if not it may be because you wand an hard sphere solute but no hardsphere solvent, in that case
!radius=0, if you want an hard sphere solvent there is a bug
if ( .not. allocated ( radius ) ) then
allocate(radius(1))
print*,'!!!!!!!!!!!!WARNING!!!!!!!!!!!!!!!!!!WARNING!!!!!!!!!!!!!!!!!!WARNING!!!!!!!!!!!!!!!!WARNING!!!!!!!!!!!!!!WARNING!!'
print*,'WARNING solvent radius is not defined in comput_vext_hard_sphere it is now equal to 0, which is ok ONLY if you want '
print*,'an HS solute but no HS Solvent'
print*,'!!!!!!!!!!!!WARNING!!!!!!!!!!!!!!!!!!WARNING!!!!!!!!!!!!!!!!!!WARNING!!!!!!!!!!!!!!!!WARNING!!!!!!!!!!!!!!WARNING!!'
radius(1)=0.00_dp
END IF
write(*,*) 'radius' , radius
! compute the sum of solute and solvent radius
do solute_site = 1 , nb_solute_sites
  do species = 1 , nb_species
    sum_of_solute_and_solvent_radius = ( hard_sphere_solute_radius + radius ( species ) ) ** 2
    do k = 1 , nfft3
      z_grid = real( k - 1 , dp ) * deltaz
      z_nm2 = ( z_grid - z_mol ( solute_site ) ) ** 2
      do j = 1 , nfft2
        y_grid = real ( j - 1 , dp ) * deltay
        y_nm2 = ( y_grid - y_mol ( solute_site ) ) ** 2
        do i = 1 , nfft1
          x_grid = real ( i - 1 , dp ) * deltax
          x_nm2 = ( x_grid - x_mol ( solute_site ) ) ** 2
          r_nm2 = x_nm2 + y_nm2 + z_nm2
          ! test if distance is smaller than sum of radius of solvent and solute
          if ( r_nm2 <= sum_of_solute_and_solvent_radius ) then
            ! we're inside the hard potential. vext is very high
            Vext_total ( i , j , k , : , : , species ) = huge ( 1.0_dp ) ! no dependancy over Omega here
            ! if not , Vext_total is unchanged
          END IF
        END DO
      END DO
    END DO
  END DO
END DO
END SUBROUTINE compute_Vext_hard_sphere

! compute the external potential as created by a hard cylinder solute
! hard cylinder is along axe z=0, x=Lx/2, y=Ly/2
subroutine compute_Vext_hard_cylinder
use precision_kinds , only : i2b , dp
use input , only : input_line, input_dp
use system , only : nfft1 , nfft2 , nfft3 , x_mol , y_mol , nb_solute_sites , Lx , Ly , radius , nb_omega , nb_species, nb_psi
use external_potential , only : Vext_total
implicit none
integer(i2b):: i , j , n ! dummy
real(dp):: x_nm2 , y_nm2 , r_nm2 ! distance**2 between solute and grid point
real(dp):: hard_cylinder_radius ! radius of the hard cylinder solute
real(dp):: deltax , deltay ! == Lx / nfft1 , Ly / nfft2 , Lz / nfft3
real(dp):: sum_rad2 ! sum of solute cylinder radius and radius of solvent hard sphere
integer(i2b):: species ! dummy between 1 and nb_species
! init variables
deltax = Lx / real ( nfft1 , dp )
deltay = Ly / real ( nfft2 , dp )
! be sure Vext_total is allocated
if ( .not. allocated ( Vext_total ) ) then
  allocate ( Vext_total ( nfft1 , nfft2 , nfft3 , nb_omega , nb_psi, nb_species ) )
  Vext_total = 0.0_dp
end if
! tell user
write(*,*)'>>> Compute HS Vext Daniel for cylinder in z dir in hs fluid'
! look for tag 'hard_cylinder_radius' in dft.in
hard_cylinder_radius=input_dp('radius_of_hard_cylinder')
write ( * , * ) 'the hard cylinder solute has a radius of ' , hard_cylinder_radius
write ( * , * ) 'this cylinder is along Z. position in x and y are extracted from input/solute.in x_mol and y_mol'
write ( * , * ) 'there are ' , nb_solute_sites , ' cylinders in the system. They all have same radius'
! compute if grid point is or not inside the hard cylinder
! hard cylinder is along axe x=Lx/2, y=Ly/2
do n = 1 , nb_solute_sites
  do species = 1 , nb_species
    sum_rad2 = ( hard_cylinder_radius + radius ( species ) ) ** 2
    do j = 1 , nfft2
      y_nm2 = ( real(j-1,dp) * deltay - y_mol (n) )**2
      do i = 1 , nfft1
        x_nm2 = ( real(i-1,dp) * deltax - x_mol (n) )**2
        ! relative distances between solute site and solvent site
        r_nm2 = x_nm2 + y_nm2
        ! if we are inside , vext is huge, else vext is unchanged
        if ( r_nm2 <= sum_rad2 ) then
          Vext_total ( i , j , : , : , : , species ) = huge ( 1.0_dp )
        end if
      end do
    end do
  end do
end do
end subroutine compute_Vext_hard_cylinder

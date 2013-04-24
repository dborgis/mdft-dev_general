! This subroutine reads the positions of the charges of the solute sites. It extrapolates them to the grid points in order to 
! get a charge density rho_c ( i , j , k )

! Written by Maximilien Levesque, 2011 12 11 in Daniel Borgis's group, Ecole Normale Superieure, Paris.




subroutine charge_density_from_point_charge_positions

! Precision_kinds defines the precisions of each type
use precision_kinds , only : i2b , dp
! i2b = integer simple precision
! dp = real double precision

use system , only : nfft1 , nfft2 , nfft3 , chg_mol , id_mol , Lx , Ly , Lz , rho_c , nb_solute_sites , deltax , deltay , deltaz , &
                    x_mol , y_mol , z_mol
! nffti = nb of grid nodes in each direction
! chg_mol = charge of point charges
! Lx = box size
! rho_c = charge density
! nb_solute_sites = number of sites of the solute
! deltax = Lx / nfft1




implicit none

integer ( kind = i2b ) :: solute ! dummy

real ( kind = dp ) :: xq , yq, zq ! coordinates of the charge in indicial coordinates

integer ( kind = i2b ) :: im , jm , km , ip , jp , kp ! indices of corner in indicial coordinates

real ( kind = dp ) :: wim , wjm , wkm , wip , wjp , wkp ! weight associated to each index

character ( 50 ) :: filename




! We know the number of grid points since "allocate from input". We may allocate the charge density here

allocate ( rho_c ( nfft1 , nfft2 , nfft3 ) )

rho_c = 0.0_dp




!print*, x_mol , y_mol , z_mol

! extrapolate each solute point charge to grid nodes

do solute = 1 , nb_solute_sites

  ! if the solute does not have charge, go to next solute

  if ( chg_mol ( id_mol ( solute ) ) == 0.0_dp ) cycle


  ! transform cartesian coordinates 0 <= x < Lx in 'indicial coordinates' 0 <= xq < nfft1

  xq = modulo ( x_mol ( solute ) , Lx ) / deltax
  yq = modulo ( y_mol ( solute ) , Ly ) / deltay
  zq = modulo ( z_mol ( solute ) , Lz ) / deltaz


  ! get coordinates of grid node juste below (corner of the cube with smallest indices)
  ! +1 is because indexation does not begin to 0 but to 1. Thus, when coordinate xq=0, it corresponds to index 1.

  im = int ( xq ) + 1
  jm = int ( yq ) + 1
  km = int ( zq ) + 1


  ! grid node juste above (corner of the cube with highest indices)

  ip = im + 1
  jp = jm + 1
  kp = km + 1


  ! this corner should never have coordinates higher than nfft

  if ( ip == nfft1 + 1 ) ip = 1
  if ( jp == nfft2 + 1 ) jp = 1
  if ( kp == nfft3 + 1 ) kp = 1


  ! define weights associated with each corner

  wim = (    1.0_dp - (   xq - real(int(xq,kind=i2b),kind=dp)   )    )
  wjm = (    1.0_dp - (   yq - real(int(yq,kind=i2b),kind=dp)   )    )
  wkm = (    1.0_dp - (   zq - real(int(zq,kind=i2b),kind=dp)   )    )
  wip = (             (   xq - real(int(xq,kind=i2b),kind=dp)   )    )
  wjp = (             (   yq - real(int(yq,kind=i2b),kind=dp)   )    )
  wkp = (             (   zq - real(int(zq,kind=i2b),kind=dp)   )    )


  ! increase density accordingly

  rho_c ( im , jm , km ) = rho_c ( im , jm , km ) + chg_mol ( id_mol ( solute ) ) * wim * wjm * wkm
  rho_c ( ip , jm , km ) = rho_c ( ip , jm , km ) + chg_mol ( id_mol ( solute ) ) * wip * wjm * wkm
  rho_c ( im , jp , km ) = rho_c ( im , jp , km ) + chg_mol ( id_mol ( solute ) ) * wim * wjp * wkm
  rho_c ( im , jm , kp ) = rho_c ( im , jm , kp ) + chg_mol ( id_mol ( solute ) ) * wim * wjm * wkp
  rho_c ( ip , jp , km ) = rho_c ( ip , jp , km ) + chg_mol ( id_mol ( solute ) ) * wip * wjp * wkm
  rho_c ( ip , jm , kp ) = rho_c ( ip , jm , kp ) + chg_mol ( id_mol ( solute ) ) * wip * wjm * wkp
  rho_c ( im , jp , kp ) = rho_c ( im , jp , kp ) + chg_mol ( id_mol ( solute ) ) * wim * wjp * wkp
  rho_c ( ip , jp , kp ) = rho_c ( ip , jp , kp ) + chg_mol ( id_mol ( solute ) ) * wip * wjp * wkp

end do





! charge density is in charge per unit volume

rho_c = rho_c * real ( nfft1 * nfft2 * nfft3 , kind = dp ) / ( Lx * Ly * Lz )

! Write charge density in .cube file

filename = 'output/charge_density.cube'
call write_to_cube_file ( rho_c , filename )



end subroutine charge_density_from_point_charge_positions

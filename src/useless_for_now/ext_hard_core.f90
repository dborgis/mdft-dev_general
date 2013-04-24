! This subroutine creates a hard core potential at each atom. It has several interests.
! First is not to let point charges to get too close (potential catastrophe)
! Second is to do pure hard calculations
! The hard cores have radius of the van der Walls radius

subroutine ext_hard_core


! contains real etc kinds related stuff
use precision_kinds , only : dp , i2b
!dp : real double precision

! contains all informations given in periodic table
use periodic_table , only : ptable
! ptable = for instance ptable(1)%name gives the name of the element of Z=1

use constants , only : infty

! contains all informpations concerning the system we're studying
use system , only : nfft1 , nfft2 , nfft3 , nb_species , nb_solute_sites , x_mol , y_mol , z_mol , deltax , deltay , deltaz ,&
                    nb_omega , radius
! nffti = nbr of grid nodes
! nb_species = number of total species
! nb_solute_sites = nb of sites in the solute
! x_mol , y_mol , z_mol = coordinates of the sites of the solute

! contains informations relative to the external potential
use external_potential , only : vext_hard_core
! vext_hard_core = external potential due to what we are calculation


! contains input informations
use input , only : input_line, input_log
! input_line = table containing all 



implicit none


integer ( kind = i2b ) :: species ! dummy

integer ( kind = i2b ) :: solute ! dummy

real ( kind = dp ) :: d_hard ! vdw_radius of the solute site + radius of the solvent

integer ( kind = i2b ) :: i , j , k ! dummy

real ( kind = dp ) :: xs , ys , zs , x , y , z ! local dummy

real ( kind = dp ) :: time0 , time1 ! timers



stop 'Vexthardcore'



! time

call cpu_time ( time0 )




! test if one wants 
if (.not. input_log('vdw_hard_core')) then 
return
end if






! for now this routine only works for hard sphere fluids, not for anything else, overall if there are orientations
if (.not. input_log ('hard_sphere_fluid')) then
print *, 'problem in compute_final_vext.f90 implementation only for nb_omega = 1 and pure hard sphere fluid'
stop
end if

if ( nb_omega /= 1 ) then
  print *, 'problem in compute_final_vext.f90 implementation only for nb_omega = 1 and pure hard sphere fluid'
  stop
end if







! init

if ( .not. allocated ( vext_hard_core ) ) allocate ( vext_hard_core ( nfft1 , nfft2 , nfft3 , nb_species ) )

vext_hard_core = 0.0_dp





! check if radius has been allocated correctly

if ( .not. allocated ( radius ) ) then
  print *, '**************************************************************'
  print *, 'radius is not allocated in ext_hard_core.f90 where is computed the vdw hard core repulsion'
  print *, '**************************************************************'
  stop
end if




! for each solvent species

do solute = 1 , nb_solute_sites

  xs = x_mol (solute)
  ys = y_mol (solute)
  zs = z_mol (solute)

  do species = 1 , nb_species

    d_hard = ptable(solute) % vdw_radius + radius ( species ) ! TODO we may improve this by using the fact that if a species is too big to be at a point, all species with bigger radius won't come neither at this point.
    d_hard = d_hard ** 2

    do k = 1 , nfft3

      z = real(k,kind=dp) * deltaz - zs
      z = z ** 2

      do j = 1 , nfft2

        y = real(j,kind=dp) * deltay - ys
        y = y ** 2

        do i = 1 , nfft1

          x = real(i,kind=dp) * deltax - xs
          x = x ** 2

          if ( x+y+z <= d_hard ) vext_hard_core (i,j,k,species) = infty

        end do ! k

      end do ! j

    end do ! i

  end do ! species

end do ! solute






! time

call cpu_time ( time1 )


print *, 'ext_hard_core took (sec) ' , time1 - time0






end subroutine ext_hard_core

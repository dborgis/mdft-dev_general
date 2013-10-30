! this subroutine gets the partial vext ( hard sphere + hard wall + hard cylinder + purely repulsive)
! and adds to it Vext_lj and Vext_q.
! then it gives a upper value (100 kJ/mol) to vext_total.
subroutine vext_total_sum
use precision_kinds , only : dp , i2b
use system , only : nfft1 , nfft2 , nfft3 , nb_species
use quadrature, only : angGrid, molRotGrid
use constants , only : fourpi
use external_potential , only : Vext_total , Vext_lj , Vext_q , vext_hard_core
! vext_total = total external potential used for minimization
! vext_lj = lennard jones part
! vext_q = electrostatic part
! vext_hard_core = vdw hard core repulsion
use quadrature , only : sym_order
implicit none
character ( len = 50 ) :: filename ! dummy
real(dp), dimension ( nfft1 , nfft2 , nfft3 ) :: temparray ! dummy
! be sure Vext_total is allocated
if ( .not. allocated ( Vext_total ) ) then
  print*,'V_ext_total is allocated in vext_total_sum'
  allocate ( Vext_total ( nfft1 , nfft2 , nfft3 , angGrid%n_angles , molRotGrid%n_angles, nb_species ) )
  Vext_total = 0.0_dp
else 
  print*, 'Warning V_ext_total has been allocated before vext_total_sum if you do not want to use HS you may have a problem'
end if
!    Vext_total = 0.0_dp
! vext is the sum over all external potentials
! note that purely repulsive and hard potentials are already included in vext
if ( allocated ( Vext_q ) ) Vext_total = Vext_total + Vext_q
if ( allocated ( Vext_lj        ) ) Vext_total             = Vext_total             + Vext_lj
if ( allocated ( vext_hard_core ) ) vext_total (:,:,:,1,1,:) = vext_total (:,:,:,1,1,:) + vext_hard_core ! TODO generalize
! limit Vext to 100 kJ/mol
print*, 'Vext is limited to 100 kJ/mol' ! which is enormous == 10^28 Kelvin
where ( Vext_total > 100.0_dp ) Vext_total = 100.0_dp
if ( allocated ( vext_lj        ) ) where ( Vext_lj        > 100.0_dp ) Vext_lj        = 100.0_dp
if ( allocated ( vext_hard_core ) ) where ( vext_hard_core > 100.0_dp ) vext_hard_core = 100.0_dp
! give vext extrema to user for visual debugging
print *, minval ( Vext_total ) , ' < Vext_total < ' , maxval ( Vext_total )
print *, minval ( Vext_lj ) , ' < Vext_lj < ' , maxval ( Vext_lj )
if ( allocated (Vext_q) ) then
print *, minval ( Vext_q ) , ' < Vext_q < ' , maxval ( Vext_q )
end if
if ( allocated (Vext_hard_core) ) then
print *, minval ( Vext_hard_core ) , ' < Vext_hard_core < ' , maxval ( Vext_hard_core )
end if
! mean over orientations and print
call mean_over_orientations ( Vext_total ( : , : , : , : , : , 1 ) , temparray )
temparray = temparray / (fourpi**2/(2.0_dp*sym_order))
filename = 'output/Vext.cube' ! care when HS or multispec
call write_to_cube_file ( temparray , filename )
filename = 'output/Vext_along-z.dat' ! care when HS or multispec'
call compute_z_density ( temparray , filename )
!! Only Vext_total and Vext_q are used in the following. Vext_q only in the initiation of the density, and Vext_total in the functional
!! We may deallocate Vext_lj.
!if ( allocated ( Vext_lj ) ) deallocate ( Vext_lj )
end subroutine vext_total_sum

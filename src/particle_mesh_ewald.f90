subroutine particle_mesh_ewald


use system , only : nfft1 => nfftx , nfft2 => nffty , nfft3 => nfftz , Lx , Ly, Lz

use precision_kinds , only : dp

implicit none


include 'sizes.i'
include 'cutoff.i'
include 'warp.i'
include 'atoms.i'
include 'bound.i'
include 'boxes.i'
include 'cell.i'
include 'charge.i'
include 'chgpot.i'
include 'couple.i'
include 'energi.i'
include 'ewald.i'
include 'group.i'
include 'math.i'
include 'shunt.i'
include 'usage.i'
include 'pme.i'


! here sort of interface for preparing echarge

use_smooth = .false. ! warp.i

use_ewald = .true. ! cutoff.i


! n is the total number of atoms (sites)

n = nb_mol + nb_solv


! nion is the total number of partial charges in the system

call compute_total_number_of_partial_charges_in_system ( nion )


! put partial charges in table pchg ( nion )

pchg ( 1 : nb_mol ) = chg_mol ( 1 : nb_mol )

pchg ( nb_mol + 1 : nb_solv ) = chg_solv ( 1 : nb_solv )


! put coordinates in table x y z ( nion )

x ( 1 : nb_mol ) = x_mol ( 1 : nb_mol )

x ( nb_mol + 1 : nb_solv ) = x_solv ( 1 : nb_solv )

y ( 1 : nb_mol ) = y_mol ( 1 : nb_mol )

y ( nb_mol + 1 : nb_solv ) = y_solv ( 1 : nb_solv )

z ( 1 : nb_mol ) = z_mol ( 1 : nb_mol )

z ( nb_mol + 1 : nb_solv ) = z_solv ( 1 : nb_solv )


! flag to use partitioning of system into groups

use_group = .false.


!flag to mark presence of infinite polymer

use_polymer = .false.


! fft with other tool

nfft1 = nfftx
nfft2 = nffty
nfft3 = nfftz


! volume in Ang**3 of the periodic box

volbox = Lx * Ly * Lz

xbox = Lx

ybox = Ly

zbox = Lz

xbox2 = 0.5_dp * Lx

ybox2 = 0.5_dp * Ly

zbox2 = 0.5_dp * Lz


! use_bounds    flag to use periodic boundary conditions

use_bounds = .true.


!assigns particle mesh Ewald parameters and options for a periodic system
call kewald

! everything is ready so call the charge-charge energy computer
call echarge




CONTAINS





!===================================================================================================================================
!  compute the total number of partial charges in system: nion 
!===================================================================================================================================
subroutine compute_total_number_of_partial_charges_in_system ( nion )
use precision_kinds , only : i2b
use system , only : nb_solute_sites , chg_mol , nb_solvent_sites , chg_solv
implicit none
integer(i2b) , intent(out) :: nion ! total number of partial charges in system
integer(i2b) :: i ! dummy
! init
nion = 0
! count in solute
do i = 1 , nb_solute_sites
  if ( chg_mol ( i ) /= 0.0_dp ) nion = nion + 1
end do
! count in solvent
do i = 1 , nb_solvent_sites
  if ( chg_solv ( i ) /= 0.0_dp ) nion = nion + 1
end do
! tell user
print *, 'Total number of partial charges in system = ' , nion 
end subroutine compute_total_number_of_partial_charges_in_system


end subroutine particle_mesh_ewald

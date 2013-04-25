!> read solvent atomic positions, charge, and lennard jones values in solvent.in
!! charge in electron units, sigma in Angstroms, epsilon in KJ/mol.
subroutine read_solvent
use precision_kinds , only : i2b , dp
use system , only : nb_solvent_sites , x_solv , y_solv , z_solv , chg_solv , sig_solv , eps_solv , id_solv
implicit none
integer(i2b):: n ! dummy
integer(i2b):: stat ! status du fichier ouvert
integer(i2b):: nb_id_solv ! number of types of solvent sites
real(dp):: total_charge ! charge of the whole system (should be zero in almost every case)
! open input/solvent.in and check if it is readable
open ( 5 , file = 'input/solvent.in' , status = 'old' , iostat = stat )
! if not readable ( stat /= 0 ) then stop simu and tell user
if ( stat /= 0 ) then
  write (*,*) 'solvent.in cannot be opened ! => STOP !'
  stop
end if
! pass first line of comments
read (5,*)
!> Second line is the total number of atom sites of the solvent AND the total number of different types of atoms
read (5,*) nb_solvent_sites , nb_id_solv
! allocate accordingly
allocate ( x_solv ( nb_solvent_sites ) )
allocate ( y_solv ( nb_solvent_sites ) )
allocate ( z_solv ( nb_solvent_sites ) )
allocate ( id_solv ( nb_solvent_sites ) ) ! from solvent_site to id for instance id_solv(1)=1, id_solv(2)=2 and id_solv(3)=2 for OH2
allocate ( chg_solv ( nb_id_solv ) )
allocate ( sig_solv ( nb_id_solv ) )
allocate ( eps_solv ( nb_id_solv ) )
! read rest of the file
! check the total charge of the solvent to WARN USER if it is not neutral
! init total charge of the solvent to zero
total_charge = 0.0_dp
read(5,*) ! comment line
do n = 1 , nb_solvent_sites
  read (5,*) id_solv ( n ) , chg_solv ( id_solv ( n ) ) , sig_solv ( id_solv ( n ) ) , &
             eps_solv ( id_solv ( n ) ) , x_solv ( n ) , y_solv ( n ) , z_solv ( n )
  total_charge = total_charge + chg_solv ( id_solv ( n ) )
end do
! close input/solvent.in
close(5)
end subroutine read_solvent

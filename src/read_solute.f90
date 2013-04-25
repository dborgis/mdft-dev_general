!> read solute atomic positions, charge, and lennard jones values in solute.in
!! charge in electron units, sigma in Angstroms, epsilon in KJ/mol.
subroutine read_solute

use precision_kinds , only : i2b,dp

use system , only : nb_solute_sites , x_mol , y_mol , z_mol , chg_mol , sig_mol , eps_mol , atomic_nbr , id_mol , Lx , Ly , Lz&
&, lambda1_mol , lambda2_mol

use input , only : input_line

use periodic_table , only : init_periodic_table , ptable




implicit none

integer(i2b):: n ! dummy

integer(i2b):: stat ! status du fichier ouvert

integer(i2b):: nb_id_mol ! number of different kinds of site (ie two LJ sites with same epsilon and sigma, but only diff positions)

real(dp):: system_charge ! charge of the whole system (should be zero in almost every case)







! init module periodic_table so that all informations are available

call init_periodic_table

!print *, ptable ( 1 ) % name






!> open and test if input/solute.in is ok

open ( 5 , file='input/solute.in', status='old', iostat=stat)

if (stat /= 0) then

  write(*,*) 'solute.in cannot be opened ! => STOP !'

  stop

end if

!> first line is a comment
read ( 5 , * )

!> Second line is the total number of atom sites of the solute AND the total number of different types of atoms
read (5,*) nb_solute_sites, nb_id_mol





!> Allocate accordingly

allocate(x_mol(nb_solute_sites))
allocate(y_mol(nb_solute_sites))
allocate(z_mol(nb_solute_sites))
allocate(id_mol(nb_solute_sites)) ! from solute_site to id for instance id_mol(1)=1, id_mol(2)=2 and id_mol(3)=2 for OH2
allocate(atomic_nbr(nb_solute_sites))
allocate(chg_mol(nb_id_mol))
allocate(sig_mol(nb_id_mol))
allocate(eps_mol(nb_id_mol))
!allocate lambda1_mol et lambda2_mol
allocate(lambda1_mol(nb_solute_sites))
allocate(lambda2_mol(nb_solute_sites))


! Read rest of the file and test the total chage

read(5,*) ! comment line

do n = 1 , nb_solute_sites

  read(5,*) id_mol(n) , chg_mol(id_mol(n)) , sig_mol(id_mol(n)) , eps_mol(id_mol(n)) ,lambda1_mol(n) , lambda2_mol(n) &
& , x_mol(n) , y_mol(n) , z_mol(n) ,atomic_nbr(n)
 

end do




! sum all charges

system_charge = 0.0_dp

do n = 1 , nb_solute_sites

  system_charge = system_charge + chg_mol ( id_mol ( n ) )

end do





! warn user if total charge is not zero

if ( abs ( system_charge ) >= 0.001 ) then

  write(*,*)
  write(*,*)'*******************************************************************'
  write(*,*)'          TOTAL CHARGE IS NOT ZERO.'
  write(*,*)'          IT IS ' , system_charge
  write(*,*)'          CONTINUE AT YOUR OWN RISK'
  write(*,*)'*******************************************************************'
  write(*,*)

end if




! close input/solute.in

close(5)


! ask user if he wants all the sites to be translated to the center of the box, ie by Lx/2, Ly/2, Lz/2

call translate_to_center_of_supercell_if_needed



! Print periodic XSF file to be read by VMD or equivalent

call print_supercell_xsf

! check if cartesian coordinates read in input/solute.in are in the supercell

call check_if_coordinates_are_in_supercell

!write(*,*), x_mol(1) , y_mol(1) , z_mol(1)
!write(*,*), x_mol(2) , y_mol(2) , z_mol(2)
!write(*,*), x_mol(3) , y_mol(3) , z_mol(3)

! x_mol(1) =  15.000000000000000_dp     
! y_mol(1) = 15.000999999999999_dp     
! z_mol(1) = 15.000000000000000_dp    
  
! x_mol(2) =  15.816495000000000_dp 
! y_mol(2) =     15.000999999999999_dp    
! z_mol(2) =   15.577352500000000_dp     
 
! x_mol(3)  = 14.183505000000000_dp    
! y_mol(3)  =  15.000999999999999_dp  
! z_mol(3) =  15.577352500000000_dp    

!write(*,*), x_mol(1) , y_mol(1) , z_mol(1)
!write(*,*), x_mol(2) , y_mol(2) , z_mol(2)
!write(*,*), x_mol(3) , y_mol(3) , z_mol(3)
!write(*,*), x_mol(4) , y_mol(4) , z_mol(4)
!write(*,*), x_mol(5) , y_mol(5) , z_mol(5)
!write(*,*), x_mol(6) , y_mol(6) , z_mol(6)


contains




! if user asks for it (tag 'translate_solute_to_center'), add Lx/2, Ly/2, Lz/2 to all solute coordinates
subroutine translate_to_center_of_supercell_if_needed
use precision_kinds , only : i2b,dp
use input , only : input_line,input_log
use system , only : Lx , Ly , Lz , x_mol , y_mol , z_mol
implicit none
integer(i2b) :: i , j
logical :: translate_solute_to_center

! then do the translation
!print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
!print*,'ATTENTION J ai Changé Translate Solute to Center Pour Test ne marche que selon z!!!!'
!print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
if (input_log( 'translate_solute_to_center' )) then
  x_mol = x_mol + Lx / 2.0_dp
  y_mol = y_mol + Ly / 2.0_dp
  z_mol = z_mol + Lz / 2.0_dp
end if
end subroutine translate_to_center_of_supercell_if_needed










! this subroutine checks if coordinates are in supercell, and if not, add i*Lx to the cartesian coordinate to put back in supercell
subroutine check_if_coordinates_are_in_supercell
use precision_kinds , only : i2b , dp
use system , only : nb_solute_sites , Lx , Ly , Lz , x_mol , y_mol , z_mol
implicit none
integer (i2b) :: i 
! check if some positions are out of the supercell

!j is a test tag. We loop over this test until every atom is in the box.
! This allows for instance, if a site is two boxes too far to still be ok.

forall ( i = 1 : nb_solute_sites )

  x_mol ( i ) = modulo ( x_mol ( i ) , Lx )
  y_mol ( i ) = modulo ( y_mol ( i ) , Ly )
  z_mol ( i ) = modulo ( z_mol ( i ) , Lz )

end forall

end subroutine check_if_coordinates_are_in_supercell




end subroutine read_solute

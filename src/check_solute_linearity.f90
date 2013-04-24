!> Check solute linearity
subroutine check_solute_linearity(islinear)

use precision_kinds, only: dp
use system, only: nb_solute_sites,x_mol,y_mol,z_mol

implicit none
logical,intent(out):: islinear
integer:: n,sim1,sim2,sim3 !> dummy


!> One site cannot be linear, while two are always, else ...
if (nb_solute_sites==1) then
  islinear=.false.
else if (nb_solute_sites==2) then
  islinear=.true.
else

!> Are there two common coordinates to all sites ?
   sim1=0
   sim2=0
   sim3=0
   do n=1,nb_solute_sites
     if (x_mol(n)==y_mol(n)) sim1=sim1+1
     if (x_mol(n)==z_mol(n)) sim2=sim2+1
     if (y_mol(n)==z_mol(n)) sim3=sim3+1
   end do
   if (sim1==nb_solute_sites .or. sim2==nb_solute_sites .or. sim3==nb_solute_sites) then
      islinear=.true.
   else
      islinear=.false.
   end if
end if

!> Inform user
if (islinear) write(*,*)'>>> Solute is found to be linear'

end subroutine

!> Check solute planarity
!> If solute is planar then is_planar = .true. else is_planar = .false.
subroutine check_solute_planarity(is_planar)

use precision_kinds, only: dp
use system, only: nb_solute_sites,x_mol,y_mol,z_mol

implicit none
logical,intent(out):: is_planar
real(dp),dimension(2:nb_solute_sites):: xvec,yvec,zvec !> vector between solute sites 1 and N (2<N<NB_SOLUTE_SITES)
integer :: n !> dummy


!> Can be planar only if number of solute sites is 3 or more
if (nb_solute_sites==1 .or. nb_solute_sites==2) then
   is_planar=.false.
else

!> Compute all vector coordinates between first site and every others
   forall (n=2:nb_solute_sites)
      xvec(n)=x_mol(1)-x_mol(n)
      yvec(n)=y_mol(1)-y_mol(n)
      zvec(n)=z_mol(1)-z_mol(n)
   end forall
   
!> Check if one coordinate is always the same (which is the case if we're in a plan)
   if ( minval(xvec)==maxval(xvec) .or. minval(yvec)==maxval(yvec) .or. minval(zvec)==maxval(zvec)) then
      is_planar=.true.
   else
      is_planar=.false.
   end if
end if

!> Inform user
if (is_planar) write(*,*)'>>> Solute is found to be planar'

end subroutine

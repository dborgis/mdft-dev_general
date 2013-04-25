! This subroutine tabulates the potential

subroutine potgen(pot_type,param1,param2,nrgrid,drgrid,pot)
use precision_kinds, only: dp,i2b
use system, only: Lx,Ly,Lz
use constants, only: qfact
implicit none
integer(i2b), intent(in) :: pot_type ! 1==lennard jones
real(dp), intent(in) :: param1,param2 ! parameters of the potential, for instance epsilon and sigma for lj
real(dp), intent(out), dimension(0:nrgrid) :: pot !tabulated potential this subroutines computes
integer(i2b), intent(in) :: nrgrid ! nb of point of the radial grid for tabulation of the potentials
real(dp), intent(in) :: drgrid ! distance between two points in radial grid (in Angstroms). should be  =(rmax-rmin)/nrgrid

integer(i2b) :: i ! dummy
real(dp) :: r_i ! radial distance


!> Lennard Jones
if (pot_type ==1) then
  ! param1==epsilon ; param2==sigma
  pot(0)=huge(1.0_dp)
  do i=1,nrgrid
    ! U = 4*e*(s/r)^6*( (s/r)^6 -1 )
    r_i=drgrid*real(i,dp)
    pot(i)=4.0_dp*param1*((param2/r_i)**12-(param2/r_i)**6)
    write(98,*)r_i,pot(i)
  end do


!> Coulomb sum
else if (pot_type ==2) then
  !param1=chg1, param2=chg2
  pot(0)=huge(1.0_dp)
  do i=1,nrgrid
    ! U=qq'/r
    r_i=drgrid*real(i,dp)
    pot(i)=qfact*param1*param2/r_i
  end do


!> test if pot_type is coherent
else
  write(*,*)'error in potgen.f90, pot_type value not expected'
  stop


end if

end subroutine potgen









!===================================================================================================================================
!===================================================================================================================================
!===================================================================================================================================

! This subroutine tabulates the potential as a function of r**2

subroutine potgensq ( pot_type , param1 , param2 , nrgrid , drgrid , pot )


use precision_kinds , only : dp , i2b
use system          , only : Lx , Ly , Lz
use constants       , only : qfact


implicit none


integer(i2b), intent(in)                   :: pot_type ! 1==lennard jones
real(dp), intent(in)                       :: param1 , param2 ! parameters of the potential, for instance epsilon and sigma for lj
integer(i2b), intent(in)                   :: nrgrid ! nb of point of the radial grid for tabulation of the potentials
real(dp), intent(in)                       :: drgrid ! distance between two points in radial grid (in Angstroms). should be  =(rmax-rmin)/nrgrid
real(dp), intent(out), dimension(0:nrgrid) :: pot !tabulated potential this subroutines computes
integer(i2b)                               :: i ! dummy
real(dp)                                   :: r_i ! radial distance
real(dp)                                   :: r_isq ! radial distance




! Lennard Jones
if ( pot_type == 1 ) then
  ! param1 == epsilon ; param2 == sigma
  pot(0) = huge ( 1.0_dp ) ! Vlj(r=0)=infty
  do i = 1 , nrgrid
    ! U = 4*e*(s/r)^6*( (s/r)^6 -1 )
    r_isq = drgrid * real( i ,dp)
    r_i   = sqrt( r_isq )
    pot(i)= 4.0_dp * param1 * ( (param2**2/r_isq)**6-(param2**2/r_isq)**3)
  end do
  write(*,*)'eps & sigma = ',param1,param2,' , Vlj min est ',minval(pot)
end if




! Coulomb 1/r
if (pot_type ==2) then
  !param1=chg1, param2=chg2
  pot(0)=huge(1.0_dp)
  do i=1,nrgrid
    ! U=qq'/r
    r_isq=drgrid*real(i,dp)
    r_i=sqrt(r_isq)
    pot(i)=qfact*param1*param2/r_i
  end do
end if



end subroutine potgensq

subroutine get_gauss_legendre_integration_roots_and_weights
use precision_kinds , only : dp , i2b
use input , only : input_line, input_int
use system , only : nb_legendre , nb_omega
use quadrature, only : w_legendre , x_legendre , Omx , Omy , Omz , weight
implicit none
real(dp),dimension(11,11) :: w_gl ! weights for gauss legendre integration of order 1 to 6
real(dp),dimension(11,11) :: x_gl ! zero for gauss legendre integration of order 1 to 6
integer(i2b):: i , j ! dummy
! get gauss legendre integration order nb_legendre from input_line
nb_legendre=input_int('order_of_quadrature')
!> Test if nb_legendre is not too high
if (nb_legendre > size(x_gl,1)) then
  write(*,*)'Problem detected in gauss_legendre_integration.f90'
  write(*,*)'nb_legendre > size(x_gl,1)'
  write(*,*)'critical. stop.'
  stop
end if
!> in case of no angular grid (nb_legendre=1) (only one orientation, the one defined in solvent.in), no need to consider psi
select case (nb_legendre)
case (1)
  nb_omega = 1
case default
  nb_omega = 2*nb_legendre**2
end select
! GET THEM FROM up to N=100
! http://www.holoborodko.com/pavel/numerical-methods/numerical-integration/#gauss_quadrature_abscissas_table
! N=1 ************************
x_gl(1,1) = 0.0_dp
w_gl(1,1) = 2.0_dp
! N=2 ************************
x_gl(1,2) = -0.5773502691896257645091488_dp
x_gl(2,2) = -x_gl(1,2)
w_gl(1,2) = 1.0_dp
w_gl(2,2) = w_gl(1,2)
! N=3 ************************
x_gl(1,3) = -0.7745966692414833770358531_dp
x_gl(2,3) = 0.0_dp
x_gl(3,3) = -x_gl(1,3)
w_gl(1,3) = 0.5555555555555555555555556_dp
w_gl(2,3) = 0.8888888888888888888888889_dp
w_gl(3,3) = w_gl(1,3)
!  N=4 *************************
x_gl(1,4) = -0.8611363115940525752239465_dp
x_gl(2,4) = -0.3399810435848562648026658_dp
x_gl(3,4) = -x_gl(2,4)
x_gl(4,4) = -x_gl(1,4)
w_gl(1,4) = 0.3478548451374538573730639_dp
w_gl(2,4) = 0.6521451548625461426269361_dp
w_gl(3,4) = w_gl(2,4)
w_gl(4,4) = w_gl(1,4)
!  N=5 ****************
x_gl(1,5) = -0.9061798459386639927976269_dp
x_gl(2,5) = -0.5384693101056830910363144_dp
x_gl(3,5) = 0.0_dp
x_gl(4,5) = -x_gl(2,5)
x_gl(5,5) = -x_gl(1,5)
w_gl(1,5) = 0.2369268850561890875142640_dp
w_gl(2,5) = 0.4786286704993664680412915_dp
w_gl(3,5) = 0.5688888888888888888888889_dp
w_gl(4,5) = w_gl(2,5)
w_gl(5,5) = w_gl(1,5)
!  N=6 *******************
x_gl(1,6) = -0.9324695142031520278123016_dp
x_gl(2,6) = -0.6612093864662645136613996_dp
x_gl(3,6) = -0.2386191860831969086305017_dp
x_gl(4,6) = -x_gl(3,6)
x_gl(5,6) = -x_gl(2,6)
x_gl(6,6) = -x_gl(1,6)
w_gl(1,6) = 0.1713244923791703450402961_dp
w_gl(2,6) = 0.3607615730481386075698335_dp
w_gl(3,6) = 0.4679139345726910473898703_dp
w_gl(4,6) = w_gl(3,6)
w_gl(5,6) = w_gl(2,6)
w_gl(6,6) = w_gl(1,6)
!  N=7 *******************
x_gl(1,7) = -0.9491079123427585245261897_dp
x_gl(2,7) = -0.7415311855993944398638648_dp
x_gl(3,7) = -0.4058451513773971669066064_dp
x_gl(4,7) = 0.0_dp
x_gl(5,7) = -x_gl(3,7)
x_gl(6,7) = -x_gl(2,7) 
x_gl(7,7) = -x_gl(1,7)
w_gl(1,7) = 0.1294849661688696932706114
w_gl(2,7) = 0.2797053914892766679014678
w_gl(3,7) = 0.3818300505051189449503698
w_gl(4,7) = 0.4179591836734693877551020
w_gl(5,7) = w_gl(3,7)
w_gl(6,7) = w_gl(2,7)
w_gl(7,7) = w_gl(1,7)
!  N=8 *******************
x_gl(1,8) = -0.9602898564975362316835609_dp
x_gl(2,8) = -0.7966664774136267395915539_dp
x_gl(3,8) = -0.5255324099163289858177390_dp
x_gl(4,8) = -0.1834346424956498049394761_dp
x_gl(5,8) = -x_gl(4,8)
x_gl(6,8) = -x_gl(3,8)
x_gl(7,8) = -x_gl(2,8)
x_gl(8,8) = -x_gl(1,8)
w_gl(1,8) = 0.1012285362903762591525314_dp
w_gl(2,8) = 0.2223810344533744705443560_dp
w_gl(3,8) = 0.3137066458778872873379622_dp
w_gl(4,8) = 0.3626837833783619829651504_dp
w_gl(5,8) = w_gl(4,8)
w_gl(6,8) = w_gl(3,8)
w_gl(7,8) = w_gl(2,8)
w_gl(8,8) = w_gl(1,8)
!  N=9 *******************
x_gl(1,9) = -0.9681602395076260898355762_dp
x_gl(2,9) = -0.8360311073266357942994298_dp
x_gl(3,9) = -0.6133714327005903973087020_dp
x_gl(4,9) = -0.3242534234038089290385380_dp
x_gl(5,9) = 0.0_dp
x_gl(6,9) = -w_gl(4,9)
x_gl(7,9) = -x_gl(3,9)
x_gl(8,9) = -x_gl(2,9)
x_gl(9,9) = -x_gl(1,9)
w_gl(1,9) = 0.0812743883615744119718922_dp
w_gl(2,9) = 0.1806481606948574040584720_dp
w_gl(3,9) = 0.2606106964029354623187429_dp
w_gl(4,9) = 0.3123470770400028400686304_dp
w_gl(5,9) = 0.3302393550012597631645251_dp
w_gl(6,9) = w_gl(4,9)
w_gl(7,9) = w_gl(3,9)
w_gl(8,9) = w_gl(2,9)
w_gl(9,9) = w_gl(1,9)
!  N=10 *******************
x_gl(1,10) = -0.9739065285171717200779640_dp
x_gl(2,10) = -0.8650633666889845107320967_dp
x_gl(3,10) = -0.6794095682990244062343274_dp
x_gl(4,10) = -0.4333953941292471907992659_dp
x_gl(5,10) = -0.1488743389816312108848260_dp
x_gl(6,10) = -x_gl(5,10)
x_gl(7,10) = -x_gl(4,10)
x_gl(8,10) = -x_gl(3,10)
x_gl(9,10) = -x_gl(2,10)
x_gl(10,10)= -x_gl(1,10)
w_gl(1,10) = 0.0666713443086881375935688_dp
w_gl(2,10) = 0.1494513491505805931457763_dp
w_gl(3,10) = 0.2190863625159820439955349_dp
w_gl(4,10) = 0.2692667193099963550912269_dp
w_gl(5,10) = 0.2955242247147528701738930_dp
w_gl(6,10) = w_gl(5,10)
w_gl(7,10) = w_gl(4,10)
w_gl(8,10) = w_gl(3,10)
w_gl(9,10) = w_gl(2,10)
w_gl(10,10)= w_gl(1,10)
!  N=11 *******************
x_gl(1,11) = -0.9782286581460569928039380_dp
x_gl(2,11) = -0.8870625997680952990751578_dp
x_gl(3,11) = -0.7301520055740493240934163_dp
x_gl(4,11) = -0.5190961292068118159257257_dp
x_gl(5,11) = -0.2695431559523449723315320_dp
x_gl(6,11) = 0.0_dp
x_gl(7,11) = -x_gl(5,11)
x_gl(8,11) = -x_gl(4,11)
x_gl(9,11) = -x_gl(3,11)
x_gl(10,11)= -x_gl(2,11)
x_gl(11,11)= -x_gl(1,11)
w_gl(1,11) = 0.0556685671161736664827537_dp
w_gl(2,11) = 0.1255803694649046246346943_dp
w_gl(3,11) = 0.1862902109277342514260976_dp
w_gl(4,11) = 0.2331937645919904799185237_dp
w_gl(5,11) = 0.2628045445102466621806889_dp
w_gl(6,11) = 0.2729250867779006307144835_dp
w_gl(7,11) = w_gl(5,11)
w_gl(8,11) = w_gl(4,11)
w_gl(9,11) = w_gl(3,11)
w_gl(10,11)= w_gl(2,11) 
w_gl(11,11)= w_gl(1,11)
allocate (w_legendre(nb_legendre,nb_legendre))
w_legendre=w_gl(1:nb_legendre,1:nb_legendre)
allocate (x_legendre(nb_legendre,nb_legendre))
x_legendre=x_gl(1:nb_legendre,1:nb_legendre)
end subroutine get_gauss_legendre_integration_roots_and_weights

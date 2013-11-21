! in this SUBROUTINE is evaluated the Weeks Chandler Anderson diameter at which the Lennard Jones potential is cutted into
! a reference potential (the Hard sphere potential) and the attractive potential)
! here we use only the fit by Verlet et Weiss but more robust methods exist and are coded in a SUBROUTINE called compute_optimal_hs_diameter.f90
! be carefull as d_wca is in reduced units. for instance if you use it for lennard jones perturbations, the hard sphere diameter is d_wca calculated here x sigmalj
! 201109121800 Maximilien Levesque   Creation
! 201109122158 Maximilien Levesque   Better notations, equations, explicit interface for function, more coherent declaration of kb
! 201109182054 Maximilien Levesque   Comment about reduced units
SUBROUTINE compute_wca_diameter ( n_0 , temp , sig , eps , d_wca )
use precision_kinds , only : dp , i2b
use constants , only : Boltz , Navo
IMPLICIT NONE
real(dp), intent(in) :: temp ! temperature
real(dp), intent(in) :: n_0 ! density
real(dp), intent(in) :: sig ! sigma of LJ potential
real(dp), intent(in) :: eps ! epsilon of LJ potential
real(dp), intent(out) :: d_wca ! Weeks Chandler Anderson cutoff diameter between the reference potential and the WCA LJ potential
real(dp), parameter :: kb = Boltz * Navo / 1000.0_dp ! \approx 8.3145e-3, Boltzman constant in kJ/mol
real(dp):: beta ! 1/(kb.T)
real(dp):: temp_reduite ! reduced temperature
real(dp):: db ! Baker Anderson parameter
real(dp), parameter :: a1 = 0.3837_dp , a2 = 0.4293_dp , a3 = 1.068_dp , a4 = 210.31_dp , a5 = 404.6_dp ! parameter from Verlet & Weis
real(dp):: delta_b ! parameter from Verlet & Weis
real(dp):: d_old ! backup value of d_wca during convergence
! init
beta = 1.0_dp / ( kb * temp )
! Reduced temperature
temp_reduite = 1.0_dp / ( beta * eps )
write(*,*)'T*= ',temp_reduite
! Estimate the Baker Anderson parameter using the fit by Verlet and Weis, Phys. Rev. A, 1972
! one could do the integration numericaly (you find it in another routine) but the fit is
! really accurate (see the article and reviews of it). The error is estimated under 10e-4
db = ( a1 + a3 * beta ) / ( a2 + beta )
! delta_b par estimation fit
delta_b = 1.0_dp / ( a4 + a5 * beta )
! d_wca is determined by convergence. It is initiated as db (which is not dependant on sigma nor epsilon LJ) and corrected up to
! some convergence criteria
! Initiate d to db
d_wca = db
! backup newest d_wca in d_old
d_old = d_wca
! compute new d_wca
d_wca = db * ( 1.0_dp + delta_b * sig1_over_2sig0 ( n_0 , d_wca ) )
do while ( abs ( d_wca - d_old ) > 1.0e-10_dp )
  d_old = d_wca
  d_wca = db * ( 1.0_dp + delta_b * sig1_over_2sig0 ( n_0 , d_wca ) )
END DO
! go back to real units
d_wca = d_wca * sig
write(*,*)'WCA calculated hard sphere diameter is d = ',d_wca
contains
  function sig1_over_2sig0 ( n_0 , d_wca )
  use precision_kinds , only : dp
  use constants , only : pi
  IMPLICIT NONE
  real(dp):: sig1_over_2sig0
  real(dp), intent(in) :: n_0 ! fluid density
  real(dp), intent(in) :: d_wca ! iterative d
  real(dp), parameter :: s1 = -17.0_dp / 4.0_dp , s2 = 1.362_dp , s3 = -0.8751_dp
  real(dp):: eta ! packing fraction
  real(dp):: eta_w ! corrected packing fraction
  eta = pi / 6.0_dp * n_0 * d_wca ** 3
  eta_w = eta - eta ** 2 / 16.0_dp
  sig1_over_2sig0 = ( 1.0_dp + s1 * eta_w + s2 * eta_w ** 2 + s3 * eta_w ** 3 ) / ( 1.0_dp - eta_w ) ** 2
  end function sig1_over_2sig0
END SUBROUTINE compute_wca_diameter
	

! This module defines the fundamental constants
! It has to be the most complete possible and very well documented in order to be easily trackable and modified and shareable
! The best website to my knowledge to get this kind of information is the National Institute of Standards and Technology (NIST)
! http://physics.nist.gov/cuu/Constants/
! In order to be easily sharable and trackable, please add date and short description of your modification below.
! 20110922 16h08 Maximilien Levesque         creation from a basis by Daniel Borgis
! 20111120 20h00 Maximilien Levesque         addition sqpi ! square root of pi
! 20111216 08h36 Maximilien Levesque         addition of infty = huge(1.0_dp)
module constants ! fondamental_constants
USE precision_kinds , only : dp
IMPLICIT NONE
real(dp)   , parameter :: ln2           = log ( 2.0_dp )
real(dp)   , parameter :: pi            = acos ( -1.0_dp )
real(dp)   , parameter :: twopi         = 2.0_dp * pi
real(dp)   , parameter :: fourpi        = 4.0_dp * pi
real(dp)   , parameter :: eightpi       = 8.0_dp * pi
real(dp)   , parameter :: sqpi          = sqrt ( pi ) ! square root of pi
real(dp)   , parameter :: eps0          = 8.854187817e-12_dp ! Electric constant ie Vacuum permittivity in Farad per meter == Coulomb^2 per (Joule x meter)
real(dp)   , parameter :: qunit         = 1.602176565e-19_dp ! Elementary charge in Coulomb ; [C]
real(dp)   , parameter :: Navo          = 6.02214129e23_dp ! Avogadro constant in mol^(-1)
real(dp)   , parameter :: qfact         = qunit ** 2 * 1.0e-3_dp * Navo / ( fourpi * eps0 * 1.0e-10_dp ) ! electrostatic potential unit so that QFACT*q*q/r is kJ/mol
real(dp)   , parameter :: Boltz         = 1.3806488e-23_dp ! Boltzmann constant in Joule per Kelvin, [J].[K]^{-1}
complex(dp), parameter :: i_complex     = ( 0.0_dp , 1.0_dp )
real(dp)   , parameter :: infty         = huge ( 1.0_dp )
real(dp)   , parameter :: zero = 0._dp, one = 1.0_dp, two = 2.0_dp, three = 3.0_dp, four = 4.0_dp, five = 5.0_dp
end module constants

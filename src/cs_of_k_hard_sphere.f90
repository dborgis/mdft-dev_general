! This SUBROUTINE computes the radial distribution function c_s (k) 
! For now it only works for 1 SPECIES. BEWARE !!!
SUBROUTINE cs_of_k_hard_sphere
USE precision_kinds , only : i2b , dp
! i2b for integer simple precision
! dp for real double precision
use system , only : radius , n_0_multispec , c_s_hs , nb_species
! radius = radius(nb_species) = radius of each species
! n_0_multispec (nb_species) = reference bulk density (0.0332891 molecule per angsrom^3 for water, etc).
! c_s_hs (nb_k) = tabulated values of cs in kspace it is the goal of this routine to tabulate this.
use input , only : input_line
! input_line is an array containing dft.in
use constants , only : fourpi , pi
! fourpi = 4pi
! pi = 3.14159...
IMPLICIT NONE
real (dp) :: phi ! excess free energy density
real (dp) :: n0, n1, n2, n3 ! weighted densities
real (dp) :: w0, w1, w2, w3 ! weight functions
real (dp) :: R ! ==radius (for easier equations)
real (dp) :: FourPiR ! FourPi*radius
real (dp) :: kR ! k*R
real (dp) :: coskR, sinkR ! cos(k*R), sin(k*R)
real (dp), dimension (0:3,0:3) :: d2phi ! second partial derivative of phi wrt ni nj
real (dp) :: kmin, kmax ! min and max values of k
integer (i2b) :: i , j ! dummy
integer (i2b) :: ios ! input output status
integer (i2b) :: nb_k ! number of k points
real (dp) :: value1, value2 ! dummy
real (dp) :: k ! k-space vector norm
logical :: PY , CS ! the one which is true is the right equation of state
! Look at hard sphere radius
if ( .not. allocated ( radius ) ) call read_hard_sphere_radius_and_allocate_if_necessary
! Works only for one implicit species
call is_it_only_one_species
! Check if you want to use Perkus-Yevick or Carnahan Starling equation of state
call do_we_use_cs_or_py_eos ( PY , CS )
! Here we could generate as many points as wanted. In order to be coherent, we will use the same number of points as in cs.in
!> read the total number of lines in input/cs.in (which is the same as in input/cd.in and input/cdelta.in
open(11,file='input/cs.in')
nb_k=0
do while(.true.)
  read(11,*,iostat=ios)
  if (ios>0) then
    write(*,*)'Error in compute_ck_dipolar.f90'
    write(*,*)'something went wrong during the computation of the total number of lines in cs.in. stop'
    stop
  ELSE IF (ios<0) then
    ! end of file reached
    exit
  ELSE
    nb_k=nb_k+1
  END IF
END DO
close(11)
print*, nb_k , '#########@@@@@@@@@@@################'
! allocate accordingly
allocate ( c_s_hs ( nb_k ) )
! Same reasoning for kmax
! read the distance between two k points in input/cs.in  (which is the same as in input/cd.in and input/cdelta.in
open(11,file='input/cs.in')    !TODO AUTOMIZATION
read(11,*)value1, k
read(11,*)value2, k
close(11)
kmin = 0.0_dp
kmax = ( value2 - value1 ) * real ( nb_k - 1 ,dp)
kloop : do i = 0, nb_k-1
  k = real(i,dp)*(kmax-kmin)/real(nb_k-1,dp) ! pay attention to k=0
  ! weight functions
  R = radius ( 1 )
  kR = k*R
  sinkR = sin(kR)
  coskR = cos(kR)
  FourPiR = FourPi*R
  
  if (i == 0) then ! k=0
    w0 = 1.d0
    w1 = R
    w2 = FourPi *R**2
    w3 = FourPi/3.d0 *R**3
    ! weighted densities
    n0 = n_0_multispec ( 1 ) * w0
    n1 = n_0_multispec ( 1 ) * w1
    n2 = n_0_multispec ( 1 ) * w2
    n3 = n_0_multispec ( 1 ) * w3
  ELSE ! k/=0
    w0 = coskR + .5d0*kR*sinkR
    w1 = (sinkR + kR*coskR) / (2.d0*k)
    w2 = (FourPiR*sinkR) /k
    w3 = FourPi*(sinkR-kR*coskR) / k**3
  END IF
  
  ! expression of phi_exc depends obviously of the choice of the excess functional
  if ( PY ) then ! PY
    ! first partial derivative of phi wrt n_\alpha, case PY
    !dphi(0) = -Log(1.0d0-n3)
    !dphi(1) = n2/(1.0d0-n3)
    !dphi(2) = n1/(1.0d0-n3) + n2**2/(8.d0*(1.d0-n3)**2*Pi)
    !dphi(3) = (n1*n2)/(1.d0 - n3)**2 + n0/(1.d0 - n3) + n2**3/(12.d0*(1.d0 - n3)**3*Pi)
    
    ! second partial derivative of phi wrt n_\alpha n_\beta
    d2phi(0,0) = 0.d0
    d2phi(0,1) = 0.d0
    d2phi(0,2) = 0.d0
    d2phi(0,3) = 1.d0/(1.d0 - n3)
    
    d2phi(1,0) = 0.d0
    d2phi(1,1) = 0.d0
    d2phi(1,2) = 1.d0/(1.d0 - n3)
    d2phi(1,3) = n2/(1.d0 - n3)**2
    
    d2phi(2,0) = 0.d0
    d2phi(2,1) = 1.d0/(1.d0 - n3)
    d2phi(2,2) = n2/(4.d0*(1.d0 - n3)**2*Pi)
    d2phi(2,3) = n1/(1.d0 - n3)**2 + n2**2/(4.d0*(1.d0 - n3)**3*Pi)
  
    d2phi(3,0) = 1.d0/(1.d0 - n3)
    d2phi(3,1) = n2/(1.d0 - n3)**2
    d2phi(3,2) = d2phi(2,3)
    d2phi(3,3) = (2.d0*n1*n2)/(1.d0 - n3)**3 + n0/(1.d0 - n3)**2 + n2**3/(4.d0*(1.d0 - n3)**4*Pi)
  ELSE IF ( CS ) then ! CS
    ! first partial derivative of phi wrt n_\alpha, case PY
    !dphi(0) = -Log(1.0d0-n3)
    !dphi(1) = n2/(1.0d0-n3)
    !dphi(2) = n1/(1 - n3) + n2**2/(12.*(1 - n3)**2*n3*Pi) + (n2**2*Log(1 - n3))/(12.*n3**2*Pi)
    !dphi(3) = (n1*n2)/(1 - n3)**2 - (-n0 + n2**3/(36.*n3**2*Pi))/(1 - n3) - n2**3/(36.*(1 - n3)**2*n3**2*Pi) + n2**3/(18.*(1 - n3)**3*n3*Pi) - (n2**3*Log(1 - n3))/(18.*n3**3*Pi)
    
    ! second partial derivative of phi wrt n_\alpha n_\beta
    d2phi(0,0) = 0.d0
    d2phi(0,1) = 0.d0
    d2phi(0,2) = 0.d0
    d2phi(0,3) = 1.d0/(1.d0 - n3)
    
    d2phi(1,0) = 0.d0
    d2phi(1,1) = 0.d0
    d2phi(1,2) = 1.d0/(1.d0 - n3)
    d2phi(1,3) = n2/(1.d0 - n3)**2
    
    d2phi(2,0) = 0.d0
    d2phi(2,1) = 1.d0/(1.d0 - n3)
    d2phi(2,2) = n2/(6.d0*(1.d0 - n3)**2*n3*Pi) + (n2*Log(1.d0 - n3))/(6.d0*n3**2*Pi) ! simplified by Mathematica
    d2phi(2,3) = (n3*(n2**2*(2.d0 - 5.d0*n3 + n3**2) + 12.d0*n1*(-1.d0 + n3)*n3**2*Pi) &
                 - 2.d0*n2**2*(-1.d0 + n3)**3*Log(1.d0 - n3))/(12.d0*(-1.d0 + n3)**3*n3**3*Pi)  ! simplified by Mathematica
  
    d2phi(3,0) = 1.d0/(1.d0 - n3)
    d2phi(3,1) = n2/(1.d0 - n3)**2
    d2phi(3,2) = d2phi(2,3)
    d2phi(3,3) = (n3*(n2**3*(6.d0 - 21.d0*n3 + 26.d0*n3**2 - 5.d0*n3**3) - 72.d0*n1*n2*(-1.d0 + n3)*n3**3*Pi &
                 + 36.d0*n0*(-1.d0 + n3)**2*n3**3*Pi) + 6.d0*n2**3*(-1.d0 + n3)**4*Log(1.d0 - n3))/(36.d0*(-1.d0 + n3)**4*n3**4*Pi) ! simplified by Mathematica
  END IF ! only PY or CS for now
  
  ! direct correlation function
  c_s_hs ( i + 1 ) =( d2phi(0,0)*w0*w0 &
                    + d2phi(0,1)*w0*w1 &
                    + d2phi(0,2)*w0*w2 &
                    + d2phi(0,3)*w0*w3 &
                    + d2phi(1,0)*w1*w0 &
                    + d2phi(1,1)*w1*w1 &
                    + d2phi(1,2)*w1*w2 &
                    + d2phi(1,3)*w1*w3 &
                    + d2phi(2,0)*w2*w0 &
                    + d2phi(2,1)*w2*w1 &
                    + d2phi(2,2)*w2*w2 &
                    + d2phi(2,3)*w2*w3 &
                    + d2phi(3,0)*w3*w0 &
                    + d2phi(3,1)*w3*w1 &
                    + d2phi(3,2)*w3*w2 &
                    + d2phi(3,3)*w3*w3 ) * (-1.d0)
  write(99,*) k, c_s_hs ( i + 1 )
END DO kloop
contains
SUBROUTINE do_we_use_cs_or_py_eos ( PY , CS )
USE precision_kinds , only : i2b
use input , only : input_line
IMPLICIT NONE
logical , intent( out ) :: CS , PY
integer(i2b) :: i , j! dummy
CS = .false.
PY = .false.
do i = 1 , size ( input_line )
  j = len ( 'hs_functional' )
  if ( input_line (i) (1:j) == 'hs_functional' .and. input_line (i) (j+4:j+5) == 'CS' ) CS = .true.
  if ( input_line (i) (1:j) == 'hs_functional' .and. input_line (i) (j+4:j+5) == 'PY' ) PY = .true.
END DO
END SUBROUTINE do_we_use_cs_or_py_eos
! This SUBROUTINE checks if the number of implicit species is not different from 1
SUBROUTINE is_it_only_one_species
use system , only : nb_species
IMPLICIT NONE
if ( nb_species /= 1 ) then
  print *, 'SUBROUTINE cs_of_k_hard_sphere.f90 which is used to generate cs(k) for the hard sphere is only written for one species.'
  print *, 'sorry.'
  stop
END IF
END SUBROUTINE is_it_only_one_species
END SUBROUTINE cs_of_k_hard_sphere

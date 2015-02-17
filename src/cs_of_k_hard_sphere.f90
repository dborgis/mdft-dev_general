! This subroutine computes the direct correlation function of a hard sphere fluid

SUBROUTINE cs_of_k_hard_sphere

USE precision_kinds, ONLY: i2b, dp
USE system, ONLY: nb_species, solvent !@GUILLAUME c_s_hs should go into MODULE DCF
USE input, ONLY: input_line, n_linesInFile, verbose
USE constants, ONLY: fourpi, pi
USE dcf, ONLY: c_s, chi_l, c_s_hs
USE hardspheres ,ONLY: hs
use mathematica, only: spline, splint

IMPLICIT NONE

real (dp) :: phi ! excess free energy density
real (dp) :: n0, n1, n2, n3 ! weighted densities
real (dp) :: w0, w1, w2, w3 ! weight functions
real (dp) :: R ! ==radius
real (dp) :: FourPiR ! FourPi*radius
real (dp) :: kR ! k*R
real (dp) :: coskR, sinkR ! cos(k*R), sin(k*R)
real (dp), dimension (0:3,0:3) :: d2phi ! second partial derivative of phi wrt ni nj
real (dp) :: kmin, kmax ! min and max values of k
integer (i2b) :: i , j, nb_k, ios
real (dp) :: value1, value2 ! dummy
real (dp) :: k ! k-space vector norm
logical :: PY , CS ! the one which is true is the right equation of state

if ( .not. allocated ( hs ) ) call read_hard_sphere_radius_and_allocate_if_necessary
write(*,'(A,F12.4)')'I will compute direct correlation function for hard sphere fluid with radius',hs(1)%R
call is_it_only_one_species ! Works only for one implicit species
call do_we_use_cs_or_py_eos ( PY , CS )! Check if you want to use Perkus-Yevick or Carnahan Starling equation of state

! Here we could generate as many points as wanted. In order to be coherent, we will use the same number of points as in cs.in
! read the total number of lines in input/cs.in (which is the same as in input/cd.in and input/cdelta.in

! IF (.NOT. ((ALLOCATED(c_s_hs%x)) .OR. (ALLOCATED(chi_l)))) THEN
!     PRINT*, "YOU WANT TO COMPUTE correlation function for HS, i.e you want to compute HSB, but you have not included any other"
!     PRINT*, " excess FUNCTIONAL. THIS IS NON-SENSE"
!     STOP
! END IF

nb_k = size(c_s%x)
allocate( c_s_hs%x(nb_k), source=0._dp) ! k
allocate( c_s_hs%y(nb_k), source=0._dp) ! c(k) hard sphere
allocate( c_s_hs%y2(nb_k), source=0._dp)! dc(k)/dk
do i = 1, nb_k
  k = c_s%x(i)
  c_s_hs%x(i) = k
  ! weight functions
  R = hs(1)%R
  kR = k*R
  sinkR = sin(kR)
  coskR = cos(kR)
  FourPiR = FourPi*R

  if (abs(k)<=epsilon(1._dp)) then ! k=0
    w0 = 1.d0
    w1 = R
    w2 = FourPi *R**2
    w3 = FourPi/3.d0 *R**3
  ELSE ! k > 0
    w0 = coskR + .5d0*kR*sinkR
    w1 = (sinkR + kR*coskR) / (2.d0*k)
    w2 = (FourPiR*sinkR) /k
    w3 = FourPi*(sinkR-kR*coskR) / k**3
  END IF

  ! weighted densities
  n0 = solvent(1)%n0 * w0
  n1 = solvent(1)%n0 * w1
  n2 = solvent(1)%n0 * w2
  n3 = solvent(1)%n0 * w3

  ! expression of phi_exc depends obviously of the choice of the excess functional
  if ( PY ) then ! PY
    ! first partial derivative of phi wrt n_\alpha, case PY
    !dphi(0) = -Log(1.0d0-n3)
    !dphi(1) = n2/(1.0d0-n3)
    !dphi(2) = n1/(1.0d0-n3) + n2**2/(8.d0*(1.d0-n3)**2*Pi)
    !dphi(3) = (n1*n2)/(1.d0 - n3)**2 + n0/(1.d0 - n3) + n2**3/(12.d0*(1.d0 - n3)**3*Pi)

    ! second partial derivative of phi wrt n_\alpha n_\beta
    d2phi(0,0) = 0.d0 ; d2phi(0,1) = 0.d0 ; d2phi(0,2) = 0.d0             ; d2phi(0,3) = 1.d0/(1.d0 - n3)
    d2phi(1,0) = 0.d0 ; d2phi(1,1) = 0.d0 ; d2phi(1,2) = 1.d0/(1.d0 - n3) ; d2phi(1,3) = n2/(1.d0 - n3)**2
    d2phi(2,0) = 0.d0 ; d2phi(2,1) = 1.d0/(1.d0 - n3) ;
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
  c_s_hs%y(i) =    -( d2phi(0,0)*w0*w0 &
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
                    + d2phi(3,3)*w3*w3 )
end do

open(99,file="output/cs_hs_tabulated.dat")
  do i=1,size(c_s_hs%x)
    write(99,*) c_s_hs%x(i), c_s_hs%y(i)
  end do
close(99)

! Prepare the future use of splines by calculating the second order derivative (y2) at each point
call spline( x=c_s_hs%x, y=c_s_hs%y, n=size(c_s_hs%x), yp1=huge(1._dp), ypn=huge(1._dp), y2=c_s_hs%y2)


! print the splined version of the direct correlation function to test both the spline function
open(14,file="output/cs_hs_spline.dat")
block
real(dp) :: x_loc, y_loc
do i=0,1000
  x_loc = i*0.1
  call splint( xa=c_s_hs%x, ya=c_s_hs%y, y2a=c_s_hs%y2, n=size(c_s_hs%y), x=x_loc, y=y_loc)
  write(14,*) x_loc, y_loc
end do
end block
close(14)



! cs is replaced by cs-cshs
block
real(dp) :: x_loc,y_loc
do i=1,size(c_s%x)
  x_loc=c_s%x(i)
  call splint( xa=c_s_hs%x, ya=c_s_hs%y, y2a=c_s_hs%y2, n=size(c_s_hs%y), x=x_loc, y=y_loc)
  c_s%y(i)= c_s%y(i) -y_loc
end do
call spline( x=c_s%x, y=c_s%y, n=size(c_s%x), yp1=huge(1._dp), ypn=huge(1._dp), y2=c_s%y2)
end block



!
! open(14,file="output/cs_analytic_PY_wertheim.dat")
!   real(dp) :: e, R, n
!   n = solvent(1)%n0
!   R = hs(1)%R*2 ! diameter
!   e = acos(-1._dp)/6. * R**3 * n
! close(14)

contains



SUBROUTINE do_we_use_cs_or_py_eos ( PY , CS )
    USE precision_kinds,only : i2b
    use input,only : input_line
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
    use system,only : nb_species
    IMPLICIT NONE
    if ( nb_species /= 1 ) then
        print *, 'SUBROUTINE cs_of_k_hard_sphere.f90 used to generate cs(k) for HS is only written for one solvant species.'
        print *, 'sorry.'
        stop
    END IF
END SUBROUTINE is_it_only_one_species




END SUBROUTINE cs_of_k_hard_sphere

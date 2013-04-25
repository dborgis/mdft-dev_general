! This program has been written by Maximilien Levesque
! for educational purpose at Ecole Normale Superieure, Paris
! on Friday the 27th, 2011.
! The idea here is to develop the easiest to read, simplest steepest descent
! code one can write. It is written in order everyone to be able
! to take it, modify it, rewritte one by him/herself.
! Please, use it and modify it as much as you want. Remarks and
! thanks :) are always welcome @ maximilien.levesque@gmail.com
! If you redistribute this file, please let at least my name and or the
! adress of my website in it.

!Steepest descent is something very common, very easy to understand,
!but everything I could find was written in fortran 77.
!For example I found this very nice website by XXXX
!www
!Fortran 77 is such a bad language to read (that's my feeling don't argue :))
!that I felt like obliged to write one in Fortran 90+. By 90+ I mean that
!you'll find below use of intrinsic functions up to fortran 2008 (norm2()).
!That's easy to read thanks to the use of intrinsic functions and implicit work
!on arrays which were dramaticaly absent in fortran 77. It's also easy 
!to modify because arrays are allocated dynamically. They are of dimension
!L so that L can be whatever you want.
!This is certainly not the most intelligent steepest descent you'll find.
!This is certainly not the best fortran programming you'll read.
!But it's to my knowledge the easiest to read and understand for
!physicists and chemists using fortran.
!You'll find a lot of comments and variables with intelligible names.
!You could find them too long, I find them clear.


! Program to find a multi dimensional minimum for non-mathematicians :P
! Steepest descent method
! Written by Maximilien Levesque, Ecole Normale Superieure, Paris
! timetag : 201105271742

!Let's begin with an example Find a local minimum of
!Y(x,y,z)=sin(x)+2*cos(y)-sin(z)

!This example is very common on the net, that's why I used it.
!In order to begin the reading of the program and its understanding,
!you just have to get that we're looking for the variables X(1),X(2),X(i)...
!so that F(X(1),X(2),...,X(i)...) is maximum.
!In our example, x is X(1), y is X(2) and z is X(3).
!Y is a scalar number. In the program you'll find Y to be an array of
!dimension 3 but that's because Y(2) and Y(1) are the value of Y at
!step before and step before before :). It's nothing more a 3 step memory.
!The last important object is dY which is the partial derivative of Y
!with respect to each X(i). It is thus a vector of the same dimension
!as X.

!number of dimensions: 3
!convergence criterion: 0.000000001
!maximum number of iterations: 50
!starting constant: 1
!input starting point X(1)=X(2)=X(3)=1
!RESULT is X(1)=1.5707963
!X(2)=-0.0000435
!X(3)=-1.5707963
!local maximum = 4.000000
!total number of steps: 30

subroutine steepest_descent_main(xdim,conv_criteria,itermax,X)
use precision_kinds,only: dp,i2b
implicit none
real(dp),dimension(xdim),intent(inout) :: X ! variables of Y
real(dp),dimension(xdim) :: XOLD ! backup of X
real(dp),dimension(3) :: Y ! function of X at 3 last iterations
real(dp),dimension(xdim) :: dY ! partial derivative of Y with respect to Xi
real(dp),intent(in) :: conv_criteria ! convergence criterion
real(dp) :: stepsize ! initial guess for step size (how fast we follow dY)
integer(i2b),intent(in) :: itermax ! maximum number of iterations
integer(i2b) :: n !dummy
integer(i2b),intent(in) :: xdim ! nb of dimensions of X

!> What is the initial guess for the variables ?
stepsize=10.0_dp

!> steepest descent
call steepds(xdim,conv_criteria,itermax,stepsize,X,XOLD,Y,dY,n)

!> Give user the results
write(*,*)'----- results -----'
!write(*,*)'Maximum found = ',eval_Y_and_dY(xdim,X,dY)
write(*,*)'Number of iterations = ',n

end subroutine steepest_descent_main








subroutine steepds(xdim,conv_criteria,itermax,stepsize,X,XOLD,Y,dY,n)
use precision_kinds,only: dp,i2b
implicit none
real(dp), parameter :: macheps=epsilon(1.0d0) ! machine precision given by fortran function epsilon
real(dp),intent(in) :: conv_criteria !convergence criteria
real(dp),intent(inout) :: stepsize! initial guess for step size (how fast we follow dY)
real(dp),dimension(xdim),intent(inout) :: X !variable vector 
real(dp),dimension(xdim),intent(inout) :: XOLD !old X
real(dp),dimension(xdim),intent(inout) :: dY(3) !partial derivatives
real(dp),dimension(3) :: Y(3) ! Y(X) at current step, step before, and step before before
real(dp) :: eval_Y_and_dY !eval_Y_and_dYuation of the function Y(X) and its partial derivatives dY(X)
integer(i2b) :: j!dummy
integer(i2b),intent(in) :: itermax
integer(i2b),intent(out) :: n
integer(i2b),intent(in) :: xdim!dimension of X

!Start initial probe
do j=1,3
  !obtain yy and dY(i)
  Y(j)=eval_Y_and_dY(xdim,X,dY)
  !update X(i)
  call updateX(xdim,stepsize,dY,X,XOLD)
end do
!We now have a history to base the subsequent search on

!Open file to write the convergence
open(10,file='output/iterations.dat',form='formatted')
write(10,*)'# iteration     F       ΔF      stepsize'

!At maximum do m iterations, if convergence is reached exit the loop
do n=1,itermax

  !Accelerate search if approach is monotonic
  if (abs(Y(2)-Y(1))>macheps .and. (Y(3)-Y(2))/(Y(2)-Y(1))>0.0d0) then
    if (stepsize <1000.d0) stepsize=stepsize*5.0d0
  end if

  !If heading the wrong way (here we're searching for the minimum)
  if(Y(3)>Y(2)) then
    !decelerate search
    stepsize=stepsize/100.d0
    !restore values
    X=XOLD
    !recall values of Y and dY for X=XOLD
    Y(3)=eval_Y_and_dY(xdim,X,dY)
  else !if head the right way
    Y(1)=Y(2)
    Y(2)=Y(3)
  end if
  
  !Update X(i)
  call updateX(xdim,stepsize,dY,X,XOLD)

  !Get Y(Xnew) and its derivatives
  Y(3)=eval_Y_and_dY(xdim,X,dY)

  write(*,*)'at iteration ',n,' ΔF =',Y(3)-Y(2),' stepsize = ',stepsize
  write(10,*)n,Y(3),Y(3)-Y(2),stepsize

  !Check for convergence : if converged, exit loop
  if (abs(Y(3)-Y(2))<conv_criteria) exit

end do

!close iteration output
close(10)
end subroutine steepds






!this function is the one to replace with any compute energy and gradient you wish for
function eval_Y_and_dY(xdim,X,dY)
use precision_kinds,only: dp,i2b
use cg, only: cg_vect,FF,dF
implicit none
real(dp) :: eval_Y_and_dY ! value of Y(X) and its partial derivatives dY
integer(i2b),intent(in) :: xdim !dimension of X
real(dp),dimension(xdim),intent(in) :: X !variables of Y
real(dp),dimension(xdim),intent(out) :: dY ! partial derivatives of Y
 cg_vect=X
!call compute_energy_and_gradients_DCF
call energy_and_gradient
eval_Y_and_dY=FF
dY=dF
end function eval_Y_and_dY








subroutine updateX(xdim,stepsize,dY,X,XOLD)
use precision_kinds,only: dp,i2b
integer(i2b),intent(in) :: xdim!dimension of variables X
real(dp),intent(in) :: stepsize ! initial guess for step size (how fast we follow dY)
real(dp),dimension(xdim),intent(inout) :: X!variables of Y
real(dp),dimension(xdim),intent(out) :: XOLD!backup of X
real(dp),dimension(xdim),intent(inout) :: dY!partial derivatives of Y with respect to Xi
integer(i2b) :: i!dummy
real(dp) :: dYnorm ! norm of the gradient
!Find the magnitude of the gradient (norm2 is fortran2008 language)
!dYnorm=norm2(dY)
dYnorm=0.0_dp
do i=1,xdim
  dYnorm=dYnorm+dY(xdim)**2
end do
dYnorm=sqrt(dYnorm)
!Backup the X(i)
XOLD=X
!Update the X(i)
X=X-stepsize*dY/dYnorm
end subroutine updateX

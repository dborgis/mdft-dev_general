!> This subroutine should minimize FF by finding the ground state of CGVECT array
! two methods for now: LBFGS and Conjugate Gradient
! LBFGS and CG+ seem equivalent in terms of velocity but
! LBFGS recquires much more memory than CG+ (around 8 times)
subroutine find_equilibrium_density
use precision_kinds , only : i2b , dp
use input , only : input_line
use cg
implicit none
integer(i2b):: iter ! step number ...
integer(i2b):: i , j ! dummy
real(dp):: time1 , time2
real(dp):: energy_before ! value of total energy (FF) at step actual and before for comparison
! get from input the type of minimization one wants to do
do i = 1 , size ( input_line )
  j = len ( 'minimizer' )
  if ( input_line (i) (1:j) == 'minimizer' ) minimizer_type = input_line (i) (j+4:j+8)
end do
! call the correct minimization technique
! if minimization technique is L-BFGS
if ( minimizer_type ( 1 : 4 ) == 'bfgs' ) then
  ! init
  iter = 0 ! id of the iteration. 0 for first iteration
  task = 'START' ! task is modified by lbfgs itself to communicate the state of the convergence
  ! setuplb finds the next step minima using LBFGS-B
  111 continue
  call cpu_time ( time1 ) ! start timer to get the time used for one loop
  call setulb ( ncg , mcg , cg_vect , ll , uu , nbd , FF , dF , factr , pgtol , wa , iwa , task , iprint , csave , lsave &
              , isave , dsave )
  if ( task (1:2) == 'FG' ) then ! continue one more step
    iter = iter + 1
    write (*,*)
    write (*,*)
    write (*,*) 'start L-BFGS iteration ' , iter
    if( iter > iterMAX ) go to 999 ! get out of minimizer.
    ! store energy at step i-1
    energy_before = FF
    ! compute energy at step i
    call energy_and_gradient
    call cpu_time(time2)
    ! tell user the difference in energy in between steps i-1 and i
    write (*,*) 'F(n)-F(n-1) = ' , FF - energy_before
    ! tell user the norm of the gradient vector
    write (*,*) 'norm2(dF) = ' , norm2 ( dF )
    ! tell user the time needed to make the iteration
    write (*,*) 'in (sec) =' , time2 - time1
    write (*,*)
    write (*,*)
    write (*,*)
    ! L-BFGS knows when it is converged but I prefer to control it here by myself ...
    if ( abs ( FF - energy_before ) <= epsg ) goto 999
    goto 111 ! minimization continues ...
  else if ( task (1:5) == 'NEW_X' ) then
    goto 111
  else
    if ( iprint <= -1 .and. task (1:4) /= 'STOP' ) write (6,*) task
  end if
  ! loops to find the minima is ended for one reason or the other
  999 continue
  ! close properly lbfgs
  deallocate (nbd)
  deallocate (iwa)
  deallocate (ll)
  deallocate (uu)
  deallocate (wa)
  deallocate (dF)
end if
! test if minimizer is cg+
!if ( minimizer_type ( 1 : 2 ) == 'cg' ) then
!  call driver ( ncg , cg_vect , FF , dF , epsg )
!end if
!! if minimizer is steepest descent
!if ( minimizer_type ( 1 : 2 ) == 'sd' ) then
!  call steepest_descent_main(ncg,epsg,itermax,cg_vect)
!end if 
end subroutine find_equilibrium_density
!====================================================================
subroutine interface_compute_energy_and_gradients ( n , x , f , g , stopouencore )
use precision_kinds , only : dp , i2b
use cg , only : ncg , cg_vect , FF , dF , minimizer_iter , itermax
implicit none
integer(i2b):: n
real(dp), dimension ( ncg ) , intent ( inout ) :: x
real(dp), intent ( inout ) :: f
real(dp), dimension ( ncg ) , intent ( inout ) :: g
character ( len = 10 ) , intent(out) :: stopouencore
minimizer_iter = minimizer_iter + 1
if ( minimizer_iter <= itermax ) then
  call energy_and_gradient
  n = ncg
  x = cg_vect
  f = FF
  g = dF
  stopouencore = 'ENCORE'
else if ( minimizer_iter > itermax ) then
  stopouencore = 'STOP'
end if
end subroutine interface_compute_energy_and_gradients

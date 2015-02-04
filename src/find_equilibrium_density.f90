! This SUBROUTINE should minimize FF by finding the ground state of CGVECT array
! two methods for now: LBFGS and Conjugate Gradient
! LBFGS and CG+ seem equivalent in terms of velocity but
! LBFGS recquires much more memory than CG+ (around 8 times)

SUBROUTINE find_equilibrium_density

  USE precision_kinds,only : i2b , dp
  USE minimizer, ONLY: ncg, mcg, cg_vect, ll, uu, nbd, FF, dF, factr, pgtol, wa, iwa, task, iprint, csave, lsave, isave , dsave,&
    minimizer_type, itermax, epsg, finalizeMinimizer, iter
  USE bfgs, ONLY: startBFGS => setulb
  use input, only: input_dp

  IMPLICIT NONE

  REAL(dp):: time1, time2
  REAL(dp):: energy_before, dfoverf, FFdiff

  if( minimizer_type(1:4)=='bfgs' ) then
    iter = 0
    task = 'START'
    111 continue
    call CPU_TIME ( time1 )
    call startBFGS (ncg, mcg, cg_vect, ll, uu, nbd, FF, dF, factr, pgtol, wa, iwa, task, iprint, csave, lsave, isave , dsave )

    if ( task (1:2) == 'FG' ) then ! continue one more step
      iter = iter + 1
      if ( iter > iterMAX ) go to 999 ! get out of minimizer.
      energy_before = FF
      call energy_and_gradient (iter)
      call cpu_time(time2)
      dfoverf = (FF-energy_before)/FF*100._dp
      ! L-BFGS knows when it is converged but I prefer to control it here by myself ...
      if ( abs ( FF - energy_before ) <= epsg ) goto 999
      goto 111 ! minimization continues ...
    else if ( task (1:5) == 'NEW_X' ) then
      goto 111
    else
      if ( iprint <= -1 .and. task (1:4) /= 'STOP' ) write (6,*) task
    end if



  else if(minimizer_type(1:2)=='SD'.or.minimizer_type(1:2)=='sd') then
    epsg = input_dp("epsg") ! we want a variation in energy between two consecutive steps to be less than epsg
    FF = 0._dp
    dF = 0._dp
    FFdiff=huge(1._dp) ! FFdiff can only be positive or zero
    iter = 0
    do while(FFdiff>epsg .and. iter<iterMAX)
      iter = iter+1
      FFdiff = FF
      call energy_and_gradient(iter) ! après ca on a le nouveau FF et le nouveau dF et le nouveau cg_vect
      cg_vect = cg_vect - 0.5_dp*dF
      FFdiff = abs(FFdiff-FF) ! old - new
    end do

  else
    stop "You ask for a minimizer that is not implemented."
  end if

  999 continue ! loops to find the minima is ended for one reason or the other

  print*,
  write(*,'(A,F12.2,A)') "FF before correction", FF," kJ/mol"

END SUBROUTINE find_equilibrium_density

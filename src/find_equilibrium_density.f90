!> This SUBROUTINE should minimize FF by finding the ground state of CGVECT array
! two methods for now: LBFGS and Conjugate Gradient
! LBFGS and CG+ seem equivalent in terms of velocity but
! LBFGS recquires much more memory than CG+ (around 8 times)

SUBROUTINE find_equilibrium_density

    USE precision_kinds,only : i2b , dp
    USE minimizer, ONLY: ncg, mcg, cg_vect, ll, uu, nbd, FF, dF, factr, pgtol, wa, iwa, task, iprint, csave, lsave, isave , dsave,&
                        minimizer_type, itermax, epsg, finalizeMinimizer
    USE bfgs, ONLY: startBFGS => setulb

    IMPLICIT NONE

    INTEGER(i2b):: iter
    REAL(dp):: time1, time2
    REAL(dp):: energy_before, dfoverf

  if (minimizer_type(1:4) /= 'bfgs') stop "You ask for a minimizer that is not implemented."

        iter = 0
        task = 'START'
        111 continue
        CALL CPU_TIME ( time1 )
        CALL startBFGS (ncg, mcg, cg_vect, ll, uu, nbd, FF, dF, factr, pgtol, wa, iwa, task, iprint, csave, lsave, isave , dsave )

        IF ( task (1:2) == 'FG' ) then ! continue one more step
            iter = iter + 1
            IF ( iter > iterMAX ) go to 999 ! get out of minimizer.
            energy_before = FF
            call energy_and_gradient (iter)
            call cpu_time(time2)
            dfoverf = (FF-energy_before)/FF*100._dp
            ! L-BFGS knows when it is converged but I prefer to control it here by myself ...
            if ( abs ( FF - energy_before ) <= epsg ) goto 999
            goto 111 ! minimization continues ...
        ELSE IF ( task (1:5) == 'NEW_X' ) then
            goto 111
        ELSE
            if ( iprint <= -1 .and. task (1:4) /= 'STOP' ) write (6,*) task
        END IF

        999 continue ! loops to find the minima is ended for one reason or the other
        print*,
        PRINT*,"FREE ENERGY =",FF,"kJ/mol"
        call adhoc_corrections_to_gsolv ! ad hoc corrections

        CALL finalizeMinimizer
END SUBROUTINE find_equilibrium_density

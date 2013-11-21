!> This SUBROUTINE should minimize FF by finding the ground state of CGVECT array
! two methods for now: LBFGS and Conjugate Gradient
! LBFGS and CG+ seem equivalent in terms of velocity but
! LBFGS recquires much more memory than CG+ (around 8 times)

SUBROUTINE find_equilibrium_density

    use precision_kinds , only : i2b , dp
    use input , only : input_line
    USE cg
    use bfgs, ONLY: startBFGS => setulb
    
    IMPLICIT NONE
    
    integer(i2b):: iter ! step number ...
    integer(i2b):: i , j ! dummy
    real(dp):: time1 , time2
    real(dp):: energy_before ! value of total energy (FF) at step actual and before for comparison
    
    ! get from input the type of minimization one wants to do
    DO i = 1 , size ( input_line )
        j = len ( 'minimizer' )
        IF ( input_line (i) (1:j) == 'minimizer' ) THEN
            minimizer_type = input_line (i) (j+4:j+8)
            EXIT
        END IF
    END DO
    
    ! call the correct minimization technique
    if ( minimizer_type ( 1 : 4 ) == 'bfgs' ) then
        iter = 0
        task = 'START' ! task is modified by lbfgs itself to communicate the state of the convergence
        111 continue
        CALL CPU_TIME ( time1 )
        CALL startBFGS (ncg, mcg, cg_vect, ll, uu, nbd, FF, dF, factr, pgtol, wa, iwa, task, iprint, csave, lsave, isave , dsave )

        if ( task (1:2) == 'FG' ) then ! continue one more step
            iter = iter + 1
            write (*,*)
            write (*,*)
            write (*,*)
            write (*,'(a,I4)') '# L-BFGS iteration ' , iter
            IF ( iter > iterMAX ) go to 999 ! get out of minimizer.
            energy_before = FF
            call energy_and_gradient
            call cpu_time(time2)

            write (*,'(''At iter '',i3,'' Fsolv= '',f9.3,'' kJ/mol converged of '',f11.6,'' w/ |dF|='',f11.6,'' in '',i5,'' s'')')&
                iter, FF, FF-energy_before, NORM2(dF), nint(time2-time1)
            
            ! L-BFGS knows when it is converged but I prefer to control it here by myself ...
            if ( abs ( FF - energy_before ) <= epsg ) goto 999
            goto 111 ! minimization continues ...
        ELSE IF ( task (1:5) == 'NEW_X' ) then
            goto 111
        ELSE
            if ( iprint <= -1 .and. task (1:4) /= 'STOP' ) write (6,*) task
        END IF
        ! loops to find the minima is ended for one reason or the other
        999 continue
        CALL finalizeMinimizer
    ELSE
        STOP "The minimizer you asked is not implemented yet."
    END IF
        ! test if minimizer is cg+
        !if ( minimizer_type ( 1 : 2 ) == 'cg' ) then
        !  call driver ( ncg , cg_vect , FF , dF , epsg )
        !END IF
        !! if minimizer is steepest descent
        !if ( minimizer_type ( 1 : 2 ) == 'sd' ) then
        !  call steepest_descent_main(ncg,epsg,itermax,cg_vect)
        !END IF 
END SUBROUTINE find_equilibrium_density


SUBROUTINE interface_compute_energy_and_gradients ( n , x , f , g , stopouencore )

    use precision_kinds , only : dp , i2b
    USE cg, ONLY: ncg , cg_vect , FF , dF , minimizer_iter , itermax

    IMPLICIT NONE
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
    ELSE IF ( minimizer_iter > itermax ) then
    stopouencore = 'STOP'
    END IF
    
END SUBROUTINE interface_compute_energy_and_gradients

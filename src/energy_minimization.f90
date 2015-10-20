subroutine energy_minimization

    use precision_kinds, only : i2b, dp
    use module_input, only: getinput
    use module_minimizer, only: lbfgsb
    use module_grid, only: grid
use module_solvent, only: solvent
    use module_minimizer, only: init_minimizer => init

    implicit none

    real(dp) :: f ! functional to minimize
    real(dp) :: df (grid%nx, grid%ny, grid%nz, grid%no, species(1)%nspec )

    !     Init minimization process. From allocations to optimizer parameters to guess solution
    call init_minimizer
    !     We start the iteration by initializing task.
    lbfgsb%task = 'START'
    !     The beginning of the loop
    iter = 0
    do while(lbfgsb%task(1:2).eq.'FG'.or.lbfgsb%task.eq.'NEW_X'.or.lbfgsb%task.eq.'START'.or.iter>itermax)
        iter = iter +1
        !     This is the call to the L-BFGS-B code.
        call setulb ( lbfgsb%n, lbfgsb%m, lbfgsb%x, lbfgsb%l, lbfgsb%u, lbfgsb%nbd, f, lbfgsb%g, lbfgsb%factr, &
        lbfgsb%pgtol, lbfgsb%wa, lbfgsb%iwa, lbfgsb%task, lbfgsb%iprint, lbfgsb%csave, lbfgsb%lsave, lbfgsb%isave,&
        lbfgsb%dsave )

        if (lbfgsb%task(1:2) == "FG") then
            call energy_and_gradient(isave(34), f, df) ! f and df are intent(out) of energy_and_gradient. see below for explanations of isave(34)
            !   change x and df to vector format
            i=0
            do is=1,species(1)%nspec; do io=1,grid%no; do iz=1,grid%nz; do iy=1,grid%ny; do ix=1,grid%nx
                i=i+1
                lbfgsb%g(i)=df(ix,iy,iz,io,is)
            end do;                   end do;          end do;          end do;          end do
            i=0
            do is=1,species(1)%nspec; do io=1,grid%no; do iz=1,grid%nz; do iy=1,grid%ny; do ix=1,grid%nx
                i=i+1
                lbfgsb%x(i)=solvent(is)%density(ix,iy,iz,io)
            end do;                   end do;          end do;          end do;          end do

        else if (lbfgsb%task(1:5) == 'NEW_X') then
            !       the minimization routine has returned with a new iterate.
            !       At this point have the opportunity of stopping the iteration
            !       or observing the values of certain parameters
            !
            !       First are two examples of stopping tests.

            !       Note: task(1:4) must be assigned the value 'STOP' to terminate
            !         the iteration and ensure that the final results are
            !         printed in the default format. The rest of the character
            !         string TASK may be used to store other information.

            !       1) Terminate if the total number of f and g evaluations
            !            exceeds 99.

            if (isave(34) >= itermax) task='STOP: TOTAL NO. of f AND g EVALUATIONS EXCEEDS LIMIT'

            !       2) Terminate if  |proj g|/(1+|f|) < 1.0d-10, where
            !          "proj g" denoted the projected gradient

            if (dsave(13) .le. 1.d-10*(1.0d0 + abs(f))) &
                task='STOP: THE PROJECTED GRADIENT IS SUFFICIENTLY SMALL'
        end if
    end do

    print*,
    write(*,'(A,F12.2,A)') "FF before correction", FF," kJ/mol"

end subroutine energy_minimization

subroutine energy_minimization

    use precision_kinds, only : dp
    use module_minimizer, only: lbfgsb, init_lbfgsb
    use module_grid, only: grid
    use module_solvent, only: solvent
    use module_energy_and_gradient, only: energy_and_gradient

    implicit none

    real(dp) :: f ! functional to minimize
    real(dp), target :: df (grid%no, grid%nx, grid%ny, grid%nz, solvent(1)%nspec )
    integer :: is, io, ix, iy, iz, i
    real :: time(1:10)

    print*,
    print*,
    print*, "===== Energy minimization ====="

    ! Init minimization process. From allocations to optimizer parameters to guess solution
    call init_lbfgsb

    do while(lbfgsb%task(1:2).eq.'FG'.or.lbfgsb%task.eq.'NEW_X'.or.lbfgsb%task.eq.'START')
time(:)=0
call cpu_time(time(1))
        !
        !     This is the call to the L-BFGS-B code. It updates the one-column vector lbfgsb%x.
        !
        call setulb ( lbfgsb%n, lbfgsb%m, lbfgsb%x, lbfgsb%l, lbfgsb%u, lbfgsb%nbd, f, lbfgsb%g, lbfgsb%factr, &
        lbfgsb%pgtol, lbfgsb%wa, lbfgsb%iwa, lbfgsb%task, lbfgsb%iprint, lbfgsb%csave, lbfgsb%lsave, lbfgsb%isave,&
        lbfgsb%dsave )

call cpu_time(time(2))

        if (lbfgsb%task(1:2) == "FG") then
            !
            ! lbfgs propose un nouveau vecteur colonne à tester, lbfgsb%x(:), qui correspond à un reshape de notre solvent(:)%density(:,:,:,:)
            ! Passons ce x dnas solvent%density qui va ensuite servir de base à un nouveau calcul de f et df
            ! Note pour plus tard : on pourrait utiliser la fonction RESHAPE de fortran
            !
call cpu_time(time(3))

            !
            ! l-bfgs-b%x is a vector (1dim array). We have to reshape it to fit into density(o,x,y,z)
            !
            do is=1,solvent(1)%nspec
                solvent(is)%density = reshape( lbfgsb%x, [grid%no, grid%nx, grid%ny, grid%nz] )
            end do

call cpu_time(time(4))

            !
            ! Given a density of orientation, xyz, return the value of the functional at this point F[rho] and the gradient dF/drho
            !
            call energy_and_gradient(f, df) ! f and df are intent(out) of energy_and_gradient. see below for explanations of isave(34)

call cpu_time(time(5))

            !
            !   Energy_and_gradient does not change solvent%density but updates df and gives the value of the functional, f
            !   lbfgs%g is one-column
            !
            lbfgsb%g = pack(df, .true.)


call cpu_time(time(6))
        end if
        call cpu_time(time(7))

        if (lbfgsb%task(1:2) == "FG") then
            print*, "]=><=[ Time to setulb: ", time(2)-time(1)
            print*, "]=><=[ Time to mv lbfgs%x into density(:): ", time(4)-time(3)
            print*, "]=><=[ Time to energy_and_gradient: ", time(5)-time(4)
            print*, "]=><=[ Time to mv new gradient to lbfgs%g: ", time(6)-time(5)
            print*, "]=><=[ Time to full cycle (evaluation + bfgs stuff): ", time(7)-time(1)
            print*, "==============================================================================="
            print*,
            print*,
        end if
    end do

end subroutine energy_minimization

subroutine energy_minimization

    use precision_kinds, only : dp
    use module_grid, only: grid
    use module_solvent, only: solvent
    use module_energy_and_gradient, only: energy_and_gradient
    use module_lbfgs_nocedal_mdft, only: lbfgsb

    implicit none

    real(dp) :: f ! functional to minimize
    real(dp) :: df (grid%no, grid%nx, grid%ny, grid%nz, solvent(1)%nspec )
    real :: time(1:10)
    integer :: i

    print*,
    print*,
    print*, "===== Energy minimization ====="

    ! Init minimization process. From allocations to optimizer parameters to guess solution
    call lbfgsb%init

    do while( (lbfgsb%task(1:2).eq.'FG'.or.lbfgsb%task.eq.'NEW_X'.or.lbfgsb%task.eq.'START') &
      .and. (lbfgsb%isave(30)<lbfgsb%itermax) )

time(:)=0
call cpu_time(time(1))
error stop "tralala l28"
        !
        !     This is the call to the L-BFGS-B code. It updates the one-column vector lbfgsb%x.
        !
        ! call setulb ( lbfgsb%n, lbfgsb%m, lbfgsb%x, lbfgsb%l, lbfgsb%u, lbfgsb%nbd, f, lbfgsb%g, lbfgsb%factr, &
        ! lbfgsb%pgtol, lbfgsb%wa, lbfgsb%iwa, lbfgsb%task, lbfgsb%iprint, lbfgsb%csave, lbfgsb%lsave, lbfgsb%isave,&
        ! lbfgsb%dsave )

        ! ATTENTION ici astuce de malade. Je passe le array (solvent%xi(:,:,:,:)) qui est un tableau à 4 colonnes.
        ! MAIS lbfgsb a besoin d'un vecteur, c'est à dire d'une seule colonne.
        ! MAIS il est con et mon array est contigu en mémoire, donc je lui passe mon tableau à 4 colonnes dans le bon sens
        ! et il ne s'en rend pas compte: c'est pour lui une seule colonne 4 fois plus grande. Joix du fortran dégueux :)

        call lbfgsb%setulb ( lbfgsb%n, lbfgsb%m, solvent(1)%xi, f, df, &
        lbfgsb%factr, lbfgsb%pgtol, lbfgsb%wa, lbfgsb%iwa, lbfgsb%task, lbfgsb%iprint, lbfgsb%csave, lbfgsb%lsave, lbfgsb%isave,&
        lbfgsb%dsave )

ERROR STOP "POOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOpokazdpoakzdpokazdpokazdpokazd========================"
call cpu_time(time(2))

        if (lbfgsb%task(1:2) == "FG") then
            !
            ! lbfgs propose un nouveau vecteur colonne à tester, lbfgsb%x(:), qui correspond à un reshape de notre solvent(:)%xi(:,:,:,:)
            ! Passons ce x dnas solvent%xi qui va ensuite servir de base à un nouveau calcul de f et df
            ! Note pour plus tard : on pourrait utiliser la fonction RESHAPE de fortran
            !
call cpu_time(time(3))

            !
            ! Given a density of orientation, xyz, return the value of the functional at this point F[xi] and the gradient dF/dxi
            !
            call energy_and_gradient(f, df) ! f and df are intent(out) of energy_and_gradient. see below for explanations of isave(34)

call cpu_time(time(4))

        end if


call cpu_time(time(5))

        if (lbfgsb%task(1:2) == "FG") then
            print*, "]=><=[ Time to setulb: ", time(2)-time(1)
            print*, "]=><=[ Time to energy_and_gradient: ", time(4)-time(3)
            print*, "]=><=[ Time to full cycle (evaluation + bfgs stuff): ", time(5)-time(1)
            print*, "==============================================================================="
            print*,
            print*,
        end if
    end do

end subroutine energy_minimization

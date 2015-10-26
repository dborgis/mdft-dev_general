subroutine energy_minimization

    use precision_kinds, only : i2b, dp
    use module_input, only: getinput
    use module_minimizer, only: lbfgsb
    use module_grid, only: grid
    use module_solvent, only: solvent
    use module_minimizer, only: init_minimizer => init

    implicit none

    real(dp) :: f ! functional to minimize
    real(dp) :: df (grid%nx, grid%ny, grid%nz, grid%no, solvent(1)%nspec )
    integer :: is, io, ix, iy, iz, i

    !     Init minimization process. From allocations to optimizer parameters to guess solution
    call init_minimizer
    !     We start the iteration by initializing task.
    lbfgsb%task = 'START'
    !     The beginning of the loop

    do while(lbfgsb%task(1:2).eq.'FG'.or.lbfgsb%task.eq.'NEW_X'.or.lbfgsb%task.eq.'START')
        !
        !     This is the call to the L-BFGS-B code. It updates the one-column vector lbfgsb%x.
        !
        call setulb ( lbfgsb%n, lbfgsb%m, lbfgsb%x, lbfgsb%l, lbfgsb%u, lbfgsb%nbd, f, lbfgsb%g, lbfgsb%factr, &
        lbfgsb%pgtol, lbfgsb%wa, lbfgsb%iwa, lbfgsb%task, lbfgsb%iprint, lbfgsb%csave, lbfgsb%lsave, lbfgsb%isave,&
        lbfgsb%dsave )

        if (lbfgsb%task(1:2) == "FG") then
            !
            ! lbfgs propose un nouveau vecteur colonne à tester, lbfgsb%x(:), qui correspond à un reshape de notre solvent(:)%density(:,:,:,:)
            ! Passons ce x dnas solvent%density qui va ensuite servir de base à un nouveau calcul de f et df
            !
            i=0
            do is=1,solvent(1)%nspec
                do io=1,grid%no
                    do iz=1,grid%nz
                        do iy=1,grid%ny
                            do ix=1,grid%nx
                                i=i+1
                                solvent(is)%density(ix,iy,iz,io) = lbfgsb%x(i)
                            end do
                        end do
                    end do
                end do
            end do

            call energy_and_gradient(f, df) ! f and df are intent(out) of energy_and_gradient. see below for explanations of isave(34)

            !   Energy_and_gradient does not change solvent%density but updates df
            !   Passons df au format lbfgs%g one-column
            i=0
            do is=1,solvent(1)%nspec
                do io=1,grid%no
                    do iz=1,grid%nz
                        do iy=1,grid%ny
                            do ix=1,grid%nx
                                i=i+1
                                lbfgsb%g(i)=df(ix,iy,iz,io,is)
                            end do
                        end do
                    end do
                end do
            end do

        end if
    end do

    print*,
    write(*,'(A,F12.2,A)') "optimal value found by L-BFGS-B =", f," kJ/mol"

end subroutine energy_minimization

subroutine energy_minimization

    use precision_kinds, only : dp
    use module_minimizer, only: lbfgsb, init_lbfgsb
    use module_grid, only: grid
    use module_solvent, only: solvent
    use module_energy_and_gradient, only: energy_and_gradient

    implicit none

    real(dp) :: f ! functional to minimize
    real(dp) :: df (grid%nx, grid%ny, grid%nz, grid%no, solvent(1)%nspec )
    integer :: is, io, ix, iy, iz, i

    print*,
    print*,
    print*, "===== Energy minimization ====="

    ! Init minimization process. From allocations to optimizer parameters to guess solution
    call init_lbfgsb

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
            ! Note pour plus tard : on pourrait utiliser la fonction RESHAPE de fortran
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
            !   Note pour plus tard : on pourrait utiliser la fonction PACK de fortran

            i=0
            do is=1,solvent(1)%nspec
                do io=1,grid%no
                    do iz=1,grid%nz
                        do iy=1,grid%ny
                            do ix=1,grid%nx
                                i=i+1
                                lbfgsb%g(i) = df(ix,iy,iz,io,is)
                            end do
                        end do
                    end do
                end do
            end do

        end if
    end do

    print*,
    print*,
    print*, "===== Energy minimization ====="
    print*,
    print*,

end subroutine energy_minimization

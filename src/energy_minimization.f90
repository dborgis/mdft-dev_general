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
    logical, parameter :: lets_minimize_with_bfgs=.FALSE., lets_minimize_with_steepest_descent=.TRUE.

    print*
    print*
    print*, "===== Energy minimization ====="

    if( lets_minimize_with_bfgs) then
      call minimization_using_lbfgs()
    else if(lets_minimize_with_steepest_descent) then
      call minimization_using_steepest_descent()
    end if

contains

subroutine minimization_using_steepest_descent
  use module_grid, only: grid
  use module_input, only: getinput
  implicit none
  integer :: ix,iy,iz,io,itermax
  real(dp),parameter :: stepsize=10.
  real(dp) :: fold
  real(dp), parameter :: factr=0.0001_dp/epsilon(1._dp)
  itermax = getinput%int("maximum_iteration_nbr", defaultvalue=50, assert=">0")
  i=0
  f=0._dp
  do while(i<itermax)
    print*,""
    print*
    print*
    print*,"ITERATION", i
    fold = f
    call energy_and_gradient(f, df)
    solvent(1)%xi = solvent(1)%xi - stepsize*df(:,:,:,:,1)
    i=i+1
    if( abs(f-fold) > factr ) exit
  end do
end subroutine minimization_using_steepest_descent


subroutine minimization_using_lbfgs
    ! Init minimization process. From allocations to optimizer parameters to guess solution
    call lbfgsb%init

    do while( (lbfgsb%task(1:2).eq.'FG'.or.lbfgsb%task.eq.'NEW_X'.or.lbfgsb%task.eq.'START') &
      .and. (lbfgsb%isave(30)<lbfgsb%itermax) )

time(:)=0
call cpu_time(time(1))

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
            print*, "Timing of whole cycle (evaluation + bfgs stuff): ", time(5)-time(1)
            print*
            print*, "                  -----------------"
            print*
            print*
            print*
        end if
    end do
end subroutine minimization_using_lbfgs

end subroutine energy_minimization

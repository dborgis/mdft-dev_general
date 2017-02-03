subroutine energy_minimization

  use precision_kinds, only : dp
  use module_grid, only: grid
  use module_solvent, only: solvent
  use module_input, only: getinput
  use module_energy_and_gradient, only: energy_and_gradient
  use module_lbfgs_nocedal_mdft, only: lbfgsb

  implicit none

  real(dp) :: f ! functional to minimize
  real(dp) :: df (grid%no, grid%nx, grid%ny, grid%nz, solvent(1)%nspec )
  real :: time(1:10)
  integer :: i
  character(80) :: minimizerName

  
  minimizerName = getinput%char( "minimizer", defaultvalue = "lbfgs", validValues = ["sd       ",&
                                                                                     "lbfgs    ",&
                                                                                     "benchmark"] )
  
  ! keep the old tag "minimize_with_steepest_descent" valid
  if( getinput%log( "minimize_with_steepest_descent", defaultvalue=.false. ) ) then
      minimizerName = "sd"
  end if

  select case( trim(adjustl(minimizerName)) )
  case( "sd" )
      print*, "===== Functional minimization by steepest descent ====="
      call minimization_using_steepest_descent()
  case( "benchmark" )
      print*, "===== Functional minimization canceled. We're benchmarking MDFT. Loop and don't minimize ====="
      call minimization_using_benchmark()
  case default
      print*, "===== Functional minimization by L-BFGS-B ====="
      call minimization_using_lbfgs()
  end select
      
contains

subroutine minimization_using_benchmark()
    ! In this mode, we don't try to minimize our functional.
    ! It is dedicated to making itermax evaluations of the functional.
    ! As such, this mode is intended to benchmarking.
    use module_input, only: getinput
    implicit none
    integer :: itermax, i
    itermax = getinput%int( "maximum_iteration_nbr", defaultvalue=30, assert=">0" )
    do i = 1, itermax
        call energy_and_gradient( f, df )
    end do
end subroutine

  subroutine minimization_using_steepest_descent
    use module_grid, only: grid
    use module_input, only: getinput
    implicit none
    integer :: ix,iy,iz,io,itermax,i,j,k
    real(dp) :: stepsize,stepsize_giving_minimum_f
    real(dp) :: fold,fmin
    real(dp), parameter :: factr=0.00001_dp
    logical :: ich_continue
    integer :: stepsize_n(0:200)
    stepsize_n=0.5
    itermax = getinput%int("maximum_iteration_nbr", defaultvalue=50, assert=">0")
    i=0
    f=0._dp
    fold=huge(1._dp)
    stepsize=0.1
    j=1
    open(12,file="output/iterate.dat")
    do while(i<itermax)
      print*,""
      print*
      print*
      print*,"ITERATION", j
      if(i>0) fold=f
      call energy_and_gradient(f,df)
      fmin=f
      ich_continue=.true.
      k=0
      do while(ich_continue)
        if(k==0.and.i==0) then
          stepsize=0.5
        else if(k==0.and.i>0) then
          stepsize=stepsize_n(i-1)
        end if
        solvent(1)%xi=solvent(1)%xi-stepsize*df(:,:,:,:,1)
        call energy_and_gradient(f)
        solvent(1)%xi=solvent(1)%xi+stepsize*df(:,:,:,:,1)
        if(f<fmin) then
          stepsize_giving_minimum_f=stepsize
          stepsize=stepsize*5.
          fmin=f
          ich_continue=.true.
        else
          if(k==0) then
            ich_continue=.true.
            stepsize=0.5
          else
            ich_continue=.false.
          end if
        end if
        stepsize_n(i)=stepsize_giving_minimum_f
        PRINT*,"AT ITERATION ",i,"BEST STEPSIZE TO DATE=",stepsize_giving_minimum_f,"F=",f,"FMIN=",fmin
        k=k+1
      end do
      PRINT*,"AT ITERATION ",i,"I WILL USE STEPSIZE",stepsize_giving_minimum_f
      PRINT*
      PRINT*
      solvent(1)%xi=solvent(1)%xi-stepsize_giving_minimum_f*df(:,:,:,:,1)
      i=i+1
      j=j+1
      write(12,*) i,fmin
      if( abs(fmin-fold)/max(abs(fmin),abs(fold),1._dp) < factr ) exit
    end do
    close(12)
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

    !   if (lbfgsb%task(1:2) == "FG") then
    !     print*, "Timing of whole cycle (evaluation + bfgs stuff): ", time(5)-time(1)
    !     print*
    !     print*, "                  -----------------"
    !     print*
    !     print*
    !     print*
    !   end if
    end do
  end subroutine minimization_using_lbfgs

end subroutine energy_minimization

! Minimize FF by finding the ground state of CGVECT array
! We use the Low Memory BFGS algorithm, provided by the NLOpt library.
subroutine find_equilibrium_density

  use precision_kinds,only : i2b, dp
  use input, only: input_dp, input_int, input_char
  use minimizer, only: ncg, cg_vect_new, FF, dF_new, minimizer_type, epsg, iter

  implicit none

  external myfunc

  real(dp):: time1, time2
  real(dp):: energy_before, dfoverf, FFdiff
  integer*8 :: opt
  integer :: xsize, ires, M, maxeval, algorithm
  real(dp) :: minf ! minimum value

!  include 'nlopt.f'

  iter = 0

  !minimizer_type = input_char("minimizer")
!  minimizer_type = "bfgs" ! TODO input_char can't now deal with default values for characters
!  select case (minimizer_type(1:3))
!  case ("bfg")
!    algorithm = NLOPT_LD_LBFGS
!  case ("var")
!    algorithm = NLOPT_LD_VAR1
!  case ("ptn")
!    algorithm = NLOPT_LD_TNEWTON_PRECOND_RESTART
!  case default
!    stop "The minimization algorithm you ask for is incorrect."
!  end select
!  print*,"Minimization algorithm: ",trim(adjustl(minimizer_type))
!
!  call nlo_create( opt, algorithm, ncg)
!  if (opt==0) error stop "STOP: problem in find_equilibrium_density.f90 at l.65: while creating the nlopt_opt object"
!
!  call nlo_set_min_objective( ires, opt, myfunc, 0)

  epsg = input_dp("epsg", 0.001_dp ) ! in kJ/mol
!  call nlo_set_ftol_abs( ires, opt, epsg)
! removethis  call nlo_set_ftol_rel( ires, opt, epsg)

!  call nlo_set_lower_bounds1(ires, opt, -huge(1._dp))
!  call nlo_set_upper_bounds1(ires, opt,  huge(1._dp))

  maxeval = input_int("maximum_iteration_nbr", defaultValue=100)
!  call nlo_set_maxeval( ires, opt, maxeval)

!  M = 1 ! default 3
!  call nlo_set_vector_storage( ires, opt, M)

!  call nlo_optimize( ires, opt, cg_vect_new, minf)

BLOCK

    ! ceci est un steepest descent de compet <= selon moi ^^
    use system, only: nb_species
    use module_grid, only: grid
    real(dp) :: FFold, dFF, best_cg_vect_new(GRID%nx,GRID%ny,GRID%nz,GRID%no,nb_species), alpha,&
                                 best_dF_new(GRID%nx,GRID%ny,GRID%nz,GRID%no,nb_species), best_FF
    real(dp), parameter :: hugedp = huge(1._dp)

    dF_new = 0
    alpha = 10._dp
    FF = hugedp
    best_FF = hugedp

    do iter = 1, maxeval
        call energy_and_gradient( iter) ! updates FF and dF_new

!        print*,"=========> at iter",iter,"FF,best_FF,normdF_new",FF,best_FF,norm2(dF_new)
        if( FF < best_FF ) then
!            print*,"FF < best_FF"
            if( abs(FF-best_FF) <= epsg ) then ! convergence reached
                print*,"CONVERGENCE REACHED"
                exit
            end if
            best_cg_vect_new = cg_vect_new
            best_FF = FF
            best_dF_new = dF_new
            cg_vect_new = cg_vect_new - alpha*dF_new
!            print*,"alpha appliqué à cg_vect_new =",alpha
            alpha = alpha*1.1
!            print*,"alpha futur",alpha
        else if( FF >= best_FF) then
!            print*,"FF >= best_FF"
!            print*,"alpha that should have been applied =",alpha
            alpha = alpha/1.27
            cg_vect_new = best_cg_vect_new - alpha*best_dF_new
!            print*,"alpha futur =",alpha
        end if
    end do

    if( iter == maxeval+1 ) print*,"MAXIMUM ITERATION REACHED. CONVERGENCE NOT REACHED."

END BLOCK
!!!!!!!!!!!!


!  select case(ires)
!  case(1)
!    print*,"NLOPT success: Generic success return value."
!  case(2)
!    print*,"NLOPT success: Optimization stopped because stopval (above) was reached."
!  case(3)
!    print*,"NLOPT success: Optimization stopped because ftol_rel or ftol_abs (above) was reached."
!  case(4)
!    print*,"NLOPT success: Optimization stopped because xtol_rel or xtol_abs (above) was reached."
!  case(5)
!    print*,"NLOPT success: Optimization stopped because maxeval (above) was reached."
!  case(6)
!    print*,"NLOPT success: Optimization stopped because maxtime (above) was reached."
!  case(-1)
!    print*,"NLOPT failure: Generic failure code."
!  case(-2)
!    print*,"NLOPT failure: Invalid arguments (e.g. lower bounds are bigger than upper bounds, an unknown algorithm was specified, etcetera)."
!  case(-3)
!    print*,"NLOPT failure: Ran out of memory."
!  case(-4)
!    print*,"NLOPT failure: Halted because roundoff errors limited progress. (In this case, the optimization still typically returns a useful result.)"
!  case(-5)
!    print*,"NLOPT failure: Halted because of a forced termination: the user called nlopt_force_stop(opt) on the optimization’s nlopt_opt object opt from the user’s objective function or constraints."
!  case default
!    error stop "NLOPT failure: dont understand the return value!"
!  end select
!
!  call nlo_destroy(opt)

  print*,
  write(*,'(A,F12.2,A)') "FF before correction", FF," kJ/mol"

end subroutine find_equilibrium_density

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!subroutine myfunc(result, n, x, grad, need_gradient, f_data)
!
!  use precision_kinds, only: dp
!  use minimizer      , only: cg_vect_new, iter, dF_new, FF
!
!  implicit none
!
!  integer, intent(in) :: n ! number of parameters
!  integer :: need_gradient ! different from 0 if one wants to use gradient-based algorithms
!  real(dp), intent(out) :: result, f_data ! the return value of our function we want to optimize
!  real(dp), intent(inout) :: x(n), grad(n)
!
!  iter = iter + 1
!  cg_vect_new = x
!  dF_new = 0._dp
!
!  PRINT*,"====================================> BEFOR ITERATION",iter
!  call energy_and_gradient(iter)
!
!  if( need_gradient/=0 ) then
!    grad = dF_new
!  else
!    stop "don't you need the gradient? something's wrong in myfunc"
!  end if
!
!  result = FF
!end subroutine myfunc

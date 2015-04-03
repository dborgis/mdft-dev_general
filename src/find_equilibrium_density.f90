! Minimize FF by finding the ground state of CGVECT array
! We use the Low Memory BFGS algorithm, provided by the NLOpt library.
subroutine find_equilibrium_density

  use precision_kinds,only : i2b, dp
  use input, only: input_dp, input_int, input_char
  use minimizer, only: ncg, cg_vect, FF, dF, minimizer_type, epsg, iter

  implicit none

  external myfunc

  real(dp):: time1, time2
  real(dp):: energy_before, dfoverf, FFdiff
  integer*8 :: opt
  integer :: xsize, ires, M, maxeval, algorithm
  real(dp) :: minf ! minimum value

  include 'nlopt.f'

  iter = 0

  !minimizer_type = input_char("minimizer")
  minimizer_type = "bfgs" ! TODO input_char can't now deal with default values
  select case (minimizer_type(1:3))
  case ("bfg")
    algorithm = NLOPT_LD_LBFGS
  case ("var")
    algorithm = NLOPT_LD_VAR1
  case ("ptn")
    algorithm = NLOPT_LD_TNEWTON_PRECOND_RESTART
  case default
    stop "The minimization algorithm you ask for is incorrect."
  end select

  call nlo_create( opt, algorithm, ncg)
  if (opt==0) stop "STOP: problem in find_equilibrium_density.f90 at l.65: while creating the nlopt_opt object"

  call nlo_set_min_objective( ires, opt, myfunc, 0)

  epsg = input_dp("epsg", 0.0001_dp ) ! in kJ/mol
  call nlo_set_ftol_abs( ires, opt, epsg)

  maxeval = input_int("maximum_iteration_nbr", defaultValue=30)
  call nlo_set_maxeval( ires, opt, maxeval)

  M = 3
  call nlo_set_vector_storage( ires, opt, M)

  call nlo_optimize( ires, opt, cg_vect, minf)
  if (ires<0) then
    write(*,*) 'minimization failed!'
  else
    write(*,*) 'minimization succeeded'
  end if

  call nlo_destroy(opt)

  print*,
  write(*,'(A,F12.2,A)') "FF before correction", FF," kJ/mol"

end subroutine find_equilibrium_density

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine myfunc(result, n, x, grad, need_gradient, f_data)

  use precision_kinds, only: dp
  use minimizer      , only: cg_vect, iter, dF, FF

  implicit none

  integer, intent(in) :: n ! number of parameters
  integer :: need_gradient ! different from 0 if one wants to use gradient-based algorithms
  real(dp), intent(out) :: result, f_data ! the return value of our function we want to optimize
  real(dp), intent(inout) :: x(n), grad(n)

  iter = iter + 1
  cg_vect = x

  call energy_and_gradient(iter)

  if( need_gradient/=0 ) then
    grad = dF
  else
    stop "don't you need the gradient? something's wrong in myfunc"
  end if

  result = FF

end subroutine myfunc

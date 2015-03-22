! This SUBROUTINE should minimize FF by finding the ground state of CGVECT array
! We use the Low Memory BFGS algorithm, provided by the NLOpt library.

SUBROUTINE find_equilibrium_density

  USE precision_kinds,only : i2b, dp
  use input, only: input_dp, input_int
  use minimizer, only: ncg, cg_vect, FF, dF, minimizer_type, epsg, iter

  IMPLICIT NONE
  external myfunc
  real(dp):: time1, time2
  real(dp):: energy_before, dfoverf, FFdiff
  integer*8 :: opt
  integer :: xsize, ires, M, maxeval
  real(dp) :: minf ! minimum value
  include 'nlopt.f'
  maxeval = input_int("maximum_iteration_nbr")
  epsg = input_dp("epsg", 0.0001_dp ) ! in kJ/mol
  iter = 0
  opt = 0
  if (minimizer_type(1:5)=='bfgs') then
    call nlo_create(opt, NLOPT_LD_LBFGS, ncg) ! Low-storage BFGS
  else if (minimizer_type(1:3)=='var') then
    call nlo_create(opt, NLOPT_LD_VAR1, ncg) ! Shifted limited-memory variable-metric
  else if (minimizer_type(1:3)=='ptn') then
    call nlo_create(opt, NLOPT_LD_TNEWTON_PRECOND_RESTART, ncg) ! Preconditioned truncated Newton
  else
    stop "You ask for a minimizer that is not implemented."
  end if
  if (opt==0) stop "STOP: problem in find_equilibrium_density.f90 at l.65: while creating the nlopt_opt object"
  call nlo_set_min_objective( ires, opt, myfunc, 0)
  call nlo_set_ftol_abs( ires, opt, epsg)
  call nlo_set_maxeval( ires, opt, maxeval)
  call nlo_optimize( ires, opt, cg_vect, minf)
  if (ires<0) then
    write(*,*) 'minimization failed!'
  else
    write(*,*) 'minimization succeeded'
  end if
  call nlo_destroy(opt)
  print*,
  write(*,'(A,F12.2,A)') "FF before correction", FF," kJ/mol"

END SUBROUTINE find_equilibrium_density

subroutine myfunc(result, n, x, grad, need_gradient, f_data)
  use precision_kinds, only: dp
  use minimizer, only: cg_vect, iter, dF, FF
  implicit none
  integer, intent(in)          :: n ! number of parameters
  integer          :: need_gradient ! different from 0 if one wants to use gradient-based algorithms
  real(dp), intent(out) :: result, f_data ! the return value of our function we want to optimize
  real(dp), intent(inout) :: x(n), grad(n)
  iter = iter +1
  cg_vect = x
  call energy_and_gradient(iter)
  if(need_gradient/=0) then
    grad = dF
  end if
  result = FF
  end subroutine myfunc

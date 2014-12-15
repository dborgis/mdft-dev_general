module time
! src: http://stackoverflow.com/questions/5083051/timing-a-fortran-multithreaded-program
contains

  subroutine tick(t)
    implicit none
    integer, intent(out) :: t
    call system_clock(t)
  end subroutine tick

  real function tock(t)
    integer, intent(in) :: t
    integer :: now, clock_rate
    call system_clock(now, clock_rate)
    tock = real(now - t)/real(clock_rate)
  end function tock

end module time


! exemple on using this:
! call tick(calc)
! ! do big calculation
! calctime = tock(calc)
! print *,'Timing summary'
! print *,'Calc: ', calctime

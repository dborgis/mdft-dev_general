!> double is the module defining precision variables.
!! This allows architecture independant programming.
module precision_kinds
  use iso_c_binding, only: C_FLOAT, C_DOUBLE
  implicit none
  integer, parameter :: dp = C_DOUBLE   ! usual double precision float
  integer, parameter :: sp = C_FLOAT       ! usual simple precision float

  integer, parameter :: i2b = KIND(1)      ! usual integer
  integer, parameter :: i4b = 2_i2b * i2b  ! usual double precision integer
end module precision_kinds

! defining the precision kind this way we say to the computer : use the default single precision and double precision.
! remember that this default value may change depending on the architecture of your CPU.
! but one might one day wish to change that. It would be such a mess do change every number in the program !
! One would just have to change here the definition of single and double precision.
! For a very nice lesson about this : http://www.owlnet.rice.edu/~ceng303/manuals/fortran/FOR2_5.html

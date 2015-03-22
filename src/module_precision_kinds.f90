!> double is the module defining precision variables.
!! This allows architecture independant programming.
MODULE precision_kinds

    IMPLICIT NONE
        
    INTEGER( KIND=KIND(1) ), PARAMETER :: i2b = KIND(1) ! usual integer 
    INTEGER(i2b), PARAMETER :: dp = KIND(0.0d0)         ! usual double precision float
    INTEGER(i2b), PARAMETER :: sp = KIND(0.0)           ! usual simple precision float
    INTEGER(i2b), PARAMETER :: i4b = 2_i2b * i2b        ! usual double precision integer

END MODULE precision_kinds
! defining the precision kind this way we say to the computer : use the default single precision and double precision.
! remember that this default value may change depending on the architecture of your CPU.
! but one might one day wish to change that. It would be such a mess do change every number in the program !
! One would just have to change here the definition of single and double precision.
! For a very nice lesson about this : http://www.owlnet.rice.edu/~ceng303/manuals/fortran/FOR2_5.html

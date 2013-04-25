! Module for angular grid and Gauss-Legendre integration.

module quadrature

use precision_kinds , only : dp , i2b ! the module defining precision kinds

implicit none


real(dp), allocatable , dimension ( : , : ) :: w_legendre , x_legendre ! w(i,L) (weights) and x(i,L) (roots) for order L integration

real(dp), allocatable , dimension ( : ) :: Omx , Omy , Omz , weight  ! unit vector for orientation OMEGA and associated weight

real (kind = dp ), allocatable, dimension ( : )  :: weight_psi , x_psi

real(dp) , allocatable, dimension(:) :: x_leb, y_leb , z_leb , weight_leb 

integer(i2b) :: sym_order

contains

  subroutine deallocate_everything_gauss_legendre

    implicit none

    if ( allocated (weight_psi ) ) deallocate ( weight_psi)

    if ( allocated (x_psi ) ) deallocate ( x_psi)

    if ( allocated ( w_legendre ) ) deallocate ( w_legendre )

    if ( allocated ( x_legendre ) ) deallocate ( x_legendre ) 

    if ( allocated ( Omx ) ) deallocate ( Omx )

    if ( allocated ( Omy ) ) deallocate ( Omy )

    if ( allocated ( Omz ) ) deallocate ( Omz )

    if ( allocated ( x_leb ) ) deallocate ( y_leb)

    if ( allocated ( y_leb ) ) deallocate ( x_leb)

    if ( allocated ( z_leb ) ) deallocate ( z_leb)

    if ( allocated ( weight ) ) deallocate ( weight)

    if ( allocated ( weight_leb ) ) deallocate ( weight_leb)

  end subroutine deallocate_everything_gauss_legendre

end module quadrature

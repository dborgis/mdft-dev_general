! !This SUBROUTINE find the index icg in the minimizer for a vector ( s, i, j, k, o, p)
! integer (i2b) function convert_coordinate_into_icg ( s, i, j, k, io)
!     USE precision_kinds,only : i2b , dp
!     use quadrature, only: angGrid, molRotGrid
!     use module_grid, only: grid
!     IMPLICIT NONE
!     integer (i2b) , intent (in) :: i , j , k , io, s
!     INTEGER(i2b) :: nx, ny, nz
!     nx= grid%n_nodes(1)
!     ny= grid%n_nodes(2)
!     nz= grid%n_nodes(3)
!
! convert_coordinate_into_icg = (s-1)*nx*ny*nz*angGrid%n_angles*molRotGrid%n_angles &
!     + ( i-1) *ny*nz*angGrid%n_angles*molRotGrid%n_angles+&
!     &( j-1)*nz*angGrid%n_angles*molRotGrid%n_angles + (k -1)*angGrid%n_angles*molRotGrid%n_angles + (o-1)*molRotGrid%n_angles +p
! end function

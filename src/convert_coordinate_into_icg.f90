!This SUBROUTINE find the index icg in the minimizer for a vector ( species, i, j, k, o, p)
integer (i2b) function convert_coordinate_into_icg ( species, i, j, k, o, p)
USE precision_kinds , only : i2b , dp
use system , only : nfft1 , nfft2 ,nfft3
use quadrature, only: angGrid, molRotGrid
IMPLICIT NONE
integer (i2b) , intent (in) :: i , j , k , o , p, species
convert_coordinate_into_icg = (species-1)*nfft1*nfft2*nfft3*angGrid%n_angles*molRotGrid%n_angles &
    + ( i-1) *nfft2*nfft3*angGrid%n_angles*molRotGrid%n_angles+&
    &( j-1)*nfft3*angGrid%n_angles*molRotGrid%n_angles + (k -1)*angGrid%n_angles*molRotGrid%n_angles + (o-1)*molRotGrid%n_angles +p
end function

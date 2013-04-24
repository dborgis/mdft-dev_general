!This subroutine find the index icg in the minimizer for a vector ( species, i, j, k, o, p)
integer (kind=i2b) function convert_coordinate_into_icg ( species, i, j, k, o, p)

use precision_kinds , only : i2b , dp
use system , only : nfft1 , nfft2 ,nfft3 , nb_omega , nb_psi


implicit none

integer (kind=i2b) , intent (in) :: i , j , k , o , p, species


convert_coordinate_into_icg = (species-1)*nfft1*nfft2*nfft3*nb_omega*nb_psi + ( i-1) *nfft2*nfft3*nb_omega*nb_psi + &
&( j-1)*nfft3*nb_omega*nb_psi + (k -1)*nb_omega*nb_psi + (o-1)*nb_psi + p

end function

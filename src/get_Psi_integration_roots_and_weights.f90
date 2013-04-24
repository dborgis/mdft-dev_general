subroutine get_psi_integration_roots_and_weights



use system, only : nb_psi 

use constants , only : twopi, pi

use precision_kinds , only : i2b , dp

use input , only : input_line

use quadrature, only : x_psi , weight_psi, sym_order



implicit none




integer (kind = i2b ) :: n 




allocate (weight_psi ( nb_psi ) )
allocate (x_psi ( nb_psi ) )

do n=1 , nb_psi
   
  weight_psi(n) =twopi/real(nb_psi*sym_order,kind=dp) !Remplacer 1 par l' ordre de Symetrie de la molécule.
  x_psi(n) = twopi*real((n-1),kind=dp)/real((nb_psi*sym_order),kind=dp)  

end do



end subroutine


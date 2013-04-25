!> Gets the final density from the last minimizer step.

subroutine get_final_polarization ( Px , Py , Pz )




use precision_kinds , only: dp , i2b

use system , only : nfft1 , nfft2 , nfft3 , nb_omega , nb_species

use constants , only : fourpi

! cg contains everything related to the minimizer
use cg , only : CG_vect
! cg_vect = density for at each grid point for each angle and each species

use quadrature , only : weight , Omx , Omy , Omz




implicit none




real(dp), dimension ( nfft1 , nfft2 , nfft3 , nb_species ) , intent(out) :: Px , Py , Pz ! equilibrium polarization(position)

integer(i2b):: i , j , k , omega , icg , species ! dummy

real(dp):: rho_over_fourpi !> = CG_vect(i)**2/fourpi

real(dp):: local_Px , local_Py , local_Pz ! dummy for speeding loops

real(dp), allocatable , dimension ( : ) :: weight_omx , weight_omy , weight_omz ! weight(:)*omx(:) for speeding up





! init outputs

Px = 0.0_dp

Py = 0.0_dp

Pz = 0.0_dp

allocate ( weight_omx ( nb_omega ) )

allocate ( weight_omy ( nb_omega ) )

allocate ( weight_omz ( nb_omega ) )

weight_omx = weight * Omx

weight_omy = weight * Omy

weight_omz = weight * Omz






! read cg_vect and get density and polarization from it

! note again that rho is the density per angle so that 

icg = 0

do species = 1 , nb_species

  do i = 1 , nfft1

    do j = 1 , nfft2

      do k = 1 , nfft3

        local_Px = 0.0_dp

        local_Py = 0.0_dp

        local_Pz = 0.0_dp

        do omega = 1 , nb_omega

          icg = icg + 1

          rho_over_fourpi = cg_vect ( icg ) ** 2 / fourpi

          local_Px = local_Px + weight_Omx ( omega ) * rho_over_fourpi

          local_Py = local_Py + weight_Omy ( omega ) * rho_over_fourpi

          local_Pz = local_Pz + weight_Omz ( omega ) * rho_over_fourpi

        end do

        Px ( i , j , k , species  ) = local_Px

        Py ( i , j , k , species  ) = local_Py

        Pz ( i , j , k , species  ) = local_Pz

      end do

    end do

  end do

end do



! deallocate

deallocate ( weight_omx )

deallocate ( weight_omy )

deallocate ( weight_omz )



end subroutine get_final_polarization

! this subroutine compute the external part of the free energy functional

subroutine energy_external


use precision_kinds , only : dp , i2b

use system , only : DeltaV , rho_0 , nfft1 , nfft2 , nfft3 , nb_omega , nb_species , rho_0_multispec , nb_psi

use input , only : input_line


use quadrature , only : weight, weight_psi

use cg , only : CG_vect , FF , dF

use external_potential , only : Vext_total, Vext_q


implicit none

real(dp):: Fext , Fext_q! external part of the total energy

integer(i2b):: icg , i , j , k , o , p! dummy for loops

real(dp):: psi !dummy

real(dp):: rho ! dummy = psi**2

real(dp):: wdfve 

real(dp):: time0,time1! timers

integer(i2b):: species  ! dummy between 1 and nb_species





call cpu_time ( time0 )

! init external potential contribution and its gradient

Fext = 0.0_dp ! scalar

Fext_q = 0.0_dp

! the following is a dummy variable to improve loop speed by decreasing the total number of calculations



! Compute Fext and dFext

! F_{ext}[\rho(\vec{r},\vec{\Omega})]=\int d \vec{r} d \vec{\Omega} V_{ext}(\vec{r},\vec{\Omega})\rho(\vec{r},\vec{\Omega})

icg = 0

do species = 1 , nb_species

  do i = 1 , nfft1

    do j = 1 , nfft2

      do k = 1 , nfft3

        do o = 1 , nb_omega

          do p=1 , nb_psi            

            icg = icg + 1

            psi = cg_vect ( icg )

            rho = psi ** 2

            wdfve = weight ( o )*weight_psi(p) * DeltaV * rho_0_multispec ( species ) * Vext_total ( i , j , k , o , p , species )

            Fext = Fext + rho * wdfve


            dF ( icg ) = dF ( icg ) + 2.0_dp * psi * wdfve

          end do

        end do

      end do

    end do

  end do

end do



! conclude

FF = FF + Fext ! scalar


! stop timer

call cpu_time(time1)

! tell user

write (*,*) 'external    = ' , Fext , 'computed in (sec)' , time1 - time0



end subroutine energy_external

! Compute total energy and gradients using direct correlation functions c_s_hs of a hard sphere fluid

subroutine energy_cs_hard_sphere


use precision_kinds, only: i2b,dp

use system, only: nfft1 , nfft2 , nfft3 , Lx , Ly , Lz , nb_omega , c_s_hs , kBT , nb_k , delta_k , deltaV , rho_0_multispec ,&
                  nb_species, nb_psi

use quadrature, only: weight, weight_psi,sym_order

use cg, only: cg_vect , FF , dF

use constants, only: fourpi , pi , twopi

use fft, only: in_forward,in_backward,out_forward,out_backward,plan_forward,plan_backward , norm_k



implicit none


integer(kind=i2b) :: i, j, k, l, m, n, o, icg, species,p !> Dummy

integer(kind=i2b) :: k_index

real(kind=dp) :: Nk !> Total number of k points = nfft1*nfft2*nfft3

real(kind=dp) :: Fint !> Internal part of the free energy

real(kind=dp) :: Vint !> Dummy for calculation of Vint

real(kind=dp) :: fact !> facteur d'integration

real(kind=dp) :: psi ! Dummy

real(kind=dp), allocatable, dimension(:,:,:) :: Delta_rho

complex(kind=dp), allocatable, dimension(:,:,:) :: rho_k , Vpolarization_k

real(kind=dp) :: time1, time0

real(kind=dp) , allocatable , dimension ( : , : , : ) :: Vpolarization



!> Initiate

call cpu_time ( time0 )

Nk = real(nfft1*nfft2*nfft3,kind=dp) ! nombre de points k

allocate ( Delta_rho ( nfft1 , nfft2 , nfft3 ) )

Delta_rho = 0.0_dp ! density



!> Put density of last minimization step in delta_rho

icg=0
do i=1,nfft1
  do j=1,nfft2
    do k=1,nfft3
      do o = 1, nb_omega
        do p=1, nb_psi
        icg=icg+1

        Delta_rho(i,j,k) = Delta_rho(i,j,k) + weight(o) * cg_vect(icg)**2*weight_psi(p)
       end do
      end do
    end do
  end do
end do
Delta_rho = Delta_rho-(twopi*fourpi)/real(sym_order,kind=dp)


!> Next FFT sequences can be done on multiple threads
!> Compute rho in k-space

in_forward = Delta_rho

call dfftw_execute ( plan_forward )

allocate ( rho_k ( nfft1 / 2 + 1 , nfft2 , nfft3 ) )

rho_k = out_forward




! Compute polarisation in k-space

allocate ( Vpolarization_k ( nfft1 / 2 + 1 , nfft2 , nfft3 ) )

! V(k)=cs(k)*rho(k)

do n = 1 , nfft3

  do m = 1 , nfft2

    do l = 1 , nfft1 / 2 + 1

      k_index = int ( norm_k ( l , m , n ) / delta_k ) + 1

      ! Here it happens that k_index gets higher than the highest c_k index.
      ! In this case one imposes k_index = k_index_max

      if ( k_index > nb_k ) k_index = nb_k




      ! V(k)=cs(k)*rho(k)

      Vpolarization_k ( l , m , n ) = rho_k ( l , m , n ) * c_s_hs ( k_index )

    end do

  end do

end do

! since rho(k) is now useless, deallocate associated array

deallocate ( rho_k )





! FFT-1

in_backward = Vpolarization_k

deallocate (Vpolarization_k)

call dfftw_execute (plan_backward)

allocate ( Vpolarization ( nfft1 , nfft2 , nfft3 ) )

Vpolarization = out_backward / Nk 




! compute excess energy and its gradient

Fint = 0.0_dp ! excess energy

icg = 0 ! index of cg_vect

do species = 1 , nb_species

  fact = DeltaV * rho_0_multispec ( species ) !> facteur d'integration

  do i = 1 , nfft1

    do j = 1 , nfft2

      do k = 1 , nfft3

        Vint   = -kBT * rho_0_multispec ( species ) * Vpolarization(i,j,k)

        do o = 1 , nb_omega

          do p=1, nb_psi

          icg = icg + 1

          psi = CG_vect ( icg )

          Fint   = Fint   + weight(o)*weight_psi(p) * fact * 0.5_dp * ( psi ** 2 - 1.0_dp) * Vint

!         dF (icg) = dF ( icg ) + 2.0_dp * psi * weight(o) * fact * Vint ! in case of bridge calculation, one deduces the pair contribution of hard spheres + => - and FF=FF-Fint
          dF (icg) = dF ( icg ) - 2.0_dp * psi * weight(o) *weight_psi(p)* fact * Vint

          end do

        end do

      end do

    end do

  end do

end do ! species

deallocate(Vpolarization)



! conclude

FF = FF - Fint

call cpu_time(time1)

write(*,*) 'Fexc c_hs   = ' , Fint , 'computed in (sec)' , time1 - time0


 

end subroutine energy_cs_hard_sphere

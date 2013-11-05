!> Compute total energy and gradients using direct correlation functions c_s
subroutine cs_from_dcf
use precision_kinds, only: i2b,dp
use system, only: nfft1 , nfft2 , nfft3 , Lx , Ly , Lz , c_s , kBT , nb_k , delta_k , deltaV , rho_0_multispec ,&
                  nb_species
use quadrature, only: sym_order, angGrid, molRotGrid
use cg, only: cg_vect , FF , dF
use constants, only: fourpi , pi , twopi
use fft, only: fftw3, norm_k
implicit none
integer(i2b) :: i, j, k, l, m, n, o, p , icg, species !> Dummy
integer(i2b) :: k_index
real(dp) :: Nk !> Total number of k points = nfft1*nfft2*nfft3
real(dp) :: Fint !> Internal part of the free energy
real(dp) :: Vint !> Dummy for calculation of Vint
real(dp) :: fact !> facteur d'integration
real(dp) :: psi ! Dummy
real(dp), allocatable, dimension(:,:,:) :: Delta_rho
real(dp) :: delta_rho_ijk ! dummy delta_rho(i,j,k) @ given ijk
complex(dp), allocatable, dimension(:,:,:) :: rho_k , Vpair_k
real(dp) :: time1, time0
real(dp) , allocatable , dimension ( : , : , : ) :: Vpair
!> Initiate
call cpu_time ( time0 )
Nk = real(nfft1*nfft2*nfft3,dp) ! nombre de points k
allocate ( Delta_rho ( nfft1 , nfft2 , nfft3 ) )
Delta_rho = 0.0_dp ! density
!> Put density of last minimization step in delta_rho
icg=0
do i=1,nfft1
  do j=1,nfft2
    do k=1,nfft3
      delta_rho_ijk = 0.0_dp
      do o = 1, angGrid%n_angles
        do p=1 , molRotGrid%n_angles
          icg=icg+1
          delta_rho_ijk = delta_rho_ijk + angGrid%weight(o) * cg_vect(icg)**2*molRotGrid%weight(p)
        end do
      end do
      delta_rho ( i , j , k ) = delta_rho_ijk
    end do
  end do
end do
Delta_rho = Delta_rho-real(2.0_dp*twopi**2/sym_order, dp)
!> Next FFT sequences can be done on multiple threads
!> Compute rho in k-space
fftw3%in_forward = Delta_rho
call dfftw_execute ( fftw3%plan_forward )
allocate ( rho_k ( nfft1 / 2 + 1 , nfft2 , nfft3 ) )
rho_k = fftw3%out_forward

! Compute polarisation in k-space
allocate ( Vpair_k ( nfft1 / 2 + 1 , nfft2 , nfft3 ) )
! V(k)=cs(k)*rho(k)
do n = 1 , nfft3
  do m = 1 , nfft2
    do l = 1 , nfft1 / 2 + 1
      k_index = int ( norm_k ( l , m , n ) / delta_k ) + 1
      ! Here it happens that k_index gets higher than the highest c_k index.
      ! In this case one imposes k_index = k_index_max
      if ( k_index > nb_k ) k_index = nb_k
      ! V(k)=cs(k)*rho(k)
      Vpair_k ( l , m , n ) = rho_k ( l , m , n ) * c_s ( k_index )
    end do
  end do
end do
! since rho(k) is now useless, deallocate associated array
deallocate ( rho_k )
! FFT-1
fftw3%in_backward = Vpair_k
deallocate (Vpair_k)
call dfftw_execute (fftw3%plan_backward)
allocate ( Vpair ( nfft1 , nfft2 , nfft3 ) )
Vpair = fftw3%out_backward / Nk
! compute excess energy and its gradient
Fint = 0.0_dp ! excess energy
icg = 0 ! index of cg_vect
do species = 1 , nb_species
  fact = DeltaV * rho_0_multispec ( species ) !> facteur d'integration
  do i = 1 , nfft1
    do j = 1 , nfft2
      do k = 1 , nfft3
        Vint   = -kBT * rho_0_multispec ( species ) * Vpair(i,j,k)
        do o = 1 , angGrid%n_angles
          do p=1, molRotGrid%n_angles
          icg = icg + 1
          psi = CG_vect ( icg )
          Fint   = Fint   + angGrid%weight(o) * fact * 0.5_dp * ( psi ** 2 - 1.0_dp) * Vint*molRotGrid%weight(p)
          dF (icg) = dF ( icg ) + 2.0_dp * psi * angGrid%weight(o) * fact * Vint*molRotGrid%weight(p)
         end do
   
        end do
      end do
    end do
  end do
end do ! species
deallocate(Vpair)
! conclude
FF = FF + Fint
call cpu_time(time1)
write(*,*) 'Fexc(rad)   = ' , Fint , 'computed in (sec)' , time1 - time0
 
end subroutine cs_from_dcf

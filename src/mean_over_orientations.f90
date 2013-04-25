!> This subroutine do the legendre integration over all orientations of any array (density, vext, ...)
subroutine mean_over_orientations(arrayin,arrayout)
use precision_kinds, only: dp,i2b
use system, only: nb_omega,nfft1,nfft2,nfft3, nb_psi
use quadrature, only:weight, weight_psi
implicit none
real(dp), dimension(nfft1,nfft2,nfft3,nb_omega,nb_psi), intent(in) :: arrayin !>@var input array
real(dp), dimension(nfft1,nfft2,nfft3), intent(out) :: arrayout !>@var output array
integer(i2b) :: n , p

arrayout=0.0_dp

do n=1,nb_omega
  do p=1,nb_psi
  arrayout(:,:,:)=arrayout(:,:,:)+arrayin(:,:,:,n,p)*weight(n)*weight_psi(p)
  end do
end do


end subroutine mean_over_orientations

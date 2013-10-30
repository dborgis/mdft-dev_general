!> This subroutine do the legendre integration over all orientations of any array (density, vext, ...)
subroutine mean_over_orientations(arrayin,arrayout)
use precision_kinds, only: dp,i2b
use system, only: nfft1,nfft2,nfft3
use quadrature, only: angGrid, molRotGrid
implicit none
real(dp), dimension(nfft1,nfft2,nfft3,angGrid%n_angles,molRotGrid%n_angles), intent(in) :: arrayin !>@var input array
real(dp), dimension(nfft1,nfft2,nfft3), intent(out) :: arrayout !>@var output array
integer(i2b) :: n , p
arrayout=0.0_dp
do n=1,angGrid%n_angles
    do p=1,molRotGrid%n_angles
        arrayout(:,:,:)=arrayout(:,:,:)+arrayin(:,:,:,n,p)*angGrid%weight(n)*molRotGrid%weight(p)
    end do
end do
end subroutine mean_over_orientations

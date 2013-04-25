!> Compute mean density perpendicular to plan xy
subroutine compute_z_density(array,filename)
use precision_kinds , only: dp,i2b
use system, only: Lz , nfft1 , nfft2 , nfft3 , deltaz
implicit none
real(dp) :: mean_density !>@var mean density in plan xy for each z
real(dp),dimension(nfft1,nfft2,nfft3), intent(in) :: array
integer(i2b) :: k !>@var dummy
real(dp) :: z
 character(50),intent(in):: filename
!> Open file "filename"
open(10,file=filename)
!> Intro file
write(10,*)'#z   <rho>_x,y'
!> Compute mean density over x and y
mean_density=0.0_dp
do k=1,nfft3
  z = real(k-1,dp) * deltaz
  mean_density = sum(array(:,:,k)) / real(nfft1*nfft2,dp)
  write(10,*)z,mean_density
end do
write(*,*)filename,' written'
end subroutine compute_z_density

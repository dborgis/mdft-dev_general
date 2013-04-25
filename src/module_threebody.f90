module threebody

use precision_kinds, only:dp,i2b

implicit none

real(dp), dimension(nfft1,nfft2,nfft3) :: Axx,Ayy,Azz,Axy,Axz,Ayz,Ax,Ay,Az,A0
complex(dp), dimension(nfft1/2+1,nfft2,nfft3) :: Axx_k,Ayy_k,Azz_k,Axy_k,Axz_k,Ayz_k,Ax_k,Ay_k,Az_k,A0_k
real(dp), dimension(nfft1,nfft2,nfft3) :: Gxx,Gyy,Gzz,Gxy,Gxz,Gyz,Gx,Gy,Gz,G0 

end module

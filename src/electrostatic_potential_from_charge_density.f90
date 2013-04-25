! This subroutine uses the Gauss theorem associated to the expression of the derivation in Fourier transforms to get V(r).
! Field E(r)=-grad(V(r))
! Local expression of Gauss theorem : div(E(r))=rho_c(r)/eps0
! => laplacien(V(r))=-rho_c(r)/eps0
! => V(k)=rho_c(k)/(esp0*k^2)
! FFT(V(k)) = V(r)
subroutine electrostatic_potential_from_charge_density
use precision_kinds , only : dp , i2b
use system , only : nfft1 , nfft2 , nfft3 , rho_c
! nfft = number of FFT grid nodes in each direction
! rho_c (nfft1,nfft2,nfft3) = discrete charge density per unit volume
use fft , only : in_forward , in_backward , out_forward , out_backward , plan_forward , plan_backward , norm_k , k2 , k2_nocoef
use constants , only : fourpi , twopi
use external_potential , only : V_c
! V_c = electrostatic potential from charge density and poisson equation
implicit none
complex(dp), dimension ( nfft1 / 2 + 1 , nfft2 , nfft3 ) :: rho_c_k
complex(dp), dimension ( nfft1 / 2 + 1 , nfft2 , nfft3 ) :: V_c_k
integer (i2b) :: i !dummy vraiable
! check if rho_c exists and is allocated
if ( .not. allocated ( rho_c ) ) then
  print *, 'rho_c is not allocated in electrostatic_potential_from_charge_density.f90'
  stop
end if
if (maxval(abs(rho_c)) < tiny(1.0_dp)) then
allocate ( V_c ( nfft1 , nfft2 , nfft3 ) )
v_c = 0.0_dp
return
end if
! FFT of rho_c
in_forward = rho_c
call dfftw_execute ( plan_forward )
rho_c_k = out_forward ! It is verified that at this point, FFT-1(rho_c_k)/ (nfft1*nfft2*nfft3) = rho_c
! FFT(Laplacian(V(r))) = FFT( - 4Pi charge density(r) ) in elecUnits = (ik)^2 V(k) = -4pi rho(k)
! V(k) = 4Pi rho(k) / k^2
where ( k2 /= 0.0_dp )
  V_c_k = rho_c_k * fourpi / k2 ! in electrostatic units : V=-4pi rho
elsewhere
  V_c_k = 0.0_dp ! arbitrary choice equivalent to an infinite dielectric material at infinity
end where
! get real space potential V(r)
if ( .not. allocated ( V_c ) ) allocate ( V_c ( nfft1 , nfft2 , nfft3 ) )
in_backward = V_c_k
call dfftw_execute ( plan_backward )
V_c = out_backward / real ( nfft1 * nfft2 * nfft3 , dp )
open(11,file='output/V_cmax.dat')
do i=1, nfft1
write(11,*), i , V_c(i, nfft2/2, nfft3/2 )
 
end do
close(11)
end subroutine electrostatic_potential_from_charge_density

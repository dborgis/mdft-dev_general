subroutine get_charge_density_k ( Rotxx,Rotxy,Rotxz,Rotyx,Rotyy,Rotyz,Rotzx,Rotzy,Rotzz ) 
!This routine compute : -The solvent molecular charge density, which can be used into Vcoul_from_solvent_charge_density.f90 to 
!evaluate the electrostatic potential.
!                       -The solvent molecular polarization (from Ranieriet al : J. Chem. Phys. 98 (11) 1993) which ca be used
!into energy_polarization_myway.f90 to compute the (multipolar) polarization Free energy.
use precision_kinds, only : i2b, dp
use constants, only : i_complex, twopi, fourpi
use system, only : chg_solv, x_solv, y_solv, z_solv, nfft1, nfft2, nfft3, Lx, Ly, Lz,nb_solvent_sites, id_solv&
, sigma_k,molec_polarx_k, molec_polary_k, molec_polarz_k,nb_species, deltaV,deltax
use external_potential, only : x_charge, y_charge, z_charge, q_charge, nb_of_interpolation
use cg , only : cg_vect
use quadrature, only : Omx , Omy , Omz, angGrid, molRotGrid
use fft , only : kx, ky, kz, k2

implicit none
integer (i2b) :: i, j, k, o , p , n,species!dummy
real(dp), dimension(angGrid%n_angles,molRotGrid%n_angles), intent(in) :: Rotxx,Rotxy,Rotxz,Rotyx,Rotyy,Rotyz,Rotzx,Rotzy,Rotzz
integer (i2b) :: nf1
real (dp) :: xmod, ymod, zmod
real (dp) :: deltaVk, Rc
real (dp), dimension(nfft1,nfft2,nfft3)::molecpolarx,molecpolary,molecpolarz
!            ====================================================
!            !    	Initialization				!
!            !							!
!            ====================================================
Rc=0.5_dp
if (Rc/=0.0_dp) then
print*, 'WARNING: you convolute molecular Charge Density and POLARIZATION with a Gaussian be sure that is what you want!!!'
end if
deltaVk=twopi**3/(Lx*Ly*Lz)
nf1=nfft1/2
allocate(sigma_k(nf1+1, nfft2, nfft3, angGrid%n_angles, molRotGrid%n_angles,nb_species))
sigma_k=(0.0_dp,0.0_dp)
allocate (molec_polarx_k (nf1+1 , nfft2 , nfft3 , angGrid%n_angles , molRotGrid%n_angles,nb_species))
molec_polarx_k= ( 0.0_dp , 0.0_dp )
allocate (molec_polary_k (nf1+1 , nfft2 , nfft3 , angGrid%n_angles , molRotGrid%n_angles,nb_species))
molec_polary_k= ( 0.0_dp , 0.0_dp )
allocate (molec_polarz_k (nf1+1 , nfft2 , nfft3 , angGrid%n_angles , molRotGrid%n_angles,nb_species))
molec_polarz_k= ( 0.0_dp , 0.0_dp )
!            ====================================================
!            !    Compute sigma and molecular polarization	!
!            !							!
!            ====================================================
do species=1,nb_species
do i = 1 , nf1 + 1
  
  do j = 1 , nfft2
    do k = 1 , nfft3
        do o=1, angGrid%n_angles
           do p=1, molRotGrid%n_angles
             do n=1, nb_solvent_sites
             xmod= Rotxx(o,p)*x_solv(n) + Rotxy(o,p)*y_solv(n) + Rotxz(o,p)*z_solv(n)
             ymod= Rotyx(o,p)*x_solv(n) + Rotyy(o,p)*y_solv(n) + Rotyz(o,p)*z_solv(n)   
             zmod= Rotzx(o,p)*x_solv(n) + Rotzy(o,p)*y_solv(n) + Rotzz(o,p)*z_solv(n)  
!print*, xmod,ymod,zmod,chg_solv(id_solv(n))
               if (xmod==Lx) then
                xmod=0.0_dp
               end if
               if (ymod==Ly) then
               ymod=0.0_dp
               end if
               if (zmod==Lz) then
               zmod=0.0_dp
               end if
  
               sigma_k (i , j, k ,o,p,species) = sigma_k (i , j , k ,o,p,species) +chg_solv(id_solv(n))*&
               exp(-i_complex*(kx(i)*xmod+ky(j)*ymod+ kz(k)*zmod ))*exp(-(Rc**2*k2(i,j,k))/2)
               if ((xmod*kx(i)+ymod*ky(j)+zmod*kz(k))==0.0_dp ) then
            
               molec_polarx_k (i,j,k,o,p,species)=molec_polarx_k(i,j,k,o,p,species) + chg_solv(id_solv(n))*xmod
               molec_polary_k (i,j,k,o,p,species)=molec_polary_k(i,j,k,o,p,species) + chg_solv(id_solv(n))*ymod
               molec_polarz_k (i,j,k,o,p,species)=molec_polarz_k(i,j,k,o,p,species) + chg_solv(id_solv(n))*zmod
              else
               molec_polarx_k (i,j,k,o,p,species)=molec_polarx_k(i,j,k,o,p,species) + chg_solv(id_solv(n))*xmod/&
              (xmod*kx(i)+ymod*ky(j)+zmod*kz(k))*(exp(i_complex*(xmod*kx(i)+ymod*ky(j)+zmod*kz(k)))-1)*(-i_complex)&
               *exp(-(Rc**2*k2(i,j,k))/2)
               molec_polary_k (i,j,k,o,p,species)=molec_polary_k(i,j,k,o,p,species) + chg_solv(id_solv(n))*ymod/&
               (xmod*kx(i)+ymod*ky(j)+zmod*kz(k))*(exp(i_complex*(xmod*kx(i)+ymod*ky(j)+zmod*kz(k)))-1)*(-i_complex)&
               *exp(-(Rc**2*k2(i,j,k))/2)
               molec_polarz_k (i,j,k,o,p,species)=molec_polarz_k(i,j,k,o,p,species) + chg_solv(id_solv(n))*zmod/&
               (xmod*kx(i)+ymod*ky(j)+zmod*kz(k))*(exp(i_complex*(xmod*kx(i)+ymod*ky(j)+zmod*kz(k)))-1)*(-i_complex)&
               *exp(-(Rc**2*k2(i,j,k))/2)
               end if
            end do  !n
           end do  !p
        end do  !o
     end do  !k
   end do  !j
end do  !i
end do!species
!in_backward=molec_polarx_k (:,:,:,1,1,1)
!call dfftw_execute(plan_backward)
!molecpolarx=out_backward*deltaVk/(twopi)**3
!in_backward=molec_polary_k (:,:,:,1,1,1)
!call dfftw_execute(plan_backward)
!molecpolary=out_backward*deltaVk/(twopi)**3
!in_backward=molec_polarz_k (:,:,:,1,1,1)
!call dfftw_execute(plan_backward)
!molecpolarz=out_backward*deltaVk/(twopi)**3
!open(11,file='output/mux')
!do i =1,nfft1
!write(11,*) i*deltax, molecpolarx(i,1,1)
!end do
!close(11)
!open(11,file='output/muy')
!do i =1,nfft1
!write(11,*) i*deltax, molecpolary(i,1,1)
!end do
!close(11)
!open(11,file='output/muz')
!do i =1,nfft1
!write(11,*) i*deltax, molecpolarz(i,1,1)
!end do
!close(11)
end subroutine

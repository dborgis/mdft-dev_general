subroutine get_charge_density_k ( Rotxx,Rotxy,Rotxz,Rotyx,Rotyy,Rotyz,Rotzx,Rotzy,Rotzz ) 

!This routine compute : -The solvent molecular charge density, which can be used into Vcoul_from_solvent_charge_density.f90 to 
!evaluate the electrostatic potential.
!                       -The solvent molecular polarization (from Ranieriet al : J. Chem. Phys. 98 (11) 1993) which ca be used
!into energy_polarization_myway.f90 to compute the (multipolar) polarization Free energy.

use precision_kinds, only : i2b, dp

use constants, only : i_complex, twopi, fourpi


use system, only : chg_solv, x_solv, y_solv, z_solv, nfft1, nfft2, nfft3, Lx, Ly, Lz,nb_solvent_sites,  nb_omega, nb_psi, id_solv&
, sigma_k,molec_polarx_k, molec_polary_k, molec_polarz_k,nb_species, deltaV, deltax, deltay,deltaz

use external_potential, only : x_charge, y_charge, z_charge, q_charge, nb_of_interpolation

use cg , only : cg_vect

use gauss_legendre, only : weight, weight_psi, Omx , Omy , Omz

use fft , only : kx, ky, kz, k2,in_forward , in_backward , out_forward , out_backward , plan_forward , plan_backward

implicit none

integer (kind=i2b) :: i, j, k, o , p , n,species, ii , jj , kk !dummy

real(kind=dp), dimension(nb_omega,nb_psi), intent(in) :: Rotxx,Rotxy,Rotxz,Rotyx,Rotyy,Rotyz,Rotzx,Rotzy,Rotzz

integer (kind=i2b) :: nf1

real (kind=dp) :: xmod, ymod, zmod

real (kind=dp) :: deltaVk

integer ( kind = i2b ) :: im , jm , km , ip , jp , kp ! indices of corner in indicial coordinates

real ( kind = dp ), dimension(2) :: wim , wjm , wkm ! weight associated to each index

real ( kind = dp ) :: xq , yq, zq ! coordinates of the atoms in indicial coordinates

real (kind=dp) , dimension (2) ::  xm, ym, zm  !coordinates of indexes where atomes are extrapolated , 1 is for the lower index, 2 for the upper
!            ====================================================
!            !    	Initialization				!
!            !							!
!            ====================================================



deltaVk=twopi**3/(Lx*Ly*Lz)


nf1=nfft1/2



allocate(sigma_k(nf1+1, nfft2, nfft3, nb_omega, nb_psi,nb_species))

sigma_k=(0.0_dp,0.0_dp)

allocate (molec_polarx_k (nf1+1 , nfft2 , nfft3 , nb_omega , nb_psi,nb_species))

molec_polarx_k= ( 0.0_dp , 0.0_dp )

allocate (molec_polary_k (nf1+1 , nfft2 , nfft3 , nb_omega , nb_psi,nb_species))

molec_polary_k= ( 0.0_dp , 0.0_dp )

allocate (molec_polarz_k (nf1+1 , nfft2 , nfft3 , nb_omega , nb_psi,nb_species))

molec_polarz_k= ( 0.0_dp , 0.0_dp )

!            ====================================================
!            !    Compute sigma and molecular polarization	!
!            !							!
!            ====================================================


do species=1,nb_species

do i = 1 , nf1 + 1
  
  do j = 1 , nfft2

    do k = 1 , nfft3

        do o=1, nb_omega

           do p=1, nb_psi

             do n=1, nb_solvent_sites


             xmod= Rotxx(o,p)*x_solv(n) + Rotxy(o,p)*y_solv(n) + Rotxz(o,p)*z_solv(n)
             ymod= Rotyx(o,p)*x_solv(n) + Rotyy(o,p)*y_solv(n) + Rotyz(o,p)*z_solv(n)   
             zmod= Rotzx(o,p)*x_solv(n) + Rotzy(o,p)*y_solv(n) + Rotzz(o,p)*z_solv(n)  

               if (xmod==Lx) then
               xmod=0.0_dp
               end if

               if (ymod==Ly) then
               ymod=0.0_dp
               end if


               if (zmod==Lz) then
               zmod=0.0_dp
               end if
!            ====================================================
!            !   Find indexes and weights of corners surrounding!
!            !		the charge				!
!            ====================================================

               xq =  xmod  / deltax
               yq =  ymod  / deltay
               zq =  zmod  / deltaz



               ! get coordinates of grid node juste below (corner of the cube with smallest indices)
               ! +1 is because indexation does not begin to 0 but to 1. Thus, when coordinate xq=0, it corresponds to index 1.

               im = int ( xq ) + 1
               jm = int ( yq ) + 1
               km = int ( zq ) + 1


               ! grid node juste above (corner of the cube with highest indices)

               ip = im + 1
               jp = jm + 1
               kp = km + 1


               ! this corner should never have coordinates higher than nfft

               if ( ip == nfft1 + 1 ) ip = 1
               if ( jp == nfft2 + 1 ) jp = 1
               if ( kp == nfft3 + 1 ) kp = 1

               !position of corners in real space

               xm(1)=(im-1)*deltax
               ym(1)=(jm-1)*deltay
               zm(1)=(km-1)*deltaz
               xm(2)=(ip-1)*deltax
               ym(2)=(jp-1)*deltay
               zm(2)=(kp-1)*deltaz

               ! define weights associated with each corner

               wim(1) = (    1.0_dp - (   xq - real(int(xq,kind=i2b),kind=dp)   )    )
               wjm(1) = (    1.0_dp - (   yq - real(int(yq,kind=i2b),kind=dp)   )    )
               wkm(1) = (    1.0_dp - (   zq - real(int(zq,kind=i2b),kind=dp)   )    )
               wim(2) = (             (   xq - real(int(xq,kind=i2b),kind=dp)   )    )
               wjm(2) = (             (   yq - real(int(yq,kind=i2b),kind=dp)   )    )
               wkm(2) = (             (   zq - real(int(zq,kind=i2b),kind=dp)   )    )


                 do ii=1,2
                   do jj=1,2
                     do kk=1,2
  


                      sigma_k (i , j, k ,o,p,species) = sigma_k (i , j , k ,o,p,species) +chg_solv(id_solv(n))*&
                      wim(ii)*wjm(jj)*wkm(kk)*exp(-i_complex*(kx(i)*xm(ii)+ky(j)*ym(jj)+ kz(k)*zm(kk) ) )
              

                        if ((xm(ii)*kx(i)+ym(jj)*ky(j)+zm(kk)*kz(k))==0.0_dp ) then


                        molec_polarx_k (i,j,k,o,p,species)=molec_polarx_k(i,j,k,o,p,species) +&
                        chg_solv(id_solv(n))*xm(ii)*wim(ii)*wjm(jj)*wkm(kk)

                        molec_polary_k (i,j,k,o,p,species)=molec_polary_k(i,j,k,o,p,species) + &
                        chg_solv(id_solv(n))*ym(jj)*wim(ii)*wjm(jj)*wkm(kk)

                        molec_polarz_k (i,j,k,o,p,species)=molec_polarz_k(i,j,k,o,p,species) + &
                        chg_solv(id_solv(n))*zm(kk)*wim(ii)*wjm(jj)*wkm(kk)


                        else

                        molec_polarx_k (i,j,k,o,p,species)=molec_polarx_k(i,j,k,o,p,species) + chg_solv(id_solv(n))*(-i_complex)*&
                        wim(ii)*wjm(jj)*wkm(kk)*xm(ii)/(xm(ii)*kx(i)+ym(jj)*ky(j)+zm(kk)*kz(k))&
                        *(exp(i_complex*(xm(ii)*kx(i)+ym(jj)*ky(j)+zm(kk)*kz(k)))-1)
              
                        molec_polary_k (i,j,k,o,p,species)=molec_polary_k(i,j,k,o,p,species) + chg_solv(id_solv(n))*(-i_complex)*&
                        wim(ii)*wjm(jj)*wkm(kk)*ym(jj)/(xm(ii)*kx(i)+ym(jj)*ky(j)+zm(kk)*kz(k))&
                        *(exp(i_complex*(xm(ii)*kx(i)+ym(jj)*ky(j)+zm(kk)*kz(k)))-1)


                        molec_polarz_k (i,j,k,o,p,species)=molec_polarz_k(i,j,k,o,p,species) + chg_solv(id_solv(n))*(-i_complex)*&
                        wim(ii)*wjm(jj)*wkm(kk)*zm(kk)/(xm(ii)*kx(i)+ym(jj)*ky(j)+zm(kk)*kz(k))&
                        *(exp(i_complex*(xm(ii)*kx(i)+ym(jj)*ky(j)+zm(kk)*kz(k)))-1)

                        end if

                    end do !ii

                  end do  !jj

                end do !kk

              end do  !n

           end do  !p

        end do  !o

     end do  !k

   end do  !j

end do  !i

end do!species

end subroutine

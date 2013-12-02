subroutine energy_polarization_with_nccoupling
use precision_kinds , only : i2b , dp
use system , only : nfft1 , nfft2 , nfft3 , nb_omega , Lx , Ly , Lz , c_delta , c_d , kBT , rho_0 , delta_k , nb_k ,&
                   deltav, nb_psi, molec_polarx_k,molec_polary_k, molec_polarz_k,delta_k,nb_k,kBT,&
                   rho_0_multispec, nb_species,pola_tot_x_k , pola_tot_y_k , pola_tot_z_k, deltax, rho_c_k_myway, chi_l, chi_t,&
                   n_0, beta,deltax,deltay,deltaz, sigma_k
use quadrature , only : weight , Omx , Omy , Omz, sym_order , weight_psi
use cg , only : cg_vect , FF , dF
use constants , only : twopi, i_complex, fourpi, eps0,qunit,Navo, qfact
use fft , only : in_forward , in_backward , out_forward , out_backward , plan_forward , plan_backward, kx, ky, kz, k2, norm_k
use input , only : input_line,input_log, input_char
implicit none
real (dp) , dimension (nfft1, nfft2, nfft3,nb_omega,nb_psi,nb_species) ::dF_pol_long , dF_pol_trans, dF_pol_tot
complex (dp) , dimension (nfft1/2+1, nfft2, nfft3, nb_omega, nb_psi, nb_species) :: rho_k, dF_pol_long_k , dF_pol_trans_k,&
 dF_pol_tot_k
real (dp) , dimension (nfft1, nfft2, nfft3, nb_omega, nb_psi, nb_species) ::rho
real (dp) , dimension (nfft1, nfft2, nfft3, nb_species) ::Px,Py,Pz,pola_tot_x,pola_tot_y,pola_tot_z,P_long_x,P_long_y,P_long_z
integer (i2b) :: icg ,nb_k_in_c  !dummy counter for cg_vect
real(dp) ::  deltaVk, F_pol_long, F_pol_trans , F_pol ,F_pol_tot  !Longitudinal , transverse and total Polarization free energy
complex(dp), allocatable, dimension(:,:,:,:) :: P_trans_x_k,P_trans_y_k,P_trans_z_k,P_long_x_k,P_long_y_k,P_long_z_k,&
pxk,pyk,pzk  !transverse part of polarization in Fourier space .
 
integer (i2b) :: i , j , k, o, p , n ,m,l, species, k_index , m1, m2, m3!dummy
real(dp) :: mu_SPCE, facsym,Pxt,Pyt,Pzt
complex(dp) :: k_tens_k_Px,k_tens_k_Py,k_tens_k_Pz     
real(dp):: time1 , time0 ,time2 , time3 ,rhot! timestamps
complex(dp) ::  F_pol_long_k , F_pol_trans_k , F_pol_tot_k 
real(dp), allocatable , dimension ( : ) :: weight_omx , weight_omy , weight_omz ! dummy
real(dp), allocatable, dimension(:) :: Ccc, Cnc, norm_k_in_c , Cnn
integer(i2b) :: ios ! iostat of the read statement: 
integer(i2b) ::kindex_in_C
real(dp) :: delta_k_in_C, rhon, fact,psi,Vint,Fint
real(dp), dimension(nfft1,nfft2,nfft3) :: rho_n  !one-particle number density
complex(dp), dimension(nfft1/2+1, nfft2, nfft3) :: rho_n_k
real(dp) ,  dimension ( nfft1 , nfft2 , nfft3 ) :: Vpair
complex(dp),  dimension(nfft1/2+1,nfft2,nfft3) :: Vpair_k
if (nb_species/=1) then
print*, 'transv_and_longi_polarization_micro IS NOT WORKING FOR MULTISPECIES'
stop
end if
! look for tag polarization in input
if (.not. input_log('polarization')) then
return
end if
!Check if you want to compute Polarization from a macroscopic point of view
!do i = 1 , size ( input_line )
!  j = len ( 'evaluate_polarization' )
!  if ( input_line (i) (1:j) == 'evaluate_polarization' .and. input_line (i) (j+4:j+8) == 'Macro' ) return ! exit this routine to get back to energy calculation skeleton.
!  
!end do
if (trim(adjustl(input_char('evaluate_polarization')))== 'dipol') return
if (.not. input_log('include_nc_coupling')) then
 return
end if

! init timer
call cpu_time(time0)
!            ====================================================
!            !    	Initialization				!
!            !							!
!            ====================================================
deltaVk=twopi**3/(Lx*Ly*Lz)
allocate (Pxk (nfft1/2+1, nfft2, nfft3, nb_species))
allocate (Pyk (nfft1/2+1, nfft2, nfft3, nb_species))
allocate (Pzk (nfft1/2+1, nfft2, nfft3, nb_species))
allocate (P_trans_x_k (nfft1/2+1, nfft2, nfft3, nb_species))
allocate (P_trans_y_k (nfft1/2+1, nfft2, nfft3, nb_species))
allocate (P_trans_z_k (nfft1/2+1, nfft2, nfft3, nb_species))
allocate (P_long_x_k (nfft1/2+1, nfft2, nfft3, nb_species))
allocate (P_long_y_k (nfft1/2+1, nfft2, nfft3, nb_species))
allocate (P_long_z_k (nfft1/2+1, nfft2, nfft3, nb_species))
P_trans_x_k=(0.0_dp,0.0_dp)
P_trans_y_k=(0.0_dp,0.0_dp)
P_trans_z_k=(0.0_dp,0.0_dp)
P_long_x_k=(0.0_dp,0.0_dp)
P_long_y_k=(0.0_dp,0.0_dp)
P_long_z_k=(0.0_dp,0.0_dp)
F_pol_long=0.0_dp
F_pol_trans=0.0_dp
F_pol_tot=0.0_dp
F_pol_long_k=(0.0_dp,0.0_dp)
F_pol_trans_k=(0.0_dp,0.0_dp)
F_pol_tot_k=(0.0_dp,0.0_dp)
F_pol=0.0_dp
dF_pol_long=0.0_dp
dF_pol_trans=0.0_dp
dF_pol_tot=0.0_dp
dF_pol_long_k=0.0_dp
dF_pol_trans_k=0.0_dp
dF_pol_tot_k=0.0_dp
icg=0
rho=0.0_dp
rho_k=(0.0_dp,0.0_dp)
if (.not. allocated (pola_tot_x_k ) )  allocate (pola_tot_x_k ( nfft1/2+1,nfft2,nfft3,nb_species) )
if (.not. allocated (pola_tot_y_k ) ) allocate (pola_tot_y_k ( nfft1/2+1,nfft2,nfft3,nb_species) )
if (.not. allocated (pola_tot_z_k ) ) allocate (pola_tot_z_k ( nfft1/2+1,nfft2,nfft3,nb_species) )
pola_tot_x_k=(0.0_dp,0.0_dp)
pola_tot_y_k=(0.0_dp,0.0_dp)
pola_tot_z_k=(0.0_dp,0.0_dp)
mu_SPCE=0.4238_dp*0.5773525_dp*2.0_dp !dipolar moment of SPCE water molecule in e.Angstromm
if (.not. allocated (rho_c_k_myway) ) allocate (rho_c_k_myway(nfft1/2+1, nfft2, nfft3))
rho_c_k_myway=(0.0_dp, 0.0_dp)
icg=0

print*,"This is the routine with nc coupling $$$$$$$$$$$"

!======================================================================================================
!            ====================================================
!            !    	Read Cnn(k), Cnc(k) and Ccc(k) 		!
             !       						!
             !							!
!            ====================================================

open(10,file='input/direct_correlation_functions/water/C_functions/Cnn.dat')

  nb_k_in_C=0
  do while(.true.)
    read(10,*,iostat=ios)
    if (ios>0) then
      write(*,*)'Error in compute_ck_dipolar.f90'
      write(*,*)'something went wrong during the computation of the total number of lines in cs.in. stop'
      stop
    else if (ios<0) then
      ! end of file reached
      exit
    else
      nb_k_in_C=nb_k_in_C+1
    end if
  end do

close(10)

allocate(norm_k_in_c(nb_k_in_c))
allocate(Ccc(nb_k_in_c))
allocate(Cnc(nb_k_in_c))
allocate(Cnn(nb_k_in_c))

open(10,file='input/direct_correlation_functions/water/C_functions/Ccc.dat')
open(11,file='input/direct_correlation_functions/water/C_functions/Cnc.dat')
open(12,file='input/direct_correlation_functions/water/C_functions/Cnn.dat')

do i=1, nb_k_in_c
  read(10,*) norm_k_in_c(i), Ccc(i)
  read(11,*) norm_k_in_c(i), Cnc(i)
  read(12,*) norm_k_in_c(i), Cnn(i)
end do
delta_k_in_c=norm_k_in_c(2)-norm_k_in_c(1)
print*, delta_k_in_c

close(10)
close(11)
close(12)
!======================================================================================================
!            ====================================================
!            !    	Compute density in real and 		!
             !       		Fourier space			!
             !							!
!            ====================================================
!Compute rho_k
do species=1,nb_species
  do i=1,nfft1
    do j=1,nfft2
      do k=1,nfft3
        rhon=0.0_dp
        do o=1,nb_omega
          do p=1,nb_psi
            icg = icg + 1
            rhon=rhon+cg_vect(icg)**2*weight(o)*weight_psi(p)!/(fourpi**2/(sym_order*2.0_dp))

            rho(i,j,k,o,p,species) = cg_vect ( icg ) ** 2
          end do
        end do
           rho_n(i,j,k)=rhon-real(2.0_dp*twopi**2/sym_order, dp)
      end do
    end do
  end do
end do
call cpu_time(time2)
!======================================================================================================
!            ====================================================
!            !    		get Density 			!
!            !			in Fourier Space		!
!            ====================================================
do species=1, nb_species
   do o=1, nb_omega
      do p=1, nb_psi
      in_forward=rho(:,:,:,o,p,species)
      call dfftw_execute (plan_forward)
      rho_k(:,:,:,o,p,species)=out_forward*deltaV
  
      end do
   end do
end do
    
      in_forward=rho_n
      call dfftw_execute (plan_forward)
      rho_n_k=out_forward*deltaV
       

call cpu_time(time3)
!do o=1,nb_omega
!     molec_polarx_k(:,:,:,o,:,:) = mu_SPCE*Omx(o)
!    molec_polary_k(:,:,:,o,:,:) = mu_SPCE*Omy(o)
!    molec_polarz_k(:,:,:,o,:,:) = mu_SPCE*Omz(o)
!end do

!======================================================================================================



!            ====================================================
!            !    		Compute 			!
!            !		part due to density/density coupling	!
!	     ! with cnn including a nc coupling part		!
!            ====================================================


Vpair=0.0_dp
Vpair_k=(0.0_dp,0.0_dp)
do n = 1 , nfft3
  do m = 1 , nfft2
    do l = 1 , nfft1 / 2 + 1
       kindex_in_c=int ( norm_k(l,m,n) / delta_k_in_C ) + 1
      ! Here it happens that k_index gets higher than the highest c_k index.
      ! In this case one imposes k_index = k_index_max
       if ( kindex_in_c > nb_k_in_c ) kindex_in_c = nb_k_in_c 

      Vpair_k ( l , m , n ) = rho_n_k ( l , m , n ) * Cnn ( kindex_in_c )   !cnn is cs if there is no coupling
    end do
  end do
end do

in_backward = Vpair_k

call dfftw_execute (plan_backward)

Vpair = out_backward *deltaVk/(twopi**3)
! compute excess energy and its gradient
Fint = 0.0_dp ! excess energy
icg = 0 ! index of cg_vect
do species = 1 , nb_species
  fact = DeltaV * rho_0_multispec ( species ) !> facteur d'integration
  do i = 1 , nfft1
    do j = 1 , nfft2
      do k = 1 , nfft3
        Vint   = -kBT * rho_0_multispec ( species ) * Vpair(i,j,k)
	!Fint=Fint-0.5* rho_n(i,j,k)**2*deltaV*rho_0_multispec ( species )**2/n_0*kBT
        do o = 1 , nb_omega
          do p=1, nb_psi
          icg = icg + 1
          psi = CG_vect ( icg )
          Fint   = Fint   + weight(o) * fact * 0.5_dp * ( psi ** 2 - 1.0_dp) * Vint*weight_psi(p)
			

          dF (icg) = dF ( icg ) + 2.0_dp * psi * weight(o) * fact * Vint*weight_psi(p)!&
			!-0.5_dp*2.0_dp*rho_0_multispec ( species )**2/n_0*2.0*psi*weight(o)*weight_psi(p)*rho_n(i,j,k)*kBT*deltaV
         end do
   
        end do
      end do
    end do
  end do
end do ! species

print*, 'Fnn=', Fint


!======================================================================================================



!            ====================================================
!            !    		Compute 			!
!            !		Total Polarization			!
!            ====================================================

!stop

do species=1, nb_species
  do i=1, nfft1/2+1
     do j=1, nfft2
        do k=1, nfft3
           do o=1, nb_omega
              do p=1, nb_psi
              pola_tot_x_k(i,j,k,species)=pola_tot_x_k(i,j,k,species)+weight(o)*weight_psi(p)*&
              molec_polarx_k(i,j,k,o,p,species)*rho_k(i,j,k,o,p,species)
              pola_tot_y_k(i,j,k,species)=pola_tot_y_k(i,j,k,species)+weight(o)*weight_psi(p)*&
              molec_polary_k(i,j,k,o,p,species)*rho_k(i,j,k,o,p,species)
              pola_tot_z_k(i,j,k,species)=pola_tot_z_k(i,j,k,species)+weight(o)*weight_psi(p)*&
              molec_polarz_k(i,j,k,o,p,species)*rho_k(i,j,k,o,p,species)
              end do
            end do
         end do
      
      end do
   end do
end do
do species=1,nb_species
in_backward= pola_tot_x_k(:,:,:,species)
call dfftw_execute(plan_backward)
pola_tot_x(:,:,:,species)=out_backward*deltaVk/(twopi)**3
in_backward= pola_tot_y_k(:,:,:,species)
call dfftw_execute(plan_backward)
pola_tot_y(:,:,:,species)=out_backward*deltaVk/(twopi)**3
in_backward= pola_tot_z_k(:,:,:,species)
call dfftw_execute(plan_backward)
pola_tot_z(:,:,:,species)=out_backward*deltaVk/(twopi)**3
end do

!!            ====================================================
!!            !    	Compute 				!
!!            !	Transverse and longitudinal Polarization	!
!!            ====================================================
do n=1, nb_species
  do i=1, nfft1/2+1
    do j=1, nfft2
      do k=1, nfft3
      P_long_x_k(i,j,k,n)=(pola_tot_x_k(i,j,k,n)*kx(i)+pola_tot_y_k(i,j,k,n)*ky(j)+pola_tot_z_k(i,j,k,n)*kz(k))*kx(i)/k2(i,j,k)
      P_long_y_k(i,j,k,n)=(pola_tot_x_k(i,j,k,n)*kx(i)+pola_tot_y_k(i,j,k,n)*ky(j)+pola_tot_z_k(i,j,k,n)*kz(k))*ky(j)/k2(i,j,k)
      P_long_z_k(i,j,k,n)=(pola_tot_x_k(i,j,k,n)*kx(i)+pola_tot_y_k(i,j,k,n)*ky(j)+pola_tot_z_k(i,j,k,n)*kz(k))*kz(k)/k2(i,j,k)
      
      if (k2(i,j,k)==0.0_dp) then
      P_long_x_k(i,j,k,n)=0.0_dp
      P_long_y_k(i,j,k,n)=0.0_dp
      P_long_z_k(i,j,k,n)=0.0_dp
      end if
      P_trans_x_k(i,j,k,n)=pola_tot_x_k(i,j,k,n)- P_long_x_k(i,j,k,n)
      P_trans_y_k(i,j,k,n)=pola_tot_y_k(i,j,k,n)- P_long_y_k(i,j,k,n)
      P_trans_z_k(i,j,k,n)=pola_tot_z_k(i,j,k,n)- P_long_z_k(i,j,k,n)
      end do
    end do
  end do
end do
do species=1,nb_species
in_backward= P_long_x_k(:,:,:,species)
call dfftw_execute(plan_backward)
P_long_x(:,:,:,species)=out_backward*deltaVk/(twopi)**3
in_backward= P_long_y_k(:,:,:,species)
call dfftw_execute(plan_backward)
P_long_y(:,:,:,species)=out_backward*deltaVk/(twopi)**3
in_backward= P_long_z_k(:,:,:,species)
call dfftw_execute(plan_backward)
P_long_z(:,:,:,species)=out_backward*deltaVk/(twopi)**3
end do
open (11, file='output/pola_tot_x')
open (12, file='output/pola_tot_y')
open (13, file='output/pola_tot_z')
do i=1,nfft1
write(11,*) , i*deltax, pola_tot_x(i,nfft2/2+1,nfft3/2+1,1), Px(i,nfft2/2+1,nfft3/2+1,1),P_long_x(i,nfft2/2+1,nfft3/2+1,1) 
write(12,*) , i*deltay, pola_tot_y(nfft1/2+1,i,nfft3/2+1,1), Py(nfft1/2+1,i,nfft3/2+1,1),P_long_y(nfft1/2+1,i,nfft3/2+1,1) 
write(13,*) , i*deltaz, pola_tot_z(nfft1/2+1,nfft2/2+1,i,1), Pz(nfft1/2+1,nfft1/2+1,i,1),P_long_z(nfft1/2+1,nfft1/2+1,i,1)
end do
close(11)
close(12)
close(13)

!            ====================================================
            !    	Compute rho_c_k	in k space		!
            !							!
!            ====================================================

do species = 1, nb_species
  do i = 1 , nfft1/2 + 1
    do j = 1 , nfft2
      do k = 1 , nfft3
        do o = 1 , nb_omega
            do p=1 , nb_psi
              rho_c_k_myway(i,j,k)=rho_c_k_myway(i,j,k)+weight(o)*weight_psi(p)*sigma_k(i,j,k,o,p,species)*&
              rho_k (i , j , k , o , p , species)*rho_0!*deltavK         
            end do  !p
          end do   !o
       end do  !k
     end do  !j
  end do  !i
end do

!            ====================================================
!            !    	Compute Free energy due to		!
!            !	Transverse and longitudinal Polarization	!
!            ====================================================
do species=1, nb_species 
  do i=1, nfft1/2+1
  
    if (i>1 .and. i<nfft1/2+1) then
    facsym=2.0_dp
    else
    facsym=1.0_dp
    end if
    do j=1, nfft2
      do k=1, nfft3
       k_index = int ( norm_k(i,j,k) / delta_k ) + 1
       kindex_in_c=int ( norm_k(i,j,k) / delta_k_in_C ) + 1
       ! Here it happens that k_index gets higher than the highest c_k index.
       ! In this case one imposes k_index = k_index_max
       if ( kindex_in_c > nb_k_in_c ) kindex_in_c = nb_k_in_c 
       if ( k_index > nb_k ) k_index = nb_k
       if(k2(i,j,k)/=0.0_dp) then
       F_pol_long_k=F_pol_long_k+deltaVk*kBT*(P_long_x_k(i,j,k,species)*conjg(p_long_x_k(i,j,k,species))+&
       p_long_y_k(i,j,k,species)*conjg(p_long_y_k(i,j,k,species))+p_long_z_k(i,j,k,species)*conjg(p_long_z_k(i,j,k,species)))&
       *Ccc(kindex_in_c)*0.5_dp*rho_0**2*facsym/(twopi**3)+&
        (deltaVk*(rho_c_k_myway(i,j,k)*conjg(rho_n_k(i,j,k)))&
       *rho_0*kBT*facsym*Cnc(kindex_in_c))/(twopi)**3/norm_k(i,j,k)
     end if

       F_pol_trans_k=F_pol_trans_k+deltaVk*(P_trans_x_k(i,j,k,species)*conjg(p_trans_x_k(i,j,k,species))+p_trans_y_k(i,j,k,species)&
       *conjg(p_trans_y_k(i,j,k,species))+p_trans_z_k(i,j,k,species)*conjg(p_trans_z_k(i,j,k,species)))/chi_t(k_index)&
       *0.5_dp*qfact*rho_0**2*facsym/(twopi**3)

       F_pol_tot_k=F_pol_tot_k-(kBT*3/(2*mu_SPCE**2*n_0)*deltaVk*facsym/(twopi**3)*rho_0**2*&
       (pola_tot_x_k(i,j,k,species)*conjg(pola_tot_x_k(i,j,k,species))+pola_tot_y_k(i,j,k,species)*&
       conjg(pola_tot_y_k(i,j,k,species))+pola_tot_z_k(i,j,k,species)*conjg(pola_tot_z_k(i,j,k,species))))
        do o=1,nb_omega
           do p=1,nb_psi

!!========================================================================================================================
!!Evaluate gradient
!!========================================================================================================================
k_tens_k_Px=conjg(kx(i)*kx(i)*molec_polarx_k(i,j,k,o,p,species)+kx(i)*ky(j)*molec_polary_k(i,j,k,o,p,species)+kx(i)*kz(k)*&
molec_polarz_k(i,j,k,o,p,species))
k_tens_k_Py=conjg(kx(i)*ky(j)*molec_polarx_k(i,j,k,o,p,species)+ky(j)*ky(j)*molec_polary_k(i,j,k,o,p,species)+ky(j)*kz(k)*&
molec_polarz_k(i,j,k,o,p,species))
k_tens_k_Pz=conjg(kz(k)*kx(i)*molec_polarx_k(i,j,k,o,p,species)+kz(k)*ky(j)*molec_polary_k(i,j,k,o,p,species)+kz(k)*kz(k)*&
molec_polarz_k(i,j,k,o,p,species))
   
  if (k2(i,j,k)==0.0_dp) then
     dF_pol_trans_k(i,j,k,o,p,species)=rho_0*0.5_dp*qfact/chi_t(k_index)*&
     (P_trans_x_k(i,j,k,species)*conjg(molec_polarx_k(i,j,k,o,p,species))&
      +P_trans_y_k(i,j,k,species)*conjg(molec_polary_k(i,j,k,o,p,species))&
     +P_trans_z_k(i,j,k,species)*conjg(molec_polarz_k(i,j,k,o,p,species)))*weight(o)*weight_psi(p)*2.0_dp
     dF_pol_long_k(i,j,k,o,p,species)=0.0_dp

!             ===============================================================================
else
      dF_pol_trans_k(i,j,k,o,p,species) = rho_0*0.5_dp*qfact/chi_t(k_index)*(P_trans_x_k(i,j,k,species)*&
      (conjg(molec_polarx_k(i,j,k,o,p,species))-k_tens_k_Px/k2(i,j,k))&
+P_trans_y_k(i,j,k,species)*(conjg(molec_polary_k(i,j,k,o,p,species))-k_tens_k_Py/k2(i,j,k))&
+P_trans_z_k(i,j,k,species)*(conjg(molec_polarz_k(i,j,k,o,p,species))-k_tens_k_Pz/k2(i,j,k)) )*weight(o)*weight_psi(p)*2.0_dp
     end if
!             ===============================================================================
dF_pol_tot_k(i,j,k,o,p,species)=-kBT*3*rho_0/(2*mu_SPCE**2*n_0)*(pola_tot_x_k(i,j,k,species)*&
conjg(molec_polarx_k(i,j,k,o,p,species))+pola_tot_y_k(i,j,k,species)*conjg(molec_polary_k(i,j,k,o,p,species))+pola_tot_z_k&  
(i,j,k,species)*conjg(molec_polarz_k(i,j,k,o,p,species)))*weight(o)*weight_psi(p)*2.0_dp

    if(k2(i,j,k)/=0.0_dp) then
       dF_pol_long_k(i,j,k,o,p,species) = rho_0*0.5_dp*Ccc(kindex_in_c)*kBT*weight(o)*weight_psi(p)*&
      (P_long_x_k(i,j,k,species)*k_tens_k_Px+P_long_y_k(i,j,k,species)*k_tens_k_Py+P_long_z_k(i,j,k,species)&
      *k_tens_k_Pz)/k2(i,j,k)*2.0_dp!+rho_0*2.0*Cnc(kindex_in_c)*(rho_n_k(i,j,k)*conjg(sigma_k(i,j,k,o,p,species))+&
   !   rho_c_k_myway(i,j,k)/rho_0)*weight(o)*weight_psi(p)/norm_k(i,j,k)
      end if

          end do  !psi
        end do   !omega
     
      end do    !k
    end do    !j
  end do    !i
end do     !species

!========================================================================================================================
!						Get gradient in real space
!========================================================================================================================
do species=1 , nb_species
  do o=1,nb_omega
    do p=1, nb_psi
    in_backward= (dF_pol_trans_k (:,:,:,o,p,species)+dF_pol_long_k (:,:,:,o,p,species)+dF_pol_tot_k (:,:,:,o,p,species))
    call dfftw_execute (plan_backward)
    dF_pol_tot (:,:,:,o,p,species)=out_backward*deltaVk/(twopi)**3
    end do
  end do
end do
icg=0
!!========================================================================================================================
!!						Allocate it for minimizing
!!========================================================================================================================
do species=1, nb_species
  do i=1,nfft1
    do j=1, nfft2
      do k=1,nfft3
    
        do o=1,nb_omega
          do p=1, nb_psi
          icg=icg+1
          dF(icg)=dF(icg)+  dF_pol_tot (i,j,k,o,p,species)*cg_vect(icg)*rho_0*2.0_dp * deltav!* weight(o) * weight_psi(p) 
          end do
 
        end do
      end do
    end do
  end do
end do
!!========================================================================================================================
!!					Check if Polarization Free energy is Real
!!========================================================================================================================
if (aimag(F_pol_tot_k+F_pol_long_k+F_pol_trans_k)<tiny(0.0_dp)) then 
F_pol_tot=real( F_pol_tot_k , dp)
F_pol_long=real( F_pol_long_k,  dp)
F_pol_trans=real( F_pol_trans_k, dp)
else
print*, 'Error in energy_polarization_withcoupling Free energy is not Real'
print*,aimag(F_pol_tot_k+F_pol_long_k+F_pol_trans_k)
print*,F_pol_tot_k+F_pol_long_k+F_pol_trans_k
stop
end if
!!!========================================================================================================================
F_pol=F_pol_tot+F_pol_long+F_pol_trans
FF=FF+Fint!+F_pol
! stop timer
call cpu_time ( time1 )
print*, 'F_polarization =' , F_pol  , 'computed in (sec)' , time1 - time0
end subroutine

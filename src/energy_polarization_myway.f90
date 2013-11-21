SUBROUTINE energy_polarization_myway (F_pol)

    USE precision_kinds,only : i2b , dp
    use system,only : nfft1 , nfft2 , nfft3 , Lx , Ly , Lz , c_delta , c_d , kBT , rho_0 , delta_k , nb_k ,&
                    deltav, molec_polarx_k,molec_polary_k, molec_polarz_k,delta_k,nb_k,kBT,&
                    rho_0_multispec, nb_species,pola_tot_x_k , pola_tot_y_k , pola_tot_z_k, deltax, rho_c_k_myway, chi_l, chi_t,&
                    n_0, beta,deltax,deltay,deltaz
    use quadrature,only : Omx , Omy , Omz, sym_order , angGrid, molRotGrid
    USE minimizer, ONLY: cg_vect , FF , dF
    use constants,only : twopi, i_complex, fourpi, eps0,qunit,Navo, qfact
    use fft,only : fftw3, kx, ky, kz, k2, norm_k
    use input,only : input_line,input_log, input_char, verbose

IMPLICIT NONE
    
    REAL(dp), INTENT(OUT) :: F_pol
real(dp), dimension (nfft1, nfft2, nfft3,angGrid%n_angles,molRotGrid%n_angles,nb_species) ::dF_pol_long,dF_pol_trans,dF_pol_tot
complex (dp) , dimension (nfft1/2+1, nfft2, nfft3, angGrid%n_angles, molRotGrid%n_angles, nb_species) :: rho_k, dF_pol_long_k ,&
    dF_pol_trans_k, dF_pol_tot_k
real (dp) , dimension (nfft1, nfft2, nfft3, angGrid%n_angles, molRotGrid%n_angles, nb_species) ::rho
real (dp) , dimension (nfft1, nfft2, nfft3, nb_species) ::Px,Py,Pz,pola_tot_x,pola_tot_y,pola_tot_z,P_long_x,P_long_y,P_long_z
integer (i2b) :: icg   !dummy counter for cg_vect
real(dp) ::  deltaVk, F_pol_long, F_pol_trans , F_pol_tot  !Longitudinal , transverse and total Polarization free energy
complex(dp), allocatable, dimension(:,:,:,:) :: P_trans_x_k,P_trans_y_k,P_trans_z_k,P_long_x_k,P_long_y_k,P_long_z_k,&
pxk,pyk,pzk  !transverse part of polarization in Fourier space .
 
integer (i2b) :: i , j , k, o, p , n , species, k_index , m1, m2, m3!dummy
real(dp) :: mu_SPCE, facsym,Pxt,Pyt,Pzt
complex(dp) :: k_tens_k_Px,k_tens_k_Py,k_tens_k_Pz     
real(dp):: time1 , time0 ,time2 , time3 ,rhot! timestamps
complex(dp) ::  F_pol_long_k , F_pol_trans_k , F_pol_tot_k 
real(dp), allocatable , dimension ( : ) :: weight_omx , weight_omy , weight_omz ! dummy

if (nb_species/=1) then
    print*, 'transv_and_longi_polarization_micro IS NOT WORKING FOR MULTISPECIES'
    stop
END IF

! look for tag polarization in input
if (.not. input_log('polarization')) return

IF (.NOT. ALLOCATED (rho_c_k_myway)) THEN
    ALLOCATE (rho_c_k_myway(nfft1/2+1, nfft2, nfft3) )
    rho_c_k_myway = (0._dp,0._dp)
END IF

!Check if you want to compute Polarization from a macroscopic point of view
!do i = 1 , size ( input_line )
!  j = len ( 'evaluate_polarization' )
!  if ( input_line (i) (1:j) == 'evaluate_polarization' .and. input_line (i) (j+4:j+8) == 'Macro' ) return ! exit this routine to get back to energy calculation skeleton.
!  
!END DO
if (trim(adjustl(input_char('evaluate_polarization')))== 'dipol') return
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
!======================================================================================================
!            ====================================================
!            !    	Compute density in real and 		!
             !       		Fourier space			!
             !							!
!            ====================================================
!Compute rho_k
do species =1 , nb_species
   do i=1,nfft1
     do j=1, nfft2
       do k=1, nfft3
         do o = 1 , angGrid%n_angles
            do p=1, molRotGrid%n_angles
            icg = icg + 1
            rho(i,j,k,o,p,species) = cg_vect ( icg ) ** 2
          END DO
        END DO
      END DO
    END DO
  END DO
END DO
call cpu_time(time2)
!======================================================================================================
!            ====================================================
!            !    		get Density 			!
!            !			in Fourier Space		!
!            ====================================================
DO species = 1, nb_species
    DO o = 1, angGrid%n_angles
        DO p = 1, molRotGrid%n_angles
            fftw3%in_forward=rho(:,:,:,o,p,species)
            CALL dfftw_execute (fftw3%plan_forward)
            rho_k(:,:,:,o,p,species)=fftw3%out_forward*deltaV
        END DO
    END DO
END DO
call cpu_time(time3)
!do o=1,angGrid%n_angles
!    molec_polarx_k(:,:,:,o,:,:) = mu_SPCE*Omx(o)
!    molec_polary_k(:,:,:,o,:,:) = mu_SPCE*Omy(o)
!    molec_polarz_k(:,:,:,o,:,:) = mu_SPCE*Omz(o)
!END DO
!======================================================================================================
!            ====================================================
!            !    		Compute 			!
!            !		Total Polarization			!
!            ====================================================
do species=1, nb_species
  do i=1, nfft1/2+1
     do j=1, nfft2
        do k=1, nfft3
           do o=1, angGrid%n_angles
              do p=1, molRotGrid%n_angles
                pola_tot_x_k(i,j,k,species)=pola_tot_x_k(i,j,k,species)+angGrid%weight(o)*molRotGrid%weight(p)*&
                molec_polarx_k(i,j,k,o,p,species)*rho_k(i,j,k,o,p,species)
                pola_tot_y_k(i,j,k,species)=pola_tot_y_k(i,j,k,species)+angGrid%weight(o)*molRotGrid%weight(p)*&
                molec_polary_k(i,j,k,o,p,species)*rho_k(i,j,k,o,p,species)
                pola_tot_z_k(i,j,k,species)=pola_tot_z_k(i,j,k,species)+angGrid%weight(o)*molRotGrid%weight(p)*&
                molec_polarz_k(i,j,k,o,p,species)*rho_k(i,j,k,o,p,species)
              END DO
            END DO
         END DO
      END DO
   END DO
END DO

DO species=1,nb_species
    fftw3%in_backward= pola_tot_x_k(:,:,:,species)
    call dfftw_execute(fftw3%plan_backward)
    pola_tot_x(:,:,:,species)=fftw3%out_backward*deltaVk/(twopi)**3
    fftw3%in_backward= pola_tot_y_k(:,:,:,species)
    call dfftw_execute(fftw3%plan_backward)
    pola_tot_y(:,:,:,species)=fftw3%out_backward*deltaVk/(twopi)**3
    fftw3%in_backward= pola_tot_z_k(:,:,:,species)
    call dfftw_execute(fftw3%plan_backward)
    pola_tot_z(:,:,:,species)=fftw3%out_backward*deltaVk/(twopi)**3
END DO

!allocate ( weight_omx ( angGrid%n_angles ) )
!allocate ( weight_omy ( angGrid%n_angles ) )
!allocate ( weight_omz ( angGrid%n_angles ) )
!weight_omx = angGrid%weight * Omx
!weight_omy = angGrid%weight * Omy
!weight_omz = angGrid%weight * Omz
!icg=0
!do species=1,nb_species
!do i = 1 , nfft1
!    
!    m1=i-1
!    if (i> nfft1/2) m1=i-1-nfft1
!  do j = 1 , nfft2
!    m2=j-1
!    if (i>nfft2/2) m2=j-1-nfft2
!    do k = 1 , nfft3    
!    m3=k-1
!    if (k>nfft3/2) m3=k-1-nfft3
!      ! init dummy variables tpx , tpy and tpz in order not to loop directly over big arrays
!      pxt = 0.0_dp
!      pyt = 0.0_dp
!      pzt = 0.0_dp
!      do o = 1 , angGrid%n_angles
!        do p=1, molRotGrid%n_angles
!          icg = icg + 1
!          rhot = cg_vect (icg) ** 2
!          pxt = pxt + weight_Omx ( o ) * molRotGrid%weight(p) * rhot
!          pyt = pyt + weight_Omy ( o ) * molRotGrid%weight(p) * rhot
!          pzt = pzt + weight_Omz ( o ) * molRotGrid%weight(p) * rhot
!        END DO
!      END DO
!      Px ( i , j , k,species ) = pxt
!      Py ( i , j , k ,species) = pyt
!      Pz ( i , j , k ,species) = pzt
!    END DO
!  END DO
!END DO
!END DO
!======================================================================================================
!            ====================================================
!            !    	Compute 				!
!            !	FFT Polarization	!
!            ====================================================
!do species=1,nb_species
!print*, shape(Pz(:,:,:,:)),shape(in_forward), shape(pola_tot_z_k(:,:,:,:)) , shape(out_forward)
!in_forward=pola_tot_x(:,:,:,species)
!call dfftw_execute ( plan_forward )
!pola_tot_x_k(:,:,:,species)=out_forward*deltaV*mu_SPCE
!in_forward=pola_tot_x_k(:,:,:,species)
!call dfftw_execute ( plan_forward )
!pola_tot_y_k(:,:,:,species)=out_forward*deltaV*mu_SPCE
!in_forward=Pz(:,:,:,species)
!call dfftw_execute ( plan_forward )
!pola_tot_z_k(:,:,:,species)=out_forward*deltaV*mu_SPCE
!END DO
!!======================================================================================================
!            ====================================================
!            !    	Compute 				!
!            !	Transverse and longitudinal Polarization	!
!            ====================================================
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
      END IF
      P_trans_x_k(i,j,k,n)=pola_tot_x_k(i,j,k,n)- P_long_x_k(i,j,k,n)
      P_trans_y_k(i,j,k,n)=pola_tot_y_k(i,j,k,n)- P_long_y_k(i,j,k,n)
      P_trans_z_k(i,j,k,n)=pola_tot_z_k(i,j,k,n)- P_long_z_k(i,j,k,n)
      END DO
    END DO
  END DO
END DO
do species=1,nb_species
fftw3%in_backward= P_long_x_k(:,:,:,species)
call dfftw_execute(fftw3%plan_backward)
P_long_x(:,:,:,species)=fftw3%out_backward*deltaVk/(twopi)**3
fftw3%in_backward= P_long_y_k(:,:,:,species)
call dfftw_execute(fftw3%plan_backward)
P_long_y(:,:,:,species)=fftw3%out_backward*deltaVk/(twopi)**3
fftw3%in_backward= P_long_z_k(:,:,:,species)
call dfftw_execute(fftw3%plan_backward)
P_long_z(:,:,:,species)=fftw3%out_backward*deltaVk/(twopi)**3
END DO

open (11, file='output/pola_tot_x')
open (12, file='output/pola_tot_y')
open (13, file='output/pola_tot_z')
    do i=1,nfft1
        write(11,*) , i*deltax, pola_tot_x(i,nfft2/2+1,nfft3/2+1,1), Px(i,nfft2/2+1,nfft3/2+1,1),P_long_x(i,nfft2/2+1,nfft3/2+1,1) 
        write(12,*) , i*deltay, pola_tot_y(nfft1/2+1,i,nfft3/2+1,1), Py(nfft1/2+1,i,nfft3/2+1,1),P_long_y(nfft1/2+1,i,nfft3/2+1,1) 
        write(13,*) , i*deltaz, pola_tot_z(nfft1/2+1,nfft2/2+1,i,1), Pz(nfft1/2+1,nfft1/2+1,i,1),P_long_z(nfft1/2+1,nfft1/2+1,i,1)
    END DO
close(11)
close(12)
close(13)
!            ====================================================
!            !    	Compute Free energy due to		!
!            !	Transverse and longitudinal Polarization	!
!            ====================================================
do species=1, nb_species 
  do i=1, nfft1/2+1
  
    if (i>1 .and. i<nfft1/2+1) then
    facsym=2.0_dp
    ELSE
    facsym=1.0_dp
    END IF
    do j=1, nfft2
      do k=1, nfft3
       k_index = int ( norm_k(i,j,k) / delta_k ) + 1
       ! Here it happens that k_index gets higher than the highest c_k index.
       ! In this case one imposes k_index = k_index_max
       if ( k_index > nb_k ) k_index = nb_k 
       F_pol_long_k=F_pol_long_k+deltaVk*fourpi*(P_long_x_k(i,j,k,species)*conjg(p_long_x_k(i,j,k,species))+&
       p_long_y_k(i,j,k,species)*conjg(p_long_y_k(i,j,k,species))+p_long_z_k(i,j,k,species)*conjg(p_long_z_k(i,j,k,species)))&
       /chi_l(k_index)*0.5_dp*qfact*rho_0**2*facsym/(twopi**3)
       F_pol_trans_k=F_pol_trans_k+deltaVk*(P_trans_x_k(i,j,k,species)*conjg(p_trans_x_k(i,j,k,species))+p_trans_y_k(i,j,k,species)&
       *conjg(p_trans_y_k(i,j,k,species))+p_trans_z_k(i,j,k,species)*conjg(p_trans_z_k(i,j,k,species)))/chi_t(k_index)&
       *0.5_dp*qfact*rho_0**2*facsym/(twopi**3)
       F_pol_tot_k=F_pol_tot_k-kBT*3/(2*mu_SPCE**2*n_0)*deltaVk*facsym/(twopi**3)*rho_0**2*&
       (pola_tot_x_k(i,j,k,species)*conjg(pola_tot_x_k(i,j,k,species))+pola_tot_y_k(i,j,k,species)*&
       conjg(pola_tot_y_k(i,j,k,species))+pola_tot_z_k(i,j,k,species)*conjg(pola_tot_z_k(i,j,k,species)))
        do o=1,angGrid%n_angles
           do p=1,molRotGrid%n_angles 
!========================================================================================================================
!Evaluate gradient
!========================================================================================================================
k_tens_k_Px=conjg(kx(i)*kx(i)*molec_polarx_k(i,j,k,o,p,species)+kx(i)*ky(j)*molec_polary_k(i,j,k,o,p,species)+kx(i)*kz(k)*&
molec_polarz_k(i,j,k,o,p,species))
k_tens_k_Py=conjg(kx(i)*ky(j)*molec_polarx_k(i,j,k,o,p,species)+ky(j)*ky(j)*molec_polary_k(i,j,k,o,p,species)+ky(j)*kz(k)*&
molec_polarz_k(i,j,k,o,p,species))
k_tens_k_Pz=conjg(kz(k)*kx(i)*molec_polarx_k(i,j,k,o,p,species)+kz(k)*ky(j)*molec_polary_k(i,j,k,o,p,species)+kz(k)*kz(k)*&
molec_polarz_k(i,j,k,o,p,species))
   
  if (k2(i,j,k)==0.0_dp) then
     dF_pol_long_k(i,j,k,o,p,species)=0.0_dp
     dF_pol_trans_k(i,j,k,o,p,species)=rho_0*0.5_dp*qfact/chi_t(k_index)*&
     (P_trans_x_k(i,j,k,species)*conjg(molec_polarx_k(i,j,k,o,p,species))&
      +P_trans_y_k(i,j,k,species)*conjg(molec_polary_k(i,j,k,o,p,species))&
     +P_trans_z_k(i,j,k,species)*conjg(molec_polarz_k(i,j,k,o,p,species)))*angGrid%weight(o)*molRotGrid%weight(p)*2.0_dp
     ELSE
      dF_pol_long_k(i,j,k,o,p,species) = rho_0*0.5_dp*qfact/chi_l(k_index)*fourpi*angGrid%weight(o)*molRotGrid%weight(p)*&
      (P_long_x_k(i,j,k,species)*k_tens_k_Px+P_long_y_k(i,j,k,species)*k_tens_k_Py+P_long_z_k(i,j,k,species)&
      *k_tens_k_Pz)/k2(i,j,k)*2.0_dp
!             ===============================================================================
      dF_pol_trans_k(i,j,k,o,p,species) = rho_0*0.5_dp*qfact/chi_t(k_index)*(P_trans_x_k(i,j,k,species)*&
      (conjg(molec_polarx_k(i,j,k,o,p,species))-k_tens_k_Px/k2(i,j,k))&
+P_trans_y_k(i,j,k,species)*(conjg(molec_polary_k(i,j,k,o,p,species))-k_tens_k_Py/k2(i,j,k))&
+P_trans_z_k(i,j,k,species)*(conjg(molec_polarz_k(i,j,k,o,p,species))-k_tens_k_Pz/k2(i,j,k)) )&
*angGrid%weight(o)*molRotGrid%weight(p)*2.0_dp
     END IF
!             ===============================================================================
dF_pol_tot_k(i,j,k,o,p,species)=-kBT*3*rho_0/(2*mu_SPCE**2*n_0)*(pola_tot_x_k(i,j,k,species)*&
conjg(molec_polarx_k(i,j,k,o,p,species))+pola_tot_y_k(i,j,k,species)*conjg(molec_polary_k(i,j,k,o,p,species))+pola_tot_z_k&  
(i,j,k,species)*conjg(molec_polarz_k(i,j,k,o,p,species)))*angGrid%weight(o)*molRotGrid%weight(p)*2.0_dp
          END DO  !psi
        END DO   !omega
     
      END DO    !k
    END DO    !j
  END DO    !i
END DO     !species
!========================================================================================================================
!						Get gradient in real space
!========================================================================================================================
do species=1 , nb_species
  do o=1,angGrid%n_angles
    do p=1, molRotGrid%n_angles
    fftw3%in_backward= (dF_pol_trans_k (:,:,:,o,p,species)+dF_pol_long_k (:,:,:,o,p,species)+dF_pol_tot_k (:,:,:,o,p,species))
    call dfftw_execute (fftw3%plan_backward)
    dF_pol_tot (:,:,:,o,p,species)=fftw3%out_backward*deltaVk/(twopi)**3
    END DO
  END DO
END DO
icg=0
!========================================================================================================================
!						Allocate it for minimizing
!========================================================================================================================
do species=1, nb_species
  do i=1,nfft1
    do j=1, nfft2
      do k=1,nfft3
    
        do o=1,angGrid%n_angles
          do p=1, molRotGrid%n_angles
          icg=icg+1
          dF(icg)=dF(icg)+  dF_pol_tot (i,j,k,o,p,species)*cg_vect(icg)*rho_0*2.0_dp * deltav!* angGrid%weight(o) * molRotGrid%weight(p) 
          END DO
 
        END DO
      END DO
    END DO
  END DO
END DO
!========================================================================================================================
!					Check if Polarization Free energy is Real
!========================================================================================================================
if (aimag(F_pol_tot_k+F_pol_long_k+F_pol_trans_k)<tiny(0.0_dp)) then 
    F_pol_tot=real( F_pol_tot_k , dp)
    F_pol_long=real( F_pol_long_k,  dp)
    F_pol_trans=real( F_pol_trans_k, dp)
ELSE
    print*, 'Error in energy_polarization_myway Free energy is not Real'
    print*,aimag(F_pol_tot_k+F_pol_long_k+F_pol_trans_k)
    print*,F_pol_tot_k+F_pol_long_k+F_pol_trans_k
    stop
END IF
!========================================================================================================================
F_pol=F_pol_tot+F_pol_long+F_pol_trans
FF=FF+F_pol
! stop timer
call cpu_time ( time1 )
IF (verbose) PRINT*, 'F_polarization =' , F_pol  , 'computed in (sec)' , time1 - time0
END SUBROUTINE

SUBROUTINE energy_polarization_multi (F_pol)

    USE precision_kinds, ONLY: i2b , dp
    USE system, ONLY: nfft1 , nfft2 , nfft3 , Lx , Ly , Lz , kBT , rho_0 , &
                    deltav, molec_polarx_k,molec_polary_k, molec_polarz_k, kBT,&
                    rho_0_multispec, nb_species, deltax, n_0, beta,deltax,deltay,deltaz, spaceGrid
    USE dcf, ONLY: chi_l, chi_t, c_s, nb_k, delta_k
    use quadrature,only : Omx , Omy , Omz, sym_order , angGrid, molRotGrid
    USE minimizer, ONLY: cg_vect , FF , dF
    use constants,only : twopi, i_complex, fourpi, eps0,qunit,Navo, qfact
    use fft,only : fftw3, kx, ky, kz, k2, norm_k
    use input,only : input_line,input_log, input_char, verbose

    IMPLICIT NONE

    REAL(dp), INTENT(OUT) :: F_pol
    REAL(dp), DIMENSION (nfft1, nfft2, nfft3,angGrid%n_angles,molRotGrid%n_angles,nb_species) ::dF_pol_long,dF_pol_trans,dF_pol_tot
    COMPLEX(dp), DIMENSION (nfft1/2+1, nfft2, nfft3, angGrid%n_angles, molRotGrid%n_angles, nb_species) :: rho_k, dF_pol_long_k ,&
                        dF_pol_trans_k, dF_pol_tot_k
    REAL(dp), DIMENSION (nfft1, nfft2, nfft3, angGrid%n_angles, molRotGrid%n_angles, nb_species) ::rho
    REAL(dp), DIMENSION (nfft1, nfft2, nfft3, nb_species) ::Px,Py,Pz,pola_tot_x,pola_tot_y,pola_tot_z,P_long_x,P_long_y,P_long_z
    INTEGER(i2b) :: icg   !dummy counter for cg_vect
    REAL(dp) ::  deltaVk, F_pol_long, F_pol_trans , F_pol_tot  !Longitudinal , transverse and total Polarization free energy
    COMPLEX(dp), ALLOCATABLE, DIMENSION (:,:,:,:) :: P_trans_x_k,P_trans_y_k,P_trans_z_k,P_long_x_k,P_long_y_k,P_long_z_k
    INTEGER(i2b) :: i , j , k, o, p , n , species, k_index , m1, m2, m3!dummy
    REAL(dp) :: mu_SPCE, facsym
    COMPLEX(dp) :: k_tens_k_Px,k_tens_k_Py,k_tens_k_Pz     
    REAL(dp) :: time1 , time0 ,time2 , time3 ,rhot! timestamps
    COMPLEX(dp) ::  F_pol_long_k , F_pol_trans_k , F_pol_tot_k 
    COMPLEX(dp), ALLOCATABLE, DIMENSION (:,:,:,:) :: pola_tot_x_k , pola_tot_y_k , pola_tot_z_k
    COMPLEX(dp), DIMENSION (nfft1/2+1, nfft2, nfft3) :: rho_c_k_myway

if (nb_species/=1) then
    PRINT*, 'transv_and_longi_polarization_micro IS NOT WORKING FOR MULTISPECIES'
    STOP
END IF


IF ( (SIZE(c_s) /= SIZE(chi_l)) .OR. (SIZE(c_s)/=SIZE(chi_t))) THEN
    WRITE(*,*)"c_s, chi_l and chi_t should have the same number of points, at least for now"
    STOP
END IF

    rho_c_k_myway = (0._dp,0._dp)


call cpu_time(time0)
!            ====================================================
!            !    	Initialization				!
!            !							!
!            ====================================================
deltaVk=twopi**3/(Lx*Ly*Lz)
allocate (P_trans_x_k (nfft1/2+1, nfft2, nfft3, nb_species), SOURCE=(0.0_dp,0.0_dp) )
allocate (P_trans_y_k (nfft1/2+1, nfft2, nfft3, nb_species), SOURCE=(0.0_dp,0.0_dp) )
allocate (P_trans_z_k (nfft1/2+1, nfft2, nfft3, nb_species), SOURCE=(0.0_dp,0.0_dp) )
allocate (P_long_x_k (nfft1/2+1, nfft2, nfft3, nb_species), SOURCE=(0.0_dp,0.0_dp) )
allocate (P_long_y_k (nfft1/2+1, nfft2, nfft3, nb_species), SOURCE=(0.0_dp,0.0_dp) )
allocate (P_long_z_k (nfft1/2+1, nfft2, nfft3, nb_species), SOURCE=(0.0_dp,0.0_dp) )
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
dF_pol_long_k=(0.0_dp,0.0_dp)
dF_pol_trans_k=(0.0_dp,0.0_dp)
dF_pol_tot_k=(0.0_dp,0.0_dp)
icg=0
rho=0.0_dp
rho_k=(0.0_dp,0.0_dp)

ALLOCATE (pola_tot_x_k ( nfft1/2+1,nfft2,nfft3,nb_species), SOURCE=(0.0_dp,0.0_dp) )
ALLOCATE (pola_tot_y_k ( nfft1/2+1,nfft2,nfft3,nb_species), SOURCE=(0.0_dp,0.0_dp) )
ALLOCATE (pola_tot_z_k ( nfft1/2+1,nfft2,nfft3,nb_species), SOURCE=(0.0_dp,0.0_dp) )
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
      P_long_x_k(i,j,k,n)=(0.0_dp,0.0_dp)
      P_long_y_k(i,j,k,n)=(0.0_dp,0.0_dp)
      P_long_z_k(i,j,k,n)=(0.0_dp,0.0_dp)
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
    DO species=1, nb_species
        DO i=1,nfft1
            DO j=1, nfft2
                DO k=1,nfft3
                    DO o=1,angGrid%n_angles
                        DO p=1, molRotGrid%n_angles
                            icg = icg + 1
                            dF(icg)=dF(icg)+  dF_pol_tot (i,j,k,o,p,species)*cg_vect(icg)*rho_0*2.0_dp * deltav!* angGrid%weight(o) * molRotGrid%weight(p) 
                        END DO
                    END DO
                END DO
            END DO
        END DO
    END DO
!===================================================================================================================================
!					Check if Polarization Free energy is Real
!===================================================================================================================================
    IF ( AIMAG(F_pol_tot_k)==0._dp .AND. AIMAG(F_pol_long_k)==0._dp .AND. AIMAG(F_pol_trans_k)==0._dp ) THEN
        F_pol_tot   = REAL( F_pol_tot_k , dp)
        F_pol_long  = REAL( F_pol_long_k,  dp)
        F_pol_trans = REAL( F_pol_trans_k, dp)
    ELSE
        WRITE(*,*) 'Error in energy_polarization_myway Free energy is not Real'
        WRITE(*,*) 'Imaginary part of F_pol_tot_k, F_pol_long_k and F_pol_trans_k:'
        WRITE(*,*) AIMAG(F_pol_tot_k), AIMAG(F_pol_long_k), AIMAG(F_pol_trans_k)
        STOP
    END IF
!===================================================================================================================================

    F_pol = F_pol_tot + F_pol_long + F_pol_trans
    FF = FF + F_pol

    CALL CPU_TIME ( time1 )
    IF (verbose) PRINT*, 'F_polarization =' , F_pol  , 'computed in (sec)' , time1 - time0

END SUBROUTINE energy_polarization_multi

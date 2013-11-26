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
    REAL(dp) ::  deltaVk, F_pol_long, F_pol_trans , F_pol_tot  !Longitudinal , transverse and total Polarization free energy
    REAL(dp) :: mu_SPCE, facsym
    REAL(dp) :: time1 , time0 ,time2 , time3 ,rhot! timestamps
    REAL(dp), DIMENSION (nfft1, nfft2, nfft3,angGrid%n_angles,molRotGrid%n_angles,nb_species) ::dF_pol_long,dF_pol_trans,dF_pol_tot

    COMPLEX(dp), DIMENSION (nfft1/2+1, nfft2, nfft3, angGrid%n_angles, molRotGrid%n_angles, nb_species) :: rho_k, dF_pol_long_k ,&
                        dF_pol_trans_k, dF_pol_tot_k
    REAL(dp), DIMENSION (nfft1, nfft2, nfft3, angGrid%n_angles, molRotGrid%n_angles, nb_species) ::rho
    REAL(dp), DIMENSION (nfft1, nfft2, nfft3, nb_species) ::Px,Py,Pz,pola_tot_x,pola_tot_y,pola_tot_z,P_long_x,P_long_y,P_long_z
    INTEGER(i2b) :: icg   !dummy counter for cg_vect

    COMPLEX(dp), ALLOCATABLE, DIMENSION (:,:,:,:) :: P_trans_x_k,P_trans_y_k,P_trans_z_k,P_long_x_k,P_long_y_k,P_long_z_k
    INTEGER(i2b) :: i , j , k, o, p , n , s, k_index , m1, m2, m3!dummy

    COMPLEX(dp) :: k_tens_k_Px, k_tens_k_Py, k_tens_k_Pz, toto
    COMPLEX(dp) :: F_pol_long_k, F_pol_trans_k, F_pol_tot_k 
    COMPLEX(dp), ALLOCATABLE, DIMENSION (:,:,:,:) :: pola_tot_x_k , pola_tot_y_k , pola_tot_z_k
    COMPLEX(dp), DIMENSION (nfft1/2+1, nfft2, nfft3) :: rho_c_k_myway

    IF (nb_species/=1) STOP 'transv_and_longi_polarization_micro IS NOT WORKING FOR MULTISPECIES'

    IF ( (SIZE(c_s) /= SIZE(chi_l)) .OR. (SIZE(c_s)/=SIZE(chi_t))) THEN
        WRITE(*,*)"c_s, chi_l and chi_t should have the same number of points, at least for now"
        STOP
    END IF

    rho_c_k_myway = (0._dp,0._dp)

    CALL CPU_TIME (time0)
!===================================================================================================================================
!                	Initialization				
!            							
!===================================================================================================================================
    deltaVk=twopi**3/(Lx*Ly*Lz)
    ALLOCATE (P_trans_x_k (nfft1/2+1, nfft2, nfft3, nb_species), SOURCE=(0.0_dp,0.0_dp) )
    ALLOCATE (P_trans_y_k (nfft1/2+1, nfft2, nfft3, nb_species), SOURCE=(0.0_dp,0.0_dp) )
    ALLOCATE (P_trans_z_k (nfft1/2+1, nfft2, nfft3, nb_species), SOURCE=(0.0_dp,0.0_dp) )
    ALLOCATE (P_long_x_k  (nfft1/2+1, nfft2, nfft3, nb_species), SOURCE=(0.0_dp,0.0_dp) )
    ALLOCATE (P_long_y_k  (nfft1/2+1, nfft2, nfft3, nb_species), SOURCE=(0.0_dp,0.0_dp) )
    ALLOCATE (P_long_z_k  (nfft1/2+1, nfft2, nfft3, nb_species), SOURCE=(0.0_dp,0.0_dp) )
    ALLOCATE (pola_tot_x_k (nfft1/2+1, nfft2, nfft3, nb_species), SOURCE=(0.0_dp,0.0_dp) )
    ALLOCATE (pola_tot_y_k (nfft1/2+1, nfft2, nfft3, nb_species), SOURCE=(0.0_dp,0.0_dp) )
    ALLOCATE (pola_tot_z_k (nfft1/2+1, nfft2, nfft3, nb_species), SOURCE=(0.0_dp,0.0_dp) )
    F_pol_long = 0.0_dp
    F_pol_trans = 0.0_dp
    F_pol_tot = 0.0_dp
    F_pol_long_k = (0.0_dp,0.0_dp)
    F_pol_trans_k = (0.0_dp,0.0_dp)
    F_pol_tot_k = (0.0_dp,0.0_dp)
    F_pol = 0.0_dp
    dF_pol_long = 0.0_dp
    dF_pol_trans = 0.0_dp
    dF_pol_tot = 0.0_dp
    dF_pol_long_k = (0.0_dp,0.0_dp)
    dF_pol_trans_k = (0.0_dp,0.0_dp)
    dF_pol_tot_k = (0.0_dp,0.0_dp)
    rho = 0.0_dp
    rho_k = (0.0_dp,0.0_dp)
    rho_c_k_myway = 0._dp

    mu_SPCE=0.4238_dp*0.5773525_dp*2.0_dp !dipolar moment of SPCE water molecule in e.Angstromm
!======================================================================================================
!            ====================================================
!            !    	Compute density in real and 		!
             !       		Fourier space			!
             !							!
!            ====================================================
!Compute rho_k
    icg=0
    DO s =1 , nb_species
        DO i=1,nfft1
            DO j=1, nfft2
                DO k=1, nfft3
                    DO o = 1 , angGrid%n_angles
                        DO p=1, molRotGrid%n_angles
                            icg = icg + 1
                            rho(i,j,k,o,p,s) = cg_vect(icg)**2
                        END DO
                    END DO
                END DO
            END DO
        END DO
    END DO

    CALL CPU_TIME (time2)
!===================================================================================================================================
!            !    		get Density 			!
!            !			in Fourier Space		!
!===================================================================================================================================
    DO s = 1, nb_species
        DO o = 1, angGrid%n_angles
            DO p = 1, molRotGrid%n_angles
                fftw3%in_forward = rho(:,:,:,o,p,s)
                CALL dfftw_execute (fftw3%plan_forward)
                rho_k(:,:,:,o,p,s) = fftw3%out_forward*deltaV
            END DO
        END DO
    END DO

    CALL CPU_TIME (time3)
!===================================================================================================================================
!            !    		Compute 			
!            !		Total Polarization			
!===================================================================================================================================
    
    DO CONCURRENT ( s=1:nb_species, i=1:nfft1/2+1, j=1:nfft2, k=1:nfft3, o=1:angGrid%n_angles, p=1:molRotGrid%n_angles )
        pola_tot_x_k(i,j,k,s) = &
            pola_tot_x_k(i,j,k,s)+angGrid%weight(o)*molRotGrid%weight(p) * molec_polarx_k(i,j,k,o,p,s)*rho_k(i,j,k,o,p,s)
        pola_tot_y_k(i,j,k,s) = &
            pola_tot_y_k(i,j,k,s)+angGrid%weight(o)*molRotGrid%weight(p) * molec_polary_k(i,j,k,o,p,s)*rho_k(i,j,k,o,p,s)
        pola_tot_z_k(i,j,k,s) = &
            pola_tot_z_k(i,j,k,s)+angGrid%weight(o)*molRotGrid%weight(p) * molec_polarz_k(i,j,k,o,p,s)*rho_k(i,j,k,o,p,s)
    END DO


    DO s=1,nb_species
        fftw3%in_backward= pola_tot_x_k(:,:,:,s)
        call dfftw_execute(fftw3%plan_backward)
        pola_tot_x(:,:,:,s)=fftw3%out_backward*deltaVk/(twopi)**3
    
        fftw3%in_backward= pola_tot_y_k(:,:,:,s)
        call dfftw_execute(fftw3%plan_backward)
        pola_tot_y(:,:,:,s)=fftw3%out_backward*deltaVk/(twopi)**3
    
        fftw3%in_backward= pola_tot_z_k(:,:,:,s)
        call dfftw_execute(fftw3%plan_backward)
        pola_tot_z(:,:,:,s)=fftw3%out_backward*deltaVk/(twopi)**3
    END DO

!            ====================================================
!            !    	Compute 				!
!            !	Transverse and longitudinal Polarization	!
!            ====================================================
    DO CONCURRENT ( s=1:nb_species, i=1:nfft1/2+1, j=1:nfft2, k=1:nfft3 )
        IF (k2(i,j,k)==0.0_dp) THEN
            P_long_x_k(i,j,k,s) = (0.0_dp,0.0_dp)
            P_long_y_k(i,j,k,s) = (0.0_dp,0.0_dp)
            P_long_z_k(i,j,k,s) = (0.0_dp,0.0_dp)
        ELSE
            toto = (pola_tot_x_k(i,j,k,s)*kx(i)+pola_tot_y_k(i,j,k,s)*ky(j)+pola_tot_z_k(i,j,k,s)*kz(k)) / k2(i,j,k)
            P_long_x_k(i,j,k,s) = toto * kx(i)
            P_long_y_k(i,j,k,s) = toto * ky(j)
            P_long_z_k(i,j,k,s) = toto * kz(k)
        END IF
    END DO
    P_trans_x_k = pola_tot_x_k - P_long_x_k
    P_trans_y_k = pola_tot_y_k - P_long_y_k
    P_trans_z_k = pola_tot_z_k - P_long_z_k



    DO s=1,nb_species
        fftw3%in_backward= P_long_x_k(:,:,:,s)
        CALL dfftw_execute(fftw3%plan_backward)
        P_long_x(:,:,:,s)=fftw3%out_backward*deltaVk/(twopi)**3
        
        fftw3%in_backward= P_long_y_k(:,:,:,s)
        CALL dfftw_execute(fftw3%plan_backward)
        P_long_y(:,:,:,s)=fftw3%out_backward*deltaVk/(twopi)**3
        
        fftw3%in_backward= P_long_z_k(:,:,:,s)
        CALL dfftw_execute(fftw3%plan_backward)
        P_long_z(:,:,:,s)=fftw3%out_backward*deltaVk/(twopi)**3
    END DO

!            ====================================================
!            !    	Compute Free energy due to		!
!            !	Transverse and longitudinal Polarization	!
!            ====================================================
    DO s=1, nb_species 
        DO i=1, nfft1/2+1
            
            IF (i>1 .and. i<nfft1/2+1) THEN
                facsym=2.0_dp
            ELSE
                facsym=1.0_dp
            END IF
            
            DO j=1, nfft2
                DO k=1, nfft3
                    
                    k_index = INT ( norm_k(i,j,k) / delta_k ) + 1
                    IF ( k_index > nb_k ) k_index = nb_k ! Here it happens that k_index gets higher than the highest c_k index. In this case one imposes k_index = k_index_max
       
                    F_pol_long_k = &
                        F_pol_long_k+deltaVk*fourpi*(P_long_x_k(i,j,k,s)*CONJG(p_long_x_k(i,j,k,s))+&
                        p_long_y_k(i,j,k,s)*CONJG(p_long_y_k(i,j,k,s))+p_long_z_k(i,j,k,s)*CONJG(p_long_z_k(i,j,k,s)))&
                        /chi_l(k_index)*0.5_dp*qfact*rho_0**2*facsym/(twopi**3)

                    F_pol_trans_k = &
                        F_pol_trans_k+deltaVk*(P_trans_x_k(i,j,k,s)*CONJG(p_trans_x_k(i,j,k,s))+p_trans_y_k(i,j,k,s)&
                        *CONJG(p_trans_y_k(i,j,k,s))+p_trans_z_k(i,j,k,s)*CONJG(p_trans_z_k(i,j,k,s)))/chi_t(k_index)&
                        *0.5_dp*qfact*rho_0**2*facsym/(twopi**3)
       
                    F_pol_tot_k = &
                        F_pol_tot_k-kBT*3/(2*mu_SPCE**2*n_0)*deltaVk*facsym/(twopi**3)*rho_0**2*&
                        (pola_tot_x_k(i,j,k,s)*CONJG(pola_tot_x_k(i,j,k,s))+pola_tot_y_k(i,j,k,s)*&
                        CONJG(pola_tot_y_k(i,j,k,s))+pola_tot_z_k(i,j,k,s)*CONJG(pola_tot_z_k(i,j,k,s)))

                    DO o=1,angGrid%n_angles
                        DO p=1,molRotGrid%n_angles 
!========================================================================================================================
!Evaluate gradient
!========================================================================================================================
                            toto = CONJG(kx(i)*molec_polarx_k(i,j,k,o,p,s)&
                                +ky(j)*molec_polary_k(i,j,k,o,p,s)+kz(k)*molec_polarz_k(i,j,k,o,p,s))
                            k_tens_k_Px = kx(i) * toto
                            k_tens_k_Py = ky(j) * toto
                            k_tens_k_Pz = kz(k) * toto
   
                            IF ( k2(i,j,k)==0.0_dp ) THEN
                                dF_pol_long_k(i,j,k,o,p,s)=0.0_dp
                                dF_pol_trans_k(i,j,k,o,p,s)=rho_0*0.5_dp*qfact/chi_t(k_index)*&
                                    (P_trans_x_k(i,j,k,s)*CONJG(molec_polarx_k(i,j,k,o,p,s))&
                                    +P_trans_y_k(i,j,k,s)*CONJG(molec_polary_k(i,j,k,o,p,s))&
                                    +P_trans_z_k(i,j,k,s)*CONJG(molec_polarz_k(i,j,k,o,p,s)))&
                                    *angGrid%weight(o)*molRotGrid%weight(p)*2.0_dp
                            ELSE
                                dF_pol_long_k(i,j,k,o,p,s) = rho_0*0.5_dp*&
                                    qfact/chi_l(k_index)*fourpi*angGrid%weight(o)*molRotGrid%weight(p)*&
                                    (P_long_x_k(i,j,k,s)*k_tens_k_Px+P_long_y_k(i,j,k,s)*k_tens_k_Py+P_long_z_k(i,j,k,s)&
                                    *k_tens_k_Pz)/k2(i,j,k)*2.0_dp
                                dF_pol_trans_k(i,j,k,o,p,s) = rho_0*0.5_dp*qfact/chi_t(k_index)*(P_trans_x_k(i,j,k,s)*&
                                    (CONJG(molec_polarx_k(i,j,k,o,p,s))-k_tens_k_Px/k2(i,j,k))&
                                    +P_trans_y_k(i,j,k,s)*(CONJG(molec_polary_k(i,j,k,o,p,s))-k_tens_k_Py/k2(i,j,k))&
                                    +P_trans_z_k(i,j,k,s)*(CONJG(molec_polarz_k(i,j,k,o,p,s))-k_tens_k_Pz/k2(i,j,k)) )&
                                    *angGrid%weight(o)*molRotGrid%weight(p)*2.0_dp
                            END IF
                            
                            dF_pol_tot_k(i,j,k,o,p,s)=-kBT*3*rho_0/(2*mu_SPCE**2*n_0)*(pola_tot_x_k(i,j,k,s)*&
                                CONJG(molec_polarx_k(i,j,k,o,p,s))+pola_tot_y_k(i,j,k,s)*CONJG(molec_polary_k(i,j,k,o,p,s))&
                                +pola_tot_z_k(i,j,k,s)*CONJG(molec_polarz_k(i,j,k,o,p,s)))&
                                *angGrid%weight(o)*molRotGrid%weight(p)*2.0_dp
                        END DO
                    END DO
                END DO
            END DO
        END DO
    END DO
!========================================================================================================================
!						Get gradient in real space
!========================================================================================================================
    DO s=1,nb_species
        DO o=1,angGrid%n_angles
            DO p=1, molRotGrid%n_angles
                fftw3%in_backward= (dF_pol_trans_k (:,:,:,o,p,s)+dF_pol_long_k (:,:,:,o,p,s)+dF_pol_tot_k (:,:,:,o,p,s))
                call dfftw_execute (fftw3%plan_backward)
                dF_pol_tot (:,:,:,o,p,s)=fftw3%out_backward*deltaVk/(twopi)**3
            END DO
        END DO
    END DO
!========================================================================================================================
!						Allocate it for minimizing
!========================================================================================================================
    icg = 0
    DO s=1, nb_species
        DO i=1,nfft1
            DO j=1, nfft2
                DO k=1,nfft3
                    DO o=1,angGrid%n_angles
                        DO p=1, molRotGrid%n_angles
                            icg = icg + 1
                            dF(icg)=dF(icg)+ dF_pol_tot (i,j,k,o,p,s)*cg_vect(icg)*rho_0*2.0_dp * deltav!* angGrid%weight(o) * molRotGrid%weight(p) 
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

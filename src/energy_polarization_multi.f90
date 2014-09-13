SUBROUTINE energy_polarization_multi (F_pol)

    USE precision_kinds, ONLY: i2b, dp
    USE system, ONLY: thermocond, molec_polarx_k, molec_polary_k, molec_polarz_k, nb_species, spaceGrid, solvent
    USE dcf, ONLY: chi_l, chi_t, c_s, nb_k, delta_k
    USE quadrature,ONLY : angGrid, molRotGrid
    USE minimizer, ONLY: cg_vect, FF, dF
    USE constants, ONLY : twopi, fourpi, qfact
    USE fft, ONLY: fftw3, kx, ky, kz, k2, norm_k
    USE input, ONLY: verbose

    IMPLICIT NONE

    REAL(dp), INTENT(OUT) :: F_pol
    REAL(dp) :: F_pol_long, F_pol_trans , F_pol_tot  !Longitudinal , transverse and total Polarization free energy
    REAL(dp) :: mu_SPCE, facsym, deltaVk, Lweight
    REAL(dp), ALLOCATABLE, DIMENSION (:,:,:,:,:,:) :: rho, dF_pol_tot
    COMPLEX(dp), ALLOCATABLE, DIMENSION (:,:,:,:,:,:) :: rho_k, dF_pol_tot_k,  dF_pol_long_k,   dF_pol_trans_k
    COMPLEX(dp), ALLOCATABLE, DIMENSION (:,:,:,:) :: &
        P_trans_x_k, P_trans_y_k, P_trans_z_k, P_long_x_k, P_long_y_k, P_long_z_k, pola_tot_x_k, pola_tot_y_k, pola_tot_z_k
    INTEGER(i2b) :: icg, i, j, k, o, p, n, s, k_index
    COMPLEX(dp) :: k_tens_k_Px, k_tens_k_Py, k_tens_k_Pz, toto, t_
    COMPLEX(dp), PARAMETER :: zeroC = (0.0_dp,0.0_dp)
    INTEGER(i2b) :: nfft1, nfft2, nfft3

    nfft1= spaceGrid%n_nodes(1)
    nfft2= spaceGrid%n_nodes(2)
    nfft3= spaceGrid%n_nodes(3)

    IF (nb_species/=1) STOP 'transv_and_longi_polarization_micro IS NOT WORKING FOR MULTISPECIES'

    IF ( (SIZE(c_s) /= SIZE(chi_l)) .OR. (SIZE(c_s)/=SIZE(chi_t))) THEN
        WRITE(*,*)"c_s, chi_l and chi_t should have the same number of points, at least for now"
        PRINT*, SIZE(c_s), SIZE(chi_l), SIZE(chi_t)
        STOP
    END IF

!===================================================================================================================================
!                	Initialization				
!            							
!===================================================================================================================================
    deltaVk = twopi**3/PRODUCT(spaceGrid%length)
    F_pol = 0.0_dp
    F_pol_long = 0.0_dp
    F_pol_trans = 0.0_dp
    F_pol_tot = 0.0_dp
    mu_SPCE = 0.4238_dp * 0.5773525_dp * 2.0_dp !dipolar moment of SPCE water molecule in e.Angstromm
!======================================================================================================
!            ====================================================
!                 	Compute density in real and 		        !
!       		            Fourier space		            	!
!							                                    !
!            ====================================================
    ALLOCATE ( rho (spaceGrid%n_nodes(1), spaceGrid%n_nodes(2), spaceGrid%n_nodes(3),&
                        angGrid%n_angles, molRotGrid%n_angles, nb_species), SOURCE=0._dp )
    icg=0
    DO s=1, nb_species
        DO i=1, nfft1
            DO j=1, nfft2
                DO k=1, nfft3
                    DO o=1, angGrid%n_angles
                        DO p=1, molRotGrid%n_angles
                            icg = icg + 1
                            rho(i,j,k,o,p,s) = cg_vect(icg)**2
                        END DO
                    END DO
                END DO
            END DO
        END DO
    END DO

!===================================================================================================================================
!            !    		get Density 			!
!            !			in Fourier Space		!
!===================================================================================================================================
    ALLOCATE (rho_k (nfft1/2+1, nfft2, nfft3, angGrid%n_angles, molRotGrid%n_angles, nb_species), SOURCE=(0._dp,0._dp))
    DO s = 1, nb_species
        DO p = 1, molRotGrid%n_angles
            DO o = 1, angGrid%n_angles
                fftw3%in_forward = rho(:,:,:,o,p,s)
                CALL dfftw_execute (fftw3%plan_forward)
                rho_k(:,:,:,o,p,s) = fftw3%out_forward *spacegrid%dv
            END DO
        END DO
    END DO
    DEALLOCATE (rho)

!===================================================================================================================================
!            !    		Compute 			
!            !		Total Polarization			
!===================================================================================================================================
    ALLOCATE (pola_tot_x_k (nfft1/2+1, nfft2, nfft3, nb_species), SOURCE=zeroC )
    ALLOCATE (pola_tot_y_k (nfft1/2+1, nfft2, nfft3, nb_species), SOURCE=zeroC )
    ALLOCATE (pola_tot_z_k (nfft1/2+1, nfft2, nfft3, nb_species), SOURCE=zeroC )    
    DO CONCURRENT ( i=1:nfft1/2+1, j=1:nfft2, k=1:nfft3, o=1:angGrid%n_angles, p=1:molRotGrid%n_angles, s=1:nb_species )
        t_ = angGrid%weight(o)*molRotGrid%weight(p)
        toto = rho_k(i,j,k,o,p,s)
        pola_tot_x_k(i,j,k,s) = pola_tot_x_k(i,j,k,s)+t_ * molec_polarx_k(i,j,k,o,p,s)*toto
        pola_tot_y_k(i,j,k,s) = pola_tot_y_k(i,j,k,s)+t_ * molec_polary_k(i,j,k,o,p,s)*toto
        pola_tot_z_k(i,j,k,s) = pola_tot_z_k(i,j,k,s)+t_ * molec_polarz_k(i,j,k,o,p,s)*toto
    END DO
    DEALLOCATE (rho_k)

!            ====================================================
!            !    	Compute 				!
!            !	Transverse and longitudinal Polarization	!
!            ====================================================
    ALLOCATE (P_long_x_k  (nfft1/2+1, nfft2, nfft3, nb_species), SOURCE=zeroC )
    ALLOCATE (P_long_y_k  (nfft1/2+1, nfft2, nfft3, nb_species), SOURCE=zeroC )
    ALLOCATE (P_long_z_k  (nfft1/2+1, nfft2, nfft3, nb_species), SOURCE=zeroC )
    DO CONCURRENT ( i=1:nfft1/2+1, j=1:nfft2, k=1:nfft3, s=1:nb_species )
        IF (k2(i,j,k)==0.0_dp) THEN
            P_long_x_k(i,j,k,s) = zeroC
            P_long_y_k(i,j,k,s) = zeroC
            P_long_z_k(i,j,k,s) = zeroC
        ELSE
            toto = (pola_tot_x_k(i,j,k,s)*kx(i)+pola_tot_y_k(i,j,k,s)*ky(j)+pola_tot_z_k(i,j,k,s)*kz(k)) / k2(i,j,k)
            P_long_x_k(i,j,k,s) = toto * kx(i)
            P_long_y_k(i,j,k,s) = toto * ky(j)
            P_long_z_k(i,j,k,s) = toto * kz(k)
        END IF
    END DO
    ALLOCATE (P_trans_x_k (nfft1/2+1, nfft2, nfft3, nb_species), SOURCE=(pola_tot_x_k - P_long_x_k) )
    ALLOCATE (P_trans_y_k (nfft1/2+1, nfft2, nfft3, nb_species), SOURCE=(pola_tot_y_k - P_long_y_k) )
    ALLOCATE (P_trans_z_k (nfft1/2+1, nfft2, nfft3, nb_species), SOURCE=(pola_tot_z_k - P_long_z_k) )

!            ====================================================
!            !    	Compute Free energy due to		!
!            !	Transverse and longitudinal Polarization	!
!            ====================================================
    ALLOCATE( dF_pol_tot (nfft1, nfft2, nfft3, angGrid%n_angles, molRotGrid%n_angles, nb_species), source=0._dp)
    ALLOCATE( dF_pol_tot_k (nfft1/2+1, nfft2, nfft3, angGrid%n_angles, molRotGrid%n_angles, nb_species), source=zeroC)
    ALLOCATE( dF_pol_long_k (nfft1/2+1, nfft2, nfft3, angGrid%n_angles, molRotGrid%n_angles, nb_species), source=zeroC)
    ALLOCATE( dF_pol_trans_k (nfft1/2+1, nfft2, nfft3, angGrid%n_angles, molRotGrid%n_angles, nb_species), source=zeroC)
    DO CONCURRENT ( i=1:nfft1/2+1, j=1:nfft2, k=1:nfft3, s=1:nb_species )
         
        IF (i>1 .AND. i<nfft1/2+1) THEN
            facsym=2.0_dp
        ELSE
            facsym=1.0_dp
        END IF
        
        k_index = INT ( norm_k(i,j,k) / delta_k ) + 1
        IF ( k_index > nb_k ) k_index = nb_k ! Here it happens that k_index gets higher than the highest c_k index. In this case one imposes k_index = k_index_max
        
        F_pol_long = &
            F_pol_long+deltaVk*fourpi*REAL(   P_long_x_k(i,j,k,s)*CONJG(P_long_x_k(i,j,k,s))&
                                            + P_long_y_k(i,j,k,s)*CONJG(P_long_y_k(i,j,k,s))&
                                            + P_long_z_k(i,j,k,s)*CONJG(P_long_z_k(i,j,k,s)))&
                                        /chi_l(k_index)*0.5_dp*qfact*solvent(s)%rho0**2*facsym/(twopi**3)
        
        F_pol_trans = &
            F_pol_trans+deltaVk*REAL(  P_trans_x_k(i,j,k,s)*CONJG(P_trans_x_k(i,j,k,s))&
                                      +P_trans_y_k(i,j,k,s)*CONJG(P_trans_y_k(i,j,k,s))&
                                      +P_trans_z_k(i,j,k,s)*CONJG(P_trans_z_k(i,j,k,s)))&
                                        /chi_t(k_index)*0.5_dp*qfact*solvent(s)%rho0**2*facsym/(twopi**3)
        
        F_pol_tot = &
            F_pol_tot-thermocond%kbT*3/(2*mu_SPCE**2*solvent(1)%n0)*deltaVk*facsym/(twopi**3)*solvent(s)%rho0**2*&
            REAL(pola_tot_x_k(i,j,k,s)*CONJG(pola_tot_x_k(i,j,k,s))&
                +pola_tot_y_k(i,j,k,s)*CONJG(pola_tot_y_k(i,j,k,s))&
                +pola_tot_z_k(i,j,k,s)*CONJG(pola_tot_z_k(i,j,k,s)))
        
        DO CONCURRENT ( o=1:angGrid%n_angles, p=1:molRotGrid%n_angles )
            Lweight = angGrid%weight(o) * molRotGrid%weight(p)
!========================================================================================================================
!Evaluate gradient
!========================================================================================================================
            IF ( k2(i,j,k)<=tiny(1.0_dp) ) THEN
                dF_pol_long_k(i,j,k,o,p,s)=0.0_dp
                dF_pol_trans_k(i,j,k,o,p,s)=solvent(s)%rho0*0.5_dp*qfact/chi_t(k_index)*&
                    (P_trans_x_k(i,j,k,s)*CONJG(molec_polarx_k(i,j,k,o,p,s))&
                    +P_trans_y_k(i,j,k,s)*CONJG(molec_polary_k(i,j,k,o,p,s))&
                    +P_trans_z_k(i,j,k,s)*CONJG(molec_polarz_k(i,j,k,o,p,s)))&
                    *Lweight*2.0_dp
            ELSE
                toto = CONJG(kx(i)*molec_polarx_k(i,j,k,o,p,s)&
                            +ky(j)*molec_polary_k(i,j,k,o,p,s)&
                            +kz(k)*molec_polarz_k(i,j,k,o,p,s))
                k_tens_k_Px = kx(i) * toto
                k_tens_k_Py = ky(j) * toto
                k_tens_k_Pz = kz(k) * toto

                dF_pol_long_k(i,j,k,o,p,s) = solvent(s)%rho0*0.5_dp*&
                    qfact/chi_l(k_index)*fourpi*Lweight*&
                    (P_long_x_k(i,j,k,s)*k_tens_k_Px&
                    +P_long_y_k(i,j,k,s)*k_tens_k_Py&
                    +P_long_z_k(i,j,k,s)*k_tens_k_Pz)/k2(i,j,k)*2.0_dp

                dF_pol_trans_k(i,j,k,o,p,s) = solvent(s)%rho0*0.5_dp*qfact/chi_t(k_index)*&
                    (P_trans_x_k(i,j,k,s)*(CONJG(molec_polarx_k(i,j,k,o,p,s))-k_tens_k_Px/k2(i,j,k))&
                    +P_trans_y_k(i,j,k,s)*(CONJG(molec_polary_k(i,j,k,o,p,s))-k_tens_k_Py/k2(i,j,k))&
                    +P_trans_z_k(i,j,k,s)*(CONJG(molec_polarz_k(i,j,k,o,p,s))-k_tens_k_Pz/k2(i,j,k)) )&
                    *Lweight*2.0_dp
            END IF
            
            dF_pol_tot_k(i,j,k,o,p,s)=-thermocond%kbT*3._dp*solvent(s)%rho0/(2._dp*mu_SPCE**2*solvent(1)%n0)*(&
                 pola_tot_x_k(i,j,k,s)*CONJG(molec_polarx_k(i,j,k,o,p,s))&
                +pola_tot_y_k(i,j,k,s)*CONJG(molec_polary_k(i,j,k,o,p,s))&
                +pola_tot_z_k(i,j,k,s)*CONJG(molec_polarz_k(i,j,k,o,p,s)))&
                *Lweight*2.0_dp
        END DO
    END DO

!========================================================================================================================
!						Get gradient in real space
!========================================================================================================================
    DO s=1,nb_species
        DO p=1, molRotGrid%n_angles
            DO o=1,angGrid%n_angles
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
        DO i=1, nfft1
            DO j=1, nfft2
                DO k=1,nfft3
                    DO o=1,angGrid%n_angles
                        DO p=1, molRotGrid%n_angles
                            icg = icg + 1
                            dF(icg) = dF(icg) + dF_pol_tot (i,j,k,o,p,s) * cg_vect(icg) * solvent(s)%rho0 * 2.0_dp * spaceGrid%dv
                        END DO
                    END DO
                END DO
            END DO
        END DO
    END DO

    F_pol = F_pol_tot + F_pol_long + F_pol_trans
    FF = FF + F_pol

    DEALLOCATE (dF_pol_tot, dF_pol_tot_k, dF_pol_long_k, dF_pol_trans_k)
END SUBROUTINE energy_polarization_multi

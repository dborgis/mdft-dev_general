SUBROUTINE energy_polarization_multi (F_pol)

    USE precision_kinds, ONLY: i2b, dp
    USE system, ONLY: thermocond, nb_species, spaceGrid, solvent
    USE dcf, ONLY: chi_l, chi_t, c_s, nb_k, delta_k
    USE quadrature,ONLY : angGrid, molRotGrid
    USE minimizer, ONLY: cg_vect, FF, dF
    USE constants, ONLY : twopi, fourpi, qfact, zeroC
    USE fft, ONLY: fftw3, kx, ky, kz, k2, norm_k
    USE input, ONLY: verbose

    IMPLICIT NONE

    REAL(dp), INTENT(OUT) :: F_pol
    REAL(dp) :: F_pol_long, F_pol_trans , F_pol_tot  ! Longitudinal, transverse and total Polarization free energy
    REAL(dp) :: mu_SPCE, facsym, deltaVkn, Lweight, kvec(3)
    REAL(dp), ALLOCATABLE, DIMENSION (:,:,:,:,:,:) :: rho, dF_pol_tot
    COMPLEX(dp), ALLOCATABLE, DIMENSION (:,:,:,:,:,:) :: rho_k, dF_pol_tot_k,  dF_pol_long_k,   dF_pol_trans_k
    COMPLEX(dp), ALLOCATABLE, DIMENSION (:,:,:,:,:) :: P_trans_k, P_long_k, P_tot_k
    INTEGER(i2b) :: icg, i, j, k, o, p, n, s, k_index
    COMPLEX(dp) :: k_tens_k_P(3), toto, temp, P_trans_k_loc(3), P_long_k_loc(3), P_tot_k_loc(3), molec_polar_k_loc(3)
    integer(i2b), pointer :: nfft1 => spacegrid%n_nodes(1), nfft2 => spacegrid%n_nodes(2), nfft3 => spacegrid%n_nodes(3)

    if (size(solvent)/=1) STOP 'energy_polarization_multi.f90 not working for multisolvent species'

    IF ( (SIZE(c_s) /= SIZE(chi_l)) .OR. (SIZE(c_s)/=SIZE(chi_t))) THEN
        WRITE(*,*)"c_s, chi_l and chi_t should have the same number of points, at least for now"
        PRINT*, SIZE(c_s), SIZE(chi_l), SIZE(chi_t)
        STOP
    END IF

!===================================================================================================================================
!                	Initialization				
!            							
!===================================================================================================================================
    deltaVkn = 1._dp/PRODUCT(spaceGrid%length)
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
    allocate (P_tot_k (3,nfft1/2+1, nfft2, nfft3, nb_species), SOURCE=zeroC )
    do concurrent ( i=1:nfft1/2+1, j=1:nfft2, k=1:nfft3, o=1:angGrid%n_angles, p=1:molRotGrid%n_angles, s=1:nb_species )
        P_tot_k(:,i,j,k,s) = P_tot_k(:,i,j,k,s) &
                    + rho_k(i,j,k,o,p,s) * angGrid%weight(o) * molRotGrid%weight(p) * solvent(s)%molec_polar_k(:,i,j,k,o,p)
    end do
    deallocate(rho_k)


!            ====================================================
!            !    	Compute 				!
!            !	Transverse and longitudinal Polarization	!
!            ====================================================
    allocate (P_long_k (3,nfft1/2+1, nfft2, nfft3, nb_species), SOURCE=zeroC )
    do concurrent ( i=1:nfft1/2+1, j=1:nfft2, k=1:nfft3, s=1:nb_species )
        if (abs(k2(i,j,k))<=epsilon(1._dp)) THEN
            P_long_k(:,i,j,k,s) = zeroC
        else
            kvec = [kx(i),ky(j),kz(k)]
            P_long_k(:,i,j,k,s) = kvec * sum(kvec * P_tot_k(:,i,j,k,s)) / k2(i,j,k)
        end if
    end do

    allocate (P_trans_k (3,nfft1/2+1, nfft2, nfft3, nb_species), SOURCE=(P_tot_k - P_long_k) )

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
        
        P_trans_k_loc = P_trans_k(:,i,j,k,s)
        P_long_k_loc  = P_long_k(:,i,j,k,s)
        P_tot_k_loc   = P_tot_k(:,i,j,k,s)
        
        k_index = INT ( norm_k(i,j,k) / delta_k ) + 1
        IF ( k_index > nb_k ) k_index = nb_k ! Here it happens that k_index gets higher than the highest c_k index. In this case one imposes k_index = k_index_max
        
        toto = deltaVkn*0.5_dp*qfact*solvent(s)%rho0**2*facsym
        
        F_pol_long = F_pol_long + fourpi * toto * dot_product( P_long_k_loc, P_long_k_loc)    /chi_l(k_index)
        
        F_pol_trans = F_pol_trans +        toto * dot_product( P_trans_k_loc, P_trans_k_loc)  /chi_t(k_index)
        
        F_pol_tot = F_pol_tot &
                   - thermocond%kbT*3/(2*mu_SPCE**2*solvent(1)%n0) *deltaVkn*facsym*solvent(s)%rho0**2&
                        * dot_product( P_tot_k_loc , P_tot_k_loc )

        DO CONCURRENT ( o=1:angGrid%n_angles, p=1:molRotGrid%n_angles )
            Lweight = angGrid%weight(o) * molRotGrid%weight(p)
            molec_polar_k_loc = solvent(s)%molec_polar_k(:,i,j,k,o,p)
!========================================================================================================================
!Evaluate gradient
!========================================================================================================================
            IF ( k2(i,j,k)<=epsilon(1.0_dp) ) THEN
                dF_pol_long_k(i,j,k,o,p,s) = zeroC
                dF_pol_trans_k(i,j,k,o,p,s) = solvent(s)%rho0*0.5_dp*qfact/chi_t(k_index)*&
                    dot_product( molec_polar_k_loc , P_trans_k_loc) *Lweight*2.0_dp
            ELSE
                kvec = [kx(i),ky(j),kz(k)]
                k_tens_k_P = kvec * dot_product(    molec_polar_k_loc   ,kvec)

                dF_pol_long_k(i,j,k,o,p,s) = solvent(s)%rho0*0.5_dp*&
                    qfact/chi_l(k_index)*fourpi*Lweight*&
                    sum( k_tens_k_P * P_long_k_loc ) /k2(i,j,k)*2.0_dp

                dF_pol_trans_k(i,j,k,o,p,s) = 0.5_dp*solvent(s)%rho0 *qfact /chi_t(k_index) *Lweight *2.0_dp *&
                    sum(   P_trans_k_loc * (conjg(molec_polar_k_loc) - k_tens_k_P/k2(i,j,k))      )
            END IF
            
            dF_pol_tot_k(i,j,k,o,p,s)=-thermocond%kbT*3._dp*solvent(s)%rho0/(2._dp*mu_SPCE**2*solvent(1)%n0)*Lweight*2.0_dp*&
                dot_product(   molec_polar_k_loc     ,     P_tot_k_loc   )
                
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
                dF_pol_tot (:,:,:,o,p,s)=fftw3%out_backward*deltaVkn
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

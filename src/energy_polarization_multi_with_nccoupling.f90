SUBROUTINE energy_polarization_multi_with_nccoupling(F_Pol)

    USE precision_kinds, ONLY : i2b , dp
    USE system, ONLY : nfft1 , nfft2 , nfft3, Lx, Ly, Lz, kBT, rho_0,&
                   deltav, molec_polarx_k, molec_polary_k, molec_polarz_k,&
                   rho_0_multispec, nb_species, n_0, beta, deltax, deltay, deltaz, sigma_k, spaceGrid
    USE dcf, ONLY: Cnn, Cnc, Ccc, chi_t, c_s, nb_k, delta_k, delta_k_in_C, nb_k_in_c
    USE quadrature, ONLY: angGrid, molRotGrid,sym_order
    USE minimizer, ONLY: cg_vect , FF , dF
    USE constants, ONLY : twopi, i_complex, fourpi, eps0,qunit,Navo, qfact
    USE fft, ONLY: fftw3, kx, ky, kz, k2, norm_k
    USE input, ONLY: input_line, input_log, input_char, verbose

    IMPLICIT NONE
    
    REAL(dp), DIMENSION(nfft1, nfft2, nfft3,angGrid%n_angles,molRotGrid%n_angles,nb_species)::dF_pol_long , dF_pol_trans, dF_pol_tot
    COMPLEX(dp), DIMENSION(nfft1/2+1, nfft2, nfft3, angGrid%n_angles,molRotGrid%n_angles,nb_species) :: rho_k, dF_pol_long_k,&
     dF_pol_trans_k,dF_pol_tot_k
    REAL(dp), ALLOCATABLE, DIMENSION (:,:,:,:,:,:) :: rho
    REAL(dp), DIMENSION (nfft1, nfft2, nfft3, nb_species) ::Px,Py,Pz,P_long_x,P_long_y,P_long_z
    INTEGER(i2b) :: icg  !dummy counter for cg_vect
    REAL(dp) ::  deltaVk, F_pol_long, F_pol_trans , F_pol ,F_pol_tot  !Longitudinal , transverse and total Polarization free energy
    COMPLEX(dp), ALLOCATABLE, DIMENSION (:,:,:,:) :: rho_c_k_myway  !transverse part of polarization in Fourier space .
    INTEGER(i2b) :: i , j , k, o, p , n ,m,l, species, k_index , m1, m2, m3!dummy
    REAL(dp) :: mu_SPCE, facsym,Pxt,Pyt,Pzt
    COMPLEX(dp) :: k_tens_k_Px,k_tens_k_Py,k_tens_k_Pz     
    REAL(dp):: time1 , time0 ,time2 , time3 ,rhot! timestamps
    COMPLEX(dp) ::  F_pol_long_k , F_pol_trans_k , F_pol_tot_k, molec_polarx_k_local, molec_polary_k_local, molec_polarz_k_local
    REAL(dp), ALLOCATABLE, DIMENSION(:) :: weight_omx , weight_omy , weight_omz ! dummy
    INTEGER(i2b) :: ios ! iostat of the read statement: 
    INTEGER(i2b) :: kindex_in_C
    REAL(dp) ::  rhon, fact,psi,Vint,Fint
    REAL(dp), ALLOCATABLE, DIMENSION (:,:,:) :: rho_n  !one-particle number density
    COMPLEX(dp), ALLOCATABLE, DIMENSION(:,:,:) :: rho_n_k
    REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) :: Vpair
    COMPLEX(dp), ALLOCATABLE, DIMENSION(:,:,:) :: Vpair_k
    COMPLEX(dp), ALLOCATABLE, DIMENSION (:,:,:,:) :: &
        P_trans_x_k, P_trans_y_k, P_trans_z_k, P_long_x_k, P_long_y_k, P_long_z_k, pola_tot_x_k, pola_tot_y_k, pola_tot_z_k
    
    IF (nb_species/=1) THEN
        PRINT*, 'transv_and_longi_polarization_micro IS NOT WORKING FOR MULTISPECIES'
        STOP
    END IF
    CALL CPU_TIME(time0)
!===================================================================================================================================
!                	Initialization
!===================================================================================================================================
    deltaVk=twopi**3/PRODUCT(spaceGrid%length)
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
    rho_k=(0.0_dp,0.0_dp)
    IF (.NOT. ALLOCATED (pola_tot_x_k) ) ALLOCATE (pola_tot_x_k ( nfft1/2+1,nfft2,nfft3,nb_species), SOURCE=(0.0_dp,0.0_dp) )
    IF (.NOT. ALLOCATED (pola_tot_y_k) ) ALLOCATE (pola_tot_y_k ( nfft1/2+1,nfft2,nfft3,nb_species), SOURCE=(0.0_dp,0.0_dp) )
    IF (.NOT. ALLOCATED (pola_tot_z_k) ) ALLOCATE (pola_tot_z_k ( nfft1/2+1,nfft2,nfft3,nb_species), SOURCE=(0.0_dp,0.0_dp) )
    mu_SPCE=0.4238_dp*0.5773525_dp*2.0_dp !dipolar moment of SPCE water molecule in e.Angstromm
    IF (.NOT. ALLOCATED (rho_c_k_myway) ) ALLOCATE (rho_c_k_myway(nfft1/2+1, nfft2, nfft3, nb_species), SOURCE=(0.0_dp, 0.0_dp))
    

    !======================================================================================================
    !            ====================================================
    !            !    	Read Cnn(k), Cnc(k) and Ccc(k) 		!
    !            ====================================================
    
!    OPEN(10,file='input/direct_correlation_functions/water/Ccc.dat')
!    
!      nb_k_in_C=0
!      DO WHILE(.true.)
!        READ(10,*,iostat=ios)
!        IF (ios>0) THEN
!          WRITE(*,*)'Error in compute_ck_dipolar.f90'
!          WRITE(*,*)'something went wrong during the computation of the total number of lines in cs.in. stop'
!          STOP
!        ELSE IF (ios<0) THEN
!          ! end of file reached
!          EXIT
!        ELSE
!          nb_k_in_C=nb_k_in_C+1
!        END IF
!      END DO
!    
!    CLOSE(10)
!    
!    ALLOCATE(norm_k_in_c(nb_k_in_c))
!    ALLOCATE(Ccc(nb_k_in_c))
!    ALLOCATE(Cnc(nb_k_in_c))
!    ALLOCATE(Cnn(nb_k_in_c))
!    
!    OPEN(10,file='input/direct_correlation_functions/water/Ccc.dat')
!    OPEN(11,file='input/direct_correlation_functions/water/Cnc.dat')
!    OPEN(12,file='input/direct_correlation_functions/water/Cnn.dat')
!    
!    DO i=1, nb_k_in_c
!      READ(10,*) norm_k_in_c(i), Ccc(i)
!      READ(11,*) norm_k_in_c(i), Cnc(i)
!      READ(12,*) norm_k_in_c(i), Cnn(i)
!    END DO
!    delta_k_in_c = norm_k_in_c(2) - norm_k_in_c(1)
!    
!    CLOSE(10)
!    CLOSE(11)
!    CLOSE(12)
    !======================================================================================================
    !            ====================================================
    !            !    	Compute density in real and 		!
                 !       		Fourier space			!
                 !							!
    !            ====================================================
    !Compute rho_k
    ALLOCATE ( rho_n (nfft1,nfft2,nfft3), SOURCE=0.0_dp )
    ALLOCATE ( rho (nfft1, nfft2, nfft3, angGrid%n_angles, molRotGrid%n_angles, nb_species), SOURCE=0.0_dp )
    icg = 0
    DO species=1,nb_species
        DO i=1,nfft1
            DO j=1,nfft2
                DO k=1,nfft3
                    rhon=0.0_dp
                    DO o=1,angGrid%n_angles
                        DO p=1,molRotGrid%n_angles
                          icg = icg + 1
                          rhon=rhon+cg_vect(icg)**2*angGrid%weight(o)*molRotGrid%weight(p)!/(fourpi**2/(sym_order*2.0_dp))
                          rho(i,j,k,o,p,species) = cg_vect ( icg ) ** 2
                        END DO
                    END DO
                    rho_n(i,j,k)=rhon-REAL(2.0_dp*twopi**2/sym_order, dp)
                END DO
            END DO
        END DO
    END DO

    fftw3%in_forward=rho_n
    DEALLOCATE (rho_n)
    CALL dfftw_execute (fftw3%plan_forward)
    ALLOCATE ( rho_n_k (nfft1/2+1, nfft2, nfft3), SOURCE=(fftw3%out_forward *deltaV) )

    CALL cpu_time(time2)
    !======================================================================================================
    !            ====================================================
    !            !    		get Density 			!
    !            !			in Fourier Space		!
    !            ====================================================
    DO species=1, nb_species
        DO o=1, angGrid%n_angles
            DO p=1, molRotGrid%n_angles
                fftw3%in_forward = rho(:,:,:,o,p,species)
                CALL dfftw_execute (fftw3%plan_forward)
                rho_k(:,:,:,o,p,species)=fftw3%out_forward*deltaV
            END DO
        END DO
    END DO
    DEALLOCATE (rho)
       
    CALL cpu_time(time3)
    
    !======================================================================================================
    !            ====================================================
    !            !    		Compute 	                       		!
    !            !		part due to density/density coupling	    !
    !	         ! with cnn including a nc coupling part		    !
    !            ====================================================
    !======================================================================================================
    
    
    ALLOCATE( Vpair_k (nfft1/2+1,nfft2,nfft3), SOURCE=(0.0_dp,0.0_dp))
    DO n = 1 , nfft3
        DO m = 1 , nfft2
            DO l = 1 , nfft1 / 2 + 1
                kindex_in_c=int ( norm_k(l,m,n) / delta_k_in_C ) + 1
                ! Here it happens that k_index gets higher than the highest c_k index.
                ! In this case one imposes k_index = k_index_max
                IF ( kindex_in_c > nb_k_in_c ) kindex_in_c = nb_k_in_c     
                Vpair_k ( l , m , n ) = rho_n_k ( l , m , n ) * Cnn ( kindex_in_c )   !cnn is cs if there is no coupling
            END DO
        END DO
    END DO
    fftw3%in_backward = Vpair_k
    DEALLOCATE (Vpair_k)
    CALL dfftw_execute (fftw3%plan_backward)
    ALLOCATE ( Vpair (nfft1, nfft2, nfft3), SOURCE=(fftw3%out_backward *deltaVk/(twopi**3)))

    ! compute excess energy and its gradient
    Fint = 0.0_dp ! excess energy
    icg = 0 ! index of cg_vect
    DO species = 1 , nb_species
        fact = DeltaV * rho_0_multispec ( species ) !> facteur d'integration
        DO i = 1 , nfft1
            DO j = 1 , nfft2
                DO k = 1 , nfft3
                Vint   = -kBT * rho_0_multispec ( species ) * Vpair(i,j,k)
                    DO o = 1 , angGrid%n_angles
                        DO p=1, molRotGrid%n_angles
                            icg=icg + 1
                            psi=CG_vect ( icg )
                            Fint=Fint+angGrid%weight(o)*fact*0.5_dp*(psi**2-1.0_dp)*Vint*molRotGrid%weight(p)
                            dF (icg)=dF(icg)+2.0_dp*psi*angGrid%weight(o)*fact* Vint*molRotGrid%weight(p)!&
                        END DO       
                    END DO
                END DO
            END DO
        END DO
    END DO
    DEALLOCATE (Vpair)

    !======================================================================================================
    !            ====================================================
    !            !    		Compute 		                      	!
    !            !		Total Polarization		                	!
    !            ====================================================
    !======================================================================================================
    DO species=1, nb_species
        DO i=1, nfft1/2+1
            DO j=1, nfft2
                DO k=1, nfft3
                    DO o=1, angGrid%n_angles
                        DO p=1, molRotGrid%n_angles
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

    
    !!            ====================================================
    !!            !    	Compute 		                      	   	 !
    !!            !	Transverse and longitudinal Polarization    	 !
    !!            ====================================================
    ALLOCATE (P_trans_x_k (nfft1/2+1, nfft2, nfft3, nb_species), SOURCE=(0._dp,0._dp))
    ALLOCATE (P_trans_y_k (nfft1/2+1, nfft2, nfft3, nb_species), SOURCE=(0._dp,0._dp))
    ALLOCATE (P_trans_z_k (nfft1/2+1, nfft2, nfft3, nb_species), SOURCE=(0._dp,0._dp))
    ALLOCATE (P_long_x_k (nfft1/2+1, nfft2, nfft3, nb_species), SOURCE=(0._dp,0._dp))
    ALLOCATE (P_long_y_k (nfft1/2+1, nfft2, nfft3, nb_species), SOURCE=(0._dp,0._dp))
    ALLOCATE (P_long_z_k (nfft1/2+1, nfft2, nfft3, nb_species), SOURCE=(0._dp,0._dp))

    BLOCK
        COMPLEX(dp) :: riri, fifi, kiki, toto
        DO n=1, nb_species
            DO i=1, nfft1/2+1
                DO j=1, nfft2
                    DO k=1, nfft3
                        IF ( k2(i,j,k) /= 0._dp ) THEN
                            riri = pola_tot_x_k(i,j,k,n)
                            fifi = pola_tot_y_k(i,j,k,n)
                            kiki = pola_tot_z_k(i,j,k,n)
                            toto = ( riri*kx(i) + fifi*ky(j) + kiki*kz(k) ) / k2(i,j,k)
                            P_long_x_k(i,j,k,n)= toto * kx(i)
                            P_long_y_k(i,j,k,n)= toto * ky(j)
                            P_long_z_k(i,j,k,n)= toto * kz(k)
                        ELSE
                            P_long_x_k(i,j,k,n)=0.0_dp
                            P_long_y_k(i,j,k,n)=0.0_dp
                            P_long_z_k(i,j,k,n)=0.0_dp
                        END IF
                        P_trans_x_k(i,j,k,n) = pola_tot_x_k(i,j,k,n)- P_long_x_k(i,j,k,n)
                        P_trans_y_k(i,j,k,n) = pola_tot_y_k(i,j,k,n)- P_long_y_k(i,j,k,n)
                        P_trans_z_k(i,j,k,n) = pola_tot_z_k(i,j,k,n)- P_long_z_k(i,j,k,n)
                    END DO
                END DO
            END DO
        END DO
    END BLOCK
    
    DO species=1,nb_species
        fftw3%in_backward= P_long_x_k(:,:,:,species)
        CALL dfftw_execute(fftw3%plan_backward)
        P_long_x(:,:,:,species)=fftw3%out_backward*deltaVk/(twopi)**3
        fftw3%in_backward= P_long_y_k(:,:,:,species)
        CALL dfftw_execute(fftw3%plan_backward)
        P_long_y(:,:,:,species)=fftw3%out_backward*deltaVk/(twopi)**3
        fftw3%in_backward= P_long_z_k(:,:,:,species)
        CALL dfftw_execute(fftw3%plan_backward)
        P_long_z(:,:,:,species)=fftw3%out_backward*deltaVk/(twopi)**3
    END DO

    IF (verbose) THEN
        BLOCK
            REAL(dp), DIMENSION(nfft1, nfft2, nfft3, nb_species) :: pola_tot_x, pola_tot_y, pola_tot_z
            DO species=1,nb_species
                fftw3%in_backward= pola_tot_x_k(:,:,:,species)
                CALL dfftw_execute(fftw3%plan_backward)
                pola_tot_x(:,:,:,species)=fftw3%out_backward*deltaVk/(twopi)**3
                
                fftw3%in_backward= pola_tot_y_k(:,:,:,species)
                CALL dfftw_execute(fftw3%plan_backward)
                pola_tot_y(:,:,:,species)=fftw3%out_backward*deltaVk/(twopi)**3
                
                fftw3%in_backward= pola_tot_z_k(:,:,:,species)
                CALL dfftw_execute(fftw3%plan_backward)
                pola_tot_z(:,:,:,species)=fftw3%out_backward*deltaVk/(twopi)**3
            END DO
            OPEN (11, file='output/pola_tot_x')
            OPEN (12, file='output/pola_tot_y')
            OPEN (13, file='output/pola_tot_z')
            DO i=1,nfft1
                WRITE(11,*) i*deltax, pola_tot_x(i,nfft2/2+1,nfft3/2+1,1), Px(i,nfft2/2+1,nfft3/2+1,1), &
                            P_long_x(i,nfft2/2+1,nfft3/2+1,1) 
                WRITE(12,*) i*deltay, pola_tot_y(nfft1/2+1,i,nfft3/2+1,1), Py(nfft1/2+1,i,nfft3/2+1,1), &
                            P_long_y(nfft1/2+1,i,nfft3/2+1,1) 
                WRITE(13,*) i*deltaz, pola_tot_z(nfft1/2+1,nfft2/2+1,i,1), Pz(nfft1/2+1,nfft1/2+1,i,1), &
                            P_long_z(nfft1/2+1,nfft1/2+1,i,1)
            END DO
            CLOSE(11)
            CLOSE(12)
            CLOSE(13)
        END BLOCK
    END IF
    
    !            ====================================================
                !    	Compute rho_c_k	in k space		            !
                !						                         	!
    !            ====================================================
    
    DO species = 1, nb_species
        DO i = 1 , nfft1/2 + 1
            DO j = 1 , nfft2
                DO k = 1 , nfft3
                    DO o = 1 , angGrid%n_angles
                        DO p=1 , molRotGrid%n_angles
                            rho_c_k_myway(i,j,k,species)=rho_c_k_myway(i,j,k,species)+sigma_k(i,j,k,o,p,species)*&
                            rho_k (i , j , k , o , p , species)*rho_0*angGrid%weight(o)*molRotGrid%weight(p)!*deltavK         
                        END DO  !P
                    END DO   !O
                END DO  !K
            END DO  !J
        END DO  !I
    END DO
    
    !            ====================================================
    !            !    	Compute Free energy due to		            !
    !            !	Transverse and longitudinal Polarization	    !
    !            ====================================================

    DO species=1, nb_species 
        DO i=1, nfft1/2+1
            IF (i>1 .and. i<nfft1/2+1) THEN
                facsym=2.0_dp
            ELSE
                facsym=1.0_dp
            END IF
            DO j=1, nfft2
                DO k=1, nfft3
                     k_index = INT ( norm_k(i,j,k) / delta_k ) + 1
                     kindex_in_c=INT ( norm_k(i,j,k) / delta_k_in_C ) + 1
                     ! Here it happens that k_index gets higher than the highest c_k index.
                     ! In this case one imposes k_index = k_index_max
                     IF ( kindex_in_c > nb_k_in_c ) kindex_in_c = nb_k_in_c 
                     IF ( k_index > nb_k ) k_index = nb_k
                     !  if(k2(i,j,k)/=0.0_dp) then
                     F_pol_long_k=F_pol_long_k+deltaVk*fourpi*(P_long_x_k(i,j,k,species)*CONJG(p_long_x_k(i,j,k,species))+&
                     p_long_y_k(i,j,k,species)*CONJG(p_long_y_k(i,j,k,species))+p_long_z_k(i,j,k,species)&
                     *CONJG(p_long_z_k(i,j,k,species)))/Ccc(kindex_in_c)*0.5_dp*qfact*rho_0**2*facsym/(twopi**3)+&
                     (deltaVk*(rho_c_k_myway(i,j,k,species)*CONJG(rho_n_k(i,j,k))+CONJG(rho_c_k_myway(i,j,k,species))&
                     *rho_n_k(i,j,k))*rho_0*0.5_dp*facsym*Cnc(kindex_in_c))/(twopi)**3
    
                     F_pol_trans_k=F_pol_trans_k+deltaVk*(P_trans_x_k(i,j,k,species)*CONJG(p_trans_x_k(i,j,k,species))&
                     +p_trans_y_k(i,j,k,species)*CONJG(p_trans_y_k(i,j,k,species))+p_trans_z_k(i,j,k,species)*&
                     CONJG(p_trans_z_k(i,j,k,species)))/chi_t(k_index)*0.5_dp*qfact*rho_0**2*facsym/(twopi**3)
    
                     F_pol_tot_k=F_pol_tot_k-(kBT*3/(2*mu_SPCE**2*n_0)*deltaVk*facsym/(twopi**3)*rho_0**2*&
                     (pola_tot_x_k(i,j,k,species)*CONJG(pola_tot_x_k(i,j,k,species))+pola_tot_y_k(i,j,k,species)*&
                     CONJG(pola_tot_y_k(i,j,k,species))+pola_tot_z_k(i,j,k,species)*CONJG(pola_tot_z_k(i,j,k,species))))
                     DO o=1,angGrid%n_angles
                         DO p=1,molRotGrid%n_angles 
    !!========================================================================================================================
    !!Evaluate gradient
    !!========================================================================================================================
                                  molec_polarx_k_local = molec_polarx_k(i,j,k,o,p,species)
                                  molec_polary_k_local = molec_polary_k(i,j,k,o,p,species)
                                  molec_polarz_k_local = molec_polarz_k(i,j,k,o,p,species)
                                  k_tens_k_Px=CONJG(kx(i)*kx(i)*molec_polarx_k_local+&
                                  kx(i)*ky(j)*molec_polary_k_local+kx(i)*kz(k)*molec_polarz_k_local)
                                  k_tens_k_Py=CONJG(kx(i)*ky(j)*molec_polarx_k_local+&
                                  ky(j)*ky(j)*molec_polary_k_local+ky(j)*kz(k)*molec_polarz_k_local)
                                  k_tens_k_Pz=CONJG(kz(k)*kx(i)*molec_polarx_k_local+&
                                  kz(k)*ky(j)*molec_polary_k_local+kz(k)*kz(k)*molec_polarz_k_local)
       
                             IF (k2(i,j,k)==0.0_dp) THEN
                                    dF_pol_trans_k(i,j,k,o,p,species)=&
                                    rho_0*qfact/chi_t(k_index)*angGrid%weight(o)*molRotGrid%weight(p)*&
                                    (P_trans_x_k(i,j,k,species)*CONJG(molec_polarx_k_local)&
                                    +P_trans_y_k(i,j,k,species)*CONJG(molec_polary_k_local)&
                                    +P_trans_z_k(i,j,k,species)*CONJG(molec_polarz_k_local))
                                    dF_pol_long_k(i,j,k,o,p,species)=0.0_dp
    !             ===============================================================================
                             ELSE
                                    dF_pol_trans_k(i,j,k,o,p,species)=&
                                    rho_0*qfact/chi_t(k_index)*angGrid%weight(o)*molRotGrid%weight(p)*&
                                    (P_trans_x_k(i,j,k,species)*(CONJG(molec_polarx_k_local)-k_tens_k_Px/k2(i,j,k))&
                                    +P_trans_y_k(i,j,k,species)*(CONJG(molec_polary_k_local)-k_tens_k_Py/k2(i,j,k))&
                                    +P_trans_z_k(i,j,k,species)*(CONJG(molec_polarz_k_local)-k_tens_k_Pz/k2(i,j,k)))
                             END IF
    !             ===============================================================================
                                    dF_pol_tot_k(i,j,k,o,p,species)=-kBT*3.0_dp*rho_0/(mu_SPCE**2*n_0)*molRotGrid%weight(p)*&
                                    (pola_tot_x_k(i,j,k,species)*CONJG(molec_polarx_k_local)&
                                    +pola_tot_y_k(i,j,k,species)*CONJG(molec_polary_k_local)&
                                    +pola_tot_z_k(i,j,k,species)*CONJG(molec_polarz_k_local))*angGrid%weight(o)
            
                             IF(k2(i,j,k)/=0.0_dp) THEN
                                    dF_pol_long_k(i,j,k,o,p,species)=rho_0*qfact/Ccc(kindex_in_c)*fourpi*angGrid%weight(o)&
                                    *molRotGrid%weight(p)*(P_long_x_k(i,j,k,species)*k_tens_k_Px+P_long_y_k(i,j,k,species)&
                                    *k_tens_k_Py+P_long_z_k(i,j,k,species)*k_tens_k_Pz)/k2(i,j,k)+&
                                    rho_0*Cnc(kindex_in_c)*angGrid%weight(o)*molRotGrid%weight(p)*&
                                    (rho_n_k(i,j,k)*CONJG(sigma_k(i,j,k,o,p,species))+rho_c_k_myway(i,j,k,species)/rho_0)
                             END IF
                        END DO  !PSI
                    END DO   !OMEGA
                END DO    !K
            END DO    !J
        END DO    !I
    END DO     !SPECIES
    DEALLOCATE ( rho_n_k )
    DEALLOCATE ( P_trans_x_k, P_trans_y_k, P_trans_z_k )
    !!========================================================================================================================
    !!						Get gradient in real space
    !!========================================================================================================================
    DO species=1 , nb_species
        DO o=1,angGrid%n_angles
            DO p=1, molRotGrid%n_angles
                fftw3%in_backward= (dF_pol_trans_k (:,:,:,o,p,species)+dF_pol_long_k (:,:,:,o,p,species)+&
                dF_pol_tot_k (:,:,:,o,p,species))
                CALL dfftw_execute (fftw3%plan_backward)
                dF_pol_tot (:,:,:,o,p,species)=fftw3%out_backward*deltaVk/(twopi)**3
            END DO
        END DO
    END DO
    icg=0
    !!========================================================================================================================
    !!						Allocate it for minimizing
    !!========================================================================================================================
    DO species=1, nb_species
        DO i=1,nfft1
            DO j=1, nfft2
                DO k=1,nfft3
                    DO o=1,angGrid%n_angles
                        DO p=1, molRotGrid%n_angles
                            icg=icg+1
                            dF(icg)=dF(icg)+dF_pol_tot(i,j,k,o,p,species)*cg_vect(icg)*rho_0*2.0_dp*deltaV!* weight(o) * molRotGrid%weight(p) 
                        END DO
                    END DO
                END DO
            END DO
        END DO
    END DO
    !!========================================================================================================================
    !!					Check if Polarization Free energy is Real
    !!========================================================================================================================
!    IF (AIMAG(F_pol_tot_k+F_pol_long_k+F_pol_trans_k)<tiny(0.0_dp)) THEN 
!       F_pol_tot=REAL( F_pol_tot_k , dp)
!       F_pol_long=REAL( F_pol_long_k,  dp)
!       F_pol_trans=REAL( F_pol_trans_k, dp)
!    ELSE
!        PRINT*, 'Error in energy_polarization_myway Free energy is not Real'
!        PRINT*,AIMAG(F_pol_tot_k+F_pol_long_k+F_pol_trans_k)
!        PRINT*,F_pol_tot_k+F_pol_long_k+F_pol_trans_k
!        STOP
!    END IF

    !!!========================================================================================================================
    F_pol=F_pol_tot+F_pol_long+F_pol_trans+Fint
    FF=FF+F_pol
!    Print*, 'Fpol=', F_pol, 'Fint=', Fint
    ! stop timer
    CALL cpu_time ( time1 )

END SUBROUTINE energy_polarization_multi_with_nccoupling

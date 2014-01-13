SUBROUTINE energy_polarization_multi_with_nccoupling(F_pol)

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
    
    REAL(dp), DIMENSION(nfft1, nfft2, nfft3,angGrid%n_angles,molRotGrid%n_angles,nb_species)::dF_pol_tot!, dF_pol_long , dF_pol_trans
    COMPLEX(dp), DIMENSION(nfft1/2+1, nfft2, nfft3, angGrid%n_angles,molRotGrid%n_angles,nb_species) :: dF_pol_long_k,&
                                                                                            dF_pol_trans_k, dF_pol_tot_k
    REAL(dp), ALLOCATABLE, DIMENSION (:,:,:,:,:,:) :: rho
    REAL(dp) ::  deltaVk, F_pol_long, F_pol_trans , F_pol ,F_pol_tot  !Longitudinal , transverse and total Polarization free energy
    COMPLEX(dp), ALLOCATABLE, DIMENSION (:,:,:,:) :: rho_c_k_myway  !transverse part of polarization in Fourier space .
    INTEGER(i2b) :: icg,i,j,k,o,p,n,m,l,s,species, k_index , m1, m2, m3!dummy
    REAL(dp) :: mu_SPCE, facsym,Pxt,Pyt,Pzt
    INTEGER(i2b) :: ios,kindex_in_C
    REAL(dp) ::  rhon, fact,psi,Vint,Fint
    REAL(dp), ALLOCATABLE, DIMENSION(:) :: weight_omx , weight_omy , weight_omz ! dummy
    REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) :: rho_n, Vpair
    COMPLEX(dp), PARAMETER :: zeroC=(0._dp,0._dp)
    COMPLEX(dp) :: F_pol_long_k, F_pol_trans_k, F_pol_tot_k, molec_polarx_k_local, molec_polary_k_local, molec_polarz_k_local
    COMPLEX(dp), ALLOCATABLE, DIMENSION(:,:,:) :: rho_n_k
    COMPLEX(dp), ALLOCATABLE, DIMENSION(:,:,:) :: Vpair_k
    COMPLEX(dp), ALLOCATABLE, DIMENSION(:,:,:,:) :: P_tot_x_k, P_tot_y_k, P_tot_z_k, &
                                                    P_long_x_k, P_long_y_k, P_long_z_k, &
                                                    P_trans_x_k, P_trans_y_k, P_trans_z_k
    COMPLEX(dp), ALLOCATABLE, DIMENSION(:,:,:,:,:,:) :: rho_k
    
    IF (nb_species/=1) THEN
        PRINT*, 'transv_and_longi_polarization_micro IS NOT WORKING FOR MULTISPECIES'
        STOP
    END IF
!===================================================================================================================================
!                	Initialization
!===================================================================================================================================
    deltaVk=twopi**3/PRODUCT(spaceGrid%length)

    dF_pol_tot=0.0_dp
    dF_pol_long_k=0.0_dp
    dF_pol_trans_k=0.0_dp
    dF_pol_tot_k=0.0_dp
    mu_SPCE=0.4238_dp*0.5773525_dp*2.0_dp !dipolar moment of SPCE water molecule in e.Angstromm
    


    ALLOCATE ( rho_n (nfft1,nfft2,nfft3), SOURCE=0.0_dp )
    ALLOCATE ( rho (nfft1, nfft2, nfft3, angGrid%n_angles, molRotGrid%n_angles, nb_species), SOURCE=0.0_dp )
    icg = 0
    DO s=1,nb_species
        DO i=1,nfft1
            DO j=1,nfft2
                DO k=1,nfft3
                    rhon=0.0_dp
                    DO o=1,angGrid%n_angles
                        DO p=1,molRotGrid%n_angles
                          icg = icg + 1
                          rhon=rhon+cg_vect(icg)**2*angGrid%weight(o)*molRotGrid%weight(p)!/(fourpi**2/(sym_order*2.0_dp))
                          rho(i,j,k,o,p,s) = cg_vect ( icg ) ** 2
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

    !======================================================================================================
    !            ====================================================
    !            !    		get Density 			!
    !            !			in Fourier Space		!
    !            ====================================================
    ALLOCATE( rho_k(nfft1/2+1,nfft2,nfft3,angGrid%n_angles,molRotGrid%n_angles,nb_species), SOURCE=zeroC)
    DO s=1, nb_species
        DO p=1, molRotGrid%n_angles
            DO o=1, angGrid%n_angles
                fftw3%in_forward = rho(:,:,:,o,p,s)
                CALL dfftw_execute (fftw3%plan_forward)
                rho_k(:,:,:,o,p,s)=fftw3%out_forward*deltaV
            END DO
        END DO
    END DO
    DEALLOCATE (rho)
       
    !======================================================================================================
    !            ====================================================
    !            !    		Compute 	                       		!
    !            !		part due to density/density coupling	    !
    !	         ! with cnn including a nc coupling part		    !
    !            ====================================================
    !======================================================================================================
    
    
    ALLOCATE( Vpair_k (nfft1/2+1,nfft2,nfft3), SOURCE=zeroC)
    DO n = 1 , nfft3
        DO m = 1 , nfft2
            DO l = 1 , nfft1 / 2 + 1
                kindex_in_c= MIN( INT( norm_k(l,m,n) / delta_k_in_C ) + 1, nb_k_in_c)
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
    DO s = 1 , nb_species
        fact = DeltaV * rho_0_multispec ( s ) !> facteur d'integration
        DO i = 1 , nfft1
            DO j = 1 , nfft2
                DO k = 1 , nfft3
                Vint   = -kBT * rho_0_multispec ( s ) * Vpair(i,j,k)
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
    ALLOCATE (P_tot_x_k ( nfft1/2+1,nfft2,nfft3,nb_species), SOURCE=zeroC )
    ALLOCATE (P_tot_y_k ( nfft1/2+1,nfft2,nfft3,nb_species), SOURCE=zeroC )
    ALLOCATE (P_tot_z_k ( nfft1/2+1,nfft2,nfft3,nb_species), SOURCE=zeroC )
    BLOCK
        COMPLEX(dp) :: rho_k_tmp
        DO CONCURRENT( i=1:nfft1/2+1, j=1:nfft2, k=1:nfft3, s=1:nb_species, o=1:angGrid%n_angles, p=1:molRotGrid%n_angles )
            rho_k_tmp = rho_k(i,j,k,o,p,s)*angGrid%weight(o)*molRotGrid%weight(p)
            P_tot_x_k(i,j,k,s)=P_tot_x_k(i,j,k,s)+molec_polarx_k(i,j,k,o,p,s)*rho_k_tmp
            P_tot_y_k(i,j,k,s)=P_tot_y_k(i,j,k,s)+molec_polary_k(i,j,k,o,p,s)*rho_k_tmp
            P_tot_z_k(i,j,k,s)=P_tot_z_k(i,j,k,s)+molec_polarz_k(i,j,k,o,p,s)*rho_k_tmp
        END DO
    END BLOCK

    
    !!            ====================================================
    !!            !    	Compute 		                      	   	 !
    !!            !	Transverse and longitudinal Polarization    	 !
    !!            ====================================================
    ALLOCATE (P_long_x_k (nfft1/2+1, nfft2, nfft3, nb_species),  SOURCE=zeroC)
    ALLOCATE (P_long_y_k (nfft1/2+1, nfft2, nfft3, nb_species),  SOURCE=zeroC)
    ALLOCATE (P_long_z_k (nfft1/2+1, nfft2, nfft3, nb_species),  SOURCE=zeroC)
    BLOCK
        COMPLEX(dp) :: riri, fifi, kiki, toto
        DO s=1, nb_species
            DO k=1, nfft3
                DO j=1, nfft2
                    DO i=1, nfft1/2+1
                        riri = P_tot_x_k(i,j,k,s)
                        fifi = P_tot_y_k(i,j,k,s)
                        kiki = P_tot_z_k(i,j,k,s)
                        IF ( k2(i,j,k) /= 0._dp ) THEN
                            toto = ( riri*kx(i) + fifi*ky(j) + kiki*kz(k) ) / k2(i,j,k)
                            P_long_x_k(i,j,k,s)= toto * kx(i)
                            P_long_y_k(i,j,k,s)= toto * ky(j)
                            P_long_z_k(i,j,k,s)= toto * kz(k)
                        ELSE
                            P_long_x_k(i,j,k,s)=zeroC
                            P_long_y_k(i,j,k,s)=zeroC
                            P_long_z_k(i,j,k,s)=zeroC
                        END IF
                    END DO
                END DO
            END DO
        END DO
    END BLOCK
    ALLOCATE (P_trans_x_k (nfft1/2+1, nfft2, nfft3, nb_species), SOURCE=(P_tot_x_k-P_long_x_k))
    ALLOCATE (P_trans_y_k (nfft1/2+1, nfft2, nfft3, nb_species), SOURCE=(P_tot_y_k-P_long_y_k))
    ALLOCATE (P_trans_z_k (nfft1/2+1, nfft2, nfft3, nb_species), SOURCE=(P_tot_z_k-P_long_z_k))

    IF (verbose) THEN
        BLOCK
            REAL(dp), DIMENSION(nfft1, nfft2, nfft3, nb_species) :: P_tot_x, P_tot_y, P_tot_z, P_long_x, P_long_y, P_long_z

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

                fftw3%in_backward= P_tot_x_k(:,:,:,s)
                CALL dfftw_execute(fftw3%plan_backward)
                P_tot_x(:,:,:,s)=fftw3%out_backward*deltaVk/(twopi)**3
                
                fftw3%in_backward= P_tot_y_k(:,:,:,s)
                CALL dfftw_execute(fftw3%plan_backward)
                P_tot_y(:,:,:,s)=fftw3%out_backward*deltaVk/(twopi)**3
                
                fftw3%in_backward= P_tot_z_k(:,:,:,s)
                CALL dfftw_execute(fftw3%plan_backward)
                P_tot_z(:,:,:,s)=fftw3%out_backward*deltaVk/(twopi)**3
            END DO
            OPEN (11, file='output/P_tot_x')
            OPEN (12, file='output/P_tot_y')
            OPEN (13, file='output/P_tot_z')
            DO i=1,nfft1
                WRITE(11,*) i*deltax, P_tot_x(i,nfft2/2+1,nfft3/2+1,1), P_long_x(i,nfft2/2+1,nfft3/2+1,1) 
                WRITE(12,*) i*deltay, P_tot_y(nfft1/2+1,i,nfft3/2+1,1), P_long_y(nfft1/2+1,i,nfft3/2+1,1) 
                WRITE(13,*) i*deltaz, P_tot_z(nfft1/2+1,nfft2/2+1,i,1), P_long_z(nfft1/2+1,nfft1/2+1,i,1)
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
    IF (.NOT. ALLOCATED (rho_c_k_myway) ) ALLOCATE (rho_c_k_myway(nfft1/2+1, nfft2, nfft3, nb_species), SOURCE=zeroC)
    DO s =1, nb_species
        DO p =1, molRotGrid%n_angles
            DO o =1, angGrid%n_angles
                DO k =1, nfft3
                    DO j =1, nfft2
                        DO i =1, nfft1/2+1
                            rho_c_k_myway(i,j,k,s)=rho_c_k_myway(i,j,k,s)+sigma_k(i,j,k,o,p,s)*&
                            rho_k(i,j,k,o,p,s)*rho_0*angGrid%weight(o)*molRotGrid%weight(p)
                        END DO
                    END DO
                END DO
            END DO
        END DO
    END DO
    DEALLOCATE(rho_k)


    !            ====================================================
    !            !    	Compute Free energy due to		            !
    !            !	Transverse and longitudinal Polarization	    !
    !            ====================================================
    BLOCK
        COMPLEX(dp) :: k_tens_k_Px, k_tens_k_Py, k_tens_k_Pz
        
        F_pol_long_k = zeroC
        F_pol_trans_k = zeroC
        F_pol_tot_k = zeroC
        
        DO CONCURRENT (i=1:nfft1/2+1, j=1:nfft2, k=1:nfft3, s=1:nb_species)
            IF (i>1 .and. i<nfft1/2+1) THEN; facsym=2.0_dp; ELSE; facsym=1.0_dp; END IF
            k_index = MIN( INT( norm_k(i,j,k) / delta_k ) + 1, nb_k)
            kindex_in_c= MIN( INT( norm_k(i,j,k) / delta_k_in_C ) + 1, nb_k_in_c)
   
            F_pol_long_k =F_pol_long_k +deltaVk*fourpi*(&
                P_long_x_k(i,j,k,s)*CONJG(P_long_x_k(i,j,k,s))+&
                P_long_y_k(i,j,k,s)*CONJG(P_long_y_k(i,j,k,s))+&
                P_long_z_k(i,j,k,s)*CONJG(P_long_z_k(i,j,k,s)))&
                /Ccc(kindex_in_c)*0.5_dp*qfact*rho_0**2*facsym/(twopi**3)+&
                (deltaVk*(rho_c_k_myway(i,j,k,s)*CONJG(rho_n_k(i,j,k))+CONJG(rho_c_k_myway(i,j,k,s))&
                *rho_n_k(i,j,k))*rho_0*0.5_dp*facsym*Cnc(kindex_in_c))/(twopi)**3

            F_pol_trans_k =F_pol_trans_k +deltaVk*(&
                P_trans_x_k(i,j,k,s)*CONJG(P_trans_x_k(i,j,k,s))+&
                P_trans_y_k(i,j,k,s)*CONJG(P_trans_y_k(i,j,k,s))+&
                P_trans_z_k(i,j,k,s)*CONJG(P_trans_z_k(i,j,k,s)))&
                /chi_t(k_index)*0.5_dp*qfact*rho_0**2*facsym/(twopi**3)

            F_pol_tot_k =F_pol_tot_k-(kBT*3/(2*mu_SPCE**2*n_0)*deltaVk*facsym/(twopi**3)*rho_0**2*&
                (P_tot_x_k(i,j,k,s)*CONJG(P_tot_x_k(i,j,k,s))&
                +P_tot_y_k(i,j,k,s)*CONJG(P_tot_y_k(i,j,k,s))&
                +P_tot_z_k(i,j,k,s)*CONJG(P_tot_z_k(i,j,k,s))))

            DO CONCURRENT (o=1:angGrid%n_angles, p=1:molRotGrid%n_angles)
        !!========================================================================================================================
        !!Evaluate gradient
        !!========================================================================================================================
                molec_polarx_k_local = molec_polarx_k(i,j,k,o,p,s)
                molec_polary_k_local = molec_polary_k(i,j,k,o,p,s)
                molec_polarz_k_local = molec_polarz_k(i,j,k,o,p,s)
        
                IF (k2(i,j,k)==0.0_dp) THEN
                    dF_pol_trans_k(i,j,k,o,p,s)=&
                        rho_0*qfact/chi_t(k_index)*angGrid%weight(o)*molRotGrid%weight(p)*&
                        (P_trans_x_k(i,j,k,s)*CONJG(molec_polarx_k_local)&
                        +P_trans_y_k(i,j,k,s)*CONJG(molec_polary_k_local)&
                        +P_trans_z_k(i,j,k,s)*CONJG(molec_polarz_k_local))
                    
                    dF_pol_long_k(i,j,k,o,p,s)=0.0_dp
                ELSE
                    k_tens_k_Px = CONJG(kx(i)*kx(i)*molec_polarx_k_local+&
                                        kx(i)*ky(j)*molec_polary_k_local+&
                                        kx(i)*kz(k)*molec_polarz_k_local)/k2(i,j,k)
                    k_tens_k_Py = CONJG(ky(j)*kx(i)*molec_polarx_k_local+&
                                        ky(j)*ky(j)*molec_polary_k_local+&
                                        ky(j)*kz(k)*molec_polarz_k_local)/k2(i,j,k)
                    k_tens_k_Pz = CONJG(kz(k)*kx(i)*molec_polarx_k_local+&
                                        kz(k)*ky(j)*molec_polary_k_local+&
                                        kz(k)*kz(k)*molec_polarz_k_local)/k2(i,j,k)
                                        
                    dF_pol_trans_k(i,j,k,o,p,s)=&
                        rho_0*qfact/chi_t(k_index)*angGrid%weight(o)*molRotGrid%weight(p)*&
                        (P_trans_x_k(i,j,k,s)*(CONJG(molec_polarx_k_local)-k_tens_k_Px)&
                        +P_trans_y_k(i,j,k,s)*(CONJG(molec_polary_k_local)-k_tens_k_Py)&
                        +P_trans_z_k(i,j,k,s)*(CONJG(molec_polarz_k_local)-k_tens_k_Pz))
                    
                    dF_pol_long_k(i,j,k,o,p,s)=rho_0*qfact/Ccc(kindex_in_c)*fourpi*angGrid%weight(o)&
                        *molRotGrid%weight(p)*(&
                        P_long_x_k(i,j,k,s)*k_tens_k_Px+&
                        P_long_y_k(i,j,k,s)*k_tens_k_Py+&
                        P_long_z_k(i,j,k,s)*k_tens_k_Pz)+&
                        rho_0*Cnc(kindex_in_c)*angGrid%weight(o)*molRotGrid%weight(p)*&
                        (rho_n_k(i,j,k)*CONJG(sigma_k(i,j,k,o,p,s))+rho_c_k_myway(i,j,k,s)/rho_0)
                END IF
    
                dF_pol_tot_k(i,j,k,o,p,s)=-kBT*3.0_dp*rho_0/(mu_SPCE**2*n_0)*molRotGrid%weight(p)*&
                    (P_tot_x_k(i,j,k,s)*CONJG(molec_polarx_k_local)&
                    +P_tot_y_k(i,j,k,s)*CONJG(molec_polary_k_local)&
                    +P_tot_z_k(i,j,k,s)*CONJG(molec_polarz_k_local))*angGrid%weight(o)
            END DO
        END DO
    END BLOCK
    DEALLOCATE ( rho_n_k, P_trans_x_k, P_trans_y_k, P_trans_z_k, rho_c_k_myway )

    !!========================================================================================================================
    !!						Get gradient in real space
    !!========================================================================================================================
    DO s=1 , nb_species
        DO p=1, molRotGrid%n_angles
            DO o=1,angGrid%n_angles
                fftw3%in_backward = dF_pol_trans_k(:,:,:,o,p,s) + dF_pol_long_k(:,:,:,o,p,s) + dF_pol_tot_k(:,:,:,o,p,s)
                CALL dfftw_execute (fftw3%plan_backward)
                dF_pol_tot (:,:,:,o,p,s)= fftw3%out_backward*deltaVk/(twopi)**3
            END DO
        END DO
    END DO
    
    !!========================================================================================================================
    !!						Allocate it for minimizing
    !!========================================================================================================================
    icg=0
    DO s=1, nb_species
        DO i=1,nfft1
            DO j=1, nfft2
                DO k=1,nfft3
                    DO o=1,angGrid%n_angles
                        DO p=1, molRotGrid%n_angles
                            icg=icg+1
                            dF(icg)=dF(icg)+dF_pol_tot(i,j,k,o,p,s)*cg_vect(icg)*rho_0*2.0_dp*deltaV
                        END DO
                    END DO
                END DO
            END DO
        END DO
    END DO
    !!========================================================================================================================
    !!					Check if Polarization Free energy is Real
    !!========================================================================================================================
       F_pol_tot   =REAL( F_pol_tot_k, dp)
       F_pol_long  =REAL( F_pol_long_k, dp)
       F_pol_trans =REAL( F_pol_trans_k, dp)

    !!!========================================================================================================================
    F_pol = F_pol_tot+F_pol_long+F_pol_trans+Fint
    FF = FF +F_pol

END SUBROUTINE energy_polarization_multi_with_nccoupling

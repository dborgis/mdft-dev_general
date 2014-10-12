SUBROUTINE energy_polarization_multi_with_nccoupling(F_pol)

    USE precision_kinds, ONLY: i2b, dp
    USE system,          ONLY: spaceGrid, thermocond, nb_species, solvent
    USE dcf,             ONLY: Cnn, Cnc, Ccc, chi_t, nb_k, delta_k, delta_k_in_C, nb_k_in_c
    USE quadrature,      ONLY: angGrid, molRotGrid, molRotSymOrder
    USE minimizer,       ONLY: cg_vect, FF, dF
    USE constants,       ONLY: twopi, fourpi, qfact
    USE fft,             ONLY: fftw3, kx, ky, kz, k2, norm_k
    USE input,           ONLY: verbose

    IMPLICIT NONE
    
    INTEGER(i2b) :: icg,i,j,k,o,p,n,m,l,s,k_index,ios,kindex_in_C,d
    INTEGER(i2b), POINTER :: nfft1=>spaceGrid%n_nodes(1),&
                             nfft2=>spaceGrid%n_nodes(2),&
                             nfft3=>spaceGrid%n_nodes(3)
    REAL(dp) :: facsym,rhon,fact,psi,Vint,Fint,deltaVk,F_pol_long,F_pol_trans,F_pol,F_pol_tot
    REAL(dp), POINTER :: Lx => spaceGrid%length(1), Ly => spaceGrid%length(2), Lz => spaceGrid%length(3)
    REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) :: rho_n, Vpair
    REAL(dp), ALLOCATABLE, DIMENSION(:,:,:,:,:,:) :: rho,dF_pol
    COMPLEX(dp) :: F_pol_long_k, F_pol_trans_k, F_pol_tot_k, molec_polar_k_local(3)
    COMPLEX(dp), PARAMETER :: zeroC=(0._dp,0._dp)
    COMPLEX(dp), allocatable :: rho_n_k(:,:,:), Vpair_k(:,:,:)
    COMPLEX(dp), ALLOCATABLE, DIMENSION(:,:,:,:) :: rho_c_k_myway
    complex(dp), allocatable, dimension(:,:,:,:,:) :: P_tot_k, P_long_k, P_trans_k
    COMPLEX(dp), ALLOCATABLE, DIMENSION(:,:,:,:,:,:) :: rho_k, dF_pol_k
    real(dp) :: kvec(3)

    IF (size(solvent)/=1) STOP 'energy_polarization_...withnccoupling not working for several solvent species'

!===================================================================================================================================
!                	Initialization
!===================================================================================================================================
    deltaVk =twopi**3/PRODUCT(spaceGrid%length)
    !TODO HERE SHOULD BE A CHECK IF SOLVENT IS SPCE. IF NOT THEN INCONSISTENCY WITH MU_SPCE

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
                          rhon=rhon+cg_vect(icg)**2*angGrid%weight(o)*molRotGrid%weight(p)
                          rho(i,j,k,o,p,s) = cg_vect ( icg ) ** 2
                        END DO
                    END DO
                    rho_n(i,j,k)=rhon-REAL(2.0_dp*twopi**2/molRotSymOrder, dp)
                END DO
            END DO
        END DO
    END DO

    fftw3%in_forward=rho_n
    DEALLOCATE (rho_n)
    CALL dfftw_execute (fftw3%plan_forward)
    ALLOCATE ( rho_n_k (nfft1/2+1, nfft2, nfft3), SOURCE=(fftw3%out_forward *spaceGrid%dv) )

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
                rho_k(:,:,:,o,p,s)=fftw3%out_forward*spaceGrid%dv
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
        fact = spaceGrid%dv * solvent(s)%rho0 ! facteur d'integration
        DO i = 1 , nfft1
            DO j = 1 , nfft2
                DO k = 1 , nfft3
                Vint   = -thermocond%kbT * solvent(s)%rho0 * Vpair(i,j,k)
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
    allocate (P_tot_k (3,nfft1/2+1,nfft2,nfft3,nb_species), SOURCE=zeroC )
    BLOCK
        COMPLEX(dp) :: rho_k_tmp
        DO CONCURRENT( i=1:nfft1/2+1, j=1:nfft2, k=1:nfft3, s=1:nb_species, o=1:angGrid%n_angles, p=1:molRotGrid%n_angles )
            rho_k_tmp = rho_k(i,j,k,o,p,s) * angGrid%weight(o) * molRotGrid%weight(p)
            P_tot_k (:,i,j,k,s) = P_tot_k (:,i,j,k,s) + rho_k_tmp * solvent(s)%molec_polar_k(:,i,j,k,o,p)
        END DO
    END BLOCK

    
    !!            ====================================================
    !!            !    	Compute 		                      	   	 !
    !!            !	Transverse and longitudinal Polarization    	 !
    !!            ====================================================
    allocate (P_long_k (3,nfft1/2+1, nfft2, nfft3, nb_species),  SOURCE=zeroC)
    DO s=1, nb_species
        DO k=1, nfft3
            DO j=1, nfft2
                DO i=1, nfft1/2+1
                    kvec = [kx(i),ky(j),kz(k)]
                    IF ( abs(k2(i,j,k)) >= epsilon(1.0_dp) ) THEN
                        P_long_k(:,i,j,k,s) = kvec * sum( P_tot_k(:,i,j,k,s) * kvec ) / k2(i,j,k)
                    ELSE
                        P_long_k(:,i,j,k,s) = zeroC
                    END IF
                END DO
            END DO
        END DO
    END DO
    
    allocate (P_trans_k (3,nfft1/2+1, nfft2, nfft3, nb_species), SOURCE=(P_tot_k-P_long_k))

    IF (verbose) THEN
        BLOCK
            REAL(dp), DIMENSION(3,nfft1,nfft2,nfft3,nb_species) :: P_tot, P_long
            DO s=1,nb_species
                do d=1,3
                    fftw3%in_backward= P_long_k(d,:,:,:,s)
                    call dfftw_execute(fftw3%plan_backward)
                    P_long(d,:,:,:,s)=fftw3%out_backward*deltaVk/(twopi)**3
    
                    fftw3%in_backward= P_tot_k(d,:,:,:,s)
                    call dfftw_execute(fftw3%plan_backward)
                    P_tot(d,:,:,:,s)=fftw3%out_backward*deltaVk/(twopi)**3
                end do
            END DO
            OPEN (11, file='output/P_tot_x')
            DO i=1,nfft1
                WRITE(11,*) (i-1)*spaceGrid%dl(1), P_tot(1,i,nfft2/2+1,nfft3/2+1,1), P_long(1,i,nfft2/2+1,nfft3/2+1,1) 
            END DO
            close(11)
            OPEN (12, file='output/P_tot_y')
            do i=1,nfft2
                WRITE(12,*) (i-1)*spaceGrid%dl(2), P_tot(2,nfft1/2+1,i,nfft3/2+1,1), P_long(2,nfft1/2+1,i,nfft3/2+1,1) 
            end do
            close(12)
            OPEN (13, file='output/P_tot_z')
            do i=1,nfft3
                WRITE(13,*) (i-1)*spaceGrid%dl(3), P_tot(3,nfft1/2+1,nfft2/2+1,i,1), P_long(3,nfft1/2+1,nfft1/2+1,i,1)
            end do
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
                            rho_c_k_myway(i,j,k,s)=rho_c_k_myway(i,j,k,s)+solvent(s)%sigma_k(i,j,k,o,p)*&
                            rho_k(i,j,k,o,p,s)*solvent(s)%rho0*angGrid%weight(o)*molRotGrid%weight(p)
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
    ALLOCATE (dF_pol_k(nfft1/2+1,nfft2,nfft3,angGrid%n_angles,molRotGrid%n_angles,nb_species), SOURCE=zeroC)
    BLOCK
        COMPLEX(dp) :: k_tens_k_P(3), dF_pol_trans_k_tmp, dF_pol_long_k_tmp, dF_pol_tot_k_tmp
        
        F_pol_long_k = zeroC
        F_pol_trans_k = zeroC
        F_pol_tot_k = zeroC
        
        DO CONCURRENT (i=1:nfft1/2+1, j=1:nfft2, k=1:nfft3, s=1:nb_species)
            IF (i>1 .and. i<nfft1/2+1) THEN; facsym=2.0_dp; ELSE; facsym=1.0_dp; END IF
            k_index = MIN( INT( norm_k(i,j,k) / delta_k ) + 1, nb_k)
            kindex_in_c= MIN( INT( norm_k(i,j,k) / delta_k_in_C ) + 1, nb_k_in_c)
   
            F_pol_long_k =F_pol_long_k +deltaVk*fourpi*&
                dot_product( P_long_k(:,i,j,k,s) , P_long_k(:,i,j,k,s) ) &
                /Ccc(kindex_in_c)*0.5_dp*qfact*solvent(s)%rho0**2*facsym/(twopi**3)+&
                (deltaVk*(rho_c_k_myway(i,j,k,s)*CONJG(rho_n_k(i,j,k))+CONJG(rho_c_k_myway(i,j,k,s))&
                *rho_n_k(i,j,k))*solvent(s)%rho0*0.5_dp*facsym*Cnc(kindex_in_c))/(twopi)**3

            F_pol_trans_k =F_pol_trans_k +deltaVk* &
                dot_product( P_trans_k(:,i,j,k,s) , P_trans_k(:,i,j,k,s)  ) &
                /chi_t(k_index)*0.5_dp*qfact*solvent(s)%rho0**2*facsym/(twopi**3)

            F_pol_tot_k =F_pol_tot_k-(thermocond%kbT*3/(2*norm2(solvent(s)%dipole)**2&
                *solvent(1)%n0)*deltaVk*facsym/(twopi**3)*solvent(s)%rho0**2*&
                dot_product( P_tot_k(:,i,j,k,s) , P_tot_k(:,i,j,k,s))  )

            DO CONCURRENT (o=1:angGrid%n_angles, p=1:molRotGrid%n_angles)
        !!========================================================================================================================
        !!Evaluate gradient
        !!========================================================================================================================
        
                molec_polar_k_local = solvent(s)%molec_polar_k(:,i,j,k,o,p)
        
                IF (abs(k2(i,j,k))<=epsilon(1.0_dp)) THEN
                    dF_pol_trans_k_tmp =&
                        solvent(s)%rho0*qfact/chi_t(k_index)*angGrid%weight(o)*molRotGrid%weight(p)*&
                        dot_product( molec_polar_k_local, P_trans_k(:,i,j,k,s) )
                    dF_pol_long_k_tmp = cmplx(0.0_dp,0.0_dp)
                ELSE
                    kvec = [kx(i),ky(j),kz(k)]
                    k_tens_k_P(1) = conjg( kx(i) * sum(kvec * molec_polar_k_local) ) /k2(i,j,k)
                    k_tens_k_P(2) = conjg( ky(j) * sum(kvec * molec_polar_k_local) ) /k2(i,j,k)
                    k_tens_k_P(3) = conjg( kz(k) * sum(kvec * molec_polar_k_local) ) /k2(i,j,k)
                
                    dF_pol_trans_k_tmp =solvent(s)%rho0*qfact/chi_t(k_index)*angGrid%weight(o)*molRotGrid%weight(p)*&
                            sum( P_trans_k(:,i,j,k,s) * (conjg(molec_polar_k_local)-k_tens_k_P) )
                    
                    dF_pol_long_k_tmp =solvent(s)%rho0*qfact/Ccc(kindex_in_c)*fourpi*angGrid%weight(o) *molRotGrid%weight(p)*&
                        sum( P_long_k(:,i,j,k,s) * k_tens_k_P ) &
                        +solvent(s)%rho0*Cnc(kindex_in_c)*angGrid%weight(o)*molRotGrid%weight(p)*&
                        (rho_n_k(i,j,k)*CONJG(solvent(s)%sigma_k(i,j,k,o,p))&
                        +rho_c_k_myway(i,j,k,s)/solvent(s)%rho0)
                END IF
    
                dF_pol_tot_k_tmp =-thermocond%kbT*3.0_dp*solvent(s)%rho0/(norm2(solvent(s)%dipole)**2&
                    *solvent(1)%n0)*molRotGrid%weight(p)*&
                    dot_product(molec_polar_k_local , P_tot_k(:,i,j,k,s)) *angGrid%weight(o)

                dF_pol_k(i,j,k,o,p,s) = dF_pol_trans_k_tmp + dF_pol_long_k_tmp + dF_pol_tot_k_tmp
            END DO
        END DO
    END BLOCK
    deallocate ( rho_n_k, P_trans_k, P_long_k, P_tot_k, rho_c_k_myway )

    !!========================================================================================================================
    !!						Get gradient in real space
    !!========================================================================================================================
    ALLOCATE( dF_pol(nfft1, nfft2, nfft3,angGrid%n_angles,molRotGrid%n_angles,nb_species), SOURCE=0._dp)
    DO s=1 , nb_species
        DO p=1, molRotGrid%n_angles
            DO o=1,angGrid%n_angles
                fftw3%in_backward = dF_pol_k(:,:,:,o,p,s)
                CALL dfftw_execute (fftw3%plan_backward)
                dF_pol (:,:,:,o,p,s) =fftw3%out_backward*deltaVk/(twopi)**3
            END DO
        END DO
    END DO
    DEALLOCATE( dF_pol_k )
    
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
                            dF(icg)=dF(icg)+dF_pol(i,j,k,o,p,s)*cg_vect(icg)*solvent(s)%rho0*2.0_dp*spaceGrid%dv
                        END DO
                    END DO
                END DO
            END DO
        END DO
    END DO
    !!========================================================================================================================
    !!					Check if Polarization Free energy is Real
    !!========================================================================================================================
       if ( abs(imag(F_pol_tot_k))>=epsilon(1.0_dp) ) stop "F_pol_tot has non-nul imaginary part in energy...nccoupling.f90"
       if ( abs(imag(F_pol_long_k))>=epsilon(1.0_dp) ) stop "F_pol_long has non-nul imaginary part in energy...nccoupling.f90"
       if ( abs(imag(F_pol_trans_k))>=epsilon(1.0_dp) ) stop "F_pol_trans has non-nul imaginary part in energy...nccoupling.f90"
       F_pol_tot   =REAL( F_pol_tot_k, dp)
       F_pol_long  =REAL( F_pol_long_k, dp)
       F_pol_trans =REAL( F_pol_trans_k, dp)

    !!!========================================================================================================================
    F_pol = F_pol_tot+F_pol_long+F_pol_trans+Fint
    FF = FF +F_pol

END SUBROUTINE energy_polarization_multi_with_nccoupling

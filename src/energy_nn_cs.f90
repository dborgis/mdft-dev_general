! Compute total energy and gradients using direct correlation functions c_s
SUBROUTINE energy_nn_cs (Fint)

    USE precision_kinds ,ONLY: i2b,dp
    USE system          ,ONLY: kBT, rho_0_multispec, nb_species, spaceGrid
    USE dcf             ,ONLY: c_s, nb_k, delta_k
    USE quadrature      ,ONLY: molRotSymOrder, angGrid, molRotGrid
    USE minimizer       ,ONLY: cg_vect, FF, dF
    USE constants       ,ONLY: fourpi, pi, twopi, zeroC
    USE fft             ,ONLY: fftw3, norm_k

    IMPLICIT NONE

    INTEGER(i2b)             :: i,j,k,l,m,n,o,p,icg,species,nfft1,nfft2,nfft3,k_index,s
    REAL(dp), INTENT(OUT)    :: Fint ! Internal part of the free energy
    REAL(dp)                 :: Vint, fact, psi, time1, time0
    REAL(dp)   , ALLOCATABLE :: delta_rho(:,:,:), Vpair(:,:,:)
    COMPLEX(dp), ALLOCATABLE :: delta_rho_k(:,:,:), Vpair_k(:,:,:)

    nfft1 =spaceGrid%n_nodes(1)
    nfft2 =spaceGrid%n_nodes(2)
    nfft3 =spaceGrid%n_nodes(3)

    CALL CPU_TIME ( time0 )

    ALLOCATE ( delta_rho(nfft1,nfft2,nfft3) ,SOURCE=0._dp )

    ! Put density of last minimization step in delta_rho
    icg=0
    DO i=1,nfft1
        DO j=1,nfft2
            DO k=1,nfft3
                psi = 0.0_dp
                DO o = 1, angGrid%n_angles
                    DO p=1 , molRotGrid%n_angles
                        icg=icg+1
                        psi = psi + angGrid%weight(o) * cg_vect(icg)**2*molRotGrid%weight(p)
                    END DO
                END DO
                delta_rho(i,j,k) = psi
            END DO
        END DO
    END DO
    delta_rho = delta_rho -REAL(2.0_dp*twopi**2/molRotSymOrder, dp)
    
    ! Next FFT sequences can be done on multiple threads
    ! Compute rho in k-space
    fftw3%in_forward = delta_rho
    DEALLOCATE (delta_rho)
    CALL dfftw_execute ( fftw3%plan_forward )
    ALLOCATE ( delta_rho_k (nfft1/2+1,nfft2,nfft3) ,SOURCE=fftw3%out_forward )
    
    ! Compute polarisation in k-space
    ALLOCATE ( Vpair_k (nfft1/2+1,nfft2,nfft3) ,SOURCE=zeroC)
    ! V(k)=cs(k)*rho(k)
    DO CONCURRENT (l=1:nfft1/2+1, m=1:nfft2, n=1:nfft3)
        k_index =MIN(  INT(norm_k(l,m,n)/delta_k)+1  ,  nb_k  )
        Vpair_k(l,m,n) = delta_rho_k(l,m,n) * c_s(k_index)
    END DO
    DEALLOCATE ( delta_rho_k )

    fftw3%in_backward = Vpair_k
    DEALLOCATE ( Vpair_k )
    CALL dfftw_execute (fftw3%plan_backward)
    ALLOCATE ( Vpair (nfft1,nfft2,nfft3) ,SOURCE=fftw3%out_backward/REAL(nfft1*nfft2*nfft3, dp))

    ! compute excess energy and its gradient
    Fint = 0.0_dp ! excess energy
    icg = 0 ! index of cg_vect
    DO s =1,nb_species
        fact = -kBT * rho_0_multispec(s)**2 * spaceGrid%dV
        DO i =1,nfft1
            DO j =1,nfft2
                DO k =1,nfft3
                    Vint = Vpair(i,j,k) * fact
                    DO o =1,angGrid%n_angles
                        DO p =1,molRotGrid%n_angles
                            icg = icg + 1
                            psi = CG_vect(icg)
                            Fint = Fint + 0.5_dp * ( psi**2 - 1.0_dp) * angGrid%weight(o) * molRotGrid%weight(p) * Vint
                            dF(icg) = dF(icg) + 2.0_dp * psi          * angGrid%weight(o) * molRotGrid%weight(p) * Vint
                        END DO
                    END DO
                END DO
            END DO
        END DO
    END DO
    DEALLOCATE (Vpair)

    FF = FF + Fint
    CALL CPU_TIME (time1)
 
END SUBROUTINE energy_nn_cs

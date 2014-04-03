! This subroutine computes the excess free energy in the HNC approximation of the homogeneous reference fluid approximation of MDFT.
! It uses the rotational invariant c_s(r), which is an input of the code.

SUBROUTINE energy_nn_cs (Fint)

    USE precision_kinds ,ONLY: i2b, dp
    USE system          ,ONLY: kBT, rho_0_multispec, spaceGrid
    USE dcf             ,ONLY: c_s, nb_k, delta_k
    USE quadrature      ,ONLY: molRotSymOrder, angGrid, molRotGrid
    USE minimizer       ,ONLY: cg_vect, FF, dF
    USE constants       ,ONLY: twopi, zeroC
    USE fft             ,ONLY: fftw3, norm_k

    IMPLICIT NONE

    REAL(dp), INTENT(OUT)    :: Fint
    INTEGER(i2b)             :: i,j,k,l,m,n,o,p,icg,species,nfft1,nfft2,nfft3,k_index,s,nspec
    REAL(dp)                 :: Vint, fact, psi, time1, time0
    REAL(dp)   , ALLOCATABLE :: delta_rho(:,:,:), gamma(:,:,:)
    COMPLEX(dp), ALLOCATABLE :: delta_rho_k(:,:,:), gamma_k(:,:,:)
    LOGICAL    , SAVE        :: c_s_isOK =.FALSE.

    nfft1 =spaceGrid%n_nodes(1)
    nfft2 =spaceGrid%n_nodes(2)
    nfft3 =spaceGrid%n_nodes(3)

    CALL CPU_TIME ( time0 )

    ALLOCATE ( delta_rho(nfft1,nfft2,nfft3) ,SOURCE=0._dp )

    CALL fill_delta_rho_from_cgvect
    
    fftw3%in_forward = delta_rho
    DEALLOCATE (delta_rho)
    CALL dfftw_execute ( fftw3%plan_forward )
    ALLOCATE ( delta_rho_k (nfft1/2+1,nfft2,nfft3) ,SOURCE=fftw3%out_forward ) ! deltarho in k-space
    
    CALL is_c_s_OK
    
    ! Polarisation in k-space
    ALLOCATE ( gamma_k (nfft1/2+1,nfft2,nfft3) ,SOURCE=zeroC)
    CALL fill_gamma_in_kspace_ie_deltarho_times_c_s
    DEALLOCATE ( delta_rho_k )

    fftw3%in_backward = gamma_k
    DEALLOCATE ( gamma_k )
    CALL dfftw_execute (fftw3%plan_backward)
    ALLOCATE ( gamma (nfft1,nfft2,nfft3) ,SOURCE=fftw3%out_backward/REAL(nfft1*nfft2*nfft3, dp))

    ! excess free energy and its gradient
    Fint = 0.0_dp ! excess free energy
    icg = 0 ! index of cg_vect
    nspec = SIZE( rho_0_multispec ) ! number of implicit solvant species
    DO s =1,nspec
        fact = -kBT * rho_0_multispec(s)**2 * spaceGrid%dV
        DO i =1,nfft1
            DO j =1,nfft2
                DO k =1,nfft3
                    Vint = gamma(i,j,k) * fact
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
    DEALLOCATE (gamma)

    FF = FF + Fint
    CALL CPU_TIME (time1)

 
    CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        SUBROUTINE is_c_s_OK
            IF (.NOT. c_s_isOK) THEN
                IF (.NOT. ALLOCATED(c_s) .OR. ALL(c_s==0._dp)) STOP "I find c_s(k) to be unconsistent in energy_nn_cs.f90"
                c_s_isOK = .TRUE.
            END IF
        END SUBROUTINE is_c_s_OK
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
        SUBROUTINE fill_delta_rho_from_cgvect
            IMPLICIT NONE
            REAL(dp)        :: psi
            INTEGER(i2b)    :: icg,i,j,k,o,p
            icg=0
            DO i=1,nfft1
                DO j=1,nfft2
                    DO k=1,nfft3
                        psi = 0.0_dp
                        DO o=1,angGrid%n_angles
                            DO p=1,molRotGrid%n_angles
                                icg = icg+1
                                psi = psi + angGrid%weight(o) * cg_vect(icg)**2*molRotGrid%weight(p)
                            END DO
                        END DO
                        delta_rho(i,j,k) = psi
                    END DO
                END DO
            END DO
            delta_rho = delta_rho -REAL(2.0_dp*twopi**2/molRotSymOrder, dp)
        END SUBROUTINE fill_delta_rho_from_cgvect
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
        SUBROUTINE fill_gamma_in_kspace_ie_deltarho_times_c_s
            INTEGER(i2b) :: l,m,n
            DO CONCURRENT (l=1:nfft1/2+1, m=1:nfft2, n=1:nfft3)
                k_index =MIN(  INT(norm_k(l,m,n)/delta_k)+1  ,  nb_k  )
                gamma_k(l,m,n) = delta_rho_k(l,m,n) * c_s(k_index)     ! V(k)=cs(k)*rho(k)
            END DO
        END SUBROUTINE fill_gamma_in_kspace_ie_deltarho_times_c_s

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
END SUBROUTINE energy_nn_cs

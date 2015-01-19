! Compute total energy and gradients using direct correlation functions c_s_hs of a hard sphere fluid
SUBROUTINE energy_cs_hard_sphere (Fint)

    USE precision_kinds ,ONLY: i2b, dp
    USE system          ,ONLY: thermocond, nb_species, spaceGrid, solvent
    USE quadrature      ,ONLY: molRotSymOrder, angGrid, molRotGrid
    USE minimizer       ,ONLY: cg_vect, FF, dF
    USE constants       ,ONLY: fourpi, twopi, zero, zeroC
    USE fft             ,ONLY: fftw3, norm_k
    USE input           ,ONLY: verbose
    use dcf             ,only: c_s_hs
    use mathematica     ,only: splint

    IMPLICIT NONE
    INTEGER(i2b) :: i, j, k, l, m, n, o, p, s, icg, species, nfft1, nfft2, nfft3
    REAL(dp), INTENT(OUT) :: Fint ! Internal part of the free energy
    REAL(dp) :: Vint, fact, psi, lx, ly, lz, time1, time0
    REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) :: delta_rho, gamma
    COMPLEX(dp), ALLOCATABLE, DIMENSION(:,:,:) :: delta_rho_k


    CALL CPU_TIME ( time0 )

    nfft1 =spaceGrid%n_nodes(1)
    nfft2 =spaceGrid%n_nodes(2)
    nfft3 =spaceGrid%n_nodes(3)
    lx =spaceGrid%length(1)
    ly =spaceGrid%length(2)
    lz =spaceGrid%length(3)

    CALL build_delta_rho_from_last_minimizer_step
    CALL build_delta_rho_in_Fourier_space
    CALL build_gamma_is_delta_rho_k_dot_cs_k


    ! compute excess energy and its gradient
    Fint = 0.0_dp
    icg = 0
    DO s = 1 , nb_species
        fact = PRODUCT(spaceGrid%dl) * solvent(s)%rho0
        DO i = 1 , nfft1
            DO j = 1 , nfft2
                DO k = 1 , nfft3
                    Vint = -thermocond%kbt * solvent(s)%rho0 * gamma(i,j,k)
                    DO o = 1 , angGrid%n_angles
                        DO p=1, molRotGrid%n_angles
                            icg = icg + 1
                            psi = CG_vect ( icg )
                            Fint = Fint   + angGrid%weight(o)*molRotGrid%weight(p) * fact * 0.5_dp * ( psi ** 2 - 1.0_dp) * Vint
                  !         dF (icg) = dF ( icg ) + 2.0_dp * psi * angGrid%weight(o) * fact * Vint ! in case of bridge calculation, one deduces the pair contribution of hard spheres + => - and FF=FF-Fint
                            dF (icg) = dF ( icg ) - 2.0_dp * psi * angGrid%weight(o) *molRotGrid%weight(p)* fact * Vint
                        END DO
                    END DO
                END DO
            END DO
        END DO
    END DO
    DEALLOCATE(gamma)

    FF = FF - Fint

    CALL CPU_TIME(time1)
    IF (verbose) WRITE(*,*) 'Fexc c_hs   = ' , Fint , 'computed in (sec)' , time1 - time0


    CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        SUBROUTINE build_delta_rho_from_last_minimizer_step
            ALLOCATE( Delta_rho (nfft1,nfft2,nfft3) ,SOURCE=zero)
            icg=0
            DO i=1,nfft1
                DO j=1,nfft2
                    DO k=1,nfft3
                        DO o = 1, angGrid%n_angles
                            DO p=1, molRotGrid%n_angles
                                icg=icg+1
                                Delta_rho(i,j,k) = Delta_rho(i,j,k) + cg_vect(icg)**2 *molRotGrid%weight(p)*angGrid%weight(o)
                            END DO
                        END DO
                    END DO
                END DO
            END DO
            delta_rho = delta_rho-(twopi*fourpi)/REAL(molRotSymOrder,dp)
        END SUBROUTINE build_delta_rho_from_last_minimizer_step

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        SUBROUTINE build_delta_rho_in_Fourier_space
            ! Compute rho in k-space
            fftw3%in_forward = delta_rho
            DEALLOCATE (delta_rho)
            CALL dfftw_execute ( fftw3%plan_forward )
            ALLOCATE ( delta_rho_k (nfft1/2+1,nfft2,nfft3) ,SOURCE=fftw3%out_forward)
        END SUBROUTINE build_delta_rho_in_Fourier_space

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        SUBROUTINE build_gamma_is_delta_rho_k_dot_cs_k
            USE dcf ,ONLY: nb_k, delta_k
            use fft, only: norm_k
            implicit none
            real(dp) :: c_s_hs_loc
            ! gamma(k)=cs(k)*rho(k)
            ! gamma is named fftw3%in_backward not to have useless big table
            DO CONCURRENT (l=1:nfft1/2+1, m=1:nfft2, n=1:nfft3)
              call splint( xa=c_s_hs%x, ya=c_s_hs%y, y2a=c_s_hs%y2, n=size(c_s_hs%y), x=norm_k(l,m,n), y=c_s_hs_loc)
              fftw3%in_backward(l,m,n) = delta_rho_k(l,m,n) * c_s_hs_loc
            END DO
            DEALLOCATE ( delta_rho_k )
            CALL dfftw_execute (fftw3%plan_backward)
            ALLOCATE ( gamma(nfft1,nfft2,nfft3) ,SOURCE=fftw3%out_backward/REAL(nfft1*nfft2*nfft3,dp))
        END SUBROUTINE build_gamma_is_delta_rho_k_dot_cs_k

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE energy_cs_hard_sphere

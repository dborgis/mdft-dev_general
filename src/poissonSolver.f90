! This SUBROUTINE uses the Gauss theorem associated to the expression of the derivation in Fourier transforms to get V(r).
! Field E(r)=-grad(V(r))
! Local expression of Gauss theorem : div(E(r))=soluteChargeDensity(r)/eps0
! => laplacien(V(r))=-soluteChargeDensity(r)/eps0
! => V(k)=soluteChargeDensity(k)/(esp0*k^2)
! FFT(V(k)) = V(r)

SUBROUTINE poissonSolver (soluteChargeDensity)

    USE precision_kinds,    ONLY: dp, i2b
    USE system,             ONLY: nfft1 , nfft2 , nfft3, spaceGrid
    USE fft,                ONLY: fftw3 , norm_k , k2
    USE constants,          ONLY: fourpi , twopi
    USE external_potential, ONLY: V_c
    USE input,              ONLY: verbose
    ! V_c = electrostatic potential from charge density and poisson equation

    IMPLICIT NONE
    REAL(dp), DIMENSION (spaceGrid%n_nodes(1),spaceGrid%n_nodes(2),spaceGrid%n_nodes(3)), INTENT(IN) :: soluteChargeDensity
    COMPLEX(dp), DIMENSION (nfft1/2+1, nfft2, nfft3 ) :: soluteChargeDensity_k, V_c_k
    INTEGER (i2b) :: i,j,k

    IF ( MAXVAL(ABS(soluteChargeDensity)) < TINY(1.0_dp)) THEN
        ALLOCATE ( V_c ( nfft1 , nfft2 , nfft3 ), SOURCE=0._dp ) !~ v_c = 0.0_dp
    END IF

    ! FFT of soluteChargeDensity
    fftw3%in_forward = soluteChargeDensity
    CALL dfftw_execute ( fftw3%plan_forward )
    soluteChargeDensity_k = fftw3%out_forward ! It is verified that at this point, FFT-1(soluteChargeDensity_k)/ (nfft1*nfft2*nfft3) = soluteChargeDensity
    ! FFT(Laplacian(V(r))) = FFT( - 4Pi charge density(r) ) in elecUnits = (ik)^2 V(k) = -4pi rho(k)
    ! V(k) = 4Pi rho(k) / k^2

    DO CONCURRENT ( i=1:nfft1/2+1, j=1:nfft2, k=1:nfft3 )
        IF ( k2(i,j,k) /= 0._dp ) THEN
            V_c_k(i,j,k) = soluteChargeDensity_k(i,j,k) * fourpi/k2(i,j,k) ! in electrostatic units : V=-4pi rho
        ELSE
            V_c_k(i,j,k) = 0._dp
        END IF
    END DO

    ! get real space potential V(r)
    IF ( .NOT. ALLOCATED ( V_c ) ) ALLOCATE ( V_c ( nfft1 , nfft2 , nfft3 ) )
    fftw3%in_backward = V_c_k
    CALL dfftw_execute ( fftw3%plan_backward )
    V_c = fftw3%out_backward / REAL ( nfft1 * nfft2 * nfft3 , dp )

    IF (verbose) THEN
        OPEN(11,FILE='output/V_cmax.dat')
            DO i=1, nfft1
                WRITE(11,*), i , V_c(i, nfft2/2, nfft3/2 )
            END DO
        CLOSE(11)
    END IF

END SUBROUTINE poissonSolver

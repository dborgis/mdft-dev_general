! This SUBROUTINE uses the Gauss theorem associated to the expression of the derivation in Fourier transforms to get V(r).
! Field E(r)=-grad(V(r))
! Local expression of Gauss theorem : div(E(r))=soluteChargeDensity(r)/eps0
! => laplacien(V(r))=-soluteChargeDensity(r)/eps0
! => V(k)=soluteChargeDensity(k)/(esp0*k^2)
! FFT(V(k)) = V(r)

SUBROUTINE poissonSolver (gridnode, soluteChargeDensity, Vpoisson)

    USE precision_kinds,    ONLY: dp, i2b
    USE system,             ONLY: nfft1 , nfft2 , nfft3, spaceGrid
    USE fft,                ONLY: fftw3 , norm_k , k2
    USE constants,          ONLY: fourpi , twopi, zeroC
    USE input,              ONLY: verbose
    ! Vpoisson = electrostatic potential from charge density and poisson equation

    IMPLICIT NONE
    
    INTEGER(i2b), INTENT(IN) :: gridnode(3)
    REAL(dp), DIMENSION(gridnode(1),gridnode(2),gridnode(3)), INTENT(IN) :: soluteChargeDensity
    REAL(dp), DIMENSION(gridnode(1),gridnode(2),gridnode(3)), INTENT(OUT) :: Vpoisson
    COMPLEX(dp), ALLOCATABLE, DIMENSION(:,:,:) :: soluteChargeDensity_k,Vpoisson_k
    INTEGER (i2b) :: i,j,k


    IF ( ALL(soluteChargeDensity==0._dp) ) THEN
        PRINT*,'The solute wears no charge.'
        Vpoisson =0._dp
        RETURN
   
    ELSE
        ALLOCATE( soluteChargeDensity_k(nfft1/2+1,nfft2,nfft3), SOURCE=zeroC )
        ALLOCATE( Vpoisson_k(nfft1/2+1,nfft2,nfft3), SOURCE=zeroC )

        ! FFT of soluteChargeDensity
        fftw3%in_forward = soluteChargeDensity
        CALL dfftw_execute ( fftw3%plan_forward )
        soluteChargeDensity_k = fftw3%out_forward ! It is verified that at this point, FFT-1(soluteChargeDensity_k)/ (nfft1*nfft2*nfft3) = soluteChargeDensity
        ! FFT(Laplacian(V(r))) = FFT( - 4Pi charge density(r) ) in elecUnits = (ik)^2 V(k) = -4pi rho(k)
        ! V(k) = 4Pi rho(k) / k^2

        DO CONCURRENT ( i=1:nfft1/2+1, j=1:nfft2, k=1:nfft3 )
            IF ( k2(i,j,k) /= 0._dp ) THEN
                Vpoisson_k(i,j,k) = soluteChargeDensity_k(i,j,k) * fourpi/k2(i,j,k) ! in electrostatic units : V=-4pi rho
            ELSE
                Vpoisson_k(i,j,k) = 0._dp
            END IF
        END DO
    
        ! get real space potential V(r)
        fftw3%in_backward = Vpoisson_k
        CALL dfftw_execute ( fftw3%plan_backward )
        Vpoisson = fftw3%out_backward / REAL(nfft1*nfft2*nfft3,dp)
    
    END IF

END SUBROUTINE poissonSolver

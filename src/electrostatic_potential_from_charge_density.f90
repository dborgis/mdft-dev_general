! This subroutine uses the Gauss theorem associated to the expression of the derivation in Fourier transforms to get V(r).
! Field E(r)=-grad(V(r))
! Local expression of Gauss theorem : div(E(r))=rho_c(r)/eps0
! => laplacien(V(r))=-rho_c(r)/eps0
! => V(k)=rho_c(k)/(esp0*k^2)
! FFT(V(k)) = V(r)
SUBROUTINE electrostatic_potential_from_charge_density

    USE precision_kinds , ONLY : dp , i2b
    USE system , ONLY : nfft1 , nfft2 , nfft3 , rho_c
    ! nfft = number of FFT grid nodes in each direction
    ! rho_c (nfft1,nfft2,nfft3) = discrete charge density per unit volume
    USE fft , ONLY : in_forward , in_backward , out_forward , out_backward , plan_forward , plan_backward , norm_k , k2
    USE constants , ONLY : fourpi , twopi
    USE external_potential , ONLY : V_c
    ! V_c = electrostatic potential from charge density and poisson equation
    
    IMPLICIT NONE
    COMPLEX(dp), DIMENSION ( nfft1 / 2 + 1 , nfft2 , nfft3 ) :: rho_c_k
    COMPLEX(dp), DIMENSION ( nfft1 / 2 + 1 , nfft2 , nfft3 ) :: V_c_k
    INTEGER (i2b) :: i,j,k

    ! check if rho_c exists and is allocated
    IF ( .NOT. ALLOCATED ( rho_c ) ) THEN
        PRINT*, 'rho_c is not allocated in electrostatic_potential_from_charge_density.f90'
        STOP
    END IF

    IF ( MAXVAL(ABS(rho_c)) < TINY(1.0_dp)) THEN
        ALLOCATE ( V_c ( nfft1 , nfft2 , nfft3 ), SOURCE=0._dp ) !~ v_c = 0.0_dp
    END IF

    ! FFT of rho_c
    in_forward = rho_c
    CALL dfftw_execute ( plan_forward )
    rho_c_k = out_forward ! It is verified that at this point, FFT-1(rho_c_k)/ (nfft1*nfft2*nfft3) = rho_c
    ! FFT(Laplacian(V(r))) = FFT( - 4Pi charge density(r) ) in elecUnits = (ik)^2 V(k) = -4pi rho(k)
    ! V(k) = 4Pi rho(k) / k^2

    DO CONCURRENT ( i=1:nfft1/2+1, j=1:nfft2, k=1:nfft3 )
        IF ( K2(i,j,k) /= 0._dp ) THEN
            V_c_k(i,j,k) = rho_c_k(i,j,k) * fourpi/k2(i,j,k) ! in electrostatic units : V=-4pi rho
        ELSE
            V_c_k(i,j,k) = 0._dp
        END IF
    END DO

    ! get real space potential V(r)
    IF ( .NOT. ALLOCATED ( V_c ) ) ALLOCATE ( V_c ( nfft1 , nfft2 , nfft3 ) )
    in_backward = V_c_k
    CALL dfftw_execute ( plan_backward )
    V_c = out_backward / REAL ( nfft1 * nfft2 * nfft3 , dp )
    
    OPEN(11,FILE='output/V_cmax.dat')
        DO i=1, nfft1
            WRITE(11,*), i , V_c(i, nfft2/2, nfft3/2 )
        END DO
    CLOSE(11)

END SUBROUTINE electrostatic_potential_from_charge_density

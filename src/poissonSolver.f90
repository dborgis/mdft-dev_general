! This SUBROUTINE uses the Gauss theorem associated to the expression of the derivation in Fourier transforms to get V(r).
! Field E(r)=-grad(V(r))
! Local expression of Gauss theorem : div(E(r))=soluteChargeDensity(r)/eps0
! => laplacien(V(r))=-soluteChargeDensity(r)/eps0
! => V(k)=soluteChargeDensity(k)/(esp0*k^2)
! FFT(V(k)) = V(r)

SUBROUTINE poissonSolver (gridnode, gridlen, soluteChargeDensity, Vpoisson)

    USE precision_kinds,    ONLY: dp, i2b, i4b
    USE constants,          ONLY: fourpi , zeroC, twopi
    ! Vpoisson = electrostatic potential from charge density and poisson equation

    IMPLICIT NONE
    
    INTEGER(i2b), INTENT(IN) :: gridnode(3)
    REAL(dp), DIMENSION(gridnode(1),gridnode(2),gridnode(3)), INTENT(IN) :: soluteChargeDensity
    REAL(dp), INTENT(IN) :: gridlen(3)
    REAL(dp), DIMENSION(gridnode(1),gridnode(2),gridnode(3)), INTENT(OUT) :: Vpoisson
    COMPLEX(dp), ALLOCATABLE, DIMENSION(:,:,:) :: soluteChargeDensity_k,Vpoisson_k
    INTEGER (i2b) :: i,j,k,m1,m2,m3
    REAL(dp) :: k2
    REAL(dp), ALLOCATABLE :: fftw3InForward(:,:,:), fftw3OutBackward(:,:,:)
    COMPLEX(dp), ALLOCATABLE :: fftw3OutForward(:,:,:), fftw3InBackward(:,:,:)
    INTEGER(i4b) :: fpspf, fpspb ! fast Poisson Solver Plan Forward or Backward
    INCLUDE "fftw3.f"


    IF ( ALL(soluteChargeDensity==0._dp) ) THEN
        PRINT*,'The solute wears no charge.'
        Vpoisson =0._dp
        RETURN
   
    ELSE
        ALLOCATE( soluteChargeDensity_k (gridnode(1)/2+1,gridnode(2),gridnode(3)) ,SOURCE=zeroC)
        ALLOCATE( Vpoisson_k (gridnode(1)/2+1,gridnode(2),gridnode(3)) ,SOURCE=zeroC)

        CALL prepare_fftw3_for_PoissonGrid

        ! Fourier transform of the solute charge density
        fftw3InForward = soluteChargeDensity 
        CALL dfftw_execute (fpspf)
        soluteChargeDensity_k = fftw3OutForward ! It is verified that at this point, FFT-1(soluteChargeDensity_k)/ (nfft1*nfft2*nfft3) = soluteChargeDensity
        ! FFT(Laplacian(V(r))) = FFT( - 4Pi charge density(r) ) in elecUnits = (ik)^2 V(k) = -4pi rho(k)
        ! V(k) = 4Pi rho(k) / k^2

        DO k = 1, gridnode(3)
            DO j = 1, gridnode(2)
                DO i = 1, gridnode(1)/2+1
            
                    IF ( i<=gridnode(1)/2 ) THEN
                        m1 = i-1
                    ELSE
                        m1 = i-1-gridnode(1)
                    END IF
                    
                    IF ( j<=gridnode(2)/2 ) THEN
                        m2 = j-1
                    ELSE
                        m2 = j-1-gridnode(2)
                    END IF

                    IF ( k<=gridnode(3)/2 ) THEN
                        m3 = k-1
                    ELSE
                        m3 = k-1-gridnode(3)
                    END IF
                    
                    k2 = (twopi/gridlen(1)*REAL(m1,dp))**2 + (twopi/gridlen(2)*REAL(m2,dp))**2 + (twopi/gridlen(3)*REAL(m3,dp))**2

                    IF ( k2 /= 0._dp ) THEN
                        Vpoisson_k(i,j,k) = soluteChargeDensity_k(i,j,k) * fourpi/k2 ! in electrostatic units : V=-4pi rho
                    ELSE
                        Vpoisson_k(i,j,k) = CMPLX(0._dp,0._dp)
                    END IF
                
                END DO
            END DO
        END DO

        ! get real space potential V(r)
        fftw3InBackward = Vpoisson_k
        CALL dfftw_execute (fpspb)
        Vpoisson = fftw3OutBackward / REAL(PRODUCT(gridnode),dp)
    
    END IF

    DEALLOCATE (soluteChargeDensity_k, Vpoisson_k, fftw3InForward, fftw3OutForward, fftw3OutBackward, fftw3InBackward)

    CONTAINS

        !===========================================================================================================================
        SUBROUTINE prepare_fftw3_for_poissongrid
        !===========================================================================================================================
            ! allocate the arrays needed as input for FFT (in_forward) or output for FFT (out_forward)
            ! or needed as input for inverse FFT (in_backward) etc.
            ALLOCATE ( fftw3InForward   ( gridnode(1)      , gridnode(2) , gridnode(3) ) )
            ALLOCATE ( fftw3OutForward  ( gridnode(1)/2 +1 , gridnode(2) , gridnode(3) ) )
            ALLOCATE ( fftw3OutBackward ( gridnode(1)      , gridnode(2) , gridnode(3) ) )
            ALLOCATE ( fftw3InBackward  ( gridnode(1)/2 +1 , gridnode(2) , gridnode(3) ) )
            ! prepare plans needed by fftw3
            CALL dfftw_plan_dft_r2c_3d &
                                ( fpspf, gridnode(1), gridnode(2), gridnode(3), fftw3InForward, fftw3OutForward, FFTW_ESTIMATE )
            CALL dfftw_plan_dft_c2r_3d &
                                ( fpspb, gridnode(1), gridnode(2), gridnode(3), fftw3InBackward, fftw3OutBackward, FFTW_ESTIMATE )
            ! Note that since the fast Poisson solver implies only 1 FFT in each direct, it is useless to use FFTW_MEASURE or even 
            ! more rigorous planning-flags. See http://www.fftw.org/doc/Planner-Flags.html
        END SUBROUTINE prepare_fftw3_for_poissongrid
        !===========================================================================================================================


END SUBROUTINE poissonSolver

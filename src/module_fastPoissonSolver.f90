!===================================================================================================================================
MODULE fastPoissonSolver
!===================================================================================================================================
! This module contains everything related to the fast poisson solver(s) of MDFT.

    USE precision_kinds     ,ONLY: dp,i2b
    USE system              ,ONLY: spaceGrid
    
    IMPLICIT NONE
    
    INTEGER(i2b), PRIVATE :: i
    CHARACTER(180), PRIVATE :: j
    REAL(dp), ALLOCATABLE, DIMENSION(:,:,:), PRIVATE :: soluteChargeDensity, Vpoisson ! TODO THIS ON FINE GRID
    TYPE :: PoissonGrid
        INTEGER(i2b) :: nnod(3)
        REAL(dp) :: len(3)
    END TYPE
    TYPE(PoissonGrid), PRIVATE :: pgrid
    PRIVATE
    PUBLIC :: init

    CONTAINS
    
        !===========================================================================================================================
        SUBROUTINE init
        !===========================================================================================================================
            IMPLICIT NONE
            pgrid%nnod = spaceGrid%n_nodes
            pgrid%len = spaceGrid%length
            IF ( ANY( pgrid%nnod /= spaceGrid%n_nodes )) STOP "Remember that FFTW3 plans are for now only for MDFT's grid"
            ALLOCATE ( soluteChargeDensity (pgrid%nnod(1),pgrid%nnod(2),pgrid%nnod(3)) ,SOURCE=0._dp, STAT=i, ERRMSG=j)
                IF (i/=0) THEN; PRINT*,j; STOP "This problem arises in module fastPoissonSolver"; END IF
            ALLOCATE ( Vpoisson            (pgrid%nnod(1),pgrid%nnod(2),pgrid%nnod(3)) ,SOURCE=0._dp, STAT=i, ERRMSG=j)
                IF (i/=0) THEN; PRINT*,j; STOP "This problem arises in module fastPoissonSolver"; END IF

            CALL soluteChargeDensityFromSoluteChargeCoordinates (pgrid%nnod, pgrid%len, soluteChargeDensity)
            CALL poissonSolver (pgrid%nnod, pgrid%len, soluteChargeDensity, Vpoisson)
            CALL vext_q_from_v_c (pgrid%nnod, pgrid%len, Vpoisson)
        END SUBROUTINE init
        !===========================================================================================================================

        
END MODULE fastPoissonSolver
!===================================================================================================================================

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE fastPoissonSolver
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This module contains everything related to the fast poisson solver(s) of MDFT.
    
    USE precision_kinds     ,ONLY: dp,i2b
    USE system              ,ONLY: spaceGrid
    
    IMPLICIT NONE
    
    INTEGER(i2b), PRIVATE :: i
    CHARACTER(180), PRIVATE :: j
    REAL(dp), ALLOCATABLE, DIMENSION(:,:,:), PRIVATE :: soluteChargeDensity, Vpoisson ! TODO THIS ON FINE GRID


    CONTAINS
    
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE init
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ALLOCATE ( soluteChargeDensity (spaceGrid%n_nodes(1),spaceGrid%n_nodes(2),spaceGrid%n_nodes(3)) ,&
                        SOURCE=0._dp, STAT=i, ERRMSG=j)
                IF (i/=0) THEN; PRINT*,j; STOP "This problem arises in module fastPoissonSolver"; END IF
            ALLOCATE ( Vpoisson            (spaceGrid%n_nodes(1),spaceGrid%n_nodes(2),spaceGrid%n_nodes(3)) ,&
                        SOURCE=0._dp, STAT=i, ERRMSG=j)
                IF (i/=0) THEN; PRINT*,j; STOP "This problem arises in module fastPoissonSolver"; END IF

            CALL soluteChargeDensityFromSoluteChargeCoordinates (soluteChargeDensity)
            CALL poissonSolver (soluteChargeDensity, Vpoisson)
            CALL vext_q_from_v_c (Vpoisson)
        END SUBROUTINE init
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE fastPoissonSolver
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

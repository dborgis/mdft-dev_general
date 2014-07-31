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
    INTEGER(i2b), PRIVATE :: gridnode(3)


    CONTAINS
    
        !===========================================================================================================================
        SUBROUTINE init
        !===========================================================================================================================
            gridnode = spaceGrid%n_nodes
            ALLOCATE ( soluteChargeDensity (gridnode(1),gridnode(2),gridnode(3)) ,SOURCE=0._dp, STAT=i, ERRMSG=j)
                IF (i/=0) THEN; PRINT*,j; STOP "This problem arises in module fastPoissonSolver"; END IF
            ALLOCATE ( Vpoisson            (gridnode(1),gridnode(2),gridnode(3)) ,SOURCE=0._dp, STAT=i, ERRMSG=j)
                IF (i/=0) THEN; PRINT*,j; STOP "This problem arises in module fastPoissonSolver"; END IF

            CALL soluteChargeDensityFromSoluteChargeCoordinates (gridnode, soluteChargeDensity)
            CALL poissonSolver (gridnode, soluteChargeDensity, Vpoisson)
            CALL vext_q_from_v_c (gridnode, Vpoisson)
        END SUBROUTINE init
        !===========================================================================================================================

END MODULE fastPoissonSolver
!===================================================================================================================================

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
            CALL UTest_floorNode
            pgrid%nnod = spaceGrid%n_nodes
            ALLOCATE ( soluteChargeDensity (pgrid%nnod(1),pgrid%nnod(2),pgrid%nnod(3)) ,SOURCE=0._dp, STAT=i, ERRMSG=j)
                IF (i/=0) THEN; PRINT*,j; STOP "This problem arises in module fastPoissonSolver"; END IF
            ALLOCATE ( Vpoisson            (pgrid%nnod(1),pgrid%nnod(2),pgrid%nnod(3)) ,SOURCE=0._dp, STAT=i, ERRMSG=j)
                IF (i/=0) THEN; PRINT*,j; STOP "This problem arises in module fastPoissonSolver"; END IF

            CALL soluteChargeDensityFromSoluteChargeCoordinates (pgrid%nnod, pgrid%len, soluteChargeDensity)
            CALL poissonSolver (pgrid%nnod, pgrid%len, soluteChargeDensity, Vpoisson)
            CALL vext_q_from_v_c (pgrid%nnod, pgrid%len, Vpoisson)
        END SUBROUTINE init
        !===========================================================================================================================


        !===========================================================================================================================
        PURE FUNCTION floorNode(gridnode,gridlen,x,pbc)
        !===========================================================================================================================
            IMPLICIT NONE
            INTEGER(i2b) :: floorNode(3)
            INTEGER(i2b), INTENT(IN) :: gridnode(3)
            REAL(dp), INTENT(IN) :: gridlen(3), x(3)
            LOGICAL, INTENT(IN) :: pbc ! periodic boundary counditions
            REAL(dp) :: dx(3)
            dx = gridlen/REAL(gridnode)
            floorNode = FLOOR(MODULO(x,gridlen)/dx) +1
        END FUNCTION floorNode
        !===========================================================================================================================
        
        
        !===========================================================================================================================
        SUBROUTINE UTest_floorNode
        !===========================================================================================================================
            IMPLICIT NONE
            REAL(dp), PARAMETER :: z=0._dp, o=1.0_dp
            INTEGER(i2b) :: gridnode(3)
            REAL(dp) :: gridlen(3), x(3)
            CALL RANDOM_NUMBER(gridlen)
            IF( ANY( floorNode([1,1,1],gridlen*100,[z,z,z],.TRUE.) /=[1,1,1]) ) STOP "problem in UTest_floorNode"
            IF( ANY( floorNode(INT(gridlen*1000)*10,gridlen*100,[z,z,z],.TRUE.) /=[1,1,1]) ) STOP "problem in UTest_floorNode"
            IF( ANY( floorNode([100,1,1],[50._dp,o,o],[51._dp,z,z],.TRUE.) /=[3,1,1]) ) STOP "problem in UTest_floorNode"
            IF( ANY( floorNode([100,1,1],[50._dp,o,o],[49.999_dp,z,z],.TRUE.) /=[100,1,1]) ) STOP "problem in UTest_floorNode"
        END SUBROUTINE UTest_floorNode
        !===========================================================================================================================
        
        
END MODULE fastPoissonSolver
!===================================================================================================================================

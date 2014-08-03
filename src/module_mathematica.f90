!===================================================================================================================================
MODULE mathematica
!===================================================================================================================================
! This module implements several usefull functions of Mathematica

    IMPLICIT NONE
    PRIVATE
    PUBLIC :: chop, TriLinearInterpolation, UTest_TrilinearInterpolation, floorNode, UTest_floorNode
    
    CONTAINS
    
    !===============================================================================================================================
    PURE FUNCTION chop(x,delta)
    !===============================================================================================================================
    ! see http://reference.wolfram.com/mathematica/ref/Chop.html
    ! It replaces numbers smaller in absolute magnitude than delta by 0.
    ! chop uses a default tolerance of 10._dp**(-10)
        USE precision_kinds ,ONLY: dp
        IMPLICIT NONE
        REAL(dp) :: chop
        REAL(dp), INTENT(IN) :: x
        REAL(dp), OPTIONAL, INTENT(IN) :: delta
        REAL(dp), PARAMETER :: defaultdelta=10._dp**(-10)
        REAL(dp) :: d
        IF (PRESENT(delta)) THEN
            d=delta
        ELSE
            d=defaultdelta
        END IF
        IF (x<=d) THEN
            chop=0._dp
        ELSE
            chop=x
        END IF
    END FUNCTION chop
    !===============================================================================================================================
    
    
    !===============================================================================================================================
    PURE FUNCTION TriLinearInterpolation (cube,x)
    !===============================================================================================================================
    ! Returns the value at position x(1:3) within the cube. The value is known at each corner of the cube.
    ! It is thus an interpolation of value at the corners the cube to a point inside the cube.
        USE precision_kinds ,ONLY: dp
        IMPLICIT NONE
        REAL(dp), INTENT(IN) :: cube(0:1,0:1,0:1), x(1:3)
        REAL(dp) :: TriLinearInterpolation
        IF( ALL(cube==cube(0,0,0)) ) THEN ! homogeneous case
            TrilinearInterpolation = cube(0,0,0)
        ELSE
            TriLinearInterpolation = cube(0,0,0) * (1._dp-x(1)) * (1._dp-x(2)) * (1._dp-x(3)) &
                                    +cube(1,0,0) * x(1) * (1._dp-x(2)) * (1._dp-x(3)) &
                                    +cube(0,1,0) * (1._dp-x(1)) * x(2) * (1._dp-x(3)) &
                                    +cube(0,0,1) * (1._dp-x(1)) * (1._dp-x(2)) * x(3) &
                                    +cube(1,0,1) * x(1) * (1._dp-x(2)) * x(3) &
                                    +cube(0,1,1) * (1._dp-x(1)) * x(2) * x(3) &
                                    +cube(1,1,0) * x(1) * x(2) * (1._dp-x(3)) &
                                    +cube(1,1,1) * x(1) * x(2) * x(3)
        END IF
    END FUNCTION TriLinearInterpolation
    !===============================================================================================================================
    
    
    !===============================================================================================================================
    SUBROUTINE UTest_TrilinearInterpolation
    !===============================================================================================================================
    ! Tests the pure function TriLinearInterpolation where result is known:
    ! - if the point is one of the corners
    ! - if it is on the center of the cube
    ! Then it tests that no answer is higher or lower than the maximum or minimum value of any corner.
        USE precision_kinds     ,ONLY:dp
        IMPLICIT NONE
        REAL(dp) :: A(0:1,0:1,0:1), x(1:3)
        REAL(dp), PARAMETER :: z=0._dp, o=1.0_dp
        LOGICAL, SAVE :: alreadydone=.FALSE.
        REAL(dp) :: cube(0:1,0:1,0:1), t0, t1
        IF (alreadydone) RETURN
        CALL CPU_TIME(t0)
        CALL RANDOM_NUMBER(cube)
        cube = cube * 1000._dp
        IF( TriLinearInterpolation(cube,[z,z,z]) /= cube(0,0,0) .OR.&
            TriLinearInterpolation(cube,[o,z,z]) /= cube(1,0,0) .OR.&
            TriLinearInterpolation(cube,[z,o,z]) /= cube(0,1,0) .OR.&
            TriLinearInterpolation(cube,[z,z,o]) /= cube(0,0,1) .OR.&
            TriLinearInterpolation(cube,[o,o,z]) /= cube(1,1,0) .OR.&
            TriLinearInterpolation(cube,[o,z,o]) /= cube(1,0,1) .OR.&
            TriLinearInterpolation(cube,[z,o,o]) /= cube(0,1,1) .OR.&
            TriLinearInterpolation(cube,[o,o,o]) /= cube(1,1,1) .OR.&
            ABS(TriLinearInterpolation(cube,[o,o,o]/2._dp) - SUM(cube)/8._dp )>EPSILON(1.0_dp) ) THEN
                STOP "Problem detected in UTest_TriLinearInterpolation"
        END IF
        CALL CPU_TIME(t1)
        DO WHILE(t1-t0<0.005_dp) ! Test for 5 ms
            CALL RANDOM_NUMBER(x)
            CALL RANDOM_NUMBER(cube)
            IF ( TriLinearInterpolation(cube,x) < MINVAL(cube) &
                .OR. TriLinearInterpolation(cube,x) > MAXVAL(cube) ) STOP "Problem detected in UTest_TriLinearInterpolation"
            CALL CPU_TIME(t1)
        END DO
        alreadydone = .TRUE.
    END SUBROUTINE UTest_TrilinearInterpolation
    !===============================================================================================================================


    !===========================================================================================================================
    PURE FUNCTION floorNode(gridnode,gridlen,x,pbc)
    !===========================================================================================================================
        USE precision_kinds, ONLY: i2b, dp
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
        USE precision_kinds, ONLY: i2b, dp
        IMPLICIT NONE
        LOGICAL, SAVE :: alreadydone=.FALSE.
        REAL(dp), PARAMETER :: z=0._dp, o=1.0_dp
        INTEGER(i2b) :: gridnode(3)
        REAL(dp) :: gridlen(3), x(3)
        IF (alreadydone) RETURN
        CALL RANDOM_NUMBER(gridlen)
        IF( ANY( floorNode([1,1,1],gridlen*100,[z,z,z],.TRUE.) /=[1,1,1]) ) STOP "problem in UTest_floorNode"
        IF( ANY( floorNode(INT(gridlen*1000)*10,gridlen*100,[z,z,z],.TRUE.) /=[1,1,1]) ) STOP "problem in UTest_floorNode"
        IF( ANY( floorNode([100,1,1],[50._dp,o,o],[51._dp,z,z],.TRUE.) /=[3,1,1]) ) STOP "problem in UTest_floorNode"
        IF( ANY( floorNode([100,1,1],[50._dp,o,o],[49.999_dp,z,z],.TRUE.) /=[100,1,1]) ) STOP "problem in UTest_floorNode"
        alreadydone=.TRUE.
    END SUBROUTINE UTest_floorNode
    !===========================================================================================================================
        

END MODULE
!===================================================================================================================================

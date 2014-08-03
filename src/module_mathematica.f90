!===================================================================================================================================
MODULE mathematica
!===================================================================================================================================
! This module implements several usefull functions of Mathematica

    IMPLICIT NONE
    PRIVATE
    PUBLIC :: chop, TriLinearInterpolation, UTest_TrilinearInterpolation
    
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
        USE precision_kinds     ,ONLY:dp,i2b
        IMPLICIT NONE
        REAL(dp) :: A(0:1,0:1,0:1), x(1:3)
        REAL(dp), PARAMETER :: z=0._dp, o=1.0_dp
        LOGICAL, SAVE :: alreadydone=.FALSE.
        REAL(dp) :: cube(0:1,0:1,0:1)
        INTEGER(i2b) :: i
        IF (alreadydone) RETURN
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
        DO i=1,10000 ! 10000 has been chosen so that the whole execution time of UTest_TrilinearInterpolation is 1ms on my laptop.
            CALL RANDOM_NUMBER(x)
            IF ( TriLinearInterpolation(cube,x) < MINVAL(cube) &
                .OR. TriLinearInterpolation(cube,x) > MAXVAL(cube) ) STOP "Problem detected in UTest_TriLinearInterpolation"
        END DO
        alreadydone = .TRUE.
    END SUBROUTINE UTest_TrilinearInterpolation
    !===============================================================================================================================

END MODULE
!===================================================================================================================================

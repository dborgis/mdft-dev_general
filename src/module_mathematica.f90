!===================================================================================================================================
MODULE mathematica
!===================================================================================================================================
! This module implements several usefull functions of Mathematica

    IMPLICIT NONE
    PRIVATE
    PUBLIC :: chop, TriLinearInterpolation
    
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

END MODULE
!===================================================================================================================================

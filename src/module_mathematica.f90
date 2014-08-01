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
    PURE FUNCTION TriLinearInterpolation (cube,r)
    !===============================================================================================================================
    ! Returns the value at position r(1:3) within the cube. The value is known at each corner of the cube.
    ! It is thus an interpolation of value at the corners the cube to a point inside the cube.
        USE precision_kinds ,ONLY: dp
        IMPLICIT NONE
        REAL(dp), INTENT(IN) :: cube(0:1,0:1,0:1), r(1:3)
        REAL(dp) :: TriLinearInterpolation
        TriLinearInterpolation = cube(0,0,0) * (1._dp-r(1)) * (1._dp-r(2)) * (1._dp-r(3)) &
                                +cube(1,0,0) * r(1) * (1._dp-r(2)) * (1._dp-r(3)) &
                                +cube(0,1,0) * (1._dp-r(1)) * r(2) * (1._dp-r(3)) &
                                +cube(0,0,1) * (1._dp-r(1)) * (1._dp-r(2)) * r(3) &
                                +cube(1,0,1) * r(1) * (1._dp-r(2)) * r(3) &
                                +cube(0,1,1) * (1._dp-r(1)) * r(2) * r(3) &
                                +cube(1,1,0) * r(1) * r(2) * (1._dp-r(3)) &
                                +cube(1,1,1) * r(1) * r(2) * r(3)
    END FUNCTION TriLinearInterpolation
    !===============================================================================================================================

END MODULE
!===================================================================================================================================

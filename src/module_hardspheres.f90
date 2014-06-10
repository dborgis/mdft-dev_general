MODULE hardspheres

    USE precision_kinds  ,ONLY: dp
    IMPLICIT NONE
    
    CONTAINS
    
    
        ! packing fraction : eta = 4/3 * pi * R^3 * solvent density of the constituant ( /= total solvent density)
        PURE REAL(dp) FUNCTION packfrac (n,R)
            REAL(dp), INTENT(IN) :: n ! density of the constituant (/= total solvant density)
            REAL(dp), INTENT(IN) :: R ! radius of the constituant
            packfrac = 4._dp/3._dp * ACOS(-1._dp) * R**3 * n
        END FUNCTION packfrac

END MODULE hardspheres

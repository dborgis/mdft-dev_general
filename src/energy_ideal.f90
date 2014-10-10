! This SUBROUTINE computes the ideal part of the free energy functional.
SUBROUTINE energy_ideal (Fideal)

    USE precision_kinds, ONLY: i2b, dp
    USE minimizer, ONLY: cg_vect, FF, dF
    USE system, ONLY: thermocond, spaceGrid, solvent
    USE quadrature, ONLY: molRotSymOrder, angGrid, molRotGrid

    IMPLICIT NONE
    
    REAL(dp), INTENT(OUT) :: Fideal
    INTEGER(i2b) :: icg , i , j , k , o , p, s! dummy for loops
    REAL(dp) :: psi ! dummy for cg_vext(i)
    REAL(dp) :: rho, rhon ! local density
    REAL(dp) :: logrho ! dummy for log(rho)
    REAL(dp) :: time0, time1
    
    CALL CPU_TIME (time0) ! init timer

    Fideal = 0.0_dp! init Fideal to zero and its gradient

    icg = 0
    do s = 1 , size(solvent)
        do i = 1 , spaceGrid%n_nodes(1)
            do j = 1 , spaceGrid%n_nodes(2)
                do k = 1 , spaceGrid%n_nodes(3)
                    do o = 1 , angGrid%n_angles
                        do p = 1 , molRotGrid%n_angles
                            icg = icg + 1
                            psi = CG_vect (icg)
                            rho = psi**2
                            Fideal = Fideal + Fideal_local (o,p,s,rho)
                            dF (icg) = dF (icg) + dFideal_local (o,p,s,psi,0._dp)
                        END DO
                    END DO
                END DO
            END DO
        END DO
    END DO

    Fideal = Fideal * thermocond%kbT * spaceGrid%dv ! integration factor
    FF = FF + Fideal
    
    CALL CPU_TIME (time1)

    CONTAINS

!===================================================================================================================================
    
    PURE FUNCTION dFideal_local (o,p,s,psi,toadd)
        INTEGER(i2b), INTENT(IN) :: o,p,s
        REAL(dp), INTENT(IN) :: psi, toadd
        REAL(dp) :: dFideal_local
        IF (abs(psi) > epsilon(1._dp)) THEN
            dFideal_local = 2.0_dp * psi * prefactor(o,p,s) * spaceGrid%dv * ( thermocond%kbT*LOG(psi**2) + toadd )
        ELSE
            dFideal_local = 0._dp
        END IF
    END FUNCTION dFideal_local

!===================================================================================================================================
    
    PURE FUNCTION Fideal_local (o,p,s,rho)
        INTEGER(i2b), INTENT(IN) :: o,p,s
        REAL(dp), INTENT(IN) :: rho
        REAL(dp) :: Fideal_local
        IF (abs(rho) > epsilon(1.0_dp) ) THEN
            Fideal_local = prefactor(o,p,s) * (rho*LOG(rho)-rho+1.0_dp)
        ELSE
            Fideal_local = prefactor(o,p,s)
        END IF
    END FUNCTION Fideal_local

!===================================================================================================================================

    PURE FUNCTION prefactor (o,p,s)
        INTEGER(i2b), INTENT(IN) :: o,p,s
        REAL(dp) :: prefactor
        prefactor = angGrid%weight(o) * molRotGrid%weight(p) * solvent(s)%rho0
    END FUNCTION

!===================================================================================================================================

END SUBROUTINE energy_ideal

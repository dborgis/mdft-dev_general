! This SUBROUTINE computes the ideal part of the free energy functional.
SUBROUTINE energy_ideal (Fideal)

    USE precision_kinds, ONLY: dp
    USE minimizer, ONLY: cg_vect_new, FF, dF_new
    USE system, ONLY: thermocond, spaceGrid, solvent, nb_species
    USE quadrature, ONLY: molRotSymOrder, angGrid, molRotGrid

    IMPLICIT NONE

    REAL(dp), INTENT(OUT) :: Fideal
    INTEGER :: icg , i , j , k , o , p, s, nfft1, nfft2, nfft3
    REAL(dp) :: psi ! dummy for cg_vext(i)
    REAL(dp) :: rho, rhon ! local density
    REAL(dp) :: logrho ! dummy for log(rho)
    REAL(dp) :: time0, time1

    nfft1 = spacegrid%n_nodes(1)
    nfft2 = spacegrid%n_nodes(2)
    nfft3 = spacegrid%n_nodes(3)

    CALL CPU_TIME (time0) ! init timer

    Fideal = 0.0_dp! init Fideal to zero and its gradient
    do concurrent( i=1:nfft1, j=1:nfft2, k=1:nfft3, o=1:anggrid%n_angles, p=1:molrotgrid%n_angles, s=1:nb_species )
      psi = cg_vect_new(i,j,k,o,p,s)
      rho = psi**2
      Fideal = Fideal + Fideal_local (o,p,s,rho)
      dF_new(i,j,k,o,p,s) = dF_new(i,j,k,o,p,s) + dFideal_local (o,p,s,psi,0._dp)
    end do

    Fideal = Fideal * thermocond%kbT * spaceGrid%dv ! integration factor
    FF = FF + Fideal

    CALL CPU_TIME (time1)
    print*, "IDEAL TERM BENCH:", time1-time0

    CONTAINS

!===================================================================================================================================

    PURE FUNCTION dFideal_local (o,p,s,psi,toadd)
        INTEGER, INTENT(IN) :: o,p,s
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
        INTEGER, INTENT(IN) :: o,p,s
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
        INTEGER, INTENT(IN) :: o,p,s
        REAL(dp) :: prefactor
        prefactor = angGrid%weight(o) * molRotGrid%weight(p) * solvent(s)%rho0
    END FUNCTION

!===================================================================================================================================

END SUBROUTINE energy_ideal

! This SUBROUTINE computes the ideal part of the free energy functional.
SUBROUTINE energy_ideal (Fideal)

    use precision_kinds, ONLY: dp
    use minimizer, ONLY: cg_vect_new, FF, dF_new
    use system, ONLY: thermocond, solvent
    use module_grid, only: grid

    IMPLICIT NONE

    REAL(dp), INTENT(OUT) :: Fideal
    INTEGER :: icg, i, j, k, o, p, s, nx, ny, nz, io
    REAL(dp) :: psi ! dummy for cg_vext(i)
    REAL(dp) :: rho, rhon ! local density
    REAL(dp) :: logrho ! dummy for log(rho)
    REAL(dp) :: time0, time1

    nx = grid%nx
    ny = grid%ny
    nz = grid%nz

    CALL CPU_TIME (time0) ! init timer

    Fideal = 0.0_dp! init Fideal to zero and its gradient
    do concurrent( i=1:nx, j=1:ny, k=1:nz, io=1:grid%no, s=1:solvent(1)%nspec )
      psi = cg_vect_new(i,j,k,io,s)
      rho = psi**2
      Fideal = Fideal + Fideal_local (io,s,rho)
      dF_new(i,j,k,io,s) = dF_new(i,j,k,io,s) + dFideal_local (io,s,psi)
    end do

    Fideal = Fideal * thermocond%kbT * grid%dv ! integration factor
    FF = FF + Fideal

    CALL CPU_TIME (time1)
    print*, "IDEAL TERM BENCH:", time1-time0

    CONTAINS

!===================================================================================================================================

    PURE FUNCTION dFideal_local (io,s,psi)
        INTEGER, INTENT(IN) :: io, s
        REAL(dp), INTENT(IN) :: psi
        REAL(dp) :: dFideal_local
        IF (abs(psi) > epsilon(1._dp)) THEN
            dFideal_local = 2.0_dp * psi * prefactor(io,s) * grid%dv * ( thermocond%kbT*LOG(psi**2) )
        ELSE
            dFideal_local = 0._dp
        END IF
    END FUNCTION dFideal_local

!===================================================================================================================================

    PURE FUNCTION Fideal_local (io,s,rho)
        INTEGER, INTENT(IN) :: io, s
        REAL(dp), INTENT(IN) :: rho
        REAL(dp) :: Fideal_local
        IF (abs(rho) > epsilon(1.0_dp) ) THEN
            Fideal_local = prefactor(io,s) * (rho*LOG(rho)-rho+1.0_dp)
        ELSE
            Fideal_local = prefactor(io,s)
        END IF
    END FUNCTION Fideal_local

!===================================================================================================================================

    PURE FUNCTION prefactor (io,s)
        INTEGER, INTENT(IN) :: io,s
        REAL(dp) :: prefactor
        prefactor = grid%w(io) * solvent(s)%rho0
    END FUNCTION

!===================================================================================================================================

END SUBROUTINE energy_ideal

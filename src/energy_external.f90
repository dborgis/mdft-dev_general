! this SUBROUTINE compute the external part of the free energy functional
SUBROUTINE energy_external (Fext)

    USE precision_kinds, ONLY: dp , i2b
    USE system, ONLY: rhoBulk=>rho_0_multispec, nb_species, spaceGrid
    USE quadrature, ONLY: angGrid, molRotGrid
    USE minimizer, ONLY: CG_vect , FF , dF
    USE external_potential, ONLY: Vext_total
    
    IMPLICIT NONE

    REAL(dp), INTENT(OUT) :: Fext
    INTEGER(i2b) :: icg , i , j , k , o , p! dummy for loops
    REAL(dp) :: psi
    REAL(dp) :: wdfve
    REAL(dp) :: time0, time1
    INTEGER(i2b) :: spec
    
    CALL CPU_TIME ( time0 )

    Fext = 0.0_dp

    ! F_{ext}[\rho(\vec{r},\vec{\Omega})]=\int d \vec{r} d \vec{\Omega} V_{ext}(\vec{r},\vec{\Omega})\rho(\vec{r},\vec{\Omega})

    icg = 0
    DO spec = 1 , nb_species
        DO i = 1 , spaceGrid%n_nodes(1)
        DO j = 1 , spaceGrid%n_nodes(2)
        DO k = 1 , spaceGrid%n_nodes(3)
            DO o = 1 , angGrid%n_angles
                DO p=1 , molRotGrid%n_angles            
                    icg = icg + 1
                    psi = cg_vect(icg)
                    wdfve = angGrid%weight(o) * molRotGrid%weight(p) * spaceGrid%dv * rhoBulk(spec) * Vext_total(i,j,k,o,p,spec)
                    Fext = Fext + psi**2 * wdfve
                    dF(icg) = dF(icg) + 2.0_dp*psi*wdfve
                END DO
            END DO
        END DO
        END DO
        END DO
    END DO

    FF = FF + Fext

    CALL CPU_TIME (time1)
END SUBROUTINE energy_external

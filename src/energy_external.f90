! this SUBROUTINE compute the external part of the free energy functional
SUBROUTINE energy_external (Fext)

    USE precision_kinds, ONLY: dp , i2b
    USE system, ONLY: nb_species, spaceGrid, solvent
    USE quadrature, ONLY: angGrid, molRotGrid
    USE minimizer, ONLY: CG_vect , FF , dF
    USE external_potential, ONLY: Vext_total
    USE input, ONLY: input_dp

    IMPLICIT NONE

    REAL(dp), INTENT(OUT) :: Fext
    INTEGER(i2b) :: icg , i , j , k , o , p,s
    REAL(dp) :: psi, wdfve, imposedChemPot

    Fext = 0.0_dp

    ! Impose a chemical potential
    imposedChemPot = input_dp('imposed_chempot',0._dp) ! 0 means that you don't impose any chemical potential (delta_chempot is built this way)
    IF( nb_species/=1 .AND. imposedChemPot/=0._dp) STOP "Imposing a chemical potential is valid only for single-species solvent"

    ! F_{ext}[\rho(\vec{r},\vec{\Omega})]=\int d \vec{r} d \vec{\Omega} V_{ext}(\vec{r},\vec{\Omega})\rho(\vec{r},\vec{\Omega})

    icg = 0
    DO s = 1 , nb_species
        DO i = 1 , spaceGrid%n_nodes(1)
        DO j = 1 , spaceGrid%n_nodes(2)
        DO k = 1 , spaceGrid%n_nodes(3)
            DO o = 1 , angGrid%n_angles
                DO p=1 , molRotGrid%n_angles
                    icg = icg + 1
                    psi = cg_vect(icg)
                    wdfve = angGrid%weight(o) * molRotGrid%weight(p) * spaceGrid%dv * solvent(s)%rho0 &
                            * (Vext_total(i,j,k,o,p,s) - imposedChemPot)
                    Fext = Fext + psi**2 * wdfve
                    dF(icg) = dF(icg) + 2.0_dp*psi*wdfve
                END DO
            END DO
        END DO
        END DO
        END DO
    END DO

    FF = FF + Fext

END SUBROUTINE energy_external

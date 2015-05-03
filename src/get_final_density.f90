! Gets the last density from the minimizer
PURE SUBROUTINE get_final_density ( neq , solventspecies)

    USE precision_kinds, ONLY: dp , i2b
    USE system,          ONLY: spaceGrid
    USE constants,       ONLY: twopi
    USE minimizer,       ONLY: CG_vect_new
    USE quadrature,      ONLY: molRotSymOrder, angGrid, molRotGrid
    USE fft,             ONLY: fftw3 , timesExpPrefactork2
    USE input,           ONLY: verbose
    
    IMPLICIT NONE

    INTEGER(i2b) :: i,j,k,o,p,s,icg
    integer(i2b), intent(in) :: solventspecies ! the solvent species of which we want the number density
    REAL(dp), INTENT(OUT) :: neq (spaceGrid%n_nodes(1),spaceGrid%n_nodes(2),spaceGrid%n_nodes(3)) ! equilibrium density(position)
    REAL(dp) :: rho,local_density
    
    neq = 0._dp
    icg = 0
    DO s =1, solventspecies
        DO i =1, spacegrid%n_nodes(1)
            DO j =1, spacegrid%n_nodes(2)
                DO k =1, spacegrid%n_nodes(3)
                    local_density = 0.0_dp
                    DO o =1,angGrid%n_angles
                        DO p =1,molRotGrid%n_angles
                            icg = icg + 1
                            rho = cg_vect_new(i,j,k,o,p,s) ** 2 / (twopi**2 * 2.0_dp/molRotSymOrder)
                            local_density = local_density + angGrid%weight(o) * molRotGrid%weight(p)* rho ! integral of rho over all orientations ie 'n'
                        END DO
                    END DO
                    if (s == solventspecies) neq(i,j,k) = local_density
                END DO
            END DO
        END DO
    END DO

END SUBROUTINE get_final_density

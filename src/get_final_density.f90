! Gets the last density from the minimizer
PURE SUBROUTINE get_final_density ( neq , solventspecies)

    USE precision_kinds, ONLY: dp , i2b
    USE system,          ONLY: spaceGrid, nb_species
    USE minimizer,       ONLY: CG_vect_new
    USE fft,             ONLY: fftw3 , timesExpPrefactork2
    USE input,           ONLY: verbose

    IMPLICIT NONE

    INTEGER(i2b) :: i,j,k,s,icg,io
    integer(i2b), intent(in) :: solventspecies ! the solvent species of which we want the number density
    REAL(dp), INTENT(OUT) :: neq (spaceGrid%n_nodes(1),spaceGrid%n_nodes(2),spaceGrid%n_nodes(3)) ! equilibrium density(position)
    REAL(dp) :: rho,local_density
    real(dp), parameter :: twopi = 2._dp*acos(-1._dp), twopi2 = twopi**2

    neq = 0._dp
    icg = 0
    DO s =1, nb_species
        DO i =1, spacegrid%nx
            DO j =1, spacegrid%ny
                DO k =1, spacegrid%nz
                    local_density = 0.0_dp
                    DO io =1,spacegrid%no
                        icg = icg + 1
                        rho = cg_vect_new(i,j,k,io,s)**2 / (twopi2 * 2.0_dp/spacegrid%molRotSymOrder)
                        local_density = local_density + spacegrid%w(io)* rho ! integral of rho over all orientations ie 'n'
                    END DO
                    if (s == solventspecies) neq(i,j,k) = local_density
                END DO
            END DO
        END DO
    END DO

END SUBROUTINE get_final_density

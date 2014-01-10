! this subroutine gets the partial vext ( hard sphere + hard wall + hard cylinder + purely repulsive)
! and adds to it Vext_lj and Vext_q.
! then it gives a upper value (100 kJ/mol) to vext_total.
SUBROUTINE vext_total_sum

    USE precision_kinds,    ONLY: dp, i2b
    USE system,             ONLY: nfft1, nfft2, nfft3, nb_species
    USE quadrature,         ONLY: angGrid, molRotGrid
    USE constants,          ONLY: fourpi, zero
    USE external_potential, ONLY: Vext_total, Vext_lj, Vext_q, vext_hard_core
    USE quadrature,         ONLY: sym_order
    USE input,              ONLY: verbose

    IMPLICIT NONE
    REAL(dp), PARAMETER :: vmax = 100._dp
    
   
    IF ( .NOT. ALLOCATED ( Vext_total ) ) THEN ! be sure Vext_total is allocated
        ALLOCATE ( Vext_total ( nfft1 , nfft2 , nfft3 , angGrid%n_angles , molRotGrid%n_angles, nb_species ), SOURCE=zero )
    END IF

    ! Vext_total = 0.0_dp
    ! vext is the sum over all external potentials
    ! note that purely repulsive and hard potentials are already included in vext
    IF ( ALLOCATED ( Vext_q ) )         Vext_total = Vext_total + Vext_q
    IF ( ALLOCATED ( Vext_lj ) )        Vext_total = Vext_total + Vext_lj
    IF ( ALLOCATED ( vext_hard_core ) ) vext_total (:,:,:,1,1,:) = vext_total (:,:,:,1,1,:) + vext_hard_core ! TODO generalize

    IF ( ALL(vext_total==zero) ) STOP "The external potential is zero everywhere. Something's wrong in input files"

    WHERE ( Vext_total > vmax ) Vext_total = vmax
    IF ( ALLOCATED ( vext_lj        ) ) WHERE ( Vext_lj        > vmax ) Vext_lj        = vmax
    IF ( ALLOCATED ( vext_hard_core ) ) WHERE ( vext_hard_core > vmax ) vext_hard_core = vmax

    IF (verbose) THEN
        BLOCK
            CHARACTER(50) :: filename ! dummy
            REAL(dp), DIMENSION (nfft1,nfft2,nfft3) :: temparray ! dummy
            PRINT*, MINVAL ( Vext_total ) , ' < Vext_total < ' , MAXVAL(Vext_total)! give vext extrema to user for visual debugging
            PRINT*, MINVAL ( Vext_lj ) , ' < Vext_lj < ' , MAXVAL( Vext_lj )
            IF ( ALLOCATED (Vext_q) ) PRINT*, MINVAL(Vext_q), ' < Vext_q < ' , MAXVAL(Vext_q)
            IF ( ALLOCATED (Vext_hard_core) ) PRINT*, MINVAL(Vext_hard_core), ' < Vext_hard_core < ', MAXVAL(Vext_hard_core)
            ! mean over orientations and print
            CALL mean_over_orientations ( Vext_total ( : , : , : , : , : , 1 ) , temparray )
            temparray = temparray / (fourpi**2/(2.0_dp*sym_order))
            filename = 'output/Vext.cube' ! care when HS or multispec
            CALL write_to_cube_file ( temparray , filename )
            filename = 'output/Vext_along-z.dat' ! care when HS or multispec'
            CALL compute_z_density ( temparray , filename )
        END BLOCK
    END IF

END SUBROUTINE vext_total_sum

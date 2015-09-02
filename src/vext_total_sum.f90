! this subroutine gets the partial vext ( hard sphere + hard wall + hard cylinder + purely repulsive)
! and adds to it Vext_lj and Vext_q.
! then it gives a upper value (100 kJ/mol) to vext_total.
SUBROUTINE vext_total_sum

    use precision_kinds,    ONLY: dp, i2b
    use system,             ONLY: solvent
    use module_grid, only: grid
    use external_potential, ONLY: Vext_total, Vext_lj, Vext_q, vext_hard_core
    use module_input,              ONLY: verbose

    IMPLICIT NONE

    real(dp), parameter :: vmax = HUGE(1.0_dp)
    real(dp), parameter :: fourpi=4._dp*acos(-1._dp), zero=0._dp
    integer :: nx, ny, nz, no, ns


    nx = grid%nx
    ny = grid%ny
    nz = grid%nz
    no = grid%no
    ns = solvent(1)%nspec


    IF ( .NOT. ALLOCATED ( Vext_total ) ) THEN ! be sure Vext_total is allocated
        ALLOCATE ( Vext_total (nx,ny,nz,no,ns), SOURCE=zero )
    END IF

    ! Vext_total = 0.0_dp
    ! vext is the sum over all external potentials
    ! note that purely repulsive and hard potentials are already included in vext
    IF ( ALLOCATED ( Vext_q ) )         Vext_total = Vext_total + Vext_q
    IF ( ALLOCATED ( Vext_lj ) )        Vext_total = Vext_total + Vext_lj
    IF ( ALLOCATED ( vext_hard_core ) ) stop "vext_hard_core not ok"!vext_total (:,:,:,1,:) = vext_total (:,:,:,1,:) + vext_hard_core ! TODO generalize

    IF ( ALL(vext_total==zero) ) PRINT*, "WARNING: The external potential is zero everywhere. Is it what you want?"

    IF ( ANY(vext_total/=vext_total) ) STOP "There is a NaN somewhere in vext_total."

    IF ( ANY(ABS(vext_total)>HUGE(1.0_dp)) ) STOP "There is an Infinity somewhere in vext_total."

    WHERE ( Vext_total > vmax ) Vext_total = vmax
    IF ( ALLOCATED ( vext_lj        ) ) THEN
        WHERE ( Vext_lj >= vmax )
            Vext_lj = vmax
            Vext_total = vmax
        END WHERE
    END IF
    IF ( ALLOCATED ( vext_hard_core ) ) WHERE ( vext_hard_core > vmax ) vext_hard_core = vmax

    ! IF (verbose) THEN
    !     BLOCK
    !         CHARACTER(50) :: filename ! dummy
    !         REAL(dp), DIMENSION (nx,ny,nz) :: temparray ! dummy
    !         PRINT*, MINVAL ( Vext_total ) , ' < Vext_total < ' , MAXVAL(Vext_total)! give vext extrema to user for visual debugging
    !         PRINT*, MINVAL ( Vext_lj ) , ' < Vext_lj < ' , MAXVAL( Vext_lj )
    !         IF ( ALLOCATED (Vext_q) ) PRINT*, MINVAL(Vext_q), ' < Vext_q < ' , MAXVAL(Vext_q)
    !         IF ( ALLOCATED (Vext_hard_core) ) PRINT*, MINVAL(Vext_hard_core), ' < Vext_hard_core < ', MAXVAL(Vext_hard_core)
    !         ! mean over orientations and print
    !         CALL mean_over_orientations ( Vext_total ( : , : , : , : , : , 1 ) , temparray )
    !         temparray = temparray / (fourpi**2/(2.0_dp*molRotSymOrder))
    !         filename = 'output/Vext.cube' ! care when HS or multispec
    !         CALL write_to_cube_file ( temparray , filename )
    !         filename = 'output/z_Vext.dat' ! care when HS or multispec'
    !         CALL compute_z_density ( temparray , filename )
    !     END BLOCK
    ! END IF

END SUBROUTINE vext_total_sum

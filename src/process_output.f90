! This subroutine computes every output asked by user.

SUBROUTINE process_output

    ! defines precision of reals and intergers
    USE precision_kinds,    ONLY: dp, i2b
    USE system,             ONLY: nb_species, spaceGrid, solvent
    USE input,              ONLY: verbose, input_log, input_char
    USE solute_geometry,    ONLY: soluteIsPlanar => isPlanar, soluteIsLinear => isLinear
    USE constants,          ONLY: zerodp
    USE minimizer,          ONLY : finalizeminimizer

    IMPLICIT NONE

    CHARACTER(50):: filename
    REAL(dp), ALLOCATABLE , DIMENSION (:,:,:,:) :: neq, Px, Py, Pz ! equilibrium density, ie rho(r), and Pi polarization(r)
    INTEGER(i2b) ,POINTER :: nfft1=>spaceGrid%n_nodes(1), nfft2=>spaceGrid%n_nodes(2), nfft3=>spaceGrid%n_nodes(3)
    INTEGER(i2b) :: s

    CALL print_cg_vect ! print output/density.bin that contains cg_vect

    allocate ( neq (nfft1,nfft2,nfft3,nb_species) ,SOURCE=zerodp)
    do s=1,size(solvent)
        call get_final_density (neq,s)
    end do

    DO s=1,nb_species
        PRINT*,"NUMBER OF PARTICLES OF SPECIES",s,"IN SUPERCELL =",SUM(neq)*spaceGrid%dv *solvent(s)%n0
    END DO


    IF (verbose) THEN
        filename = 'output/density.cube'
        CALL write_to_cube_file (neq(:,:,:,1), filename) ! TODO for now only write for the first species
        IF (input_log('polarization').EQV. .TRUE.) THEN
            ALLOCATE ( Px (nfft1,nfft2,nfft3,nb_species) ,SOURCE=zerodp)
            ALLOCATE ( Py (nfft1,nfft2,nfft3,nb_species) ,SOURCE=zerodp)
            ALLOCATE ( Pz (nfft1,nfft2,nfft3,nb_species) ,SOURCE=zerodp)
            CALL get_final_polarization (Px,Py,Pz)
            filename='output/normP.cube' ; CALL write_to_cube_file ( (SQRT(Px(:,:,:,1)**2+Py(:,:,:,1)**2+Pz(:,:,:,1)**2)), filename) ! TODO for now only write for the first species
            filename='output/Px.cube' ; CALL write_to_cube_file(Px(:,:,:,1), filename)
            filename='output/z_Px.dat'; CALL compute_z_density(Px(:,:,:,1) , filename)
            DEALLOCATE (Px)
            filename='output/Py.cube' ; CALL write_to_cube_file(Py(:,:,:,1), filename)
            filename='output/z_Py.dat'; CALL compute_z_density(Py(:,:,:,1) , filename)
            DEALLOCATE (Py)
            filename='output/Pz.cube' ; CALL write_to_cube_file(Pz(:,:,:,1), filename)
            filename='output/z_Pz.dat'; CALL compute_z_density(Pz(:,:,:,1) , filename)
            DEALLOCATE (Pz)
        END IF


        ! If calculation is for hard sphere fluid in presence of a hard wall compute profile perp wall
        ! TODO: DONT HAVE TIME TO WRITE THE TEST TODAY
        filename = 'output/z_density.out'
        CALL compute_z_density ( neq (:,:,:,1) , filename ) ! TODO for now only write for the first species

        IF (soluteIsLinear()) THEN
            ! nothing for now
        END IF

        IF( soluteIsPlanar() ) THEN
            PRINT*,'This solute has planar symetry'
            filename = 'output/planardensity.out'
            CALL compute_planar_density ( neq (:,:,:,1) , filename ) ! TODO for now only write for the first species
        END IF

        IF ( input_char('other_predefined_vext')=='vextdef0' ) THEN
            filename = 'output/molecular_density_in_xy_plane.out'
            PRINT*,"I am writing file ",filename
            OPEN(378,FILE=filename)
                BLOCK
                    INTEGER(i2b) :: i,j
                    DO i=1,SIZE(neq,1)
                        DO j=1,SIZE(neq,2)
                            WRITE(378,*)[i,j]*spaceGrid%length(1:2)/spaceGrid%n_nodes(1:2),neq(i,j,1,1)
                        END DO
                        WRITE(378,*)
                    END DO
                END BLOCK
            CLOSE(378)
        END IF
    END IF


    filename = 'output/rdf.out'
    CALL compute_rdf ( neq(:,:,:,1) , filename ) ! Get radial distribution functions

    filename = 'output/g-r-theta.dat'
    print*,"TODO CHANGE THE NAME OF OUTPUT. SEE PROCESS_OUTPUT.F90"
    call compute_gOfRandCosTheta

    CALL adhoc_corrections_to_gsolv
    CALL finalizeMinimizer


    CONTAINS

        SUBROUTINE print_cg_vect
            USE minimizer, ONLY: cg_vect
            if ( .not. allocated ( cg_vect ) ) then
                print *, 'cg_vect is not allocated in SUBROUTINE print_cg_vect in process_output.f90. STOP.'
                stop
            END IF
            OPEN (10, file = 'output/density.bin.out' , form = 'unformatted' )
                write ( 10 ) cg_vect
            CLOSE (10)
        END SUBROUTINE print_cg_vect

END SUBROUTINE process_output

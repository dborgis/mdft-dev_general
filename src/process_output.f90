! This subroutine computes every output asked by user.

SUBROUTINE process_output

    ! defines precision of reals and intergers
    USE precision_kinds,    ONLY: dp, i2b
    USE system,             ONLY: nfft1, nfft2, nfft3, nb_species
    USE input,              ONLY: input_line, verbose
    USE solute_geometry,    ONLY: soluteIsPlanar => isPlanar, soluteIsLinear => isLinear
    
    IMPLICIT NONE
    
    real(dp), dimension( nfft1 , nfft2 , nfft3 , nb_species ) :: neq ! equilibrium density(position)
    real(dp), dimension( nfft1 , nfft2 , nfft3 , nb_species ) :: Px , Py , Pz ! equilibrium polarization(position)
    character(50):: filename
    real(dp), allocatable , dimension ( : , : , : , : ) :: temparray
    integer(i2b):: i , j ! dummy

    call print_cg_vect ! print output/density.bin which has the content of cg_vect 
    call get_final_density ( neq ) ! Get the final density (position) from the last minimizer step.

    IF (verbose) THEN
        filename = 'output/density.cube'
        CALL write_to_cube_file (neq(:,:,:,1), filename) ! TODO for now only write for the first species
        DO i =1,SIZE(input_line)
            j =LEN('polarization ')
            IF ( input_line(i)(1:j)=='polarization ' .AND. input_line(i)(j+3:j+3)=='T' ) THEN
                CALL get_final_polarization (Px,Py,Pz)
                filename ='output/P.cube'
                ALLOCATE (temparray (nfft1,nfft2,nfft3,nb_species) )
                temparray =SQRT(Px**2+Py**2+Pz**2)
                CALL write_to_cube_file ( temparray(:,:,:,1), filename ) ! TODO for now only write for the first species
                filename='output/Px.cube' ; CALL write_to_cube_file(Px(:,:,:,1),filename)
                filename='output/Py.cube' ; CALL write_to_cube_file(Py(:,:,:,1),filename)
                filename='output/Pz.cube' ; CALL write_to_cube_file(Pz(:,:,:,1),filename)
                filename='output/z_Px.dat'; CALL compute_z_density(Px(:,:,:,1),filename)
                filename='output/z_Py.dat'; CALL compute_z_density(Py(:,:,:,1),filename)
                filename='output/z_Pz.dat'; CALL compute_z_density(Pz(:,:,:,1),filename)
                DEALLOCATE(temparray)
                EXIT
            END IF
        END DO
        
        ! If calculation is for hard sphere fluid in presence of a hard wall compute profile perp wall
        ! TODO: DONT HAVE TIME TO WRITE THE TEST TODAY
        filename = 'output/z_density.dat'
        CALL compute_z_density ( neq (:,:,:,1) , filename ) ! TODO for now only write for the first species
    
        IF (soluteIsLinear() ) THEN
            ! nothing for now
        END IF
        
        IF( soluteIsPlanar() ) THEN
            filename = 'output/planardensity.out'
            CALL compute_planar_density ( neq (:,:,:,1) , filename ) ! TODO for now only write for the first species
        END IF
    END IF


    filename = 'output/rdf.out'
    CALL compute_rdf ( neq , filename ) ! Get radial distribution functions


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

! This subroutine computes every output asked by user.

SUBROUTINE process_output

    ! defines precision of reals and intergers
    use precision_kinds,    ONLY: dp, i2b
    use system,             ONLY: solvent, thermocond
    use module_input,              ONLY: verbose, getinput
    use solute_geometry,    ONLY: soluteIsPlanar => isPlanar, soluteIsLinear => isLinear
    use constants,          ONLY: zerodp
    use hardspheres,        only: hs
    use module_grid, only: grid

    IMPLICIT NONE

    CHARACTER(50):: filename
    REAL(dp), ALLOCATABLE , DIMENSION (:,:,:,:) :: neq, Px, Py, Pz ! equilibrium density, ie rho(r), and Pi polarization(r)
    INTEGER(i2b) :: nfft1, nfft2, nfft3
    INTEGER(i2b) :: s

    nfft1 = grid%nx
    nfft2 = grid%ny
    nfft3 = grid%nz

    CALL print_cg_vect_new ! print output/density.bin that contains cg_vect_new

    allocate ( neq (nfft1,nfft2,nfft3,solvent(1)%nspec) ,SOURCE=zerodp)
    do s=1,size(solvent)
        call get_final_density (neq,s)
    end do

    DO s=1,solvent(1)%nspec
      write(*,'(A,F12.2)') "Solvent molecules in supercell", SUM(neq)*grid%dv *solvent(s)%n0
    END DO


    IF (verbose) THEN
        filename = 'output/density.cube'
        CALL write_to_cube_file (neq(:,:,:,1), filename) ! TODO for now only write for the first species
        IF (getinput%log('polarization').EQV. .TRUE.) THEN
            ALLOCATE ( Px (nfft1,nfft2,nfft3,solvent(1)%nspec) ,SOURCE=zerodp)
            ALLOCATE ( Py (nfft1,nfft2,nfft3,solvent(1)%nspec) ,SOURCE=zerodp)
            ALLOCATE ( Pz (nfft1,nfft2,nfft3,solvent(1)%nspec) ,SOURCE=zerodp)
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

        IF ( getinput%char('other_predefined_vext')=='vextdef0' ) THEN
            filename = 'output/molecular_density_in_xy_plane.out'
            PRINT*,"I am writing file ",filename
            OPEN(378,FILE=filename)
                BLOCK
                    INTEGER(i2b) :: i,j
                    DO i=1,SIZE(neq,1)
                        DO j=1,SIZE(neq,2)
                            WRITE(378,*)[i,j]*grid%length(1:2)/grid%n_nodes(1:2),neq(i,j,1,1)
                        END DO
                        WRITE(378,*)
                    END DO
                END BLOCK
            CLOSE(378)
        END IF
    END IF


    filename = 'output/rdf.out'
    CALL output_rdf ( neq(:,:,:,1) , filename ) ! Get radial distribution functions
    call output_gsitesite
    call output_gOfRandCosTheta

    CALL adhoc_corrections_to_gsolv


    write(*,'(A,F7.2,A)') "T       ", thermocond%T,    " K"
    write(*,'(A,F7.2,A)') "kT      ", thermocond%kbT,  " kJ/mol"
    write(*,'(A,F7.2,A)') "β=(kT)⁻¹", thermocond%beta, " (kJ/mol)⁻¹"
    if( allocated(hs) ) then
      block
        real(dp)::x
        x=hs(1)%pf
        write(*,'(A,F7.2)') "packing fraction η =",x
        write(*,'(A,F7.2)') "βP/n PY by pressure route        = (1+2η+3η²)           /(1-η)² =",(1+2*x+3*x**2)/(1-x)**2
        write(*,'(A,F7.2)') "βP/n PY by compressibility route = (1+ η+ η²)           /(1-η)³ =",(1+x+x**2)/(1-x)**3
        write(*,'(A,F7.2)') "βP/n CS                          = (1+ η+ η²-η³)        /(1-η)³ =",(1+x+x**2-x**3)/(1-x)**3
        write(*,'(A,F7.2)') "βP/n CSK                         = (1+ η+ η²-2(1+η)η³/3)/(1-η)³ =",&
          (1+x+x**2-(2./3.)*(1+x)*x**3)/((1.-x)**3)
      end block
    end if


    CONTAINS

        SUBROUTINE print_cg_vect_new
            use module_minimizer, ONLY: cg_vect_new
            if ( .not. allocated ( cg_vect_new ) ) then
                print *, 'cg_vect_new is not allocated in SUBROUTINE print_cg_vect_new in process_output.f90. STOP.'
                stop
            END IF
            OPEN (10, file = 'output/density.bin.out' , form = 'unformatted' )
                write ( 10 ) cg_vect_new
            CLOSE (10)
        END SUBROUTINE print_cg_vect_new

END SUBROUTINE process_output

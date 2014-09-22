!> read solute atomic positions, charge, and lennard jones values in solute.in
!! charge in electron units, sigma in Angstroms, epsilon in KJ/mol.
SUBROUTINE read_solute

    USE precision_kinds, ONLY: i2b,dp
    USE system,          ONLY: nb_solute_sites, solute, nb_species, spaceGrid
    USE input,           ONLY: input_line, input_log, input_dp
    USE periodic_table,  ONLY: init_periodic_table, ptable

    IMPLICIT NONE

    INTEGER(i2b) :: n,i,stat
    REAL(dp) :: lx, ly, lz
    
    lx= spaceGrid%length(1)
    ly= spacegrid%length(2)
    lz= spaceGrid%length(3)

    CALL init_periodic_table
    ! print *, ptable ( 1 ) % name
    ! open and test if input/solute.in is ok

    OPEN (5, FILE='input/solute.in', STATUS='old', IOSTAT=stat)
        IF (stat /= 0) THEN
            PRINT*,'solute.in cannot be opened ! => STOP !'
            STOP
        END IF

        READ (5,*) ! comment line
        READ (5,*) nb_solute_sites ! total number of atom sites of the solute AND the total number of different types of atoms
        ALLOCATE(solute%site(nb_solute_sites))

        READ (5,*)
        DO n = 1, SIZE(solute%site)
            READ(5,*) i, solute%site(n)%q, solute%site(n)%sig, solute%site(n)%eps, &
                    solute%site(n)%lambda1, solute%site(n)%lambda2, solute%site(n)%r, solute%site(n)%Z
        END DO
    CLOSE (5)
    solute%site%q = solute%site%q * input_dp('solute_charges_scale_factor')

    CALL translate_to_center_of_supercell_if_needed ! if user wants all the sites to be translated to the center of the box, ie by Lx/2, Ly/2, Lz/2
    CALL assure_coo_inside_cell ! check if cartesian coordinates read in input/solute.in are in the supercell
    CALL print_supercell_xsf ! Print periodic XSF file to be read by VMD or equivalent


    CONTAINS

    ! if user asks for it (tag 'translate_solute_to_center'), add Lx/2, Ly/2, Lz/2 to all solute coordinates
    SUBROUTINE translate_to_center_of_supercell_if_needed
        USE input, ONLY: input_log
        USE system, ONLY: spaceGrid
        IF (input_log( 'translate_solute_to_center' )) THEN
            solute%site%r(1) = solute%site%r(1) + spaceGrid%length(1)/2.0_dp
            solute%site%r(2) = solute%site%r(2) + spaceGrid%length(2)/2.0_dp
            solute%site%r(3) = solute%site%r(3) + spaceGrid%length(3)/2.0_dp
        END IF
    END SUBROUTINE translate_to_center_of_supercell_if_needed



    SUBROUTINE assure_coo_inside_cell
        USE SYSTEM, ONLY: spaceGrid
        INTEGER(i2b) :: i 
        ! check if some positions are out of the supercell
        !j is a test tag. We loop over this test until every atom is in the box.
        ! This allows for instance, if a site is two boxes too far to still be ok.
        DO CONCURRENT ( i=1:SIZE(solute%site) )
            solute%site(i)%r(1) = MODULO ( solute%site(i)%r(1) , spaceGrid%length(1) )
            solute%site(i)%r(2) = MODULO ( solute%site(i)%r(2) , spaceGrid%length(2) )
            solute%site(i)%r(3) = MODULO ( solute%site(i)%r(3) , spaceGrid%length(3) )
        END DO
    END SUBROUTINE assure_coo_inside_cell


END SUBROUTINE read_solute

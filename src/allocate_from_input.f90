! Here are allocated variables declared in modules
SUBROUTINE allocate_from_input

    USE precision_kinds     ,ONLY: i2b , dp
    USE input               ,ONLY: input_line, input_int, input_dp, input_log, verbose
    USE system              ,ONLY: thermoCond, n_0, rho_0, n_0_multispec, rho_0_multispec, nb_species, mole_fraction, spaceGrid
    USE constants           ,ONLY: fourpi, boltz, navo, twopi
    USE quadrature          ,ONLY: molRotSymOrder

    IMPLICIT NONE
    INTEGER(i2b):: i, j, s

    verbose = input_log('verbose')
    
    molRotSymOrder = input_int('molRotSymOrder') !Get the order of the main symmetry axis of the solvent
    IF (molRotSymOrder < 1) THEN
        PRINT*,'order of main symetric axe cannot be less than 1. molRotSymOrder is declared as ',molRotSymOrder
        STOP 'CRITICAL STOP. CHANGE molRotSymOrder IN INPUT'
    END IF
    
    spaceGrid%n_nodes = [ input_int('nfft1'), input_int('nfft2'), input_int('nfft3') ] ! number of grid nodes in each direction
    IF (ANY( spaceGrid%n_nodes  <= 0) ) THEN
        PRINT*,'The space is divided into grid nodes. For each direction, you ask', spaceGrid%n_nodes,'node.'
        STOP 'This is unphysical.'
    END IF

    spaceGrid%length = [ input_dp('Lx'), input_dp('Ly'), input_dp('Lz') ]
    IF (ANY( spaceGrid%length  <= 0._dp ) ) THEN
        PRINT*,'The supercell cannot have negative length.'
        PRINT*,'Here are your Lx, Ly and Lz as defined in input/dft.in :',spaceGrid%length
        STOP 'CRITICAL STOP BECAUSE OF NON-PHYSICAL INPUT'
    END IF
    spaceGrid%dl = spaceGrid%length/REAL(spaceGrid%n_nodes,dp)
    spaceGrid%dv = product(spaceGrid%dl)

    thermoCond%T = input_dp('temperature') ! look for temperature in input
    IF (thermoCond%T <= 0 ) THEN
        PRINT*,'CRITICAL STOP. NEGATIVE TEMPERATURE IN INPUT FILE tag temperature :',thermoCond%T
        STOP
    END IF
    thermoCond%kbT = Boltz * Navo * thermoCond%T * 1.0e-3_dp
    thermoCond%beta = 1.0_dp / thermocond%kbT
    
    nb_species = input_int('nb_implicit_species') ! get the number of implicit solvant species
    IF (nb_species < 1 ) THEN
        PRINT*,'Solvent is without species. This is not what you want. Check nb_implicit_species in dft.in :', nb_species
        STOP
    END IF
    
    ! look for bulk density of the reference solvent fluid. for instance 0.0332891 for H2O and 0.0289 for Stockmayer  
    ALLOCATE ( n_0_multispec ( nb_species ) )
    do i = 1 , size ( input_line )
        if ( input_line ( i ) ( 1 : len ( 'ref_bulk_density' ) ) == 'ref_bulk_density' ) THEN
            do s = 1 , nb_species
                read ( input_line ( i + s ) , * ) n_0_multispec ( s )
            END DO
            exit
        END IF
    END DO
    IF (ANY( n_0_multispec <= 0._dp) ) THEN
        PRINT*,'you ask for several species as solvent, but some have negative bulk density'
        do i=1, size(n_0_multispec)
            PRINT*,'species',i,'bulk density',n_0_multispec(i)
        END DO
        STOP 'UNPHYSICAL INPUT CHECK ref_bulk_density'
    END IF
    ALLOCATE ( rho_0_multispec ( nb_species ) )
    rho_0_multispec = molRotSymOrder*n_0_multispec / (2.0_dp*twopi**2) ! for single specie compatibility while not fully complete :
    n_0 = n_0_multispec(1) ! for single specie compatibility while not fully complete : 
    rho_0 = rho_0_multispec(1) ! rho is the density per angle

    if (ALLOCATEd(mole_fraction)) THEN
        write(*,*)'something is not under control with respect to mole fraction reading and definition in ALLOCATE_from_input.f90'
        write(*,*)'this was done, in a previous version of the program, in compute_hard_sphere_parameters.f90'
        STOP
    ELSE
        ALLOCATE ( mole_fraction ( nb_species ) ) ! molar fraction of each species.
        call read_mole_fractions ( nb_species , mole_fraction )
    END IF


    
    
    contains  
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to read the mole fractions in dft.in.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! This SUBROUTINE open the array input_line which contains every line of input/dft.in
    ! It then reads every line of input_line and looks for the tag "mole_fractions"
    ! Then, it reads, one line after the other, the mole fractions of every constituant.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE read_mole_fractions ( nb_species , mole_fraction )
        USE precision_kinds,only : dp , i2b
        use input,only : input_line
        IMPLICIT NONE
        integer(i2b), intent(in) :: nb_species
        real(dp), dimension ( nb_species ) , intent ( inout ) :: mole_fraction
        integer(i2b):: i, j, s
        do i = 1 , size ( input_line )
        j = len ( 'mole_fractions' )
        if ( input_line (i) (1:j) == 'mole_fractions' ) THEN
            do s = 1 , nb_species
                read ( input_line ( i + s ) , * ) mole_fraction ( s )
            END DO
            exit ! loop over i
        END IF
        END DO
    
        ! check error in mole fraction : sum of all mole fractions should be equal to 1
        call check_error_in_mole_fraction ( mole_fraction )
    END SUBROUTINE read_mole_fractions


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! check_error_in_mole_fraction checks if there are errors in mole fractions
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! It checks that the sum of all mole fractions is 1, and that every mole fractions are between 0 and 1.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE check_error_in_mole_fraction ( mole_fraction )
        USE precision_kinds,only : i2b , dp
        use system,only : nb_species
        IMPLICIT NONE
        real(dp), dimension ( nb_species ) , intent(in) :: mole_fraction
        integer(i2b):: s
    
        ! sum of all mole fractions should be equal to 1
        if ( sum ( mole_fraction ) /= 1.0_dp ) THEN
            write (*,*) 'Critial error. Sum of all mole fraction should be equal to one.'
            write (*,*) 'here are the number of the species and its associated mole fraction'
            do s = 1 , nb_species
                PRINT*, s , mole_fraction(s)
            END DO
            write (*,*) 'STOP'
            STOP
        END IF
        ! a mole fraction should be between 0 and 1
        do s = 1 , nb_species
            if ( mole_fraction ( s ) < 0.0_dp .or. mole_fraction ( s ) > 1.0_dp ) THEN
                write (*,*) 'Critical errror in ALLOCATE_from_input.f90. Mole fractions should be between 0 and 1'
                write (*,*) 'here are the number of the species and its associated mole fraction'
                write (*,*) 'species number ' , s , ' has mole fraction ' , mole_fraction ( s )
                write (*,*) 'STOP'
                STOP
            END IF
        END DO
    END SUBROUTINE check_error_in_mole_fraction



END SUBROUTINE ALLOCATE_from_input

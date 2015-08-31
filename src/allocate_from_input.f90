! Here are allocated variables declared in modules
SUBROUTINE allocate_from_input

    USE precision_kinds     ,ONLY: i2b , dp
    USE input               ,ONLY: input_line, input_int, input_dp, input_log, verbose, input_dp3, input_int3
    USE system              ,ONLY: thermoCond, nb_species, mole_fraction, gr=>spacegrid, solvent, solute, nb_solute_sites
    USE constants           ,ONLY: eightpiSQ, boltz, navo
    USE quadrature          ,ONLY: molRotSymOrder

    IMPLICIT NONE
    INTEGER(i2b):: i, s

    verbose = input_log('verbose',defaultvalue=.false.)

    molRotSymOrder = input_int('molRotSymOrder', defaultvalue=1) !Get the order of the main symmetry axis of the solvent

    if (molRotSymOrder < 1) THEN
      print*, 'order of main symetric axe cannot be less than 1. molRotSymOrder is declared as ',molRotSymOrder
      stop    'CRITICAL STOP. CHANGE molRotSymOrder IN INPUT'
    else if (molRotSymOrder > 2) then
      print*, "I am surprised your molrotsymorder >2. Certainly a problem somewhere."
      stop
    end if

    gr%length = input_dp3( "lxlylz" , defaultvalue= gr%length )
    if (ANY( gr%length  <= 0._dp ) ) THEN
        PRINT*,'The supercell cannot have negative length.'
        PRINT*,'Here are your Lx, Ly and Lz as defined in input/dft.in :',gr%length
        STOP
    end if
    gr%l = gr%length
    gr%lx = gr%length(1)
    gr%ly = gr%length(2)
    gr%lz = gr%length(3)

    gr%n_nodes = input_int3( "nxnynz" , defaultvalue= nint(gr%length/0.3_dp) )
    if ( any(gr%n_nodes <= 0) ) then
        print*, 'The space is divided into grid nodes. For each direction, you ask', gr%n_nodes,'node.'
        error stop
    end if
    gr%n = gr%n_nodes
    gr%nx = gr%n_nodes(1)
    gr%ny = gr%n_nodes(2)
    gr%nz = gr%n_nodes(3)

    gr%dl = gr%length / real(gr%n_nodes,dp)
    gr%dx = gr%dl(1)
    gr%dy = gr%dl(2)
    gr%dz = gr%dl(3)

    gr%v = product(gr%length)
    gr%dv = product(gr%dl)

    ! We now have a full description of the space grid
    print*,
    print*, "====GRID==========="
    print*, "Box Length :", gr%length
    print*, "nodes      :", gr%n_nodes
    print*, "dx, dy, dz :", gr%dl
    print*, "====/GRID==========="
    print*,


    thermoCond%T = input_dp('temperature', defaultvalue=300._dp) ! look for temperature in input
    IF (thermoCond%T <= 0 ) THEN
        PRINT*,'CRITICAL STOP. NEGATIVE TEMPERATURE IN INPUT FILE tag temperature :',thermoCond%T
        STOP
    end if
    thermoCond%kbT = Boltz * Navo * thermoCond%T * 1.0e-3_dp
    thermoCond%beta = 1.0_dp / thermocond%kbT



    nb_species = input_int('nb_implicit_species',defaultvalue=1) ! get the number of implicit solvant species
    if (nb_species < 1) then
        print*,"STOP. nb_species =",nb_species,". The number of solvent species must be >0."
        stop
    end if
    allocate (solvent(nb_species), stat=i)
    if (i /= 0) stop "Problem during allocation of type solvent in allocate from input."

    ! look for bulk density of the reference solvent fluid. for instance 0.0332891 for H2O and 0.0289 for Stockmayer
    do i = 1, size(input_line)
        if (input_line(i)(1:len('ref_bulk_density')) == 'ref_bulk_density') then
            do s = 1, nb_species
                read (input_line(i+s),*) solvent(s)%n0
            end do
            exit
        end if
    end do

    if (any (solvent%n0 <= 0._dp) ) then
        print *,"You ask for negative densities!"
        do s =1, nb_species
            print *,"For species",s,"you want density (molecule/Ang^3):",solvent(s)%n0
        end do
        stop
    end if
    solvent%rho0 = solvent%n0 / (eightpiSQ/molrotsymorder)

    if (nb_species > 1) stop "molRotSymOrder must be solvent specific. See github issu #60"

    if (ALLOCATEd(mole_fraction)) THEN
        write(*,*)'something is not under control with respect to mole fraction reading and definition in ALLOCATE_from_input.f90'
        write(*,*)'this was done, in a previous version of the program, in compute_hard_sphere_parameters.f90'
        STOP
    ELSE
        ALLOCATE ( mole_fraction ( nb_species ) ) ! molar fraction of each species.
        call read_mole_fractions ( nb_species , mole_fraction )
    end if


    CALL mv_solute_to_center ! if user wants all the sites to be translated to the center of the box, ie by Lx/2, Ly/2, Lz/2
    CALL assert_solute_inside ! check if cartesian coordinates read in input/solute.in are in the supercell
    CALL print_supercell_xsf ! Print periodic XSF file to be read by VMD or equivalent

    CONTAINS

    ! if user asks for it (tag 'translate_solute_to_center'), add Lx/2, Ly/2, Lz/2 to all solute coordinates
    subroutine mv_solute_to_center
        use input  ,only: input_log
        use system ,only: gr=>spacegrid, solute
        implicit none
        if( input_log( 'translate_solute_to_center', defaultvalue=.true. )) then
            solute%site%r(1) = solute%site%r(1) + gr%length(1)/2.0_dp
            solute%site%r(2) = solute%site%r(2) + gr%length(2)/2.0_dp
            solute%site%r(3) = solute%site%r(3) + gr%length(3)/2.0_dp
        end if
    end subroutine mv_solute_to_center

    subroutine assert_solute_inside
        use system, only: gr=>spacegrid, solute
        implicit none
        integer :: i
        ! check if some positions are out of the supercell
        !j is a test tag. We loop over this test until every atom is in the box.
        ! This allows for instance, if a site is two boxes too far to still be ok.
        do concurrent( i=1:nb_solute_sites )
            solute%site(i)%r(1) = MODULO ( solute%site(i)%r(1) , gr%length(1) )
            solute%site(i)%r(2) = MODULO ( solute%site(i)%r(2) , gr%length(2) )
            solute%site(i)%r(3) = MODULO ( solute%site(i)%r(3) , gr%length(3) )
        end do
    end subroutine assert_solute_inside

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
        select case (nb_species)
        case (1)
          mole_fraction = 1._dp
        case default
          do i = 1, size( input_line )
            j = len ( 'mole_fractions' )
            if ( input_line(i)(1:j) == 'mole_fractions' ) then
              do s = 1, nb_species
                  read( input_line(i+s),*) mole_fraction(s)
              end do
              exit ! loop over i
            end if
          end do
        end select
        call check_error_in_mole_fraction ( mole_fraction ) ! check mole fractions are physical
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
              print*, s , mole_fraction(s)
            end do
            write (*,*) 'STOP'
            stop
        end if
        ! a mole fraction should be between 0 and 1
        do s = 1 , nb_species
          if ( mole_fraction ( s ) < 0.0_dp .or. mole_fraction ( s ) > 1.0_dp ) THEN
            write (*,*) 'Critical errror in ALLOCATE_from_input.f90. Mole fractions should be between 0 and 1'
            write (*,*) 'here are the number of the species and its associated mole fraction'
            write (*,*) 'species number ' , s , ' has mole fraction ' , mole_fraction ( s )
            write (*,*) 'STOP'
            stop
          end if
        end do
    END SUBROUTINE check_error_in_mole_fraction

END SUBROUTINE ALLOCATE_from_input

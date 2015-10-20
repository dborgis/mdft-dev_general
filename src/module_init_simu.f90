module module_init_simu

    use quadrature, only: prepare_quadrature => init
    use fft, only: init_fft => init
    use dcf, only: init_dcf => init
    use module_input, only: getinput
    use hardspheres, only: compute_hard_spheres_parameters
    use module_density, only: init_density
    use module_solvent, only: init_solvent => read_solvent
    use module_solute, only: read_solute

    implicit none

    private
    public :: init_simu

contains

    !
    !   Init the simulation
    !
    subroutine init_simu

        implicit none

        call print_header ! package name, day & time
        call read_solute ! AWFUL TODO trick to default lx
        call allocate_from_input ! TODO this should be removed and done only when needed ! read dft.in, solute.in and solvent.in
        call prepare_quadrature
        call init_fft
        call init_solvent ! Read solvent atomic positions, charge and Lennard-Jones param
        call init_dcf
        if (getinput%log('hard_sphere_fluid', defaultvalue=.false.)) call compute_hard_spheres_parameters ! If calculation based on Fundamental Measure Theory read FMT parameters and compute weight functions etc
        call init_external_potential
        call init_density
    end subroutine init_simu

    !
    !   Print a header to stdout to welcome the user
    !.
    subroutine print_header
        character(8)  :: date
        character(10) :: time
        call date_and_time ( DATE=date,TIME=time)
        print*,' ******************************************************************'
        print*,' ******************************************************************'
        print*,' **             *************************   ',date(1:4),'/',date(5:6),'/',date(7:8),'   **********'
        print*,' **   M D F T   *************************                **********'
        print*,' **             *************************    ',time(1:2),'h',time(3:4),'m',time(5:6),'    **********'
        print*,' ******************************************************************'
        print*,' ******************************************************************'
    end subroutine print_header

    ! Here are allocated variables declared in modules
    SUBROUTINE allocate_from_input

        use precision_kinds     ,ONLY: i2b , dp
        use module_input        ,ONLY: input_line, getinput, verbose
        use system              ,ONLY: thermoCond, mole_fraction
        use module_solvent, only: solvent
        use constants           ,ONLY: eightpiSQ, boltz, navo
        use module_grid, only: grid, init_grid

        IMPLICIT NONE
        integer:: i, s

        verbose = getinput%log( "verbose", defaultvalue=.false.)

        print*,
        print*, "[Conditions]====="
        thermoCond%T = getinput%dp('temperature', defaultvalue=300._dp) ! look for temperature in input
        IF (thermoCond%T <= 0 ) THEN
            PRINT*,'CRITICAL STOP. NEGATIVE TEMPERATURE IN INPUT FILE tag temperature :',thermoCond%T
            STOP
        end if
        thermoCond%kbT = Boltz * Navo * thermoCond%T * 1.0e-3_dp
        thermoCond%beta = 1.0_dp / thermocond%kbT
        print*, "Temperature (K)  :", thermocond%t
        print*, "kT               :", thermocond%kbt
        print*, "\beta = 1/(kT)   :", thermocond%beta
        print*, "[/Conditions]====="
        print*,

        call init_grid


        s = getinput%int('nb_implicit_species', defaultvalue=1, assert=">0") ! get the number of implicit solvant species
        allocate( solvent(s) )
        solvent(:)%nspec = s
        solvent(1)%nspec = s



        ! look for bulk density of the reference solvent fluid. for instance 0.0332891 for H2O and 0.0289 for Stockmayer
        do i = 1, size(input_line)
            if (input_line(i)(1:len('ref_bulk_density')) == 'ref_bulk_density') then
                do s = 1, solvent(1)%nspec
                    read (input_line(i+s),*) solvent(s)%n0
                end do
                exit
            end if
        end do

        if (any (solvent%n0 <= 0._dp) ) then
            print *,"You ask for negative densities!"
            do s =1, solvent(1)%nspec
                print *,"For species",s,"you want density (molecule/Ang^3):",solvent(s)%n0
            end do
            stop
        end if
        solvent%rho0 = solvent%n0 / (eightpiSQ/grid%molrotsymorder)


        if (ALLOCATEd(mole_fraction)) THEN
            write(*,*)'something is not under control with respect to mole fraction reading and definition in ALLOCATE_from_input.f90'
            write(*,*)'this was done, in a previous version of the program, in compute_hard_sphere_parameters.f90'
            STOP
        ELSE
            ALLOCATE ( mole_fraction(solvent(1)%nspec)) ! molar fraction of each species.
            call read_mole_fractions(solvent(1)%nspec, mole_fraction )
        end if


        CALL mv_solute_to_center ! if user wants all the sites to be translated to the center of the box, ie by Lx/2, Ly/2, Lz/2
        CALL assert_solute_inside ! check if cartesian coordinates read in input/solute.in are in the supercell
        CALL print_supercell_xsf ! Print periodic XSF file to be read by VMD or equivalent

    END SUBROUTINE ALLOCATE_from_input




    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to read the mole fractions in dft.in.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! This SUBROUTINE open the array input_line which contains every line of input/dft.in
    ! It then reads every line of input_line and looks for the tag "mole_fractions"
    ! Then, it reads, one line after the other, the mole fractions of every constituant.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE read_mole_fractions ( nspec , mole_fraction )
        use precision_kinds,only : dp , i2b
        IMPLICIT NONE
        integer(i2b), intent(in) :: nspec
        real(dp), dimension ( nspec ) , intent ( inout ) :: mole_fraction
        integer(i2b):: i, j, s
        select case (nspec)
        case (1)
            mole_fraction = 1._dp
        case default
            do i = 1, size( input_line )
                j = len ( 'mole_fractions' )
                if ( input_line(i)(1:j) == 'mole_fractions' ) then
                    do s = 1, nspec
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
        use precision_kinds, only: i2b , dp
        use module_solvent, only: solvent
        IMPLICIT NONE
        real(dp), dimension ( solvent(1)%nspec ) , intent(in) :: mole_fraction
        integer(i2b):: s

        ! sum of all mole fractions should be equal to 1
        if ( sum ( mole_fraction ) /= 1.0_dp ) THEN
            write (*,*) 'Critial error. Sum of all mole fraction should be equal to one.'
            write (*,*) 'here are the number of the species and its associated mole fraction'
            do s = 1 , solvent(1)%nspec
                print*, s , mole_fraction(s)
            end do
            write (*,*) 'STOP'
            stop
        end if
        ! a mole fraction should be between 0 and 1
        do s = 1 , solvent(1)%nspec
            if ( mole_fraction ( s ) < 0.0_dp .or. mole_fraction ( s ) > 1.0_dp ) THEN
                write (*,*) 'Critical errror in ALLOCATE_from_input.f90. Mole fractions should be between 0 and 1'
                write (*,*) 'here are the number of the species and its associated mole fraction'
                write (*,*) 'species number ' , s , ' has mole fraction ' , mole_fraction ( s )
                write (*,*) 'STOP'
                stop
            end if
        end do
    END SUBROUTINE check_error_in_mole_fraction


end module module_init_simu

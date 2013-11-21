! Here are allocated variables declared in modules
subroutine allocate_from_input

    use precision_kinds , only : i2b , dp
    use input , only : input_line, input_int, input_dp, input_log
    USE system, ONLY: n_0, rho_0, temp, beta, kbT, Rc, n_0_multispec, rho_0_multispec, nb_species, mole_fraction, spaceGrid
    use constants , only : fourpi , boltz, navo , twopi
    use quadrature , only : sym_order

    implicit none

    integer(i2b):: i , j ! dummy
    integer(i2b):: species ! dummy between 1 and nb_species
    
    sym_order = input_int('sym_order') !Get the order of the main symmetry axis of the solvent
    if( sym_order <= 1 ) then
        print*,'order of main symetric axe cannot be less than 1. sym_order is declared as ',sym_order
        stop 'CRITICAL STOP. CHANGE sym_order IN INPUT'
    end if
    
    spaceGrid%n_nodes = [ input_int('nfft1'), input_int('nfft2'), input_int('nfft3') ] ! number of grid nodes in each direction
    if( any( spaceGrid%n_nodes  <= 0) ) then
        print*,'The space is divided into grid nodes. The number of nodes must be greater than 1.'
        print*,'You have nfft1, nfft2, nfft3 nodes per direction : ',spaceGrid%n_nodes
        stop 'CRITICAL STOP BECAUSE OF NON-PHYSICAL INPUT.'
    end if

    spaceGrid%length = [ input_dp('Lx'), input_dp('Ly'), input_dp('Lz') ]
    if( any( spaceGrid%length  <= 0._dp ) ) then
        print*,'The supercell cannot have negative length.'
        print*,'Here are your Lx, Ly and Lz as defined in input/dft.in :',spaceGrid%length
        stop 'CRITICAL STOP BECAUSE OF NON-PHYSICAL INPUT'
    end if
    spaceGrid%dl = spaceGrid%length/REAL(spaceGrid%n_nodes,dp)
    spaceGrid%dv = product(spaceGrid%dl)

    temp = input_dp('temperature') ! look for temperature in input
    if( temp <= 0 ) then
        print*,'CRITICAL STOP. NEGATIVE TEMPERATURE IN INPUT FILE tag temperature :',temp
        stop
    end if
    kBT = Boltz * Navo * temp * 1.0e-3_dp
    beta = 1.0_dp / kBT
    
    nb_species=input_int('nb_implicit_species') ! get the number of implicit solvant species
    if( nb_species < 1 ) then
        print*,'Solvent is without species. This is not what you want. Check nb_implicit_species in dft.in :',nb_species
        stop
    end if
    
    ! look for bulk density of the reference solvent fluid. for instance 0.0332891 for H2O and 0.0289 for Stockmayer  
    allocate ( n_0_multispec ( nb_species ) )
    do i = 1 , size ( input_line )
        if ( input_line ( i ) ( 1 : len ( 'ref_bulk_density' ) ) == 'ref_bulk_density' ) then
            do species = 1 , nb_species
                read ( input_line ( i + species ) , * ) n_0_multispec ( species )
            end do
            exit
        end if
    end do
    if( any( n_0_multispec <= 0._dp) ) then
        print*,'you ask for several species as solvent, but some have negative bulk density'
        do i=1, size(n_0_multispec)
            print*,'species',i,'bulk density',n_0_multispec(i)
        end do
        STOP 'UNPHYSICAL INPUT CHECK ref_bulk_density'
    end if
    allocate ( rho_0_multispec ( nb_species ) )
    rho_0_multispec = sym_order*n_0_multispec / (2.0_dp*twopi**2) ! for single specie compatibility while not fully complete :
    n_0 = n_0_multispec ( 1 ) ! for single specie compatibility while not fully complete : 
    rho_0 = sym_order*n_0 / ( twopi**2*2.0_dp) ! rho is the density per angle

    if (allocated(mole_fraction)) then
        write(*,*)'something is not under control with respect to mole fraction reading and definition in allocate_from_input.f90'
        write(*,*)'this was done, in a previous version of the program, in compute_hard_sphere_parameters.f90'
        stop
    else
        allocate ( mole_fraction ( nb_species ) ) ! molar fraction of each species.
        call read_mole_fractions ( nb_species , mole_fraction )
    end if

    ! look for Rc
    Rc=input_dp('Rc')
    
    
    
    contains  
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to read the mole fractions in dft.in.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! This subroutine open the array input_line which contains every line of input/dft.in
    ! It then reads every line of input_line and looks for the tag "mole_fractions"
    ! Then, it reads, one line after the other, the mole fractions of every constituant.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine read_mole_fractions ( nb_species , mole_fraction )
        use precision_kinds , only : dp , i2b
        use input , only : input_line
        implicit none
        integer(i2b), intent(in) :: nb_species
        real(dp), dimension ( nb_species ) , intent ( inout ) :: mole_fraction
        integer(i2b):: i , j , species
        do i = 1 , size ( input_line )
        j = len ( 'mole_fractions' )
        if ( input_line (i) (1:j) == 'mole_fractions' ) then
            do species = 1 , nb_species
                read ( input_line ( i + species ) , * ) mole_fraction ( species )
            end do
            exit ! loop over i
        end if
        end do
    
        ! check error in mole fraction : sum of all mole fractions should be equal to 1
        call check_error_in_mole_fraction ( mole_fraction )
    end subroutine read_mole_fractions


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! check_error_in_mole_fraction checks if there are errors in mole fractions
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! It checks that the sum of all mole fractions is 1, and that every mole fractions are between 0 and 1.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine check_error_in_mole_fraction ( mole_fraction )
        use precision_kinds , only : i2b , dp
        use system , only : nb_species
        implicit none
        real(dp), dimension ( nb_species ) , intent(in) :: mole_fraction
        integer(i2b):: species
    
        ! sum of all mole fractions should be equal to 1
        if ( sum ( mole_fraction ) /= 1.0_dp ) then
            write (*,*) 'Critial error. Sum of all mole fraction should be equal to one.'
            write (*,*) 'here are the number of the species and its associated mole fraction'
            do species = 1 , nb_species
                write (*,*) species , mole_fraction ( species )
            end do
            write (*,*) 'stop'
            stop
        end if
        ! a mole fraction should be between 0 and 1
        do species = 1 , nb_species
            if ( mole_fraction ( species ) < 0.0_dp .or. mole_fraction ( species ) > 1.0_dp ) then
                write (*,*) 'Critical errror in allocate_from_input.f90. Mole fractions should be between 0 and 1'
                write (*,*) 'here are the number of the species and its associated mole fraction'
                write (*,*) 'species number ' , species , ' has mole fraction ' , mole_fraction ( species )
                write (*,*) 'stop'
                stop
            end if
        end do
    end subroutine check_error_in_mole_fraction



end subroutine allocate_from_input

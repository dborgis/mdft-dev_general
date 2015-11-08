module module_solvent

    use precision_kinds, only: dp
    use system, only: site_type
    use module_input, only: getinput
    use module_grid, only: grid

    implicit none
    private

    type :: do_type
        logical ::  id_and_ext = .false.,&
        exc_cs = .false.,&
        exc_cdeltacd = .false.,&
        exc_cproj = .false.,&
        exc_fmt = .false.,&
        exc_wca = .false.,&
        exc_3b = .false.,&
        exc_multipolar_without_coupling_to_density = .false.,&
        exc_multipolar_with_coupling_to_density = .false.,&
        exc_hydro = .false.,&
        exc_nn_cs_plus_nbar = .false.
    end type

    type :: correlationfunction_type
        character(150) :: filename
        real(dp), allocatable :: x(:), y(:)
        logical :: isok = .false.
    end type

    type :: solvent_type
        logical :: is_initiated=.false.
        character(130) :: name
        integer :: molrotsymorder
        integer :: nsite ! number of site of the solvent molecule
        integer :: nspec ! number of solvent species
        real(dp) :: monopole, dipole(3), quadrupole(3,3), octupole(3,3,3), hexadecapole(3,3,3,3)
        real(dp) :: diameter ! hard sphere diameter, for instance
        real(dp), allocatable :: density(:,:,:,:) ! ix, iy, iz, io
        type (site_type), allocatable :: site(:)
        real(dp)              :: n0        ! number density of the homogeneous reference fluid in molecules per Angstrom^3, e.g., 0.033291 molecule.A**-3 for water
        real(dp)              :: rho0      ! number density per orientation of the homogeneous reference fluid in molecules per Angstrom^3 per orientation
        complex(dp), allocatable :: sigma_k(:,:,:,:) ! charge factor
        complex(dp), allocatable :: molec_polar_k(:,:,:,:,:) ! molecule polarization factor
        real(dp), allocatable :: vext(:,:,:,:), vextq(:,:,:,:)
        real(dp) :: vext_threeshold = 100._dp!36.04_dp ! 36.something is the maximum value of v so that exp(-beta.v) does not return underflow at 300 K
        type(do_type) :: do
        real(dp) :: mole_fraction = 1._dp
        type(correlationfunction_type) :: cs
        type(correlationfunction_type) :: cdelta
        type(correlationfunction_type) :: cd
    contains
        procedure, nopass :: init => read_solvent
    end type
    type (solvent_type), allocatable :: solvent(:)

    public :: solvent, read_solvent, print_solvent_not_allocated

contains

    subroutine print_solvent_not_allocated (message)
        implicit none
        character(*), intent(in) :: message
        print*, "solvent% is not allocated where it should"
        print*, message
        error stop
    end subroutine

    subroutine functional_decision_tree
        use module_input, only: getinput
        implicit none
        integer :: s
        !
        !   Later, we can imagine having different "functional trees" for each solvent.
        !   That makes sense for instance in water + ions, with HNC for water and MSA for ions.
        !   For now all species have the same functional.
        !
        do s=1,solvent(1)%nspec
            solvent(s)%do%id_and_ext = .true.
            solvent(s)%do%exc_fmt = getinput%log ('hard_sphere_fluid', defaultvalue=.false.)
            solvent(s)%do%exc_wca = getinput%log ('wca', defaultvalue=.false.)
            solvent(s)%do%exc_3b = getinput%log ('threebody', defaultvalue=.false. )

            select case (grid%mmax)
            case (0)
                solvent(s)%do%exc_cs=.true.
                solvent(s)%do%exc_cproj=.true.
            case (1:5)
                solvent(s)%do%exc_cs=.true.
                solvent(s)%do%exc_cdeltacd=.true.
                solvent(s)%do%exc_cproj=.true.
            case default
                print*, "see module_solvent > functional decision tree"
                print*, "mmax is trop grand"
                error stop
            end select
            !
            ! select case (getinput%char("polarization", defaultvalue="no"))
            ! case("no","none")
            ! case("cdeltacd")
            !     solvent(s)%do%exc_cdeltacd = .true.
            ! case("multipolar_without_coupling_to_density")
            !     solvent(s)%do%exc_multipolar_without_coupling_to_density = .true.
            ! case("multipolar_with_coupling_to_density")
            !     solvent(s)%do%exc_multipolar_with_coupling_to_density = .true.
            ! case default
            !     print*, "The tag 'polarization' in input reads ", getinput%char("polarization", defaultvalue="no")&
            !     ,". This is not correct"
            !     stop "in energy_and_gradient"
            ! end select
            !
            !
            ! if (getinput%log('readDensityDensityCorrelationFunction', defaultvalue=.true.)) THEN
            !     if (getinput%log('hydrophobicity', defaultvalue=.false.)) THEN
            !         select case (getinput%char('treatment_of_hydro'))
            !         case ('C')
            !             solvent(s)%do%exc_nn_cs_plus_nbar = .true.
            !         case ('VdW')
            !             if (getinput%log('bridge_hard_sphere', defaultvalue=.false.) ) THEN
            !                 print*, 'You are using HSB and VdW so you are withdrawing twice the HS second order term'
            !                 stop
            !             end if
            !             solvent(s)%do%exc_hydro = .true.
            !         case default
            !             stop "Hydrophobicity TRUE can only be associated to treatment_of_hydro == C or VdW"
            !         end select
            !     else
            !         if( .not. getinput%log('include_nc_coupling', defaultvalue=.false.) ) then
            !             solvent(s)%do%exc_cs = .true.
            !         end if
            !     end if
            ! end if

        end do
        !
        ! !
        ! ! bridge calculation: F(FMT)-F(c2hs)+F(c2H2O)
        ! !
        ! IF (getinput%log('bridge_hard_sphere', defaultvalue=.false.) .AND. .NOT. getinput%log('hard_sphere_fluid',defaultvalue=.false.)) THEN
        !     STOP 'bridge_hard_sphere and hard_sphere_fluid should be both turned TRUE for a calculation with bridge'
        ! END IF

    end subroutine functional_decision_tree


    !===================================================================================================================================
    subroutine read_solvent
        !===================================================================================================================================
        ! Read solvent atomic positions, charge, and lennard jones values in solvent.in
        ! charge in electron units, sigma in Angstroms, epsilon in KJ/mol.
        use mathematica, only: chop
        use module_input, only: getinput
        implicit none
        integer :: n, ios, i, j, k, l, s
        character(180) :: polarization

        if (allocated(solvent)) then
            print*, "bug dans read_solvent. Le 21 octobre 2015, j'ai transféré l'allocation de allocate_from_input a read_solvent"
            print*, "ca a bcp plus de sens de le mettre dans read_solvent que dans le fourre-tout"
            print*, "cependant il est alloue alors qu'il ne devrait pas l'etre"
            stop "dans module_solvent/read_solvent"
        end if


        s = getinput%int('nb_implicit_species', defaultvalue=1, assert=">0") ! get the number of implicit solvant species
        allocate( solvent(s) )
        solvent(:)%nspec = s


        OPEN(5, FILE= 'input/solvent.in', STATUS= 'old', IOSTAT= ios )! open input/solvent.in and check if it is readable
        IF ( ios/=0 ) STOP 'ERROR: solvent.in can not be opened.'
        READ (5,*) solvent(1)%name
        READ (5,*) solvent(1)%nsite!, solvent(1)%molrotsymorder
        i = solvent(1)%nsite
        allocate (solvent(1)%site(i), stat=ios)
        if (ios /= 0) stop "ERROR: wrong allocate of solvent%site in read_solvent.f90"
        READ(5,*) ! comment line
        DO n = 1 , size(solvent(1)%site)
            READ(5,*) i, solvent(1)%site(n)%q, solvent(1)%site(n)%sig, solvent(1)%site(n)%eps, solvent(1)%site(n)%r
            if (i/=n) then
                print*, "in solvent.in, index in first column is very strange for site number", n
                stop "have a look at read_solvent.f90"
            end if
        END DO
        CLOSE(5)

        call read_mole_fractions

        !... compute monopole, dipole, quadrupole, octupole and hexadecapole of each solvent species
        !... 1 Debye (D)  = 3.33564095 x10-30 C·m (= -0.20819435 e-·Å)
        do concurrent (s=1:size(solvent))
            !... monopole = net charge
            solvent(s)%monopole = chop(sum( solvent(s)%site%q ))

            !... dipole
            do concurrent (i=1:3)
                solvent(s)%dipole(i) = chop(sum( solvent(s)%site%q * solvent(s)%site%r(i) ))
            end do

            !... quadrupole
            do concurrent (i=1:3, j=1:3)
                solvent(s)%quadrupole(i,j) = chop(sum( solvent(s)%site%q * solvent(s)%site%r(i) * solvent(s)%site%r(j) ))
            end do

            !... octupole
            do concurrent (i=1:3, j=1:3, k=1:3)
                solvent(s)%octupole(i,j,k) = chop(&
                sum( solvent(s)%site%q * solvent(s)%site%r(i) * solvent(s)%site%r(j) * solvent(s)%site%r(k) ))
            end do

            !... hexadecapole
            do concurrent (i=1:3, j=1:3, k=1:3, l=1:3)
                solvent(s)%hexadecapole(i,j,k,l) = chop( sum( solvent(s)%site%q * &
                solvent(s)%site%r(i) * solvent(s)%site%r(j) * solvent(s)%site%r(k) * solvent(s)%site%r(l) ))
            end do
        end do

        ! Compute the charge density of a single solvant molecule in Fourier-space, and electrostatic potential generated by a such distribution
        polarization = getinput%char("polarization", defaultvalue="no")
        select case (polarization)
        case("no","none")
        case default
            call chargeDensityAndMolecularPolarizationOfASolventMoleculeAtOrigin
        end select

        ! look for bulk density of the reference solvent fluid. for instance 0.0332891 for H2O and 0.0289 for Stockmayer
        select case (solvent(1)%nspec)
        case (1)
            solvent(1)%n0= getinput%dp("bulk_density", assert=">0")
        case (2)
            solvent(1:2)%n0= getinput%dp2("bulk_density", assert=">0")
        case (3)
            solvent(1:3)%n0= getinput%dp3("bulk_density", assert=">0")
        case default
            print*, "In module solvent you are looking for the bulk density for nspec>3. Not implemented yet"
            error stop
        end select

        if (any (solvent%n0 <= 0._dp) ) then
            print *,"You ask for negative densities!"
            do s =1, solvent(1)%nspec
                print *,"For species",s,"you want density (molecule/Ang^3):",solvent(s)%n0
            end do
            stop
        end if

        solvent%rho0 = solvent%n0 / (8._dp*acos(-1._dp)**2/grid%molrotsymorder)

        call read_mole_fractions
        call functional_decision_tree

        solvent%is_initiated = .true.
    end subroutine read_solvent

    !This routine compute : -The solvent molecular charge density, which can be used into Vcoul_from_solvent_charge_density.f90 to
    !evaluate the electrostatic potential.
    !                       -The solvent molecular polarization (from Ranieriet al : J. Chem. Phys. 98 (11) 1993) that can be used
    !into energy_polarization_..._.f90 to compute the (multipolar) polarization Free energy.

    subroutine chargeDensityAndMolecularPolarizationOfASolventMoleculeAtOrigin
        implicit none
        integer :: nx, ny, nz, no, ns
        integer :: i, j, k, n, s, io, d
        real(dp)     :: r(3), kr, kvec(3)
        complex(dp)  :: fac, X
        complex(dp), parameter :: zeroc = complex(0._dp,0._dp), ic = complex(0._dp,1._dp)
        real(dp), parameter :: epsdp = epsilon(1._dp)
        type :: smoother_type
            real(dp) :: radius = 0.5_dp ! dramaticaly important
            real(dp) :: factor
        end type smoother_type
        type (smoother_type) :: smoother

        nx = grid%nx
        ny = grid%ny
        nz = grid%nz
        no = grid%no
        ns = solvent(1)%nspec

        ! sigma_k is the Fourier transformed charge density of a single solvent molecule in the reference frame defined by solvent.in
        ! molec_polar_k is the Fourier transformed molecular polarization
        do concurrent( s=1:ns , sum(abs(solvent(s)%site%q))>0) ! mask elimitates solvent molecules without point charges
            allocate( solvent(s)%sigma_k       (   nx/2+1, ny, nz, no), SOURCE=zeroC )
            allocate( solvent(s)%molec_polar_k (3, nx/2+1, ny, nz, no), SOURCE=zeroC )
        end do

        do concurrent ( i=1:nx/2+1, j=1:ny, k=1:nz, s=1:ns  , sum(abs(solvent(s)%site%q))>epsdp) ! mask elimitates solvent molecules without point charges)

            kvec = [ grid%kx(i), grid%ky(j), grid%kz(k) ]
            smoother%factor =  exp(-smoother%radius**2 * sum( kvec**2 )/2._dp)

            do concurrent ( io=1:grid%no, n=1:SIZE(solvent(s)%site), abs(solvent(s)%site(n)%q)>epsdp )
                r(1) = dot_product(   [grid%Rotxx(io),grid%Rotxy(io),grid%Rotxz(io)]  ,  solvent(s)%site(n)%r  )
                r(2) = dot_product(   [grid%Rotyx(io),grid%Rotyy(io),grid%Rotyz(io)]  ,  solvent(s)%site(n)%r  )
                r(3) = dot_product(   [grid%Rotzx(io),grid%Rotzy(io),grid%Rotzz(io)]  ,  solvent(s)%site(n)%r  )
                kr = dot_product( kvec, r )
                X = -iC*kr
                solvent(s)%sigma_k(i,j,k,io) = solvent(s)%sigma_k(i,j,k,io) + solvent(s)%site(n)%q *exp(X) *smoother%factor ! exact
                ! solvent(s)%sigma_k(i,j,k,o,p) = solvent(s)%sigma_k(i,j,k,o,p) +solvent(s)%site(n)%q* sum([(X**i/factorial(i), i=0,4)])&
                ! * smoother%factor ! Series expansion of exp(x) at 0 => multipole expansion of Vcoul(x). i=4 :: hexadecapole (16)
                if ( abs(kr)<=epsdp ) then
                    solvent(s)%molec_polar_k(:,i,j,k,io) = solvent(s)%molec_polar_k(:,i,j,k,io) + solvent(s)%site(n)%q *r
                else
                    fac = -iC*(exp(iC*kr)-1._dp)/kr *smoother%factor
                    solvent(s)%molec_polar_k(:,i,j,k,io) = solvent(s)%molec_polar_k(:,i,j,k,io) + fac*solvent(s)%site(n)%q *r
                end if
            end do
        end do

        !
        ! Substract the trace of the molecular polarization tensor
        !
        do concurrent (i=1:nx/2+1, j=1:ny, k=1:nz, s=1:ns, d=1:3)
            solvent(s)%molec_polar_k(d,i,j,k,:) = solvent(s)%molec_polar_k(d,i,j,k,:)  &
            -sum( solvent(s)%molec_polar_k(d,i,j,k,:) ) /real(grid%no,dp)
        end do
    end subroutine chargeDensityAndMolecularPolarizationOfASolventMoleculeAtOrigin


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to read the mole fractions in dft.in.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! This SUBROUTINE open the array input_line which contains every line of input/dft.in
    ! It then reads every line of input_line and looks for the tag "mole_fractions"
    ! Then, it reads, one line after the other, the mole fractions of every constituant.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE read_mole_fractions
        use precision_kinds, only: dp
        use module_input, only: input_line
        implicit none
        integer :: i, j, s
        if (.not.allocated(solvent)) stop "in read_mole_fractions, solvent is not allocated"
        select case (solvent(1)%nspec)
        case (1)
            solvent(1)%mole_fraction = 1._dp
            return
        case default
            do i = 1, size( input_line )
                j = len ( 'mole_fractions' )
                if ( input_line(i)(1:j) == 'mole_fractions' ) then
                    do s = 1, solvent(1)%nspec
                        read( input_line(i+s),*) solvent(s)%mole_fraction
                    end do
                    exit ! loop over i
                end if
            end do
        end select
        if (sum(solvent(:)%mole_fraction)/=1._dp) then
            write (*,*) 'Critial error. Sum of all mole fraction should be equal to one.'
            write (*,*) 'here are the number of the species and its associated mole fraction'
        end if
        if ( any(solvent(:)%mole_fraction<0._dp) .or. any(solvent(:)%mole_fraction>1._dp) ) THEN
            write (*,*) 'Critical errror in ALLOCATE_from_input.f90. Mole fractions should be between 0 and 1'
            write (*,*) 'here are the number of the species and its associated mole fraction'
            write (*,*) 'STOP'
            stop
        end if
    END SUBROUTINE read_mole_fractions


end module module_solvent

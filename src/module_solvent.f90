module module_solvent

    use precision_kinds
    use system, only: site_type
    use module_input, only: getinput
    use module_grid, only: grid
    use constants, only: qfact
!    use module_postprocessing !! for test only. REMOVE

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
                    exc_nn_cs_plus_nbar = .false.,&
                    exc_ck_angular = .false.
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
        real(dp) :: hs_radius !hard sphere radius for HS bridge
        real(dp) :: diameter ! hard sphere diameter, for instance
        real(dp), allocatable :: xi(:,:,:,:) ! io, ix, iy, iz    xi**2=rho/rho0
        type (site_type), allocatable :: site(:)
        real(dp)              :: n0        ! number density of the homogeneous reference fluid in molecules per Angstrom^3, e.g., 0.033291 molecule.A**-3 for water
        real(dp)              :: rho0      ! number density per orientation of the homogeneous reference fluid in molecules per Angstrom^3 per orientation
        character(80) :: quantum_description
        complex(dp), allocatable :: pseudo_charge_density_k(:,:,:,:) ! solvent molecule pseudo-charge density in k-space for QM/MM calculations
        complex(dp), allocatable :: charge_density_k(:,:,:,:) ! solvent molecule charge density, same as sigma_k with diffrent normalization
        complex(dp), allocatable :: sigma_k(:,:,:,:)! solvent molecule charge density
        complex(dp), allocatable :: molec_polar_k(:,:,:,:,:) ! solvent molecule polarization density
        real(dp), allocatable :: vext(:,:,:,:), vextq(:,:,:,:)
        real(dp) :: vext_threeshold = qfact/2.0_dp !100._dp!36.04_dp ! 36.something is the maximum value of v so that exp(-beta.v) does not return underflow at 300 K
        type(do_type) :: do
        real(dp) :: mole_fraction = 1._dp
        type(correlationfunction_type) :: cs
        type(correlationfunction_type) :: cdelta
        type(correlationfunction_type) :: cd
        complex(dp), allocatable :: ck_angular(:,:,:,:,:,:) ! TODO REMOVE THIS IS FOR TESTING PURPOSE ONLY!
        real(dp) :: relativePermittivity ! relative permittivity == static dielectric constant = dielectric constant = coonstante diélectrique
        integer:: npluc(0:6)
        integer:: n_line_cfile
    contains
        procedure, nopass :: init => read_solvent
        procedure, nopass :: init_chargedensity_molecularpolarization => &
            chargeDensityAndMolecularPolarizationOfASolventMoleculeAtOrigin
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

            !
            ! Always compute the ideal and external contributions
            !
            solvent(s)%do%id_and_ext = .true.

            !
            ! You may want to minimize only with ideal and external contributions, for instance to restart later with a not-so-bad guess
            !
            if(.not.getinput%log('minimize_wrt_vext_only',defaultvalue=.false.)) then

                solvent(s)%do%exc_fmt = getinput%log ('hard_sphere_fluid', defaultvalue=.false.)
                solvent(s)%do%exc_wca = getinput%log ('wca', defaultvalue=.false.)
                solvent(s)%do%exc_3b  = getinput%log ('threebody', defaultvalue=.false. )
                solvent(s)%do%exc_cs         = .false.
                solvent(s)%do%exc_cdeltacd   = .false.
                solvent(s)%do%exc_cproj      = .true.
                solvent(s)%do%exc_ck_angular = .false.

            else
                write(*,*) "MINIMIZING WRT VEXT ONLY *************************"
                solvent(s)%do%exc_fmt        = .false.
                solvent(s)%do%exc_wca        = .false.
                solvent(s)%do%exc_3b         = .false.
                solvent(s)%do%exc_cs         = .false.
                solvent(s)%do%exc_cdeltacd   = .false.
                solvent(s)%do%exc_cproj      = .false.
                solvent(s)%do%exc_ck_angular = .false.
            end if

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
        use module_solute, only: solute
        implicit none
        integer :: i, j, k, l, s, ncomma
        character(80) :: polarization, dummychar,solvent_quantum_description


        if (allocated(solvent)) then
            print*, "bug dans read_solvent. Le 21 octobre 2015, j'ai transféré l'allocation de allocate_from_input a read_solvent"
            print*, "ca a bcp plus de sens de le mettre dans read_solvent que dans le fourre-tout"
            print*, "cependant il est alloue alors qu'il ne devrait pas l'etre"
            stop "dans module_solvent/read_solvent"
        end if

        !
        ! How many solvent species are there ?
        !
        s= getinput%int('nb_solvent', defaultvalue=1)
        allocate( solvent(s) )
        solvent(:)%nspec = s
        if (s==1) then
          solvent(1)%name = getinput%char('solvent', defaultvalue="spce")
        else
          dummychar= getinput%char_multiple('solvent', defaultvalue="spce spce spce spce spce spce")
          read(dummychar,*) solvent(:)%name
          do i=1,s
            solvent(i)%name=Trim(adjustl(solvent(i)%name))
          end do
        end if
        
        call read_mole_fractions
        !
        ! Get the information about the solvent
        !
        !solvent(1)%name = getinput%char('solvent', defaultvalue="spce") ! This wont be valid anymore when several solvents will be used.
        Do s=1,size(solvent)
          write(*,'(A,I3,A,A,A,F7.3)') "Solvent number  " , s, " is " , solvent(s)%name," with a molecular fraction of", solvent(s)%mole_fraction
          select case (solvent(s)%name)
          case ("spce")
              solvent(s)%hs_radius = 0.0_dp
              solvent(s)%nsite = 3
              solvent(s)%molrotsymorder = 2
              allocate( solvent(s)%site(3) )
              solvent(s)%site(1:3)%q = [-0.8476, 0.4238, 0.4238]
              solvent(s)%site(1:3)%sig = [3.166, 0., 0.]
              solvent(s)%site(1:3)%eps = [0.65, 0., 0.]
              solvent(s)%site(1)%r = [0., 0., 0.]
              solvent(s)%site(2)%r = [0.816495, 0.0, 0.5773525]
              solvent(s)%site(3)%r = [-0.816495, 0.0, 0.5773525]
              solvent(s)%site(1:3)%Z = [8, 1, 1]
              solvent(s)%n0 = 0.0332891*solvent(s)%mole_fraction
              solvent(s)%rho0 = solvent(s)%n0 / (8._dp*acos(-1._dp)**2/solvent(s)%molrotsymorder)
              solvent(s)%relativePermittivity = 71._dp
              solvent(s)%npluc(0:5)=[1,6,75,252,877,2002]
              solvent(s)%n_line_cfile=1024
              if( grid%mmax>5 .or. grid%mmax<0) error stop "solvent spce only avail with mmax between 0 and 5"
          case ("spce-m")
              solvent(s)%name='spce'
              solvent(s)%hs_radius = 0.0_dp
              solvent(s)%nsite = 3
              solvent(s)%molrotsymorder = 2
              allocate( solvent(s)%site(3) )
              solvent(s)%site(1:3)%q = [-1.0, 0.5, 0.5]
              solvent(s)%site(1:3)%sig = [3.166, 0., 0.]
              solvent(s)%site(1:3)%eps = [0.65, 0., 0.]
              solvent(s)%site(1)%r = [0., 0., 0.]
              solvent(s)%site(2)%r = [0.816495, 0.0, 0.5773525]
              solvent(s)%site(3)%r = [-0.816495, 0.0, 0.5773525]
              solvent(s)%site(1:3)%Z = [8, 1, 1]
              solvent(s)%n0 = 0.0332891*solvent(s)%mole_fraction
              solvent(s)%rho0 = solvent(s)%n0 / (8._dp*acos(-1._dp)**2/solvent(s)%molrotsymorder)
              solvent(s)%relativePermittivity = 71._dp
              solvent(s)%npluc(0:5)=[1,6,75,252,877,2002]
              solvent(s)%n_line_cfile=1024
              if( grid%mmax>5 .or. grid%mmax<0) error stop "solvent spce only avail with mmax between 0 and 5"
          case ("tip3p")
              ! cf
              solvent(s)%nsite = 3
              solvent(s)%molrotsymorder = 2
              allocate( solvent(s)%site(3) )
              solvent(s)%site(1:3)%q = [-0.834, 0.417, 0.417]
              solvent(s)%site(1:3)%sig = [3.15061, 0., 0.]
              solvent(s)%site(1:3)%eps = [0.636386, 0., 0.]
              solvent(s)%site(1)%r = [0., 0., 0.]
              solvent(s)%site(2)%r = [0.756950, 0.0, 0.585882]
              solvent(s)%site(3)%r = [-0.756950, 0.0, 0.585882] 
              solvent(s)%site(1:3)%Z = [8, 1, 1]
              solvent(s)%n0 = 0.03349459*solvent(s)%mole_fraction
              solvent(s)%rho0 = solvent(s)%n0 / (8._dp*acos(-1._dp)**2/solvent(s)%molrotsymorder)
              solvent(s)%relativePermittivity = 91._dp ! cf mail de Luc du 16/12/2016 :
              solvent(s)%npluc(0:5)=[1,6,75,252,877,2002]
              solvent(s)%n_line_cfile=1024
              if( grid%mmax>5 .or. grid%mmax<0) error stop "solvent tip3p only avail with mmax between 0 and 5"
              ! Je connais ce site. C'est bizarre, la ref.3 pour epsilon(tip3p) n'a pas fait tip3p!
              ! Il y aussi J.Chem.Phys.108, 10220 (1998) qui donne 82, 94, 86 suivant N et paramètres de réaction field.
              ! Ma simulation rapide N=100 donne 100, et MC/HNC résultant donne 91.
              ! Luc
         case ("tip3p-m")
              solvent(s)%name='tip3p'
              solvent(s)%nsite = 3
              solvent(s)%molrotsymorder = 2
              allocate( solvent(s)%site(3) )
              solvent(s)%site(1:3)%q = [-0.95, 0.475, 0.475]
              solvent(s)%site(1:3)%sig = [3.15061, 0., 0.]
              solvent(s)%site(1:3)%eps = [0.636386, 0., 0.]
              solvent(s)%site(1)%r = [0., 0., 0.]
              solvent(s)%site(2)%r = [0.756950, 0.0, 0.585882]
              solvent(s)%site(3)%r = [-0.756950, 0.0, 0.585882]
              solvent(s)%site(1:3)%Z = [8, 1, 1]
              solvent(s)%n0 = 0.03349459*solvent(s)%mole_fraction
              solvent(s)%rho0 = solvent(s)%n0 / (8._dp*acos(-1._dp)**2/solvent(s)%molrotsymorder)
              solvent(s)%relativePermittivity = 91._dp ! cf mail de Luc du 16/12/2016 :
              solvent(s)%npluc(0:5)=[1,6,75,252,877,2002]
              solvent(s)%n_line_cfile=1024
              if( grid%mmax>5 .or. grid%mmax<0) error stop "solvent tip3p only avail with mmax between 0 and 5"
          case ("acetonitrile")
              ! Reference: Edwards, Madden and McDonald, doi:10.1080/00268978400100731
              solvent(s)%nsite = 3 ! z<---Me---C--N--
              solvent(s)%molrotsymorder = 1000
              allocate( solvent(s)%site(3) )
              solvent(s)%site(1:3)%q = [0.269, 0.129, -0.398]
              solvent(s)%site(1:3)%sig = [3.6, 3.4, 3.3]
              solvent(s)%site(1:3)%eps = [1.59, 0.416, 0.416]
              solvent(s)%site(1)%r = [0., 0., 1.46]  !C CH3
              solvent(s)%site(2)%r = [0., 0., 0.]    !C N
              solvent(s)%site(3)%r = [0., 0., -1.17] !N
              solvent(s)%site(1:3)%Z = [9, 6, 7]
              !solvent(1)%n0 = 0.0289
              solvent(s)%n0 = 0.012044*solvent(s)%mole_fraction
              solvent(s)%rho0 = solvent(s)%n0/ (8._dp*acos(-1._dp)**2/solvent(s)%molrotsymorder)
              !solvent(1)%rho0 = solvent(1)%n0/ (8._dp*acos(-1._dp)**2)
              solvent(s)%relativePermittivity = 31.6_dp ! TODO TO BE CHECKED AND INCLUDED.
              solvent(s)%npluc(0:6)=[1,6, 19, 44, 85, 146, 231]
              solvent(s)%n_line_cfile=500
              if( grid%mmax>6 .or. grid%mmax<0) error stop "mmax is not between 0 and 6"
          case ("acn6_hnc")
             ! Reference: Bohm, Mc Donald, Madden
              solvent(s)%nsite = 6 ! z<---CH3---C--N--
              solvent(s)%molrotsymorder = 3
              allocate( solvent(s)%site(6) )
              solvent(s)%site(1:6)%q = [-0.577, 0.488, -0.514, 0.201, 0.201, 0.201]
              solvent(s)%site(1:6)%sig = [3.0, 3.4, 3.3, 2.2, 2.2, 2.2]
              solvent(s)%site(1:6)%eps = [0.4177, 0.4177, 0.4177, 0.0835, 0.0835, 0.0835]
              solvent(s)%site(1)%r = [0., 0., 1.46]  !C CH3
              solvent(s)%site(2)%r = [0., 0., 0.]    !C N
              solvent(s)%site(3)%r = [0., 0., -1.17] !N
              solvent(s)%site(4)%r = [1.02274, 0.0, 1.8282]  !H1
              solvent(s)%site(5)%r = [-0.51137,  0.88572, 1.8282] ! H2
              solvent(s)%site(6)%r = [-0.51137, -0.88572, 1.8282] ! H3
              solvent(s)%site(1:6)%Z = [9, 6, 7, 1, 1, 1]
              solvent(s)%n0 = 0.011694635*solvent(s)%mole_fraction
              solvent(s)%rho0 = solvent(s)%n0/ (8._dp*acos(-1._dp)**2/solvent(s)%molrotsymorder)
              solvent(s)%relativePermittivity = 26.38_dp ! from Luc
              solvent(s)%npluc(0:6)=[1, 6, 19, 136, 397, 244, 571]
              solvent(s)%n_line_cfile=512
              if( grid%mmax>4 .or. grid%mmax<0) error stop "mmax is not between 0 and 4"
          case default
               print*,  "the solvent you want to use: ", trim(solvent(s)%name), " is not available, you can only use spce, tip3p or acetonitrile for now, sorry :("
               stop
          end select
        end do
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

! Compute solvent molecule charge-density in reciprocal space (solvent(1)%sigma_k

        call Get_solvent_molecule_reciprocal_charge_density
  !      call chargeDensityAndMolecularPolarizationOfASolventMoleculeAtOrigin


! Compute solvent molecule pseudo charge-density for QM calculations (added by daniel, 7-12-2018)
COMPUTE_SOLVENT_MOLECULE_RECIPROCAL_PSEUDO_CHARGE_DENSITY: BLOCK
     allocate( solvent(1)%pseudo_charge_density_k(grid%nx/2+1, grid%ny, grid%nz, grid%no), SOURCE=(0._dp,0._dp) )

     if( solute%nature == "QM") then

        select case (solute%solvent_coupling)
        case ("pseudo")
              call Get_solvent_molecule_reciprocal_pseudo_charge_density
            print*, 'solvent_quantum_representation: pseudo'
        case ("classical_charges")
              solvent(1)%pseudo_charge_density_k(:,:,:,:) = solvent(1)%sigma_k(:,:,:,:)
             print*, 'solvent_quantum_representation: classical_charges'
        end select

     else
           solvent(1)%pseudo_charge_density_k(:,:,:,:) = solvent(1)%sigma_k(:,:,:,:)
!           print*, 'if classical solute: solvent_quantum_representation= classical_charges'

     end if

END BLOCK COMPUTE_SOLVENT_MOLECULE_RECIPROCAL_PSEUDO_CHARGE_DENSITY




        call functional_decision_tree

        solvent%is_initiated = .true.
        if (any(solvent%molrotsymorder/=grid%molrotsymorder)) then
          print*, "########################################################"
          print*, "WARNING WARNING WARNING WARNING WARNING WARNING WARNING"
          print*, "the grid molrotsymorder is", grid%molrotsymorder, "which differs from at least one solvent molrotsymorder that are"
          print*, solvent(:)%molrotsymorder, "be very carefull I did not tested what you are trying to do" 
          print*, "WARNING WARNING WARNING WARNING WARNING WARNING WARNING"
          print*, "########################################################"
        end if
    end subroutine read_solvent

    !This routine compute : -The solvent molecular charge density, which can be used into Vcoul_from_solvent_charge_density.f90 to
    !evaluate the electrostatic potential.
    !                       -The solvent molecular polarization (from Ranieriet al : J. Chem. Phys. 98 (11) 1993) that can be used
    !into energy_polarization_..._.f90 to compute the (multipolar) polarization Free energy.

    subroutine chargeDensityAndMolecularPolarizationOfASolventMoleculeAtOrigin
        use module_grid, only: grid
        implicit none
        integer :: nx, ny, nz, no, ns
        integer :: i, j, k, n, s, io, d
        real(dp)     :: r(3), kr, kvec(3)
        complex(dp)  :: fac, X
        complex(dp), parameter :: zeroc = (0._dp,0._dp), ic = (0._dp,1._dp)
        real(dp), parameter :: epsdp = epsilon(1._dp)
        real(dp) :: smootherfactor
        real(dp) :: smootherradius = 0.2_dp ! dramaticaly important
        nx = grid%nx
        ny = grid%ny
        nz = grid%nz
        no = grid%no
        ns = size(solvent) ! Count of solvent species

    ! sigma_k is the Fourier transformed charge density of a single solvent molecule in the reference frame defined by solvent.in
        ! molec_polar_k is the Fourier transformed molecular polarization
        do s = 1, ns
            if( sum( abs( solvent(s)%site%q )) > epsdp ) then
                allocate( solvent(s)%sigma_k       (   nx/2+1, ny, nz, no), SOURCE=zeroC )
                allocate( solvent(s)%molec_polar_k (3, nx/2+1, ny, nz, no), SOURCE=zeroC )
            end if
        end do
        
        do s=1,ns  
           if( sum( abs( solvent(s)%site%q )) <= epsdp ) cycle ! Don't compute the polarization of a solvent molecule that has no point charge
           !$omp parallel private(i, j, k, kvec, smootherfactor, r, kr, X, fac)
           !$omp do
           do k = 1, nz
              do j = 1, ny
                 do i = 1, nx/2+1
                    kvec = [ grid%kx(i), grid%ky(j), grid%kz(k) ]
                    smootherfactor =  exp(-smootherradius**2 * sum( kvec**2 )/2._dp)
                    
                    do n = 1, SIZE(solvent(s)%site)
                       do io = 1, grid%no
                          if ( abs(solvent(s)%site(n)%q) > epsdp ) then
                             r(1) = dot_product(   [grid%Rotxx(io),grid%Rotxy(io),grid%Rotxz(io)]  ,  solvent(s)%site(n)%r  )
                             r(2) = dot_product(   [grid%Rotyx(io),grid%Rotyy(io),grid%Rotyz(io)]  ,  solvent(s)%site(n)%r  )
                             r(3) = dot_product(   [grid%Rotzx(io),grid%Rotzy(io),grid%Rotzz(io)]  ,  solvent(s)%site(n)%r  )
                             kr = dot_product( kvec, r )
                             X = -iC*kr
                             solvent(s)%sigma_k(i,j,k,io) = solvent(s)%sigma_k(i,j,k,io) + solvent(s)%site(n)%q *exp(X) *smootherfactor ! exact
                             ! solvent(s)%sigma_k(i,j,k,o,p) = solvent(s)%sigma_k(i,j,k,o,p) +solvent(s)%site(n)%q* sum([(X**i/factorial(i), i=0,4)])&
                             ! * smootherfactor ! Series expansion of exp(x) at 0 => multipole expansion of Vcoul(x). i=4 :: hexadecapole (16)
                             if ( abs(kr)<=epsdp ) then
                                solvent(s)%molec_polar_k(:,i,j,k,io) = solvent(s)%molec_polar_k(:,i,j,k,io) + solvent(s)%site(n)%q *r
                             else
                                fac = -iC*(exp(iC*kr)-1._dp)/kr *smootherfactor
                                solvent(s)%molec_polar_k(:,i,j,k,io) = solvent(s)%molec_polar_k(:,i,j,k,io) + fac*solvent(s)%site(n)%q *r
                             end if
                          end if
                       end do
                    end do

                 end do
              end do
           end do
           !$omp end do
           !$omp end parallel
        end do
        !
        ! Substract the trace of the molecular polarization tensor
        !
        do concurrent (i=1:nx/2+1, j=1:ny, k=1:nz, s=1:ns, d=1:3)
            if( sum( abs( solvent(s)%site%q )) > epsdp ) cycle
            solvent(s)%molec_polar_k(d,i,j,k,:) = solvent(s)%molec_polar_k(d,i,j,k,:)  &
                -sum( solvent(s)%molec_polar_k(d,i,j,k,:) ) /real(grid%no,dp)
        end do
      end subroutine chargeDensityAndMolecularPolarizationOfASolventMoleculeAtOrigin

subroutine Get_solvent_molecule_reciprocal_pseudo_charge_density
use module_grid, only: grid
implicit none
integer :: nx, ny, nz, no, ns
integer :: i, j, k, n, s, io, d
real(dp)     :: r(3), kr, kvec(3)
complex(dp)  :: fac, X
real(dp) :: GaussianFactor
complex(dp), parameter :: zeroc = (0._dp,0._dp), ic = (0._dp,1._dp)
real(dp), parameter :: AngtoBohr = 1.889725989_dp, BohrtoAng = 1.0_dp/AngtoBohr
real(dp), parameter :: q_h = 0.4238_dp, A_1h = 0.75_dp/BohrtoAng, B_1h = 0.150_dp, B_2h = 0.5_dp/BohrtoAng, B_3h = 0.35_dp/BohrtoAng
real(dp), parameter :: q_o = -0.8476_dp, A_1o = 0.575_dp/BohrtoAng, B_1o = 0.620_dp, B_2o = 1.0_dp/BohrtoAng, B_3o = 0.4_dp/BohrtoAng
!real(dp), parameter :: q_o = -0.8476_dp, A_1o = 0.575_dp/BohrtoAng, B_1o = 3.0*0.620_dp, B_2o = 1.0_dp/BohrtoAng, B_3o = 0.4_dp/BohrtoAng
real(dp), parameter :: C_1o = 0.03195, C_2o = 0.48129_dp/BohrtoAng, C_3o = 0.12092_dp/BohrtoAng

nx = grid%nx
ny = grid%ny
nz = grid%nz
no = grid%no
ns = size(solvent) ! Count of solvent species

!print*, 'passed in subroutine get_solvent_molecule_pseudo_charge_density, created by Daniel on 7-12-2018 for introducing QM/MM electron-water pseudopotential'

! At this stage: pseudo_charge_density is the Fourier transformed charge density of a single water molecule in the reference frame defined by solvent.in

if( ns /= 1) stop 'init_solvent_molecule_pseudo_charge_density only works for ns = 1'
if( solvent(1)%name /= 'spce' ) stop 'init_solvent_molecule_pseudo_charge_density only work for spc or spce'


!$omp parallel private(i, j, k, kvec, smootherfactor, r, kr, X, fac)
!$omp do

  do io = 1, no

    do k = 1, nz
      do j = 1, ny
        do i = 1, nx/2+1
           kvec = [ grid%kx(i), grid%ky(j), grid%kz(k) ]

           n= 1 ! oxygen site

           GaussianFactor =  q_o*exp( -sum(kvec**2 )/(4._dp*A_1o**2)) &
                     - B_1o*( exp(-sum( kvec**2 )/(4._dp*B_2o**2)) - exp(-sum( kvec**2 )/(4._dp*B_3o**2)) ) &
                     + C_1o*( exp(-sum( kvec**2 )/(4._dp*C_2o**2)) - exp(-sum( kvec**2 )/(4._dp*C_3o**2)) )

           r(1) = dot_product(   [grid%Rotxx(io),grid%Rotxy(io),grid%Rotxz(io)]  ,  solvent(1)%site(n)%r  )
           r(2) = dot_product(   [grid%Rotyx(io),grid%Rotyy(io),grid%Rotyz(io)]  ,  solvent(1)%site(n)%r  )
           r(3) = dot_product(   [grid%Rotzx(io),grid%Rotzy(io),grid%Rotzz(io)]  ,  solvent(1)%site(n)%r  )
           kr = dot_product( kvec, r )
           X = -iC*kr
           solvent(1)%pseudo_charge_density_k(i,j,k,io) = exp(X) *GaussianFactor!! exact

        do n= 2, 3 ! sum of hydrogen sites site

           GaussianFactor =  q_h*exp(-sum( kvec**2 )/(4._dp*A_1h**2))&
                - B_1h*( exp(-sum( kvec**2 )/(4._dp*B_2h**2)) - exp(-sum( kvec**2 )/(4._dp*B_3h**2)) )

            r(1) = dot_product(   [grid%Rotxx(io),grid%Rotxy(io),grid%Rotxz(io)]  ,  solvent(1)%site(n)%r  )
            r(2) = dot_product(   [grid%Rotyx(io),grid%Rotyy(io),grid%Rotyz(io)]  ,  solvent(1)%site(n)%r  )
            r(3) = dot_product(   [grid%Rotzx(io),grid%Rotzy(io),grid%Rotzz(io)]  ,  solvent(1)%site(n)%r  )
            kr = dot_product( kvec, r )
            X = -iC*kr
            solvent(1)%pseudo_charge_density_k(i,j,k,io) =  solvent(1)%pseudo_charge_density_k(i,j,k,io) &
                                 + exp(X) *GaussianFactor! exact

           end do ! loop over hydrogen sites


        end do !loop nx
      end do ! loop ny
    end do !loop nz

  end do !loop no


!$omp end do
!$omp end parallel
end subroutine Get_solvent_molecule_reciprocal_pseudo_charge_density

subroutine Get_solvent_molecule_reciprocal_charge_density  !not useful for now, same as chargeDensityAndMolecularPolarizationOfASolventMoleculeAtOrigin exept for normalization by grid%dv. charge_density_k is sigma_k
use module_grid, only: grid
implicit none
integer :: nx, ny, nz, no, ns
integer :: i, j, k, n, s, io, d
real(dp)     :: r(3), kr, kvec(3)
complex(dp)  :: fac, X
complex(dp), parameter :: zeroc = (0._dp,0._dp), ic = (0._dp,1._dp)
real(dp), parameter :: epsdp = epsilon(1._dp)
real(dp) :: smootherfactor
real(dp) :: smootherradius = 0.5_dp ! dramaticaly important
nx = grid%nx
ny = grid%ny
nz = grid%nz
no = grid%no
ns = size(solvent) ! Count of solvent species

!print*, 'passed in subroutine init_solvent_molecule_charge_density, created by Daniel on 7-12-2018, same as chargeDensityAndMolecularPolarizationOfASolventMoleculeAtOrigin'

!charge_density is the Fourier transformed "classical" charge density of a single water molecule in the reference frame defined by solvent.in

if( ns /= 1) stop 'init_solvent_molecule_pseudo_charge_density only works for ns = 1'

allocate( solvent(1)%sigma_k(nx/2+1, ny, nz, no), SOURCE=zeroC )

!$omp parallel private(i, j, k, kvec, smootherfactor, r, kr, X, fac)
!$omp do
do io = 1, no

do k = 1, nz
do j = 1, ny
do i = 1, nx/2+1

do n=1, solvent(1)%nsite

kvec = [ grid%kx(i), grid%ky(j), grid%kz(k) ]
smootherfactor =  exp(-smootherradius**2 * sum( kvec**2 )/2._dp)

r(1) = dot_product(   [grid%Rotxx(io),grid%Rotxy(io),grid%Rotxz(io)]  ,  solvent(1)%site(n)%r  )   !  - grid%length(1)/2._dp
r(2) = dot_product(   [grid%Rotyx(io),grid%Rotyy(io),grid%Rotyz(io)]  ,  solvent(1)%site(n)%r  )   !  - grid%length(2)/2._dp
r(3) = dot_product(   [grid%Rotzx(io),grid%Rotzy(io),grid%Rotzz(io)]  ,  solvent(1)%site(n)%r  )   !  - grid%length(3)/2._dp
kr = dot_product( kvec, r )
X = -iC*kr
solvent(1)%sigma_k(i,j,k,io) = solvent(1)%sigma_k(i,j,k,io) &
                         + solvent(1)%site(n)%q *exp(X) *smootherfactor ! exact
end do ! loop over solvent sites


end do !loop nx
end do ! loop ny
end do !loop nz

end do !loop no
!$omp end do
!$omp end parallel

end subroutine Get_solvent_molecule_reciprocal_charge_density


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
        logical :: tagfound
        tagfound=.false.
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
                    tagfound=.true.
                    exit ! loop over i
                end if
            end do
            if (.not. tagfound) then
              print*, "you are using several solvent and you are not specifying the molefraction, I will assume it is equimolar i.e each are", 1.0/real(solvent(1)%nspec,dp) 
              solvent(:)%mole_fraction=1.0/real(solvent(1)%nspec,dp)
            end if
        end select
        if (sum(solvent(:)%mole_fraction)/=1._dp) then
            write (*,*) 'Critial error. Sum of all mole fraction should be equal to one.'
            write (*,*) 'here are the number of the species and its associated mole fraction'
            stop
        end if
        if ( any(solvent(:)%mole_fraction<0._dp) .or. any(solvent(:)%mole_fraction>1._dp) ) THEN
            write (*,*) 'Critical errror in ALLOCATE_from_input.f90. Mole fractions should be between 0 and 1'
            write (*,*) 'here are the number of the species and its associated mole fraction'
            write (*,*) 'STOP'
            stop
        end if
    END SUBROUTINE read_mole_fractions


end module module_solvent

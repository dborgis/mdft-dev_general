! Variable relative to the system being studied
MODULE system

    USE precision_kinds , ONLY: i2b,dp

    IMPLICIT NONE

    INTEGER(i2b) :: nb_species ! number of solvents in the species, e.g. 2 if the solvent is a mixture of water and acetone
    INTEGER(i2b) :: nb_solute_sites, nb_solvent_sites  ! nombre de site pour le solute et pour le solvent

    TYPE :: sites
        CHARACTER(100) :: name
        INTEGER(i2b) :: type, n_sites
        REAL(dp), DIMENSION(3) :: r
        REAL(dp) :: q, sig, eps, lambda1, lambda2
        INTEGER(i2b) :: Z ! atomic number
    END TYPE sites
    
    TYPE (sites), ALLOCATABLE, DIMENSION(:) :: soluteSite
    
    TYPE :: solventType
        type (sites), allocatable :: site(:)
        real(dp), allocatable :: n(:,:,:)  ! number density
        real(dp)              :: n0        ! number density of the homogeneous reference fluid in molecules per Angstrom^3, e.g., 0.033291 molecule.A**-3 for water
        real(dp), allocatable :: Dn(:,:,:) ! Dn = n - n0
        real(dp), allocatable :: rho(:,:,:,:,:)! number density per orientation = n0/(8pi²/molrotsymorder)
        real(dp)              :: rho0      ! number density per orientation of the homogeneous reference fluid in molecules per Angstrom^3 per orientation
        real(dp), allocatable :: Drho(:,:,:,:,:) ! Drho = rho - rho0
    END TYPE
    TYPE (solventType), ALLOCATABLE :: solvent(:)
    
    TYPE :: thermoCondType
        REAL(dp) :: T ! temperature
        REAL(dp) :: kbT ! temperature energy unit
        REAL(dp) :: beta ! 1/kbT
    END TYPE
    TYPE (thermoCondType) :: thermoCond ! everything related to the thermodynamic conditions
    
    TYPE :: spaceGridType
        INTEGER(i2b), DIMENSION(3) :: n_nodes ! number of grid nodes in direction x, y and z
        REAL(dp), DIMENSION(3) :: length ! total length in direction x, y and z
        REAL(dp), DIMENSION(3) :: dl ! elemental distance between two nodes in direction x, y and z
        REAL(dp) :: dv ! elemental volume
    END TYPE spaceGridType
    TYPE( spaceGridType ), TARGET :: spaceGrid
    REAL(dp), POINTER :: Lx => spaceGrid%length(1), Ly => spaceGrid%length(2), Lz => spaceGrid%length(3)
    REAL(dp), POINTER :: DeltaV => spaceGrid%dv
    REAL(dp), POINTER :: deltax => spaceGrid%dl(1), deltay => spaceGrid%dl(2), deltaz => spaceGrid%dl(3)
    INTEGER(i2b), POINTER :: nfft1 => spaceGrid%n_nodes(1), nfft2 => spaceGrid%n_nodes(2), nfft3 => spaceGrid%n_nodes(3) ! deprecated. Should be removed at some point

    REAL(dp), ALLOCATABLE, DIMENSION (:) :: c_s_hs ! c(2)(k) of a hard sphere
    COMPLEX(dp),ALLOCATABLE, DIMENSION(:,:,:) :: Vk !>@var perturabtion in kspace


    ! Electrostatics
    ! charge factor & molecule polarization factor
    COMPLEX(dp), ALLOCATABLE, DIMENSION (:,:,:,:,:,:) :: sigma_k
    !Polarization
    COMPLEX(dp), ALLOCATABLE, DIMENSION (:,:,:,:,:,:) ::  molec_polarx_k, molec_polary_k, molec_polarz_k 
    
    !> 3body parameters
    !NEW
    !!a=100.0_dp
    !!gam=2.0_dp/3.0_dp
    !!rmin1=1.5_dp
    !!rsw1=2.0_dp
    !!rmax1=0.5_dp*(sig_mol(1)+sig_solv(1))+1.9_dp
    !!rmin2=2.25_dp
    !!rsw2=2.5_dp
    !!rmax2=4.75_dp
    !!lambda_fs=100.0_dp
    !OLD
    !gamma  :  2.9
    !rsw1  :  2.0
    !rsw2  :  2.75
    !rmin1  :  1.5
    !rmax1  :  5.0
    !rmin2  :  2.25
    !rmax2  :  5.00
    !lambda  :  125.00
    !lambda_fs  :  0.0
    !a3  :  0.0
    !b3  :  3.0
    !REAL(dp), ALLOCATABLE , DIMENSION (:,:,:, : ) :: rho_n ! density per angle (=n/4pi) of each species at each node (nfft1,nfft2,nfft3,nb_species)
    ! for maximum efficiency, rho_species(nfft1,nfft2,nfft3,nb_species) should be used in loops where species is the outer loop (varies slowliest)
    !do species = 1 , nb_species
    !do z
    !do y 
    !do x
    !density(x,y,z,species)
    !> Hard spheres
    ! TODO this should be put in another module
    logical(4) :: HS
    COMPLEX(dp), ALLOCATABLE , DIMENSION (:,:,:) :: v_perturbation_k ! fourier transform of the lennard jones perturbation (WCA)
    REAL(dp), ALLOCATABLE , DIMENSION (:) :: mole_fraction ! mole fraction of each species "x_i"
    !> for xsf printing


END MODULE system

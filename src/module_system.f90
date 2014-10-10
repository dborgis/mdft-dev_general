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
    
    TYPE :: solventType
        type (sites), allocatable :: site(:)
        real(dp), allocatable :: n(:,:,:)  ! number density
        real(dp)              :: n0        ! number density of the homogeneous reference fluid in molecules per Angstrom^3, e.g., 0.033291 molecule.A**-3 for water
        real(dp), allocatable :: Dn(:,:,:) ! Dn = n - n0
        real(dp), allocatable :: rho(:,:,:,:,:)! number density per orientation = n0/(8pi²/molrotsymorder)
        real(dp)              :: rho0      ! number density per orientation of the homogeneous reference fluid in molecules per Angstrom^3 per orientation
        real(dp), allocatable :: Drho(:,:,:,:,:) ! Drho = rho - rho0
        complex(dp), allocatable :: sigma_k(:,:,:,:,:) ! charge factor
        complex(dp), allocatable, dimension(:,:,:,:,:,:) :: molec_polar_k ! molecule polarization factor
    END TYPE
    TYPE (solventType), ALLOCATABLE :: solvent(:)
    TYPE (solventType) :: solute
    
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

    REAL(dp), ALLOCATABLE, DIMENSION (:) :: c_s_hs ! c(2)(k) of a hard sphere
    COMPLEX(dp),ALLOCATABLE, DIMENSION(:,:,:) :: Vk !>@var perturabtion in kspace

    ! Electrostatics
    ! charge factor & molecule polarization factor
    COMPLEX(dp), ALLOCATABLE, DIMENSION (:,:,:,:,:,:) :: sigma_k
    !Polarization
    COMPLEX(dp), ALLOCATABLE, DIMENSION (:,:,:,:,:,:) ::  molec_polarx_k, molec_polary_k, molec_polarz_k 
    
    logical(4) :: HS
    COMPLEX(dp), ALLOCATABLE , DIMENSION (:,:,:) :: v_perturbation_k ! fourier transform of the lennard jones perturbation (WCA)
    REAL(dp), ALLOCATABLE , DIMENSION (:) :: mole_fraction ! mole fraction of each species "x_i"
    !> for xsf printing

    type :: enertype
        real(dp) :: tot, ext, id, fmt, fmtcs, excnn, excpol, b31, b32, exc_ck_angular
    end type enertype
    type (enertype) :: FF

END MODULE system

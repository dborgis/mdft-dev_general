! Variable relative to the system being studied
MODULE system

    USE precision_kinds , ONLY: i2b,dp

    IMPLICIT NONE

    INTEGER(i2b) :: nb_solute_sites, nb_solvent_sites  ! nombre de site pour le solute et pour le solvent

    TYPE :: sites
        CHARACTER(100) :: name
        INTEGER(i2b) :: type, n_sites
        REAL(dp), DIMENSION(3) :: r
        REAL(dp) :: q, sig, eps, lambda1, lambda2
        INTEGER(i2b) :: Z ! atomic number
    END TYPE sites

    type :: vextType
        real(dp) :: tot
        real(dp) :: lj
        real(dp) :: q, qold
        real(dp) :: h
    end type vextType

    TYPE :: solventType
        integer :: nspec ! number of solvent species
        real(dp) :: monopole, dipole(3), quadrupole(3,3), octupole(3,3,3), hexadecapole(3,3,3,3)
        real(dp) :: diameter ! hard sphere diameter, for instance
        type (sites), allocatable :: site(:)
        real(dp), allocatable :: n(:,:,:)  ! number density
        real(dp)              :: n0        ! number density of the homogeneous reference fluid in molecules per Angstrom^3, e.g., 0.033291 molecule.A**-3 for water
        real(dp), allocatable :: Dn(:,:,:) ! Dn = n - n0
        real(dp), allocatable :: rho(:,:,:,:)! number density per orientation = n0/(8pi²/molrotsymorder)
        real(dp)              :: rho0      ! number density per orientation of the homogeneous reference fluid in molecules per Angstrom^3 per orientation
        real(dp), allocatable :: Drho(:,:,:,:) ! Drho = rho - rho0
        complex(dp), allocatable :: sigma_k(:,:,:,:) ! charge factor
        complex(dp), allocatable :: molec_polar_k(:,:,:,:,:) ! molecule polarization factor
        type(vextType), allocatable :: vext(:,:,:,:) ! nfft1,nfft2,nfft3,orientation
    END TYPE

    TYPE (solventType), ALLOCATABLE :: solvent(:)
    TYPE (solventType) :: solute

    TYPE :: thermoCondType
        REAL(dp) :: T ! temperature
        REAL(dp) :: kbT ! temperature energy unit
        REAL(dp) :: beta ! 1/kbT
    END TYPE
    TYPE (thermoCondType) :: thermoCond ! everything related to the thermodynamic conditions


    ! type :: somegrid
    !     !
    !     ! Spatial grid
    !     !
    !     integer, dimension(3) :: n_nodes, n ! number of grid nodes in direction x, y and z
    !     integer :: nx, ny, nz
    !     real(dp), dimension(3) :: length, l ! total length in direction x, y and z
    !     real(dp) :: lx, ly, lz
    !     real(dp), dimension(3) :: dl ! elemental distance between two nodes in direction x, y and z
    !     real(dp) :: dx, dy, dz
    !     real(dp) :: dv ! elemental volume
    !     real(dp) :: v
    !     real(dp) :: buffer_length ! length of free space between the extremam of the solute.
    !     !
    !     ! Angular grid .. angular quadrature
    !     !
    !     integer :: molrotsymorder, mmax, ntheta, nphi, npsi, no
    !     real(dp) :: dphi, dpsi
    !     real(dp), allocatable :: theta(:), phi(:), psi(:), wtheta(:), wphi(:), wpsi(:), w(:)
    !     integer, allocatable :: tio(:,:,:) ! table of index of orientations
    !     real(dp), allocatable, dimension(:) :: rotxx, rotxy, rotxz, rotyx, rotyy, rotyz, rotzx, rotzy, rotzz
    !     real(dp), allocatable, dimension(:) :: OMx, OMy, OMz
    ! end type somegrid
    ! type(somegrid), target, public :: grid ! TODO remove target. Was used for retrocompatibility reason



    COMPLEX(dp),ALLOCATABLE, DIMENSION(:,:,:) :: Vk !>@var perturabtion in kspace

    logical(4) :: HS
    COMPLEX(dp), ALLOCATABLE , DIMENSION (:,:,:) :: v_perturbation_k ! fourier transform of the lennard jones perturbation (WCA)
    REAL(dp), ALLOCATABLE , DIMENSION (:) :: mole_fraction ! mole fraction of each species "x_i"
    !> for xsf printing

    type :: enertype
        real(dp) :: tot, ext, id, fmt, fmtcs, excnn, excpol, b31, b32, exc_ck_angular
    end type enertype
    type (enertype) :: FF

END MODULE system

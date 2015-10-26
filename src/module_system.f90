! Variable relative to the system being studied
MODULE system
    use precision_kinds , ONLY: i2b,dp
    implicit none
    type :: site_type
        CHARACTER(100) :: name
        INTEGER(i2b) :: type, n_sites
        REAL(dp), DIMENSION(3) :: r
        REAL(dp) :: q, sig, eps, lambda1, lambda2
        INTEGER(i2b) :: Z ! atomic number
    end type site_type
    type :: vextType
        real(dp) :: tot
        real(dp) :: lj
        real(dp) :: q, qold
        real(dp) :: h
    end type vextType
    type :: thermoCondType
        REAL(dp) :: T ! temperature
        REAL(dp) :: kbT ! temperature energy unit
        REAL(dp) :: beta ! 1/kbT
    END TYPE
    type (thermoCondType) :: thermoCond ! everything related to the thermodynamic conditions
    private
    public :: site_type, thermocond, mole_fraction


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

    ! type :: enertype
    !     real(dp) :: tot, ext, id, fmt, fmtcs, excnn, excpol, b31, b32, exc_ck_angular
    ! end type enertype
    ! type (enertype) :: FF

END MODULE system

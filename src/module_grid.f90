module module_grid
    use precision_kinds, only: dp
    implicit none
    private
    type :: somegrid
        !
        ! Spatial grid
        !
        integer, dimension(3) :: n_nodes, n ! number of grid nodes in direction x, y and z
        integer :: nx, ny, nz
        real(dp), dimension(3) :: length, l ! total length in direction x, y and z
        real(dp) :: lx, ly, lz
        real(dp), dimension(3) :: dl ! elemental distance between two nodes in direction x, y and z
        real(dp) :: dx, dy, dz
        real(dp) :: dv ! elemental volume
        real(dp) :: v
        real(dp) :: buffer_length ! length of free space between the extremam of the solute.
        !
        ! Angular grid .. angular quadrature
        !
        integer :: molrotsymorder, mmax, ntheta, nphi, npsi, no
        real(dp) :: dphi, dpsi
        real(dp), allocatable :: theta(:), phi(:), psi(:), wtheta(:), wphi(:), wpsi(:), w(:)
        integer, allocatable :: tio(:,:,:) ! table of index of orientations
        real(dp), allocatable, dimension(:) :: rotxx, rotxy, rotxz, rotyx, rotyy, rotyz, rotzx, rotzy, rotzz
        real(dp), allocatable, dimension(:) :: OMx, OMy, OMz
    end type somegrid
    type(somegrid), target, public :: grid ! TODO remove target. Was used for retrocompatibility reason
end module module_grid

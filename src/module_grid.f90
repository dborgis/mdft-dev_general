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
    type(somegrid), public :: grid ! TODO remove target. Was used for retrocompatibility reason

    public :: init_grid

contains
    subroutine init_grid
        use input, only: input_int, input_int3, input_dp, input_dp3
        grid%molRotSymOrder = input_int('molRotSymOrder', defaultvalue=1) !Get the order of the main symmetry axis of the solvent
        grid%length = input_dp3( "lxlylz" , defaultvalue= grid%length )
        if (ANY( grid%length  <= 0._dp ) ) THEN
            PRINT*,'The supercell cannot have negative length.'
            PRINT*,'Here are your Lx, Ly and Lz as defined in input/dft.in :',grid%length
            STOP
        end if
        grid%l = grid%length
        grid%lx = grid%length(1)
        grid%ly = grid%length(2)
        grid%lz = grid%length(3)

        grid%n_nodes = input_int3( "nxnynz" , defaultvalue= nint(grid%length/0.3_dp) )
        if ( any(grid%n_nodes <= 0) ) then
            print*, 'The space is divided into grid nodes. For each direction, you ask', grid%n_nodes,'node.'
            error stop
        end if
        grid%n = grid%n_nodes
        grid%nx = grid%n(1)
        grid%ny = grid%n(2)
        grid%nz = grid%n(3)

        grid%dl = grid%length / real(grid%n,dp) !
        grid%dx = grid%dl(1)
        grid%dy = grid%dl(2)
        grid%dz = grid%dl(3)

        grid%v = product(grid%length)
        grid%dv = product(grid%dl)

        ! We now have a full description of the space grid
        print*,
        print*, "[GRID]====="
        print*, "Box Length :", grid%length
        print*, "nodes      :", grid%n_nodes
        print*, "dx, dy, dz :", grid%dl
        print*, "[/GRID]===="
        print*,

    end subroutine init_grid
end module module_grid

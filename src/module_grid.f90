module module_grid
    use precision_kinds, only: dp
    implicit none
    private
    type :: somegrid
        logical :: isinitiated = .false. ! Once the grid is initiated, change this value to .true.
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
        real(dp), allocatable, dimension(:) :: kx, ky, kz
        !
        ! Angular grid .. angular quadrature
        !
        integer :: molrotsymorder, mmax, ntheta, nphi, npsi, no
        real(dp) :: dphi, dpsi
        real(dp), allocatable :: theta(:), phi(:), psi(:), wtheta(:), wphi(:), wpsi(:), w(:)
        integer, allocatable :: tio(:,:,:) ! table of index of orientations
        real(dp), allocatable, dimension(:) :: rotxx, rotxy, rotxz, rotyx, rotyy, rotyz, rotzx, rotzy, rotzz
        real(dp), allocatable, dimension(:) :: OMx, OMy, OMz
    contains
        procedure, nopass :: init
    end type somegrid
    type(somegrid), public :: grid

    public :: norm_k, timesExpPrefactork2, k2!,init_grid

contains

    subroutine init
        use module_input, only: getinput
        implicit none
        if (grid%isinitiated) then
            print*, "Dans init_grid, c'est bizarre. On veut initialiser le type derivé grid mais il semble deja initialisé"
            stop "dans module_grid/init_grid "
        end if
        grid%molRotSymOrder = getinput%int('molRotSymOrder', defaultvalue=1, assert=">0") !Get the order of the main symmetry axis of the solvent
        grid%length = getinput%dp3( "boxlen" , defaultvalue=[128._dp,128._dp,128._dp], assert=">0" )
        if (ANY( grid%length  <= 0._dp ) ) THEN
            PRINT*,'The supercell cannot have negative length.'
            PRINT*,'Here are your Lx, Ly and Lz as defined in input/dft.in :',grid%length
            STOP "in module_grid> init_grid"
        end if
        grid%l(1:3) = grid%length(1:3)
        grid%lx = grid%length(1)
        grid%ly = grid%length(2)
        grid%lz = grid%length(3)

        grid%n_nodes = getinput%int3( "boxnod" , defaultvalue= nint(grid%length/0.3_dp), assert=">0" )
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
        print*, "[grid]====="
        print*, "Box Length :", grid%length
        print*, "nodes      :", grid%n_nodes
        print*, "dx, dy, dz :", grid%dl
        print*, "[/grid]===="
        print*,

        call tabulate_kx_ky_kz
        grid%isinitiated = .true.

    end subroutine init

    ! subroutine init_grid
    !     use module_input, only: getinput
    !     use module_quadrature, only: prepare_quadrature => init
    !     implicit none
    !     if (grid%isinitiated) then
    !         print*, "Dans init_grid, c'est bizarre. On veut initialiser le type derivé grid mais il semble deja initialisé"
    !         stop "dans module_grid/init_grid "
    !     end if
    !     grid%molRotSymOrder = getinput%int('molRotSymOrder', defaultvalue=1, assert=">0") !Get the order of the main symmetry axis of the solvent
    !     grid%length = getinput%dp3( "lxlylz" , defaultvalue= grid%length )
    !     if (ANY( grid%length  <= 0._dp ) ) THEN
    !         PRINT*,'The supercell cannot have negative length.'
    !         PRINT*,'Here are your Lx, Ly and Lz as defined in input/dft.in :',grid%length
    !         STOP "in module_grid> init_grid"
    !     end if
    !     grid%l(1:3) = grid%length(1:3)
    !     grid%lx = grid%length(1)
    !     grid%ly = grid%length(2)
    !     grid%lz = grid%length(3)
    !
    !     grid%n_nodes = getinput%int3( "nxnynz" , defaultvalue= nint(grid%length/0.3_dp) )
    !     if ( any(grid%n_nodes <= 0) ) then
    !         print*, 'The space is divided into grid nodes. For each direction, you ask', grid%n_nodes,'node.'
    !         error stop
    !     end if
    !     grid%n = grid%n_nodes
    !     grid%nx = grid%n(1)
    !     grid%ny = grid%n(2)
    !     grid%nz = grid%n(3)
    !
    !     grid%dl = grid%length / real(grid%n,dp) !
    !     grid%dx = grid%dl(1)
    !     grid%dy = grid%dl(2)
    !     grid%dz = grid%dl(3)
    !
    !     grid%v = product(grid%length)
    !     grid%dv = product(grid%dl)
    !
    !     ! We now have a full description of the space grid
    !     print*,
    !     print*, "[grid]====="
    !     print*, "Box Length :", grid%length
    !     print*, "nodes      :", grid%n_nodes
    !     print*, "dx, dy, dz :", grid%dl
    !     print*, "[/grid]===="
    !     print*,
    !
    !     call prepare_quadrature
    !     call tabulate_kx_ky_kz
    !     grid%isinitiated = .true.
    !
    ! end subroutine init_grid


    SUBROUTINE tabulate_kx_ky_kz
        integer :: l
        integer :: nx, ny, nz
        nx = grid%nx
        ny = grid%ny
        nz = grid%nz
        allocate ( grid%kx(nx/2+1), source=0._dp)
        allocate ( grid%ky(ny)    , source=0._dp)
        allocate ( grid%kz(nz)    , source=0._dp)
        do concurrent ( l=1:nx/2+1 )
            grid%kx(l) = kproj(1,l)
        end do
        do concurrent ( l=1:ny )
            grid%ky(l) = kproj(2,l)
        end do
        do concurrent ( l=1:nz )
            grid%kz(l) = kproj(3,l)
        end do
    END SUBROUTINE


    PURE FUNCTION k2 (l,m,n) ! UTILE ????
        integer, INTENT(IN) :: l,m,n
        REAL(dp) :: k2
        k2 = grid%kx(l)**2 + grid%ky(m)**2 + grid%kz(n)**2
    END FUNCTION k2


    PURE FUNCTION norm_k (l,m,n)
        integer, intent(in) :: l,m,n
        real(dp) :: norm_k
        norm_k = sqrt(k2(l,m,n))
    END FUNCTION norm_k


    PURE FUNCTION kproj (dir,l)
        ! note the special ordering for negative values. See FFTW (FFTW3) documentation
        ! http://www.fftw.org/doc/Real_002ddata-DFT-Array-Format.html#Real_002ddata-DFT-Array-Format
        ! http://www.fftw.org/doc/The-1d-Discrete-Fourier-Transform-_0028DFT_0029.html#The-1d-Discrete-Fourier-Transform-_0028DFT_0029
        integer, INTENT(IN) :: dir, l ! dir is 1 for x, 2 for y, 3 for z
        REAL(dp) :: kproj
        integer :: m1
        real(dp), parameter :: twopi=2._dp*acos(-1._dp)
        IF ( l <= grid%n_nodes(dir)/2 ) THEN
            m1 = l - 1
        ELSE
            m1 = l - 1 - grid%n_nodes(dir)
        END IF
        kproj = twopi/grid%length(dir)*REAL(m1,dp)
    END FUNCTION


    PURE FUNCTION kvec (l,m,n)
        integer, INTENT(IN) :: l,m,n
        REAL(dp), DIMENSION(3) :: kvec
        kvec(1:3) = [ kproj(1,l), kproj(2,l), kproj(3,l) ]
    END FUNCTION kvec


    PURE FUNCTION timesExpPrefactork2 (array3D, prefactor)
        COMPLEX(dp), DIMENSION(:,:,:), INTENT(IN) :: array3D
        COMPLEX(dp), DIMENSION(SIZE(array3D,1),SIZE(array3D,2),SIZE(array3D,3)) :: timesExpPrefactork2
        REAL(dp), INTENT(IN) :: prefactor
        integer :: i,j,k,imax,jmax,kmax
        imax = SIZE(array3D,1)
        jmax = SIZE(array3D,2)
        kmax = SIZE(array3D,3)
        DO CONCURRENT ( i=1:imax, j=1:jmax, k=1:kmax )
            timesExpPrefactork2 (i,j,k) = array3D (i,j,k) * EXP( prefactor* k2 (i,j,k) )
        END DO
    END FUNCTION

end module module_grid

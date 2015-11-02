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
        integer :: molrotsymorder, mmax, ntheta, nphi, npsi, no, np
        real(dp) :: dphi, dpsi
        real(dp), allocatable :: theta(:), phi(:), psi(:), wtheta(:), wphi(:), wpsi(:), w(:), thetaofitheta(:)
        integer, allocatable :: indo(:,:,:) ! table of index of orientations
        real(dp), allocatable, dimension(:) :: rotxx, rotxy, rotxz, rotyx, rotyy, rotyz, rotzx, rotzy, rotzz
        real(dp), allocatable, dimension(:) :: OMx, OMy, OMz
    contains
        procedure, nopass :: init, integrate_over_orientations
    end type somegrid
    type(somegrid), protected :: grid
    real(dp), parameter, private :: eightpisq=8._dp*acos(-1._dp)**2
    real(dp), parameter, private :: quadrature_norm=eightpisq

    public :: norm_k, timesExpPrefactork2, k2, mean_over_orientations, grid


contains

    subroutine init
        use module_input, only: getinput
        implicit none
        real(dp), parameter :: twopi=2._dp*acos(-1._dp)
        integer :: io, m, mup, mu

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
        print*, "===== Grid ====="
        print*, "   lx, ly, lz :", real(grid%length)
        print*, "   nx, ny, nz :", grid%n_nodes
        print*, "   dx, dy, dz :", real(grid%dl)

        grid%mmax = getinput%int("mmax", defaultvalue=0, assert=">=0")
        grid%ntheta = grid%mmax+1
        grid%nphi = 2*grid%mmax+1
        grid%npsi = 2*(grid%mmax/grid%molrotsymorder)+1
        grid%no = grid%ntheta*grid%nphi*grid%npsi   ! nombez d'orientations dans la representation Euler
        grid%dphi = twopi/real(grid%nphi,dp)
        grid%dpsi = twopi/real(grid%npsi*grid%molrotsymorder,dp)
        grid%np=sum( [( [( [( 1 ,mu=0,m,grid%molrotsymorder)], mup=-m,m)], m=0,grid%mmax)] ) ! number of projections
        allocate( grid%theta(grid%no) , source=0._dp) ! io => theta
        allocate( grid%phi(grid%no) , source=0._dp)
        allocate( grid%psi(grid%no) , source=0._dp)
        allocate( grid%wtheta(grid%no) , source=0._dp)
        allocate( grid%wphi(grid%no) , source=0._dp)
        allocate( grid%wpsi(grid%no) , source=0._dp)
        allocate( grid%w(grid%no) , source=0._dp)
        allocate( grid%thetaofitheta(grid%ntheta), source=0._dp) ! itheta => theta
        block
            real(dp) :: thetaofitheta(grid%ntheta), wthetaofitheta(grid%ntheta)
            integer :: err, io, i, itheta, iphi, ipsi
            real(dp) :: phiofiphi(grid%nphi), psiofipsi(grid%npsi)
            real(dp) :: wphiofiphi(grid%nphi)
            real(dp) :: wpsiofipsi(grid%npsi)
            wphiofiphi = 1._dp/real(grid%nphi,dp)
            wpsiofipsi = 1._dp/real(grid%npsi,dp)
            psiofipsi = [(   real(i-1,dp)*grid%dpsi   , i=1,grid%npsi )]
            phiofiphi = [(   real(i-1,dp)*grid%dphi   , i=1,grid%nphi )]
            call gauss_legendre( grid%ntheta, thetaofitheta, wthetaofitheta, err)
            if (err /= 0) error stop "problem in gauss_legendre"
            thetaofitheta = acos(thetaofitheta)
            allocate( grid%indo(grid%ntheta,grid%nphi,grid%npsi), source=-huge(1) )
            io = 0
            do itheta = 1, grid%ntheta
                grid%thetaofitheta(itheta) = thetaofitheta(itheta)
                do iphi = 1, grid%nphi
                    do ipsi = 1, grid%npsi
                        io = io+1
                        grid%theta(io) = thetaofitheta(itheta)
                        grid%phi(io) = phiofiphi(iphi)
                        grid%psi(io) = psiofipsi(ipsi)
                        grid%indo(itheta,iphi,ipsi) = io
                        grid%wtheta(io) = wthetaofitheta(itheta)
                        grid%wphi(io) = wphiofiphi(iphi)
                        grid%wpsi(io) = wpsiofipsi(ipsi)
                        grid%w(io) = grid%wtheta(io) * grid%wphi(io) * grid%wpsi(io) *quadrature_norm
                    end do
                end do
            end do
        end block


        if (abs(sum(grid%w)-quadrature_norm)>1000*epsilon(1._dp)) then
            print*, "In module_grid.f90, I am checking the normalization of the angular quadrature"
            print*, "The sum of all quadrature weights is", sum(grid%w)
            print*, "It should be ", quadrature_norm
            print*, "deviation to reference value is", abs(sum(grid%w)-quadrature_norm)
            print*, "I fixed the error acceptance to", 1000*epsilon(1._dp)
            error stop
        end if

        open(56, file="output/su3-roots-weights.dat")
        write(56,*) "theta, weight_theta, phi, weight_phi, psi, weight_psi, weight_total_of_theta_phi_psi"
        do io=1,grid%no
            write(56,*) grid%theta(io), grid%wtheta(io), grid%phi(io), grid%wphi(io), grid%psi(io), grid%wpsi(io), grid%w(io)
        end do
        close(56)

        print*, "   mmax =",int(grid%mmax,1),"and molrotsymorder =",int(grid%molrotsymorder,1)
        print*, "=>",grid%no,"orientations and",grid%np,"projections of the density"
        print*, "   θ φ ψ:",int([grid%ntheta,grid%nphi,grid%npsi],1)
        print*, "===== Grid ====="

        allocate( grid%rotxx(grid%no) , grid%rotxy(grid%no), grid%rotxz(grid%no), source=0._dp)
        allocate( grid%rotyx(grid%no) , grid%rotyy(grid%no), grid%rotyz(grid%no), source=0._dp)
        allocate( grid%rotzx(grid%no) , grid%rotzy(grid%no), grid%rotzz(grid%no), source=0._dp)
        allocate( grid%OMx(grid%no), grid%OMy(grid%no), grid%OMz(grid%no) , source=0._dp)

        select case (grid%no)
        case (1)
            grid%Rotxx = 1.0_dp ; grid%Rotxy = 0.0_dp ; grid%Rotxz = 0.0_dp
            grid%Rotyx = 0.0_dp ; grid%Rotyy = 1.0_dp ; grid%Rotyz = 0.0_dp
            grid%Rotzx = 0.0_dp ; grid%Rotzy = 0.0_dp ; grid%Rotzz = 1.0_dp
            grid%OMx = 0._dp ; grid%OMy = 0._dp ; grid%OMz = 1._dp ! ATTENTION "completely" arbitrary decision to put OMEGA along z.
        case default
            block
                integer :: itheta, io
                real(dp) :: cos_theta, sin_theta, cos_phi, sin_phi, cos_psi, sin_psi
                do io = 1, grid%no
                    cos_theta = cos(grid%theta(io))
                    sin_theta = sin(grid%theta(io))
                    cos_phi = cos(grid%phi(io))
                    sin_phi = sin(grid%phi(io))
                    grid%OMx(io) = sin_theta * cos_phi
                    grid%OMy(io) = sin_theta * sin_phi
                    grid%OMz(io) = cos_theta
                    cos_psi = cos(grid%psi(io))
                    sin_psi = sin(grid%psi(io))
                    grid%rotxx(io) = cos_theta*cos_phi*cos_psi-sin_phi*sin_psi
                    grid%rotxy(io) = -cos_theta*cos_phi*sin_psi-sin_phi*cos_psi
                    grid%rotxz(io) = sin_theta*cos_phi
                    grid%rotyx(io) = cos_theta*sin_phi*cos_psi+cos_phi*sin_psi
                    grid%rotyy(io) = -cos_theta*sin_phi*sin_psi+cos_phi*cos_psi
                    grid%rotyz(io) = sin_theta*sin_phi
                    grid%rotzx(io) = -sin_theta*cos_psi
                    grid%rotzy(io) = sin_theta*sin_psi
                    grid%rotzz(io) = cos_theta
                end do
            end block
        end select

        call tabulate_kx_ky_kz
        grid%isinitiated = .true.
    end subroutine init

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
    END SUBROUTINE tabulate_kx_ky_kz

    !> This SUBROUTINE do the legendre integration over all orientations of any array (density, vext, ...)
    SUBROUTINE mean_over_orientations ( arrayin , arrayout )
        use precision_kinds, only: dp
        implicit none
        real(dp), intent(in) :: arrayin(:,:,:,:)  ! x, y, z, orientation
        real(dp), intent(out) :: arrayout(:,:,:)   ! <x, y, z>_orientations
        ! real(dp), dimension(grid%n_nodes(1),grid%n_nodes(2),grid%n_nodes(3),angGrid%n_angles,molRotGrid%n_angles),&
        !         intent(in) :: arrayin ! input array
        ! real(dp), dimension(grid%n_nodes(1),grid%n_nodes(2),grid%n_nodes(3)), intent(out) :: arrayout ! output array
        integer :: io

        arrayout = 0.0_dp
        do io = 1, grid%no
            arrayout = arrayout + arrayin(:,:,:,io)*grid%w(io)
        end do

        ! arrayout = sum( arrayin * grid%w(:))
        ! do n=1,angGrid%n_angles
        !     do p=1,molRotGrid%n_angles
        !         arrayout(:,:,:)=arrayout(:,:,:)+arrayin(:,:,:,n,p)*angGrid%weight(n)*molRotGrid%weight(p)
        !     END DO
        ! END DO

    END SUBROUTINE mean_over_orientations

    subroutine gauss_legendre( n, x, w, exitstatus) ! copy paste from Luc's subroutine Luc74p85
        !
        ! Returns the n roots (x) and associated weights(w) of a gauss legendre quadrature of order n
        ! The roots are the Cos(theta) so that if you need theta, don't forget to acos(x)
        !
        implicit none
        integer, intent(in) :: n
        real(dp), intent(out) :: x(n), w(n)
        integer, optional, intent(out) :: exitstatus
        integer :: m, i, j
        real(dp), parameter :: pi=acos(-1._dp)
        real(dp) :: xi, p1, p2, p3, pp, deltaxi
        exitstatus = 0
        if( n <= 0 ) then
            exitstatus = -1
            return
        else
            m = (n+1)/2
            do i = 1,m          ! on s'interesse au ième zero du polynome pn(x) de legendre
                xi = cos(pi*(i-0.25_dp)/(n+0.5_dp))     ! estimation de départ qu'on va raffiner par nr
                deltaxi = 1._dp
                do while (abs(deltaxi)>epsilon(1._dp))!1.d-15)
                  p1 = 1._dp
                  p2 = 0._dp
                  do j = 1,n
                    p3 = p2
                    p2 = p1
                    p1 = (real(2*j-1,dp)*xi*p2-real(j-1,dp)*p3)/real(j,dp)         ! relation de récurrence entre les pj
                  end do
                  pp = n*(xi*p1-p2)/(xi**2-1._dp)                   ! donne pn' en fonction de pn et pn-1
                  deltaxi = -p1/pp                                  ! nr
                  xi = xi + deltaxi
                end do
                x(i) = xi
                w(i) = 1._dp/((1._dp-xi**2)*pp**2)                  ! poids normalise a 1
                x(n+1-i) = -xi
                w(n+1-i) = w(i)
            end do
        end if
    end subroutine gauss_legendre



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
    END FUNCTION kproj


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
    END FUNCTION timesExpPrefactork2

    subroutine integrate_over_orientations (array, integrated_array)
        use precision_kinds, only: dp
        implicit none
        real(dp), intent(in) :: array(:,:,:,:)
        real(dp), allocatable, intent(out) :: integrated_array(:,:,:)
        integer :: i, il, iu, j, jl, ju, k, kl, ku
        if (size(array,4)/=size(grid%w)) then
            print*, "In module_grid > integrate_over_orientations"
            print*, "You want to integrate an array(:,:,:,:) whose last dimension is not the same as grid%w(:)"
            error stop
        end if
        il=lbound(array,1)
        iu=ubound(array,1)
        jl=lbound(array,2)
        ju=ubound(array,2)
        kl=lbound(array,3)
        ku=ubound(array,3)
        if (.not. allocated(integrated_array)) then
            allocate (integrated_array(il:iu,jl:ju,kl:ku))
        end if
        if (lbound(integrated_array,1)/=il) stop "pb of bounds in integrate_over_orientations"
        if (ubound(integrated_array,1)/=iu) stop "pb of bounds in integrate_over_orientations"
        if (lbound(integrated_array,2)/=jl) stop "pb of bounds in integrate_over_orientations"
        if (ubound(integrated_array,2)/=ju) stop "pb of bounds in integrate_over_orientations"
        if (lbound(integrated_array,3)/=kl) stop "pb of bounds in integrate_over_orientations"
        if (ubound(integrated_array,3)/=ku) stop "pb of bounds in integrate_over_orientations"

        integrated_array=0._dp
        do k=lbound(array,3),ubound(array,3)
            do j=lbound(array,2),ubound(array,2)
                do i=lbound(array,1),ubound(array,1)
                    integrated_array(i,j,k)=sum(array(i,j,k,:)*grid%w(:))
                end do
            end do
        end do
    end subroutine integrate_over_orientations

end module module_grid

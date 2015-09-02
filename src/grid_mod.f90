module grid_mod
    use iso_c_binding, only: dp => c_double
    implicit none
    private
    type, public :: grid
        integer :: nx, ny, nz, ptperangx, ptperangy, ptperangz
        integer :: mmax, molrotsymorder, nproj, no, ntheta, nphi, npsi, dphi, dpsi
        integer, allocatable :: tproj(:,:,:), tmup(:), tmu(:), tm(:)
        real(dp) :: lx, ly, lz, dx, dy, dz, dv, volume
        real(dp), allocatable :: tx(:), ty(:), tz(:)
        real(dp), allocatable :: w(:)
        real(dp), allocatable :: wtheta(:), theta(:), wphi(:), phi(:), wpsi(:), psi(:)
        logical :: is_built = .false.
        integer, allocatable :: tio(:,:,:) ! index of orientation corresponding to itheta,iphi,ipsi
        contains
        procedure :: build
    end type
contains
    subroutine build( gr )
        use input, only: getinput%dp, getinput%int
        implicit none
        class(grid), intent(inout) :: gr
        integer :: err, i, m, mup, mu, itheta, iphi, ipsi, io
        real(dp) :: weight_of_each_phi, weight_of_each_psi
        real(dp), parameter :: twopi = acos(-1._dp)*2._dp
        real(dp), allocatable :: thetaofitheta(:), phiofiphi(:), psiofipsi(:), wthetaofitheta(:), wphiofiphi(:), wpsiofipsi(:)
        if (grid%is_built) return
        stop "I dont want to use this yet. PB IN RESOX INT HERE DP IN ALLOCATE_FROM_INPUT"
        grid%ptperangx = getinput%int("resox", defaultvalue=3)
        grid%ptperangy = getinput%int("resoy", defaultvalue=3)
        grid%ptperangz = getinput%int("resoz", defaultvalue=3)
        print*, "grid%ptperangx=",grid%ptperangx
        print*, "grid%ptperangy=",grid%ptperangy
        print*, "grid%ptperangz=",grid%ptperangz
        stop "ptperang"
        grid%lx = getinput%dp("Lx", defaultvalue=10._dp)
        grid%ly = getinput%dp("Ly", defaultvalue=10._dp)
        grid%lz = getinput%dp("Lz", defaultvalue=10._dp)
        grid%nx = nint(grid%ptperangx * grid%lx)
        grid%ny = nint(grid%ptperangy * grid%ly)
        grid%nz = nint(grid%ptperangz * grid%lz)
        grid%dx = grid%lx/real(grid%nx,dp)
        grid%dy = grid%ly/real(grid%ny,dp)
        grid%dz = grid%lz/real(grid%nz,dp)
        grid%dv = grid%dx * grid%dy * grid%dz
        grid%volume = grid%lx * grid%ly * grid%lz
        allocate(grid%tx(grid%nx), source=[( (i-1)*grid%dx ,i=1,grid%nx )], stat=err)
        if (err /= 0) error stop "grid%tx: Allocation request denied"
        allocate(grid%ty(grid%ny), source=[( (i-1)*grid%dy ,i=1,grid%ny )], stat=err)
        if (err /= 0) error stop "grid%ty: Allocation request denied"
        allocate(grid%tz(grid%nz), source=[( (i-1)*grid%dz ,i=1,grid%nz )], stat=err)
        if (err /= 0) error stop "grid%tz: Allocation request denied"
        grid%molrotsymorder = getinput%int("molrotsymorder", defaultvalue=1)
        grid%mmax = getinput%int("mmax", defaultvalue=4)
        grid%ntheta = grid%mmax+1
        grid%nphi = 2*grid%mmax+1
        grid%npsi = 2*(grid%mmax/grid%molrotsymorder)+1
        grid%no = grid%ntheta*grid%nphi*grid%npsi
        ! grid%nproj = SUM( [( [( [( 1 ,mu=0,m,molrotsymorder)], mup=-m,m)], m=0,mmax)] )
        ! allocate (grid%tproj(0:mmax,-mmax:mmax,0:mmax), source=-huge(kind(iproj)) ) ! alpha(m,mup,mu)
        ! allocate (grid%tm(grid%nproj), source=-huge(1))
        ! allocate (grid%tmup(grid%nproj), source=-huge(1))
        ! allocate (grid%tmu(grid%nproj), source=-huge(1))
        allocate (grid%w(grid%no), source=1._dp, stat=err)
        if (err /= 0) error stop "grid%w: Allocation request denied"
        !
        ! Build the Euler representation suitable for Gauss-Legendre quadrature
        !
        grid%dphi = twopi/real(grid%nphi,dp)
        grid%dpsi = twopi/real(grid%npsi*grid%molrotsymorder,dp)
        allocate( thetaofitheta(grid%ntheta) ,source=huge(1._dp) )
        allocate( phiofiphi(grid%nphi)   ,source=[(   real(i-1,dp)*grid%dphi   , i=1,grid%nphi )] )
        allocate( psiofipsi(grid%npsi)   ,source=[(   real(i-1,dp)*grid%dpsi   , i=1,grid%npsi )] )
        allocate( wthetaofitheta(grid%ntheta))
        allocate( wphiofiphi(grid%nphi) ,source=1._dp/grid%nphi)
        allocate( wpsiofipsi(grid%npsi) ,source=1._dp/grid%npsi)
        allocate( grid%theta(grid%no) )
        allocate( grid%phi(grid%no) )
        allocate( grid%psi(grid%no) )
        allocate( grid%wtheta(grid%no))
        allocate( grid%wphi(grid%no))
        allocate( grid%wpsi(grid%no))
        call gauss_legendre( grid%ntheta, thetaofitheta, wthetaofitheta, err )
        if (err /= 0) error stop "problem in gauss_legendre"
        thetaofitheta = acos(thetaofitheta)
        allocate( grid%tio(grid%ntheta,grid%nphi,grid%npsi), source=-huge(1) )
        io = 0
        do itheta = 1, grid%ntheta
            do iphi = 1, grid%nphi
                do ipsi = 1, grid%npsi
                    io = io+1
                    grid%theta(io) = thetaofitheta(itheta)
                    grid%phi(io) = phiofiphi(iphi)
                    grid%psi(io) = psiofipsi(ipsi)
                    grid%tio(itheta,iphi,ipsi) = io
                    grid%wtheta(io) = wthetaofitheta(itheta)
                    grid%wphi(io) = wphiofiphi(iphi)
                    grid%wpsi(io) = wpsiofipsi(ipsi)
                    grid%w(io) = grid%wtheta(io) * grid%wphi(io) * grid%wpsi(io)
                end do
            end do
        end do
        grid%is_built = .true.
    end subroutine
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
        real(dp), PARAMETER :: pi=acos(-1._dp)
        real(dp) :: xi, p1, p2, p3, pp, deltaxi
        exitstatus = 0
        if( n <= 0 ) then
            exitstatus = -1
            return
        else
            m = (n+1)/2
            do i = 1,m          ! on s'interesse au ième zero du polynome pn(x) de legendre
                xi = cos(pi*(i-0.25)/(n+0.5))     ! estimation de départ qu'on va raffiner par nr
                deltaxi = 1
                do while (abs(deltaxi)>1.d-13)
                  p1 = 1.
                  p2 = 0.
                  do j = 1,n
                    p3 = p2
                    p2 = p1
                    p1 = ((2*j-1.)*xi*p2-(j-1.)*p3)/j         ! relation de récurrence entre les pj
                  end do
                  pp = n*(xi*p1-p2)/(xi**2-1.)                   ! donne pn' en fonction de pn et pn-1
                  deltaxi = -p1/pp                                  ! nr
                  xi = xi + deltaxi
                end do
                x(i) = xi
                w(i) = 1./((1.-xi**2)*pp**2)                  ! poids normalise a 1
                x(n+1-i) = -xi
                w(n+1-i) = w(i)
            end do
        end if
    end subroutine
end module

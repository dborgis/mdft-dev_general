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
        use input, only: input_dp, input_int
        implicit none
        class(grid), intent(inout) :: gr
        integer :: err, i, m, mup, mu, itheta, iphi, ipsi, io
        real(dp) :: weight_of_each_phi, weight_of_each_psi
        real(dp), parameter :: twopi = acos(-1._dp)*2._dp
        real(dp), allocatable :: thetaofitheta(:), phiofiphi(:), psiofipsi(:), wthetaofitheta(:), wphiofiphi(:), wpsiofipsi(:)
        if (gr%is_built) return
        stop "I dont want to use this yet. PB IN RESOX INT HERE DP IN ALLOCATE_FROM_INPUT"
        gr%ptperangx = input_int("resox", defaultvalue=3)
        gr%ptperangy = input_int("resoy", defaultvalue=3)
        gr%ptperangz = input_int("resoz", defaultvalue=3)
        print*, "gr%ptperangx=",gr%ptperangx
        print*, "gr%ptperangy=",gr%ptperangy
        print*, "gr%ptperangz=",gr%ptperangz
        stop "ptperang"
        gr%lx = input_dp("Lx", defaultvalue=10._dp)
        gr%ly = input_dp("Ly", defaultvalue=10._dp)
        gr%lz = input_dp("Lz", defaultvalue=10._dp)
        gr%nx = nint(gr%ptperangx * gr%lx)
        gr%ny = nint(gr%ptperangy * gr%ly)
        gr%nz = nint(gr%ptperangz * gr%lz)
        gr%dx = gr%lx/real(gr%nx,dp)
        gr%dy = gr%ly/real(gr%ny,dp)
        gr%dz = gr%lz/real(gr%nz,dp)
        gr%dv = gr%dx * gr%dy * gr%dz
        gr%volume = gr%lx * gr%ly * gr%lz
        allocate(gr%tx(gr%nx), source=[( (i-1)*gr%dx ,i=1,gr%nx )], stat=err)
        if (err /= 0) error stop "gr%tx: Allocation request denied"
        allocate(gr%ty(gr%ny), source=[( (i-1)*gr%dy ,i=1,gr%ny )], stat=err)
        if (err /= 0) error stop "gr%ty: Allocation request denied"
        allocate(gr%tz(gr%nz), source=[( (i-1)*gr%dz ,i=1,gr%nz )], stat=err)
        if (err /= 0) error stop "gr%tz: Allocation request denied"
        gr%molrotsymorder = input_int("molrotsymorder", defaultvalue=1)
        gr%mmax = input_int("mmax", defaultvalue=4)
        gr%ntheta = gr%mmax+1
        gr%nphi = 2*gr%mmax+1
        gr%npsi = 2*(gr%mmax/gr%molrotsymorder)+1
        gr%no = gr%ntheta*gr%nphi*gr%npsi
        ! gr%nproj = SUM( [( [( [( 1 ,mu=0,m,molrotsymorder)], mup=-m,m)], m=0,mmax)] )
        ! allocate (gr%tproj(0:mmax,-mmax:mmax,0:mmax), source=-huge(kind(iproj)) ) ! alpha(m,mup,mu)
        ! allocate (gr%tm(gr%nproj), source=-huge(1))
        ! allocate (gr%tmup(gr%nproj), source=-huge(1))
        ! allocate (gr%tmu(gr%nproj), source=-huge(1))
        allocate (gr%w(gr%no), source=1._dp, stat=err)
        if (err /= 0) error stop "gr%w: Allocation request denied"
        !
        ! Build the Euler representation suitable for Gauss-Legendre quadrature
        !
        gr%dphi = twopi/real(gr%nphi,dp)
        gr%dpsi = twopi/real(gr%npsi*gr%molrotsymorder,dp)
        allocate( thetaofitheta(gr%ntheta) ,source=huge(1._dp) )
        allocate( phiofiphi(gr%nphi)   ,source=[(   real(i-1,dp)*gr%dphi   , i=1,gr%nphi )] )
        allocate( psiofipsi(gr%npsi)   ,source=[(   real(i-1,dp)*gr%dpsi   , i=1,gr%npsi )] )
        allocate( wthetaofitheta(gr%ntheta))
        allocate( wphiofiphi(gr%nphi) ,source=1._dp/gr%nphi)
        allocate( wpsiofipsi(gr%npsi) ,source=1._dp/gr%npsi)
        allocate( gr%theta(gr%no) )
        allocate( gr%phi(gr%no) )
        allocate( gr%psi(gr%no) )
        allocate( gr%wtheta(gr%no))
        allocate( gr%wphi(gr%no))
        allocate( gr%wpsi(gr%no))
        call gauss_legendre( gr%ntheta, thetaofitheta, wthetaofitheta, err )
        if (err /= 0) error stop "problem in gauss_legendre"
        thetaofitheta = acos(thetaofitheta)
        allocate( gr%tio(gr%ntheta,gr%nphi,gr%npsi), source=-huge(1) )
        io = 0
        do itheta = 1, gr%ntheta
            do iphi = 1, gr%nphi
                do ipsi = 1, gr%npsi
                    io = io+1
                    gr%theta(io) = thetaofitheta(itheta)
                    gr%phi(io) = phiofiphi(iphi)
                    gr%psi(io) = psiofipsi(ipsi)
                    gr%tio(itheta,iphi,ipsi) = io
                    gr%wtheta(io) = wthetaofitheta(itheta)
                    gr%wphi(io) = wphiofiphi(iphi)
                    gr%wpsi(io) = wpsiofipsi(ipsi)
                    gr%w(io) = gr%wtheta(io) * gr%wphi(io) * gr%wpsi(io)
                end do
            end do
        end do
        gr%is_built = .true.
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

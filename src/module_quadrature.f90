! Module for numerical integration
module quadrature

  use precision_kinds ,only: dp, i2b
  use constants       ,only: pi, twopi, fourpi, zero, epsdp
  use input           ,only: getinput
  use module_grid, only: grid

  implicit none

  real(dp), allocatable, dimension(:), private :: x_leb, y_leb , z_leb
  real(dp), allocatable, dimension(:), public :: Omx , Omy , Omz  ! unit vector for Euler angle Omega in lab frame
  real(dp), allocatable, dimension(:,:), public :: Rotxx, Rotxy, Rotxz, Rotyx, Rotyy, Rotyz, Rotzx, Rotzy, Rotzz

  type angularGrid
      integer(i2b) :: n_angles, N, o
      character(80) :: name
      real(dp), allocatable, dimension(:) :: weight, root, w, x, y, z
  end type

  type (angularGrid), public :: angGrid ! angular grid
  type (angularGrid), public :: molRotGrid ! rotation of molecule around its main axis, e.g., around C2v axis for H2O.
  type (angularGrid), public :: qsu2  ! quadrature for integration of a function over the unit sphere, often called s2 or su2

  type integrationScheme
      character(80) :: name
      integer(i2b) :: order
      real(dp), allocatable, dimension(:) :: weight, root
  end type

  type (integrationScheme), public :: intScheme
  integer(i2b), public :: molrotsymorder

  contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init
    use mathematica, only: chop
    implicit none
    integer :: N, P, i, j, no, itheta, iphi, ipsi

    print*,
    print*,"==="
    print*,"grid pour quadrature angulaire"
    grid%molrotsymorder = getinput%int("molrotsymorder", defaultvalue=1)
    molrotsymorder = grid%molrotsymorder
    grid%mmax = getinput%int("mmax", defaultvalue=4)
    grid%ntheta = grid%mmax+1
    grid%nphi = 2*grid%mmax+1
    grid%npsi = 2*(grid%mmax/grid%molrotsymorder)+1
    grid%no = grid%ntheta*grid%nphi*grid%npsi   ! nombez d'orientations dans la representation Euler
    grid%dphi = twopi/real(grid%nphi,dp)
    grid%dpsi = twopi/real(grid%npsi*grid%molrotsymorder,dp)
    n = grid%no
    allocate( grid%theta(n) , source=0._dp)
    allocate( grid%phi(n) , source=0._dp)
    allocate( grid%psi(n) , source=0._dp)
    allocate( grid%wtheta(n) , source=0._dp)
    allocate( grid%wphi(n) , source=0._dp)
    allocate( grid%wpsi(n) , source=0._dp)
    allocate( grid%w(n) , source=0._dp)
    block
        real(dp) :: thetaofitheta(grid%ntheta), wthetaofitheta(grid%ntheta)
        integer :: err, io, itheta, iphi, ipsi
        real(dp) :: phiofiphi(grid%nphi), psiofipsi(grid%npsi)
        real(dp) :: wphiofiphi(grid%nphi)
        real(dp) :: wpsiofipsi(grid%npsi)
        wphiofiphi = 1._dp/grid%nphi
        wpsiofipsi = 1._dp/grid%npsi
        psiofipsi = [(   real(i-1,dp)*grid%dpsi   , i=1,grid%npsi )]
        phiofiphi = [(   real(i-1,dp)*grid%dphi   , i=1,grid%nphi )]
        call gauss_legendre( grid%ntheta, thetaofitheta, wthetaofitheta, err)
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
    end block
    print*, "mmax =",grid%mmax
    print*, "number of orientations",grid%no
    print*, "number of theta =", grid%ntheta
    print*, "number of phi=", grid%nphi
    print*, "molrotsymorder =", grid%molrotsymorder
    print*, "number of psi=", grid%npsi
    print*, "sum of all weights =", sum(grid%w)
    allocate( grid%rotxx(n) , grid%rotxy(n), grid%rotxz(n), source=0._dp)
    allocate( grid%rotyx(n) , grid%rotyy(n), grid%rotyz(n), source=0._dp)
    allocate( grid%rotzx(n) , grid%rotzy(n), grid%rotzz(n), source=0._dp)
    allocate( grid%OMx(n), grid%OMy(n), grid%OMz(n) , source=0._dp)
    select case (grid%no)
    case (1)
        grid%Rotxx = 1.0_dp ; grid%Rotxy = 0.0_dp ; grid%Rotxz = 0.0_dp
        grid%Rotyx = 0.0_dp ; grid%Rotyy = 1.0_dp ; grid%Rotyz = 0.0_dp
        grid%Rotzx = 0.0_dp ; grid%Rotzy = 0.0_dp ; grid%Rotzz = 1.0_dp
        grid%OMx = 0._dp ; grid%OMy = 0._dp ; grid%OMz = 1._dp ! ATTENTION completely arbitrary. We decide to put it along z.
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
    print*,"=== /grid"
    print*,


    return





    ! quadrature for psi. For now a uniform grid over psi between 0 and 2pi.
    molRotGrid%n_angles = getinput%int('nb_psi', defaultvalue=4 ) ! number of nodes for the grid over psi
    N = molrotgrid%n_angles
    print*,"Quadrature for psi : uniform"
    print*,"Number of nodes for psi :",N * molrotsymorder
    print*,"Thanks to symetries, we will only use this number of nodes for psi :",N
    if ( N < 1 ) stop "in module_quadrature, nb_psi, readen from input/dft.in is unphysical. critical stop."
    allocate( molRotGrid%root(N), source=0._dp)
    do concurrent (i=1:N) ! equidistant repartition between 0 and 2Pi
      molRotGrid%root(i) = real(i-1,dp)*twopi/real(N*molrotsymorder,dp)
    end do
    allocate( molrotgrid%weight(N) , source=  twopi/real(N*molrotsymorder ,dp)   ) ! homogeneous grid between 0 and 2pi. all weights equa l

    ! NOW SU2
    ! integration scheme
    qsu2%o = getinput%int('order_of_quadrature', defaultvalue=3)

    qsu2%name = trim(adjustl(   getinput%char('quadrature')   ))
    select case (qsu2%name(1:2))
    case ('LE')
      print*,"Quadrature for SU2 : Lebedev"
    case ('SD')
      print*,"Quadrature for SU2 : Spherical Design"
    case default
      print*, "To ingrate over the unit sphere, you want a quadrature called ", qsu2%name,"lala"
      stop "It is not implemented"
    end select
    print*,"Order of quadrature for SU2 :",qsu2%o

    ! number of quadrature nodes
    qsu2%n = huge(1) ! start by non-physical value
    select case (qsu2%name)
    case ("LE")
      select case (qsu2%o)
      case(0); qsu2%n=1
      case(3); qsu2%n=6
      case(5); qsu2%n=14
      case(7); qsu2%n=26
      case(9); qsu2%n=38
      case default
        print*, "Quadrature order",qsu2%o,"is not implemented with Lebedev."
        error stop
      end select
    case ("SD")
      select case (qsu2%o)
      case(0); qsu2%n = 1
      case(1); qsu2%n = 2
      case(2); qsu2%n = 4
      case(3); qsu2%n = 6
      case(4); qsu2%n = 10
      case(5); qsu2%n = 12
      case(6); qsu2%n = 18
      case(7); qsu2%n = 22
      case(8); qsu2%n = 28
      case(9); qsu2%n = 32
      case(10); qsu2%n = 42
      case(11); qsu2%n = 48
      case(12); qsu2%n = 58
      case(13); qsu2%n = 64
      case(14); qsu2%n = 72
      case(15); qsu2%n = 82
      case(16); qsu2%n = 98
      case(17); qsu2%n = 104
      case(18); qsu2%n = 122
      case(19); qsu2%n = 130
      case(20); qsu2%n = 148
      case(21); qsu2%n = 156
      case(22); qsu2%n = 178
      case(23); qsu2%n = 186
      case(24); qsu2%n = 210
      case(25); qsu2%n = 220
      case(26); qsu2%n = 244
      case(27); qsu2%n = 254
      case(28); qsu2%n = 282
      case(29); qsu2%n = 292
      case(30); qsu2%n = 322
      case(32); qsu2%n = 364
      case(34); qsu2%n = 410
      case(35); qsu2%n = 422
      case(36); qsu2%n = 458
      case(37); qsu2%n = 472
      case(38); qsu2%n = 508
      case(39); qsu2%n = 522
      case(44); qsu2%n = 672
      case default
        print*, "You ask for a quadrature over unit sphere of order", qsu2%o," that is not compatible with Gauss type spherical designs"
        stop
      end select
    end select
    print*, "Number of grid nodes for SU2 :", qsu2%n

    N = qsu2%n
    P = molrotgrid%n_angles

    ! Rotij is the rotation matrix in Euler coordinates. See manual files.
    allocate( Rotxx(N,P), Rotxy(N,P), Rotxz(N,P), Rotyx(N,P), Rotyy(N,P), Rotyz(N,P), Rotzx(N,P), Rotzy(N,P), Rotzz(N,P) )

    ! OM is the dipole vector, or more exactly the z axis of the molecular frame
    allocate( OMx(N), OMy(N), OMz(N) )

    ! quadrature nodes
    allocate( qsu2%x(N), qsu2%y(N), qsu2%z(N), qsu2%w(N) )
    qsu2%x = huge(1._dp)
    qsu2%y = huge(1._dp)
    qsu2%z = huge(1._dp)
    qsu2%w = huge(1._dp)

    select case( qsu2%name )
    case('LE') ! LEBEDEV nodes
      select case (qsu2%o)
      case (0)
        qsu2%x( 1) = 0._dp; qsu2%y(1)=0._dp; qsu2%z(1)=1._dp; qsu2%w(1)=fourpi
      case (3)
        qsu2%x( 1) =  1.0_dp ; qsu2%y( 1) =  0.0_dp ; qsu2%z( 1) =  0.0_dp ; qsu2%w(1)= 2.0943951023931953_dp
        qsu2%x( 2) = -1.0_dp ; qsu2%y( 2) =  0.0_dp ; qsu2%z( 2) =  0.0_dp ; qsu2%w(2)= 2.0943951023931953_dp
        qsu2%x( 3) =  0.0_dp ; qsu2%y( 3) =  1.0_dp ; qsu2%z( 3) =  0.0_dp ; qsu2%w(3)= 2.0943951023931953_dp
        qsu2%x( 4) =  0.0_dp ; qsu2%y( 4) = -1.0_dp ; qsu2%z( 4) =  0.0_dp ; qsu2%w(4)= 2.0943951023931953_dp
        qsu2%x( 5) =  0.0_dp ; qsu2%y( 5) =  0.0_dp ; qsu2%z( 5) =  1.0_dp ; qsu2%w(5)= 2.0943951023931953_dp
        qsu2%x( 6) =  0.0_dp ; qsu2%y( 6) =  0.0_dp ; qsu2%z( 6) = -1.0_dp ; qsu2%w(6)= 2.0943951023931953_dp
      case (5)
        qsu2%x( 1) =   1.0000000000000000_dp ; qsu2%y( 1) =   0.0000000000000000_dp ; qsu2%z( 1) =   0.0000000000000000_dp ; qsu2%w( 1) = 0.83775804095727813_dp
        qsu2%x( 2) =  -1.0000000000000000_dp ; qsu2%y( 2) =   0.0000000000000000_dp ; qsu2%z( 2) =   0.0000000000000000_dp ; qsu2%w( 2) = 0.83775804095727813_dp
        qsu2%x( 3) =   0.0000000000000000_dp ; qsu2%y( 3) =   1.0000000000000000_dp ; qsu2%z( 3) =   0.0000000000000000_dp ; qsu2%w( 3) = 0.83775804095727813_dp
        qsu2%x( 4) =   0.0000000000000000_dp ; qsu2%y( 4) =  -1.0000000000000000_dp ; qsu2%z( 4) =   0.0000000000000000_dp ; qsu2%w( 4) = 0.83775804095727813_dp
        qsu2%x( 5) =   0.0000000000000000_dp ; qsu2%y( 5) =   0.0000000000000000_dp ; qsu2%z( 5) =   1.0000000000000000_dp ; qsu2%w( 5) = 0.83775804095727813_dp
        qsu2%x( 6) =   0.0000000000000000_dp ; qsu2%y( 6) =   0.0000000000000000_dp ; qsu2%z( 6) =  -1.0000000000000000_dp ; qsu2%w( 6) = 0.83775804095727813_dp
        qsu2%x( 7) =  0.57735026918962573_dp ; qsu2%y( 7) =  0.57735026918962573_dp ; qsu2%z( 7) =  0.57735026918962573_dp ; qsu2%w( 7) = 0.94247779607693793_dp
        qsu2%x( 8) = -0.57735026918962573_dp ; qsu2%y( 8) =  0.57735026918962573_dp ; qsu2%z( 8) =  0.57735026918962573_dp ; qsu2%w( 8) = 0.94247779607693793_dp
        qsu2%x( 9) =  0.57735026918962573_dp ; qsu2%y( 9) = -0.57735026918962573_dp ; qsu2%z( 9) =  0.57735026918962573_dp ; qsu2%w( 9) = 0.94247779607693793_dp
        qsu2%x(10) = -0.57735026918962573_dp ; qsu2%y(10) = -0.57735026918962573_dp ; qsu2%z(10) =  0.57735026918962573_dp ; qsu2%w(10) = 0.94247779607693793_dp
        qsu2%x(11) =  0.57735026918962573_dp ; qsu2%y(11) =  0.57735026918962573_dp ; qsu2%z(11) = -0.57735026918962573_dp ; qsu2%w(11) = 0.94247779607693793_dp
        qsu2%x(12) = -0.57735026918962573_dp ; qsu2%y(12) =  0.57735026918962573_dp ; qsu2%z(12) = -0.57735026918962573_dp ; qsu2%w(12) = 0.94247779607693793_dp
        qsu2%x(13) =  0.57735026918962573_dp ; qsu2%y(13) = -0.57735026918962573_dp ; qsu2%z(13) = -0.57735026918962573_dp ; qsu2%w(13) = 0.94247779607693793_dp
        qsu2%x(14) = -0.57735026918962573_dp ; qsu2%y(14) = -0.57735026918962573_dp ; qsu2%z(14) = -0.57735026918962573_dp ; qsu2%w(14) = 0.94247779607693793_dp
      case (7)
        qsu2%x( 1) =   1.0000000000000000_dp ; qsu2%y( 1) =   0.0000000000000000_dp ; qsu2%z( 1) =   0.0000000000000000_dp ; qsu2%w( 1) = 0.59839860068377004_dp
        qsu2%x( 2) =  -1.0000000000000000_dp ; qsu2%y( 2) =   0.0000000000000000_dp ; qsu2%z( 2) =   0.0000000000000000_dp ; qsu2%w( 2) = 0.59839860068377004_dp
        qsu2%x( 3) =   0.0000000000000000_dp ; qsu2%y( 3) =   1.0000000000000000_dp ; qsu2%z( 3) =   0.0000000000000000_dp ; qsu2%w( 3) = 0.59839860068377004_dp
        qsu2%x( 4) =   0.0000000000000000_dp ; qsu2%y( 4) =  -1.0000000000000000_dp ; qsu2%z( 4) =   0.0000000000000000_dp ; qsu2%w( 4) = 0.59839860068377004_dp
        qsu2%x( 5) =   0.0000000000000000_dp ; qsu2%y( 5) =   0.0000000000000000_dp ; qsu2%z( 5) =   1.0000000000000000_dp ; qsu2%w( 5) = 0.59839860068377004_dp
        qsu2%x( 6) =   0.0000000000000000_dp ; qsu2%y( 6) =   0.0000000000000000_dp ; qsu2%z( 6) =  -1.0000000000000000_dp ; qsu2%w( 6) = 0.59839860068377004_dp
        qsu2%x( 7) =   0.0000000000000000_dp ; qsu2%y( 7) =  0.70710678118654757_dp ; qsu2%z( 7) =  0.70710678118654757_dp ; qsu2%w( 7) = 0.47871888054701611_dp
        qsu2%x( 8) =   0.0000000000000000_dp ; qsu2%y( 8) = -0.70710678118654757_dp ; qsu2%z( 8) =  0.70710678118654757_dp ; qsu2%w( 8) = 0.47871888054701611_dp
        qsu2%x( 9) =   0.0000000000000000_dp ; qsu2%y( 9) =  0.70710678118654757_dp ; qsu2%z( 9) = -0.70710678118654757_dp ; qsu2%w( 9) = 0.47871888054701611_dp
        qsu2%x(10) =   0.0000000000000000_dp ; qsu2%y(10) = -0.70710678118654757_dp ; qsu2%z(10) = -0.70710678118654757_dp ; qsu2%w(10) = 0.47871888054701611_dp
        qsu2%x(11) =  0.70710678118654757_dp ; qsu2%y(11) =   0.0000000000000000_dp ; qsu2%z(11) =  0.70710678118654757_dp ; qsu2%w(11) = 0.47871888054701611_dp
        qsu2%x(12) = -0.70710678118654757_dp ; qsu2%y(12) =   0.0000000000000000_dp ; qsu2%z(12) =  0.70710678118654757_dp ; qsu2%w(12) = 0.47871888054701611_dp
        qsu2%x(13) =  0.70710678118654757_dp ; qsu2%y(13) =   0.0000000000000000_dp ; qsu2%z(13) = -0.70710678118654757_dp ; qsu2%w(13) = 0.47871888054701611_dp
        qsu2%x(14) = -0.70710678118654757_dp ; qsu2%y(14) =   0.0000000000000000_dp ; qsu2%z(14) = -0.70710678118654757_dp ; qsu2%w(14) = 0.47871888054701611_dp
        qsu2%x(15) =  0.70710678118654757_dp ; qsu2%y(15) =  0.70710678118654757_dp ; qsu2%z(15) =   0.0000000000000000_dp ; qsu2%w(15) = 0.47871888054701611_dp
        qsu2%x(16) = -0.70710678118654757_dp ; qsu2%y(16) =  0.70710678118654757_dp ; qsu2%z(16) =   0.0000000000000000_dp ; qsu2%w(16) = 0.47871888054701611_dp
        qsu2%x(17) =  0.70710678118654757_dp ; qsu2%y(17) = -0.70710678118654757_dp ; qsu2%z(17) =   0.0000000000000000_dp ; qsu2%w(17) = 0.47871888054701611_dp
        qsu2%x(18) = -0.70710678118654757_dp ; qsu2%y(18) = -0.70710678118654757_dp ; qsu2%z(18) =   0.0000000000000000_dp ; qsu2%w(18) = 0.47871888054701611_dp
        qsu2%x(19) =  0.57735026918962573_dp ; qsu2%y(19) =  0.57735026918962573_dp ; qsu2%z(19) =  0.57735026918962573_dp ; qsu2%w(19) = 0.40391905546154477_dp
        qsu2%x(20) = -0.57735026918962573_dp ; qsu2%y(20) =  0.57735026918962573_dp ; qsu2%z(20) =  0.57735026918962573_dp ; qsu2%w(20) = 0.40391905546154477_dp
        qsu2%x(21) =  0.57735026918962573_dp ; qsu2%y(21) = -0.57735026918962573_dp ; qsu2%z(21) =  0.57735026918962573_dp ; qsu2%w(21) = 0.40391905546154477_dp
        qsu2%x(22) = -0.57735026918962573_dp ; qsu2%y(22) = -0.57735026918962573_dp ; qsu2%z(22) =  0.57735026918962573_dp ; qsu2%w(22) = 0.40391905546154477_dp
        qsu2%x(23) =  0.57735026918962573_dp ; qsu2%y(23) =  0.57735026918962573_dp ; qsu2%z(23) = -0.57735026918962573_dp ; qsu2%w(23) = 0.40391905546154477_dp
        qsu2%x(24) = -0.57735026918962573_dp ; qsu2%y(24) =  0.57735026918962573_dp ; qsu2%z(24) = -0.57735026918962573_dp ; qsu2%w(24) = 0.40391905546154477_dp
        qsu2%x(25) =  0.57735026918962573_dp ; qsu2%y(25) = -0.57735026918962573_dp ; qsu2%z(25) = -0.57735026918962573_dp ; qsu2%w(25) = 0.40391905546154477_dp
        qsu2%x(26) = -0.57735026918962573_dp ; qsu2%y(26) = -0.57735026918962573_dp ; qsu2%z(26) = -0.57735026918962573_dp ; qsu2%w(26) = 0.40391905546154477_dp
      case (9)
        qsu2%x( 1) =   1.0000000000000000_dp ; qsu2%y( 1) =   0.0000000000000000_dp ; qsu2%z( 1) =   0.0000000000000000_dp ; qsu2%w( 1) = 0.11967972013675403_dp
        qsu2%x( 2) =  -1.0000000000000000_dp ; qsu2%y( 2) =   0.0000000000000000_dp ; qsu2%z( 2) =   0.0000000000000000_dp ; qsu2%w( 2) = 0.11967972013675403_dp
        qsu2%x( 3) =   0.0000000000000000_dp ; qsu2%y( 3) =   1.0000000000000000_dp ; qsu2%z( 3) =   0.0000000000000000_dp ; qsu2%w( 3) = 0.11967972013675403_dp
        qsu2%x( 4) =   0.0000000000000000_dp ; qsu2%y( 4) =  -1.0000000000000000_dp ; qsu2%z( 4) =   0.0000000000000000_dp ; qsu2%w( 4) = 0.11967972013675403_dp
        qsu2%x( 5) =   0.0000000000000000_dp ; qsu2%y( 5) =   0.0000000000000000_dp ; qsu2%z( 5) =   1.0000000000000000_dp ; qsu2%w( 5) = 0.11967972013675403_dp
        qsu2%x( 6) =   0.0000000000000000_dp ; qsu2%y( 6) =   0.0000000000000000_dp ; qsu2%z( 6) =  -1.0000000000000000_dp ; qsu2%w( 6) = 0.11967972013675403_dp
        qsu2%x( 7) =  0.57735026918962573_dp ; qsu2%y( 7) =  0.57735026918962573_dp ; qsu2%z( 7) =  0.57735026918962573_dp ; qsu2%w( 7) = 0.40391905546154477_dp
        qsu2%x( 8) = -0.57735026918962573_dp ; qsu2%y( 8) =  0.57735026918962573_dp ; qsu2%z( 8) =  0.57735026918962573_dp ; qsu2%w( 8) = 0.40391905546154477_dp
        qsu2%x( 9) =  0.57735026918962573_dp ; qsu2%y( 9) = -0.57735026918962573_dp ; qsu2%z( 9) =  0.57735026918962573_dp ; qsu2%w( 9) = 0.40391905546154477_dp
        qsu2%x(10) = -0.57735026918962573_dp ; qsu2%y(10) = -0.57735026918962573_dp ; qsu2%z(10) =  0.57735026918962573_dp ; qsu2%w(10) = 0.40391905546154477_dp
        qsu2%x(11) =  0.57735026918962573_dp ; qsu2%y(11) =  0.57735026918962573_dp ; qsu2%z(11) = -0.57735026918962573_dp ; qsu2%w(11) = 0.40391905546154477_dp
        qsu2%x(12) = -0.57735026918962573_dp ; qsu2%y(12) =  0.57735026918962573_dp ; qsu2%z(12) = -0.57735026918962573_dp ; qsu2%w(12) = 0.40391905546154477_dp
        qsu2%x(13) =  0.57735026918962573_dp ; qsu2%y(13) = -0.57735026918962573_dp ; qsu2%z(13) = -0.57735026918962573_dp ; qsu2%w(13) = 0.40391905546154477_dp
        qsu2%x(14) = -0.57735026918962573_dp ; qsu2%y(14) = -0.57735026918962573_dp ; qsu2%z(14) = -0.57735026918962573_dp ; qsu2%w(14) = 0.40391905546154477_dp
        qsu2%x(15) =  0.45970084338098310_dp ; qsu2%y(15) =  0.88807383397711526_dp ; qsu2%z(15) =   0.0000000000000000_dp ; qsu2%w(15) = 0.35903916041026207_dp
        qsu2%x(16) = -0.45970084338098310_dp ; qsu2%y(16) =  0.88807383397711526_dp ; qsu2%z(16) =   0.0000000000000000_dp ; qsu2%w(16) = 0.35903916041026207_dp
        qsu2%x(17) =  0.45970084338098310_dp ; qsu2%y(17) = -0.88807383397711526_dp ; qsu2%z(17) =   0.0000000000000000_dp ; qsu2%w(17) = 0.35903916041026207_dp
        qsu2%x(18) = -0.45970084338098310_dp ; qsu2%y(18) = -0.88807383397711526_dp ; qsu2%z(18) =   0.0000000000000000_dp ; qsu2%w(18) = 0.35903916041026207_dp
        qsu2%x(19) =  0.88807383397711526_dp ; qsu2%y(19) =  0.45970084338098310_dp ; qsu2%z(19) =   0.0000000000000000_dp ; qsu2%w(19) = 0.35903916041026207_dp
        qsu2%x(20) = -0.88807383397711526_dp ; qsu2%y(20) =  0.45970084338098310_dp ; qsu2%z(20) =   0.0000000000000000_dp ; qsu2%w(20) = 0.35903916041026207_dp
        qsu2%x(21) =  0.88807383397711526_dp ; qsu2%y(21) = -0.45970084338098310_dp ; qsu2%z(21) =   0.0000000000000000_dp ; qsu2%w(21) = 0.35903916041026207_dp
        qsu2%x(22) = -0.88807383397711526_dp ; qsu2%y(22) = -0.45970084338098310_dp ; qsu2%z(22) =   0.0000000000000000_dp ; qsu2%w(22) = 0.35903916041026207_dp
        qsu2%x(23) =  0.45970084338098310_dp ; qsu2%y(23) =   0.0000000000000000_dp ; qsu2%z(23) =  0.88807383397711526_dp ; qsu2%w(23) = 0.35903916041026207_dp
        qsu2%x(24) = -0.45970084338098310_dp ; qsu2%y(24) =   0.0000000000000000_dp ; qsu2%z(24) =  0.88807383397711526_dp ; qsu2%w(24) = 0.35903916041026207_dp
        qsu2%x(25) =  0.45970084338098310_dp ; qsu2%y(25) =   0.0000000000000000_dp ; qsu2%z(25) = -0.88807383397711526_dp ; qsu2%w(25) = 0.35903916041026207_dp
        qsu2%x(26) = -0.45970084338098310_dp ; qsu2%y(26) =   0.0000000000000000_dp ; qsu2%z(26) = -0.88807383397711526_dp ; qsu2%w(26) = 0.35903916041026207_dp
        qsu2%x(27) =  0.88807383397711526_dp ; qsu2%y(27) =   0.0000000000000000_dp ; qsu2%z(27) =  0.45970084338098310_dp ; qsu2%w(27) = 0.35903916041026207_dp
        qsu2%x(28) = -0.88807383397711526_dp ; qsu2%y(28) =   0.0000000000000000_dp ; qsu2%z(28) =  0.45970084338098310_dp ; qsu2%w(28) = 0.35903916041026207_dp
        qsu2%x(29) =  0.88807383397711526_dp ; qsu2%y(29) =   0.0000000000000000_dp ; qsu2%z(29) = -0.45970084338098310_dp ; qsu2%w(29) = 0.35903916041026207_dp
        qsu2%x(30) = -0.88807383397711526_dp ; qsu2%y(30) =   0.0000000000000000_dp ; qsu2%z(30) = -0.45970084338098310_dp ; qsu2%w(30) = 0.35903916041026207_dp
        qsu2%x(31) =   0.0000000000000000_dp ; qsu2%y(31) =  0.45970084338098310_dp ; qsu2%z(31) =  0.88807383397711526_dp ; qsu2%w(31) = 0.35903916041026207_dp
        qsu2%x(32) =   0.0000000000000000_dp ; qsu2%y(32) = -0.45970084338098310_dp ; qsu2%z(32) =  0.88807383397711526_dp ; qsu2%w(32) = 0.35903916041026207_dp
        qsu2%x(33) =   0.0000000000000000_dp ; qsu2%y(33) =  0.45970084338098310_dp ; qsu2%z(33) = -0.88807383397711526_dp ; qsu2%w(33) = 0.35903916041026207_dp
        qsu2%x(34) =   0.0000000000000000_dp ; qsu2%y(34) = -0.45970084338098310_dp ; qsu2%z(34) = -0.88807383397711526_dp ; qsu2%w(34) = 0.35903916041026207_dp
        qsu2%x(35) =   0.0000000000000000_dp ; qsu2%y(35) =  0.88807383397711526_dp ; qsu2%z(35) =  0.45970084338098310_dp ; qsu2%w(35) = 0.35903916041026207_dp
        qsu2%x(36) =   0.0000000000000000_dp ; qsu2%y(36) = -0.88807383397711526_dp ; qsu2%z(36) =  0.45970084338098310_dp ; qsu2%w(36) = 0.35903916041026207_dp
        qsu2%x(37) =   0.0000000000000000_dp ; qsu2%y(37) =  0.88807383397711526_dp ; qsu2%z(37) = -0.45970084338098310_dp ; qsu2%w(37) = 0.35903916041026207_dp
        qsu2%x(38) =   0.0000000000000000_dp ; qsu2%y(38) = -0.88807383397711526_dp ; qsu2%z(38) = -0.45970084338098310_dp ; qsu2%w(38) = 0.35903916041026207_dp
      case default
        stop "order not implemented"
      end select
    case("SD")
      select case (qsu2%o)
      case (0)
        qsu2%x( 1) = 0._dp; qsu2%y(1)=0._dp; qsu2%z(1)=1._dp; qsu2%w(1)=1._dp
      case(1)
        qsu2%x(1)=0._dp;  qsu2%y(1)= 0._dp ; qsu2%z(1)=   1._dp  ; qsu2%w(1)=  0.5_dp
        qsu2%x(2)=0._dp;  qsu2%y(2)= 0._dp ; qsu2%z(2)=  -1._dp  ; qsu2%w(2)=  0.5_dp
      case(2)
        qsu2%x(1)= 0.57735026918962573_dp  ;qsu2%y(1)= 0.57735026918962573_dp ;qsu2%z(1)= 0.57735026918962573_dp  ; qsu2%w(1)=   0.25_dp
        qsu2%x(2)=-0.57735026918962573_dp  ;qsu2%y(2)=-0.57735026918962573_dp ;qsu2%z(2)= 0.57735026918962573_dp  ; qsu2%w(2)=   0.25_dp
        qsu2%x(3)=-0.57735026918962573_dp  ;qsu2%y(3)= 0.57735026918962573_dp ;qsu2%z(3)=-0.57735026918962573_dp ; qsu2%w(3)=    0.25_dp
        qsu2%x(4)= 0.57735026918962573_dp  ;qsu2%y(4)=-0.57735026918962573_dp ;qsu2%z(4)=-0.57735026918962573_dp ; qsu2%w(4)=    0.25_dp
      case(3)
        qsu2%x(1)= 1._dp  ;qsu2%y(1)=      0._dp  ;qsu2%z(1)=      0._dp ; qsu2%w(1)=0.166666666666666667_dp
        qsu2%x(2)= 0._dp  ;qsu2%y(2)=      1._dp  ;qsu2%z(2)=      0._dp ; qsu2%w(2)=0.166666666666666667_dp
        qsu2%x(3)= 0._dp  ;qsu2%y(3)=      0._dp  ;qsu2%z(3)=      1._dp ; qsu2%w(3)=0.166666666666666667_dp
        qsu2%x(4)=-1._dp ;qsu2%y(4)=      0._dp  ;qsu2%z(4)=      0._dp ; qsu2%w(4)=0.166666666666666667_dp
        qsu2%x(5)= 0._dp  ;qsu2%y(5)=     -1._dp  ;qsu2%z(5)=      0._dp ; qsu2%w(5)=0.166666666666666667_dp
        qsu2%x(6)= 0._dp  ;qsu2%y(6)=      0._dp  ;qsu2%z(6)=     -1._dp ; qsu2%w(6)=0.166666666666666667_dp
      case(4)
         qsu2%x(1)=0._dp  ;qsu2%y(1)= 0._dp ;qsu2%z(1)= 1._dp ;qsu2%w(1)=0.0833333333333333333_dp
         qsu2%x(2)=0._dp  ;qsu2%y(2)= 0._dp ;qsu2%z(2)=-1._dp ;qsu2%w(2)=0.0833333333333333333_dp
         qsu2%x(3)=0.047060422787  ;qsu2%y(3)=  0.89318828732  ;qsu2%z(3)= 0.4472135955 ;qsu2%w(3)=      0.10416666666666667_dp
         qsu2%x(4)=0.66485623892   ;qsu2%y(4)= 0.59830275076   ;qsu2%z(4)=-0.4472135955 ;qsu2%w(4)=      0.10416666666666667_dp
         qsu2%x(5)=-0.047060422787 ;qsu2%y(5)= -0.89318828732  ;qsu2%z(5)=0.4472135955  ;qsu2%w(5)=      0.10416666666666667_dp
         qsu2%x(6)=-0.66485623892  ;qsu2%y(6)= -0.59830275076  ;qsu2%z(6)= -0.4472135955;qsu2%w(6)=      0.10416666666666667_dp
         qsu2%x(7)=-0.89318828732  ;qsu2%y(7)= 0.047060422787  ;qsu2%z(7)=0.4472135955  ;qsu2%w(7)=      0.10416666666666667_dp
         qsu2%x(8)=-0.59830275076  ;qsu2%y(8)= 0.66485623892   ;qsu2%z(8)=-0.4472135955 ;qsu2%w(8)=      0.10416666666666667_dp
         qsu2%x(9)=0.89318828732   ;qsu2%y(9)= -0.047060422787 ;qsu2%z(9)=0.4472135955  ;qsu2%w(9)=      0.10416666666666667_dp
         qsu2%x(10)=0.59830275076  ;qsu2%y(10)= -0.66485623892  ;qsu2%z(10)=-0.4472135955 ;qsu2%w(10)=   0.10416666666666667_dp
      case(5)
        qsu2%x(1 )=  0.85065080835   ; qsu2%y(1 )=  0.52573111212   ; qsu2%z(1 )=  0.00000000000  ; qsu2%w(1 )=0.083333333333
        qsu2%x(2 )=  0.85065080835   ; qsu2%y(2 )=  -0.52573111212  ; qsu2%z(2 )=  0.00000000000  ; qsu2%w(2 )=0.083333333333
        qsu2%x(3 )=  -0.85065080835  ; qsu2%y(3 )=  0.52573111212   ; qsu2%z(3 )=  0.00000000000  ; qsu2%w(3 )=0.083333333333
        qsu2%x(4 )=  -0.85065080835  ; qsu2%y(4 )=  -0.52573111212  ; qsu2%z(4 )=  0.00000000000  ; qsu2%w(4 )=0.083333333333
        qsu2%x(5 )=  0.00000000000   ; qsu2%y(5 )=  0.85065080835   ; qsu2%z(5 )=  0.52573111212  ; qsu2%w(5 )=0.083333333333
        qsu2%x(6 )=  0.00000000000   ; qsu2%y(6 )=  0.85065080835   ; qsu2%z(6 )=  -0.52573111212 ; qsu2%w(6 )=0.083333333333
        qsu2%x(7 )=  0.00000000000   ; qsu2%y(7 )=  -0.85065080835  ; qsu2%z(7 )=  0.52573111212  ; qsu2%w(7 )=0.083333333333
        qsu2%x(8 )=  0.00000000000   ; qsu2%y(8 )=  -0.85065080835  ; qsu2%z(8 )=  -0.52573111212 ; qsu2%w(8 )=0.083333333333
        qsu2%x(9 )=  0.52573111212   ; qsu2%y(9 )=  0.00000000000   ; qsu2%z(9 )=  0.85065080835  ; qsu2%w(9 )=0.083333333333
        qsu2%x(10)=  -0.52573111212  ; qsu2%y(10)=  0.00000000000   ; qsu2%z(10)=  0.85065080835  ; qsu2%w(10)=0.083333333333
        qsu2%x(11)=  0.52573111212   ; qsu2%y(11)=  0.00000000000   ; qsu2%z(11)=  -0.85065080835 ; qsu2%w(11)=0.083333333333
        qsu2%x(12)=  -0.52573111212  ; qsu2%y(12)=  0.00000000000   ; qsu2%z(12)=  -0.85065080835 ; qsu2%w(12)=0.083333333333
      case default
        stop "Spherical design quadratures of order >5 not implemented"
      end select
      qsu2%w = qsu2%w * fourpi ! tabulates weights for SD are normalized to 1
    case default
      stop "your quadrature is not implemented."
    end select

    ! Correspondance between Euler angles and quadrature node
    block
      real(dp) :: x,y,z,w,theta,r,phi,psi,cos_phi,sin_phi,cos_theta,sin_theta,sin_psi,cos_psi
      do i=1,N ! loop over all quadrature nodes
        x = qsu2%x(i)
        y = qsu2%y(i)
        z = qsu2%z(i)
        w = qsu2%w(i)
        theta = acos(z)
        r = sqrt(x**2+y**2)

        ! phi
        if (r <= epsdp) then
          phi = 0
        else if (y >= 0) then
          phi = acos(x/r)
        else if (y < 0) then
          phi = twopi-acos(x/r)
        else
          error stop "something is wrong with phi, this else statement should not appear"
        end if

        cos_phi = cos(phi)
        sin_phi = sin(phi)
        cos_theta = cos(theta)
        sin_theta = sin(theta)

        OMx(i) = chop( cos_phi*sin_theta )
        OMy(i) = chop( sin_phi*sin_theta )
        OMz(i) = chop( cos_theta         )

        P = molrotgrid%n_angles
        if (P<=0) stop "problem sur P in module_quadrature"
        do concurrent (j=1:P)
          cos_psi = cos(molrotgrid%root(j))
          sin_psi = sin(molrotgrid%root(j))
          Rotxx(i,j) =  cos_theta*cos_phi*cos_psi-sin_phi*sin_psi
          Rotxy(i,j) = -cos_theta*cos_phi*sin_psi-sin_phi*cos_psi
          Rotxz(i,j) =  sin_theta*cos_phi
          Rotyx(i,j) =  cos_theta*sin_phi*cos_psi+cos_phi*sin_psi
          Rotyy(i,j) = -cos_theta*sin_phi*sin_psi+cos_phi*cos_psi
          Rotyz(i,j) =  sin_theta*sin_phi
          Rotzx(i,j) = -sin_theta*cos_psi
          Rotzy(i,j) =  sin_theta*sin_psi
          Rotzz(i,j) =  cos_theta
        end do
      end do
    end block

    ! Assert sum of all weights is 2pi for psi and 4pi for phi theta
    if( abs( twopi - molrotsymorder*sum( molrotgrid%weight )) > 10*epsdp  ) then
        print*,"sum of all weights over psi is",sum(molrotgrid%weight) *real(molrotsymorder,dp)
        print*,"should be 2pi=",twopi
        print*,"the difference is",twopi - sum(molrotgrid%weight) *real(molrotsymorder,dp)
        print*,"epsilon machine=",epsdp
        error stop
    end if

    if( abs( fourpi -sum(qsu2%w) ) > 10*epsdp ) then
        print*,"sum of all weights over phi theta is",sum(qsu2%w)
        print*,"should be 4pi=",fourpi
        print*,"the difference is",fourpi-sum(qsu2%w)
        print*,"epsilon machine=",epsdp
        error stop
    end if

    N = qsu2%n
    allocate( anggrid%weight(N) ,source=qsu2%w )
    angGrid%n_angles = N
    anggrid%n = anggrid%n_angles ! just a shortcut

    !!!
!    print*,"anggrid"
!    print*,"n_angles=",anggrid%n_angles
!    print*,"weights=",anggrid%weight
!    error stop 'subroutine init module quadrature'
  end subroutine init

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
  end subroutine gauss_legendre


end module quadrature

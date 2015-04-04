! Module for numerical integration
module quadrature

  use precision_kinds ,only: dp, i2b
  use constants       ,only: pi, twopi, fourpi, zero
  use input           ,only: input_log, input_char, input_int

  implicit none

  real(dp), allocatable, dimension(:), private :: x_leb, y_leb , z_leb
  real(dp), allocatable, dimension(:), public :: Omx , Omy , Omz  ! unit vector for Euler angle Omega in lab frame
  real(dp), allocatable, dimension(:,:), public :: Rotxx, Rotxy, Rotxz, Rotyx, Rotyy, Rotyz, Rotzx, Rotzy, Rotzz

  type angularGrid
      integer(i2b) :: n_angles, N, o
      real(dp), allocatable, dimension(:) :: weight, root, w, x
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
  integer(i2b), public :: molRotSymOrder

  contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init

    implicit none

    molRotSymOrder = input_int('molRotSymOrder', defaultvalue=1)
    if ( molRotSymOrder <= 0 ) then
        print*,"in module_quadrature, molRotSymOrder, readen from input/dft.in, is unphysical. critical stop"
        stop
    end if

    ! molecular rotation grid
    call get_psi_integration_roots_and_weights (molRotGrid, molRotSymOrder)
    call check_weights_psi(molRotGrid%weight)

    ! integration scheme
    intScheme%name = input_char('quadrature')

    select case (intscheme%name)
    case ('L')
      print*,"Quadrature for SU2 : Lebedev"
    case ('GL')
      print*,"Quadrature for SU2 : Gauss-Legendre"
    case ('SD')
      print*,"Quadrature for SU2 : Spherical Design"
    case default
      print*, "To ingrate over the unit sphere, you want a quadrature called", intscheme%name
      stop "It is not implemented"
    end select

    call read_order_of_quadrature (intScheme%order)

    allocate (intScheme%weight(intScheme%order), source=zero)
    allocate (intScheme%root(intScheme%order), source=zero)


!    call qsu2_get (qsu2%N, qsu2%w, qsu2%x)
!    ask order of quadrature (qsu2%o), in the sense the maximum order of the polynome that can be exactly integrated by the quadrature
!    you will get back number of roots, their coordinates and weights

    select case (intscheme%name)
    case ("L")
      select case (qsu2%o)
      case (3)
        qsu2%n = 6
      case (5)
        qsu2%n = 14
      case (7)
        qsu2%n = 26
      case (9)
        qsu2%n = 38
      case default
        print*, "You ask for a quadrature over unit sphere of order",qsu2%o,"that is not implemented with Lebedev."
      end select
    case ("GL")
      if ( (qsu2%o/2)*2 == qsu2%o ) then ! i.e. if order is odd
        qsu2%n = (qsu2%o +1)**2/2
      else
        print*, "You ask for a quadrature over unit sphere of order",qsu2%o,"that is not compatible with Gauss-Legendre quadratures"
        stop
      end if
    case ("SDG")
      select case (qsu2%o)
      case (1)
        qsu2%n = 2
      case (2)
        qsu2%n = 4
      case (3)
        qsu2%n = 6
      case (4)
        qsu2%n = 10
      case (5)
        qsu2%n = 12
      case (6)
        qsu2%n = 18
      case (7)
        qsu2%n = 22
      case (8)
        qsu2%n = 28
      case (9)
        qsu2%n = 32
      case (10)
        qsu2%n = 42
      case (11)
        qsu2%n = 48
      case (12)
        qsu2%n = 58
      case (13)
        qsu2%n = 64
      case (14)
        qsu2%n = 72
      case (15)
        qsu2%n = 82
      case (16)
        qsu2%n = 98
      case (17)
        qsu2%n = 104
      case (18)
        qsu2%n = 122
      case (19)
        qsu2%n = 130
      case (20)
        qsu2%n = 148
      case (21)
        qsu2%n = 156
      case (22)
        qsu2%n = 178
      case (23)
        qsu2%n = 186
      case (24)
        qsu2%n = 210
      case (25)
        qsu2%n = 220
      case (26)
        qsu2%n = 244
      case (27)
        qsu2%n = 254
      case (28)
        qsu2%n = 282
      case (29)
        qsu2%n = 292
      case (30)
        qsu2%n = 322
      case (32)
        qsu2%n = 364
      case (34)
        qsu2%n = 410
      case (35)
        qsu2%n = 422
      case (36)
        qsu2%n = 458
      case (37)
        qsu2%n = 472
      case (38)
        qsu2%n = 508
      case (39)
        qsu2%n = 522
      case (44)
        qsu2%n = 672
      case default
        print*, "You ask for a quadrature over unit sphere of order", qsu2%o," that is not compatible with Gauss type spherical designs"
      end select
    end select


    ! LEBEDEV
    if (intScheme%name=='L') then
        angGrid%n_angles = intScheme%order
        call allocate_Rotij (angGrid%n_angles,molRotGrid%n_angles,Rotxx,Rotxy,Rotxz,Rotyx,Rotyy,Rotyz,Rotzx,Rotzy,Rotzz)
        allocate( x_leb(intScheme%order), y_leb(intScheme%order), z_leb(intScheme%order) )
        allocate (angGrid%weight(angGrid%n_angles), source=zero)
        call lebedev_integration_roots_and_weights (intScheme%order, x_leb ,y_leb , z_leb, intScheme%weight)
        call lebedev (angGrid, intScheme, molRotGrid, Rotxx, Rotxy, Rotxz, Rotyx, Rotyy, Rotyz, Rotzx, Rotzy, Rotzz)
        deallocate( x_leb, y_leb, z_leb)

    ! GAUSS-LEGENDRE
    else if (intScheme%name=='GL') then
        select case (intScheme%order)
        case (1)
            angGrid%n_angles = 1
        case default
            angGrid%n_angles = 2*(intScheme%order**2)
        end select
        call allocate_Rotij (angGrid%n_angles,molRotGrid%n_angles,Rotxx,Rotxy,Rotxz,Rotyx,Rotyy,Rotyz,Rotzx,Rotzy,Rotzz)
        call gauss_legendre_integration_roots_and_weights (intScheme%order, intScheme%weight , intScheme%root)
        allocate (angGrid%weight(angGrid%n_angles), source=0._dp)
        call gauss (angGrid, intScheme, molRotGrid, Rotxx, Rotxy, Rotxz, Rotyx, Rotyy, Rotyz, Rotzx, Rotzy, Rotzz)
    end if
print*,anggrid%n_angles
print*,intscheme%root;stop
  end subroutine init

  ! GAUSS
  subroutine gauss (angGrid, intScheme, molRotGrid, Rotxx, Rotxy, Rotxz, Rotyx, Rotyy, Rotyz, Rotzx, Rotzy, Rotzz)
      type (angularGrid), intent(inout) :: angGrid
      type (angularGrid), intent(in) :: molRotGrid
      type (integrationScheme), intent(in) :: intScheme
      real(dp), dimension(:,:), intent(out) :: Rotxx,Rotxy,Rotxz,Rotyx,Rotyy,Rotyz,Rotzx,Rotzy,Rotzz
      integer(i2b) :: omega, n_psi, theta, n_phi
      real(dp) :: phi, cos_theta, sin_theta, cos_phi, sin_phi, cos_psi, sin_psi

      allocate(OMx(angGrid%n_angles) , SOURCE=0._dp)
      allocate(OMy(angGrid%n_angles) , SOURCE=0._dp)
      allocate(OMz(angGrid%n_angles) , SOURCE=0._dp)

      select case (intScheme%order)
      case (1)
          angGrid%weight (1)  = fourpi
          Rotxx = 1.0_dp ; Rotxy = 0.0_dp ; Rotxz = 0.0_dp
          Rotyx = 0.0_dp ; Rotyy = 1.0_dp ; Rotyz = 0.0_dp
          Rotzx = 0.0_dp ; Rotzy = 0.0_dp ; Rotzz = 1.0_dp
          OMx = 0._dp ; OMy = 0._dp ; OMz = 1._dp ! ATTENTION completely arbitrary. We decide to put it along z.
      case default
          do  theta = 1, intScheme%order
              cos_theta = intScheme%root(theta)
              sin_theta = sqrt ( 1.0_dp - cos_theta ** 2 )
              do  n_phi = 1, 2*intScheme%order
                  omega = (theta-1)*2*intScheme%order + n_phi
                  phi = real(n_phi-1,dp) * twopi / real ( 2 * intScheme%order , dp ) !0,pi
                  cos_phi = cos ( phi ) ; sin_phi = sin ( phi )
                  OMx ( omega ) = sin_theta * cos_phi
                  OMy ( omega ) = sin_theta * sin_phi
                  OMz ( omega ) = cos_theta
                  angGrid%weight(omega) = intScheme%weight(theta) *pi /real(intScheme%order,dp)
                  do n_psi = 1 , molRotGrid%n_angles
                      cos_psi = cos(  molRotGrid%root(n_psi)  )
                      sin_psi = sin(  molRotGrid%root(n_psi)  )
                      Rotxx(omega,n_psi) =  cos_theta*cos_phi*cos_psi-sin_phi*sin_psi
                      Rotxy(omega,n_psi) = -cos_theta*cos_phi*sin_psi-sin_phi*cos_psi
                      Rotxz(omega,n_psi) =  sin_theta*cos_phi
                      Rotyx(omega,n_psi) =  cos_theta*sin_phi*cos_psi+cos_phi*sin_psi
                      Rotyy(omega,n_psi) = -cos_theta*sin_phi*sin_psi+cos_phi*cos_psi
                      Rotyz(omega,n_psi) =  sin_theta*sin_phi
                      Rotzx(omega,n_psi) = -sin_theta*cos_psi
                      Rotzy(omega,n_psi) =  sin_theta*sin_psi
                      Rotzz(omega,n_psi) =  cos_theta
                  end DO
              end DO
          end DO
      end select

      call check_weights(angGrid%weight)

  end subroutine gauss


  ! LEBEDEV
  subroutine lebedev (angGrid, intScheme, molRotGrid, Rotxx, Rotxy, Rotxz, Rotyx, Rotyy, Rotyz, Rotzx, Rotzy, Rotzz)

      type (angularGrid), intent(inout) :: angGrid
      type (angularGrid), intent(in) :: molRotGrid
      type (integrationScheme), intent(in) :: intScheme
      real(dp), dimension(:,:), intent(out) :: Rotxx,Rotxy,Rotxz,Rotyx,Rotyy,Rotyz,Rotzx,Rotzy,Rotzz
      integer(i2b) ::  n, psi
      real(dp) :: phi, theta, cos_theta, sin_theta, cos_phi, sin_phi, cos_psi, sin_psi

      allocate(OMx(angGrid%n_angles))
      allocate(OMy(angGrid%n_angles))
      allocate(OMz(angGrid%n_angles))

      angGrid%weight(1:angGrid%n_angles) = intScheme%weight(1:angGrid%n_angles)
      call check_weights(angGrid%weight)

      do concurrent (n=1:angGrid%n_angles)

          if ( abs(x_leb(n))<=epsilon(1.0_dp) .and. abs(y_leb(n))<=epsilon(1._dp) ) then
              phi = zero
          else if (y_leb(n)>=zero) then
              phi = acos(x_leb(n)/(sqrt ( x_leb(n)**2 + y_leb(n)**2 ) ) )
          else
              phi = twopi - acos(x_leb(n)/(sqrt ( x_leb(n)**2 + y_leb(n)**2 ) ) )
          end if
          theta = acos(z_leb(n))
          OMx (n) = sin(theta) * cos(phi)
          OMy (n) = sin(theta) * sin(phi)
          OMz (n) = cos(theta)
      end DO

      do concurrent (n=1:angGrid%n_angles)
          theta=acos(z_leb(n))
          cos_theta=cos(theta)
          sin_theta=sin(theta)

          if    ( abs(x_leb(n)) <= epsilon(1._dp) .and. abs(y_leb(n)) <= epsilon(1._dp) ) then
              phi=0.0_dp
          else if  (y_leb(n)>=0.0_dp) then
              phi=acos(x_leb(n)/(sqrt ( x_leb(n)**2 + y_leb(n)**2 ) ) )
          else
              phi=twopi - acos(x_leb(n)/(sqrt ( x_leb(n)**2 + y_leb(n)**2 ) ) )
          end if
          cos_phi=cos(phi)
          sin_phi=sin(phi)
          do concurrent (psi=1:molRotGrid%n_angles)
              cos_psi = cos(  molRotGrid%root(psi)  )
              sin_psi = sin(  molRotGrid%root(psi)  )
              Rotxx(n,psi) =  cos_theta*cos_phi*cos_psi-sin_phi*sin_psi
              Rotxy(n,psi) = -cos_theta*cos_phi*sin_psi-sin_phi*cos_psi
              Rotxz(n,psi) =  sin_theta*cos_phi
              Rotyx(n,psi) =  cos_theta*sin_phi*cos_psi+cos_phi*sin_psi
              Rotyy(n,psi) = -cos_theta*sin_phi*sin_psi+cos_phi*cos_psi
              Rotyz(n,psi) =  sin_theta*sin_phi
              Rotzx(n,psi) = -sin_theta*cos_psi
              Rotzy(n,psi) =  sin_theta*sin_psi
              Rotzz(n,psi) =  cos_theta
          end DO
      end DO
  end subroutine lebedev

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine check_weights_psi(weight_over_molecular_rotations)
      implicit none
      real(dp), dimension(:), intent(in) :: weight_over_molecular_rotations
      if ( abs ( sum ( weight_over_molecular_rotations( : ) ) - twopi/molRotSymOrder )  > 1.0e-10_dp ) then
          print*, 'problem detected in module_quadrature.f90 :'
          print*, 'sum over omegas of molRotGrid%weight(omega) is not 2pi/molRotSymOrder, it is ',&
              sum ( weight_over_molecular_rotations )
          stop 'CRITICAL'
      end if
  end subroutine check_weights_psi

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine check_weights(weight)
      implicit none
      real(dp), dimension(:), intent(in) :: weight
      ! check if sum over all omega of weight(omega) is fourpi
      if ( abs ( sum ( weight ( : ) ) - fourpi )  > 1.0e-10_dp ) then
          print *, 'problem detected in compute_angular_grid.f90 :'
          print *, 'sum over omegas of weight(omega) is not 4pi. it is ',sum ( weight )
          stop 'CRITICAL'
      end if
  end subroutine check_weights

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine get_psi_integration_roots_and_weights (molRotGrid, molRotSymOrder)
      type (angularGrid), intent(inout) :: molRotGrid
      integer(i2b), intent(out) :: molRotSymOrder
      integer(i2b) :: i
      molRotGrid%n_angles = input_int('nb_psi')
      if ( molRotGrid%n_angles < 1 ) then
          print*,"in module_quadrature, nb_psi, readen from input/dft.in is unphysical. critical stop."
          stop
      end if
      allocate( molRotGrid%weight(molRotGrid%n_angles), source=0._dp)
      molRotGrid%weight = twopi/real(molRotGrid%n_angles*molRotSymOrder,dp) ! equidistant repartition between 0 and 2Pi => all weights are equal
      allocate( molRotGrid%root(molRotGrid%n_angles), source=0._dp)
      do concurrent (i=1:molRotGrid%n_angles) ! equidistant repartition between 0 and 2Pi
          molRotGrid%root(i) = real(i-1,dp)*twopi/real(molRotGrid%n_angles*molRotSymOrder,dp) ! roots
      end DO
  end subroutine get_psi_integration_roots_and_weights

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  pure subroutine  allocate_Rotij (n_angles, n_psi,Rotxx,Rotxy,Rotxz,Rotyx,Rotyy,Rotyz,Rotzx,Rotzy,Rotzz)
      integer(i2b), intent(in) :: n_angles, n_psi
      real(dp), allocatable, dimension(:,:), intent(out) :: Rotxx,Rotxy,Rotxz,Rotyx,Rotyy,Rotyz,Rotzx,Rotzy,Rotzz
      allocate( Rotxx(n_angles,n_psi), Rotxy(n_angles,n_psi), Rotxz(n_angles,n_psi) )
      allocate( Rotyx(n_angles,n_psi), Rotyy(n_angles,n_psi), Rotyz(n_angles,n_psi) )
      allocate( Rotzx(n_angles,n_psi), Rotzy(n_angles,n_psi), Rotzz(n_angles,n_psi) )
  end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine read_order_of_quadrature (order_of_quadrature)
    integer(i2b), intent(out) :: order_of_quadrature
    order_of_quadrature = input_int('order_of_quadrature',defaultvalue=2)
    if (order_of_quadrature < 1) then
      print*,"In input/dft.in, I read order_of_quadrature is < 1. This is not physical as you must have"
      print*,"at least one point, even without angular grid. CRITICAL"
      stop
    end if
  end subroutine read_order_of_quadrature

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine qsu2_get (N, w, x)
    implicit none
    integer, intent(in) :: N
    real(dp), allocatable, intent(out) :: w(:), x(:)
    select case (N)
    case (1)

      allocate( w(2), source=[0.5_dp,0.5_dp] )
    case (2)
    case (3)
    case (4)
    case (5)
    case default
      print*,"such high order in quadrature for su2 is not implemented"
    end select
  end subroutine qsu2_get

end module quadrature

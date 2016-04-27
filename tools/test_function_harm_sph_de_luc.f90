program test_function_wigner_small_d

  use iso_c_binding, only: dp=>c_double
  implicit none
  integer :: i, N
  double precision :: theta
  double precision, parameter :: pi = acos(-1._dp)

  !
  ! on va faire N fois le test
  !
  do N=1,20

    !
    ! genere un theta aleatoire entre 0 et 1
    !
    call random_number(theta)

    !
    ! multiplie par pi pour que ce soit entre 0 et pi
    !
    theta = theta * pi

    print*, "100", cos(theta)            - wigner_small_d(1,0,0,theta)
    print*, "111", (1+cos(theta))/2.     - wigner_small_d(1,1,1,theta)
    print*, "110", -sin(theta)/sqrt(2.)  - wigner_small_d(1,1,0,theta)
    print*, "11-1", (1-cos(theta))/2.    - wigner_small_d(1,1,-1,theta)
    print*, "222",  ((1+cos(theta))/2.)**2 - wigner_small_d(2,2,2,theta)
    print*, "221", -(1+cos(theta))/2.*sin(theta) - wigner_small_d(2,2,1,theta)
    print*, "220", sqrt(6.)/4.*sin(theta)**2 - wigner_small_d(2,2,0,theta)
    print*, "22-1", -(1-cos(theta))/2.*sin(theta) - wigner_small_d(2,2,-1,theta)
    print*, "22-2", ((1-cos(theta))/2.)**2 - wigner_small_d(2,2,-2,theta)
    print*, "211", (1+cos(theta))/2.*(2*cos(theta)-1) - wigner_small_d(2,1,1,theta)
    print*, "210", -sqrt(3./2.)*sin(theta)*cos(theta) - wigner_small_d(2,1,0,theta)
    print*, "21-1", (1-cos(theta))/2.*(2*cos(theta)+1) - wigner_small_d(2,1,-1,theta)
    print*, "200", (3./2.*cos(theta)**2 -0.5) - wigner_small_d(2,0,0,theta)

  end do


  ! OK

  contains

  pure function wigner_small_d( m, mu, mup, beta ) ! Luc's luc72p143
      !
      ! R^m_{mu,mup}(\beta)
      ! beta is the angle in radian. It is often named theta in molecular physics
      !
      use iso_c_binding, only: dp => c_double
      use mathematica, only: fact
      implicit none
      REAL(dp) :: wigner_small_d
      INTEGER, intent(in) :: m, mu, mup
      REAL(dp), intent(in) :: beta
      REAL(dp) :: beta0, x, cc, ss, pm, pm1, pm2
      INTEGER :: mu0, l, it
      if( abs(mu)>m .or. abs(mup)>m ) then
          wigner_small_d = 0
          return
      end if
      if( m==0 ) then
          wigner_small_d = 1
          return
      end if
      if( mu==0 .or. mup==0 ) then
          mu0 = mu
          beta0 = beta
          if(mu==0) then            ! je mets le 0 en second
              mu0 = mup
              beta0 = -beta
          end if
          x = 1                 ! si mu negatif, ca vaut (-1)**mu * valeur pour -mu
          if(mu0<0) then
              x = (-1)**mu0
              mu0 = -mu0
          end if
          cc = cos(beta0)           ! plutôt a partir d'une formule de recurrence stable
          pm1 = 0.               ! des polynomes de legendre associes pl,m
          pm = 1.                ! luc73p96
          do l = mu0+1,m
              pm2 = pm1
              pm1 = pm
              pm = (cc*REAL(2*l-1,dp)*pm1-REAL(l+mu0-1,dp)*pm2)/REAL(l-mu0,dp)
          end do
          wigner_small_d = x*(-1.)**mu0*sqrt(fact(m-mu0)/fact(m+mu0))&
          *fact(2*mu0)/(2.**mu0*fact(mu0))*sin(beta0)**mu0*pm
      else                        ! donc mu et mup non nuls, utiliser betement la formule de wigner
          wigner_small_d = 0
          cc = cos(0.5*beta)
          ss = sin(0.5*beta)
          do it = max(0,mu-mup),min(m+mu,m-mup)
              wigner_small_d = wigner_small_d+(-1.)**it/(fact(m+mu-it)*fact(m-mup-it)*fact(it)*&
              fact(it-mu+mup))*cc**(2*m+mu-mup-2*it)*ss**(2*it-mu+mup)
          end do
          wigner_small_d = sqrt(fact(m+mu)*fact(m-mu)*fact(m+mup)*fact(m-mup))*wigner_small_d
      end if
  end function wigner_small_d

end program

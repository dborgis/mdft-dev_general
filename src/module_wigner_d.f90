module module_wigner_d
  private
  public :: wigner_small_d, wigner_big_d, symbol_3j
contains

  pure function wigner_small_d( m, mu, mup, theta ) ! Luc's luc72p143
      !
      ! r^m_{mu,mup}(\theta)
      ! theta is the angle in radian.
      ! see http://pdg.lbl.gov/2015/reviews/rpp2014-rev-clebsch-gordan-coefs.pdf  <= we have the same normalization etc for d^m_{mu,mup}(theta)
      !
      use precision_kinds, only: dp
      use mathematica, only: fact
      implicit none
      REAL(dp) :: wigner_small_d
      INTEGER, intent(in) :: m, mu, mup
      REAL(dp), intent(in) :: theta
      REAL(dp) :: theta0, x, cc, ss, pm, pm1, pm2
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
          theta0 = theta
          if(mu==0) then            ! je mets le 0 en second
              mu0 = mup
              theta0 = -theta
          end if
          x = 1                 ! si mu negatif, ca vaut (-1)**mu * valeur pour -mu
          if(mu0<0) then
              x = (-1)**mu0
              mu0 = -mu0
          end if
          cc = cos(theta0)           ! plutôt a partir d'une formule de recurrence stable
          pm1 = 0.               ! des polynomes de legendre associes pl,m
          pm = 1.                ! luc73p96
          do l = mu0+1,m
              pm2 = pm1
              pm1 = pm
              pm = (cc*REAL(2*l-1,dp)*pm1-REAL(l+mu0-1,dp)*pm2)/REAL(l-mu0,dp)
          end do
          wigner_small_d = x*(-1.)**mu0*sqrt(fact(m-mu0)/fact(m+mu0))&
          *fact(2*mu0)/(2.**mu0*fact(mu0))*sin(theta0)**mu0*pm
      else                        ! donc mu et mup non nuls, utiliser betement la formule de wigner
          wigner_small_d = 0
          cc = cos(0.5*theta)
          ss = sin(0.5*theta)
          do it = max(0,mu-mup),min(m+mu,m-mup)
              wigner_small_d = wigner_small_d+(-1.)**it/(fact(m+mu-it)*fact(m-mup-it)*fact(it)*&
              fact(it-mu+mup))*cc**(2*m+mu-mup-2*it)*ss**(2*it-mu+mup)
          end do
          wigner_small_d = sqrt(fact(m+mu)*fact(m-mu)*fact(m+mup)*fact(m-mup))*wigner_small_d
      end if
  end function wigner_small_d

  pure function wigner_big_D( m, mu, mup, theta, phi, psi )
    use precision_kinds, only: dp
    implicit none
    real(dp), intent(in) :: theta, phi, psi
    integer, intent(in) :: m, mu, mup
    complex(dp) :: wigner_big_D
    complex(dp), parameter :: ii=complex(0._dp,1._dp)
    wigner_big_D = wigner_small_d( m, mu, mup, theta ) * exp(-ii*(mu*phi+mup*psi))
  end function

  pure function symbol_3j(m,n,l,mu,nu,lu)
    use precision_kinds, only: dp
    use mathematica, only: fac=>fact
    implicit none
    integer, intent(in) :: m, n, l, mu, nu, lu
    real(dp) :: symbol_3j, som
    integer :: it
      ! COMMON/facto/fac(0:50)
    !
    !        symbole 3j
    !        Messiah page 910 eq.21
    !
    IF(itriangle(m,n,l)==0.or.mu+nu+lu/=0.or.                                     &
    ABS(mu).gt.m.or.abs(nu).gt.n.or.abs(lu).gt.l) then
    symbol_3j=0.
  else
    som=0.
    do it=MAX(0,n-l-mu,m-l+nu),MIN(m+n-l,m-mu,n+nu)
      som=som+(-1)**it/(fac(it)*fac(l-n+it+mu)*fac(l-m+it-nu)*                          &
      fac(m+n-l-it)*fac(m-it-mu)*fac(n-it+nu))
    end do
    symbol_3j=(-1)**(m-n-lu)*SQRT(delta(m,n,l))*                                       &
    SQRT(fac(m+mu)*fac(m-mu)*fac(n+nu)*fac(n-nu)*fac(l+lu)*fac(l-lu))                 &
    *som
    IF(mu==0.and.nu==0.and.lu==0.and.2*((m+n+l)/2)/=m+n+l) symbol_3j=0.
  endif
  end
  !
  pure function itriangle(m,n,l)
    implicit none
    integer, intent(in) :: m, n, l
    integer :: itriangle
    !        nul sauf si |m-n|<l<m+n
    !        rq: ne depend pas de l'ordre des 3 entiers
    itriangle=0
    IF(l>=ABS(m-n).and.l<=m+n) itriangle=1
  end
  !
  pure function delta(m,n,l)
    use mathematica, only: fac=>fact
    implicit none
    real(8) :: delta
    integer, intent(in) :: m, n, l
    ! COMMON/facto/fac(0:50)
    delta=fac(m+n-l)*fac(n+l-m)*fac(l+m-n)/fac(m+n+l+1)
  end




end module

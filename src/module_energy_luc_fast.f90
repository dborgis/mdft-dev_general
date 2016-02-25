module tableaux_dft_3d_fast
  INTEGER::mnmax,mnmax2,nbeta,nphi,nomeg
  REAL(8),DIMENSION(:),ALLOCATABLE::beta,cosbeta,sinbeta,wb
  REAL(8),DIMENSION(:),ALLOCATABLE::phi,cosphi,sinphi
  COMPLEX(8),DIMENSION(:),ALLOCATABLE::expiphi
  COMPLEX(8),DIMENSION(:,:),ALLOCATABLE::expikhiphi
  REAL(8),DIMENSION(:),ALLOCATABLE::omeg
  COMPLEX(8),DIMENSION(:,:,:),ALLOCATABLE::harsph1_q
  COMPLEX(8),DIMENSION(:,:,:),ALLOCATABLE:: delta_rho,delta_rho_proj1,auxi_1,auxi_2
  COMPLEX(8),DIMENSION(:,:,:),ALLOCATABLE:: gamma4_proj1
  INTEGER, DIMENSION(:),ALLOCATABLE:: mm,nn,ll,mumu,nunu
  REAL(8), DIMENSION(:),ALLOCATABLE:: ck
  COMPLEX(8),DIMENSION(:,:,:,:,:,:),ALLOCATABLE:: tab5
end module tableaux_dft_3d_fast
!
!
module module_energy_luc_fast
  real(8) :: fac(0:50)
  complex(8), parameter :: zeroc=cmplx(0.d0,0.d0)
  private
  public :: energy_luc_fast
contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine energy_luc_fast (ff,df)
  use iso_c_binding
  use module_grid, only: grid
  use module_solvent, only: solvent
  use module_wigner_d, only: wigner_big_D, wigner_small_d
  implicit none
  include 'fftw3.f03'
  integer, parameter :: dp=c_double
  real(dp), intent(out) :: ff
  real(dp), intent(inout) :: df(:,:,:,:,:)
  type(c_ptr) :: plan_fft_c2c_3d_signe_plus, plan_fft_c2c_3d_signe_minus
  complex(dp), allocatable :: in(:,:,:)
  real :: time(0:20)
  integer :: io, no, nx, ny, nz, np, iqx, iqy, iqz, itheta, iphi, ipsi, mmax, mmax2, m, mup, mu, mu2, ix, iy ,iz, i, j, k, khi, m2
  complex(dp), parameter :: zeroc=cmplx(0,0,dp)
  complex(dp), allocatable :: delta_rho_k_angle(:,:,:,:), gamma_k_angle(:,:,:,:)
  complex(dp), allocatable :: delta_rho_k_proj(:,:,:)
  complex(dp), allocatable :: gamma_k_proj(:,:,:) ! m mup mu2
  complex(dp), allocatable :: gamma_angle(:,:,:,:) ! omega, x y z
  real(dp) :: q(3), rho0
  complex(dp) :: a, b
  integer :: imqx, imqy, imqz, ntheta, nphi, npsi
  real(dp) :: fm(0:grid%mmax)
  complex(dp), allocatable :: in2d(:,:), delta_rho_theta_mup_mu(:,:,:)
  type(c_ptr) :: planc2c2dplus, planc2c2dminus

  ff=0
  no = grid%no
  nx = grid%nx
  ny = grid%ny
  nz = grid%nz
  np = grid%np
  rho0 = solvent(1)%rho0
  mmax=grid%mmax
  mmax2=mmax/grid%molrotsymorder
  ntheta=grid%ntheta
  nphi=grid%nphi
  npsi=grid%npsi
  fm(0:mmax) = [( sqrt(real(2*m+1,dp))  ,m=0,mmax  )]
  call cpu_time(time(1))


allocate ( in(nx,ny,nz), source=zeroc)
call dfftw_plan_dft_3d (plan_fft_c2c_3d_signe_plus,  nx, ny, nz, in, in, FFTW_BACKWARD, FFTW_ESTIMATE) ! le tag FFTW_BACKWARD indique un signe + dans l'exponentiel
call dfftw_plan_dft_3d (plan_fft_c2c_3d_signe_minus, nx, ny, nz, in, in, FFTW_FORWARD, FFTW_ESTIMATE) ! le tag FFTW_FORWARD indique un signe - dans l'exponentiel

!
!
! 1. FOURIER TRANSFORM DELTA RHO
!
!
! call random_number( solvent(1)%xi )
! solvent(1)%xi = solvent(1)%xi * 10

allocate ( delta_rho_k_angle(no,nx,ny,nz) , source=zeroc)

do io=1,no
  in = cmplx(solvent(1)%xi(io,:,:,:)**2*rho0-rho0,0,dp)  ! Δρ(Ω,x,y,z)
  call dfftw_execute_dft( plan_fft_c2c_3d_signe_plus, in, in )
  delta_rho_k_angle(io,:,:,:) = in  ! Δρ(Ω,qx,qy,qz)
end do




allocate (delta_rho_k_proj(0:mmax,-mmax:mmax,-mmax2:mmax2) , source=zeroc)
allocate (gamma_k_angle   (no,nx,ny,nz) , source=zeroc)
allocate (gamma_k_proj    (0:mmax,-mmax:mmax,-mmax2:mmax2) , source=zeroc)


!
!
! FOR A GIVEN Q
!
!
do iqz=1,nz
  print*, "... avancement ...",real(iqz)/real(nz)*100.,"%"
  do iqy=1,ny
    do iqx=1,nx

q = [grid%kx(iqx), grid%ky(iqy), grid%kz(iqz)]

!
!
! PASSAGE EN PROJECTIONS
!
!
  if(.not.allocated(in2d)) then
    allocate( in2d(nphi,npsi), source=zeroc)
    allocate( delta_rho_theta_mup_mu(ntheta,-mmax:mmax,-mmax2:mmax2), source=zeroc)
    call dfftw_plan_dft_2d (planc2c2dplus, nphi, npsi, in2d, in2d, FFTW_BACKWARD, FFTW_ESTIMATE) ! le tag FFTW_BACKWARD indique un signe + dans l'exponentiel
    call dfftw_plan_dft_2d (planc2c2dminus, nphi, npsi, in2d, in2d, FFTW_FORWARD, FFTW_ESTIMATE) ! le tag FFTW_BACKWARD indique un signe - dans l'exponentiel
  end if
  do itheta=1,ntheta
    in2d = zeroc
    do concurrent( iphi=1:nphi, ipsi=1:npsi )
      in2d(iphi,ipsi) = delta_rho_k_angle(grid%io(itheta,iphi,ipsi), iqx,iqy,iqz) /real(nphi*npsi)
    end do
    call dfftw_execute (planc2c2dplus, in2d, in2d)
    do iphi=1,nphi
      if (iphi <= nphi/2+1) then
        mup = iphi-1
      else
        mup = iphi-1-nphi
      end if
      do ipsi=1,npsi
        if(ipsi<=npsi/2+1) then
          mu2=ipsi-1
        else
          mu2=ipsi-1-npsi
        end if
        delta_rho_theta_mup_mu(itheta,mup,mu2) = in2d(iphi,ipsi)
      end do
    end do
  end do ! sur theta
  do m=0,mmax
    m2 = m/2
    do mup=-m,m
      do mu2=-m2,m2
        mu=2*mu2
        do itheta=1,ntheta
delta_rho_k_proj(m,mup,mu2) = delta_rho_k_proj(m,mup,mu2) + delta_rho_theta_mup_mu(itheta,mup,mu2)&
         * wigner_small_d(m,mup,mu,grid%thetaofntheta(itheta)) * grid%wthetaofntheta(itheta) &
         /sum(grid%wthetaofntheta) * fm(m)
        end do
      end do
    end do
  end do



!
!
!   PASSAGE REPERE LIE A Q
! + OZ
! + REPASSAGE REPERE FIXE
!
!
! On garde delta_rho_k_proj_full et gamma_k_proj_full qui sont des tableaux énormes mais pour vérifier qu'on a bien la symétrie
! de l'équation 1.15 de Luc, c'est à dire entre delta_rho_alpha(q) et delta_rho_alpha'(-q)
!
call luc_oz (q, delta_rho_k_proj, gamma_k_proj)


!
!
! RETOUR EN ANGLES
!
!
gamma_k_angle(:,iqx,iqy,iqz) = (0,0)
do io=1,no
  do m=0,mmax
    do mup=-m,m
      do mu2=-m/2,m/2
        mu=mu2*2
        gamma_k_angle(io,iqx,iqy,iqz) = gamma_k_angle(io,iqx,iqy,iqz)   &
         + sqrt(real(2*m+1,dp)) * wigner_big_D(m,mup,mu,grid%theta(io),grid%phi(io),grid%psi(io)) * gamma_k_proj(m,mup,mu2)
      end do
    end do
  end do
end do

end do ! iqx
end do ! iqy
end do ! iqz





deallocate (delta_rho_k_proj)
deallocate (gamma_k_proj )

!
!
! INVERSE FOURIER TRANSFORM
!
!
if(.not.allocated(gamma_angle)) allocate( gamma_angle(no,nx,ny,nz))
gamma_angle = zeroc
do io=1,no
  in = delta_rho_k_angle(io,:,:,:)
  call dfftw_execute_dft( plan_fft_c2c_3d_signe_minus, in, in)
  gamma_angle(io,:,:,:) = in/real(nx*ny*nz,dp)
end do


call cpu_time (time(3))
deallocate(in)
deallocate(delta_rho_k_angle)

print*, "energy_luc took ",time(3)-time(1),"sec"
print*, "ff luc fast =", ff

end subroutine energy_luc_fast






  subroutine luc_oz (qvec, delta_rho_proj, gamma4_proj)
    USE lecture                                    ! pour lire des valeurs � la Luc
    USE tableaux_dft_3d_fast                            ! d�finit le type des tableaux de valeurs � partager
    use module_rotation, only: thetaofq, phiofq, rotation_matrix_between_complex_spherical_harmonics_lu
    use module_grid,only: grid
    use iso_c_binding, only: dp=>c_double
    implicit real(8) (a-h,o-z)
    real(8),    intent(in) :: qvec(3)
    complex(8), intent(in) :: delta_rho_proj(0:grid%mmax,-grid%mmax:grid%mmax,-grid%mmax/2:grid%mmax/2)
    complex(8), intent(out) :: gamma4_proj(0:grid%mmax,-grid%mmax:grid%mmax,-grid%mmax/2:grid%mmax/2)
    CHARACTER(50)::nomfic,texte
    CHARACTER(10) texte1
    COMPLEX(8)::xi_cmplx,auxi_cmplx,ck_cmplx,coeff_cmplx,coeff1_cmplx,coeff2_cmplx
    integer :: mmax, mmax2
    logical, save :: first_time_here = .true.
    logical :: qz_was_negatif
    fac(0)=1.
    do k=1,50
      fac(k)=fac(k-1)*REAL(k)
    end do
    mnmax=grid%mmax
    mnmax2=mnmax/2
    mmax=mnmax
    mmax2=mnmax2
    nbeta=grid%ntheta
    nphi=grid%nphi
    nomeg=grid%npsi
    pi=4.d0*ATAN(1.d0)
    xi_cmplx=(0.,1.d0)
    !
    !                           coordonnees du vecteur q (en A-1)
    !
    qx=qvec(1)
    qy=qvec(2)
    qz=qvec(3)
    qq=SQRT(qx**2+qy**2+qz**2)

    !           les elements Rllambda0(q)
    !                      les elements Rmmu'khi(q)          (psi(q)=0)


    if(.not.allocated(harsph1_q)) ALLOCATE(harsph1_q(0:mnmax,-mnmax:mnmax,-mnmax:mnmax))
    harsph1_q = rotation_matrix_between_complex_spherical_harmonics_lu (qvec,mnmax)


    !
    !                           on definit ou lit les projections de delta_rho de depart (complexes!)
    !

    !
    !                          on lit les projections mnlmunu de ck
    !
    ! PRINT*, '********************************************************************************'
    ! 130 PRINT*, 'Lecture des cmnlmunu(q) dans un fichier'
    ! PRINT*, 'ce fichier contient apres quelques lignes de baratin la ligne des m puis celle des n...'
    ! PRINT*, 'suivies des lignes de q,cmnlmunu(q) (on choisira le q le plus proche, par exces)'
    ! PRINT*, 'attention: si l pair, fonction reelle; si l impair, fonction imaginaire pure donc i implicite!'
    ! PRINT*, 'Nom du fichier --->'
    dq=0.613592315E-01
    if(.not.allocated(tab5)) then
      OPEN(7,FILE="/home/levesque/Recherche/00__BELLONI/2016__FEVRIER__OZ/ck_h2o39-3916_l.txt",STATUS='old')
      ialpmax=549
      ALLOCATE(mm(ialpmax),nn(ialpmax),ll(ialpmax),mumu(ialpmax),nunu(ialpmax),ck(ialpmax))
      do i=1,5
        READ(7,*) texte
      end do
      READ(7,*) texte1,mm
      READ(7,*) texte1,nn
      READ(7,*) texte1,ll
      READ(7,*) texte1,mumu
      READ(7,*) texte1,nunu
      iqmax = int( sqrt(maxval(grid%kx)**2+maxval(grid%ky)**2+maxval(grid%kz)**2) / dq ) +2 ! +2 is surely useless
      if(iqmax>1024) error stop "mon fichier de ck est trop restreint en valeurs de q for our nfft and L"
      ALLOCATE (tab5(0:mnmax,0:mnmax,-mnmax:mnmax,-mnmax2:mnmax2,-mnmax2:mnmax2,iqmax), source=zeroc)
      do iq=1,iqmax
        READ(7,*) q,ck(:)

          do khi=0,mnmax                ! que khi>=0 pour l'instant
            do ialp=1,ialpmax
              m=mm(ialp); n=nn(ialp); l=ll(ialp); mu=mumu(ialp); nu=nunu(ialp)
              IF(khi>MIN(m,n)) cycle
              mu2=mu/2; nu2=nu/2
              coeff_cmplx=symbol_3j(m,n,l,khi,-khi,0)*ck(ialp)
              IF((-1)**l==-1) coeff_cmplx=xi_cmplx*coeff_cmplx             ! imaginaire pur si l impair
              tab5(m,n,khi,mu2,nu2,iq)=tab5(m,n,khi,mu2,nu2,iq)+coeff_cmplx                     ! tab5 est donc complexe
              IF(mu/=0.or.nu/=0) tab5(m,n,khi,-mu2,-nu2,iq)=tab5(m,n,khi,-mu2,-nu2,iq)+(-1)**(m+n+l)*coeff_cmplx
              IF(m/=n.or.ABS(mu)/=ABS(nu)) then
                tab5(n,m,khi,nu2,mu2,iq)=tab5(n,m,khi,nu2,mu2,iq)+(-1)**(m+n)*coeff_cmplx
                IF(mu/=0.or.nu/=0) tab5(n,m,khi,-nu2,-mu2,iq)=tab5(n,m,khi,-nu2,-mu2,iq)+(-1)**(l)*coeff_cmplx
              endif
            end do                           ! fin ialp
          end do                             ! fin khi>=0
          !                            et je complete les khi<0 avec cmnmunu-khi=(-1)**(m+n)*cmnmunukhi*
          do m=0,mnmax
            m2=m/2
            do n=0,mnmax
              n2=n/2
              coeff=(-1)**(m+n)
              do khi=-MIN(m,n),-1
                tab5(m,n,khi,-m2:m2,-n2:n2,iq)=coeff*CONJG(tab5(m,n,-khi,-m2:m2,-n2:n2,iq))
              end do
            end do
          end do

      end do
      CLOSE(7)
      deallocate(ck)
    end if

    !
    !
    !            METHODE 4: projections mais en passant par le repere local
    !
    !
    !     d'abord, passer des projections delta_rho a delta_rho'
    if(.not.allocated(delta_rho_proj1)) ALLOCATE(delta_rho_proj1(0:mnmax,-mnmax:mnmax,-mnmax2:mnmax2))
    ! PRINT*, 'deltarho_proj'''
    delta_rho_proj1=0.
    do m=0,mnmax
      m2=m/2
      do khi=-m,m
        do mu2=-m2,m2
          delta_rho_proj1(m,khi,mu2)=SUM(delta_rho_proj(m,-m:m,mu2)*harsph1_q(m,-m:m,khi))
        end do
      end do
    end do
    !               rappel: cmnmunu;khi est dans tab5
    !               faire alors OZ le plus simple!
    iq=int(qq/dq +0.5) +1

    if(.not.allocated(gamma4_proj1)) ALLOCATE(gamma4_proj1(0:mnmax,-mnmax:mnmax,-mnmax2:mnmax2))
    gamma4_proj1=0.
    do khi=-mnmax,mnmax
      do mu2=-mnmax2,mnmax2
        mu=2*mu2
        do m=MAX(ABS(khi),ABS(mu)),mnmax

          do nu2=-mnmax2,mnmax2
            nu=2*nu2
            do n=MAX(ABS(khi),ABS(nu)),mnmax
              gamma4_proj1(m,khi,mu2)=gamma4_proj1(m,khi,mu2)+(-1)**khi *tab5(m,n,khi,mu2,nu2,iq) *delta_rho_proj1(n,khi,-nu2)
            end do
          end do

        end do
      end do
    end do
    !             et revenir des projections gamma' aux projections gamma
    gamma4_proj=cmplx(0.,0)
    do m=0,mnmax
      m2=m/2
      do mu1=-m,m
        do mu2=-m2,m2
          gamma4_proj(m,mu1,mu2)=SUM(gamma4_proj1(m,-m:m,mu2)*CONJG(harsph1_q(m,mu1,-m:m)))
        end do
      end do
    end do
  end subroutine










  !
  !
  !
  !
  function symbol_3j(m,n,l,mu,nu,lu)
    implicit REAL(8) (a-h,o-z)
    !
    !        symbole 3j
    !        Messiah page 910 eq.21
    !
    IF(itriangle(m,n,l).eq.0.or.mu+nu+lu.NE.0.or.                                     &
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
function itriangle(m,n,l)
  !        nul sauf si |m-n|<l<m+n
  !        rq: ne depend pas de l'ordre des 3 entiers
  itriangle=0
  IF(l.ge.ABS(m-n).and.l.le.m+n) itriangle=1
end
!
function delta(m,n,l)
  implicit REAL(8) (a-h,o-z)
  delta=fac(m+n-l)*fac(n+l-m)*fac(l+m-n)/fac(m+n+l+1)
end
!
function harm_sph(m,mu,mup,beta)
  implicit real(8) (a-h,o-z)
  !
  !          pour harmonique spherique Rm,mu,mup(Omega=omega,beta,phi)
  !                                    =exp(-i*omega)*r(m)mu,mup*exp(-i*phi)
  !          calcul l'element mu,mup de la matrice r(m) en fonction de l'angle beta
  !          formule de Wigner dans Messiah eq.72 p922 betement
  !          luc72p143
  !          si mu ou mup nul, formule de recurrence stable
  !
  IF(ABS(mu)>m.or.ABS(mup)>m) THEN; harm_sph=0.; RETURN; endif
    !
    IF(mu==0.or.mup==0) THEN                      ! mu=0 ou mup=0
      !
      mu0=mu; beta0=beta
      IF(mu==0) THEN; mu0=mup; beta0=-beta; ENDIF       ! je mets le 0 en second
        !          si mu negatif, ca vaut (-1)**mu * valeur pour -mu
        x=1.d0;
        IF(mu0<0) THEN; x=(-1)**mu0; mu0=-mu0; endif
          cc=COS(beta0)                            ! plut�t a partir d'une formule de recurrence stable
          pm1=0.                                  ! des polynomes de Legendre associes Pl,m
          pm=1.d0                                 ! luc73p96
          do l=mu0+1,m
            pm2=pm1
            pm1=pm
            pm=(cc*dble(2*l-1)*pm1-dble(l+mu0-1)*pm2)/DBLE(l-mu0)
          end do
          harm_sph=x*(-1)**mu0*SQRT(fac(m-mu0)/fac(m+mu0))*fac(2*mu0)/(2.d0**mu0*fac(mu0))*SIN(beta0)**mu0*pm
          !
        else                  !   donc mu et mup non nuls, utiliser betement la formule de Wigner
          !
          harm_sph=0.
          cc=COS(0.5d0*beta)
          ss=SIN(0.5d0*beta)
          do it=MAX(0,mu-mup),MIN(m+mu,m-mup)
            harm_sph=harm_sph+(-1)**it/(fac(m+mu-it)*fac(m-mup-it)*fac(it)*fac(it-mu+mup))*  &
            cc**(2*m+mu-mup-2*it)*ss**(2*it-mu+mup)
          end do
          harm_sph=SQRT(fac(m+mu)*fac(m-mu)*fac(m+mup)*fac(m-mup))*harm_sph
        endif
      end
      !
      subroutine gauleg(x,w,n)
        implicit REAL(8) (a-h,o-z)
        DIMENSION x(n),w(n)
        !
        !      calcule les abscisses x() (sur -1,1) et poids () de la quadrature gauss-legendre pour n points
        !      luc74p85
        !
        pi=4.d0*ATAN(1.d0)
        !
        m=(n+1)/2                       ! racines symetriques par rapport a 0
        do i=1,m                        ! on s'interesse au i�me zero du polynome Pn(x) de Legendre
          xi=COS(pi*(i-0.25)/(n+0.5))     ! estimation de d�part qu'on va raffiner par NR
          100 p1=1.d0
          p2=0.
          do j=1,n
            p3=p2
            p2=p1
            p1=((2*j-1.d0)*xi*p2-(j-1.d0)*p3)/j    ! relation de r�currence entre les Pj
          end do
          pp=n*(xi*p1-p2)/(xi**2-1.d0)            ! donne Pn' en fonction de Pn et Pn-1
          deltaxi=-p1/pp                          ! NR
          xi=xi+deltaxi
          IF(ABS(deltaxi)>1.d-13) GOTO 100
          x(i)=xi
          w(i)=1.d0/((1.d0-xi**2)*pp**2)         ! poids normalise a 1
          x(n+1-i)=-xi
          w(n+1-i)=w(i)
        end do
      end
      !


    end module

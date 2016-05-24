
!     Last change:  LB    1 Feb 2016   12:51 pm
!
module tableaux_dft_3d
  use precision_kinds, only: dp
  implicit none
  INTEGER::mnmax,mnmax2,nbeta,nphi,nomeg
  real(dp),DIMENSION(:),ALLOCATABLE::beta,cosbeta,sinbeta,wb
  real(dp),DIMENSION(:),ALLOCATABLE::phi,cosphi,sinphi
  complex(dp),DIMENSION(:),ALLOCATABLE::expiphi
  complex(dp),DIMENSION(:,:),ALLOCATABLE::expikhiphi
  real(dp),DIMENSION(:),ALLOCATABLE::omeg
  complex(dp),DIMENSION(:),ALLOCATABLE::expiomeg
  complex(dp),DIMENSION(:,:),ALLOCATABLE::expimuomeg
  real(dp),DIMENSION(:,:,:,:),ALLOCATABLE::harsph
  complex(dp),DIMENSION(:,:),ALLOCATABLE::harsph_q
  complex(dp),DIMENSION(:,:,:),ALLOCATABLE::harsph1_q
  complex(dp),DIMENSION(:,:,:),ALLOCATABLE::delta_rho_proj,delta_rho,delta_rho_proj1,auxi_1,auxi_2
  complex(dp),DIMENSION(:,:,:),ALLOCATABLE::gamma1,gamma1_proj,gamma2,gamma2_proj,gamma3_proj,gamma4_proj,gamma4_proj1
  INTEGER, DIMENSION(:),ALLOCATABLE:: mm,nn,ll,mumu,nunu
  real(dp), DIMENSION(:),ALLOCATABLE:: ck
  real(dp), allocatable :: my_ck(:,:)
  ! complex(dp), DIMENSION(:,:,:,:,:,:),ALLOCATABLE:: ck_omega_omega
  complex(dp),DIMENSION(:,:,:,:,:),ALLOCATABLE:: tab5
  ! complex(dp),DIMENSION(:,:,:,:,:),ALLOCATABLE:: tab6
  ! complex(dp),DIMENSION(:,:,:,:),ALLOCATABLE:: tab7
  ! complex(dp),DIMENSION(:,:,:,:),ALLOCATABLE:: tab4
  !complex(dp),DIMENSION(:,:,:),ALLOCATABLE:: tab3
end module tableaux_dft_3d
!
!
module module_energy_luc
  use precision_kinds, only: dp
  real(dp), allocatable :: fac(:)

contains



subroutine energy_luc (ff,df)
  use iso_c_binding, only: C_PTR, C_INT, C_INT32_T, C_INTPTR_T, C_DOUBLE_COMPLEX, C_DOUBLE, C_FUNPTR, C_SIZE_T, C_FLOAT, &
                           C_FLOAT_COMPLEX, C_CHAR
  use precision_kinds, only: dp
  use module_grid, only: grid
  use module_solvent, only: solvent
  use module_wigner_d, only: wigner_big_D
  implicit none
  include 'fftw3.f03'
  real(dp), intent(out) :: ff
  ! real(dp), intent(inout) :: df(:,:,:,:,:) ! (grid%no, grid%nx, grid%ny, grid%nz, solvent(1)%nspec)
  real(dp), intent(inout) :: df(grid%no, grid%nx, grid%ny, grid%nz, solvent(1)%nspec)
  type(c_ptr) :: plan_fft_c2c_3d_signe_plus, plan_fft_c2c_3d_signe_minus
  complex(dp), allocatable :: in(:,:,:), out(:,:,:)
  real :: time(0:20)
  integer :: io, no, nx, ny, nz, np, iqx, iqy, iqz, itheta, iphi, ipsi, mmax, mmax2, m, mup, mu, mu2, ix, iy ,iz, i, j, k, khi, m2
  complex(dp), parameter :: zeroc=cmplx(0,0,dp)
  complex(dp), allocatable :: delta_rho_k_angle(:,:,:,:), gamma_k_angle(:,:,:,:)
  complex(dp), allocatable :: gamma_full(:,:,:,:) ! gamma(Ω,x,y,z)
  complex(dp), allocatable :: delta_rho_k_proj(:,:,:)
  complex(dp), allocatable :: gamma_k_proj(:,:,:) ! m mup mu2
  complex(dp), allocatable :: gamma_angle(:,:,:,:) ! omega, x y z
  real(dp) :: q(3), rho0
  complex(dp), allocatable :: delta_rho_k_proj_full(:,:,:,:,:,:)
  complex(dp), allocatable :: gamma_k_proj_full(:,:,:,:,:,:)
  complex(dp), allocatable :: delta_rho_prime_k_proj_full(:,:,:,:,:,:)
  complex(dp) :: a, b
  integer :: imqx, imqy, imqz
  logical :: erreur_trouvee

  ff=0
  no = grid%no
  nx = grid%nx
  ny = grid%ny
  nz = grid%nz
  np = grid%np
  rho0 = solvent(1)%rho0
  mmax=grid%mmax
  mmax2=mmax/2
  call cpu_time(time(1))


  allocate ( in(nx,ny,nz), source=zeroc)
  allocate ( out(nx,ny,nz), source=zeroc)
  select case (dp)
  case(c_double)
    call dfftw_plan_dft_3d (plan_fft_c2c_3d_signe_plus,  nx, ny, nz, in, out, FFTW_BACKWARD, FFTW_ESTIMATE) ! le tag FFTW_BACKWARD indique un signe + dans l'exponentiel
    call dfftw_plan_dft_3d (plan_fft_c2c_3d_signe_minus, nx, ny, nz, in, out, FFTW_FORWARD, FFTW_ESTIMATE) ! le tag FFTW_FORWARD indique un signe - dans l'exponentiel
  case(c_float)
    call sfftw_plan_dft_3d (plan_fft_c2c_3d_signe_plus,  nx, ny, nz, in, out, FFTW_BACKWARD, FFTW_ESTIMATE) ! le tag FFTW_BACKWARD indique un signe + dans l'exponentiel
    call sfftw_plan_dft_3d (plan_fft_c2c_3d_signe_minus, nx, ny, nz, in, out, FFTW_FORWARD, FFTW_ESTIMATE) ! le tag FFTW_FORWARD indique un signe - dans l'exponentiel
  end select

!
!
! 1. FOURIER TRANSFORM DELTA RHO
!
!

allocate ( delta_rho_k_angle(no,nx,ny,nz) , source=zeroc)
allocate ( gamma_k_angle(no,nx,ny,nz) , source=zeroc)

do io=1,no
  in = cmplx(solvent(1)%xi(io,:,:,:)**2*rho0-rho0,0,dp)  ! Δρ(Ω,x,y,z)
  select case (dp)
  case(c_double)
    call dfftw_execute_dft( plan_fft_c2c_3d_signe_plus, in, out )
  case(c_float)
    call sfftw_execute_dft( plan_fft_c2c_3d_signe_plus, in, out )
  end select
  delta_rho_k_angle(io,:,:,:) = out  ! Δρ(Ω,qx,qy,qz)
end do


!
! on s'assure :
! - la symétrie hermitienne est gardée
! - surtout, on sait aller chercher q et -q dans les tableaux pour tous les q. Il n'y a pas d'exception, que nx, ny, nz soient pairs ou impairs.
! Autrement dit, ce n'est pas vraiment la symétrie hermitienne qu'on teste, mais la fonction grid%ix_mq(ixq) qui a l'indice iqx trouve l'indice ix_mq correspondant à -q(ix_q)
!
block
real(dp) :: maxdiff, diff
maxdiff=0._dp
erreur_trouvee=.false.
do concurrent (iqx=1:nx, iqy=1:ny, iqz=1:nz, io=1:no)
  imqx = grid%ix_mq(iqx)
  imqy = grid%iy_mq(iqy)
  imqz = grid%iz_mq(iqz)
  a = delta_rho_k_angle(io,iqx,iqy,iqz)
  b = conjg(delta_rho_k_angle(io,imqx,imqy,imqz))
  diff=abs(a-b)
  if (diff>1.D-10) erreur_trouvee=.true.
  if(diff>maxdiff) maxdiff=diff
end do
if(erreur_trouvee) error stop "raté pour la symétrie hermitienne de delta_rho_k_angle"
print*,"maxdiff=",maxdiff
end block





allocate (delta_rho_k_proj(0:mmax,-mmax:mmax,-mmax2:mmax2) , source=zeroc)
allocate (gamma_k_proj    (0:mmax,-mmax:mmax,-mmax2:mmax2) , source=zeroc)

allocate (delta_rho_k_proj_full(0:mmax,-mmax:mmax,-mmax2:mmax2,nx,ny,nz) , source=zeroc)
allocate (gamma_k_proj_full(0:mmax,-mmax:mmax,-mmax2:mmax2,nx,ny,nz) , source=zeroc)

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
delta_rho_k_proj = cmplx(0,0)
do m=0,mmax
  do mup=-m,m
    do mu2=-m/2,m/2
      mu=2*mu2
      delta_rho_k_proj(m,mup,mu2) = sqrt(real(2*m+1,dp)) /sum(grid%w) *sum( delta_rho_k_angle(:,iqx,iqy,iqz) *grid%w *  &
        [( conjg(wigner_big_D(m,mup,mu,grid%theta(io),grid%phi(io),grid%psi(io)))  ,io=1,no )]   )
    end do
  end do
end do
! put delta_rho_k_proj into an array that stores all projections for all q. Utile seulement pour les tests
delta_rho_k_proj_full(0:mmax,-mmax:mmax,-mmax2:mmax2,iqx,iqy,iqz) = delta_rho_k_proj(:,:,:)


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
gamma_k_proj_full(0:mmax,-mmax:mmax,-mmax2:mmax2,iqx,iqy,iqz) = gamma_k_proj(0:mmax,-mmax:mmax,-mmax2:mmax2)


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






! est ce qu'on a bien la symetrie hermitienne sur gamma_k_angle ?

erreur_trouvee = .false.
do concurrent (iqx=1:nx, iqy=1:ny, iqz=1:nz, io=1:no)
  imqx = grid%ix_mq(iqx)
  imqy = grid%iy_mq(iqy)
  imqz = grid%iz_mq(iqz)
  a = gamma_k_angle(io,iqx,iqy,iqz)
  b = conjg(gamma_k_angle(io,imqx,imqy,imqz))
  if (abs(a-b)>1.D-10) then
      erreur_trouvee = .true.
      print*, iqx,iqy,iqz,io,a,b
  end if
end do
if (erreur_trouvee)  error stop "on n'a pas la symetrie hermitienne sur gamma_k_angle"

!
!
! ! on teste 1.15
erreur_trouvee = .false.
do concurrent( iqx=1:nx, iqy=1:ny, iqz=1:nz )
  imqx = grid%ix_mq(iqx)
  imqy = grid%iy_mq(iqy)
  imqz = grid%iz_mq(iqz)
  do m=0,mmax
    do mup=-m,m
      do mu2=-m/2,m/2
        mu = 2*mu2
        a =        gamma_k_proj_full(m,-mup,-mu2,iqx,iqy,iqz)
        b =  (-1)**(mup+mu)*conjg(gamma_k_proj_full(m, mup, mu2,imqx,imqy,imqz))
        if (abs(a-b)>1.D-10) erreur_trouvee=.true.
      end do
    end do
  end do
end do
if (erreur_trouvee) error stop "1.15 pas verifiee sur gamma_k_proj_full"




! ! ! on teste 1.25 sur delta_rho_k_proj_prime
! erreur_trouvee = .false.
! do concurrent( iqx=1:nx, iqy=1:ny, iqz=1:nz )
!   imqx = grid%ix_mq(iqx)
!   imqy = grid%iy_mq(iqy)
!   imqz = grid%iz_mq(iqz)
!   do m=0,mmax
!     do khi=-m,m
!       do mu2=-m/2,m/2
!         mu = 2*mu2
!         a = gamma_k_proj_full(m,khi,-mu2,iqx,iqy,iqz)
!         b = conjg(gamma_k_proj_full(m, khi, mu2,imqx,imqy,imqz)) * (-1)**(m+khi+mu)
!         if (abs(a-b)>1.D-10) then
!           erreur_trouvee=.true.
!           print*,iqx,iqy,iqz, m,khi,mu,a,b
!         end if
!       end do
!     end do
!   end do
! end do
! if(erreur_trouvee) error stop "1.25 pas verifiee pour gamma_k_proj_full"


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
  select case(dp)
  case(c_double)
    call dfftw_execute_dft( plan_fft_c2c_3d_signe_minus, in, out)
  case(c_float)
    call sfftw_execute_dft( plan_fft_c2c_3d_signe_minus, in, out)
  end select
  gamma_angle(io,:,:,:) = out/real(nx*ny*nz,dp)
end do


! On a maintenant notre gamma(r,omega). On veut testé qu'il soit purement réel.
erreur_trouvee=.false.
do concurrent (iqx=1:nx, iqy=1:ny, iqz=1:nz, io=1:no)
  if (imag(gamma_angle(io,iqx,iqy,iqz))>1.D-10) erreur_trouvee=.true.
end do
if(erreur_trouvee) error stop "gamma(r,omega) n'est pas purement reel"


call cpu_time (time(3))
deallocate(in)
deallocate(delta_rho_k_angle)


block
use module_thermo, only: thermo
real(dp) :: vexc(no), xi, rho, dv, kT
real(dp), parameter :: fourpisq=4*acos(-1._dp)**2
dv = grid%dv
kT = thermo%kBT
do iz=1,nz
  do iy=1,ny
    do ix=1,nx
      vexc(1:no) = -kT*grid%w(1:no)*gamma_angle(:,ix,iy,iz) *fourpisq  /solvent(1)%n0
      do io=1,no
        xi = solvent(1)%xi(io,ix,iy,iz)
        rho = xi**2*rho0
        ff = ff + (rho-rho0)*0.5_dp*vexc(io)*dv
        df(io,ix,iy,iz,1) = df(io,ix,iy,iz,1) + 2._dp*vexc(io)*rho0*xi
      end do
    end do
  end do
end do
end block


print*, "energy_luc took ",time(3)-time(1),"sec"
end subroutine energy_luc
















































  subroutine luc_oz (qvec, my_delta_rho_proj, my_gamma_proj)
    use iso_c_binding, only: c_double, c_float
    USE lecture                                    ! pour lire des valeurs � la Luc
    USE tableaux_dft_3d                            ! d�finit le type des tableaux de valeurs � partager
    use module_grid,only: grid
    use precision_kinds, only: dp
    implicit real(dp) (a-h,o-z)
    real(dp),    intent(in) :: qvec(3)
    complex(dp), intent(in) :: my_delta_rho_proj(0:grid%mmax,-grid%mmax:grid%mmax,-grid%mmax/2:grid%mmax/2)
    complex(dp), intent(out) :: my_gamma_proj(0:grid%mmax,-grid%mmax:grid%mmax,-grid%mmax/2:grid%mmax/2)
    character(50)::nomfic,texte
    character(10) texte1
    complex(dp)::xi_cmplx,auxi_cmplx,ck_cmplx,coeff_cmplx,coeff1_cmplx,coeff2_cmplx
    integer :: mmax
    logical :: qz_was_negatif
    !
    !interface
    !subroutine proj_angl(f_proj,f_ang)
    !!implicit real(dp) (a-h,o-z)
    !complex(dp), DIMENSION(:,:,:):: f_proj,f_ang
    !end subroutine
    !END interface
    !
    !interface
    !subroutine angl_proj(f_ang,f_proj)
    !!implicit real(dp) (a-h,o-z)
    !complex(dp), DIMENSION(:,:,:):: f_proj,f_ang
    !end subroutine
    !END interface
    !
    if(.not.allocated(fac)) then
      select case(dp)
      case(c_double)
        allocate(fac(0:50))
        fac(0)=1.
        do k=1,50
          fac(k)=fac(k-1)*REAL(k,dp)
        end do
      case(c_float)
        allocate(fac(0:34))
        fac(0)=1.
        do k=1,34
          fac(k)=fac(k-1)*REAL(k,dp)
        end do
      end select
    end if

    !
    !               pour test DFT_3D, calcul des projections de gamma a partir de celles de deltarho   (solute-solvant)
    !               dans l'espace de fourier pour un vecteur q donne
    !               4 methodes differentes
    !               luc91p88
    !
    ! 9999 continue
    ! PRINT*, '********************************************************************************'
    ! PRINT*, 'Calcul des projections gam mmu''mu a partir des delta_rho mmu''mu par convolution angulaire avec ck'
    ! PRINT*, 'dans l''espace de fourier, pour un vecteur q particulier'
    ! PRINT*, '(ces projections sont definies dans le repere fixe, z axe principal)'
    ! PRINT*, '4 methodes differentes'
    ! PRINT*, '********************************************************************************'
    mmax=grid%mmax
    mmax2=grid%mmax/2
    mnmax=grid%mmax
    mnmax2=mnmax/2
    ! PRINT*, 'Valeur de nmax --->', mnmax
    nbeta=grid%ntheta
    nphi=grid%nphi
    nomeg=grid%npsi
    ! PRINT*, 'Nombre d''angles theta, phi et psi --->', nbeta, nphi, nomeg
    !
    !                       on calcule les harmoniques spheriques generalisees
    !
    pi=4._dp*ATAN(1._dp)
    xi_cmplx=(0.,1._dp)
    !                       theta (ou beta)
    if(.not.allocated(beta)) then
      ALLOCATE(beta(nbeta),cosbeta(nbeta),sinbeta(nbeta),wb(nbeta))
      call gauleg(cosbeta,wb,nbeta)               ! cosbeta()=racines de Pnbeta(x), wb() le poids
      beta=ACOS(cosbeta)
      sinbeta=SIN(beta)
    end if

    if(.not.allocated(phi)) then
      ALLOCATE(phi(nphi),cosphi(nphi),sinphi(nphi),expiphi(nphi))
      phi=2._dp*pi* (/(i-0.5,i=1,nphi)/) /REAL(nphi)               ! sur 0-2pi
      ! print*,"phi de luc",phi
      ! print*,"phi de max",grid%phiofnphi
      ! PHI=GRID%PHIOFNPHI
      expiphi=EXP(xi_cmplx*phi)
      cosphi=REAL(expiphi)
      sinphi=AIMAG(expiphi)
    end if

    if(.not.allocated(expikhiphi)) then
      ALLOCATE(expikhiphi(-mnmax:mnmax,nphi))
      do khi=-mnmax,mnmax
        expikhiphi(khi,:)=EXP(xi_cmplx*khi*phi(:))
      end do
    end if
    !                        psi (ou omega)
    if(.not.allocated(omeg)) then
      ALLOCATE(omeg(nomeg),expiomeg(nomeg))
      omeg=pi* (/(i-0.5,i=1,nomeg)/) /REAL(nomeg)                  ! sur 0-pi car H2O
      ! print*,"psi de luc",omeg
      ! print*,"psi de max",grid%psiofnpsi
      ! OMEG=GRID%PSIOFNPSI
      expiomeg=EXP(xi_cmplx*omeg)
    end if

    if(.not.allocated(expimuomeg)) then
      ALLOCATE(expimuomeg(-mnmax2:mnmax2,nomeg))     ! et mu pair=2*mu2
      do mu2=-mnmax2,mnmax2
        mu=2*mu2
        expimuomeg(mu2,:)=EXP(xi_cmplx*mu*omeg(:))
      end do
    end if

    !                            rmmu'mu(theta)
    if(.not.allocated(harsph)) then
      ALLOCATE(harsph(0:mnmax,-mnmax:mnmax,-mnmax2:mnmax2,nbeta))
      harsph=0.
      do khi=-mnmax,mnmax
        do mu2=-mnmax2,mnmax2
          mu=2*mu2
          do m=MAX(ABS(khi),ABS(mu)),mnmax
            do ibeta=1,nbeta
              harsph(m,khi,mu2,ibeta)=harm_sph(m,khi,mu,beta(ibeta))
            end do
          end do
        end do
      end do
    end if
    !
    !                           coordonnees du vecteur q (en A-1)
    !
    ! PRINT*, '********************************************************************************'
    qx=qvec(1)
    qy=qvec(2)
    qz=qvec(3)
    qq=SQRT(qx**2+qy**2+qz**2)
    ! PRINT*, 'Coordonnees x,y,z du vecteur q  (en A-1) --->',qx,qy,qz,'ce qui donne q --->',qq




    if(abs(qx)+abs(qy)<=epsilon(1.)) then
      theta_q = 0.
      phi_q = 0.
    else
      c_theta_q=qz/qq                        ! angles theta,phi pour q
      theta_q=ACOS(c_theta_q)
      s_theta_q=SIN(theta_q)
      c_phi_q=qx/qq/s_theta_q
      if(abs(c_phi_q-1)<epsilon(1.)) then
        phi_q=0.
      else if(abs(c_phi_q+1)<epsilon(1.)) then
        phi_q=acos(-1._dp)
      else
        phi_q=ACOS(c_phi_q)
        IF(qy<0.) phi_q=-phi_q
      end if
      s_phi_q=SIN(phi_q)
    end if

block
  use module_rotation, only: thetaofq, phiofq
    theta_q = thetaofq(qx,qy,qz)
    phi_q   = phiofq(qx,qy,qz)
end block

! block
!   use module_rotation, only: thetaofq, phiofq
!   print*, theta_q,phi_q, thetaofq(qx,qy,qz), phiofq(qx,qy,qz)
! end block


    ! if(abs(qx)+abs(qy)<=epsilon(1.)) then
    !   theta_q = 0.
    !   phi_q = 0.
    ! else
    !   c_theta_q=qz/qq                        ! angles theta,phi pour q
    !   theta_q=ACOS(c_theta_q)
    !   s_theta_q=SIN(theta_q)
    !   c_phi_q=qx/qq/s_theta_q
    !   if(abs(c_phi_q-1)<epsilon(1.)) then
    !     phi_q=0.
    !   else if(abs(c_phi_q+1)<epsilon(1.)) then
    !     phi_q=acos(-1._dp)
    !   else
    !     phi_q=ACOS(c_phi_q)
    !     IF(qy<0.) phi_q=-phi_q  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    !   end if
    !   s_phi_q=SIN(phi_q)
    ! end if


    !           les elements Rllambda0(q)
    if(.not.allocated(harsph_q)) ALLOCATE(harsph_q(0:2*mnmax,-2*mnmax:2*mnmax))
    harsph_q=0.
    do l=0,2*mnmax
      do lambda=-l,l
        harsph_q(l,lambda)=harm_sph(l,lambda,0,theta_q)*EXP(-xi_cmplx*lambda*phi_q)
      end do
    end do
    if(any(harsph_q/=harsph_q)) error stop "problem detected in harsph_q"
    !                      les elements Rmmu'khi(q)          (psi(q)=0)
    if(.not.allocated(harsph1_q)) ALLOCATE(harsph1_q(0:mnmax,-mnmax:mnmax,-mnmax:mnmax))
    harsph1_q=0.
    do m=0,mnmax
      do mu1=-m,m
        do khi=-m,m
          harsph1_q(m,mu1,khi)=harm_sph(m,mu1,khi,theta_q)*EXP(-xi_cmplx*mu1*phi_q)
        end do
      end do
    end do
    if(any(harsph1_q/=harsph1_q)) error stop "problem detected in harsph1_q"

    !          test des R*R
    ! do khi=-mnmax,mnmax
    !   do m=0,mnmax
    !     do n=0,mnmax
    !       do l=-ABS(m-n),m+n
    !         do nu1=-n,n
    !           coeff_cmplx=0.
    !           do mu1=-m,m
    !             do lambda1=-l,l
    !               print*, coeff_cmplx, symbol_3j(m,n,l,mu1,nu1,lambda1),harsph_q(l,lambda1),harsph1_q(m,mu1,khi)
    !               coeff_cmplx=coeff_cmplx+symbol_3j(m,n,l,mu1,nu1,lambda1)*harsph_q(l,lambda1)*harsph1_q(m,mu1,khi)
    !             end do
    !           end do
    !           coeff1_cmplx=(-1)**(nu1+khi)*symbol_3j(m,n,l,khi,-khi,0)*harsph1_q(n,-nu1,khi)
    !           IF(ABS(coeff1_cmplx-coeff_cmplx)>1.d-10) then
    !             PRINT*, m,n,l,nu1,khi,coeff_cmplx,coeff1_cmplx
    !             error stop "test des R*R raté"
    !           end if
    !         end do
    !       end do
    !     end do
    !   end do
    ! end do
    !          fin de test
    !
    !                           on definit ou lit les projections de delta_rho de depart (complexes!)
    !
    if(.not.allocated(delta_rho_proj)) ALLOCATE(delta_rho_proj(0:mnmax,-mnmax:mnmax,-mnmax2:mnmax2))
    delta_rho_proj=cmplx(0.,0.)

    ! PRINT*, '********************************************************************************'
    i_depart = 3
    ! 110 PRINT*, 'Projections delta_rho mmu''mu de depart 1) analytique  2) fichier 3) passé en argument --->', i_depart
    110 CONTINUE
    IF(i_depart==1) then ! ANALYTIQUE
      do khi=-mnmax,mnmax
        do mu2=-mnmax2,mnmax2
          mu=2*mu2
          do m=MAX(ABS(khi),ABS(mu)),mnmax
            delta_rho_proj(m,khi,mu2)=EXP(-m*2.3)*((1.2+2.3*xi_cmplx)*khi-(2.1+0.6*xi_cmplx)*mu)     ! exemple
          end do
        end do
      end do

    else if (i_depart==2) then ! DEPUIS UN FICHIER
      PRINT*, 'Nom du fichier --->'
      READ(*,*) nomfic
      OPEN(7,FILE=nomfic,STATUS='old',ERR=110)
      do j=1,10000
        READ(7,*,END=120) m,khi,mu,auxi_cmplx
        delta_rho_proj(m,khi,mu/2)=auxi_cmplx
      end do
      120 close(7)

    else if (i_depart==3) then
      delta_rho_proj = my_delta_rho_proj
    else
      goto 110
    endif
    !PRINT*, 'deltarho_proj',delta_rho_proj
    !
    !                          on lit les projections mnlmunu de ck
    !
    ! PRINT*, '********************************************************************************'
    ! 130 PRINT*, 'Lecture des cmnlmunu(q) dans un fichier'
    ! PRINT*, 'ce fichier contient apres quelques lignes de baratin la ligne des m puis celle des n...'
    ! PRINT*, 'suivies des lignes de q,cmnlmunu(q) (on choisira le q le plus proche, par exces)'
    ! PRINT*, 'attention: si l pair, fonction reelle; si l impair, fonction imaginaire pure donc i implicite!'
    ! PRINT*, 'Nom du fichier --->'
    OPEN(7,FILE="/home/levesque/Recherche/00__BELLONI/2016__FEVRIER__OZ/ck_h2o39-3916_l.txt",STATUS='old')
    n_baratin=5
    ! PRINT*, 'nombre de lignes de baratin --->', n_baratin
    ialpmax=549
    ! PRINT*, 'Nombre de projections --->', ialpmax
    if(.not.allocated(ck)) ALLOCATE(mm(ialpmax),nn(ialpmax),ll(ialpmax),mumu(ialpmax),nunu(ialpmax),ck(ialpmax))
    do i=1,n_baratin
      READ(7,*) texte
      ! PRINT*, texte
    end do
    READ(7,*) texte1,mm
    ! PRINT*, texte1,mm(1:10)
    READ(7,*) texte1,nn
    ! PRINT*, texte1,nn(1:10)
    READ(7,*) texte1,ll
    ! PRINT*, texte1,ll(1:10)
    READ(7,*) texte1,mumu
    ! PRINT*, texte1,mumu(1:10)
    READ(7,*) texte1,nunu
    ! PRINT*, texte1,nunu(1:10)
    do i=1,10000
      READ(7,*) q,ck
      IF(q>=qq) exit
    end do
    ! PRINT*, 'valeur de q choisie dans le tableau: ',q
    ! PRINT*, 'ck= ',ck(1:MIN(10,ialpmax))
    CLOSE(7)
    !   je construis en tableau a 5 entrees cmnlmunu pour tests ulterieurs
    ! if(.not.allocated(tab6)) ALLOCATE (tab6(0:mnmax,0:mnmax,0:2*mnmax,-mnmax2:mnmax2,-mnmax2:mnmax2))
    ! tab6=0.
    ! do ialp=1,ialpmax
    !   m=mm(ialp)
    !   n=nn(ialp)
    !   l=ll(ialp)
    !   mu=mumu(ialp)
    !   nu=nunu(ialp)
    !   mu2=mu/2
    !   nu2=nu/2
    !   coeff_cmplx=ck(ialp)
    !   IF((-1)**l==-1) coeff_cmplx=xi_cmplx*ck(ialp)
    !   tab6(m,n,l,mu2,nu2)=coeff_cmplx
    !   tab6(m,n,l,-mu2,-nu2)=(-1)**(m+n+l)*coeff_cmplx
    !   tab6(n,m,l,nu2,mu2)=(-1)**(m+n)*coeff_cmplx
    !   tab6(n,m,l,-nu2,-mu2)=(-1)**l*coeff_cmplx
    ! END do
            ! PRINT*, '********************************************************************************'
            ! !
            ! !                           on a donc tout ce qui nous faut pour calculer le produit de convolution
            ! !
            ! !
            ! !           methode 1: on integre betement en angles dans le repere fixe avec c calcule directement dans le repere fixe
            ! !
            ! PRINT*, "Methode 1: j'integre en angles dans le repere fixe et c calculé directement dans le repère fixe"
            ! ALLOCATE(auxi_1(nbeta,-mnmax:mnmax,-mnmax2:mnmax2),auxi_2(nbeta,nphi,-mnmax2:mnmax2))
            ! auxi_1=0.; auxi_2=0.
            ! ALLOCATE(delta_rho(nbeta,nphi,nomeg))
            ! ALLOCATE(gamma1(nbeta,nphi,nomeg),gamma1_proj(0:mnmax,-mnmax:mnmax,-mnmax2:mnmax2))
            ! !               	je transforme delta_rho_proj en delta_rho(Omega)
            ! PRINT*, 'deltarho(O)'
            ! call proj_angl(delta_rho_proj,delta_rho)
            ! ! test
            ! call angl_proj(delta_rho,gamma1_proj)
            ! PRINT*, 'test proj_angl_proj '
            ! IF(SUM(ABS(gamma1_proj-delta_rho_proj))>1.d-10) PRINT*, 'attention: maillage insuffisant'
            ! !GOTO 200
            ! ! 			je fabrique c(Omega,Omega')
            ! ALLOCATE(ck_omega_omega(nbeta,nphi,nomeg,nbeta,nphi,nomeg))
            ! PRINT*, 'c(O,O'')'
            ! ck_omega_omega=0.
            ! do ialp=1,ialpmax
            !   m=mm(ialp); n=nn(ialp); l=ll(ialp); mu=mumu(ialp); nu=nunu(ialp)
            !   mu2=mu/2; nu2=nu/2
            !   fmn=SQRT((2._dp*m+1._dp)*(2._dp*n+1._dp))
            !   ck_cmplx=ck(ialp)
            !   IF((-1)**l==-1) ck_cmplx=xi_cmplx*ck_cmplx
            !   do mu1=-m,m
            !     do nu1=-n,n
            !       lambda1=-mu1-nu1
            !       IF(ABS(lambda1)>l) cycle
            !       s3j=symbol_3j(m,n,l,mu1,nu1,lambda1)
            !       IF(s3j==0.) cycle
            !       coeff_cmplx=fmn*ck_cmplx*s3j*harsph_q(l,lambda1)
            !       do i=1,nbeta
            !         do j=1,nphi
            !           do k=1,nomeg
            !             coeff1_cmplx=coeff_cmplx*harsph(m,mu1,mu2,i)*expikhiphi(-mu1,j)*expimuomeg(-mu2,k)
            !             do i1=1,nbeta
            !               do j1=1,nphi
            !                 do k1=1,nomeg
            !                   coeff2_cmplx=coeff1_cmplx*harsph(n,nu1,nu2,i1)*expikhiphi(-nu1,j1)*expimuomeg(-nu2,k1)
            !                   IF(mu/=0.or.nu/=0) coeff2_cmplx=coeff2_cmplx+(-1)**l*CONJG(coeff2_cmplx)
            !                   ck_omega_omega(i,j,k,i1,j1,k1)=ck_omega_omega(i,j,k,i1,j1,k1)+coeff2_cmplx
            !                 IF(m/=n.or.ABS(mu)/=ABS(nu)) ck_omega_omega(i1,j1,k1,i,j,k)=ck_omega_omega(i1,j1,k1,i,j,k)+(-1)**l*coeff2_cmplx
            !                 end do
            !               end do
            !             end do
            !           end do
            !         end do
            !       end do
            !     end do
            !   end do
            ! end do
            ! !     test en partant de tab6 pour voir si on retrouve bien le meme ck
            ! GOTO 117                   ! je shunte les 2 tests a venir puisque reussis
            ! PRINT*, 'test en partant de tab6'
            ! do i=1,nbeta
            !   do j=1,nphi
            !     do k=1,nomeg
            !       do i1=1,nbeta
            !         do j1=1,nphi
            !           do k1=1,nomeg
            !             coeff_cmplx=0.
            !             do m=0,mnmax
            !               m2=m/2
            !               do mu2=-m2,m2
            !                 do n=0,mnmax
            !                   n2=n/2
            !                   do nu2=-n2,n2
            !                     do l=ABS(m-n),m+n
            !                       coeff1_cmplx=SQRT((2._dp*m+1._dp)*(2._dp*n+1._dp))*tab6(m,n,l,mu2,nu2)
            !                       do mu1=-m,m
            !                         do nu1=-n,n
            !                           lambda1=-mu1-nu1
            !                           IF(ABS(lambda1)>l) cycle
            !                           s3j=symbol_3j(m,n,l,mu1,nu1,lambda1)
            !                           IF(s3j==0.) cycle
            !                           coeff_cmplx=coeff_cmplx+coeff1_cmplx*s3j*				   	     &
            !                           harsph(m,mu1,mu2,i)*expikhiphi(-mu1,j)*expimuomeg(-mu2,k)*                      &
            !                           harsph(n,nu1,nu2,i1)*expikhiphi(-nu1,j1)*expimuomeg(-nu2,k1)*                   &
            !                           harsph_q(l,lambda1)
            !                         end do
            !                       end do
            !                     end do
            !                   end do
            !                 end do
            !               end do
            !             end do
            !     IF(ABS(coeff_cmplx-ck_omega_omega(i,j,k,i1,j1,k1))>1.d-10) PRINT*, i,j,k,i1,j1,k1,coeff_cmplx,ck_omega_omega(i,j,k,i1,j1,k1)
            !           end do
            !         end do
            !       end do
            !     end do
            !   end do
            ! end do
            ! !          fin de test, reussi
            ! !     test de retour
            ! print*, 'test retour de c(o,o'') aux proj'
            ! do ialp=1,ialpmax
            !   m=mm(ialp); n=nn(ialp); l=ll(ialp); mu=mumu(ialp); nu=nunu(ialp)
            !   mu2=mu/2; nu2=nu/2
            !   fmn=SQRT((2._dp*m+1._dp)*(2._dp*n+1._dp))
            !   coeff_cmplx=0.
            !   do mu1=-m,m
            !     do nu1=-n,n
            !       lambda1=-mu1-nu1
            !       IF(ABS(lambda1)>l) cycle
            !       s3j=symbol_3j(m,n,l,mu1,nu1,lambda1)
            !       IF(s3j==0.) cycle
            !       coeff1_cmplx=0.
            !       do i=1,nbeta
            !         do j=1,nphi
            !           do k=1,nomeg
            !             coeff2_cmplx=0.
            !             do i1=1,nbeta
            !               do j1=1,nphi
            !                 do k1=1,nomeg
            !                   coeff2_cmplx=coeff2_cmplx+ck_omega_omega(i,j,k,i1,j1,k1)*	&
            !                   wb(i1)*CONJG(harsph(n,nu1,nu2,i1)*expikhiphi(-nu1,j1)*expimuomeg(-nu2,k1))
            !                 end do
            !               end do
            !             end do
            !             coeff1_cmplx=coeff1_cmplx+coeff2_cmplx*wb(i)*CONJG(harsph(m,mu1,mu2,i)*expikhiphi(-mu1,j)*expimuomeg(-mu2,k))
            !           end do
            !         end do
            !       end do
            !       coeff_cmplx=coeff_cmplx+(2._dp*l+1._dp)*fmn*s3j*CONJG(harsph_q(l,lambda1))*coeff1_cmplx
            !     end do
            !   end do
            !   coeff_cmplx=coeff_cmplx/(nphi*nomeg)**2
            !   ck_cmplx=ck(ialp)
            !   IF((-1)**l==-1) ck_cmplx=xi_cmplx*ck_cmplx
            !   IF(ABS(coeff_cmplx-ck_cmplx)>1.d-10) PRINT*, m,m,l,mu,nu,ck_cmplx,coeff_cmplx
            ! end do
            ! !       fin de test de retour, reussi
            ! 117 continue
            ! !
            ! !      j'integre en omega' par gauss
            ! PRINT*, 'convolution'
            ! gamma1=0.
            ! do i=1,nbeta
            !   do j=1,nphi
            !     do k=1,nomeg
            !       do i1=1,nbeta
            !         do j1=1,nphi
            !           do k1=1,nomeg
            !             gamma1(i,j,k)=gamma1(i,j,k)+wb(i1)*ck_omega_omega(i,j,k,i1,j1,k1)*delta_rho(i1,j1,k1)
            !           end do
            !         end do
            !       end do
            !     end do
            !   end do
            ! end do
            ! gamma1=gamma1/(nphi*nomeg)
            ! !      et je projette
            ! PRINT*, 'gam_proj'
            ! call angl_proj(gamma1,gamma1_proj)
            ! !
            ! !
            ! !          methode 2: j'integre betement en angles dans le repere fixe mais c calcule via le repere local
            ! !
            ! !
            ! 200 continue
            ! PRINT*, '********************************************************************************'
            ! PRINT*, 'Methode 2: j''integre en angles dans le repere fixe MAIS c calcule via le repere local (lie a q)'
            ! !
            ! !          d'abord, je calcule les cmnmunu;khi     mis dans tab5(m,n,khi,mu,nu)
            ! !
            ! if(.not.allocated(tab3)) ALLOCATE (tab3(-mnmax:mnmax,-mnmax2:mnmax2,-mnmax2:mnmax2))
            ! if(.not.allocated(tab4)) ALLOCATE (tab4(0:mnmax,-mnmax:mnmax,-mnmax2:mnmax2,-mnmax2:mnmax2))
            if(.not.allocated(tab5)) ALLOCATE (tab5(0:mnmax,0:mnmax,-mnmax:mnmax,-mnmax2:mnmax2,-mnmax2:mnmax2))
            ! PRINT*, 'ckhi'
            !
            tab5=0._dp
            do khi=0,mnmax                ! que khi>=0 pour l'instant
                do ialp=1,ialpmax
                    m=mm(ialp)
                    n=nn(ialp)
                    l=ll(ialp)
                    mu=mumu(ialp)
                    nu=nunu(ialp)
                    IF(khi>MIN(m,n)) cycle
                    mu2=mu/2; nu2=nu/2
                    coeff_cmplx=symbol_3j(m,n,l,khi,-khi,0)*ck(ialp)
                    IF((-1)**l==-1) coeff_cmplx=xi_cmplx*coeff_cmplx             ! imaginaire pur si l impair
                    tab5(m,n,khi,mu2,nu2)=tab5(m,n,khi,mu2,nu2)+coeff_cmplx                     ! tab5 est donc complexe
                    IF(mu/=0.or.nu/=0) tab5(m,n,khi,-mu2,-nu2)=tab5(m,n,khi,-mu2,-nu2)+(-1)**(m+n+l)*coeff_cmplx
                    IF(m/=n.or.ABS(mu)/=ABS(nu)) then
                        tab5(n,m,khi,nu2,mu2)=tab5(n,m,khi,nu2,mu2)+(-1)**(m+n)*coeff_cmplx
                        IF(mu/=0.or.nu/=0) tab5(n,m,khi,-nu2,-mu2)=tab5(n,m,khi,-nu2,-mu2)+(-1)**(l)*coeff_cmplx
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
                  tab5(m,n,khi,-m2:m2,-n2:n2)=coeff*CONJG(tab5(m,n,-khi,-m2:m2,-n2:n2))
                end do
              end do
            end do

TAB5=CONJG(TAB5)
    !
    !
    !            METHODE 4: projections mais en passant par le repere local
    !
    !
    ! PRINT*, '********************************************************************************'
    ! PRINT*, 'Methode4: je combine les projections via le repere local!'
    if(.not.allocated(gamma4_proj)) ALLOCATE(gamma4_proj(0:mnmax,-mnmax:mnmax,-mnmax2:mnmax2))
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
    ! test aller-retour  reussi donc shunte
    GOTO 417
    PRINT*, 'test retour a deltarho_proj'
    gamma4_proj=0.
    do m=0,mnmax
      m2=m/2
      do mu1=-m,m
        do mu2=-m2,m2
          gamma4_proj(m,mu1,mu2)=SUM(delta_rho_proj1(m,-m:m,mu2)*CONJG(harsph1_q(m,mu1,-m:m)))
        end do
      end do
    end do
    IF(MAXVAL(ABS(gamma4_proj-delta_rho_proj))>1.d-8) PRINT*, 'AR pas bon!'
    !     fin test AR
    417 continue
    !               rappel: cmnmunu;khi est dans tab5
    !               faire alors OZ le plus simple!
    if(.not.allocated(gamma4_proj1)) ALLOCATE(gamma4_proj1(0:mnmax,-mnmax:mnmax,-mnmax2:mnmax2))
    ! PRINT*, 'OZ'
    gamma4_proj1=0.
    do khi=-mnmax,mnmax
      do mu2=-mnmax2,mnmax2
        mu=2*mu2
        do m=MAX(ABS(khi),ABS(mu)),mnmax

          do nu2=-mnmax2,mnmax2
            nu=2*nu2
            do n=MAX(ABS(khi),ABS(nu)),mnmax
              gamma4_proj1(m,khi,mu2)=gamma4_proj1(m,khi,mu2)+(-1)**khi *tab5(m,n,khi,mu2,nu2) *delta_rho_proj1(n,khi,-nu2)
            end do
          end do

        end do
      end do
    end do
    !             et revenir des projections gamma' aux projections gamma
    ! PRINT*, 'gam_proj'
    gamma4_proj=cmplx(0.,0)
    do m=0,mnmax
      m2=m/2
      do mu1=-m,m
        do mu2=-m2,m2
          gamma4_proj(m,mu1,mu2)=SUM(gamma4_proj1(m,-m:m,mu2)*CONJG(harsph1_q(m,mu1,-m:m)))
        end do
      end do
    end do
    my_gamma_proj = gamma4_proj

  end subroutine luc_oz


  !
subroutine proj_angl(f_proj,f_ang)
    USE tableaux_dft_3d, ONLY: mnmax,mnmax2,nbeta,nphi,nomeg,auxi_1,auxi_2,harsph,expikhiphi,expimuomeg
    implicit real(dp) (a-h,o-z)
    !complex(dp), DIMENSION(:,:,:):: f_proj,f_ang
    complex(dp), DIMENSION(0:mnmax,-mnmax:mnmax,-mnmax2:mnmax2), intent(in):: f_proj
    complex(dp), DIMENSION(nbeta,nphi,nomeg), intent(out):: f_ang
    !
    !             passe des projections a une fonction angulaire
    !
    !             d'abord, je transforme m en theta
    auxi_1=0.
    do mu2=-mnmax2,mnmax2
      mu=2*mu2
      do khi=-mnmax,mnmax
        do m=MAX(ABS(khi),ABS(mu)),mnmax
          auxi_1(1:nbeta,khi,mu2)=auxi_1(1:nbeta,khi,mu2)+SQRT(2._dp*m+1._dp)*f_proj(m,khi,mu2)*harsph(m,khi,mu2,1:nbeta)
        end do
      end do
    end do
    !            puis khi en phi
    auxi_2=0.
    do j=1,nphi
      do khi=-mnmax,mnmax
        auxi_2(:,j,:)=auxi_2(:,j,:)+auxi_1(:,khi,:)*expikhiphi(-khi,j)
      end do
    end do
    !            et enfin mu2 en psi
    f_ang=0.
    do k=1,nomeg
      do mu2=-mnmax2,mnmax2
        f_ang(:,:,k)=f_ang(:,:,k)+auxi_2(:,:,mu2)*expimuomeg(-mu2,k)
      end do
    end do
    !
end subroutine proj_angl


subroutine angl_proj(f_ang,f_proj)
    USE tableaux_dft_3d, ONLY: mnmax,mnmax2,nbeta,nphi,nomeg,auxi_1,auxi_2,wb,harsph,expikhiphi,expimuomeg
    implicit real(dp) (a-h,o-z)
    !complex(dp), DIMENSION(:,:,:):: f_proj,f_ang
    complex(dp), DIMENSION(nbeta,nphi,nomeg), intent(in):: f_ang
    complex(dp), DIMENSION(0:mnmax,-mnmax:mnmax,-mnmax2:mnmax2), intent(out):: f_proj
    !
    !             passe de la fonction angulaire aux projections
    !
    !             d'abord, je transforme psi en mu2
    auxi_2=0.
    do mu2=-mnmax2,mnmax2
      do k=1,nomeg
        auxi_2(:,:,mu2)=auxi_2(:,:,mu2)+f_ang(:,:,k)*expimuomeg(mu2,k)
      end do
    end do
    !            puis phi en khi
    auxi_1=0.
    do khi=-mnmax,mnmax
      do j=1,nphi
        auxi_1(:,khi,:)=auxi_1(:,khi,:)+auxi_2(:,j,:)*expikhiphi(khi,j)
      end do
    end do
    auxi_1=auxi_1/(nphi*nomeg)
    !            et enfin theta en m
    f_proj=0.
    do mu2=-mnmax2,mnmax2
      mu=2*mu2
      do khi=-mnmax,mnmax
        do m=MAX(ABS(khi),ABS(mu)),mnmax
          f_proj(m,khi,mu2)=f_proj(m,khi,mu2)+SQRT(2._dp*m+1._dp)*SUM(wb(:)*auxi_1(:,khi,mu2)*harsph(m,khi,mu2,:))
        end do
      end do
    end do
    !
end subroutine angl_proj


pure function symbol_3j(m,n,l,mu,nu,lu)
  implicit real(dp) (a-h,o-z)
  integer, intent(in) :: m, n, l, mu, nu, lu
  real(dp) :: symbol_3j
    !
    !        symbole 3j
    !        Messiah page 910 eq.21
    !
  IF(itriangle(m,n,l).eq.0.or.mu+nu+lu.NE.0.or.                                     &
    ABS(mu).gt.m.or.abs(nu).gt.n.or.abs(lu).gt.l) then
    symbol_3j=0._dp
  else
    som=0._dp
    do it=MAX(0,n-l-mu,m-l+nu),MIN(m+n-l,m-mu,n+nu)
      som=som+(-1)**it/(fac(it)*fac(l-n+it+mu)*fac(l-m+it-nu)*fac(m+n-l-it)*fac(m-it-mu)*fac(n-it+nu))
    end do
    symbol_3j=(-1)**(m-n-lu)*SQRT(delta(m,n,l))*                                       &
    SQRT(fac(m+mu)*fac(m-mu)*fac(n+nu)*fac(n-nu)*fac(l+lu)*fac(l-lu))*som
    IF(mu==0.and.nu==0.and.lu==0.and.2*((m+n+l)/2)/=m+n+l) symbol_3j=0._dp
  endif
end function symbol_3j


pure function itriangle(m,n,l)
  implicit none
  integer :: itriangle
  integer, intent(in) :: m, n, l
  !        nul sauf si |m-n|<l<m+n
  !        rq: ne depend pas de l'ordre des 3 entiers
  itriangle=0
  IF(l.ge.ABS(m-n).and.l.le.m+n) itriangle=1
end function itriangle


pure function delta(m,n,l)
  use precision_kinds, only: dp
  implicit real(dp) (a-h,o-z)
  integer, intent(in) :: m, n, l
  delta=fac(m+n-l)*fac(n+l-m)*fac(l+m-n)/fac(m+n+l+1)
end function delta


pure function harm_sph(m,mu,mup,beta)
  use precision_kinds, only: dp
  implicit real(dp) (a-h,o-z)
  complex(dp) :: harm_sph
  integer, intent(in) :: m, mu, mup
  real(dp), intent(in) :: beta

  !
  !          pour harmonique spherique Rm,mu,mup(Omega=omega,beta,phi)
  !                                    =exp(-i*omega)*r(m)mu,mup*exp(-i*phi)
  !          calcul l'element mu,mup de la matrice r(m) en fonction de l'angle beta
  !          formule de Wigner dans Messiah eq.72 p922 betement
  !          luc72p143
  !          si mu ou mup nul, formule de recurrence stable
  !
  IF(ABS(mu)>m.or.ABS(mup)>m) THEN
    harm_sph=0._dp
    RETURN
  endif
    !
  IF(mu==0.or.mup==0) THEN                      ! mu=0 ou mup=0
      !
    mu0=mu; beta0=beta
    IF(mu==0) THEN
      mu0=mup
      beta0=-beta
    ENDIF       ! je mets le 0 en second
        !          si mu negatif, ca vaut (-1)**mu * valeur pour -mu
    x=1._dp;
    IF(mu0<0) THEN
      x=(-1)**mu0
      mu0=-mu0
    endif
    cc=COS(beta0)                            ! plut�t a partir d'une formule de recurrence stable
    pm1=0._dp                                  ! des polynomes de Legendre associes Pl,m
    pm=1._dp                                 ! luc73p96
    do l=mu0+1,m
      pm2=pm1
      pm1=pm
      pm=(cc*real(2*l-1,dp)*pm1-real(l+mu0-1,dp)*pm2)/real(l-mu0,dp)
    end do
    harm_sph=x*(-1)**mu0*SQRT(fac(m-mu0)/fac(m+mu0))*fac(2*mu0)/(2._dp**mu0*fac(mu0))*SIN(beta0)**mu0*pm
          !
  else                  !   donc mu et mup non nuls, utiliser betement la formule de Wigner
          !
    harm_sph=0.
    cc=COS(0.5_dp*beta)
    ss=SIN(0.5_dp*beta)
    do it=MAX(0,mu-mup),MIN(m+mu,m-mup)
      harm_sph=harm_sph+(-1)**it/(fac(m+mu-it)*fac(m-mup-it)*fac(it)*fac(it-mu+mup))*  &
      cc**(2*m+mu-mup-2*it)*ss**(2*it-mu+mup)
    end do
    harm_sph=SQRT(fac(m+mu)*fac(m-mu)*fac(m+mup)*fac(m-mup))*harm_sph
  endif
end function harm_sph


subroutine gauleg(x,w,n)
    !
    !      calcule les abscisses x() (sur -1,1) et poids () de la quadrature gauss-legendre pour n points
    !      luc74p85
    !
    !
    use precision_kinds, only: dp
    ! implicit real(dp) (a-h,o-z)
    use iso_c_binding, only: c_double, c_float
    implicit none
    integer, intent(in) :: n
    real(dp), intent(out) :: x(n), w(n)
    real(dp), parameter :: pi=acos(-1._dp)
    integer :: m, i, j
    real(dp) :: xi, p1, p2, p3, pp, deltaxi, deltaximax
    m=(n+1)/2                       ! racines symetriques par rapport a 0
    select case(dp)
    case(c_double)
      deltaximax=10._dp**(-13)
    case(c_float)
      deltaximax=10._dp**(-7)
    end select
    do i=1,m                        ! on s'interesse au ieme zero du polynome Pn(x) de Legendre
        xi=COS(pi*(i-0.25_dp)/(n+0.5_dp))     ! estimation de depart qu'on va raffiner par NR
        deltaxi=huge(1._dp)
        do while(abs(deltaxi)>deltaximax)
            p1=1._dp
            p2=0.
            do j=1,n
                p3=p2
                p2=p1
                p1=((2*j-1._dp)*xi*p2-(j-1._dp)*p3)/j    ! relation de recurrence entre les Pj
            end do
            pp=n*(xi*p1-p2)/(xi**2-1._dp)            ! donne Pn' en fonction de Pn et Pn-1
            deltaxi=-p1/pp                          ! NR
            xi=xi+deltaxi
        end do
        x(i)=xi
        w(i)=1._dp/((1._dp-xi**2)*pp**2)         ! poids normalise a 1
        x(n+1-i)=-xi
        w(n+1-i)=w(i)
    end do
end subroutine gauleg      !


end module

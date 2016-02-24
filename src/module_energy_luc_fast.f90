module tableaux_dft_3d_fast
  INTEGER::mnmax,mnmax2,nbeta,nphi,nomeg
  REAL(8),DIMENSION(:),ALLOCATABLE::beta,cosbeta,sinbeta,wb
  REAL(8),DIMENSION(:),ALLOCATABLE::phi,cosphi,sinphi
  COMPLEX(8),DIMENSION(:),ALLOCATABLE::expiphi
  COMPLEX(8),DIMENSION(:,:),ALLOCATABLE::expikhiphi
  REAL(8),DIMENSION(:),ALLOCATABLE::omeg
  COMPLEX(8),DIMENSION(:),ALLOCATABLE::expiomeg
  COMPLEX(8),DIMENSION(:,:),ALLOCATABLE::expimuomeg
  REAL(8),DIMENSION(:,:,:,:),ALLOCATABLE::harsph
  COMPLEX(8),DIMENSION(:,:),ALLOCATABLE::harsph_q
  COMPLEX(8),DIMENSION(:,:,:),ALLOCATABLE::harsph1_q
  COMPLEX(8),DIMENSION(:,:,:),ALLOCATABLE::delta_rho_proj,delta_rho,delta_rho_proj1,auxi_1,auxi_2
  COMPLEX(8),DIMENSION(:,:,:),ALLOCATABLE::gamma1,gamma1_proj,gamma2,gamma2_proj,gamma3_proj,gamma4_proj,gamma4_proj1
  INTEGER, DIMENSION(:),ALLOCATABLE:: mm,nn,ll,mumu,nunu
  REAL(8), DIMENSION(:),ALLOCATABLE:: ck
  real(8), allocatable :: my_ck(:,:)
  COMPLEX(8), DIMENSION(:,:,:,:,:,:),ALLOCATABLE:: ck_omega_omega
  COMPLEX(8),DIMENSION(:,:,:,:,:),ALLOCATABLE:: tab5,tab6
  COMPLEX(8),DIMENSION(:,:,:,:),ALLOCATABLE:: tab4,tab7
  COMPLEX(8),DIMENSION(:,:,:),ALLOCATABLE:: tab3
end module tableaux_dft_3d_fast
!
!
module module_energy_luc_fast
  real(8) :: fac(0:50)
  complex(8), parameter :: zeroc=cmplx(0.d0,0.d0)
  private
  public :: energy_luc_fast
contains



subroutine energy_luc_fast (ff,df)
  use iso_c_binding
  use module_grid, only: grid
  use module_solvent, only: solvent
  use module_wigner_d, only: wigner_big_D
  implicit none
  include 'fftw3.f03'
  integer, parameter :: dp=c_double
  real(dp), intent(out) :: ff
  real(dp), intent(inout) :: df(:,:,:,:,:)
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
  complex(dp) :: a, b
  integer :: imqx, imqy, imqz, ntheta, nphi, npsi
  real(dp) :: fm(0:grid%mmax)

  ff=0
  no = grid%no
  nx = grid%nx
  ny = grid%ny
  nz = grid%nz
  np = grid%np
  rho0 = solvent(1)%rho0
  mmax=grid%mmax
  mmax2=mmax/2
  ntheta=grid%ntheta
  nphi=grid%nphi
  npsi=grid%npsi
  fm(0:mmax) = [( sqrt(real(2*m+1,dp))  ,m=0,mmax  )]
  call cpu_time(time(1))


allocate ( in(nx,ny,nz), source=zeroc)
allocate ( out(nx,ny,nz), source=zeroc)
call dfftw_plan_dft_3d (plan_fft_c2c_3d_signe_plus,  nx, ny, nz, in, out, FFTW_BACKWARD, FFTW_ESTIMATE) ! le tag FFTW_BACKWARD indique un signe + dans l'exponentiel
call dfftw_plan_dft_3d (plan_fft_c2c_3d_signe_minus, nx, ny, nz, in, out, FFTW_FORWARD, FFTW_ESTIMATE) ! le tag FFTW_FORWARD indique un signe - dans l'exponentiel

!
!
! 1. FOURIER TRANSFORM DELTA RHO
!
!
! call random_number( solvent(1)%xi )
! solvent(1)%xi = solvent(1)%xi * 10

allocate ( delta_rho_k_angle(no,nx,ny,nz) , source=zeroc)
allocate ( gamma_k_angle(no,nx,ny,nz) , source=zeroc)

do io=1,no
  in = cmplx(solvent(1)%xi(io,:,:,:)**2*rho0-rho0,0,dp)  ! Δρ(Ω,x,y,z)
  call dfftw_execute_dft( plan_fft_c2c_3d_signe_plus, in, out )
  delta_rho_k_angle(io,:,:,:) = out  ! Δρ(Ω,qx,qy,qz)
end do


!
! on s'assure :
! - la symétrie hermitienne est gardée
! - surtout, on sait aller chercher q et -q dans les tableaux pour tous les q. Il n'y a pas d'exception, que nx, ny, nz soient pairs ou impairs.
! Autrement dit, ce n'est pas vraiment la symétrie hermitienne qu'on teste, mais la fonction grid%ix_mq(ixq) qui a l'indice iqx trouve l'indice ix_mq correspondant à -q(ix_q)
!
goto 177
do concurrent (iqx=1:nx, iqy=1:ny, iqz=1:nz, io=1:no)
  imqx = grid%ix_mq(iqx)
  imqy = grid%iy_mq(iqy)
  imqz = grid%iz_mq(iqz)
  a = delta_rho_k_angle(io,iqx,iqy,iqz)
  b = conjg(delta_rho_k_angle(io,imqx,imqy,imqz))
  if (abs(a-b)>1.D-10)  error stop "raté pour la symétrie hermitienne de delta_rho_k_angle"
end do
177 continue







allocate (delta_rho_k_proj(0:mmax,-mmax:mmax,-mmax2:mmax2) , source=zeroc)
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
block
  use module_wigner_d, only: wigner_small_d
  integer :: itheta, iphi, ipsi, ntheta, nphi, npsi, io
  complex(dp), allocatable :: in(:,:), out(:,:), delta_rho_theta_mup_mu(:,:,:), delta_rho_p(:,:,:)
  complex(dp), parameter :: zeroc=cmplx(0,0)
  type(c_ptr) :: planc2c2dplus
  mmax   = grid%mmax
  mmax2  = mmax/2
  ntheta = grid%ntheta
  nphi   = grid%nphi
  npsi   = grid%npsi
  allocate( in(nphi,npsi), out(nphi,npsi), source=zeroc)
  allocate( delta_rho_theta_mup_mu(ntheta,-mmax:mmax,-mmax2:mmax2), source=zeroc)
  call dfftw_plan_dft_2d (planc2c2dplus, nphi, npsi, in, out, FFTW_BACKWARD, FFTW_ESTIMATE) ! le tag FFTW_BACKWARD indique un signe + dans l'exponentiel
  do itheta=1,ntheta
    in = zeroc
    do concurrent( iphi=1:nphi, ipsi=1:npsi )
      in(iphi,ipsi) = delta_rho_k_angle(grid%io(itheta,iphi,ipsi), iqx,iqy,iqz) /real(nphi*npsi)
    end do
    call dfftw_execute (planc2c2dplus, in, out)

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
        delta_rho_theta_mup_mu(itheta,mup,mu2) = out(iphi,ipsi)
      end do
    end do
  end do ! sur theta

  allocate ( delta_rho_p(0:mmax,-mmax:mmax,-mmax2:mmax2), source=zeroc)

  do m=0,mmax
    m2 = m/2
    do mup=-m,m
      do mu2=-m2,m2
        mu=2*mu2
        do itheta=1,ntheta
delta_rho_p(m,mup,mu2) = delta_rho_p(m,mup,mu2) + delta_rho_theta_mup_mu(itheta,mup,mu2)&
         * wigner_small_d(m,mup,mu,grid%thetaofntheta(itheta)) * grid%wthetaofntheta(itheta) &
         /sum(grid%wthetaofntheta) * fm(m)
        end do
      end do
    end do
  end do

! ANCIENNE METHODE VERIFIEE OK
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

! COMPARAISON
do m=0,mmax
  do mup=-m,m
    do mu2=-m/2,m/2
      mu=2*mu2
      a = delta_rho_k_proj(m,mup,mu2)
      b = delta_rho_p(m,mup,mu2)
      ! if(abs(a-b)>1.D-10) print*, m,mup,mu,a, b, "le rapide est ",real(b)/real(a), "fois plus grand que b"
      if(abs(a-b)>1.D-10 )  print*, iqx,iqy,iqz,m,mup,mu,"LucREF=",real(a),"LucFFT=",real(b)
    end do
  end do
end do
deallocate(in,out,delta_rho_p,delta_rho_theta_mup_mu)
end block
cycle
stop "block projection fini"



! delta_rho_k_proj = cmplx(0,0)
! do m=0,mmax
!   do mup=-m,m
!     do mu2=-m/2,m/2
!       mu=2*mu2
!       delta_rho_k_proj(m,mup,mu2) = sqrt(real(2*m+1,dp)) /sum(grid%w) *sum( delta_rho_k_angle(:,iqx,iqy,iqz) *grid%w *  &
!         [( conjg(wigner_big_D(m,mup,mu,grid%theta(io),grid%phi(io),grid%psi(io)))  ,io=1,no )]   )
!     end do
!   end do
! end do

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
  call dfftw_execute_dft( plan_fft_c2c_3d_signe_minus, in, out)
  gamma_angle(io,:,:,:) = out/real(nx*ny*nz,dp)
end do


call cpu_time (time(3))
deallocate(in)
deallocate(delta_rho_k_angle)

print*, "energy_luc took ",time(3)-time(1),"sec"
print*, "ff luc fast =", ff

end subroutine energy_luc_fast






  subroutine luc_oz (qvec, my_delta_rho_proj, my_gamma_proj)
    USE lecture                                    ! pour lire des valeurs � la Luc
    USE tableaux_dft_3d_fast                            ! d�finit le type des tableaux de valeurs � partager
    use module_rotation, only: thetaofq, phiofq
    use module_grid,only: grid
    implicit real(8) (a-h,o-z)
    real(8),    intent(in) :: qvec(3)
    complex(8), intent(in) :: my_delta_rho_proj(0:grid%mmax,-grid%mmax:grid%mmax,-grid%mmax/2:grid%mmax/2)
    complex(8), intent(out) :: my_gamma_proj(0:grid%mmax,-grid%mmax:grid%mmax,-grid%mmax/2:grid%mmax/2)
    CHARACTER(50)::nomfic,texte
    CHARACTER(10) texte1
    COMPLEX(8)::xi_cmplx,auxi_cmplx,ck_cmplx,coeff_cmplx,coeff1_cmplx,coeff2_cmplx
    integer :: mmax
    logical :: qz_was_negatif
    !
    !interface
    !subroutine proj_angl(f_proj,f_ang)
    !!implicit real(8) (a-h,o-z)
    !COMPLEX(8), DIMENSION(:,:,:):: f_proj,f_ang
    !end subroutine
    !END interface
    !
    !interface
    !subroutine angl_proj(f_ang,f_proj)
    !!implicit real(8) (a-h,o-z)
    !COMPLEX(8), DIMENSION(:,:,:):: f_proj,f_ang
    !end subroutine
    !END interface
    !
    fac(0)=1.
    do k=1,50
      fac(k)=fac(k-1)*REAL(k)
    end do
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
    pi=4.d0*ATAN(1.d0)
    xi_cmplx=(0.,1.d0)
    !                       theta (ou beta)
    if(.not.allocated(beta)) then
      ALLOCATE(beta(nbeta),cosbeta(nbeta),sinbeta(nbeta),wb(nbeta))
      call gauleg(cosbeta,wb,nbeta)               ! cosbeta()=racines de Pnbeta(x), wb() le poids
      beta=ACOS(cosbeta)
      sinbeta=SIN(beta)
    end if
    !                       phi
    if(.not.allocated(phi)) then
      ALLOCATE(phi(nphi),cosphi(nphi),sinphi(nphi),expiphi(nphi))
      phi=2.d0*pi* (/(i-0.5,i=1,nphi)/) /REAL(nphi)               ! sur 0-2pi
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

    theta_q = thetaofq(qx,qy,qz)
    phi_q   = phiofq(qx,qy,qz)

    ! PRINT*, 'Coordonnees x,y,z du vecteur q  (en A-1) --->',qx,qy,qz,'ce qui donne q --->',qq


    !
    !
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
    !     phi_q=acos(-1.d0)
    !   else
    !     phi_q=ACOS(c_phi_q)
    !     IF(qy<0.) phi_q=-phi_q
    !   end if
    !   s_phi_q=SIN(phi_q)
    ! end if


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
    !     phi_q=acos(-1.d0)
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
            !   fmn=SQRT((2.d0*m+1.d0)*(2.d0*n+1.d0))
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
            !                       coeff1_cmplx=SQRT((2.d0*m+1.d0)*(2.d0*n+1.d0))*tab6(m,n,l,mu2,nu2)
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
            !   fmn=SQRT((2.d0*m+1.d0)*(2.d0*n+1.d0))
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
            !       coeff_cmplx=coeff_cmplx+(2.d0*l+1.d0)*fmn*s3j*CONJG(harsph_q(l,lambda1))*coeff1_cmplx
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
            if(.not.allocated(tab3)) ALLOCATE (tab3(-mnmax:mnmax,-mnmax2:mnmax2,-mnmax2:mnmax2),       &
                                               tab4(0:mnmax,-mnmax:mnmax,-mnmax2:mnmax2,-mnmax2:mnmax2),       &
                                               tab5(0:mnmax,0:mnmax,-mnmax:mnmax,-mnmax2:mnmax2,-mnmax2:mnmax2))
            ! PRINT*, 'ckhi'
            !
            tab5=0.
            do khi=0,mnmax                ! que khi>=0 pour l'instant
              do ialp=1,ialpmax
                m=mm(ialp); n=nn(ialp); l=ll(ialp); mu=mumu(ialp); nu=nunu(ialp)
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



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! !         test a partir de tab6 , reussi
            ! GOTO 217                   ! je shunte les 2 tests a venir puisque reussis
            ! PRINT*, 'test a partir de tableaux a 5 entrees cmnmunu'
            ! do m=0,mnmax
            !   do n=0,mnmax
            !     do khi=-mnmax,mnmax
            !       do mu2=-mnmax2,mnmax2
            !         do nu2=-mnmax2,mnmax2
            !           coeff_cmplx=0.
            !           do l=0,2*mnmax
            !             coeff_cmplx=coeff_cmplx+symbol_3j(m,n,l,khi,-khi,0)*tab6(m,n,l,mu2,nu2)
            !           end do
            !       IF(ABS(tab5(m,n,khi,mu2,nu2)-coeff_cmplx)>1.d-10) PRINT*, m,n,khi,2*mu2,2*nu2,coeff_cmplx,tab5(m,n,khi,mu2,nu2)
            !         end do
            !       end do
            !     end do
            !   end do
            ! end do
            ! !     fin test
            ! !         test l-khi-l  , reussi
            ! PRINT*, 'test c l khi l'
            ! do ialp=1,ialpmax
            !   m=mm(ialp); n=nn(ialp); l=ll(ialp); mu=mumu(ialp); nu=nunu(ialp)      ! mnlmunu
            !   coeff_cmplx=0.
            !   do khi=-mnmax,mnmax
            !     coeff_cmplx=coeff_cmplx+(2.*l+1.)*symbol_3j(m,n,l,khi,-khi,0)*tab5(m,n,khi,mu/2,nu/2)
            !   end do
            !   coeff1_cmplx=ck(ialp)
            !   IF((-1)**l==-1) coeff1_cmplx=xi_cmplx*coeff1_cmplx
            !   IF(ABS(coeff_cmplx-coeff1_cmplx)>1.d-10) PRINT*, m,n,l,mu,nu,coeff1_cmplx,coeff_cmplx
            !   IF(m/=0.or.n/=0) then
            !     m=mm(ialp); n=nn(ialp); l=ll(ialp); mu=-mumu(ialp); nu=-nunu(ialp)    ! idem mnl-mu-nu
            !     coeff_cmplx=0.
            !     do khi=-mnmax,mnmax
            !       coeff_cmplx=coeff_cmplx+(2.*l+1.)*symbol_3j(m,n,l,khi,-khi,0)*tab5(m,n,khi,mu/2,nu/2)
            !     end do
            !     coeff1_cmplx=(-1)**(m+n+l)*ck(ialp)
            !     IF((-1)**l==-1) coeff1_cmplx=xi_cmplx*coeff1_cmplx
            !     IF(ABS(coeff_cmplx-coeff1_cmplx)>1.d-10) PRINT*, m,n,l,mu,nu,coeff1_cmplx,coeff_cmplx
            !   endif
            !   IF(m/=n.or.ABS(mu)/=ABS(nu)) then
            !     m=nn(ialp); n=mm(ialp); l=ll(ialp); mu=nunu(ialp); nu=mumu(ialp)      ! idem nmlnumu
            !     coeff_cmplx=0.
            !     do khi=-mnmax,mnmax
            !       coeff_cmplx=coeff_cmplx+(2.*l+1.)*symbol_3j(m,n,l,khi,-khi,0)*tab5(m,n,khi,mu/2,nu/2)
            !     end do
            !     coeff1_cmplx=(-1)**(m+n)*ck(ialp)
            !     IF((-1)**l==-1) coeff1_cmplx=xi_cmplx*coeff1_cmplx
            !     IF(ABS(coeff_cmplx-coeff1_cmplx)>1.d-10) PRINT*, m,n,l,mu,nu,coeff1_cmplx,coeff_cmplx
            !     IF(m/=0.or.n/=0) then
            !       m=nn(ialp); n=mm(ialp); l=ll(ialp); mu=-nunu(ialp); nu=-mumu(ialp)    ! idem nml-nu-mu
            !       coeff_cmplx=0.
            !       do khi=-mnmax,mnmax
            !         coeff_cmplx=coeff_cmplx+(2.*l+1.)*symbol_3j(m,n,l,khi,-khi,0)*tab5(m,n,khi,mu/2,nu/2)
            !       end do
            !       coeff1_cmplx=(-1)**l*ck(ialp)
            !       IF((-1)**l==-1) coeff1_cmplx=xi_cmplx*coeff1_cmplx
            !       IF(ABS(coeff_cmplx-coeff1_cmplx)>1.d-10) PRINT*, m,n,l,mu,nu,coeff1_cmplx,coeff_cmplx
            !     endif
            !   endif
            ! end do
            ! !      fin test
            ! 217 continue
            ! !GOTO 300
            ! !
            ! ALLOCATE(gamma2(nbeta,nphi,nomeg),gamma2_proj(0:mnmax,-mnmax:mnmax,-mnmax2:mnmax2))
            ! PRINT*, 'c(o,o'') et convolution et test'
            ! gamma2=0.
            ! !          pour chaque omega du repere fixe, transformer en Omega du repere local
            ! do i=1,nbeta
            !   do j=1,nphi
            !     cphiq=cosphi(j)*c_phi_q+sinphi(j)*s_phi_q       ! cos(phi-phiq)
            !     sphiq=sinphi(j)*c_phi_q-cosphi(j)*s_phi_q       ! sin(phi-phiq)
            !     cbet=c_theta_q*cosbeta(i)+s_theta_q*sinbeta(i)*cphiq
            !     bet=ACOS(cbet)
            !     sbet=SIN(bet)
            !     cphi=(c_theta_q*sinbeta(i)*cphiq-s_theta_q*cosbeta(i))/sbet
            !     ph=ACOS(cphi)
            !     IF(sphiq<0.) ph=-ph
            !     cpsiq=(-s_theta_q*cosbeta(i)*cphiq+c_theta_q*sinbeta(i))/sbet
            !     psiq=ACOS(cpsiq)
            !     IF(sphiq>0.) psiq=-psiq
            !     do k=1,nomeg
            !       ps=psiq+omeg(k)
            !       !      test retour o vers O
            !       cc=c_theta_q*cbet-s_theta_q*sbet*cphi
            !       bb=ACOS(cc)
            !       ss=SIN(bb)
            !       ccphiq=(c_theta_q*sbet*cphi+s_theta_q*cbet)/ss
            !       pphiq=ACOS(ccphiq)
            !       IF(ph<0.) pphiq=-pphiq
            !       pphi=pphiq+phi_q
            !       ccpsiq=(s_theta_q*cbet*cphi+c_theta_q*sbet)/ss
            !       ppsiq=ACOS(ccpsiq)
            !       IF(ph<0.) ppsiq=-ppsiq
            !       oo=ppsiq+ps
            !       !IF(ABS(bb-beta(i))+ABS(
            !       !                               Idem pour omega'
            !       do i1=1,nbeta
            !         do j1=1,nphi
            !           cphiq1=cosphi(j1)*c_phi_q+sinphi(j1)*s_phi_q       ! cos(phi-phiq)
            !           sphiq1=sinphi(j1)*c_phi_q-cosphi(j1)*s_phi_q       ! sin(phi-phiq)
            !           cbet1=c_theta_q*cosbeta(i1)+s_theta_q*sinbeta(i1)*cphiq1
            !           bet1=ACOS(cbet1)
            !           sbet1=SIN(bet1)
            !           cphi1=(c_theta_q*sinbeta(i1)*cphiq1-s_theta_q*cosbeta(i1))/sbet1
            !           ph1=ACOS(cphi1)
            !           IF(sphiq1<0.) ph1=-ph1
            !           cpsiq1=(-s_theta_q*cosbeta(i1)*cphiq1+c_theta_q*sinbeta(i1))/sbet1
            !           psiq1=ACOS(cpsiq1)
            !           IF(sphiq1>0.) psiq1=-psiq1
            !           do k1=1,nomeg
            !             ps1=psiq1+omeg(k1)
            !             !                               calculer alors c(omega,omega') dans le repere local
            !             !       transformee n-->angle beta2
            !             tab4=0.
            !             do n=0,mnmax
            !               n2=n/2
            !               coeff=SQRT(2.d0*n+1.d0)
            !               do khi=-n,n                    ! se contenter de khi>=0
            !                 do nu2=-n2,n2
            !                   nu=2*nu2
            !                   m1=ABS(khi)
            !                   tab4(m1:mnmax,khi,-mnmax2:mnmax2,nu2)=tab4(m1:mnmax,khi,-mnmax2:mnmax2,nu2)+       &
            !                   coeff*harm_sph(n,-khi,nu,bet1)*tab5(m1:mnmax,n,khi,-mnmax2:mnmax2,nu2)
            !                 end do
            !               end do
            !             end do
            !             !       transformee m-->angle beta1
            !             tab3=0.
            !             do m=0,mnmax
            !               m2=m/2
            !               coeff=SQRT(2.d0*m+1.d0)
            !               do khi=-m,m                    ! se contenter de khi>=0
            !                 do mu2=-m2,m2
            !                   mu=2*mu2
            !                   tab3(khi,mu2,-mnmax2:mnmax2)=tab3(khi,mu2,-mnmax2:mnmax2)+  &
            !                   coeff*harm_sph(m,khi,mu,bet)*tab4(m,khi,mu2,-mnmax2:mnmax2)
            !                 end do
            !               end do
            !             end do
            !             !       enfin, somme sur khi,mu,nu  qui donne ck
            !             coeff_cmplx=0.
            !             do khi=-mnmax,mnmax                       ! tjs khi>=0
            !               coeff=1.d0!; IF(khi>0) coeff=2.d0
            !               do nu2=-mnmax2,mnmax2
            !                 nu=2*nu2
            !                 do mu2=-mnmax2,mnmax2
            !                   mu=2*mu2
            !                   coeff_cmplx=coeff_cmplx+coeff*tab3(khi,mu2,nu2)*EXP(-xi_cmplx*(khi*(ph-ph1)+mu*ps+nu*ps1))   ! cos pour avoir khi et -khi! luc91p91
            !                 end do
            !               end do
            !             end do
            !             IF(ABS(coeff_cmplx-ck_omega_omega(i,j,k,i1,j1,k1))>1.d-10) PRINT*, i,j,k,i1,j1,k1,	&
            !             ck_omega_omega(i,j,k,i1,j1,k1),coeff_cmplx
            !             !
            !             ! 		finalement, convoluer
            !             gamma2(i,j,k)=gamma2(i,j,k)+wb(i1)*coeff_cmplx*delta_rho(i1,j1,k1)
            !           end do
            !         end do
            !       end do
            !     end do
            !   end do
            ! end do
            ! gamma2=gamma2/(nphi*nomeg)
            ! !      et je projette
            ! PRINT*, 'gam_proj'
            ! call angl_proj(gamma2,gamma2_proj)
            !
            !
            !               Methode 3: je n'utilise que les projections, mais je combine dans le repere fixe
            !
            !
            300 continue
            ! PRINT*, '********************************************************************************'
            ! PRINT*, 'Methode 3: je combine directement les projections dans le repere fixe'
            IF(.NOT.ALLOCATED(TAB7)) THEN
              ALLOCATE (tab7(0:mnmax,0:mnmax,-mnmax2:mnmax2,-mnmax2:mnmax2))
              ALLOCATE(gamma3_proj(0:mnmax,-mnmax:mnmax,-mnmax2:mnmax2))
            END IF
            gamma3_proj=0.
            !            pour chaque mu',nu'
            do mu1=-mnmax,mnmax
              do nu1=-mnmax,mnmax
                !
                lambda1=-mu1-nu1
                tab4=0.
                do ialp=1,ialpmax
                  m=mm(ialp); n=nn(ialp); l=ll(ialp); mu=mumu(ialp); nu=nunu(ialp)
                  mu2=mu/2; nu2=nu/2
                  coeff_cmplx=ck(ialp)*harsph_q(l,lambda1)
                  IF((-1)**l==-1) coeff_cmplx=xi_cmplx*coeff_cmplx             ! imaginaire pur si l impair
                  s3j=symbol_3j(m,n,l,mu1,nu1,lambda1)
                  tab4(m,n,mu2,nu2)=tab4(m,n,mu2,nu2)+coeff_cmplx*s3j                     ! tab5 est donc complexe
                  IF(mu/=0.or.nu/=0) tab4(m,n,-mu2,-nu2)=tab4(m,n,-mu2,-nu2)+(-1)**(m+n+l)*coeff_cmplx*s3j
                  IF(m/=n.or.ABS(mu)/=ABS(nu)) then
                    s3j=symbol_3j(n,m,l,mu1,nu1,lambda1)
                    tab4(n,m,nu2,mu2)=tab4(n,m,nu2,mu2)+(-1)**(m+n)*coeff_cmplx*s3j
                    IF(mu/=0.or.nu/=0) tab4(n,m,-nu2,-mu2)=tab4(n,m,-nu2,-mu2)+(-1)**(l)*coeff_cmplx*s3j
                  endif
                end do                           ! fin ialp
                !                  test en partant plutot de tab6  , reussi donc shunte
                ! GOTO 317
                ! tab7=0.
                ! do mu2=-mnmax2,mnmax2
                !   mu=2*mu2
                !   do m=MAX(ABS(mu1),ABS(mu)),mnmax
                !     do nu2=-mnmax2,mnmax2
                !       nu=2*nu2
                !       do n=MAX(ABS(nu1),ABS(nu)),mnmax
                !         do l=ABS(m-n),m+n
                !       tab7(m,n,mu2,nu2)=tab7(m,n,mu2,nu2)+tab6(m,n,l,mu2,nu2)*symbol_3j(m,n,l,mu1,nu1,lambda1)*harsph_q(l,lambda1)
                !         end do
                !   IF(ABS(tab7(m,n,mu2,nu2)-tab4(m,n,mu2,nu2))>1.d-10) PRINT*, mu1,nu1,m,n,mu,nu,tab7(m,n,mu2,nu2),tab4(m,n,mu2,nu2)
                !       end do
                !     end do
                !   end do
                ! end do
                ! !            fin test
                ! 317 continue
                !
                do mu2=-mnmax2,mnmax2
                  mu=2*mu2
                  do m=MAX(ABS(mu1),ABS(mu)),mnmax
                    do nu2=-mnmax2,mnmax2
                      nu=2*nu2
                      do n=MAX(ABS(nu1),ABS(nu)),mnmax
                        gamma3_proj(m,mu1,mu2)=gamma3_proj(m,mu1,mu2)+tab4(m,n,mu2,nu2)*(-1)**nu1*delta_rho_proj(n,-nu1,-nu2)
                      end do
                    end do
                  end do
                end do
                !
              end do                     ! fin nu'
            end do                      ! fin mu'
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
    !GOTO 417
    ! PRINT*, 'test retour a deltarho_proj'
    gamma4_proj=0.
    do m=0,mnmax
      m2=m/2
      do mu1=-m,m
        do mu2=-m2,m2
          gamma4_proj(m,mu1,mu2)=SUM(delta_rho_proj1(m,-m:m,mu2)*CONJG(harsph1_q(m,mu1,-m:m)))
        end do
      end do
    end do
    IF(MAXVAL(ABS(gamma4_proj-delta_rho_proj))>1.d-10) PRINT*, 'AR pas bon!'
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



    !
    !
    !                 OUF!
    !
    !
    !       resultats
    !
    ! PRINT*, '********************************************************************************'
    ! PRINT*, 'RESULTATS comparatifs:'
    do m=0,mnmax
      m2=m/2
      do mu1=-m,m
        do mu2=-m2,m2
          mu=2*mu2
          ! PRINT*, m,mu1,mu,gamma1_proj(m,mu1,mu2),gamma2_proj(m,mu1,mu2),gamma3_proj(m,mu1,mu2),gamma4_proj(m,mu1,mu2)
          if(gamma4_proj(m,mu1,mu2)/=gamma4_proj(m,mu1,mu2)) then
            print*, qx,qy,qz
            error stop "je suis triste"
          end if
        end do
      end do
    end do
! stop "RESULTATS comparatifs"


    ! PRINT*, MAXVAL(ABS(gamma1_proj-gamma4_proj)),MAXVAL(ABS(gamma2_proj-gamma4_proj)),MAXVAL(ABS(gamma3_proj-gamma4_proj))
    !
    !

    ! deallocate(beta,cosbeta,sinbeta,wb)
    ! deallocate(phi,cosphi,sinphi)
    ! deallocate(expiphi,expikhiphi)
    ! deallocate(omeg,expiomeg,expimuomeg)
    ! deallocate(harsph,harsph_q,harsph1_q)
    ! deallocate(delta_rho_proj,delta_rho,delta_rho_proj1,auxi_1,auxi_2)
    ! deallocate(gamma1,gamma1_proj,gamma2,gamma2_proj,gamma3_proj,gamma4_proj,gamma4_proj1)
    ! deallocate(mm,nn,ll,mumu,nunu)
    ! deallocate(ck,ck_omega_omega)
    ! DEALLOCATE(tab3,tab5,tab6,tab4,tab7)
    ! GOTO 9999
    !
  end subroutine


  !
  subroutine proj_angl(f_proj,f_ang)
    USE tableaux_dft_3d_fast, ONLY: mnmax,mnmax2,nbeta,nphi,nomeg,auxi_1,auxi_2,harsph,expikhiphi,expimuomeg
    implicit real(8) (a-h,o-z)
    !COMPLEX(8), DIMENSION(:,:,:):: f_proj,f_ang
    COMPLEX(8), DIMENSION(0:mnmax,-mnmax:mnmax,-mnmax2:mnmax2), intent(in):: f_proj
    COMPLEX(8), DIMENSION(nbeta,nphi,nomeg), intent(out):: f_ang
    !
    !             passe des projections a une fonction angulaire
    !
    !             d'abord, je transforme m en theta
    auxi_1=0.
    do mu2=-mnmax2,mnmax2
      mu=2*mu2
      do khi=-mnmax,mnmax
        do m=MAX(ABS(khi),ABS(mu)),mnmax
          auxi_1(1:nbeta,khi,mu2)=auxi_1(1:nbeta,khi,mu2)+SQRT(2.d0*m+1.d0)*f_proj(m,khi,mu2)*harsph(m,khi,mu2,1:nbeta)
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
  end subroutine
  !
  !
  subroutine angl_proj(f_ang,f_proj)
    USE tableaux_dft_3d_fast, ONLY: mnmax,mnmax2,nbeta,nphi,nomeg,auxi_1,auxi_2,wb,harsph,expikhiphi,expimuomeg
    implicit real(8) (a-h,o-z)
    !COMPLEX(8), DIMENSION(:,:,:):: f_proj,f_ang
    COMPLEX(8), DIMENSION(nbeta,nphi,nomeg), intent(in):: f_ang
    COMPLEX(8), DIMENSION(0:mnmax,-mnmax:mnmax,-mnmax2:mnmax2), intent(out):: f_proj
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
          f_proj(m,khi,mu2)=f_proj(m,khi,mu2)+SQRT(2.d0*m+1.d0)*SUM(wb(:)*auxi_1(:,khi,mu2)*harsph(m,khi,mu2,:))
        end do
      end do
    end do
    !
  end subroutine
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

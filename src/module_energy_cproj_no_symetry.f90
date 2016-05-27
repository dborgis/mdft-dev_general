! Je vous demande de ne pas vous préoccuper des allocations etc pour l'instant.
! L'objectif ici est d'avoir une version simple de la routine, pouvant servir de base de travail.
! Notez que je ne tiens pas compte ici de la symétrie hermicienne (du fait que rho est réel)
! Je ne tiens pas compte de la symétrie sur psi de la molécule d'eau (c'est à dire que j'ai npsi=2mmax+1).

module module_energy_cproj_no_symetry

  use iso_c_binding, only: C_PTR, C_INT, C_INT32_T, C_INTPTR_T, C_DOUBLE_COMPLEX, C_DOUBLE, C_FUNPTR, C_SIZE_T, C_FLOAT, &
                           C_FLOAT_COMPLEX, C_CHAR
  use precision_kinds, only: dp

  implicit none
  private

  !
  ! FFTW3 header - modern (fortran 2003) version. Expects iso_c_binding
  !
  include 'fftw3.f03'

  !
  ! Parameters
  !
  real(dp), parameter :: zero=0._dp
  complex(dp), parameter :: zeroc=complex(0._dp, 0._dp)


  !
  ! Direct correlation functions from Luc
  !
  type :: cq_type
    logical :: isok = .false.
    real(dp) :: q
    integer :: na ! number of projections for the direct correlation function
    integer :: nq=1024 ! number of q points in c(q)
    integer, allocatable :: a(:,:,:,:,:) ! index of m, n, mu, nu, khi => (0:mmax,0:mmax,-mmax:mmax,-mmax:mmax,-mmax:mmax) ! m n mu nu khi
    real(dp) :: dq
    real(dp), allocatable :: normq(:)
    integer, allocatable :: m(:), n(:), mu(:), nu(:), khi(:)
  end type
  type (cq_type) :: cq
  complex(dp), allocatable :: ck(:,:)


  type :: p3_type
    integer, allocatable :: p(:,:,:) ! index of the projection corresponding to m, mup, mu
    integer, allocatable :: m(:) ! m for projection 1 to np
    integer, allocatable :: mup(:) ! mup for projection 1 to np. mup corresponds to phi
    integer, allocatable :: mu(:) ! mu for projection 1 to np. mu corresponds to psi
  end type p3_type
  type (p3_type) :: p3

  public :: energy_cproj_no_symetry

contains

  subroutine energy_cproj_no_symetry (ff,df)

    use iso_c_binding, only: c_ptr
    use precision_kinds, only: dp
    use module_rotation, only: rotation_matrix_between_complex_spherical_harmonics_lu
    use module_wigner_d, only: wigner_big_d
    use module_grid, only: grid
    use module_thermo, only: thermo
    use module_solvent, only: solvent

    implicit none

    real(dp), intent(out) :: ff
    real(dp), intent(inout) :: df(:,:,:,:,:) ! grad(f) of omega, x, y, z, solventId
    real(dp) :: dv, kT
    integer :: ix, iy, iz, ix_q, iy_q, iz_q, p
    integer :: nx, ny, nz, np, no, ns, mmax, na, nq, mrso
    integer :: m, n, mu, nu, khi, mup, ip, iq, io
    real(dp) :: q(3), lx, ly, lz, rho0
    real :: time(20)
    real(dp) :: vexc(grid%no), rho, xi

    !
    ! f_m tel que définit par Luc en dessous de l'équation (1.4)
    !
    real(dp), parameter :: fm(0:7) = [( sqrt(real(2*m+1,dp)) ,m=0,7 )]


    !
    ! dans tous les tableaux qui suivent, quasi tous sont inutiles (un ou deux suffiront), mais je les garde pour faire tout un tas
    ! de test dessus au fur et à mesure
    ! On utilise toujours la première colonne pour les mettre les indices des angles (1 à no) ou les indices des projections (1 à np)
    !
    complex(dp), allocatable :: deltarho_k(:,:,:,:), & ! j'appelle _k ce qui est dans l'espace de fourier
                                deltarho_k_p(:,:,:,:),& ! j'appelle _p ce qui est en projections
                                deltarho_prime_k_p(:,:,:,:), & ! j'appelle _prime ce qui est dans le repère lié à k
                                gamma_prime_k_p(:,:,:,:), & ! gamma_prime_k_p est donc les projections de gamma dans l'espace de fourier dans le repère lié à k
                                gamma_k_p(:,:,:,:), & ! projections de gamma dans le repère fixe dans l'espace de fourier
                                gamma_k(:,:,:,:) ! gamma en orientations  (== en angles) dans l'espace de fourier (et repère fixe)
    real(dp), allocatable ::    gamma(:,:,:,:) ! gamma dans l'espace réel, repère fixe, en orientations

    complex(dp) :: cfulltest(2002,200) ! DOIT ETRE VIRE DES QUE CA FONCTIONNE TODO
    integer :: ptest(0:5,0:5,-5:5,-5:5,-5:5)  ! DOIT ETRE VIRE DES QUE CA FONCTIONNE TODO m n mu nu khi
    integer :: iqtest, mu2
    complex(dp) :: a, b
    logical :: resultat

    call cpu_time(time(12))

    mmax=grid%mmax ! lu au début du programme
    mrso=grid%molrotsymorder ! =1

    ! on definit quelques raccourcis
    lx=grid%lx
    ly=grid%ly
    lz=grid%lz
    kT=thermo%kbT
    dv=grid%dv
    nx=grid%nx
    ny=grid%ny
    nz=grid%nz
    np=grid%np
    no=grid%no
    na=cq%na
    nq=cq%nq
    ns=solvent(1)%nspec
    rho0 = solvent(1)%rho0

    !
    ! Mes routines ne sont valables que pour un nombre pair de noeuds dans chaque direction
    !
    if (mod(nx,2)/=0 .or. mod(ny,2)/=0 .or. mod(nz,2)/=0) then
      print*, "mdft-dev n'est valables que pour un nombre pair de noeuds dans chaque direction"
      print*, "or, on a nx,ny,nz = ",nx,ny,nz
      error stop
    end if



    !
    ! Toutes les projections sont stockées dans un meme vecteur de taille np
    ! p3%p(m,mup,mu) transforme le triplet (m,mup,mu) a l'indice de la projection, ip.
    ! p3%m(ip) donne, inversement à p3%p(m,mup,mu), le m (de {m,mup,mu}) qui correspondant à l'indice ip dans le tableau des projections
    ! On a donc p3%p(p3%m(p),p3%mup(p),p3%mu(p)) == p
    !
    ! On ne fait ce travail que la première fois qu'on passe dans cette routine.
    ! Les fois suivantes, p3%p reste alloué, et donc le .not.allocated(p3%p) est skippé.
    ! On utilise cette astuce de nombreuses fois dans cette routine.
    !
    if (.not.allocated(p3%p)) then
      allocate ( p3%p(0:mmax, -mmax:mmax, -mmax:mmax) ,source=-huge(1)) ! Dans p3%p, on met mu2=2*mu, pas mu
      allocate ( p3%m(np) ,source=-huge(1)) ! mettre des -huge comme valeur initiale permet de plus tard verifier que s'il reste un -huge qqpart, il y a un problème. Si on avait mis 0, on ne peut pas vérifier.
      allocate ( p3%mup(np) ,source=-huge(1))
      allocate ( p3%mu(np) ,source=-huge(1)) ! c'est bien mu qu'on met dans p3%mu(), pas mu2
      ip=0
      do m=0,mmax
        do mup=-m,m
          !
          ! We chose to have p3%p and p3%mu to contain the true value of mu, not mu/molrotsymorder.
          ! Thus, we loop over mu with steps of molrotsymorder
          ! Thus, we have holes in p3%p, but all calls to p3%p and p3%mu should be done with this true mu
          ! For instance, if molrotsymorder=2 (for water or any C2V molecule), any call to p3%p(m,mup,mu) with mu=1 will return -999
          ! Also, the array p3%mu only contains even values (des valeurs paires)
          !
          do mu=-m,m ! NOTE DANS CETTE ROUTINE ON A PAS LA SYMETRIE SUR L'EAU, LE COMMENTAIRE CI DESSUS ETAIT SURTOUT VALIDE POUR LA ROUTINE EFFICACE AVEC LES SYMETRIES.
            ip=ip+1
            if (ip > np) error stop "p > np at line 166"
            p3%p(m,mup,mu) = ip ! on met mu2 dans p3%p, pas mu
            p3%m(ip) = m
            p3%mup(ip) = mup
            p3%mu(ip) = mu ! c'est le vrai mu, pas mu2
          end do
        end do
      end do

      if (ip /= np) error stop "ip /= np in energy_cproj"

      if ( any(abs(p3%m)>mmax) .or. any(abs(p3%mup)>mmax) .or. any(abs(p3%mu)>mmax) ) then
        print*, "certaines valeurs de m, mup or mu have incorrect values"
        error stop
      end if

    end if




    !
    ! On va être carrément lourds et transférer xi (pour rappel xi²=rho/rho0) dans un tableau de complexes dont la partie imag est nulle.
    ! on va donc complètement oublier qu'on a une symétrie hermicienne, pour faire des FFT de complexe à complexe
    ! On testera ainsi deux choses : 1/ qu'on a bien en fin de cette routine un Ɣ(Ω,x,y,z) dont la partie imaginaire est nulle.
    ! ici on passe d'abord xi dans un tableau de complexes en k. On va travailler sans tenir compte de la symétrie hermitienne
    !
    allocate( deltarho_k(no,nx,ny,nz), source=zeroc )

    !
    ! Transformée de Fourier 3D de Δρ(Ω,x,y,z)
    !
    call cpu_time(time(1))
    block
      type(c_ptr) :: plan_fft_c2c_3d_signe_plus
      complex(dp), allocatable :: in(:,:,:)
      integer :: io
      allocate ( in(nx,ny,nz ))
      call dfftw_plan_dft_3d (plan_fft_c2c_3d_signe_plus, nx, ny, nz, in, in, FFTW_BACKWARD, FFTW_ESTIMATE) ! le tag FFTW_BACKWARD indique un signe + dans l'exponentiel
      do io=1,no
        in = cmplx(solvent(1)%xi(io,:,:,:)**2*rho0-rho0,0,dp)  ! solvent%xi est réel
        call dfftw_execute_dft( plan_fft_c2c_3d_signe_plus, in, in )
        deltarho_k(io,:,:,:) = in  ! Δρ(Ω,qx,qy,qz)
      end do
      deallocate (in)
    end block
    call cpu_time(time(2))
    print*, "<<<<< FFT en",time(2)-time(1),"seconds"


    do io=1,no
      call verify_hermitian_symetry (deltarho_k(io,:,:,:), resultat) ! since deltarho is real
      if (resultat.eqv..false.) error stop "Δρ(Ω,q) /= Δρ*(Ω,-q)"
    end do













    !
    ! PASSAGE EN PROJECTIONS
    !
    ! for each vector q, Δρ(Ω) -> Δρ m mup mu
    !
    call cpu_time(time(1))
    allocate( deltarho_k_p(np,nx,ny,nz) ,source=zeroc )

    block
      complex(dp), allocatable :: temp_wigner_big_D(:,:)
      integer :: o
      real(dp), allocatable :: fm_of_p(:)
      allocate (fm_of_p(np)) ! les fm_of_p sont juste les f(m) dans un tableau de taille 1 à np, et non de 0 à mmax
      do p=1,np
        fm_of_p(p) = fm(p3%m(p)) ! c'est juste un f_m de chaque projection.
      end do

      allocate (temp_wigner_big_D(no,np), source=zeroc)
      ! En fait, je tabule brutalement les harmoniques sphériques généralisées.
      ! La fonction wigner_big_D^m_µ'µ(θ) := wigner_small_d^m_µ'µ(θ) * exp(-ii.µ'.φ) * exp(-ii.µ.ψ)   Tel que définit dans Luc 1.4  (c'est un peu biblique, ça, non ? :)
      ! Savez vous que dans la bible, Luc 1.4 (Luc chapitre 1 verset 4) correspond au texte suivant : "afin que tu reconnaisses la certitude des enseignements que tu as reçus"
      ! Notez que je n'utilise pas la méthode rapide de Choi ici : seulement la fonction de Luc, harm_sph, qui correspond au calcul de wigner_small_d
!$OMP PARALLEL DO
      do p=1,np
        do o=1,no
          temp_wigner_big_D(o,p) = wigner_big_D(p3%m(p),p3%mup(p),p3%mu(p),grid%theta(o),grid%phi(o),grid%psi(o))
        end do
      end do
!$OMP END PARALLEL DO


      !
      ! On verifie les symetries sur wigner_big_D
      !
      do o=1,no
        do p=1,np
          m   = p3%m(p)
          mup = p3%mup(p)
          mu  = p3%mu(p)
          a = (-1)**(mup+mu) * conjg(temp_wigner_big_D(o,p))
          b = temp_wigner_big_D(o,p3%p(m,-mup,-mu))
          if (abs(a-b)>1.D-10) error stop "dans wigner_big_D qui ne verifie pas symetrie 1.6"
        end do
      end do
      !
      ! On prend le conjugué de la R^m_mup_mu(\Omega) car c'est bien le conjugué qui doit être utilisé pour la projection
      !
      temp_wigner_big_D = conjg(temp_wigner_big_D)
      do iz_q=1,nz
        do iy_q=1,ny
          do ix_q=1,nx
            do p=1,np
               ! les grid%w(1:no) sont les poids des quadratures angulaires. Pour rappel, la somme des poids fait 8 pi² ~ 79
              deltarho_k_p(p,ix_q,iy_q,iz_q) = fm_of_p(p) * sum( deltarho_k(:,ix_q,iy_q,iz_q) * grid%w(:) * temp_wigner_big_D(:,p) )
            end do
          end do
        end do
      end do
      deallocate (temp_wigner_big_D, fm_of_p)

    end block
    deallocate (deltarho_k)
    call cpu_time(time(2))
    print*, "<<<<< projection en r terminé en", time(2)-time(1), "seconds"


    call verify_hermitian_symetry_of_projections(deltarho_k_p ,resultat)
    if (resultat .eqv. .false.) error stop "see l. 294"




    !
    ! Vérification de la symétrie dûe a la fonction de départ purement réelle
    ! cf equations 1.15 de Luc
    !
    block
      complex(dp) :: a, b
      integer :: ix_q, iy_q, iz_q, ix_mq, iy_mq, iz_mq

      do iz_q=1,nz
        do iy_q=1,ny
          do ix_q=1,nx
            do p=1,np

              ix_mq = grid%ix_mq(ix_q)
              iy_mq = grid%iy_mq(iy_q)
              iz_mq = grid%iz_mq(iz_q)
              if (ix_mq<=0 .or. iy_mq<=0 .or. iz_mq<=0) cycle

              a = deltarho_k_p (p3%p(m, mup, mu),ix_mq,iy_mq,iz_mq)
              a = conjg(a)*(-1)**(mup+mu)
              b = deltarho_k_p (p3%p(m,-mup,-mu),ix_q,iy_q,iz_q)

              if (abs(a-b)/=0) then
                error stop "La symétrie 1.15 n'est pas vérifiée pour deltarho_k_p"
              end if

            end do

          end do
        end do
      end do
    end block




print*, "BLOCK OZ OZ OZ OZ OZ OZ OZ OZ OZ OZ OZ OZ OZ OZ OZ OZ OZ OZ OZ OZ OZ OZ OZ OZ OZ OZ OZ OZ OZ OZ FINI"

  block
    ! use module_oz, only: OZ
    logical :: foundbug

    complex(dp) :: delta_rho_proj(0:mmax,-mmax:mmax,-mmax/2:mmax/2) ! je garde les noms de variables de Luc
    delta_rho_proj = cmplx(0,0)


    do m=0,mmax
      do mup=-m,m
        do mu2=-m/2,m/2
          mu=mu2*2
          delta_rho_proj(m,mup,mu) = deltarho_k_p(p3%p(m,mup,mu),ix_q,iy_q,iz_q)
        end do
      end do
    end do


    do iz_q=1,nz
      do iy_q=1,ny
        do ix_q=1,nx
          q = [grid%kx(ix_q), grid%ky(iy_q), grid%kz(iz_q)]
          ! call OZ ( deltarho_k_p(1:np,ix_q,iy_q,iz_q), q, mmax, gamma_k_p(1:np,ix_q,iy_q,iz_q), foundbug )
          if (foundbug) error stop "l.348"
        end do
      end do
    end do
  end block

print*, "BLOCK OZ OZ OZ OZ OZ OZ OZ OZ OZ OZ OZ OZ OZ OZ OZ OZ OZ OZ OZ OZ OZ OZ OZ OZ OZ OZ OZ OZ OZ OZ FINI"
stop

    !
    ! ROTATION VERS REPAIRE LIE A Q
    !
    ! Δρ^m_µ'µ  => Δρ'^m_χ_µ
    !
    call cpu_time(time(1))
    allocate( deltarho_prime_k_p(np,nx,ny,nz) ,source=zeroc )
    block
      use module_wigner_d, only: wigner_big_D
      use module_rotation, only: thetaofq, phiofq
      complex(dp), allocatable :: R(:,:,:) ! R^m_{mup,mu} (\hat{q})
      allocate ( R(0:mmax,-mmax:mmax,-mmax:mmax) ,source=zeroc)
      do iz_q=1,nz
        do iy_q=1,ny
          do ix_q=1,nx

            q = [grid%kx(ix_q), grid%ky(iy_q), grid%kz(iz_q)]
            R = rotation_matrix_between_complex_spherical_harmonics_lu ( q, mmax ) ! Ici j'utilise Choi. C'est vérifié que c'est équivalent à wigner_big_D^m_mup_mu(q) for psi_q=0.

            do m=0,mmax
              do khi=-m,m
                do mu=-m,m
                  deltarho_prime_k_p ( p3%p(m,khi,mu) ,ix_q, iy_q, iz_q) = &
                      sum( deltarho_k_p( p3%p(m,-m:m,mu), ix_q,iy_q,iz_q) * R(m,-m:m,khi) ) ! la somme sur µ' des Δρ^m_µ'µ(qx,qy,qz)

                end do
              end do
            end do

          end do
        end do
      end do
    end block

    deallocate (deltarho_k_p)
    call cpu_time(time(2))
    print*, "<<<<< rotation vers repere moleculaire en",time(2)-time(1),"seconds"

    !
    ! Vérification de la relation de symétrie 1.25, c'est à dire la meme que 1.15 en notation χ
    !
    block
      complex(dp) :: a, b
      integer :: ix_q, iy_q, iz_q, ix_mq, iy_mq, iz_mq

      do iz_q=1,nz
        do iy_q=1,ny
          do ix_q=1,nx

            do p=1,np

              m = p3%m(p)
              khi = p3%mup(p)
              mu = p3%mu(p)

              !
              ! on cherche le l'indice correspond à -q
              !
              if (ix_q == 1) then
                ix_mq = ix_q
              else if (ix_q == nx/2 +1) then
                cycle
              else
                ix_mq = nx - ix_q + 2
              end if

              if (iy_q == 1) then
                iy_mq = iy_q
              else if (iy_q == ny/2 +1) then
                cycle
              else
                iy_mq = ny - iy_q + 2
              end if

              if (iz_q == 1) then
                iz_mq = iz_q
              else if (iz_q == nz/2 +1) then
                cycle
              else
                iz_mq = nz - iz_q + 2
              end if


              ! ATTENTION
              ! ICI LA NOTION DE Q ET -Q DEVRAIENT CHANGER PHI EN PHI+PI DANS LE CAS DE Q SELON Z
              a = deltarho_prime_k_p (p3%p(m,khi,mu),ix_mq,iy_mq,iz_mq)
              a = conjg(a)*(-1)**(m+mu+khi)
              b = deltarho_prime_k_p (p3%p(m,khi,-mu),ix_q,iy_q,iz_q)

              if (abs(a-b) >1.E-10 ) then
                ! if (ix_q==1 .and. iy_q==1 ) cycle ! VOIR ATTENTION JUSTE AU DESSUS !!!!! TODO ASK LUC DANIEL LU
                print*,
                print*, "ix_q, iy_q, iz_q",ix_q, iy_q, iz_q
                print*, "ix_mq, iy_mq, iz_mq", ix_mq, iy_mq,iz_mq
                print*, "q", [grid%kx(ix_q), grid%ky(iy_q), grid%kz(iz_q)]
                print*, "mq", [grid%kx(ix_mq), grid%ky(iy_mq), grid%kz(iz_mq)]
                print*, "m,khi,mu=",m,khi,mu
                print*, a
                print*, b
                print*, a-b
                error stop "La symétrie 1.15 n'est pas vérifiée pour deltarho_prime_k_p"
              end if

            end do

          end do
        end do
      end do
    end block


p = p3%p(1,1,1)
print*, "DELTARHO GRRR", deltarho_prime_k_p(p,2,3,4) ; stop "ok"


    !
    ! Comment on traite les deltarho(k=0)
    !
    ! if (mmax>=1) then
    !   ! on s'interesse a la projection 1 0 0 pour q=0
    !   ix_q=1
    !   iy_q=1
    !   iz_q=1
    !   deltarho_prime_k_p( p3%p(1:,0,0) ,ix_q,iy_q,iz_q) = complex(0,0)
    ! end if





    !
    ! Read Luc's c of k
    !
    if (.not.cq%isok) then
      call cpu_time(time(1))
      call read_ck_nonzero
      call cpu_time(time(2))
      print*, "<<<<< lecture des c Luc en",time(2)-time(1),"seconds"
    end if


    block
      integer :: m, n, mu, nu, khi
      print*, "opening /home/levesque/Recherche/00__BELLONI/2016__02/reconstruction_des_ck_depuis_luc_ck_nmax5/ck_nmax5_max"
      open(89, file="/home/levesque/Recherche/00__BELLONI/2016__02/reconstruction_des_ck_depuis_luc_ck_nmax5/ck_nmax5_max")
        do khi=-mmax,mmax
          do m=abs(khi),mmax
            do mu=-m,m
              if (mod(mu,2)/=0) cycle
              do n=abs(khi),mmax
                do nu=-n,n
                  if (mod(nu,2)/=0) cycle
                  read(89,*) ptest(m,n,mu,nu,khi), cfulltest(ptest(m,n,mu,nu,khi),:)
                end do
              end do
            end do
          end do
        end do
      close(89)
      print*, "closed /home/levesque/Recherche/00__BELLONI/2016__02/reconstruction_des_ck_depuis_luc_ck_nmax5/ck_nmax5_max"
    end block



    block
      use module_wigner_d, only: symbol_3j
      integer :: mmax2, khi, ialp, m, n, l, mu, nu, mu2, nu2, i, m2, n2, coeff
      integer, parameter :: ialpmax = 549
      integer, parameter :: n_baratin = 5
      complex(dp) :: coeff_cmplx
      complex(dp), parameter :: xi_cmplx=cmplx(0,1,dp)
      integer, allocatable, dimension(:) :: mm, nn, ll, mumu, nunu
      character(50) :: texte
      character(10) :: texte1
      double precision :: q
      REAL(8), allocatable :: ck(:,:)
      complex(dp), allocatable :: tab6(:,:,:,:,:   ,:) ! m n l mu2 nu2  iq
      complex(dp), allocatable :: tab5(:,:,:,:,:   ,:) ! m n khi, mu2 nu2    iq
      integer :: mnmax, mnmax2
      integer, parameter :: nq=1024
      mnmax=mmax
      mmax2=mmax/2
      mnmax2=mmax2
      !
      !                          on lit les projections mnlmunu de ck
      !
      PRINT*, '********************************************************************************'
      PRINT*, 'Lecture des cmnlmunu(q) dans un fichier'
      PRINT*, 'ce fichier contient apres quelques lignes de baratin la ligne des m puis celle des n...'
      PRINT*, 'suivies des lignes de q,cmnlmunu(q) (on choisira le q le plus proche, par exces)'
      PRINT*, 'attention: si l pair, fonction reelle; si l impair, fonction imaginaire pure donc i implicite!'
      print*, "reading /home/levesque/Recherche/00__BELLONI/2016__FEVRIER__OZ/ck_h2o39-3916_l.txt"
      OPEN(7, file="/home/levesque/Recherche/00__BELLONI/2016__FEVRIER__OZ/ck_h2o39-3916_l.txt", STATUS='old')
      ALLOCATE(mm(ialpmax),nn(ialpmax),ll(ialpmax),mumu(ialpmax),nunu(ialpmax),ck(ialpmax,nq))
      do i=1,n_baratin
        READ(7,*) texte
      end do
      READ(7,*) texte1,mm
      PRINT*, texte1,mm(1:10),"..."
      READ(7,*) texte1,nn
      PRINT*, texte1,nn(1:10),"..."
      READ(7,*) texte1,ll
      PRINT*, texte1,ll(1:10),"..."
      READ(7,*) texte1,mumu
      PRINT*, texte1,mumu(1:10),"..."
      READ(7,*) texte1,nunu
      PRINT*, texte1,nunu(1:10),"..."
      ! do i=1,10000
      do i=1,nq
        READ(7,*) q,ck(:,i)
      end do
      CLOSE(7)
      print*, "finished reading /home/levesque/Recherche/00__BELLONI/2016__FEVRIER__OZ/ck_h2o39-3916_l.txt"
      !   je construis en tableau a 5 entrees cmnlmunu pour tests ulterieurs
      ALLOCATE (tab6(0:mmax,0:mmax,0:2*mmax,-mmax2:mmax2,-mmax2:mmax2,nq))
      tab6=(0,0)
      do iq=1,nq
        do ialp=1,ialpmax
          m=mm(ialp)
          n=nn(ialp)
          l=ll(ialp)
          mu=mumu(ialp)
          nu=nunu(ialp)
          mu2=mu/2
          nu2=nu/2
          coeff_cmplx= ck(ialp,iq) ! ck est reel
          IF((-1)**l==-1) coeff_cmplx=xi_cmplx*ck(ialp,iq)
          tab6(m,n,l,mu2,nu2,iq)=coeff_cmplx
          tab6(m,n,l,-mu2,-nu2,iq)=(-1)**(m+n+l)*coeff_cmplx
          tab6(n,m,l,nu2,mu2,iq)=(-1)**(m+n)*coeff_cmplx
          tab6(n,m,l,-nu2,-mu2,iq)=(-1)**l*coeff_cmplx
        end do
      end do
    !
    !   !
    !   !                           on a donc tout ce qui nous faut pour calculer le produit de convolution
      allocate (tab5(0:mmax,0:mmax,-mmax:mmax,-mmax2:mmax2,-mmax2:mmax2 ,nq))
      tab5 = (0,0)
      do iq=1,nq
        do khi=0,mmax                ! que khi>=0 pour l'instant
          do ialp=1,ialpmax
            m=mm(ialp)
            n=nn(ialp)
            l=ll(ialp)
            mu=mumu(ialp)
            nu=nunu(ialp)
            if ( khi > min(m,n) ) cycle
            mu2=mu/2
            nu2=nu/2
            coeff_cmplx = symbol_3j(m,n,l,khi,-khi,0) * ck(ialp,iq)
            IF((-1)**l==-1) coeff_cmplx = xi_cmplx*coeff_cmplx             ! imaginaire pur si l impair
            tab5(m,n,khi,mu2,nu2,iq)=tab5(m,n,khi,mu2,nu2,iq)+coeff_cmplx                     ! tab5 est donc complexe
            IF(mu/=0.or.nu/=0) tab5(m,n,khi,-mu2,-nu2,iq)=tab5(m,n,khi,-mu2,-nu2,iq)+(-1)**(m+n+l)*coeff_cmplx
            IF(m/=n.or.ABS(mu)/=ABS(nu)) then
              tab5(n,m,khi,nu2,mu2,iq)=tab5(n,m,khi,nu2,mu2,iq)+(-1)**(m+n)*coeff_cmplx
              IF(mu/=0.or.nu/=0) tab5(n,m,khi,-nu2,-mu2,iq)=tab5(n,m,khi,-nu2,-mu2,iq)+(-1)**(l)*coeff_cmplx
            end if
          end do                           ! fin ialp
        end do                             ! fin khi>=0
      end do
      !                            et je complete les khi<0 avec cmnmunu-khi=(-1)**(m+n)*cmnmunukhi*
      do iq=1,nq
        do m=0,mmax
          m2=m/2
          do n=0,mmax
            n2=n/2
            coeff=(-1)**(m+n)
            do khi=-MIN(m,n),-1
              tab5(m,n,khi,-m2:m2,-n2:n2, iq)=coeff*CONJG(tab5(m,n,-khi,-m2:m2,-n2:n2,iq))
            end do
          end do
        end do
      end do

      ! PRINT*, 'test a partir de tableaux a 5 entrees cmnmunu'
      !
      ! do iq=1,nq
      !   do m=0,mnmax
      !     do n=0,mnmax
      !       do khi=-mnmax,mnmax
      !         do mu2=-mnmax2,mnmax2
      !           do nu2=-mnmax2,mnmax2
      !             coeff_cmplx=0.
      !             do l=0,2*mnmax
      !               coeff_cmplx=coeff_cmplx+symbol_3j(m,n,l,khi,-khi,0)*tab6(m,n,l,mu2,nu2,iq)
      !             end do
      !             IF(ABS(tab5(m,n,khi,mu2,nu2,iq)-coeff_cmplx)>1.d-10) &
      !               PRINT*, m,n,khi,2*mu2,2*nu2,coeff_cmplx,tab5(m,n,khi,mu2,nu2,iq),"ERROR"
      !           end do
      !         end do
      !       end do
      !     end do
      !   end do
      ! end do


    ! do iq=1,100
    !   do khi=-mmax,mmax
    !     do m=abs(khi),mmax
    !       do mu=-m,m
    !         do n=abs(khi),mmax
    !           do nu=-n,n
    !             if (mod(mu,2)/=0 .or. mod(nu,2)/=0) cycle
    !             a = tab5(m,n,khi,mu/2,nu/2, iq) ! c m n mu nu khi reconstruit depuis c m n l mu nu de luc
    !             b = cfulltest(ptest(m,n,mu,nu,khi),iq)
    !   if (abs(a-b)>1.D-3) then
    !     print*, iq, (iq-1)*(0.613592315E-01), m, n, mu, nu, khi, a, b, a-b
    !   end if
    !           end do
    !         end do
    !       end do
    !     end do
    !   end do
    ! end do

    end block






    !
    ! OZ
    !
    call cpu_time(time(1))
    allocate (gamma_prime_k_p(np,nx,ny,nz), source=zeroc)
    do iz_q=1,nz
      do iy_q=1,ny
        do ix_q=1,nx
          !
          ! iq = int( norm2([grid%kx(ix_q), grid%ky(iy_q), grid%kz(iz_q)]) /cq%dq +0.5) +1
          ! if (iq >cq%nq) error stop "we need larger norm of q in Luc's cq"

          iqtest = int( norm2([grid%kx(ix_q), grid%ky(iy_q), grid%kz(iz_q)]) /(0.613592315E-01) +0.5) +1
          if (iqtest>200) error stop "iqtest trop grand"

          do khi=-mmax,mmax
            do m=abs(khi),mmax
              do mu=-m,m

                p=p3%p(m,khi,mu)

                if (mod(mu,2)/=0) cycle!then
                !   gamma_prime_k_p(p3%p(m,khi,mu),ix_q,iy_q,iz_q) = complex(0,0)
                !
                ! end if

                do n=abs(khi),mmax
                  do nu=-n,n
                    if (mod(nu,2)/=0) cycle


                    ! ia = cq%a(m,n,mu,nu,khi)
                    ! if (ia<0) then
                    !   print*, "error : ia <0"
                    !   print*, "m,n,mu,nu,khi=",m,n,mu,nu,khi
                    !   error stop
                    ! end if

  ! gamma_prime_k_p(p3%p(m,khi,mu),ix_q,iy_q,iz_q) = gamma_prime_k_p(p3%p(m,khi,mu),ix_q,iy_q,iz_q) &
  !                 + (-1)**(khi+nu) * ck(ia,iq) * deltarho_prime_k_p(p3%p(n,khi,-nu),ix_q,iy_q,iz_q)

  gamma_prime_k_p(p,ix_q,iy_q,iz_q) = gamma_prime_k_p(p,ix_q,iy_q,iz_q) &
          + (-1)**(khi+nu) * cfulltest(ptest(m,n,mu,nu,khi),iqtest) * deltarho_prime_k_p(p3%p(n,khi,-nu),ix_q,iy_q,iz_q)

                  end do
                end do

              end do
            end do
          end do

        end do
      end do
    end do
    call cpu_time(time(2))
    print*, "<<<<< OZ en",time(2)-time(1),"seconds"
!
! verifie que toutes les projections à mu impairs sont nuls
!
do iz_q=1,nz
  do iy_q=1,ny
    do ix_q=1,nx
      do m=0,mmax
        do khi=-m,m
          do mu=-m,m
            a = gamma_prime_k_p(p3%p(m,khi,mu),ix_q,iy_q,iz_q)
            if (mod(mu,2)/=0 .and. abs(a)>1.D-10) then
              print*,"q", ix_q, iy_q, iz_q,"m khi mu",m,khi,mu,"gammaOZ", gamma_prime_k_p(p3%p(m,khi,mu),ix_q,iy_q,iz_q)
              error stop "OZ faux"
            end if
          end do
        end do
      end do
    end do
  end do
end do
    !
    ! Vérification de la relation de symétrie 1.25, c'est à dire la meme que 1.15 en notation χ
    !
    ! block
    !   complex(dp) :: a, b
    !   integer :: ix_q, iy_q, iz_q, ix_mq, iy_mq, iz_mq
    !
    !   do iz_q=1,nz
    !     do iy_q=1,ny
    !       do ix_q=1,nx
    !
    !         do p=1,np
    !
    !           m = p3%m(p)
    !           khi = p3%mup(p)
    !           mu = p3%mu(p)
    !
    !           !
    !           ! on cherche le l'indice correspond à -q
    !           !
    !           if (ix_q == 1) then
    !             ix_mq = ix_q
    !           else if (ix_q == nx/2 +1) then
    !             cycle
    !           else
    !             ix_mq = nx - ix_q + 2
    !           end if
    !
    !           if (iy_q == 1) then
    !             iy_mq = iy_q
    !           else if (iy_q == ny/2 +1) then
    !             cycle
    !           else
    !             iy_mq = ny - iy_q + 2
    !           end if
    !
    !           if (iz_q == 1) then
    !             iz_mq = iz_q
    !           else if (iz_q == nz/2 +1) then
    !             cycle
    !           else
    !             iz_mq = nz - iz_q + 2
    !           end if
    !
    !
    !           ! ATTENTION
    !           ! ICI LA NOTION DE Q ET -Q DEVRAIENT CHANGER PHI EN PHI+PI DANS LE CAS DE Q SELON Z
    !           a = gamma_prime_k_p (p3%p(m,khi,mu),ix_mq,iy_mq,iz_mq)
    !           a = conjg(a)*(-1)**(m+mu+khi)
    !           b = gamma_prime_k_p (p3%p(m,khi,-mu),ix_q,iy_q,iz_q)
    !
    !           ! if (abs(a-b) >1.E-10 ) then
    !             ! if (ix_q/=1 .and. iy_q/=1 .and. iz_q/=1 ) then!cycle ! VOIR ATTENTION JUSTE AU DESSUS !!!!! TODO ASK LUC DANIEL LU
    !             print*,
    !             print*, "ix_q, iy_q, iz_q",ix_q, iy_q, iz_q
    !             print*, "ix_mq, iy_mq, iz_mq", ix_mq, iy_mq,iz_mq
    !             print*, "q", [grid%kx(ix_q), grid%ky(iy_q), grid%kz(iz_q)]
    !             print*, "mq", [grid%kx(ix_mq), grid%ky(iy_mq), grid%kz(iz_mq)]
    !             print*, "m,khi,mu=",m,khi,mu
    !             print*, "conjg @ gamma_prime_k_p (p3%p(m,khi,mu),ix_mq,iy_mq,iz_mq)*(-1)**(m+mu+khi),  gamma m khi -mu"
    !             print*,a
    !             print*,b
    !             print*, "diff",a-b
    !             ! error stop "La symétrie 1.15 n'est pas vérifiée pour gamma_prime_k_p"
    !           ! end if!; end if
    !
    !         end do
    !
    !       end do
    !     end do
    !   end do
    ! end block



    do iz_q=1,nz
      do iy_q=1,ny
        do ix_q=1,nx
          do m=0,mmax
            do khi=-m,m
              do mu=-m,m
if(mod(mu,2)/=0 .and. abs(gamma_prime_k_p(p3%p(m,khi,mu),ix_q,iy_q,iz_q))>epsilon(1._dp)) then
       print*,"q", ix_q, iy_q, iz_q,"m khi mu",m,khi,mu,"gammaT", gamma_prime_k_p(p3%p(m,khi,mu),ix_q,iy_q,iz_q)
       error stop "hasoijazdui"
end if
              end do
            end do
          end do
        end do
      end do
    end do


    !
    ! Retour dans le repère fixe
    !
    call cpu_time(time(1))
    block
      complex(dp), allocatable :: R(:,:,:) ! R^m_{mup,mu} (\hat{q})
      allocate ( R(0:mmax,-mmax:mmax,-mmax:mmax) ,source=zeroc)

      do iz_q=1,nz
        do iy_q=1,ny
          do ix_q=1,nx

            q = [grid%kx(ix_q), grid%ky(iy_q), grid%kz(iz_q)]
            R = conjg(rotation_matrix_between_complex_spherical_harmonics_lu ( q, mmax ))

            do m=0,mmax
              do mup=-m,m
                do mu=-m,m

                  gamma_k_p( p3%p(m,mup,mu) ,ix_q,iy_q,iz_q) = sum(gamma_prime_k_p(p3%p(m,-m:m,mu),ix_q,iy_q,iz_q) * R(m,mup,-m:m))

                end do
              end do
            end do

          end do
        end do
      end do

    end block
    deallocate (gamma_prime_k_p)
    call cpu_time(time(2))
    print*, "<<<<< retour repère fixe en",time(2)-time(1),"seconds"


    do iz_q=1,nz
      do iy_q=1,ny
        do ix_q=1,nx
          do m=0,mmax
            do khi=-m,m
              do mu=-m,m
if(mod(mu,2)/=0 .and. abs(gamma_k_p(p3%p(m,khi,mu),ix_q,iy_q,iz_q))>epsilon(1._dp)) then
       print*,"q", ix_q, iy_q, iz_q,"m khi mu",m,khi,mu,"gammaT", gamma_k_p(p3%p(m,khi,mu),ix_q,iy_q,iz_q)
       error stop "hasoijazdui"
end if
              end do
            end do
          end do
        end do
      end do
    end do


    !
    ! Vérification de la symétrie dûe a la fonction de départ purement réelle
    ! cf equations 1.15 de Luc
    !
    block
      complex(dp) :: a, b
      integer :: ix_q, iy_q, iz_q, ix_mq, iy_mq, iz_mq

      do iz_q=1,nz
        do iy_q=1,ny
          do ix_q=1,nx

            do p=1,np

              m = p3%m(p)
              mup = p3%mup(p)
              mu = p3%mu(p)

              !
              ! on cherche le l'indice correspond à -q
              !
              if (ix_q == 1) then
                ix_mq = ix_q
              else if (ix_q == nx/2 +1) then
                cycle
              else
                ix_mq = nx - ix_q + 2
              end if

              if (iy_q == 1) then
                iy_mq = iy_q
              else if (iy_q == ny/2 +1) then
                cycle
              else
                iy_mq = ny - iy_q + 2
              end if

              if (iz_q == 1) then
                iz_mq = iz_q
              else if (iz_q == nz/2 +1) then
                cycle
              else
                iz_mq = nz - iz_q + 2
              end if

              a = gamma_k_p (p3%p(m, mup, mu),ix_mq,iy_mq,iz_mq)
              a = conjg(a)*(-1)**(mup+mu)
              b = gamma_k_p (p3%p(m,-mup,-mu),ix_q,iy_q,iz_q)

              if (abs(a-b)>1.D-10) then
                print*, "La symétrie 1.15 n'est pas vérifiée pour gamma^m_mup_mu(q)"
                print*, "ix, iy ,iz, m, mup, mu",ix_q,iy_q,iz_q,m,mup,mu
                print*, "a",a
                print*, "b",b
                print*, "a-b",a-b
              end if

            end do

          end do
        end do
      end do
    end block




    !
    ! Retour en angles
    !
    call cpu_time(time(1))
    allocate( gamma_k(no,nx,ny,nz) ,source=zeroc )
    block
      complex(dp), allocatable :: temp_wigner_big_D(:,:)
      real(dp), allocatable :: fm_of_p(:)
      integer :: o
      allocate (fm_of_p(np))
      do p=1,np
        fm_of_p(p) = fm(p3%m(p))
      end do
      allocate (temp_wigner_big_D(np,no), source=zeroc)
      do o=1,no
        do p=1,np
          temp_wigner_big_D(p,o) = wigner_big_D(p3%m(p),p3%mup(p),p3%mu(p),grid%theta(o),grid%phi(o),grid%psi(o))
        end do
      end do

      do iz=1,nz
        do iy=1,ny
          do ix=1,nx
            do o=1,no

              gamma_k(o,ix,iy,iz) = sum( fm_of_p(:) * temp_wigner_big_D(:,o) * gamma_k_p(:,ix,iy,iz) )

            end do
          end do
        end do
      end do

      deallocate (temp_wigner_big_D)
    end block
    deallocate (gamma_k_p)
    call cpu_time(time(2))
    print*, "<<<<< retour en angles en",time(2)-time(1),"seconds"


    !
    ! gamma dans l'espace reel doit être une fonction réelle. On peut donc vérifier la symétrie hermitique sur q
    ! En fait, on retrouve bien un pure réel par FFT-1 donc cette vérification est un peu inutile.
    !
    block
      complex(dp) :: a, b
      integer :: ix_q, iy_q, iz_q, ix_mq, iy_mq, iz_mq

      do iz_q=1,nz
        do iy_q=1,ny
          do ix_q=1,nx

            do io=1,no

              !
              ! on cherche le l'indice correspond à -q
              !
              if (ix_q == 1) then
                ix_mq = ix_q
              else if (ix_q == nx/2 +1) then
                cycle
              else
                ix_mq = nx - ix_q + 2
              end if

              if (iy_q == 1) then
                iy_mq = iy_q
              else if (iy_q == ny/2 +1) then
                cycle
              else
                iy_mq = ny - iy_q + 2
              end if

              if (iz_q == 1) then
                iz_mq = iz_q
              else if (iz_q == nz/2 +1) then
                cycle
              else
                iz_mq = nz - iz_q + 2
              end if

              a = conjg( gamma_k (io,ix_mq,iy_mq,iz_mq)  )
              b =        gamma_k (io,ix_q,iy_q,iz_q)

              if (abs(a-b)>1.D-10) then
                error stop "La symétrie sur q -q du fait de gamma purement réel n'est pas vérifiée"
              end if

            end do

          end do
        end do
      end do
    end block



    ! !
    ! ! Estimation de l'énergie depuis angles et q (eq. 1.33 ligne 3)
    ! !
    ! block
    !   real(dp) :: ff
    !   complex(dp) :: a
    !   ff = 0
    !   do iz_q=1,nz
    !     do iy_q=1,ny
    !       do ix_q=1,nx
    !         a = sum(grid%w * deltarho_k(:,ix_q,iy_q,iz_q)*conjg(gamma_k(:,ix_q,iy_q,iz_q)) )
    !         if( abs(imag(a)) > 1.E-10 ) then
    !           print*, "ff(rvec) devrait être purement réel mais ne l'est pas"
    !           print*, sum(grid%w * deltarho_k(:,ix_q,iy_q,iz_q)*conjg(gamma_k(:,ix_q,iy_q,iz_q)) )
    !         end if
    !         ff = ff - kT/2._dp*real(a,dp)
    !       end do
    !     end do
    !   end do
    !   block
    !     double precision, parameter :: pi=acos(-1._dp),   twopi=2*pi
    !     print*,"ff%exc_cproj calculé depuis l'espace de Fourier=",ff   /real(nx*ny*nz) /solvent(1)%n0! /twopi**3  *kT    /sqrt(real(nx*ny*nz,dp))   *2*pi**2
    !   end block
    !   print*, "dk^3=",grid%kx(2)*grid%ky(2)*grid%kz(2)
    !   print*, "dv=",dv
    ! end block

    !
    ! Transformée de Fourier 3D inverse
    !
    call cpu_time(time(1))
    allocate (gamma(no,nx,ny,nz), source=0._dp)
    block
      complex(dp), allocatable :: in(:,:,:)
      type(c_ptr) :: plan_fft_c2c_3d_signe_minus
      allocate (in(nx,ny,nz))
      call dfftw_plan_dft_3d (plan_fft_c2c_3d_signe_minus, nx, ny, nz, in, in, FFTW_FORWARD, FFTW_ESTIMATE)
      do io=1,no
        in = (0,0)
        in = gamma_k(io,:,:,:)
        call dfftw_execute_dft( plan_fft_c2c_3d_signe_minus, in, in )
        gamma_k(io,:,:,:) = in
      end do
    end block
    ! do iz=1,nz
    !   do iy=1,ny
    !     do ix=1,nx
    !       do io=1,no
    !         print*,"io ix iy iz imag(gamma_k(io,ix,iy,iz))", io,ix,iy,iz,imag(gamma_k(io,ix,iy,iz))
    !       end do
    !     end do
    !   end do
    ! end do
    !
    ! if ( any(  imag(gamma_k)   >1.D-10) ) then
    !   error stop "gamma_k has some imaginary component somewhere"
    ! end if
    gamma = real(gamma_k,dp)/real(nx*ny*nz,dp)
    deallocate (gamma_k)
    call cpu_time(time(2))
    print*, "<<<<< retour en r par FFT-1 en",time(2)-time(1),"seconds"



    !
    ! calcul de l'énergie et du gradient
    !
    call cpu_time(time(1))
    ff=0._dp
    do iz=1,nz
      do iy=1,ny
        do ix=1,nx
          vexc(1:no) = -kT*grid%w(1:no)*gamma(1:no,ix,iy,iz)   /solvent(1)%n0
          do io=1,no
            xi = solvent(1)%xi(io,ix,iy,iz)
            rho = xi**2*rho0
            ff = ff + (rho-rho0)*0.5_dp*vexc(io)*dv
            df(io,ix,iy,iz,1) = df(io,ix,iy,iz,1) + 2._dp*vexc(io)*rho0*xi
          end do
        end do
      end do
    end do
    call cpu_time(time(2))
    print*, "<<<<< calcul de ff et df en",time(2)-time(1),"seconds"


  end subroutine energy_cproj_no_symetry


! CI DESSOUS SEULEMENT LA LECTURE DU C DE LUC DANS MON FORMAT.


!!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX



    subroutine read_ck_nonzero
      use module_input, only: n_linesInFile
      use module_grid, only: grid
      implicit none
      integer :: i, iq, m, n, mu, nu, khi, ia
      character(3) :: somechar
      integer :: ufile, ios
      character(65) :: filename
      integer, parameter, dimension(0:5) :: nprojections_for_mmax = [1,6,75,252,877,2002]
      real(dp) :: qmax_effectif
      integer :: iqmax_effectif
      complex(dp), allocatable :: ck_full(:,:)


      if (cq%isok) return

      write(filename,'(a,i0,a)') "input/direct_correlation_functions/water/SPCE/ck_nonzero_nmax",grid%mmax,"_ml"

      select case (grid%mmax)
      case (0,1,2,3,4,5)
        cq%na = nprojections_for_mmax(grid%mmax)
      case default
        print*, "In module_energy_cproj, you want mmax>5 but I don't have the corresponding file."
        error stop
      end select

      ! number of q points in c_mnmunukhi(q)
      cq%nq = n_linesInFile(filename) - 17 ! there are 17 lines of "header" before the columns corresponding to cmnmunukhi(q)

      ! open the file that contains c_mnmunukhi(q)
      open(newunit=ufile, file=filename, iostat=ios, status="old", action="read")
      if ( ios /= 0 ) then
        print*, "Cant open file", filename
        error stop "in module_energy_cproj.f90"
      end if

      if (.not.allocated(cq%normq)) allocate(cq%normq(cq%nq), source=0._dp)


      ! i=(1+mmax)*(1+2*mmax)*(3+2*mmax)*(5+4*mmax*(2+mmax))/15
      !  if (cq%na/=i) then
      !      print*, "in read_ck_cproj nalpha est bizarre"
      !      print*, "cq%na=",cq%na
      !      print*, "it should be (1+mmax)*(1+2*mmax)*(3+2*mmax)*(5+4*mmax*(2+mmax))/15 =",i
      !      ! error stop
      !  end if

      ! i=sum([([([([([( 1 ,nu=-n,n)] ,mu=-m,m)], n=abs(khi),mmax)] ,m=abs(khi),mmax)] ,khi=-mmax,mmax)]  )
      !  if (cq%na/=i) then
      !      print*, "in read ck cq%na /= nalpha"
      !      print*, "cq%na=",cq%na
      !      print*, "bruteforce=",i
      !      ! error stop
      !  end if

      if (allocated(ck) .and. .not.cq%isok) then
        print*, "ck is already allocated but .not. cq%isok"
        error stop
      end if

      if (.not.allocated( cq%m)) then
        allocate ( cq%m  (cq%na) ,source=-huge(1))
        allocate ( cq%n  (cq%na) ,source=-huge(1))
        allocate ( cq%mu (cq%na) ,source=-huge(1))
        allocate ( cq%nu (cq%na) ,source=-huge(1))
        allocate ( cq%khi(cq%na) ,source=-huge(1))
        block
          integer :: mmax
          mmax = grid%mmax
          allocate ( cq%a(0:mmax,0:mmax,-mmax:mmax,-mmax:mmax,-mmax:mmax), source=-huge(1)) ! m n mu nu khi. -huge is used to spot more easily bugs that may come after
        end block
      end if
      !
      ! Skip 10 lines of comments
      !
      do i=1,10
        read(ufile,*)
      end do

      read(ufile,*) somechar
      read(ufile,*) somechar, cq%m
      read(ufile,*) somechar, cq%n
      read(ufile,*) somechar, cq%mu
      read(ufile,*) somechar, cq%nu
      read(ufile,*) somechar, cq%khi
      read(ufile,*)

      do ia=1,cq%na
        m = cq%m(ia)
        n = cq%n(ia)
        mu = cq%mu(ia)
        nu = cq%nu(ia)
        khi = cq%khi(ia)
        cq%a(m,n,mu,nu,khi) = ia
      end do


      allocate( ck_full(cq%na,cq%nq), source=zeroc)
      do iq=1,cq%nq
        read(ufile,*) cq%normq(iq), ck_full(:,iq)
      end do
      close(ufile)
      cq%dq = cq%normq(2)

      !
      ! Plutot que stocker un tableau tres grand avec toutes les valeurs de |q| dont on a besoin,
      ! on va aller voir quel est le |q| maximum dont on a besoin dans notre code (depend de nx,ny,nz)
      ! cette valuer maximum effective de |q| dans le code, on l'appelle qmax_effectif
      ! on appelle son indice : iqmax_effectif
      !
      qmax_effectif = maxval ( sqrt(grid%kx**2 + grid%ky**2 + grid%kz**2)  )
      iqmax_effectif = int( qmax_effectif / cq%dq  ) + 10 ! +10 is just to be safe. +1 suffit très certainement.
      cq%nq = iqmax_effectif
      allocate( ck(cq%na, cq%nq), source=ck_full(1:cq%na,1:cq%nq))
      deallocate( ck_full )

      deallocate (cq%m, cq%n, cq%mu, cq%nu, cq%khi)

      ! tell the world the reading and allocating of the data structure for c(q) is ok
      cq%isok=.true.

    end subroutine read_ck_nonzero



















    pure subroutine verify_hermitian_symetry (arrayofq, resultat)
    !
    ! Vérification de la symétrie hermitienne dûe a la fonction de départ réelle
    !
      implicit none
      complex(dp), intent(in) :: arrayofq(:,:,:)
      logical, intent(out) :: resultat
      complex(dp) :: a, b
      integer :: ix_q, iy_q, iz_q, ix_mq, iy_mq, iz_mq, nx, ny, nz

      resultat = .true.

      nx=ubound(arrayofq,1)
      ny=ubound(arrayofq,2)
      nz=ubound(arrayofq,1)

      do iz_q=1,nz
        do iy_q=1,ny
          do ix_q=1,nx


            if (ix_q == 1) then
              ix_mq = ix_q
            else if (ix_q == nx/2 +1) then
              cycle
            else
              ix_mq = nx - ix_q + 2
            end if

            if (iy_q == 1) then
              iy_mq = iy_q
            else if (iy_q == ny/2 +1) then
              cycle
            else
              iy_mq = ny - iy_q + 2
            end if

            if (iz_q == 1) then
              iz_mq = iz_q
            else if (iz_q == nz/2 +1) then
              cycle
            else
              iz_mq = nz - iz_q + 2
            end if

            a = arrayofq(ix_q,iy_q,iz_q)
            b = conjg(arrayofq(ix_mq,iy_mq,iz_mq))

            if (abs(a-b)>1.D-10) then
              resultat = .false.
            end if

          end do
        end do
      end do

    end subroutine verify_hermitian_symetry





    pure subroutine verify_hermitian_symetry_of_projections(arrayp ,resultat)
      !
      ! Vérification de la symétrie dûe a la fonction de départ purement réelle
      ! cf equations 1.15 de Luc
      !
      use module_grid, only: grid
      implicit none
      complex(dp), intent(in) :: arrayp(:,:,:,:)
      logical, intent(out) :: resultat
      complex(dp) :: a, b
      integer :: ix_q, iy_q, iz_q, ix_mq, iy_mq, iz_mq, nx, ny, nz, p, m, mu, mup, np

      np=ubound(arrayp,1)
      nx=ubound(arrayp,2)
      ny=ubound(arrayp,3)
      nz=ubound(arrayp,4)


      resultat=.true.

      do iz_q=1,nz
        do iy_q=1,ny
          do ix_q=1,nx

            do p=1,np

              m = p3%m(p)
              mup = p3%mup(p)
              mu = p3%mu(p)

              !
              ! on cherche le l'indice correspond à -q
              !
              ix_mq = grid%ix_mq(ix_q)
              iy_mq = grid%iy_mq(iy_q)
              iz_mq = grid%iz_mq(iz_q)
              if ( ix_mq<=0 .or. iy_mq <=0 .or. iz_mq <=0 ) cycle ! this covers the case of ix_q = nx/2 +1 which has no "-q"

              a = arrayp (p3%p(m, mup, mu),ix_mq,iy_mq,iz_mq)
              a = conjg(a)*(-1)**(mup+mu)
              b = arrayp (p3%p(m,-mup,-mu),ix_q,iy_q,iz_q)

              if (abs(a-b)/=0) then
                resultat=.false.
              end if

            end do

          end do
        end do
      end do
    end subroutine verify_hermitian_symetry_of_projections

end module module_energy_cproj_no_symetry

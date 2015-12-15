module module_energy_cproj_mrso

    use iso_c_binding
    use precision_kinds, only: dp
    use module_grid, only: grid
    use module_solvent, only: solvent
    use module_rotation

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
    complex(dp), parameter :: ii=complex(0._dp,1._dp)
    real(dp), parameter :: epsdp=epsilon(1._dp)
    real(dp), parameter :: pi=acos(-1._dp)
    real(dp), parameter :: eightpisq=8._dp*pi**2
    real(dp), parameter :: twopi=2._dp*pi



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

    real(dp), allocatable :: fm(:)

    type :: fft3d_c_type
        type(c_ptr) :: plan_sens_minus
        type(c_ptr) :: plan_sens_plus
        logical :: is_alread_planned = .false.
    end type fft3d_c_type
    type (fft3d_c_type) :: fft3d

    integer, allocatable :: tableof_ix_mq(:), tableof_iy_mq(:), tableof_iz_mq(:)
    real(dp), allocatable :: qx(:), qy(:), qz(:) ! vector q and its components tabulated

    complex(dp), allocatable, target :: deltarho_p(:,:,:,:) ! deltarho_p(np,nx,ny,nz)
    complex(dp), allocatable :: deltarho_p_q(:)
    complex(dp), allocatable :: deltarho_p_mq(:)
    complex(dp), allocatable :: gamma_p_q(:)
    complex(dp), allocatable :: gamma_p_mq(:)
    real(dp), allocatable :: gamma_o(:)

    complex(dp), allocatable :: foo_theta_mu_mup(:,:,:)


    type :: p3_type
        real(dp), allocatable :: harm_sph(:,:) ! tabulation des harmoniques sphériques r(m,mup,mu,theta) en un tableau r(itheta,p)
        integer, allocatable :: p(:,:,:) ! index of the projection corresponding to m, mup, mu
        integer, allocatable :: m(:) ! m for projection 1 to np
        integer, allocatable :: mup(:) ! mup for projection 1 to np. mup corresponds to phi
        integer, allocatable :: mu(:) ! mu for projection 1 to np. mu corresponds to psi
        complex(dp), allocatable :: foo_q(:) ! foo (:) is a temporary array of size np
        complex(dp), allocatable :: foo_mq(:) ! foo (:) is a temporary array of size np
    end type p3_type
    type (p3_type) :: p3


    type :: fft2d_r2c_type
        type(c_ptr) :: plan
        real(dp), allocatable    :: in(:,:)
        complex(dp), allocatable :: out(:,:)
        logical :: is_already_planned
    end type
    type :: fft2d_c2r_type
        type(c_ptr) :: plan
        complex(dp), allocatable :: in(:,:)
        real(dp), allocatable    :: out(:,:)
        logical :: is_already_planned
    end type
    type( fft2d_r2c_type ), save :: fft2d_r2c
    type( fft2d_c2r_type ), save :: fft2d_c2r

    complex(dp), allocatable :: R(:,:,:) ! Table of generalized spherical harmonics of m, mup, mu


    public :: energy_cproj_mrso

contains

subroutine energy_cproj_mrso (ff,df)

    use iso_c_binding, only: c_ptr, dp=>c_double
    use module_grid, only: grid
    use module_thermo, only: thermo

    implicit none

    real(dp), intent(out) :: ff
    real(dp), intent(inout) :: df(:,:,:,:,:) ! x y z o s
    real(dp) :: dv, kT
    logical :: q_eq_mq
    integer :: ix, iy, iz, ix_q, iy_q, iz_q, ix_mq, iy_mq, iz_mq, i, p, ip2
    integer :: nx, ny, nz, np, no, ns, ntheta, nphi, npsi, mmax, na, nq, mrso
    integer :: m, n, mu, nu, khi, mup, ia, ip, ipsi, iphi, itheta, iq
    complex(dp) :: gamma_m_khi_mu_q, gamma_m_khi_mu_mq
    real(dp) :: q(3), mq(3), lx, ly, lz, rho0
    real(dp) :: theta(grid%ntheta), wtheta(grid%ntheta)
    logical, allocatable :: gamma_p_isok(:,:,:)
    real :: time(20)

call cpu_time (time(1))

        ntheta=grid%ntheta
        nphi=grid%nphi
        npsi=grid%npsi

        if (.not.allocated(fft2d_r2c%in))  allocate (fft2d_r2c%in (npsi,nphi))
        if (.not.allocated(fft2d_r2c%out)) allocate (fft2d_r2c%out(npsi/2+1,nphi)) ! pay attention, for fft2d we have psi as the first index, while it is the second index everywhere else.
        if (.not.allocated(fft2d_c2r%in))  allocate (fft2d_c2r%in (npsi/2+1,nphi)) ! this is because of our choice of doing half of mu thanks to hermitian symetry
        if (.not.allocated(fft2d_c2r%out)) allocate (fft2d_c2r%out(npsi,nphi))

        mmax=grid%mmax
        mrso=grid%molrotsymorder

        if (.not. allocated( foo_theta_mu_mup) ) &
            allocate ( foo_theta_mu_mup(1:ntheta,0:mmax/mrso,-mmax:mmax))

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
        ! MDFT-dev only valid for even number of nodes in each direction
        !
        if (mod(nx,2)/=0 .or. mod(ny,2)/=0 .or. mod(nz,2)/=0) then
            print*, "mdft-dev wants even grid nodes"
            print*, "nx,ny,nz=",nx,ny,nz
            error stop
        end if

        if (.not.allocated(fm)) allocate (fm(0:mmax) ,source= [( sqrt(real(2*m+1,dp)) ,m=0,mmax  )])

        if (.not. allocated(p3%foo_q)) allocate ( p3%foo_q(1:np) )
        if (.not. allocated(p3%foo_mq)) allocate ( p3%foo_mq(1:np) )


        if (.not. allocated (deltarho_p) ) then
            allocate (deltarho_p(np,nx,ny,nz) ,source=zeroc)
            allocate (deltarho_p_q(np) ,source=zeroc)
            allocate (deltarho_p_mq(np) ,source=zeroc)
            allocate (gamma_p_q(np), source=zeroc)
            allocate (gamma_p_mq(np), source=zeroc)
        end if

        if (.not. allocated (gamma_p_isok) ) allocate (gamma_p_isok(nx,ny,nz), source=.false.)

call cpu_time (time(2))

        ! 1/ get deltarho = rho-rho0
        ! 2/ project deltarho => deltarho_p (use FFT2D-R2C)
        ! 3/ FFT3D-C2C deltarho_p(x) => deltarho_p(q)
        ! 4/ rotate to q frame
        ! 5/ OZ: deltarho_p(q) => gamma_p(q)
        ! 6/ rotate back to fixed frame
        ! 7/ FFT3D-C2D gamma_p(q) => gamma_p(r)
        ! 8/ gather projections: gamma_p(r) => gamma(r)


        !
        ! Toutes les projections sont stockées dans un meme vecteur de taille np
        ! p3%p(m,mup,mu) lie le triplet (m,mup,mu) a l'unique indice de projection ip
        ! p3%m(p) donne inversement le m correspondant à l'indice p dans le tableau des projections
        ! On a donc p3%p(p3%m(p),p3%mup(p),p3%mu(p)) == p
        !
        if (.not.allocated(p3%p)) then
            allocate ( p3%p(0:mmax,-mmax:mmax, 0:mmax/mrso) ,source=-huge(1)) ! Dans p3%p, on met mu2, pas mu
            allocate ( p3%m(np) ,source=-huge(1))
            allocate ( p3%mup(np) ,source=-huge(1))
            allocate ( p3%mu(np) ,source=-huge(1)) ! c'est bien mu, pas mu2
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
                    do mu=0,m,mrso
                        ip=ip+1
                        IF (ip > np) ERROR STOP "p > np at line 166"
                        p3%p(m,mup,mu/mrso) = ip ! on met mu2 dans p3%p, pas mu
                        p3%m(ip) = m
                        p3%mup(ip) = mup
                        p3%mu(ip) = mu ! c'est le vrai mu, pas mu2
                    end do
                end do
            end do
            if (ip /= np) error stop "ip /= np in energy_cproj"
            if ( any(abs(p3%m)>mmax) .or. any(abs(p3%mup)>mmax) .or. any(abs(p3%mu)>mmax) ) then
                print*, "tabulated m, mup or mu have incorrect values"
                error stop
            end if

            !
            ! Print all projections that will be kept in memory
            !
            open(11, file="output/nonzero-projections.out")
            write(11,*)"Non-zero projections:"
            write(11,*)"        index         m          mup          mu"
            write(11,*)"        -----        ---         ---          --"
            do ip=1,np
                write(11,*) ip, p3%m(ip), p3%mup(ip), p3%mu(ip)
            end do
            close(11)
        end if

call cpu_time (time(3))


        !
        ! Initie poids et racines de l'integration angulaire par la methode de Gauss-Legendre
        ! Les racines sont les racines d'un polynomes en cos(theta).
        ! Comme nous voulons les theta correspondant, on en prend l'arccos(cos(theta))
        !
        theta=grid%thetaofntheta
        wtheta=grid%wthetaofntheta


        !
        ! Tabulate generalized spherical harmonics in array p3%harm_sph(theta,proj)
        ! where theta can be any of the GaussLegendre integration roots for theta
        ! where proj is an index related to a tuple {m,mup,mu}
        ! Blum's notation :
        ! m is related to theta
        ! mup is related to phi
        ! mu is related to psi
        ! TOdo: a remplacer par la routine de luc, et utiliser la notation alpha plutot que m,mup,mu a ce moment
        !
        ! call test_routines_calcul_de_Rm_mup_mu_q
        if (.not. allocated(p3%harm_sph)) then
            allocate ( p3%harm_sph(1:ntheta, 1:np) ,source=0._dp)
            do p=1,np
                m = p3%m(p)
                mup = p3%mup(p)
                mu = p3%mu(p)
                do i=1,ntheta
                    p3%harm_sph(i,p) = harm_sph(m,mup,mu,theta(i))
                end do
            end do
        end if

call cpu_time (time(4))

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! PROJECTION DE Δρ(r,ω) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !
        ! On passe tout de suite en projection, dans le repère cartesien
        ! ces projections sont complexes, mais la symetrie hermitienne
        ! permet de ne garder que les mup>=0 ou les mu>=0.
        ! On choisit les mu>=0 (cf doc de Luc)
        !
        if (.not. fft2d_r2c%is_already_planned) then
            call dfftw_plan_dft_r2c_2d(  fft2d_r2c%plan, npsi, nphi,  fft2d_r2c%in,  fft2d_r2c%out, FFTW_EXHAUSTIVE ) ! npsi est en premier indice
            call dfftw_plan_dft_c2r_2d(  fft2d_c2r%plan, npsi, nphi, fft2d_c2r%in, fft2d_c2r%out, FFTW_EXHAUSTIVE )
            ! call dfftw_plan_dft_2d (fft2d_c%plan, npsi, nphi, fft2d_c%in, fft2d_c%out, FFTW_BACKWARD, FFTW_EXHAUSTIVE)
            fft2d_r2c%is_already_planned =.true.
            fft2d_c2r%is_already_planned =.true.
        end if

        if (.not. fft3d%is_alread_planned) then
            call dfftw_plan_dft_3d (fft3d%plan_sens_plus,&
                nx, ny, nz, deltarho_p(1,:,:,:), deltarho_p(1,:,:,:), FFTW_BACKWARD, FFTW_MEASURE) ! TODO CHECK ESTIMATE VS REST & IS IT WORTH CHANGEING THE PLAN FLAG FOR DIFFERENT np ? Certainly!
            call dfftw_plan_dft_3d( fft3d%plan_sens_minus,&
                nx, ny, nz, deltarho_p(1,:,:,:), deltarho_p(1,:,:,:), FFTW_FORWARD, FFTW_MEASURE )
            fft3d%is_alread_planned = .true.
        end if

call cpu_time (time(5))

        !
        ! Projection of delta rho
        !
        do iz=1,nz
            do iy=1,ny
                do ix=1,nx
                    deltarho_p(1:np,ix,iy,iz) = angl2proj( solvent(1)%density(1:no,ix,iy,iz)-solvent(1)%rho0 )
                end do
            end do
        end do
!
!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !
! ! Test que proj2angl( angl2proj( x )) = x
! !
! block
!     real(dp) :: a(no), b(no)
!     real(dp) :: c
!     integer :: io
!     do ix=1,nx
!         do iy=1,ny
!             do iz=1,nz
!                 a = proj2angl(   angl2proj(   solvent(1)%density(1:no,ix,iy,iz)  ))
!                 b = proj2angl(   angl2proj(   a))
!                 do itheta=1,ntheta
!                     do iphi=1,nphi
!                         do ipsi=1,npsi
!                             io=grid%indo(itheta,iphi,ipsi)
!                             print*,"ix,iy,iz=", ix, iy, iz,"itheta iphi ipsi=", itheta,iphi,ipsi, "on veut", a(io),"on a", b(io)
!                         end do
!                     end do
!                 end do
!             end do
!         end do
!     end do
!     stop "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
! end block
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! call cpu_time (time(6))
!
!         block
!             character(50) :: filename
!             filename = "output/rho_000.dat"
!             call output_rdf (real(deltarho_p(p3%p(0,0,0),:,:,:),dp), filename)
!             print*,"printed output/rho_000.dat"
!         end block
!         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! FFT DE Δρ_p(r) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !
        ! On a les projections sur la grille cartesienne
        ! On veut passer dans l'espace de Fourier pour calculer les convolutions spatiales
        ! On fait donc une FFT 3D pour chacune des projections.
        ! Les projections sont complexes, il s'agit donc d'une FFT3D C2C habituelle : Il n'y a pas de symétrie hermitienne.
        !
        do ip=1,np
            call dfftw_execute_dft( fft3d%plan_sens_plus, deltarho_p(ip,:,:,:), deltarho_p(ip,:,:,:) )
        end do

call cpu_time (time(7))

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ROTATION DE Δρ_p(q) VERS LE REPAIRE MOLECULAIRE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !
        ! We now want this projections in Fourier space and lab frame to the molecular frame
        !
        ! First, we need the rotation matrix to transform
        ! a REAL vector (coordinates) in lab frame to a REAL vector in molecular frame
        ! This must be done for each q vector
        !
        !
        ! We now want the projected densities in the molecular frame
        !
        ! First we tabulate the vector q (Fourier vector) corresponding to x, y and z
        ! This is specific to our FFTW's black box and how it expects indices to be ranged.
        ! See function kproj in module_grid
        ! And the inverse index that gives vector -q
        !

        call for_all_q_find_indices_of_mq

call cpu_time (time(8))

        !
        ! Read Luc's direct correlation function c^{m,n}_{mu,nu_,chi}(|q|)
        ! projected on generalized spherical harmonics
        ! in the intermolecular frame
        ! normq is norm of q, |q|, that correspond to the index iq in ck(ia,iq)
        !
        !call read_ck_nmax (ck, normq)
        !call read_ck_toutes_nmax( ck, normq)
        if (.not.cq%isok) then
            call read_ck_nonzero
            cq%isok=.true.
        end if

call cpu_time (time(9))
        if (.not. allocated(R) ) allocate ( R(0:mmax,-mmax:mmax,-mmax:mmax) ,source=zeroc)

        !
        ! For all vectors q and -q handled simultaneously
        !
        do ix_q=1,nx
            do iy_q=1,ny
                do iz_q=1,nz

                    !
                    ! cartesian coordinates of vector q  and mq (-q) in lab frame
                    !

                    !
                    ! Find ix_mq,iy_mq,iz_mq so that q(ix_mq,iy_mq,iz_mq)= -q(ix_q,iy_q,iz_q)
                    !
                    ix_mq = tableof_ix_mq(ix_q)
                    iy_mq = tableof_iy_mq(iy_q)
                    iz_mq = tableof_iz_mq(iz_q)

                    if ( gamma_p_isok(ix_q,iy_q,iz_q) .and. gamma_p_isok(ix_mq, iy_mq,iz_mq) ) cycle

                    q = [qx(ix_q), qy(iy_q), qz(iz_q)]
                    mq = [qx(ix_mq), qy(iy_mq), qz(iz_mq)]

                    if (ix_mq==ix_q .and. iy_mq==iy_q .and. iz_mq==iz_q) then ! this should only happen for ix=1 and ix=nx/2
                        q_eq_mq=.true.
                    else
                        q_eq_mq=.false.
                    end if

                    !
                    ! Move all projections for q and -q to a smaller temporary array: deltarho_p_q (for q) and deltarho_p_mq (for -q)
                    !
                    deltarho_p_q (1:np) = deltarho_p(1:np,ix_q,iy_q,iz_q)
                    deltarho_p_mq(1:np) = deltarho_p(1:np,ix_mq,iy_mq,iz_mq)

                    !
                    ! Prepare R^m_mup_khi(q)
                    !
                    R = rotation_matrix_between_complex_spherical_harmonics_lu ( q, mmax )
                    ! Eq. 1.23 We don't need to compute gshrot for -q. We do q and -q at the same time.

                    !
                    ! Rotation to molecular (q) frame
                    !
                    p3%foo_q = zeroc
                    p3%foo_mq = zeroc
                    do m=0,mmax
                        do khi=-m,m
                            do mu=0,m,mrso
                                ip=p3%p(m,khi,mu/mrso)
                                ! Equation 1.22
                                do mup=-m,m
                                    ip2=p3%p(m,mup,mu/mrso)
                                    p3%foo_q(ip) = p3%foo_q(ip)  + deltarho_p_q(ip2)*R(m,mup,khi)
                                    p3%foo_mq(ip) = p3%foo_mq(ip) + deltarho_p_mq(ip2)*(-1)**m*R(m,mup,-khi)
                                end do
                            end do
                        end do
                    end do
                    deltarho_p_q = p3%foo_q
                    deltarho_p_mq = p3%foo_mq

                    !
                    ! c^{m,n}_{mu,nu,chi}(|q|) is tabulated for 200 values of |q|.
                    ! Find the tabulated value that is closest to |q|. Its index is iq.
                    ! Note |q| = |-q| so iq is the same for both vectors.
                    !
                    iq = min( int(norm2(q)/cq%dq)+1   ,nq)

                    !
                    ! Consider the case of Pz(k=0). In cdeltacd, we impose it is zeroc because it is non defined. Here we do
                    ! the equivalent. One can show Pz(k=0) is equivalent to projection 100 of deltarho
                    !
                    if (mmax>=1 .and. ix_q==1 .and. iy_q==1 .and. iz_q==1) then
                        deltarho_p_q ( p3%p(1,0,0) ) = zeroc
                        deltarho_p_mq( p3%p(1,0,0) ) = zeroc
                    end if

                    ! if (mmax>=2 .and. ix_q==1 .and. iy_q==1 .and. iz_q==1) then
                    !     deltarho_p_q ( p3%p(2,0,0) ) = zeroc
                    !     deltarho_p_mq( p3%p(2,0,0) ) = zeroc
                    ! end if


                    !
                    ! Ornstein-Zernike in the molecular frame
                    ! We do OZ for q and -q at the same time
                    !
                    gamma_p_q(1:np) = zeroc
                    gamma_p_mq(1:np) = zeroc

                    do khi=-mmax,mmax
                        do m=abs(khi),mmax
                            do mu=0,m,mrso ! not -m,m a cause des symetries ! EST CE QU'IL Y A UNE RAISON POUR QUE LES mu IMPAIRES SOIENT NON NULS ICI ?

                                gamma_m_khi_mu_q= zeroc
                                gamma_m_khi_mu_mq = zeroc

                                do n=abs(khi),mmax
                                    do nu= -mrso*(n/mrso), mrso*(n/mrso), mrso   ! imaginons n=3, -n,n,mrso  ferait nu=-3,-1,1,3 mais en faisant /mrso puis *mrso, ça fait -2,0,2 as expected

                                        ia = cq%a(m,n,mu,nu,khi) ! the index of the projection of c(q). 1<=ia<na

                                        ip = p3%p(n,khi,abs(nu)/mrso)

                                        if (nu<0) then ! no problem with delta rho (n, khi, -nu) since -nu>0. Thus, we apply eq. 1.30 directly
                                            gamma_m_khi_mu_q  = gamma_m_khi_mu_q  + (-1)**(khi+nu) *ck(ia,iq) *deltarho_p_q(ip)
                                            gamma_m_khi_mu_mq = gamma_m_khi_mu_mq + (-1)**(khi+nu) *ck(ia,iq) *deltarho_p_mq(ip)
                                        else
                                            gamma_m_khi_mu_q = gamma_m_khi_mu_q  + (-1)**(n) *ck(ia,iq) *conjg(deltarho_p_mq(ip))
                                            gamma_m_khi_mu_mq= gamma_m_khi_mu_mq + (-1)**(n) *ck(ia,iq) *conjg(deltarho_p_q(ip))
                                        end if

                                    end do
                                end do

                                ip=p3%p(m,khi,mu/mrso)
                                gamma_p_q(ip)  = gamma_m_khi_mu_q
                                gamma_p_mq(ip) = gamma_m_khi_mu_mq

                            end do
                        end do
                    end do


                    !
                    ! Rotation from molecular frame to fix laboratory (Fourier) frame
                    !
                    R = conjg(R) ! le passage retour au repaire fixe se fait avec simplement le conjugue complexe de l'harm sph generalisee

                    p3%foo_q = zeroc
                    p3%foo_mq = zeroc
                    do m=0,mmax
                        do mup=-m,m
                            do mu=0,m,mrso
                                ip=p3%p(m,mup,mu/mrso)
                                ! Equation 1.22
                                do khi=-m,m
                                    ip2=p3%p(m,khi,mu/mrso)
                                    p3%foo_q(ip) = p3%foo_q(ip)   + gamma_p_q(ip2)  *R(m,mup,khi)
                                    p3%foo_mq(ip) = p3%foo_mq(ip) + gamma_p_mq(ip2) *(-1)**m*R(m,mup,-khi)
                                    ! p3%foo_mq(ip) = p3%foo_mq(ip) + gamma_p_mq(ip2)*gshrot_mq(m,mup,khi)
                                end do
                                !
                            end do
                        end do
                    end do
                    gamma_p_q = p3%foo_q
                    gamma_p_mq = p3%foo_mq


                    deltarho_p (1:np, ix_q, iy_q, iz_q) = gamma_p_q(1:np)
                    gamma_p_isok(ix_q,iy_q,iz_q)=.true.

                    if( q_eq_mq .and. (ix_q==nx/2+1.or.iy_q==ny/2+1.or.iz_q==nz/2+1)) then
                        deltarho_p(1:np, ix_mq, iy_mq, iz_mq) = conjg(gamma_p_mq(1:np))
                    else
                        deltarho_p(1:np, ix_mq, iy_mq, iz_mq) = gamma_p_mq(1:np)
                    end if
                    gamma_p_isok(ix_mq,iy_mq,iz_mq)=.true.

                end do
            end do
        end do

call cpu_time(time(10))


        if (.not.all(gamma_p_isok.eqv..true.)) then
            print*, "not all gamma_p(projections,ix,iy,iz) have not been computed"
            error stop
        end if

call cpu_time(time(11))

        !
        ! FFT3D from Fourier space to real space
        !
        do ip=1,np
            call dfftw_execute_dft( fft3d%plan_sens_minus, deltarho_p(ip,1:nx,1:ny,1:nz), deltarho_p(ip,1:nx,1:ny,1:nz) )
        end do
        deltarho_p=deltarho_p/real(nx*ny*nz,dp)

call cpu_time(time(12))

        !
        ! Gather projections into gamma (in fact, into the gradient, that IS gamma)
        !
        if (.not. allocated(gamma_o)) allocate (gamma_o(1:no) ,source=0._dp)
        ff=0._dp
        do iz=1,nz
            do iy=1,ny
                do ix=1,nx
                    gamma_o(1:no) = -kT*dv*proj2angl(deltarho_p(1:np,ix,iy,iz))/0.0333_dp*grid%w(:)**2*real(nphi*npsi*ntheta,dp)
                    df(:,ix,iy,iz,1) = df(:,ix,iy,iz,1) + gamma_o(:)
                    ff = ff + sum( gamma_o*(solvent(1)%density(:,ix,iy,iz)-rho0) )/2._dp
                end do
            end do
        end do

call cpu_time(time(13))

print*, "   > allocations                                                  ", time(2)-time(1),"sec"
print*, "   > print projections to file                                    ", time(3)-time(2),"sec"
print*, "   > tabulate spherical harmonics                                 ", time(4)-time(3),"sec"
print*, "   > plan all ffts                                                ", time(5)-time(4),"sec"
print*, "   > project delta rho                                            ", time(6)-time(5),"sec"
print*, "   > FFT 3D of deltarho_p                                         ", time(7)-time(6),"sec"
print*, "   > find correspondance between q and -q                         ", time(8)-time(7),"sec"
print*, "   > read ck                                                      ", time(9)-time(8),"sec"
print*, "   > rotation to molecular frame + OZ + rotation back to lab frame", time(10)-time(9),"sec"
print*, "   > check all q points have been used in OZ                      ", time(11)-time(10),"sec"
print*, "   > FFT-1 of deltarho_p                                          ", time(12)-time(11),"sec"
print*, "   > Fcproj: gather projections into gamma and sum                ", time(13)-time(12),"sec"











































































contains


    !
    ! Transform angles to projections
    !
    function angl2proj (foo_o) result (foo_p)
        implicit none
        real(dp), intent(in) :: foo_o(:) ! orientations from 1 to no
        complex(dp) :: foo_p(1:grid%np)
        integer :: itheta, iphi, ipsi, m, mup, mu, ip, p, o
        complex(dp), parameter :: zeroc=complex(0._dp,0._dp)
        complex(dp) :: xp
        real(dp) :: eightpisq_I_mrso
        eightpisq_I_mrso = 8._dp*acos(-1._dp)**2/real(mrso,dp)
        foo_theta_mu_mup = zeroc
        do itheta=1,ntheta
            do iphi=1,nphi
                do ipsi=1,npsi
                    fft2d_r2c%in(ipsi,iphi) = foo_o(grid%indo(itheta,iphi,ipsi))
                end do
            end do
            call dfftw_execute (fft2d_r2c%plan)
            foo_theta_mu_mup(itheta,0:mmax/mrso,0:mmax)   = CONJG( fft2d_r2c%out(:,1:mmax+1) )/real(nphi*npsi,dp)
            if (mmax>0) then
                foo_theta_mu_mup(itheta,0:mmax/mrso,-mmax:-1) = CONJG( fft2d_r2c%out(:,mmax+2:) )/real(nphi*npsi,dp)
            end if
        end do
        foo_p = zeroc
        do m=0,mmax
            do mup=-m,m
                do mu=0,m,mrso
                    ip = p3%p( m,mup,mu/mrso )
                    foo_p(ip)= sum(foo_theta_mu_mup(:,mu/mrso,mup)*p3%harm_sph(:,ip)*wtheta(:))*fm(m)
                end do
            end do
        end do



    ! do concurrent( m=0:mmax)
    !     do concurrent( mup=-m:m, mu=0:m:mrso)
    !         p = p3%p(m,mup,mu/mrso)
    !         foo_p(p) = zeroc
    !         do concurrent (itheta=1:ntheta, iphi=1:nphi, ipsi=1:npsi)
    !             o = grid%indo(itheta,iphi,ipsi)
    !             xp = exp(ii*mup*grid%phi(o) +ii*mu*grid%psi(o))
    !             foo_p(p) = foo_p(p) + fm(m)*p3%harm_sph(itheta,p)*grid%w(o)/sum(grid%w)*foo_o(o)*xp
    !         end do
    !     end do
    ! end do

end function angl2proj










    function proj2angl (foo_p) result (foo_o)
        implicit none
        real(dp) :: foo_o (1:grid%no)
        complex(dp), intent(in) :: foo_p(:) ! np
        integer :: o, p, m, mup, mu, ip, iphi, ipsi, itheta
        complex(dp), parameter :: ii=complex(0._dp,1._dp)
        complex(dp) :: foo_o_c
        complex(dp) :: R
        complex(dp) :: xp
        foo_theta_mu_mup = zeroc
        foo_o = 0._dp
        do itheta=1,ntheta
            do mup=-mmax,mmax
                do mu=0,mmax,mrso
                    do m= max(abs(mup),abs(mu)), mmax
                        ip=p3%p(m,mup,mu/mrso)
                        foo_theta_mu_mup(itheta,mu/mrso,mup) = foo_theta_mu_mup(itheta,mu/mrso,mup) &
                            + foo_p(ip)*p3%harm_sph(itheta,ip)
                        ! foo_theta_mu_mup(itheta,mu/mrso,mup) = foo_p(ip)*sum(p3%harm_sph(:,ip))*fm(m)
                            ! foo_theta_mu_mup(itheta,mu/mrso,mup) = foo_theta_mu_mup(itheta,mu/mrso,mup)&
                            ! +foo_p(ip)*p3%harm_sph(itheta,ip)*fm(m)
                    end do
                end do
            end do
        end do
        do itheta=1,ntheta
            fft2d_c2r%in(:,1:mmax+1) = 2*conjg( foo_theta_mu_mup(itheta,0:mmax/mrso,0:mmax)   )
            fft2d_c2r%in(:,mmax+2:)  = 2*conjg( foo_theta_mu_mup(itheta,0:mmax/mrso,-mmax:-1) )
            call dfftw_execute( fft2d_c2r%plan )
            do iphi=1,nphi
                do ipsi=1,npsi
                    o = grid%indo(itheta,iphi,ipsi)
                    foo_o(o) = fft2d_c2r%out(ipsi,iphi)
                end do
            end do
        end do


    do concurrent (itheta=1:ntheta, iphi=1:nphi, ipsi=1:npsi)
        o=grid%indo(itheta,iphi,ipsi)
        foo_o_c = zeroc
        do m=0,mmax
            do mup=-m,m
                do mu=0,m,2
                    if (mu==0) then
                        p=p3%p(m,mup,0)
                        xp = exp(-ii*mup*grid%phi(iphi))
                        foo_o_c = foo_o_c + foo_p(p) *fm(m) *xp *p3%harm_sph(itheta,p)
                    else
                        p=p3%p(m,mup,mu/mrso)
                        xp = exp(-ii*mup*grid%phi(iphi)-ii*mu*grid%psi(ipsi))
                        foo_o_c = foo_o_c + 2*foo_p(p) *fm(m)  *p3%harm_sph(itheta,p)
                    end if
                end do
            end do
        end do
        print*, foo_o(o), foo_o_c
        foo_o(o) = real(foo_o_c)
    end do
end function proj2angl



    subroutine for_all_q_find_indices_of_mq
        !
        ! for each vector q, defined by its components qx(ix_q), qy(iy_q), qz(iz_q)
        ! I am looking for the indices ix_mq, iy_mq, iz_mq for which
        ! qx(ix_mq), qy(iy_mq), qz(iz_mq) defines vector -q
        !
        implicit none
        integer :: i, ix_q, iy_q, iz_q, ix_mq, iy_mq, iz_mq
        if ( allocated(tableof_ix_mq) .and. allocated(tableof_iy_mq) .and. allocated(tableof_iz_mq)) then
            if ( size(tableof_ix_mq)==nx .and. size(tableof_iy_mq)==ny .and. size(tableof_iz_mq)==nz &
                .and. allocated(qx) .and. size(qx)==nx .and. allocated(qy) .and. size(qy)==ny .and.&
                    allocated(qz) .and. size(qz)==nz ) then
                return
            else
                print*, "problem detected in the subroutine that allocates and defines all q and mq related arrays"
                error stop
            end if
        end if

        allocate (qx(nx), qy(ny), qz(nz), source=0._dp)
        qx(1:nx) = grid%kx(1:nx)
        qy(1:ny) = grid%ky(1:ny)
        qz(1:nz) = grid%kz(1:nz)

        allocate (tableof_ix_mq(nx), tableof_iy_mq(ny), tableof_iz_mq(nz) ,source=-huge(1))
        i=0
        do ix_q=1,nx
            do ix_mq=1,nx
                if (qx(ix_q)==-qx(ix_mq) .or. qx(ix_q)+2._dp*pi*nx/lx==-qx(ix_mq) ) then
                    tableof_ix_mq(ix_q) = ix_mq
                    i=i+1
                    ! print*, "qx(ix_q), qx(ix_mq)", qx(ix_q), qx(ix_mq), "ix_q, iq_mq=", ix_q, ix_mq
                    cycle
                end if
            end do
        end do
        if (i/=nx) then
            print*, "I did not found all the -qx. Of",nx,", I found only",i
            error stop
        end if

        i=0
        do iy_q=1,ny
            do iy_mq=1,ny
                if ( qy(iy_q)==-qy(iy_mq) .or. qy(iy_q)+2._dp*pi*ny/ly==-qy(iy_mq)  ) then
                    tableof_iy_mq(iy_q) = iy_mq
                    i=i+1
                    ! print*, "qy(iy_q), qy(iy_mq)", qy(iy_q), qy(iy_mq), "iy_q, iq_mq=", iy_q, iy_mq
                    cycle
                end if
            end do
        end do
        if (i/=ny) then
            print*, "I did not found all the -qy. Of",ny,", I found only",i
            error stop
        end if

        i=0
        do iz_q=1,nz
            do iz_mq=1,nz
                if ( qz(iz_q)==-qz(iz_mq) .or. qz(iz_q)+2._dp*pi*nz/lz==-qz(iz_mq)  ) then
                    tableof_iz_mq(iz_q) = iz_mq
                    i=i+1
                    ! print*, "qz(iz_q), qz(iz_mq)", qz(iz_q), qz(iz_mq), "iz_q, iq_mq=", iz_q, iz_mq
                    cycle
                end if
            end do
        end do
        if (i/=nz) then
            print*, "I did not found all the -qz. Of",nz,", I found only",i
            error stop
        end if

        ! Here I want to check that we really have q == -q but for the point in each direction at index 0 and nx/2+1
        ! for odd numbers, please ask maximilien
        if (mod(nx,2)/=0 .or. mod(ny,2)/=0 .or. mod(nz,2)/=0) then
            print*, "nx, ny ou nz est impair"
            print*, "nx ny nz =", nx, ny, nz
            print*, "that is not compatible with our energy_cproj for now"
            error stop
        end if
        do ix_q=1,nx
            do iy_q=1,ny
                do iz_q=1,nz
                    ix_mq = tableof_ix_mq(ix_q)
                    iy_mq = tableof_iy_mq(iy_q)
                    iz_mq = tableof_iz_mq(iz_q)
                    if ( qx(ix_q)/=-qx(ix_mq) .or. qy(iy_q)/=-qy(iy_mq) .or. qz(iz_q)/=-qz(iz_mq) ) then
                        if (ix_q==1 .or. iy_q==1 .or. iz_q==1 .or. ix_q==nx/2+1 .or. iy_q==ny/2+1 .or. iz_q==nz/2+1) then
                            ! that is the expected behavior
                        else
                            print*,"q n'est pas l'opposé de mq dans la boucle de verification des q apres la recherche"
                            print*, "des indices de mq"
                            print*, "ix_q, iy_q, iz_q, q(ix_q, iy_q, iz_q)    =",ix_q, iy_q, iz_q, qx(ix_q), qy(iy_q), qz(iz_q)
                            print*, "ix_mq, iy_mq, iz_mq, q(ix_mq,iy_mq,iz_mq)=",ix_mq,iy_mq,iz_mq,qx(ix_mq),qy(iy_mq),qz(iz_mq)
                            error stop
                        end if
                    end if
                end do
            end do
        end do
        end subroutine for_all_q_find_indices_of_mq

        ! subroutine test_routines_calcul_de_Rm_mup_mu_q
        !     !
        !     ! ici test de la routine de Choi
        !     !
        !     implicit none
        !         real(dp) :: qx, qy, qz
        !         complex(dp) :: R(0:mmax,-mmax:mmax,-mmax:mmax), Rlu(0:mmax,-mmax:mmax,-mmax:mmax)
        !         complex(dp), parameter :: ii=complex(0._dp,1._dp)
        !         real(dp), parameter :: epsdp=epsilon(1._dp)
        !         integer :: ix_q, iy_q, iz_q, m, mup, mu
        !         do ix_q=1,nx
        !             do iy_q=1,ny
        !                 do iz_q=1,nz
        !                     qx=grid%kx(ix_q)
        !                     qy=grid%ky(iy_q)
        !                     qz=grid%kz(iz_q)
        !                     R= rotation_matrix_between_complex_spherical_harmonics_lu ([qx,qy,qz], mmax) ! exactement Luc
        !                     Rlu= rotation_matrix_between_complex_spherical_harmonics_lu ([qx,qy,qz], mmax) ! Lu
        !                     do m=0,mmax
        !                         do mup=-m,m
        !                             do mu=-m,m
        !                                 if (abs(&
        !                                     Rlu(m,mup,mu)-harm_sph(m,mup,mu,thetaofq(qx,qy,qz))*exp(-ii*mup*phiofq(qx,qy)))&
        !                                         >1E-10) then
        !                                     print*, "dans test_routines_calcul_de_Rm_mup_mu_q on a Choi+Lu /= Wigner"
        !                                     print*, "ix_q, iy_q, iz_q",ix_q, iy_q, iz_q
        !                                     print*, "qx qy qz", qx, qy, qz
        !                                     print*, "m, mup, mu", m, mup, mu
        !                                     print*, "Rlu(m,mup,mu)",Rlu(m,mup,mu)
        !                                     print*, "harm_sph(m,mup,mu,thetaofq(qx,qy,qz))*exp(-ii*mup*phiofq(qx,qy))",&
        !                                     harm_sph(m,mup,mu,thetaofq(qx,qy,qz))*exp(-ii*mup*phiofq(qx,qy))
        !                                     print*, "diff=",&
        !                                     abs(Rlu(m,mup,mu)-harm_sph(m,mup,mu,thetaofq(qx,qy,qz))*exp(-ii*mup*phiofq(qx,qy)))
        !                                     print*, "epsilon=",epsdp
        !                                     error stop "oihjkkhblkuyuiykjhbkjhg76"
        !                                 end if
        !                             end do
        !                         end do
        !                     end do
        !                 end do
        !             end do
        !         end do
        ! end subroutine test_routines_calcul_de_Rm_mup_mu_q

        subroutine read_ck_nonzero
            implicit none
            integer :: i, na, iq, m, n, mu, nu, khi, ia, nq
            character(3) :: somechar
            integer :: ufile, ios, mmax
            character(65) :: filename

            if (cq%isok) return
            nq=cq%nq
            if (.not.allocated(cq%normq)) allocate(cq%normq(nq), source=0._dp)
            mmax=grid%mmax

            select case (mmax)
            case (0)
                na = 1
                filename="input/direct_correlation_functions/water/SPCE/ck_nonzero_nmax0_ml"
                open(newunit=ufile, file=filename, iostat=ios, status="old", action="read")
                if ( ios /= 0 ) then
                    print*, "Cant open file", filename
                    error stop "in module_energy_cproj.f90"
                end if
            case (1)
                na = 6
                filename="input/direct_correlation_functions/water/SPCE/ck_nonzero_nmax1_ml"
                open(newunit=ufile, file=filename, iostat=ios, status="old", action="read")
                if ( ios /= 0 ) then
                    print*, "Cant open file", filename
                    error stop "in module_energy_cproj.f90"
                end if
            case (2)
                na = 75
                filename="input/direct_correlation_functions/water/SPCE/ck_nonzero_nmax2_ml"
                open(newunit=ufile, file=filename, iostat=ios, status="old", action="read")
                if ( ios /= 0 ) then
                    print*, "Cant open file", filename
                    error stop "in module_energy_cproj.f90"
                end if
            case (3)
                na = 252
                filename="input/direct_correlation_functions/water/SPCE/ck_nonzero_nmax3_ml"
                open(newunit=ufile, file=filename, iostat=ios, status="old", action="read")
                if ( ios /= 0 ) then
                    print*, "Cant open file", filename
                    error stop "in module_energy_cproj.f90"
                end if
            case (4)
                na = 877
                filename="input/direct_correlation_functions/water/SPCE/ck_nonzero_nmax4_ml"
                open(newunit=ufile, file=filename, iostat=ios, status="old", action="read")
                if ( ios /= 0 ) then
                    print*, "Cant open file", filename
                    error stop "in module_energy_cproj.f90"
                end if
            case (5)
                na = 2002
                filename="input/direct_correlation_functions/water/SPCE/ck_nonzero_nmax5_ml"
                open(newunit=ufile, file=filename, iostat=ios, status="old", action="read")
                if ( ios /= 0 ) then
                    print*, "Cant open file", filename
                    error stop "in module_energy_cproj.f90"
                end if
            case default
                print*, "In module_energy_cproj, you want mmax>5. I cant read the corresponding file."
                error stop
            end select


            ! i=(1+mmax)*(1+2*mmax)*(3+2*mmax)*(5+4*mmax*(2+mmax))/15
          !  if (na/=i) then
          !      print*, "in read_ck_cproj nalpha est bizarre"
          !      print*, "na=",na
          !      print*, "it should be (1+mmax)*(1+2*mmax)*(3+2*mmax)*(5+4*mmax*(2+mmax))/15 =",i
          !      ! error stop
          !  end if

            ! i=sum([([([([([( 1 ,nu=-n,n)] ,mu=-m,m)], n=abs(khi),mmax)] ,m=abs(khi),mmax)] ,khi=-mmax,mmax)]  )
          !  if (na/=i) then
          !      print*, "in read ck na /= nalpha"
          !      print*, "na=",na
          !      print*, "bruteforce=",i
          !      ! error stop
          !  end if

            if (allocated(ck) .and. .not.cq%isok) then
                print*, "ck is already allocated but .not. cq%isok"
                error stop
            end if

            if (.not.allocated( cq%m)) then
                allocate ( cq%m(na) ,source=-huge(1))
                allocate ( cq%n(na) ,source=-huge(1))
                allocate ( cq%mu(na) ,source=-huge(1))
                allocate ( cq%nu(na) ,source=-huge(1))
                allocate ( cq%khi(na) ,source=-huge(1))
                allocate ( cq%a(0:mmax,0:mmax,-mmax:mmax,-mmax:mmax,-mmax:mmax), source=-huge(1)) ! m n mu nu khi. -huge is used to spot more easily bugs that may come after
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


            do ia=1,na
                m = cq%m(ia)
                n = cq%n(ia)
                mu = cq%mu(ia)
                nu = cq%nu(ia)
                khi = cq%khi(ia)
                cq%a(m,n,mu,nu,khi) = ia
            end do

            allocate( ck(na,nq), source=zeroc)
            do iq=1,nq
                read(ufile,*) cq%normq(iq), ck(:,iq)
            end do
            close(ufile)
            cq%dq = cq%normq(2)

            deallocate (cq%m, cq%n, cq%mu, cq%nu, cq%khi)
            cq%isok=.true.
            cq%na=na
            cq%nq=nq
        end subroutine read_ck_nonzero


    end subroutine energy_cproj_mrso

end module module_energy_cproj_mrso

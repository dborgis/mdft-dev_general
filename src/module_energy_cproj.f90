module module_energy_cproj

    use iso_c_binding
    use precision_kinds, only: dp
    use module_grid, only: grid
    use module_solvent, only: solvent
    use module_rotation

    implicit none
    private

    !
    ! FFTW3 header - modern (f03) version. Expects iso_c_binding
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
        integer :: nq=200 ! number of q points in c(q)
        integer, allocatable :: a(:,:,:,:,:) ! index of m, n, mu, nu, khi => (0:mmax,0:mmax,-mmax:mmax,-mmax:mmax,-mmax:mmax) ! m n mu nu khi
        real(dp) :: dq
        real(dp), allocatable :: normq(:)
    end type
    type (cq_type) :: cq
    complex(dp), allocatable :: ck(:,:)

    real(dp), allocatable :: fm(:)

    type :: fft2d_c_type
        type(c_ptr) :: plan
        complex(c_double), allocatable :: in(:,:), out(:,:)
    end type fft2d_c_type
    type (fft2d_c_type) :: fft2d_c

    type :: fft3d_c_type
        type(c_ptr) :: plan_forward
        type(c_ptr) :: plan_backward
        logical :: plan_forward_ok = .false.
        logical :: plan_backward_ok = .false.
    end type fft3d_c_type
    type (fft3d_c_type) :: fft3d

    integer, allocatable :: tableof_ix_mq(:), tableof_iy_mq(:), tableof_iz_mq(:)
    real(dp), allocatable :: qx(:), qy(:), qz(:) ! vector q and its components tabulated

    complex(dp), allocatable :: deltarho_p(:,:,:,:)
    complex(dp), allocatable :: deltarho_p_q(:)
    complex(dp), allocatable :: deltarho_p_mq(:)
    complex(dp), allocatable :: gamma_p(:,:,:,:)

    type :: p3_type
        real(dp), allocatable :: harm_sph(:,:) ! tabulation des harmoniques sphériques r(m,mup,mu,theta) en un tableau r(itheta,p)
        integer, allocatable :: p(:,:,:) ! index of the projection corresponding to m, mup, mu
        integer, allocatable :: m(:) ! m for projection 1 to np
        integer, allocatable :: mup(:) ! mup for projection 1 to np. mup corresponds to phi
        integer, allocatable :: mu(:) ! mu for projection 1 to np. mu corresponds to psi
    end type p3_type
    type (p3_type) :: p3


    type :: fft2d_type
        type(c_ptr) :: plan
        real(dp), allocatable    :: in(:,:)
        complex(dp), allocatable :: out(:,:)
        logical :: isalreadyplanned
    end type
    type :: ifft2d_type
        type(c_ptr) :: plan
        complex(dp), allocatable :: in(:,:)
        real(dp), allocatable    :: out(:,:)
        logical :: isalreadyplanned
    end type
    type(  fft2d_type ), save :: fft2d
    type( ifft2d_type ), save :: ifft2d


    public :: energy_cproj

contains

    subroutine energy_cproj (ff,df)
        ! use ieee_arithmetic
        use iso_c_binding, only: c_ptr, dp=>c_double
        use module_grid, only: grid
        use module_thermo, only: thermo
        implicit none
        real(dp), intent(out) :: ff
        real(dp), intent(out) :: df(:,:,:,:,:)
        real(dp) :: dv, kT
        logical :: q_eq_mq
        integer :: ix, iy, iz, ix_q, iy_q, iz_q, ix_mq, iy_mq, iz_mq, i, p
        integer :: nx, ny, nz, np, no, ns, ntheta, nphi, npsi, mmax, na, nq, molrotsymorder
        integer :: m, n, mu, nu, khi, mup, ia, ip, ipsi, iphi, itheta, iq
        complex(dp), allocatable :: gamma_p_q(:), gamma_p_mq(:)
        complex(dp) :: gamma_m_khi_mu_q, gamma_m_khi_mu_mq
        real(dp) :: q(3), mq(3)
        real(dp) :: theta(grid%ntheta), wtheta(grid%ntheta)
        real(dp) :: lx, ly, lz, rho0
        logical, allocatable :: gamma_p_isok(:,:,:)
        integer :: points_q_considere_en_vrai


        if (.not.allocated(fft2d%in)) allocate (fft2d%in(grid%npsi,grid%nphi))
        if (.not.allocated(fft2d%out)) allocate (fft2d%out(grid%npsi/2+1,grid%nphi))
        if (.not.allocated(ifft2d%in)) allocate (ifft2d%in(grid%npsi/2+1,grid%nphi))
        if (.not.allocated(ifft2d%out)) allocate (ifft2d%out(grid%npsi,grid%nphi))

        lx=grid%lx
        ly=grid%ly
        lz=grid%lz
        mmax=grid%mmax
        molrotsymorder=grid%molrotsymorder
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
        ntheta=grid%ntheta
        nphi=grid%nphi
        npsi=grid%npsi
        rho0 = solvent(1)%rho0

        if (.not.allocated(fm)) allocate (fm(0:mmax) ,source= [( sqrt(real(2*m+1,dp))/real(nphi*npsi,dp) ,m=0,mmax  )])


        if (.not. allocated (gamma_p_q) ) allocate (gamma_p_q(np), source=zeroc)
        if (.not. allocated (gamma_p_mq) ) allocate (gamma_p_mq(np), source=zeroc)
        if (.not. allocated (gamma_p) ) allocate (gamma_p(np,nx,ny,nz) ,source=zeroc)
        if (.not. allocated (deltarho_p) ) allocate (deltarho_p(np,nx,ny,nz) ,source=zeroc)
        if (.not. allocated (deltarho_p_q) ) allocate (deltarho_p_q(np) ,source=zeroc)
        if (.not. allocated (deltarho_p_mq) ) allocate (deltarho_p_mq(np) ,source=zeroc)
        if (.not. allocated (gamma_p_isok) ) allocate (gamma_p_isok(nx,ny,nz), source=.false.)

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
            allocate ( p3%p(0:mmax,-mmax:mmax, 0:mmax) ,source=-huge(1))
            allocate ( p3%m(np) ,source=-huge(1))
            allocate ( p3%mup(np) ,source=-huge(1))
            allocate ( p3%mu(np) ,source=-huge(1))
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
                    do mu=0,m,molrotsymorder
                        ip=ip+1
                        IF (ip > np) ERROR STOP "p > np at line 166"
                        p3%p(m,mup,mu) = ip
                        p3%m(ip) = m
                        p3%mup(ip) = mup
                        p3%mu(ip) = mu
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
        call test_routines_calcul_de_Rm_mup_mu_q

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


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! PROJECTION DE Δρ(r,ω) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !
        ! On passe tout de suite en projection, dans le repère cartesien
        ! ces projections sont complexes, mais la symetrie hermitienne
        ! permet de ne garder que les mup>=0 ou les mu>=0.
        ! On choisit les mu>=0 (cf doc de Luc)
        !
        if (.not. fft2d%isalreadyplanned) then
            call dfftw_plan_dft_r2c_2d(  fft2d%plan, npsi, nphi,  fft2d%in,  fft2d%out, FFTW_EXHAUSTIVE ) ! npsi est en premier indice
            allocate (fft2d_c%in(npsi,nphi))
            allocate (fft2d_c%out(npsi,nphi))
            call dfftw_plan_dft_2d (fft2d_c%plan, npsi, nphi, fft2d_c%in, fft2d_c%out, FFTW_BACKWARD, FFTW_EXHAUSTIVE)
            fft2d%isalreadyplanned =.true.
        end if
        if (.not. ifft2d%isalreadyplanned) then
            call dfftw_plan_dft_c2r_2d(  ifft2d%plan, npsi, nphi, ifft2d%in, ifft2d%out, FFTW_EXHAUSTIVE )
            ifft2d%isalreadyplanned =.true.
        end if
        if (.not. fft3d%plan_backward_ok) then
            call dfftw_plan_dft_3d (fft3d%plan_backward,&
                                    nx, ny, nz, deltarho_p(1,:,:,:), deltarho_p(1,:,:,:), FFTW_BACKWARD, FFTW_PATIENT) ! TODO CHECK ESTIMATE VS REST & IS IT WORTH CHANGEING THE PLAN FLAG FOR DIFFERENT np ? Certainly!
            fft3d%plan_backward_ok = .true.
        end if
        if (.not. fft3d%plan_forward_ok) then
            call dfftw_plan_dft_3d( fft3d%plan_forward,&
                                    nx, ny, nz, gamma_p(1,1:nx,1:ny,1:nz), gamma_p(1,1:nx,1:ny,1:nz), FFTW_FORWARD, FFTW_PATIENT )
            fft3d%plan_forward_ok = .true.
        end if



        !
        ! Projection of delta rho
        !
        do iz=1,nz
            do iy=1,ny
                do ix=1,nx
                    deltarho_p(1:np,ix,iy,iz) = euler2proj( solvent(1)%density(ix,iy,iz,1:no)-solvent(1)%rho0 )
                    !
                    ! Test if result contains NaN or Inf
                    !
                    if ( any(deltarho_p(1:np,ix,iy,iz) /= deltarho_p(1:np,ix,iy,iz)) ) then
                        print*, "Error in euler2proj for ix,iy,iz =", ix, iy, iz
                        print*, "deltarho_p(1:np,ix,iy,iz) = ",deltarho_p(1:np,ix,iy,iz)
                        error stop "Bug found in euler2proj. NaN or Inf is hiding somewhere."
                    end if
                end do
            end do
        end do



        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! FFT DE Δρ_p(r) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !
        ! On a les projections sur la grille cartesienne
        ! On veut passer dans l'espace de Fourier pour calculer les convolutions spatiales
        ! On fait donc une FFT 3D pour chacune des projections.
        ! Les projections sont complexes, il s'agit donc d'une FFT3D C2C habituelle : Il n'y a pas de symétrie hermitienne.
        !
        ! call dfftw_plan_dft_3d( fft3d%plan, nx, ny, nz, fft3d%in, fft3d%in, FFTW_BACKWARD, FFTW_ESTIMATE ) ! in-place. BACKWARD means exp(iqr)
        ! do ip=1,np
        !     fft3d%in = deltarho_p(ip,1:nx,1:ny,1:nz)
        !     call dfftw_execute( fft3d%plan )
        !     deltarho_p(ip,1:nx,1:ny,1:nz) = fft3d%in
        ! end do

        do ip=1,np
            call dfftw_execute_dft( fft3d%plan_backward, deltarho_p(ip,:,:,:), deltarho_p(ip,:,:,:) )
        end do

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

        !
        ! Read Luc's direct correlation function c^{m,n}_{mu,nu_,chi}(|q|)
        ! projected on generalized spherical harmonics
        ! in the intermolecular frame
        ! normq is norm of q, |q|, that correspond to the index iq in ck(ia,iq)
        !
        !call read_ck_nmax (ck, normq)
        !call read_ck_toutes_nmax( ck, normq)
        if (.not.cq%isok) call read_ck_nonzero

        !
        ! For all vectors q and -q handled simultaneously
        !
        points_q_considere_en_vrai = 0

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

                    if (gamma_p_isok(ix_q,iy_q,iz_q).eqv..true.) then
                        if (gamma_p_isok(ix_mq, iy_mq,iz_mq).eqv..true.) then
                            cycle
                        else
                            print*,"in the main loop of energy_cproj, i have already done q but not mq!"
                            print*, "ix_q, iy_q, iz_q =",ix_q, iy_q, iz_q
                            print*, "ix_mq, iy_mq, iz_mq =",ix_mq, iy_mq, iz_mq
                            error stop
                        end if
                    end if

                    points_q_considere_en_vrai = points_q_considere_en_vrai +1

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

                    call rotate_to_molecular_frame ! does q and mq at the same time.


                    if(any(deltarho_p_mq/=deltarho_p_mq)) then
                        do ip=1,np
                            print*, ip, ix_q, iy_q, iz_q, deltarho_p_mq(ip)
                        end do
                        error stop "yhrgvussnehfkil"
                    end if
                    if(any(deltarho_p_q/=deltarho_p_q)) then
                        do iq=1,np
                            print*, ip, ix_q, iy_q, iz_q, deltarho_p_q(ip)
                        end do
                        error stop "opiloyukiuskerg"
                    end if


                    ! if (q_eq_mq .and. any(deltarho_p_q/=deltarho_p_mq) ) then
                    !     print*, "q = mq mais apres la rotation vers le molecular frame, on n'a plus deltarho_p_q=deltarho_p_mq"
                    !     print*, "ix_q,iy_q,iz_q=", ix_q,iy_q,iz_q
                    !     print*, "ix_mq,iy_mq,iz_mq=", ix_mq, iy_mq, iz_mq
                    !     print*, "deltarho_p_q", deltarho_p_q
                    !     print*, "deltarho_p_mq", deltarho_p_mq
                    !     error stop "oijazd"
                    ! end if

                    !
                    ! c^{m,n}_{mu,nu,chi}(|q|) is tabulated for 200 values of |q|.
                    ! Find the tabulated value that is closest to |q|. Its index is iq.
                    ! Note |q| = |-q| so iq is the same for both vectors.
                    !
                    iq = int(norm2(q)/cq%dq)+1
                    if (iq>nq) then
                        iq=nq
                    else if (iq<=0) then
                        print*, "in energy_cproj, ik <=0"
                        error stop
                    end if


                    !
                    ! Ornstein-Zernike in the molecular frame
                    ! We do OZ for q and -q at the same time
                    !
                    gamma_p_q(1:np) = zeroc
                    gamma_p_mq(1:np) = zeroc

                    do khi=-mmax,mmax
                        do m=abs(khi),mmax
                            do mu=0,m,molrotsymorder ! not -m,m a cause des symetries ! EST CE QU'IL Y A UNE RAISON POUR QUE LES mu IMPAIRES SOIENT NON NULS ICI ?
                                if (mod(mu,2) /= 0) cycle

                                gamma_m_khi_mu_q= zeroc
                                gamma_m_khi_mu_mq = zeroc

                                do n=abs(khi),mmax
                                    do nu=-n,n,molrotsymorder
                                        if (mod(nu,2) /=0) cycle ! don't threat cases where c is (0,0) for instance all nu even in water (molrotsymorder==2)

                                        ia = cq%a(m,n,mu,nu,khi) ! the index of the projection of c(q). 1<=ia<na
                                        if (ia<=0) then
                                            print*,"ia=",ia
                                            print*,"for m, n, mu, nu, khi=",m,n,mu,nu,khi
                                            error stop "ia<=0"
                                        else if (ia>cq%na) then
                                            print*,"ia=",ia
                                            print*,"for m, n, mu, nu, khi=",m,n,mu,nu,khi
                                            error stop "ia>na"
                                        end if

                                        ip = p3%p(n,khi,abs(nu))

                                        if (nu<0) then ! no problem with delta rho (n, khi, -nu) since -nu>0. Thus, we apply eq. 1.30 directly
                                            gamma_m_khi_mu_q  = gamma_m_khi_mu_q  + (-1)**(khi+nu) *ck(ia,iq) *deltarho_p_q(ip)
                                            gamma_m_khi_mu_mq = gamma_m_khi_mu_mq + (-1)**(khi+nu) *ck(ia,iq) *deltarho_p_mq(ip)
                                        else ! transform delta rho (n, khi, -nu)(q) into conjg( deltarho(n,khi,nu)(-q) )
                                            select case (mmax)
                                            case (0)
                                            gamma_m_khi_mu_q= gamma_m_khi_mu_q  + (-1)**(n) *ck(ia,iq) *conjg(deltarho_p_mq(ip))
                                            gamma_m_khi_mu_mq= gamma_m_khi_mu_mq + (-1)**(n) *ck(ia,iq) *conjg(deltarho_p_q(ip))
                                            case default
                                                if (ip/=p3%p(1,0,0) .or. ix_q/=1 .or. iy_q/=1 .or. iz_q/=1) then
                                                gamma_m_khi_mu_q= gamma_m_khi_mu_q  + (-1)**(n) *ck(ia,iq) *conjg(deltarho_p_mq(ip))
                                                gamma_m_khi_mu_mq= gamma_m_khi_mu_mq + (-1)**(n) *ck(ia,iq) *conjg(deltarho_p_q(ip))
                                                end if
                                            end select
                                        end if

                                    end do
                                end do

                                if (gamma_m_khi_mu_q/=gamma_m_khi_mu_q) error stop "jhwdlijweflkhwke"
                                if (gamma_m_khi_mu_mq/=gamma_m_khi_mu_mq) error stop "lijsfdkuglserhk"

                                ip=p3%p(m,khi,mu)
                                gamma_p_q(ip)  = gamma_m_khi_mu_q
                                gamma_p_mq(ip) = gamma_m_khi_mu_mq

                                ! if (q_eq_mq .and. any(gamma_p_q/=gamma_p_mq) ) then
                                !     print*, "we have q=mq"
                                !     print*, "for ix_q,iy_q,iz_q=",  ix_q, iy_q, iz_q
                                !     print*, "and ix_mq,iy_mq,iz_mq=", ix_mq, iy_mq, iz_mq
                                !     print*, "q=",q
                                !     print*, "mq=",mq
                                !     print*, "gamma_p_q=",gamma_p_q
                                !     print*, "gamma_p_mq=",gamma_p_mq
                                !     print*, "deltarho_p_q=",deltarho_p_q
                                !     print*, "deltarho_p_mq=",deltarho_p_mq
                                !     error stop "apzokdukhsf"
                                ! end if

                            end do
                        end do
                    end do

                    if (any(gamma_p_q/=gamma_p_q)) error stop "gammapqijosuiheofij"
                    if (any(gamma_p_mq/=gamma_p_mq)) error stop "ljsfkjhfksgrlkhh"

                    !
                    ! Rotation from molecular frame to fix laboratory (Fourier) frame
                    !
                    call rotate_to_fixed_frame ! does q and mq at the same time.

                    if (any(gamma_p_q/=gamma_p_q)) error stop "oijjiuurkgserg"
                    if (any(gamma_p_mq/=gamma_p_mq)) error stop "098123KJN"


                    ! TOdo JE NE SAIS PAS JUSTIFIER POURQUOI FAUT RAJOTUER UN CONJG SI (VOIR IF CI DESSOUS) OMG OMG OMG OMG OMG

                    if (gamma_p_isok(ix_q,iy_q,iz_q).eqv..true. .and. any( abs(gamma_p(:,ix_q,iy_q,iz_q)-gamma_p_q(:))>epsdp )) then
                        print*, "gamma(q) already been calculated for ix_q, iy_q, iz_q =",ix_q,iy_q,iz_q
                        print*, "q=",q
                        print*,'mq=',mq
                        print*, "                        ip          m          mup          mu"
                        do ip=1,np
                            print*, "old gamma_p =",ip, p3%m(ip), p3%mup(ip), p3%mu(ip), gamma_p(ip,ix_q,iy_q,iz_q)
                            print*, "new gamma_p =",ip, p3%m(ip), p3%mup(ip), p3%mu(ip), gamma_p_q(ip)
                            print*,
                        end do
                        error stop "jfhzoenfhsozk"
                    end if

        if (gamma_p_isok(ix_mq,iy_mq,iz_mq).eqv..true.  .and. any( abs(gamma_p(:,ix_mq,iy_mq,iz_mq)-gamma_p_mq(:))>epsdp)) then
                        print*, "gamma(q) already been calculated for ix_mq, iy_mq, iz_mq =",ix_mq,iy_mq,iz_mq
                        print*, "mq=",mq
                        print*,"q=",q
                        print*, "                        ip          m          mup          mu"
                        do ip=1,np
                            print*, "old gamma_mq =",ip, p3%m(ip), p3%mup(ip), p3%mu(ip), gamma_p(ip,ix_mq,iy_mq,iz_mq)
                            print*, "new gamma_mq =",ip, p3%m(ip), p3%mup(ip), p3%mu(ip), gamma_p_mq(ip)
                            print*,
                        end do
                        error stop "oihskdjhfkuhse"
                    end if

                    gamma_p(1:np, ix_q,  iy_q,  iz_q) = gamma_p_q(1:np)
                    if( q_eq_mq .and. (ix_q==nx/2+1.or.iy_q==ny/2+1.or.iz_q==nz/2+1)) then
                        gamma_p(1:np, ix_mq, iy_mq, iz_mq) = conjg(gamma_p_mq(1:np))
                    else
                        gamma_p(1:np, ix_mq, iy_mq, iz_mq) = gamma_p_mq(1:np)
                    end if

                    ! we gamma_p_isok that all gamma(:,ix, iy, iz) have been calculated
                    gamma_p_isok(ix_q,iy_q,iz_q)=.true.
                    gamma_p_isok(ix_mq,iy_mq,iz_mq)=.true.

                    if (any(gamma_p(:,ix_q,iy_q,iz_q)/=gamma_p(:,ix_q,iy_q,iz_q))) then
                        print*, "I FOUND A NAN OR SOMETHING NASTY"
                        print*, ix_q, iy_q, iz_q
                        do ip=1,np
                            print*, ip, deltarho_p(ip,ix_q,iy_q,iz_q), gamma_p(ip,ix_q,iy_q,iz_q)
                        end do
                        error stop "ihsdfuygfkhll"
                    end if

                end do
            end do
        end do


        if (.not.all(gamma_p_isok.eqv..true.)) then
            print*, "not all gamma_p(projections,ix,iy,iz) have not been computed"
            error stop
        end if


        !
        ! FFT3D from Fourier space to real space
        !
        !call dfftw_plan_dft_3d(  fft3d%plan, nx, ny, nz, fft3d%in, fft3d%in, FFTW_FORWARD, FFTW_MEASURE )
        ! do ip=1,np
        !     fft3d%in = gamma_p(ip,1:nx,1:ny,1:nz)
        !     call dfftw_execute( fft3d%plan )
        !     gamma_p(ip,:,:,:) = fft3d%in/real(nx*ny*nz,dp)
        ! end do
        do ip=1,np
            call dfftw_execute_dft( fft3d%plan_forward, gamma_p(ip,1:nx,1:ny,1:nz), gamma_p(ip,1:nx,1:ny,1:nz) )
        end do
        gamma_p=gamma_p/real(nx*ny*nz,dp)

        !
        ! Gather projections
        !
        call dfftw_plan_dft_c2r_2d( ifft2d%plan, npsi, nphi, ifft2d%in, ifft2d%out, FFTW_EXHAUSTIVE)

        df=0._dp
        ff=0._dp
        do iz=1,nz
            do iy=1,ny
                do ix=1,nx
                    ! gamma_o(1:no) = proj2euler( gamma_p(1:np,ix,iy,iz) )
                    df(ix,iy,iz,:,1) = -kT*dv*proj2euler( gamma_p(1:np,ix,iy,iz) )/0.0333_dp*grid%w(:)**2*real(nphi*npsi*ntheta,dp)
                    ff=ff +sum( df(ix,iy,iz,:,1)*(solvent(1)%density(ix,iy,iz,:)-rho0) ) /2._dp
                end do
            end do
        end do

contains


        subroutine rotate_to_fixed_frame
            implicit none
            complex(dp) :: gamma_p_q_temp(1:np)
            complex(dp) :: gamma_p_mq_temp(1:np)
            complex(dp) :: gshrot (0:grid%mmax, -grid%mmax:grid%mmax, -grid%mmax:grid%mmax)
            ! complex(dp) :: gshrot_mq (0:grid%mmax, -grid%mmax:grid%mmax, -grid%mmax:grid%mmax)
            integer :: m, khi, mu, mup, ip, ip2
            if (ix_q==1 .and. iy_q==1 .and. iz_q==1) then ! that is if |q|=0. Rotation is arbitrary: we chose the Identity and thus have nothing to do.
                !return
            else
                gamma_p_q_temp = zeroc
                gamma_p_mq_temp = zeroc
                gshrot = conjg(rotation_matrix_between_complex_spherical_harmonics_lu(q,mmax) )
                ! gshrot_mq = conjg(rotation_matrix_between_complex_spherical_harmonics_lu(mq,mmax) )
                do m=0,mmax
                    do mup=-m,m
                        do mu=0,m
                            ip=p3%p(m,mup,mu)
                            ! Equation 1.22
                            do khi=-m,m
                                ip2=p3%p(m,khi,mu)
                                gamma_p_q_temp(ip) = gamma_p_q_temp(ip)   + gamma_p_q(ip2)*gshrot(m,mup,khi)
                                gamma_p_mq_temp(ip) = gamma_p_mq_temp(ip) + gamma_p_mq(ip2)*(-1)**m*gshrot(m,mup,-khi)
                                ! gamma_p_mq_temp(ip) = gamma_p_mq_temp(ip) + gamma_p_mq(ip2)*gshrot_mq(m,mup,khi)

                                ! if (q_eq_mq .and. gamma_p_q_temp(ip)/=gamma_p_mq_temp(ip)) then
                                !     print*, "q == mq mais pourtant deltarho_p_q_temp /= deltarho_p_mq_temp"
                                !     print*, "m,khi,mu,ip,mup,ip2,gamma_p_q_temp(ip),gamma_p_mq_temp(ip)"
                                !     print*,m,khi,mu,ip,mup,ip2,gamma_p_q_temp(ip),gamma_p_mq_temp(ip)
                                !     error stop
                                ! end if

                            end do
                            !
                        end do
                    end do
                end do
                gamma_p_q = gamma_p_q_temp
                gamma_p_mq = gamma_p_mq_temp
            end if
        end subroutine rotate_to_fixed_frame

        subroutine rotate_to_molecular_frame
            implicit none
            complex(dp) :: deltarho_p_q_temp(1:np)
            complex(dp) :: deltarho_p_mq_temp(1:np)
            complex(dp) :: gshrot (0:grid%mmax, -grid%mmax:grid%mmax, -grid%mmax:grid%mmax)
            ! complex(dp) :: gshrot_mq (0:grid%mmax, -grid%mmax:grid%mmax, -grid%mmax:grid%mmax)
            integer :: m, khi, mu, mup, ip2, ip
            ! complex(dp) :: a, b, c, d, e
            !
            ! Rotate projections from laboratory frame to molecular frame
            ! If q=0 then chose the rotation you want: the identity is the simplest option.
            !
            if (ix_q==1 .and. iy_q==1 .and. iz_q==1) then ! that is if |q|=0. Rotation is arbitrary: we chose the Identity and thus have nothing to do.
                ! return
            else
                gshrot    = rotation_matrix_between_complex_spherical_harmonics_lu ( q, mmax)
                ! gshrot_mq = rotation_matrix_between_complex_spherical_harmonics_lu ( mq, mmax)

                deltarho_p_q_temp = zeroc
                deltarho_p_mq_temp = zeroc
                do m=0,mmax
                    do khi=-m,m
                        do mu=0,m
                            ip=p3%p(m,khi,mu)
                            ! Equation 1.22
                            do mup=-m,m
                                ip2=p3%p(m,mup,mu)
                                deltarho_p_q_temp(ip) = deltarho_p_q_temp(ip)  + deltarho_p_q(ip2)*gshrot(m,mup,khi)
                                ! Eq. 1.23 We don't need to compute gshrot for -q. We do q and -q at the same time.
                                deltarho_p_mq_temp(ip) = deltarho_p_mq_temp(ip) + deltarho_p_mq(ip2)*(-1)**m*gshrot(m,mup,-khi)
                        !!        ! deltarho_p_mq_temp(ip) = deltarho_p_mq_temp(ip) + deltarho_p_mq(ip2)*gshrot_mq(m,mup,khi)
                                ! if (q_eq_mq .and. deltarho_p_q_temp(ip)/=deltarho_p_mq_temp(ip)) then
                                !     print*, "q == mq mais pourtant deltarho_p_q_temp /= deltarho_p_mq_temp"
                                !     print*,m,khi,mu,ip,mup,ip2,deltarho_p_q_temp(ip),deltarho_p_mq_temp(ip)
                                !     error stop
                                ! end if

                                ! a = (-1)**m*gshrot(m,mup,-khi)
                                ! b = (-1)**(m+mup+khi)*conjg(gshrot(m,-mup,khi))
                                ! c = gshrot_mq(m,mup,khi)
                                ! d = a-b
                                ! e = b-c

                                ! if (abs(e)>epsdp .or. abs(d)>epsdp) then
                                !     print*, "ix, iy, iz =", ix_q, iy_q, iz_q,"m,mup,mu=", m, mup, khi
                                !     print*, "q=", q
                                !     print*, "mq=", mq
                                !     print*, "gshrot(q,m,mup,khi)                 =", gshrot(m,mup,khi)
                                !     print*, "gshrot(mq,m,mup,khi) par sym 1.23/1 =", (-1)**m*gshrot(m,mup,-khi)
                                !     print*, "gshrot(mq,m,mup,khi) par sym 1.23/2 =", (-1)**(m+mup+khi)*conjg(gshrot(m,-mup,khi))
                                !     print*, "gshrot(mq,m,mup,khi) brutal         =", gshrot_mq(m,mup,khi)
                                !     error stop "85643820"
                                ! end if

                            end do
                        end do
                    end do
                end do
                deltarho_p_q = deltarho_p_q_temp
                deltarho_p_mq = deltarho_p_mq_temp
            end if
        end subroutine rotate_to_molecular_frame


        subroutine read_ck_nonzero
            implicit none
            integer :: i, na, iq, m, n, mu, nu, khi, ia, nq
            character(3) :: somechar
            integer :: ufile, ios, mmax
            integer, allocatable, dimension(:) :: m_ck, n_ck, mu_ck, nu_ck, khi_ck
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

            if (.not.allocated( m_ck)) allocate( m_ck(na) )
            if (.not.allocated( n_ck)) allocate( n_ck(na) )
            if (.not.allocated( mu_ck)) allocate( mu_ck(na) )
            if (.not.allocated( nu_ck)) allocate( nu_ck(na) )
            if (.not.allocated( khi_ck)) allocate( khi_ck(na) )
            !
            ! Skip 10 lines of comments
            !
            do i=1,10
                read(ufile,*)
            end do

            read(ufile,*) somechar
            read(ufile,*) somechar, m_ck
            read(ufile,*) somechar, n_ck
            read(ufile,*) somechar, mu_ck
            read(ufile,*) somechar, nu_ck
            read(ufile,*) somechar, khi_ck
            read(ufile,*)

            allocate (cq%a(0:mmax,0:mmax,-mmax:mmax,-mmax:mmax,-mmax:mmax), source=-huge(1)) ! m n mu nu khi. -huge is used to spot more easily bugs that may come after
            do ia=1,na
                m = m_ck(ia)
                n = n_ck(ia)
                mu = mu_ck(ia)
                nu = nu_ck(ia)
                khi = khi_ck(ia)
                cq%a(m,n,mu,nu,khi) = ia
            end do

            allocate( ck(na,nq), source=zeroc)
            do iq=1,nq
                read(ufile,*) cq%normq(iq), ck(:,iq)
            end do
            close(ufile)
            cq%dq = cq%normq(2)

            deallocate (m_ck, n_ck, mu_ck, nu_ck, khi_ck)
            cq%isok=.true.
            cq%na=na
            cq%nq=nq
        end subroutine read_ck_nonzero

        FUNCTION euler2proj (foo_o) RESULT (foo_p)
            use module_rotation, only: harm_sph
            IMPLICIT NONE
            real(dp), intent(in) :: foo_o(:) ! orientations from 1 to no
            COMPLEX(dp) :: foo_theta_mup_mu(1:grid%ntheta,0:grid%mmax/grid%molrotsymorder,-grid%mmax:grid%mmax) ! itheta,mu,mup    note we changed mu(psi) to the 2nd position
            COMPLEX(dp) :: foo_p(1:grid%np)
            INTEGER :: itheta, iphi, ipsi, m, mup, mu, ip, mmax, molrotsymorder
            complex(dp), allocatable :: test_explicit(:,:,:), foo_theta_mup_mu_full(:,:,:)
            complex(dp), allocatable :: proj_m_mup_mu(:,:,:)
            mmax=grid%mmax
            molrotsymorder=grid%molrotsymorder
            foo_p = zeroc
            foo_theta_mup_mu = zeroc


            allocate (foo_theta_mup_mu_full(ntheta,-mmax:mmax,-mmax:mmax), source=zeroc)
            allocate (test_explicit(ntheta,-mmax:mmax,-mmax:mmax), source=zeroc)
            do itheta=1,ntheta
                do mu=-mmax,mmax
                    do mup=-mmax,mmax
                        do ipsi=1,grid%npsi
                            do iphi=1,grid%nphi
                                test_explicit(itheta,mu,mup) = test_explicit(itheta,mu,mup)&
                                +foo_o( grid%indo(itheta,iphi,ipsi) ) &
                                *exp(ii*mup*grid%phiofnphi(iphi)) *exp(ii*mu*grid%psiofnpsi(ipsi))
                            end do
                        end do
                    end do
                end do
            end do


            do itheta=1,ntheta
                do iphi=1,nphi
                    do ipsi=1,npsi
                        fft2d%in(ipsi,iphi) = foo_o(grid%indo(itheta,iphi,ipsi))
                    end do
                end do
                fft2d_c%in(:,:)=fft2d%in
                call dfftw_execute (fft2d%plan)
                call dfftw_execute (fft2d_c%plan)


                ! do iphi=1,nphi
                !     do ipsi=1,npsi/2+1
                !         print*,ipsi,iphi, fft2d_c%out(ipsi,iphi), fft2d%out(ipsi,iphi)
                !     end do
                !     do ipsi=npsi/2+2,npsi
                !         print*,ipsi,iphi,fft2d_c%out(ipsi,iphi)
                !     end do
                ! end do
                foo_theta_mup_mu(itheta,0:mmax/molrotsymorder,0:mmax)   = CONJG( fft2d%out(:,1:mmax+1) )
                foo_theta_mup_mu(itheta,0:mmax/molrotsymorder,-mmax:-1) = CONJG( fft2d%out(:,mmax+2:) )

                foo_theta_mup_mu_full(itheta,0:mmax,0:mmax)     = fft2d_c%out(1:npsi/2+1, 1:nphi/2+1)
                foo_theta_mup_mu_full(itheta,0:mmax,-mmax:-1)   = fft2d_c%out(1:npsi/2+1, nphi/2+2:)
                foo_theta_mup_mu_full(itheta,-mmax:-1,0:mmax)   = fft2d_c%out(npsi/2+2:, 1:nphi/2+1)
                foo_theta_mup_mu_full(itheta,-mmax:-1,-mmax:-1) = fft2d_c%out(npsi/2+2:, nphi/2+2:)

            end do
            !
            !
            ! do itheta=1,ntheta
            !     do mup=-mmax,mmax
            !         do mu=-mmax,mmax
            !             print*,"itheta, mup, mu =",itheta,mup,mu
            !             print*,"explicit                      =",test_explicit(itheta,mu,mup)
            !             print*,"foo_theta_mup_mu_full(itheta,mu,mup)=",foo_theta_mup_mu_full(itheta,mu,mup)
            !             print*,"foo_theta_mup_mu_full reconstruit   =",(-1)**(-mup-mu)*conjg(foo_theta_mup_mu_full(itheta,-mu,-mup))," WWWW"
            !             if(mu>=0) print*,"foo_theta_mup_mu(itheta,mu,mup)     =",foo_theta_mup_mu(itheta,mu,mup)
            !             if(mu<0 ) print*,"foo_theta_mup_mu(itheta,mu,mup)     =",(-1)**(-mup-mu)*conjg(foo_theta_mup_mu(itheta,-mu,-mup))," reconstruit"
            !             print*,
            !         end do
            !     end do
            ! end do
            ! print*, "loop itheta,mup,mu FINISHED"
            ! print*,
            !
            allocate (proj_m_mup_mu(0:mmax,-mmax:mmax,-mmax:mmax), source=zeroc)
            do m=0,mmax
                do mup=-m,m
                    do mu=-m,m
                        do itheta=1,ntheta
                            proj_m_mup_mu(m,mup,mu) = proj_m_mup_mu(m,mup,mu) +&
                            foo_theta_mup_mu_full(itheta,mu,mup)*wtheta(itheta)*harm_sph(m,mup,mu,theta(itheta))*fm(m)
                        end do
                    end do
                end do
            end do

            ! Check relation 1.7
            do m=0,mmax
                do mup=-m,m
                    do mu=-m,m
                        if (proj_m_mup_mu(m,-mup,-mu) - (-1)**(mup+mu)*conjg(proj_m_mup_mu(m,mup,mu)) /= zeroc) then
                            print*, "Relation 1.7 is not satisfied"
                            error stop
                        end if
                    end do
                end do
            end do

            do m=0,mmax
                do mup=-m,m
                    do mu=0,m,molrotsymorder
                        ip = p3%p( m,mup,mu )
                        foo_p(ip)= sum(foo_theta_mup_mu(:,mu,mup)*p3%harm_sph(:,ip)*wtheta(:))*fm(m)
                        if (abs(foo_p(ip)-proj_m_mup_mu(m,mup,mu))>epsilon(1._dp)) then
                            print*, "foo_p(ip)/=proj_m_mup_mu(m,mup,mu)"
                            print*, "foo_p(ip)         =",foo_p(ip)
                            print*, "proj_m_mup_mu(m,mup,mu) =",proj_m_mup_mu(m,mup,mu)
                            error stop
                        end if
                    end do
                end do
            end do
        END FUNCTION euler2proj

        function proj2euler (foo_p) result (foo_o)
            implicit none
            real(dp) :: foo_o (1:grid%no)
            complex(dp), intent(in) :: foo_p(1:grid%np)
            complex(dp) :: foo_theta_mup_mu(1:grid%ntheta,0:grid%mmax/grid%molrotsymorder,-grid%mmax:grid%mmax) ! in fact foo_theta_mu_mup
            foo_theta_mup_mu = zeroc
            do itheta=1,ntheta
                do mup=-mmax,mmax
                    do mu=0,mmax/molrotsymorder
                        do m= MAX(ABS(mup),ABS(mu)), mmax
                            ip=p3%p(m,mup,mu)
                            foo_theta_mup_mu(itheta,mu,mup) = foo_theta_mup_mu(itheta,mu,mup)&
                            +foo_p(ip)*p3%harm_sph(itheta,ip)*fm(m)
                        end do
                    end do
                end do
            end do
            do itheta=1,ntheta
                ifft2d%in(:,1:mmax+1) = CONJG( foo_theta_mup_mu(itheta,:,0:mmax)   )
                ifft2d%in(:,mmax+2:)  = CONJG( foo_theta_mup_mu(itheta,:,-mmax:-1) )
                call dfftw_execute( ifft2d%plan )
                do iphi=1,nphi
                    do ipsi=1,npsi
                        foo_o (grid%indo(itheta,iphi,ipsi)) = ifft2d%out(ipsi,iphi) *(nphi*npsi)
                    end do
                end do
            end do
        END FUNCTION proj2euler

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

        subroutine test_routines_calcul_de_Rm_mup_mu_q
            !
            ! ici test de la routine de Choi
            !
            implicit none
                real(dp) :: qx, qy, qz
                complex(dp) :: R(0:mmax,-mmax:mmax,-mmax:mmax), Rlu(0:mmax,-mmax:mmax,-mmax:mmax)
                complex(dp), parameter :: ii=complex(0._dp,1._dp)
                real(dp), parameter :: epsdp=epsilon(1._dp)
                integer :: ix_q, iy_q, iz_q, m, mup, mu
                do ix_q=1,nx
                    do iy_q=1,ny
                        do iz_q=1,nz
                            qx=grid%kx(ix_q)
                            qy=grid%ky(iy_q)
                            qz=grid%kz(iz_q)
                            R= rotation_matrix_between_complex_spherical_harmonics_lu ([qx,qy,qz], mmax) ! exactement Luc
                            Rlu= rotation_matrix_between_complex_spherical_harmonics_lu ([qx,qy,qz], mmax) ! Lu
                            do m=0,mmax
                                do mup=-m,m
                                    do mu=-m,m
                                        if (abs(&
                                            Rlu(m,mup,mu)-harm_sph(m,mup,mu,thetaofq(qx,qy,qz))*exp(-ii*mup*phiofq(qx,qy)))&
                                                >1E-10) then
                                            print*, "dans test_routines_calcul_de_Rm_mup_mu_q on a Choi+Lu /= Wigner"
                                            print*, "ix_q, iy_q, iz_q",ix_q, iy_q, iz_q
                                            print*, "qx qy qz", qx, qy, qz
                                            print*, "m, mup, mu", m, mup, mu
                                            print*, "Rlu(m,mup,mu)",Rlu(m,mup,mu)
                                            print*, "harm_sph(m,mup,mu,thetaofq(qx,qy,qz))*exp(-ii*mup*phiofq(qx,qy))",&
                                            harm_sph(m,mup,mu,thetaofq(qx,qy,qz))*exp(-ii*mup*phiofq(qx,qy))
                                            print*, "diff=",&
                                            abs(Rlu(m,mup,mu)-harm_sph(m,mup,mu,thetaofq(qx,qy,qz))*exp(-ii*mup*phiofq(qx,qy)))
                                            print*, "epsilon=",epsdp
                                            error stop "oihjkkhblkuyuiykjhbkjhg76"
                                        end if
                                    end do
                                end do
                            end do
                        end do
                    end do
                end do
        end subroutine test_routines_calcul_de_Rm_mup_mu_q


    end subroutine energy_cproj

end module module_energy_cproj

module module_energy_cproj_mrso

    use iso_c_binding, only: C_PTR, C_INT, C_INT32_T, C_INTPTR_T, C_DOUBLE_COMPLEX, C_DOUBLE, C_FUNPTR, C_SIZE_T, C_FLOAT, &
                             C_FLOAT_COMPLEX, C_CHAR
    use precision_kinds, only: dp
    use module_grid, only: grid
    use module_solvent, only: solvent
    use module_rotation, only: rotation_matrix_between_complex_spherical_harmonics_lu

    implicit none

    !
    ! Everything here but the function is private. Don't worry about all these var in the module.
    !
    private

    !
    ! FFTW3 header - modern (fortran 2003) version. Expects iso_c_binding
    !
    include 'fftw3.f03'

    !
    ! Parameters
    !
    real(dp), parameter :: zero=0._dp
    complex(dp), parameter :: zeroc=(0._dp, 0._dp)
    complex(dp), parameter :: ii=(0._dp,1._dp)
    real(dp), parameter :: epsdp=epsilon(1._dp)
    real(dp), parameter :: pi=acos(-1._dp)
    real(dp), parameter :: eightpisq=8._dp*pi**2
    real(dp), parameter :: twopi=2._dp*pi

    !
    ! Direct correlation functions from Luc
    !
    type :: c_type
        logical :: isok = .false.
        real(dp) :: q
        integer :: np ! number of projections for the direct correlation function. Called alpha  by Luc
        integer :: nq ! number of q points in c(q)
        integer, allocatable :: ip(:,:,:,:,:) ! index of m, n, mu2, nu2, khi => (0:mmax,0:mmax,-mmax/mrso:mmax/mrso,-mmax/mrso:mmax/mrso,-mmax:mmax) ! m n mu nu khi
        real(dp) :: dq
        real(dp), allocatable :: normq(:)
        integer, allocatable :: m(:), n(:), mu(:), nu(:), khi(:)
        complex(dp), allocatable :: mnmunukhi_q(:,:)
    end type c_type
    type(c_type), protected :: c

    complex(dp), allocatable, protected :: deltarho_p(:,:,:,:) ! deltarho_p(np,nx,ny,nz)
    complex(dp), allocatable, protected :: deltarho_p_q(:)
    complex(dp), allocatable, protected :: deltarho_p_mq(:)
    complex(dp), allocatable, protected :: gamma_p_q(:)
    complex(dp), allocatable, protected :: gamma_p_mq(:)

    type :: p3_type
        real(dp), allocatable :: wigner_small_d(:,:) ! tabulation des harmoniques sphériques r(m,mup,mu,theta) en un tableau r(itheta,p)
        integer, allocatable :: p(:,:,:) ! index of the projection corresponding to m, mup, mu
        integer, allocatable :: m(:) ! m for projection 1 to np
        integer, allocatable :: mup(:) ! mup for projection 1 to np. mup corresponds to phi
        integer, allocatable :: mu(:) ! mu for projection 1 to np. mu corresponds to psi
    end type p3_type
    type (p3_type), protected :: p3

    type :: fft_type
        logical :: planned = .false.
        type(c_ptr) :: plan3dp, plan3dm ! plans for 2D and 3D FFTs with sign +(p) or -(m) in the exponential
    end type fft_type
    type(fft_type), protected :: fft
    complex(dp), allocatable :: c3d(:,:,:)
    complex(dp), allocatable, protected :: R(:,:,:) ! Table of generalized spherical harmonics of m, mup, mu
    complex(dp), allocatable, protected :: Rmmupmu(:) ! same as R(m,mup,mu) with only used values of m,mup and mu

    public :: energy_cproj_mrso

contains

    subroutine energy_cproj_mrso (ff, df, print_timers)
        use omp_lib
        use precision_kinds, only: dp
        use module_grid, only: grid
        use module_thermo, only: thermo
        use module_orientation_projection_transform, only: angl2proj, proj2angl
        implicit none
        real(dp), intent(out) :: ff
        real(dp), contiguous, intent(inout), optional :: df(:,:,:,:,:) ! x y z o s
        logical, intent(in), optional :: print_timers
        real(dp) :: dv, kT
        logical :: q_eq_mq
        integer :: ix, iy, iz, ix_q, iy_q, iz_q, ix_mq, iy_mq, iz_mq, ip, ip2
        integer :: nx, ny, nz, np, no, ns, ntheta, nphi, npsi, mmax, mrso
        integer :: m, n, mu, nu, khi, mup, ia, iq, io, mu2, nu2
        complex(dp) :: gamma_m_khi_mu_q, gamma_m_khi_mu_mq
        real(dp) :: q(3), lx, ly, lz, rho0
        real(dp) :: theta(grid%ntheta), wtheta(grid%ntheta)
        logical, allocatable :: gamma_p_isok(:,:,:)
        real :: time(20)
        real(dp) :: vexc(grid%no), rho, xi
        real(dp), parameter :: fourpisq = 4._dp*acos(-1._dp)**2
        real :: total_time_in_subroutine
        complex(dp) :: R_loc(-5:5), deltarho_p_q_loc, deltarho_p_mq_loc
        integer :: ip2_loc(-5:5)

        call cpu_time (time(1))

        ntheta = grid%ntheta
        nphi   = grid%nphi
        npsi   = grid%npsi
        mmax   = grid%mmax
        ! In the case of C∞v or any symetry with ∞, the user enters 0 for ∞ in the input file in molrotsymorder.
        ! To ease the loops over all mu=-m/mrso,m/mrso, it is easier to make it large, like 100, so that m/mrso:=0.
        ! With this trick, we have the same loops as before, with only one value of mu, 0.
        if( grid%molrotsymorder == 0 ) then
            mrso = 100
        else
            mrso = grid%molrotsymorder
        end if
        lx   = grid%lx
        ly   = grid%ly
        lz   = grid%lz
        kT   = thermo%kbT
        dv   = grid%dv
        nx   = grid%nx
        ny   = grid%ny
        nz   = grid%nz
        np   = grid%np
        no   = grid%no
        ns   = solvent(1)%nspec
        rho0 = solvent(1)%rho0


        if (.not. allocated (deltarho_p) ) then
            allocate (deltarho_p(np,nx,ny,nz) ,source=zeroc)
            allocate (deltarho_p_q(np) ,source=zeroc)
            allocate (deltarho_p_mq(np) ,source=zeroc)
            allocate (gamma_p_q(np), source=zeroc)
            allocate (gamma_p_mq(np), source=zeroc)
        end if

        if (.not. allocated (gamma_p_isok) ) allocate (gamma_p_isok(nx,ny,nz), source=.false.)

        call cpu_time (time(2))

        ! 1/ get Δρ(r,ω)                   :    Δρ(r,ω)          =     ρ0·(ξ²-1)
        ! 2/ projection                    :    Δρ(r,ω)          =>    Δρ^m_mup,mu(r)
        ! 3/ FFT3D-C2C                     :    Δρ^m_mup,mu(r)   =>    Δρ^m_mup,mu(q)
        ! 4/ rotate to q frame             :    Δρ^m_mup,mu(q)   =>    Δρ'^m_mup,mu(q)
        ! 5/ OZ                            :    Δρ'^m_mup,mu(q)  =>    ɣ'^m_mup,mu(q)
        ! 6/ rotate back to fixed frame    :    ɣ'^m_mup,mu(q)   =>    ɣ^m_mup,mu(q)
        ! 7/ FFT3D-C2D                     :    ɣ^m_mup,mu(q)    =>    ɣ^m_mup,mu(r)
        ! 8/ gather projections            :    ɣ^m_mup,mu(r)    =>    ɣ(r,ω)


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
            allocate ( p3%p(0:mmax,-mmax:mmax, 0:mmax/mrso) ,source=-huge(1)) ! Dans p3%p, on met mu2=2*mu, pas mu
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
                error stop "Invalid values of m, mup or mu in energy_cproj_mrso"
            end if

            !
            ! Print all projections that will be used
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
        ! L'initialisation des racines et des poids de la quadrature angulaire a été faite beaucoup plus tot
        ! dans le module_grid.
        ! On fait attention qu'on a deux tableaux des racines et des poids:
        ! - le premier tableau, grid%theteaofntheta, est un tableau de taille  le nombre de theta (ntheta)
        ! - le deuxième, grid%theta, est de taille le nombre total d'angle. Il retourne le theta qui correspond à l'indice de chaque angle.
        ! Si on veut la liste de tous les theta, on utilisera donc le tableau grid%thetaofntheta.
        !
        theta=grid%thetaofntheta
        wtheta=grid%wthetaofntheta

        !
        ! Tabulate generalized spherical harmonics in array p3%wigner_small_d(theta,proj)
        ! where theta can be any of the GaussLegendre integration roots for theta
        ! where proj is an index related to a tuple {m,mup,mu}
        ! Blum's notation :
        ! m is related to theta
        ! mup is related to phi
        ! mu is related to psi
        ! TODO: a remplacer par la routine de luc, et utiliser la notation alpha plutot que m,mup,mu a ce moment
        !
        ! call test_routines_calcul_de_Rm_mup_mu_q
        if (.not. allocated(p3%wigner_small_d)) then
            block
                use module_wigner_d, only: wigner_small_d
                integer :: p,m,mup,mu,i
                allocate ( p3%wigner_small_d(1:ntheta, 1:np) ,source=0._dp)
                do p=1,np
                    m = p3%m(p)
                    mup = p3%mup(p)
                    mu = p3%mu(p)
                    do i=1,ntheta
                        ! Pour chaque theta, calcule la fonction de Wigner-d correspondant à toutes les projections avec la méthode de Wigner.
                        p3%wigner_small_d(i,p) = wigner_small_d(m,mup,mu,theta(i))
                    end do
                end do
            end block
        end if

        call cpu_time (time(4))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! PROJECTION DE Δρ(r,ω) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !
        ! On passe tout de suite en projection, dans le repère cartesien
        ! ces projections sont complexes, mais la symetrie hermitienne
        ! permet de ne garder que les mup>=0 ou les mu>=0.
        ! On choisit les mu>=0 (cf doc de Luc)
        !

        ! 1/ ON PREPARE LES FFT
        ! Pay attention to the definition of "forward" and "backward" for the FFTW library. Source: http://www.fftw.org/doc/Real_002ddata-DFTs.html
        ! FFTW_FORWARD = -1
        ! FFTW_BACKWARD = +1
        ! R2C == FORWARD = -1
        ! C2R == BACKWARD = +1
        if (.not. fft%planned) then
            allocate( c3d(nx,ny,nz), source=(0._dp,0._dp)  )
            select case(dp)
            case(c_double)
                call dfftw_plan_dft_3d(      fft%plan3dp, nx, ny, nz, c3d, c3d, +1, FFTW_MEASURE) ! TODO CHECK ESTIMATE VS REST & IS IT WORTH CHANGEING THE PLAN FLAG FOR DIFFERENT np ? Certainly!
                call dfftw_plan_dft_3d(      fft%plan3dm, nx, ny, nz, c3d, c3d, -1, FFTW_MEASURE )
            case(c_float)
                call sfftw_plan_dft_3d(      fft%plan3dp, nx, ny, nz, c3d, c3d, +1, FFTW_MEASURE) ! TODO CHECK ESTIMATE VS REST & IS IT WORTH CHANGEING THE PLAN FLAG FOR DIFFERENT np ? Certainly!
                call sfftw_plan_dft_3d(      fft%plan3dm, nx, ny, nz, c3d, c3d, -1, FFTW_MEASURE )
            end select
            fft%planned =.true.
        end if

        call cpu_time (time(5))

        ! 2/ ON PROJETTE delta_rho

        block
            real(dp) :: o(no)
            do iz=1,nz
                do iy=1,ny
                    do ix=1,nx
                        o = rho0*(solvent(1)%xi(:,ix,iy,iz)**2 -1._dp)
                        call angl2proj( o, deltarho_p(:,ix,iy,iz) )
                    end do
                end do
            end do
        end block

        call cpu_time (time(6))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! FFT SPATIALE [ Δρ^m_mup,mu(r) ] !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !
        ! On a les projections sur la grille cartesienne
        ! On veut passer dans l'espace de Fourier pour calculer les convolutions spatiales
        ! On fait donc une FFT 3D pour chacune des projections.
        ! Les projections sont complexes, il s'agit donc d'une FFT3D C2C habituelle : Il n'y a pas de symétrie hermitienne.
        !
        block
            integer :: ip
            do ip=1,np
                c3d = deltarho_p(ip,:,:,:)
                select case(dp);
                case(c_double)
                    call dfftw_execute_dft(fft%plan3dp, c3d, c3d)
                case(c_float)
                    call sfftw_execute_dft(fft%plan3dp, c3d, c3d)
                end select
                deltarho_p(ip,:,:,:) = c3d
            end do
        end block
        call cpu_time (time(7))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ROTATION VERS REPAIRE MOLECULAIRE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        +                       OZ !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        + RETOUR VERS REPAIRE FIXE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
        ! And the inverse index that gives vector -q: see function grid%ix_mq(ix_q) of module_grid
        !

        call cpu_time (time(8))

        !
        ! Read Luc's direct correlation function c^{m,n}_{mu,nu_,chi}(|q|)
        ! projected on generalized spherical harmonics
        ! in the intermolecular frame
        ! normq is norm of q, |q|, that correspond to the index iq in ck(ia,iq)
        !
        if (.not.c%isok) then
            block
                use module_read_c_luc, only: read_c_luc
                real(dp) :: qmaxnecessary
                qmaxnecessary = norm2([maxval(grid%kx(1:nx)), maxval(grid%ky(1:ny)), maxval(grid%kz(1:nz/2+1))])
                call read_c_luc(c%mnmunukhi_q,mmax,mrso,qmaxnecessary,c%np,c%nq,c%dq,c%m,c%n,c%mu,c%nu,c%khi,c%ip)
            end block
            c%mnmunukhi_q=conjg(c%mnmunukhi_q) ! this is strange, but certainly due to some error or misunderstanding with Luc. It does not appear in Luc's document.
            c%isok=.true.
        end if

        call cpu_time (time(9))

        if (.not. allocated(R) ) then
            allocate ( R(0:mmax,-mmax:mmax,-mmax:mmax) ,source=(0._dp,0._dp) )
            allocate( Rmmupmu(np), source=(0._dp,0._dp))
        end if

        !
        ! For all vectors q and -q handled simultaneously.
        ! ix_q,iy_q,iz_q are the coordinates of vector q, while ix_mq,iy_mq_iz_mq are those of vector -q
        !
        gamma_p_isok = .false.

        !$OMP PARALLEL DO DEFAULT(FIRSTPRIVATE) SHARED(gamma_p_isok,c,deltarho_p,grid,ck)
        do iz_q=1,nz/2+1
            iz_mq = grid%iz_mq(iz_q)
            do iy_q=1,ny
                iy_mq = grid%iy_mq(iy_q)
                do ix_q=1,nx
                    ix_mq = grid%ix_mq(ix_q)

                    !
                    ! gamma_p_isok is a logical array. If gamma(ix_q,iy_q,iz_q) has already been calculated, it is .true.
                    !
                    if ( gamma_p_isok(ix_q,iy_q,iz_q) .and. gamma_p_isok(ix_mq, iy_mq,iz_mq) ) cycle

                    !
                    ! cartesian coordinates of vector q in lab frame
                    !
                    q = [grid%kx(ix_q), grid%ky(iy_q), grid%kz(iz_q)]

                    !
                    ! pay attention to the special case(s) where q=-q
                    !
                    if (ix_mq==ix_q .and. iy_mq==iy_q .and. iz_mq==iz_q) then ! this should only happen for ix=1 and ix=nx/2
                        q_eq_mq=.true.
                    else
                        q_eq_mq=.false.
                    end if

                    !
                    ! Prepare R^m_mup_khi(q)
                    !
                    R = rotation_matrix_between_complex_spherical_harmonics_lu ( mmax, q)
                    where( abs(R)<=epsdp ) R = (0._dp,0._dp)
                    Rmmupmu(1:np) = [(    R(p3%m(ip),p3%mup(ip),p3%mu(ip))   ,ip=1,np    )]

                    ! Eq. 1.23 We don't need to compute gshrot for -q since there are symetries between R(q) and R(-q).
                    ! Thus, we do q and -q at the same time. That's the most important point in doing all q but half of mu.
                    ! Lu decided to do all mu but half of q in her code

                    !
                    ! Rotation to molecular (q) frame
                    !
                    gamma_p_q = deltarho_p(:,ix_q,iy_q,iz_q) ! this a temporary array
                    gamma_p_mq = deltarho_p(:,ix_mq,iy_mq,iz_mq) ! this a temporary array
                    ip=0
                    do m=0,mmax
                        do khi=-m,m
                            R_loc(-m:m)=R(m,-m:m,khi)
                            do mu2=0,m/mrso
                                ip=ip+1
                                ip2_loc(-m:m)=p3%p(m,-m:m,mu2)  ! Optimization added 6th of June 2016. Makes the code ugly but improves locality... :-/
                                ! Equation 1.22
                                deltarho_p_q_loc  = (0._dp,0._dp)
                                deltarho_p_mq_loc = (0._dp,0._dp)
                                do mup=-m,m
                                    ip2=ip2_loc(mup)
                                    deltarho_p_q_loc  = deltarho_p_q_loc  + gamma_p_q(ip2)  *R_loc(mup)
                                    deltarho_p_mq_loc = deltarho_p_mq_loc + gamma_p_mq(ip2) *R_loc(mup)
                                end do
                                deltarho_p_q(ip) = deltarho_p_q_loc
                                deltarho_p_mq(ip) = deltarho_p_mq_loc *(-1)**m
                            end do
                        end do
                    end do

                    !
                    ! c^{m,n}_{mu,nu,chi}(|q|) is tabulated for c%nq values of |q|.
                    ! Find the tabulated value that is closest to |q|. Its index is iq.
                    ! Note |q| = |-q| so iq is the same for both vectors.
                    !
                    iq = int( norm2(q) /c%dq +0.5) +1

                    !
                    ! Ornstein-Zernike in the molecular frame
                    ! We do OZ for q and -q at the same time
                    !
                    gamma_p_q  = zeroc
                    gamma_p_mq = zeroc !1:np

                    do khi=-mmax,mmax
                        do m=abs(khi),mmax

                            do mu=0,m,mrso ! not -m,m a cause des symetries ! EST CE QU'IL Y A UNE RAISON POUR QUE LES mu IMPAIRES SOIENT NON NULS ICI ?
                              mu2=mu/mrso

                                gamma_m_khi_mu_q= zeroc
                                gamma_m_khi_mu_mq = zeroc

                                do n=abs(khi),mmax

                                    do nu= -mrso*(n/mrso), mrso*(n/mrso), mrso   ! imaginons n=3, -n,n,mrso  ferait nu=-3,-1,1,3 mais en faisant /mrso puis *mrso, ça fait -2,0,2 as expected
                                      nu2=nu/mrso

                                        ia = c%ip(m,n,mu2,nu2,khi) ! the index of the projection of c(q). 1<=ia<na

                                        ip = p3%p(n,khi,abs(nu2))

                                        if (nu<0) then ! no problem with delta rho (n, khi, -nu) since -nu>0. Thus, we apply eq. 1.30 directly
                                            gamma_m_khi_mu_q  = gamma_m_khi_mu_q  + (-1)**(khi+nu) *c%mnmunukhi_q(ia,iq) *deltarho_p_q(ip)
                                            gamma_m_khi_mu_mq = gamma_m_khi_mu_mq + (-1)**(khi+nu) *c%mnmunukhi_q(ia,iq) *deltarho_p_mq(ip)
                                        else
                                            gamma_m_khi_mu_q = gamma_m_khi_mu_q  + (-1)**(n) *c%mnmunukhi_q(ia,iq) *conjg(deltarho_p_mq(ip))
                                            gamma_m_khi_mu_mq= gamma_m_khi_mu_mq + (-1)**(n) *c%mnmunukhi_q(ia,iq) *conjg(deltarho_p_q(ip))
                                        end if

                                    end do
                                end do

                                ip=p3%p(m,khi,mu2)

                                gamma_p_q(ip)  = gamma_m_khi_mu_q
                                gamma_p_mq(ip) = gamma_m_khi_mu_mq
                            end do
                        end do
                    end do

                    !
                    ! Rotation from molecular frame to fix frame
                    !
                    R = conjg(R) ! le passage retour au repaire fixe se fait avec simplement le conjugue complexe de l'harm sph generalisee
                    Rmmupmu =  conjg(Rmmupmu)
                    ! we use deltarho_p_q and deltarho_p_mq as temp arrays since they're not used after MOZ


                    deltarho_p_q = (0._dp,0._dp)
                    deltarho_p_mq = (0._dp,0._dp)
                    ip=0
                    do m=0,mmax
                        do mup=-m,m
                            do mu2=0,m/mrso
                                ip=ip+1
                                ! Equation 1.22
                                do khi=-m,m
                                    ip2=p3%p(m,khi,mu2)
                                    ! Rmmupmu( p3%p(m,mup,khi) )
                                    deltarho_p_q(ip)  = deltarho_p_q(ip)  + gamma_p_q (ip2) *R(m,mup,khi)
                                    deltarho_p_mq(ip) = deltarho_p_mq(ip) + gamma_p_mq(ip2) *R(m,mup,khi)
                                end do
                                deltarho_p_mq(ip) = deltarho_p_mq(ip) *(-1)**m
                            end do
                        end do
                    end do

                    !
                    ! For vector q,
                    !
                    deltarho_p (:, ix_q, iy_q, iz_q) = deltarho_p_q
                    gamma_p_isok(ix_q,iy_q,iz_q)=.true.

                    !
                    ! Again, pay attention to the singular mid-k point
                    !
                    if( q_eq_mq .and. (ix_q==nx/2+1.or.iy_q==ny/2+1.or.iz_q==nz/2+1)) then
                        deltarho_p(1:np, ix_mq, iy_mq, iz_mq) = conjg(deltarho_p_mq)
                    else
                        deltarho_p(1:np, ix_mq, iy_mq, iz_mq) = deltarho_p_mq
                    end if
                    gamma_p_isok(ix_mq,iy_mq,iz_mq)=.true.

                end do
            end do
        end do
        !$OMP END PARALLEL DO


        call cpu_time(time(10))

        if (.not.all(gamma_p_isok.eqv..true.)) then
            error stop "not all gamma_p(projections,ix,iy,iz) have not been computed"
        end if

        call cpu_time(time(11))


        !
        ! FFT3D from Fourier space to real space
        !
        block
            integer :: ip
            do ip=1,np
                c3d = deltarho_p(ip,:,:,:)
                select case(dp)
                case(c_double)
                    call dfftw_execute_dft( fft%plan3dm, c3d, c3d )
                case(c_float)
                    call sfftw_execute_dft( fft%plan3dm, c3d, c3d )
                end select
                deltarho_p(ip,:,:,:) = c3d/real(nx*ny*nz,dp)
            end do
        end block

        call cpu_time(time(12))

        !
        ! Gather projections into gamma
        ! Note that gamma==df
        !
        if(present(df)) then
            ff=0._dp
            do iz=1,nz
                do iy=1,ny
                    do ix=1,nx
                        call proj2angl( deltarho_p(:,ix,iy,iz), vexc)
                        vexc = -kT*grid%w*vexc*fourpisq  /solvent(1)%n0  ! /0.0333 vient du c de luc qui est un n0.c
                        do io=1,no
                            xi = solvent(1)%xi(io,ix,iy,iz)
                            rho = xi**2*rho0
                            ff = ff + (rho-rho0)*0.5_dp*vexc(io)*dv
                            df(io,ix,iy,iz,1) = df(io,ix,iy,iz,1) + 2._dp*vexc(io)*rho0*xi
                        end do
                    end do
                end do
            end do
        else
            ff=0._dp
            do iz=1,nz
                do iy=1,ny
                    do ix=1,nx
                        call proj2angl( deltarho_p(:,ix,iy,iz), vexc)
                        vexc = -kT*grid%w*vexc*fourpisq  /solvent(1)%n0  ! /0.0333 vient du c de luc qui est un n0.c
                        do io=1,no
                            rho = solvent(1)%xi(io,ix,iy,iz)**2*rho0
                            ff = ff + (rho-rho0)*0.5_dp*vexc(io)*dv
                        end do
                    end do
                end do
            end do
        end if
        call cpu_time(time(13))

        total_time_in_subroutine = time(13)-time(1)

        if(present(print_timers)) then
            if( print_timers) then
print*, "                  |"
print*, "                  | read ck                                      ", time(9)-time(8)  ,"sec (",nint((time(9) -time(8 ))/total_time_in_subroutine*100),"%)"
print*, "                  | tabulate spherical harmonics                 ", time(4)-time(3)  ,"sec (",nint((time(4) -time(3 ))/total_time_in_subroutine*100),"%)"
print*, "                  | plan FFTs                                    ", time(5)-time(4)  ,"sec (",nint((time(5) -time(4 ))/total_time_in_subroutine*100),"%)"
print*, "                  | angl2proj                                    ", time(6)-time(5)  ,"sec (",nint((time(6) -time(5 ))/total_time_in_subroutine*100),"%)"
print*, "                  | FFT  @Δρ^m_mup_mu(k)                         ", time(7)-time(6)  ,"sec (",nint((time(7) -time(6 ))/total_time_in_subroutine*100),"%)"
print*, "                  | rotation to molecular frame + OZ + inv rot   ", time(10)-time(9) ,"sec (",nint((time(10)-time(9 ))/total_time_in_subroutine*100),"%)"
print*, "                  | FFT⁻¹@Δρ^m_mup_mu(k)                         ", time(12)-time(11),"sec (",nint((time(12)-time(11))/total_time_in_subroutine*100),"%)"
print*, "                  | proj2angl + sum to ff and df                 ", time(13)-time(12),"sec (",nint((time(13)-time(12))/total_time_in_subroutine*100),"%)"
print*, "                  |"
            end if
        end if

    end subroutine energy_cproj_mrso

end module module_energy_cproj_mrso

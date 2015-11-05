module module_energy_cproj
    use iso_c_binding
    use module_grid, only: grid
    use module_solvent, only: solvent
    use mathematica, only: fact
    implicit none
    private
    include 'fftw3.f03'                 ! Needed by FFTW3
    integer, parameter :: dp=C_DOUBLE   ! Corresponds to double precision in usual architectures
    complex(dp), parameter :: imag=(0_dp,1_dp), zeroc=complex(0._dp,0._dp)

    integer :: m,n,mup,mu,nu,khi,i,ix,iy,iz,iix,iiy,iiz,p,itheta,ip,iq,ia,iphi,ipsi,io,no

    !
    ! Main physical parameters
    !
    LOGICAL, PARAMETER :: debug=.true.
    integer, PARAMETER :: mmax=0
    integer, PARAMETER :: nx=64, ny=nx, nz=nx
    real(dp), parameter :: lx=25._dp, ly=lx, lz=lx ! size of the supercell in angstrom
    integer, PARAMETER :: molrotsymorder=1


    integer, parameter :: nproj = SUM( [( [( [( 1 ,mu=0,m,molrotsymorder)], mup=-m,m)], m=0,mmax)] )
    integer, parameter :: ntheta =   mmax+1
    integer, parameter :: nphi   = 2*mmax+1
    integer, parameter :: npsi   = 2*(mmax/molrotsymorder)+1

    !
    ! Arrays to translate index of vector to tuple of projections, and inverse
    ! alpha(m,mup,mu), m(alpha), mup(alpha), mu(alpha)
    !
    INTEGER ::  indp(0:mmax,-mmax:mmax,0:mmax) ! alpha(m,mup,mu)
    INTEGER ::  im(nproj)                       ! m(alpha)
    INTEGER ::  imup(nproj)                    ! mup(alpha)
    INTEGER ::  imu(nproj)                        ! mu(alpha)

    !
    ! Direct correlation functions from Luc
    !
    integer, parameter :: nq=200 ! nombre de valeurs de |q| dans le fichier de c(q) de luc
    integer, parameter :: nalpha=(1+mmax)*(1+2*mmax)*(3+2*mmax)*(5+4*mmax*(2+mmax))/15 ! number of projections in file c(q)
    integer :: ialpha(0:mmax,0:mmax,-mmax:mmax,-mmax:mmax,-mmax:mmax) ! m n mu nu khi
    type :: ck_type
        logical :: isok = .false.
    end type
    type (ck_type) :: myck
    complex(dp) :: ck(nalpha,nq), ck_q(nalpha)
    real(dp) :: normq(nq) ! contient pour chaque indice de c(iq) la norme de q associée

    real(dp), parameter :: pi=acos(-1._dp), twopi=2._dp*pi
    real(dp), parameter :: dphi = twopi/real(nphi,dp)
    real(dp), parameter :: phi(nphi) = [(   (i-1)*dphi   , i=1,nphi )]
    real(dp), parameter :: wphi(nphi) = 1./nphi
    real(dp), parameter :: dpsi = twopi/real(npsi*molrotsymorder,dp)
    real(dp), parameter :: psi(npsi) = [(   (i-1)*dpsi   , i=1,npsi )]
    real(dp), parameter :: wpsi(npsi) = 1./npsi
    real(dp) :: theta(ntheta), wtheta(ntheta)

    integer :: exitstatus
    real(dp), parameter :: fm(0:mmax) = [( sqrt(real(2*m+1,dp))/real(nphi*npsi,dp) ,m=0,mmax  )]
    type :: fft2d_type
        type(c_ptr) :: plan
        real(dp)    :: in(npsi,nphi)
        complex(dp) :: out(npsi/2+1,nphi)
    end type
    type :: ifft2d_type
        type(c_ptr) :: plan
        complex(dp) :: in(npsi/2+1,nphi)
        real(dp)    :: out(npsi,nphi)
    end type
    type :: fft3d_type
        type(c_ptr) :: plan
        complex(dp) :: in(nx,ny,nz)
    end type
    type(  fft2d_type ) ::  fft2d
    type( ifft2d_type ) :: ifft2d
    type( fft3d_type ) :: fft3d
    real(dp) :: tharm_sph(ntheta,nproj), t0, t1, t2, t3, t4, t5, ta, tb, tc, t6, t7, t8
    real(dp) :: q(3), tqx(nx), tqy(ny), tqz(nz) ! vector q and its components tabulated
    integer :: tqxi(nx), tqyi(ny), tqzi(nz)
    complex(dp) :: deltarho_m_mup_mu_q   (0:mmax,-mmax:mmax,-mmax:mmax)
    complex(dp) :: deltarho_m_mup_mu_mq(0:mmax,-mmax:mmax,-mmax:mmax)
    complex(dp) :: xxx_temp_array_of_m_khi_mu_q   (0:mmax,-mmax:mmax,-mmax:mmax)
    complex(dp) :: xxx_temp_array_of_m_khi_mu_mq(0:mmax,-mmax:mmax,-mmax:mmax)
    complex(dp) :: deltarho_p_q(1:nproj), deltarho_p_mq(1:nproj)
    complex(dp) :: gamma_p_q(1:nproj), gamma_p_mq(1:nproj)
    complex(dp) :: gamma_m_khi_mu_q, gamma_m_khi_mu_mq
    integer(1), parameter :: errorgenerator=-1

    public :: energy_cproj

contains

    subroutine energy_cproj (ff,df)
        use precision_kinds, only: dp
        use module_grid, only: grid
        use module_thermo, only: thermo
        implicit none
        real(dp), parameter :: zero=0._dp
        complex(dp), parameter :: zeroc=complex(0._dp, 0._dp)
        real(dp), intent(out) :: ff
        real(dp), intent(inout) :: df(:,:,:,:,:)
        real(dp) :: dv, dq, kT
        complex(dp) :: ffc
        logical :: q_eq_mq
        integer :: ix, iy, iz, ix_q, iy_q, iz_q, ix_mq, iy_mq, iz_mq
        real(dp), allocatable :: gamma(:,:,:,:,:,:), gamma_o(:,:,:,:), deltarho(:,:,:,:,:,:)
        complex(dp), allocatable :: deltarho_p(:,:,:,:), deltarho_p_q(:), deltarho_p_mq(:), gamma_p(:,:,:,:)
        complex(dp), allocatable :: deltarho_m_mup_mu_q(:,:,:), deltarho_m_mup_mu_mq(:,:,:)


        kT=thermo%kbT
        dv=grid%dv
        no=grid%no

        allocate (deltarho(nx,ny,nz,ntheta,nphi,npsi) ,source=zero)
        allocate (gamma(nx,ny,nz,ntheta,nphi,npsi) ,source=zero)
        allocate (gamma_o(nx,ny,nz,no) ,source=zero)

        allocate (gamma_p(nproj,nx,ny,nz) ,source=zeroc)
        allocate (deltarho_p(nproj,nx,ny,nz) ,source=zeroc)
        allocate (deltarho_p_q(nproj) ,source=zeroc)
        allocate (deltarho_p_mq(nproj) ,source=zeroc)
        allocate (deltarho_m_mup_mu_q(m,mup,mu) ,source=zeroc)
        allocate (deltarho_m_mup_mu_mq(m,mup,mu) ,source=zeroc)


        ! 1/ get deltarho = rho-rho0
        ! 2/ project deltarho => deltarho_p (use FFT2D-R2C)
        ! 3/ FFT3D-C2C deltarho_p(x) => deltarho_p(q)
        ! 4/ rotate to q frame
        ! 5/ OZ: deltarho_p(q) => gamma_p(q)
        ! 6/ rotate back to fixed frame
        ! 7/ FFT3D-C2D gamma_p(q) => gamma_p(r)
        ! 8/ gather projections: gamma_p(r) => gamma(r)



        ! PRINT*,"========"
        ! print*,"mmax              =", mmax
        ! print*,"molrotsymorder    =", molrotsymorder
        ! print*,"nproj for density =", nproj
        ! print*,"nproj for c(q)    =", nalpha
        ! print*,"ntheta            =", ntheta
        ! print*,"nphi              =", nphi
        ! print*,"npsi              =", npsi
        ! print*,"nx,ny,nz          =", [nx,ny,nz]
        ! PRINT*,"lx,ly,lz          =", [lx,ly,lz]
        ! PRINT*,"========"
        ! PRINT*,

        if(mmax/=grid%mmax) stop "mmax pas bon"
        if(nx/=grid%nx) stop "nx pas bon"
        if(ny/=grid%ny) stop "ny pas bon"
        if(nz/=grid%nz) stop "nz pas bon"
        if(lx/=grid%lx) stop "lx pas bon"
        if(ly/=grid%ly) stop "ly pas bon"
        if(lz/=grid%lz) stop "lz pas bon"
        if(molrotsymorder/=grid%molrotsymorder) stop "molrotsymorder pas bon"

        call cpu_time(tc)

        !
        ! Toutes les projections sont stockées dans un meme vecteur de taille nproj
        ! indp(m,mup,mu) lie le triplet (m,mup,mu) a l'unique indice de projection ip
        ! im(p) donne inversement le m correspondant à l'indice p dans le tableau des projections
        ! On a donc indp(im(p),imup(p),imu(p)) == p
        !
        indp = -huge(1)
        im    = -huge(1)
        imup  = -huge(1)
        imu   = -huge(1)
        p=0
        do m=0,mmax
            do mup=-m,m
                !
                ! We chose to have indp and imu to contain the true value of mu, not mu/molrotsymorder.
                ! Thus, we loop over mu with steps of molrotsymorder
                ! Thus, we have holes in indp, but all calls to indp and imu should be done with this true mu
                ! For instance, if molrotsymorder=2 (for water or any C2V molecule), any call to indp(m,mup,mu) with mu=1 will return -999
                ! Also, the array imu only contains even values (des valeurs paires)
                !
                do mu=0,m,molrotsymorder
                    p=p+1
                    IF (p > nproj) ERROR STOP "p > nproj at line 166"
                    indp(m,mup,mu) = p
                    im(p) = m
                    imup(p) = mup
                    imu(p) = mu
                end do
            end do
        end do

        if (debug) then
            if (p /= nproj) error stop "p /= nproj in energy_cproj"
            if ( any(abs(im)>mmax) .or. any(abs(imup)>mmax) .or. any(abs(imu)>mmax) ) then
                print*, "tabulated m, mup or mu have incorrect values"
                error stop
            end if
        end if

        !
        ! Print all projections that will be kept in memory
        !
        if (debug) then
            open(11, file="output/nonzero-projections.out")
            write(11,*)"Non-zero projections:"
            write(11,*)"        index         m          mup          mu"
            write(11,*)"        -----        ---         ---          --"
            do p=1,nproj
                write(11,*) p, im(p), imup(p), imu(p)
            end do
            close(11)
        end if

        !
        ! Initie poids et racines de l'integration angulaire par la methode de Gauss-Legendre
        ! Les racines sont les racines d'un polynomes en cos(theta).
        ! Comme nous voulons les theta correspondant, on en prend l'arccos(cos(theta))
        !
        block
            real(dp) :: costheta(1:ntheta)
            call gauss_legendre( ntheta, costheta, wtheta, exitstatus ) ! retourne les ntheta racines (theta(:)) et poids (wtheta(:))
            theta = acos(costheta)
        end block


                !
                ! Tabulate generalized spherical harmonics in array tharm_sph(theta,proj)
                ! where theta can be any of the GaussLegendre integration roots for theta
                ! where proj is an index related to a tuple {m,mup,mu}
                ! Blum's notation :
                ! m is related to theta
                ! mup is related to phi
                ! mu is related to psi
                ! TODO: a remplacer par la routine de luc, et utiliser la notation alpha plutot que m,mup,mu a ce moment
                !
                tharm_sph = 0._dp
                do ip=1,nproj
                    m = im(ip)
                    mup = imup(ip)
                    mu = imu(ip)
                    do itheta=1,ntheta
                        tharm_sph(itheta,ip) = harm_sph(m,mup,mu,theta(itheta))
                    end do
                end do
                if (mmax==0) then
                    if (size(tharm_sph,1) /= size(tharm_sph,2) .or. tharm_sph(1,1)/=1._dp) then
                        print*, "in energy_cproj, pour mmax==0, tharm_sph ne vaut pas 1."
                        print*, "tharm_sph(itheta=1,ip=1) =", tharm_sph(1,1)
                        error stop
                    end if
                end if


        !
        ! Init la densite en representation MDFT, c'est à dire en REAL space, angles d'Euler
        !
        do ipsi=1,npsi
            do iphi=1,nphi
                do itheta=1,ntheta
                    io=grid%indo(itheta,iphi,ipsi)
                    do iz=1,nz
                        do iy=1,ny
                            do ix=1,nx
                                deltarho(ix,iy,iz,itheta,iphi,ipsi)= &
                                    solvent(1)%density(ix,iy,iz,io) - solvent(1)%rho0
                            end do
                        end do
                    end do
                end do
            end do
        end do

        !
        ! On passe tout de suite en projection, dans le repère cartesien
        ! ces projections sont complexes, mais la symetrie hermitienne
        ! permet de ne garder que les mup>=0 ou les mu>=0.
        ! On choisit les mu>=0 (cf doc de Luc)
        !
        call cpu_time(tb)
        call dfftw_plan_dft_r2c_2d(  fft2d%plan, npsi, nphi,  fft2d%in,  fft2d%out, FFTW_EXHAUSTIVE )
        call dfftw_plan_dft_c2r_2d(  ifft2d%plan, npsi, nphi, ifft2d%in, ifft2d%out, FFTW_EXHAUSTIVE )



        DO iz=1,nz
            DO iy=1,ny
                DO ix=1,nx
                    deltarho_p(1:nproj,ix,iy,iz) = euler2proj( deltarho(ix,iy,iz,1:ntheta,1:nphi,1:npsi) )
                    !
                    ! Test if result contains NaN or Inf
                    !
                    if (debug .and. any(deltarho_p(1:nproj,ix,iy,iz) /= deltarho_p(1:nproj,ix,iy,iz)) ) then
                        print*, "Error in euler2proj for ix,iy,iz =", ix, iy, iz
                        print*, "deltarho(ix,iy,iz,1:ntheta,1:nphi,1:npsi) =", deltarho(ix,iy,iz,1:ntheta,1:nphi,1:npsi)
                        print*, "deltarho_p(1:nproj,ix,iy,iz) = ",deltarho_p(1:nproj,ix,iy,iz)
                        error stop "Bug found in euler2proj. NaN or Inf is hiding somewhere."
                    end if
                END DO
            END DO
        END DO
        if (mmax==0) then
            if (complex(deltarho(1,1,1,1,1,1),0._dp) /= deltarho_p(1,1,1,1)) then
                print*, "in energy_cproj et mmax==0, deltarho /= deltarho_p avant meme la FFT3D"
                print*, "deltarho(1,1,1,1,1,1),0._dp) =", complex(deltarho(1,1,1,1,1,1),0._dp)
                print*, "deltarho_p(1,1,1,1)", deltarho_p(1,1,1,1)
                error stop
            end if
        end if


        !
        ! On a les projections sur la grille cartesienne
        ! On veut passer dans l'espace de Fourier pour calculer les convolutions spatiales
        ! On fait donc une FFT 3D pour chacune des projections.
        ! Les projections sont complexes, il s'agit donc d'une FFT3D C2C habituelle : Il n'y a pas de symétrie hermitienne.
        !
        call dfftw_plan_dft_3d(  fft3d%plan, nx, ny, nz, fft3d%in, fft3d%in, FFTW_FORWARD, FFTW_ESTIMATE ) ! in-place
        do ip=1,nproj
            ! deltarho_p in real space (is a complex number)
            fft3d%in = deltarho_p(ip,1:nx,1:ny,1:nz)
            call dfftw_execute( fft3d%plan )
            deltarho_p(ip,1:nx,1:ny,1:nz) = fft3d%in
            ! now deltarhop in fourier space
        end do
        if (mmax==0) then
            if (complex(deltarho(1,1,1,1,1,1),0._dp) /= deltarho_p(1,1,1,1)) then
                print*, "in energy_cproj et mmax=0, deltarho /= deltarho_p apres la FFT3D"
                print*, "deltarho(1,1,1,1,1,1),0._dp) =", complex(deltarho(1,1,1,1,1,1),0._dp)
                print*, "deltarho_p(1,1,1,1)", deltarho_p(1,1,1,1)
                error stop
            end if
        end if


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
        ! See function qproj
        ! And the inverse index that gives vector -q
        ! See function qproj_inv
        !
        tqx = [( qproj(ix,nx,lx), ix=1,nx )]
        tqy = [( qproj(iy,ny,ly), iy=1,ny )]
        tqz = [( qproj(iz,nz,lz), iz=1,nz )]
        tqxi = [( qproj_inv(-tqx(ix),nx,lx)  ,  ix=1,nx  )]
        tqyi = [( qproj_inv(-tqy(iy),ny,ly)  ,  iy=1,ny  )]
        tqzi = [( qproj_inv(-tqz(iz),nz,lz)  ,  iz=1,nz  )]

        !
        ! Read Luc's direct correlation function c^{m,n}_{mu,nu_,chi}(|q|)
        ! projected on generalized spherical harmonics
        ! in the intermolecular frame
        ! normq is norm of q, |q|, that correspond to the index iq in ck(ia,iq)
        !
        !call read_ck_nmax (ck, normq)
        !call read_ck_toutes_nmax( ck, normq)
        call read_ck_nonzero( ck, normq )
        dq=normq(2)

        !
        ! For all vectors q and -q handled simultaneously
        !
        do ix_q=1,nx ! /2+1 IS SUFFICIENT. THIS IS A TEST
            do iy_q=1,ny
                do iz_q=1,nz
                    !
                    ! cartesian coordinates of vector q in lab frame
                    !
                    q = [tqx(ix_q), tqy(iy_q), tqz(iz_q)]

                    ! !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
                    ! block
                    !     real(dp) :: norm2q, cos_theta, sin_theta, theta, phi
                    !     complex(dp) :: mat(0:mmax,-mmax:mmax,-mmax:mmax)
                    !     integer :: m, mup, mu
                    !     norm2q= norm2(q)
                    !     ! print*, "q=",q
                    !     if (norm2q /= 0) then
                    !         cos_theta = q(3)/norm2q
                    !         sin_theta = sqrt(1-cos_theta**2)
                    !         theta = acos(cos_theta)
                    !         phi = angle( q(1)/norm2q , q(2)/norm2q )
                    !     else
                    !         stop "norm2q==0"
                    !     end if
                    !     mat = rotation_matrix_between_complex_spherical_harmonics(q, mmax)
                    !     do m=0,mmax
                    !         do mup=-m,m
                    !             do mu=-m,m
                    !                 print*, mat(m,mup,mu)
                    !                 print*, harm_sph(m,mup,mu,theta)*exp(-imag*mu*phi)
                    !                 print*, harm_sph(m,mup,mu,theta)*exp(imag*mup*phi)
                    !                 print*,
                    !             end do
                    !         end do
                    !     end do
                    !     ! stop "ici"
                    ! end block
                    ! !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

                    !
                    ! Find ix_mq,iy_mq,iz_mq so that q(ix_mq,iy_mq,iz_mq)= -q(ix_q,iy_q,iz_q)
                    !
                    ix_mq = tqxi(ix_q)
                    iy_mq = tqyi(iy_q)
                    iz_mq = tqzi(iz_q)

                    if (ix_mq==ix_q .and. iy_mq==iy_q .and. iz_mq==iz_q) then
                        q_eq_mq=.true.
                    else
                        q_eq_mq=.false.
                    end if

                    !
                    ! On peut vérifier que
                    ! q=[grid%kx(ix_q), grid%ky(iy), grid%kz(iz)]
                    ! q=-[grid%kx(iix_q), grid%ky(iiy), grid%kz(iiz)] ! but for ix_q>nx/2+1 where kx is not allocated
                    !
                    !
                    ! Move all projections to a smaller temporary array: deltarho_p_q (for q) and deltarho_p_mq (for -q)
                    !
                    deltarho_p_q (1:nproj) = deltarho_p(1:nproj,ix_q,iy_q,iz_q)
                    deltarho_p_mq(1:nproj) = deltarho_p(1:nproj,ix_mq,iy_mq,iz_mq)


                    !
                    ! Rotate projections from laboratory frame to molecular frame
                    ! If q=0 then chose the rotation you want: the identity is the simplest option.
                    !
                    if ( any([ix_q, iy_q, iz_q]/=1) ) then  ! case 1,1,1 corresponds to |q|=0. Rotation is arbitrary: we chose the Identity and thus have nothing to do.
                        deltarho_m_mup_mu_q  = zeroc
                        deltarho_m_mup_mu_mq = zeroc
                        DO ip=1,nproj
                            m = im(ip)
                            mup = imup(ip)
                            mu = imu(ip)
                            deltarho_m_mup_mu_q(m,mup,mu) = deltarho_p_q(ip)
                            deltarho_m_mup_mu_mq(m,mup,mu) = deltarho_p_mq(ip)
                        END DO
                        xxx_temp_array_of_m_khi_mu_q = zeroc
                        xxx_temp_array_of_m_khi_mu_mq = zeroc
                        xxx_temp_array_of_m_khi_mu_q = lab2molframe  (deltarho_m_mup_mu_q, q)
                        xxx_temp_array_of_m_khi_mu_mq = lab2molframe (deltarho_m_mup_mu_mq, -q)

                        DO ip=1,nproj
                            m = im(ip)
                            khi = imup(ip)
                            mu = imu(ip)
                            deltarho_p_q(ip) = xxx_temp_array_of_m_khi_mu_q(m,khi,mu)
                            deltarho_p_mq(ip) = xxx_temp_array_of_m_khi_mu_mq(m,khi,mu)
                        END DO
                        ! now deltarho_p_q and deltarho_p_mq are in the molecular frame
                    END IF

                    !
                    ! c^{m,n}_{mu,nu,chi}(|q|) is tabulated for 200 values of |q|.
                    ! Find the tabulated value that is closest to |q|. Its index is iq.
                    ! Note |q| = |-q| so iq is the same for both vectors.
                    !
                    iq = min( int(norm2(q)/dq)+1 ,size(normq))
                    if (iq>size(normq)) then
                        iq=size(normq)
                    else if (iq<=0) then
                        print*, "in energy_cproj, ik <=0"
                        error stop
                    end if

                    !
                    ! Put ck(all projections,iq) to a smaller temp array, ck_q
                    ! c(p,q) => c_local(p)
                    !
                    ck_q(:) = ck(:,iq)

                    !
                    ! Ornstein-Zernike in the molecular frame
                    !

                    gamma_p_q(1:nproj) = zeroc
                    gamma_p_mq(1:nproj) = zeroc

                    DO khi=-mmax,mmax
                        DO m=ABS(khi),mmax
                            DO mu=0,m,molrotsymorder ! not -m,m a cause des symetries ! EST CE QU'IL Y A UNE RAISON POUR QUE LES mu IMPAIRES SOIENT NON NULS ICI ?
                                if (mod(mu,molrotsymorder) /= 0) cycle







                                gamma_m_khi_mu_q= zeroc
                                gamma_m_khi_mu_mq = zeroc
                                do n=abs(khi),mmax
                                    do nu=-n,n,molrotsymorder
                                        if (mod(nu,molrotsymorder) /=0) cycle ! don't threat cases where c is (0,0) for instance all nu even in water (molrotsymorder==2)
                                        print*, "m, n, mu, nu, khi", m, n, mu, nu, khi

                                        ia = ialpha(m,n,mu,nu,khi)
                                        if (ia==errorgenerator) stop "ia = errorgenerator omg"
                                        if (ia==0) stop "ia=0 in energy_cproj"
                                        if (nu<0) then ! no problem with delta rho (n, khi, -nu) since -nu>0. Thus, we apply eq. 1.30 directly
                                            ip = indp(n,khi,-nu)
                                            gamma_m_khi_mu_q  = gamma_m_khi_mu_q  + (-1)**(khi+nu) *ck_q(ia) *deltarho_p_q(ip)
                                            gamma_m_khi_mu_mq = gamma_m_khi_mu_mq + (-1)**(khi+nu) *ck_q(ia) *deltarho_p_mq(ip)
                                        else ! transform delta rho (n, khi, -nu)(q) into conjg( deltarho(n,khi,nu)(-q) )
                                            ip = indp(n,khi,nu)
                                            gamma_m_khi_mu_q  = gamma_m_khi_mu_q  + (-1)**(n) *ck_q(ia) *conjg(deltarho_p_mq(ip))
                                            gamma_m_khi_mu_mq = gamma_m_khi_mu_mq + (-1)**(n) *ck_q(ia) *conjg(deltarho_p_q(ip))
                                        end if
                                    end do
                                end do

                                ip=indp(m,khi,mu)
                                gamma_p_q(ip) = gamma_m_khi_mu_q
                                gamma_p_mq(ip) = gamma_m_khi_mu_mq


                            END DO
                        END DO
                    END DO

                    !
                    ! Rotation from molecular frame to fix laboratory (Fourier) frame
                    !
                    IF ( any([ix_q, iy_q, iz_q]/=1) ) then
                        xxx_temp_array_of_m_khi_mu_q=zeroc
                        xxx_temp_array_of_m_khi_mu_mq=zeroc

                        DO ip=1,nproj
                            m=im(ip)
                            khi=imup(ip)
                            mu=imu(ip)
                            xxx_temp_array_of_m_khi_mu_q(m,khi,mu) = gamma_p_q(ip)
                            xxx_temp_array_of_m_khi_mu_mq(m,khi,mu) = gamma_p_mq(ip)
                        END DO

                        deltarho_m_mup_mu_q = mol2labframe (xxx_temp_array_of_m_khi_mu_q, q)
                        deltarho_m_mup_mu_mq = mol2labframe (xxx_temp_array_of_m_khi_mu_mq, -q)

                        DO ip=1,nproj
                            m = im(ip)
                            mup = imup(ip)
                            mu = imu(ip)
                            gamma_p_q(ip) = deltarho_m_mup_mu_q(m,mup,mu)
                            gamma_p_mq(ip) = deltarho_m_mup_mu_mq(m,mup,mu)
                        END DO
                    END IF

                    gamma_p(1:nproj, ix_q,  iy_q,  iz_q) = gamma_p_q(1:nproj)
                    gamma_p(1:nproj, ix_mq, iy_mq, iz_mq) = gamma_p_mq(1:nproj)

                    if (q_eq_mq .and. any(gamma_p_q/=gamma_p_mq)) then
                        print*, "gamma_p_q", gamma_p_q
                        print*, "gamma_p_mq", gamma_p_mq
                        stop
                    end if

                end do
            end do
        end do


        call cpu_time(t6)

        !
        ! FFT3D from Fourier space to real space
        !
        call dfftw_plan_dft_3d(  fft3d%plan, nx, ny, nz, fft3d%in, fft3d%in, FFTW_BACKWARD, FFTW_MEASURE )
        DO ip= 1, nproj
            fft3d%in = gamma_p(ip,1:nx,1:ny,1:nz)
            call dfftw_execute( fft3d%plan )
            gamma_p(ip,:,:,:) = fft3d%in/real(nx*ny*nz,dp)
        END DO

        ffc=sum(gamma_p(1:nproj,1:nx,1:ny,1:nz)*deltarho_p(1:nproj,1:nx,1:ny,1:nz))/real(nx*ny*nz,dp)
        ! ffc should be real (have no imaginary part)
        if (abs(aimag(ffc))>epsilon(1._dp)) then
            print*, "aimag(ffc)/=0 in energy_cproj"
            print*, "ffc =", ffc
            error stop
        end if
        ff=-thermo%kbT/2._dp*dv*real(ffc,dp)
        print*, "ff from aimag=", ff


        !
        ! Gather projections
        !
        call dfftw_plan_dft_c2r_2d(ifft2d%plan, npsi, nphi, ifft2d%in, ifft2d%out, FFTW_EXHAUSTIVE)
        do iz=1,nz
            do iy=1,ny
                do ix=1,nx
                    gamma(ix,iy,iz,1:ntheta,1:nphi,1:npsi) = proj2euler( gamma_p(1:nproj,ix,iy,iz) )
                end do
            end do
        end do


ff=0._dp
            do ix=1,nx
                do iy=1,ny
                    do iz=1,nz
                        do itheta=1,ntheta
                            do iphi=1,nphi
                                do ipsi=1,npsi
                                    io=grid%indo(itheta,iphi,ipsi)
                                    gamma_o(ix,iy,iz,io)=gamma(ix,iy,iz,itheta,iphi,ipsi)
    ff=ff-kT/2._dp*dv*gamma_o(ix,iy,iz,io)*solvent(1)%density(ix,iy,iz,io)*grid%w(io)
                                end do
                            end do
                        end do
                    end do
                end do
            enddo



        !
        ! Timings
        !
        ! PRINT*,"============"
        ! PRINT*,"TIMINGS (sec)"
        ! PRINT*,"============"
        ! PRINT*,"Initialization               =",tb-tc
        ! PRINT*,"Projection                   =",ta-tb
        ! PRINT*,"Plan the FFT3D               =",t0-ta
        ! PRINT*,"FFT3D                        =",t1-t0
        ! PRINT*,"Rotate from lab to mol frame =",t3-t2
        ! PRINT*,"OZ                           =",t4-t3
        ! PRINT*,"Rotate from mol to lab frame =",t5-t4
        ! PRINT*,"All lab->mol + OZ + mol->lab =",t6-t1
        ! PRINT*,"FFT3D inverse                =",t7-t6
        ! PRINT*,"Gather projections           =",t8-t7
        ! PRINT*,"============"
        ! PRINT*,"TOTAL TIME                   =",t8-tc
        ! PRINT*,"============"
        df=df



        !=============================================================================================
    CONTAINS


        PURE FUNCTION harm_sph( m, mu, mup, beta ) ! Luc's luc72p143
            !
            ! R^m_{mu,mup}(\beta)
            ! beta is the angle in radian
            !
            REAL(dp) :: harm_sph
            INTEGER, intent(in) :: m, mu, mup
            REAL(dp), intent(in) :: beta
            REAL(dp) :: beta0, x, cc, ss, pm, pm1, pm2
            INTEGER :: mu0, l, it
            if( abs(mu)>m .or. abs(mup)>m ) then
                harm_sph = 0
                return
            end if
            if( m==0 ) then
                harm_sph = 1
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
                harm_sph = x*(-1.)**mu0*sqrt(fact(m-mu0)/fact(m+mu0))&
                *fact(2*mu0)/(2.**mu0*fact(mu0))*sin(beta0)**mu0*pm
            else                        ! donc mu et mup non nuls, utiliser betement la formule de wigner
                harm_sph = 0
                cc = cos(0.5*beta)
                ss = sin(0.5*beta)
                do it = max(0,mu-mup),min(m+mu,m-mup)
                    harm_sph = harm_sph+(-1.)**it/(fact(m+mu-it)*fact(m-mup-it)*fact(it)*&
                    fact(it-mu+mup))*cc**(2*m+mu-mup-2*it)*ss**(2*it-mu+mup)
                end do
                harm_sph = sqrt(fact(m+mu)*fact(m-mu)*fact(m+mup)*fact(m-mup))*harm_sph
            end if
        END FUNCTION harm_sph
        !
        !
        !
        PURE SUBROUTINE gauss_legendre( n, x, w, exitstatus) ! copy paste from Luc's subroutine Luc74p85
            !
            ! Returns the n roots (x) and associated weights(w) of a gauss legendre quadrature of order n
            ! The roots are the Cos(theta) so that if you need theta, don't forget to acos(x)
            !
            IMPLICIT NONE
            INTEGER, intent(in) :: n
            REAL(dp), intent(out) :: x(n), w(n)
            INTEGER, optional, intent(out) :: exitstatus
            INTEGER :: m, i, j
            REAL(dp), PARAMETER :: pi=acos(-1._dp)
            REAL(dp) :: xi, p1, p2, p3, pp, deltaxi
            exitstatus = 1
            if( n <= 0 ) then
                exitstatus = -1
                return
            else
                m = (n+1)/2
                do i = 1,m          ! on s'interesse au ième zero du polynome pn(x) de legendre
                    xi = cos(pi*(i-0.25)/(n+0.5))     ! estimation de départ qu'on va raffiner par nr
                    deltaxi = 1
                    do while (abs(deltaxi)>1.d-13)
                        p1 = 1.
                        p2 = 0.
                        do j = 1,n
                            p3 = p2
                            p2 = p1
                            p1 = ((2*j-1.)*xi*p2-(j-1.)*p3)/j         ! relation de récurrence entre les pj
                        end do
                        pp = n*(xi*p1-p2)/(xi**2-1.)                   ! donne pn' en fonction de pn et pn-1
                        deltaxi = -p1/pp                                  ! nr
                        xi = xi + deltaxi
                    end do
                    x(i) = xi
                    w(i) = 1./((1.-xi**2)*pp**2)                  ! poids normalise a 1
                    x(n+1-i) = -xi
                    w(n+1-i) = w(i)
                end do
            end if
        END SUBROUTINE gauss_legendre
        !
        !
        !
        PURE FUNCTION trimed(x)
            !
            ! If some REAL value is smaller than some epsilon, then give it 0 value
            ! This is usefull for printing stuff where you don't want 1E-19 to be printed everywhere.
            !
            implicit none
            REAL(dp) :: trimed
            REAL(dp), intent(in) :: x
            REAL(dp), PARAMETER :: zero = 0._dp
            if( abs(x) <= epsilon(x) ) trimed = zero
        END FUNCTION trimed
        !
        !
        !
        PURE FUNCTION cross_product(a,b)
            !
            ! Cross product of two vectors
            !
            IMPLICIT NONE
            REAL(dp) :: cross_product(3)
            REAL(dp), INTENT(IN) :: a(3), b(3)
            cross_product(1) = a(2)*b(3)-a(3)*b(2)
            cross_product(2) = a(3)*b(1)-a(1)*b(3)
            cross_product(3) = a(1)*b(2)-a(2)*b(1)
        END FUNCTION cross_product
        !
        !
        !
        PURE FUNCTION rotation_matrix_from_lab_to_q_frame(q) RESULT(R)
            !
            ! Given the coordinates q(3) in a frame (let's call it lab frame),
            ! this function produces the rotation matrix that
            ! transforms any coordinates in lab frame to some frame
            ! whose (0,0,1) vector is aligned with q.
            ! In other words, transform a vector in lab frame to a vector in intermolecular (or q) frame.
            !
            IMPLICIT NONE
            REAL(dp), INTENT(IN) :: q(3)
            REAL(dp) :: R1(3), R2(3), R3(3), R(3,3)
            R3 = q/NORM2(q) ! The rotation matrix R will induce z to be aligned to vector q
            IF( ABS(R3(1)) < 0.99 ) THEN
                R1 = cross_product( REAL([1,0,0],dp) , R3 ) ! we dont care about R2 and R1, as long as they are orthogonal to R3.
            ELSE
                R1 = cross_product( REAL([0,1,0],dp) , R3 )
            END IF
            R1 = R1/NORM2(R1) ! normalization
            R2 = cross_product( R3, R1 )
            R(1,:) = R1
            R(2,:) = R2
            R(3,:) = R3
        END FUNCTION
        !
        !
        FUNCTION rotation_matrix_between_complex_spherical_harmonics(q,lmax) RESULT (R)
            !
            ! Given a column vector of REAL values q of dimension 3,
            ! this function returns the rotation matrix between complex spherical harmonics.
            ! The algorithm is from Choi et al., J. Chem. Phys. 111, 8825 (1999)
            ! http://dx.doi.org/10.1063/1.480229
            ! see eq 6.1 and following and 7.3 and following
            ! We use Choi's notation
            !
            ! Luc Belloni wrote the original fortran code (luc85p17)
            ! Maximilien Levesque rewrote it from scratch thanks to Luc's explanations
            ! Conditions for m==-l and m==l by ML
            ! All bugs due to ML
            !
            IMPLICIT NONE
            !
            INTEGER, INTENT(IN) :: lmax
            REAL(dp), INTENT(IN) :: q(3)
            COMPLEX(dp) :: R(0:lmax,-lmax:lmax,-lmax:lmax)
            REAL(dp), DIMENSION(0:lmax,-lmax:lmax,-lmax:lmax) :: fi, gi, a, b, c, d
            REAL(dp) :: Rot(3,3)
            REAL(dp), PARAMETER :: tsqrt2=SQRT(2._dp)
            REAL(dp) :: tsqrt(-1:2*lmax+1), adenom, cdenom, anumerator, bnumerator
            INTEGER :: l,l1,m,m1,m1min
            !
            ! Init R, fi and gi to zero
            !
            fi = 0
            gi = 0
            R = CMPLX(0,0)
            !
            ! l = 0
            !
            fi(0,0,0) = 1
            gi(0,0,0) = 0
            IF( lmax == 0 ) THEN
                R = CMPLX(1,0)
                RETURN
            ELSE
                !
                ! We will need the rotation matrix from lab to q frame
                !
                Rot = rotation_matrix_from_lab_to_q_frame(q)
            END IF
            !
            ! l = 1
            !
            fi(1,-1,-1) = (Rot(2,2)+Rot(1,1))/2.
            fi(1,-1, 0) = Rot(1,3)/tsqrt2
            fi(1,-1, 1) = (Rot(2,2)-Rot(1,1))/2.
            fi(1, 0,-1) = Rot(3,1)/tsqrt2
            fi(1, 0, 0) = Rot(3,3)
            fi(1, 0, 1) = -Rot(3,1)/tsqrt2
            fi(1, 1,-1) = (Rot(2,2)-Rot(1,1))/2.
            fi(1, 1, 0) = -Rot(1,3)/tsqrt2
            fi(1, 1, 1) = (Rot(2,2)+Rot(1,1))/2.
            gi(1,-1,-1) = (Rot(2,1)-Rot(1,2))/2.
            gi(1,-1, 0) = Rot(2,3)/tsqrt2
            gi(1,-1, 1) = (-Rot(2,1)-Rot(1,2))/2.
            gi(1, 0,-1) = -Rot(3,2)/tsqrt2
            gi(1, 0, 0) = 0
            gi(1, 0, 1) = -Rot(3,2)/tsqrt2
            gi(1, 1,-1) = (Rot(2,1)+Rot(1,2))/2.
            gi(1, 1, 0) = Rot(2,3)/tsqrt2
            gi(1, 1, 1) = (Rot(1,2)-Rot(2,1))/2.
            IF( lmax == 1 ) THEN
                R = CMPLX(fi,gi,dp)
                RETURN
            ELSE
                !
                ! Tabulate sqrt
                !
                tsqrt(-1)=0 ! to avoid problems later
                tsqrt(0)=0
                tsqrt(1)=1
                tsqrt(2)=tsqrt2
                DO l=3,2*lmax+1
                    tsqrt(l) = SQRT(REAL(l,dp))
                END DO
                !
                ! Coefficients needed for the recursion
                ! see eq 6.1 to 6.5 of Choi et al.
                !
                DO l=1,lmax
                    a(l,-l,:)=0
                    b(l,-l,:)=0
                    c(l,-l,:)=0
                    d(l,-l,:)=0
                    DO m=-l+1,l
                        anumerator=tsqrt(l+m)*tsqrt(l-m)
                        bnumerator=tsqrt(l+m)*tsqrt(l+m-1)
                        DO m1=-l+1,l-1
                            adenom = 1/(tsqrt(l+m1)*tsqrt(l-m1))
                            a(l,m,m1)=anumerator*adenom
                            b(l,m,m1)=bnumerator*adenom/tsqrt2
                        END DO
                        m1=l
                        cdenom=1/(tsqrt(l+m1)*tsqrt(l+m1-1))
                        c(l,m,m1)=tsqrt2*anumerator*cdenom
                        d(l,m,m1)=       bnumerator*cdenom
                    END DO
                END DO
            END IF
            !
            ! l > 1    Use recursion
            !
            DO l=2,lmax
                l1=l-1
                DO m=-l,l
                    m1min=0
                    IF(m>0) m1min=1
                    DO m1=m1min,l-1
                        if( m==-l) then
                            fi(l,m,m1)=b(l,-m,m1)*(fi(1,-1,0)*fi(l1,m+1,m1)-gi(1,-1,0)*gi(l1,m+1,m1))
                            gi(l,m,m1)=b(l,-m,m1)*(fi(1,-1,0)*gi(l1,m+1,m1)+gi(1,-1,0)*fi(l1,m+1,m1))
                        else if( m==l) then
                            fi(l,m,m1)=b(l,m,m1)*(fi(1, 1,0)*fi(l1,m-1,m1)-gi(1, 1,0)*gi(l1,m-1,m1))
                            gi(l,m,m1)=b(l,m,m1)*(fi(1, 1,0)*gi(l1,m-1,m1)+gi(1, 1,0)*fi(l1,m-1,m1))
                        else
                            fi(l,m,m1)=a(l, m,m1)*(fi(1, 0,0)*fi(l1,m  ,m1))+   &
                            b(l, m,m1)*(fi(1, 1,0)*fi(l1,m-1,m1)-gi(1, 1,0)*gi(l1,m-1,m1))+   &
                            b(l,-m,m1)*(fi(1,-1,0)*fi(l1,m+1,m1)-gi(1,-1,0)*gi(l1,m+1,m1))
                            gi(l,m,m1)=a(l, m,m1)*(fi(1, 0,0)*gi(l1,m  ,m1))+   &
                            b(l, m,m1)*(fi(1, 1,0)*gi(l1,m-1,m1)+gi(1, 1,0)*fi(l1,m-1,m1))+   &
                            b(l,-m,m1)*(fi(1,-1,0)*gi(l1,m+1,m1)+gi(1,-1,0)*fi(l1,m+1,m1))
                        end if
                        fi(l,-m,-m1)=(-1)**(m+m1)*fi(l,m,m1)
                        gi(l,-m,-m1)=-(-1)**(m+m1)*gi(l,m,m1)
                    END DO
                    m1=l
                    if( m==-l) then
                        fi(l,m,m1)=d(l,-m,m1)*(fi(1,-1,+1)*fi(l1,m+1,m1-1)-gi(1,-1,+1)*gi(l1,m+1,m1-1))
                        gi(l,m,m1)=d(l,-m,m1)*(fi(1,-1,+1)*gi(l1,m+1,m1-1)+gi(1,-1,+1)*fi(l1,m+1,m1-1))
                    else if( m==l) then
                        fi(l,m,m1)=d(l,m,m1)*(fi(1,+1,+1)*fi(l1,m-1,m1-1)-gi(1,+1,+1)*gi(l1,m-1,m1-1))
                        gi(l,m,m1)=d(l,m,m1)*(fi(1,+1,+1)*gi(l1,m-1,m1-1)+gi(1,+1,+1)*fi(l1,m-1,m1-1))
                    else
                        fi(l,m,m1)=c(l,m,m1)*(fi(1,0,+1)*fi(l1,m,m1-1)-gi(1,0,+1)*gi(l1,m,m1-1))+         &
                        d(l,m,m1)*(fi(1,+1,+1)*fi(l1,m-1,m1-1)-gi(1,+1,+1)*gi(l1,m-1,m1-1))+   &
                        d(l,-m,m1)*(fi(1,-1,+1)*fi(l1,m+1,m1-1)-gi(1,-1,+1)*gi(l1,m+1,m1-1))
                        gi(l,m,m1)=c(l,m,m1)*(fi(1,0,+1)*gi(l1,m,m1-1)+gi(1,0,+1)*fi(l1,m,m1-1))+         &
                        d(l,m,m1)*(fi(1,+1,+1)*gi(l1,m-1,m1-1)+gi(1,+1,+1)*fi(l1,m-1,m1-1))+   &
                        d(l,-m,m1)*(fi(1,-1,+1)*gi(l1,m+1,m1-1)+gi(1,-1,+1)*fi(l1,m+1,m1-1))
                    end if
                    fi(l,-m,-m1)=(-1)**(m+m1)*fi(l,m,m1)
                    gi(l,-m,-m1)=-(-1)**(m+m1)*gi(l,m,m1)
                END DO
            END DO
            !
            ! fi and gi are the REAL and imaginary part of our rotation matrix
            !
            R = CMPLX(fi,gi,dp)
        END FUNCTION rotation_matrix_between_complex_spherical_harmonics
        !
        PURE FUNCTION qproj( i, imax, length)
            !
            ! Say you have grid nodes every 2pi/lx distance units (unit given by the unit of lx)
            ! and you have imax grid nodes, indices of which range from 1 to imax
            ! Remember you are in periodic boundary conditions.
            ! What is the coordinate of point index i in FFTW's way of representing stuff (ie with q=0 for index 1)?
            ! That is what compute qproj.
            ! For my part I was used to thinking about q vectors from negative values to positives ones. Here it starts with 0.
            ! See FFTW documentation.
            !
            ! i must be in [1,imax]
            ! imax must be >= 1
            ! lx must be > 0.
            !
            IMPLICIT NONE
            REAL(dp) :: qproj
            REAL(dp), PARAMETER :: twopi=2*ACOS(-1._dp)
            INTEGER, INTENT(IN) :: i, imax ! index the node and number of nodes in the direction you wish
            ! If you use REAL to complex transforms, and thus have in the first direction nx/2+1 nodes, you should still pass nx to this
            ! routine, not nx/2+1
            REAL(dp), INTENT(IN) :: length ! size of the supercell
            IF( i == 1 ) THEN
                qproj = 0._dp
            ELSE IF( i <= imax/2 ) THEN
                qproj = twopi*REAL(i-1,dp)/length
            ELSE
                qproj = twopi*REAL(i-1-imax,dp)/length
            END IF
        END FUNCTION qproj
        !
        PURE FUNCTION qproj_inv( qi, imax, length)
            !
            ! You know the length of vector q in one direction, and you want to know what index in the table it comes from
            ! we thus have    qproj_inv( qproj(i,imax,length) , imax, length) = i
            !
            IMPLICIT NONE
            INTEGER :: qproj_inv
            REAL(dp), INTENT(IN) :: qi ! length of vector q in the direction you want
            REAL(dp), INTENT(IN) :: length ! length of the supercell, in the direction you want
            INTEGER, INTENT(IN)  :: imax ! how many nodes in this direction
            REAL(dp), PARAMETER  :: twopi=2*ACOS(-1._dp)
            IF( qi < 0 ) THEN
                qproj_inv = NINT( qi*length/twopi +1 +imax )
            ELSE
                qproj_inv = NINT( qi*length/twopi +1 )
            END IF
        END FUNCTION qproj_inv

        SUBROUTINE read_ck_nonzero ( ck, normq )
            implicit none
            integer :: i, nproj_max, iq, m, n, mu, nu, khi, ia
            INTEGER, PARAMETER :: nalpha=SUM(  [([([([([(1,nu=-n,n)],mu=-m,m)],&
            n=abs(khi),mmax)],m=abs(khi),mmax)],khi=-mmax,mmax)]  )
            !INTEGER, PARAMETER :: nalpha=(1+mmax)*(1+2*mmax)*(3+2*mmax)*(5+4*mmax*(2+mmax))/15 ! number of projections for c^{m n}_{mu nu,khi}(q)
            complex(dp), allocatable, dimension(:,:) :: ck_max
            real(dp), intent(out) :: normq(nq)
            complex(dp), intent(out) :: ck(nalpha,nq)
            character(3) :: somechar
            integer :: ufile, ios
            integer, allocatable, dimension(:) :: alp_max ,m_max, n_max, mu_max, nu_max, khi_max
            if (myck%isok) return
            select case (mmax)
            case (0)
                nproj_max=1
                open(newunit=ufile, file="ck_nonzero_nmax0_ml", iostat=ios, status="old", action="read")
                if ( ios /= 0 ) stop "Error opening file ck_nonzero_nmax0_ml"
            case (1)
                nproj_max = 6
                open(newunit=ufile, file="ck_nonzero_nmax1_ml", iostat=ios, status="old", action="read")
                if ( ios /= 0 ) stop "Error opening file ck_nonzero_nmax1_ml"
            case (2)
                nproj_max = 75
                open(newunit=ufile, file="ck_nonzero_nmax2_ml", iostat=ios, status="old", action="read")
                if ( ios /= 0 ) stop "Error opening file ck_nonzero_nmax2_ml"
            case (3)
                nproj_max = 252
                open(newunit=ufile, file="ck_nonzero_nmax3_ml", iostat=ios, status="old", action="read")
                if ( ios /= 0 ) stop "Error opening file ck_nonzero_nmax3_ml"
            case (4)
                nproj_max = 877
                open(newunit=ufile, file="ck_nonzero_nmax4_ml", iostat=ios, status="old", action="read")
                if ( ios /= 0 ) stop "Error opening file ck_nonzero_nmax4_ml"
            case (5)
                nproj_max = 2002
                open(newunit=ufile, file="ck_nonzero_nmax5_ml", iostat=ios, status="old", action="read")
                if ( ios /= 0 ) stop "Error opening file ck_nonzero_nmax5_ml"
            case default
                error stop "In energy_cproj, read_ck_nonzero valid only for mmax=0 to 5"
            end select


            allocate( ck_max(nproj_max,nq) ,source=zeroc)
            allocate( alp_max(nproj_max) )
            allocate( m_max(nproj_max) )
            allocate( n_max(nproj_max) )
            allocate( mu_max(nproj_max) )
            allocate( nu_max(nproj_max) )
            allocate( khi_max(nproj_max) )
            !
            ! Skip 10 lines of comments
            !
            do i=1,10
                read(ufile,*)
            end do

            read(ufile,*) somechar, alp_max
            read(ufile,*) somechar, m_max
            read(ufile,*) somechar, n_max
            read(ufile,*) somechar, mu_max
            read(ufile,*) somechar, nu_max
            read(ufile,*) somechar, khi_max
            read(ufile,*)

            do iq=1,nq
                read(ufile,*) normq(iq), ck_max(:,iq)
            end do

            ck = (0,0)
            do i=1,nproj_max
                m = m_max(i)
                n = n_max(i)
                mu = mu_max(i)
                nu = nu_max(i)
                khi = khi_max(i)
                ia = alp_max(i)
                ialpha(m,n,mu,nu,khi) = ia
                ck(ia,1:nq) = ck_max(i,1:nq)
            end do

            if (allocated(ck_max)) deallocate(ck_max, stat=ios)
            if (ios /= 0) print *, "ck_max: Deallocation request denied"
            if (allocated(m_max)) deallocate(m_max, stat=ios)
            if (ios /= 0) print *, "m_max: Deallocation request denied"
            if (allocated(n_max)) deallocate(n_max, stat=ios)
            if (ios /= 0) print *, "n_max: Deallocation request denied"
            if (allocated(mu_max)) deallocate(mu_max, stat=ios)
            if (ios /= 0) print *, "mu_max: Deallocation request denied"
            if (allocated(nu_max)) deallocate(nu_max, stat=ios)
            if (ios /= 0) print *, "nu_max: Deallocation request denied"
            if (allocated(khi_max)) deallocate(khi_max, stat=ios)
            if (ios /= 0) print *, "khi_max: Deallocation request denied"
            if (allocated(alp_max)) deallocate(alp_max, stat=ios)
            if (ios /= 0) print *, "alp_max: Deallocation request denied"

            myck%isok=.true.

            close(ufile)
        end subroutine read_ck_nonzero


        FUNCTION euler2proj (eulerangles) RESULT (projections)
            IMPLICIT NONE
            REAL(dp), INTENT(IN) :: eulerangles(1:ntheta,1:nphi,1:npsi)
            COMPLEX(dp) :: proj_theta(1:ntheta,0:mmax/molrotsymorder,-mmax:mmax) ! itheta,mu,mup    note we changed mu(psi) to the 2nd position
            COMPLEX(dp) :: projections(1:nproj)
            INTEGER :: itheta, iphi, ipsi, m, mup, mu, ip
            complex(dp), parameter :: zeroc=complex(0._dp,0._dp)
            projections = zeroc
            proj_theta = zeroc
            DO itheta=1,ntheta
                DO iphi=1,nphi
                    DO ipsi=1,npsi
                        fft2d%in(ipsi,iphi) = eulerangles(itheta,iphi,ipsi)
                    END DO
                END DO
                call dfftw_execute( fft2d%plan )
                proj_theta(itheta,0:mmax/molrotsymorder,0:mmax)   = CONJG( fft2d%out(:,1:mmax+1) )
                proj_theta(itheta,0:mmax/molrotsymorder,-mmax:-1) = CONJG( fft2d%out(:,mmax+2:) )
            END DO
            DO m=0,mmax
                DO mup=-m,m
                    DO mu=0,m,molrotsymorder
                        ip = indp( m,mup,mu )
                        projections(ip)= SUM(proj_theta(:,mu,mup)*tharm_sph(:,ip)*wtheta(:))*fm(m)
                    END DO
                END DO
            END DO
        END FUNCTION euler2proj
        !
        FUNCTION proj2euler (projections) RESULT (eulerangles)
            IMPLICIT NONE
            REAL(dp) :: eulerangles(1:ntheta,1:nphi,1:npsi)
            COMPLEX(dp), INTENT(IN) :: projections(1:nproj)
            COMPLEX(dp) :: proj_theta(1:ntheta,0:mmax/molrotsymorder,-mmax:mmax)
            INTEGER :: m, mup, mu, ip, itheta, iphi, ipsi
            proj_theta = CMPLX(0,0)
            DO itheta=1,ntheta
                DO mup=-mmax,mmax
                    DO mu=0,mmax/molrotsymorder
                        DO m= MAX(ABS(mup),ABS(mu)), mmax
                            ip=indp(m,mup,mu)
                            proj_theta(itheta,mu,mup) = proj_theta(itheta,mu,mup)&
                            +projections(ip)*tharm_sph(itheta,ip)*fm(m)
                        END DO
                    END DO
                END DO
            END DO
            DO itheta=1,ntheta
                ifft2d%in(:,1:mmax+1) = CONJG( proj_theta(itheta,:,0:mmax)   )
                ifft2d%in(:,mmax+2:)  = CONJG( proj_theta(itheta,:,-mmax:-1) )
                call dfftw_execute( ifft2d%plan )
                DO iphi=1,nphi
                    DO ipsi=1,npsi
                        eulerangles(itheta,iphi,ipsi) = ifft2d%out(ipsi,iphi) *(nphi*npsi)
                    END DO
                END DO
            END DO
        END FUNCTION proj2euler

        FUNCTION lab2molframe (datalabframe, q) RESULT (datamolframe)
            IMPLICIT NONE
            COMPLEX(dp), INTENT(IN) :: datalabframe(0:mmax,-mmax:mmax,-mmax:mmax)
            COMPLEX(dp) :: datamolframe(0:mmax,-mmax:mmax,-mmax:mmax)
            REAL(dp), INTENT(IN) :: q(1:3)
            INTEGER :: m, mup, mu, khi
            COMPLEX(dp) :: GSHRot(0:mmax,-mmax:mmax,-mmax:mmax)
            GSHRot = rotation_matrix_between_complex_spherical_harmonics(q,mmax)
            datamolframe = CMPLX(0,0)
            DO m=0,mmax
                DO khi=-m,m
                    DO mu=0,m
                        !
                        DO mup=-m,m
                            datamolframe(m,khi,mu) = datamolframe(m,khi,mu)&
                            + datalabframe(m,mup,mu) *GSHRot(m,mup,khi)
                        END DO
                        !
                    END DO
                END DO
            END DO
        END FUNCTION
        !
        FUNCTION mol2labframe (datamolframe, q) RESULT (datalabframe)
            IMPLICIT NONE
            COMPLEX(dp) :: GSHRot(0:mmax,-mmax:mmax,-mmax:mmax)
            COMPLEX(dp) :: datalabframe(0:mmax,-mmax:mmax,-mmax:mmax)
            COMPLEX(dp), INTENT(IN) :: datamolframe(0:mmax,-mmax:mmax,-mmax:mmax)
            REAL(dp), INTENT(IN) :: q(1:3)
            INTEGER :: m, mup, mu, khi
            GSHRot = CONJG(rotation_matrix_between_complex_spherical_harmonics(q,mmax))
            datalabframe = CMPLX(0,0)
            DO m=0,mmax
                DO mup=-m,m
                    DO mu=0,m
                        !
                        DO khi=-m,m
                            datalabframe(m,mup,mu) = datalabframe(m,mup,mu) &
                            + datamolframe(m,khi,mu) *GSHRot(m,mup,khi)
                        END DO
                        !
                    END DO
                END DO
            END DO
        END FUNCTION mol2labframe

        PURE FUNCTION  angle(x,y)
            IMPLICIT NONE
            REAL(dp) :: angle
            REAL(dp), INTENT(IN) :: x, y
            REAL(dp) :: xx, r
            r = SQRT(x**2+y**2)
            xx = x/r
            IF (y >= 0) THEN
                angle = ACOS(xx)
            ELSE
                angle = 2*pi -ACOS(xx)
            END IF
        END FUNCTION
        !
    end subroutine energy_cproj
end module module_energy_cproj






























! FUNCTION projections_threecolumns2onecolumn (data3c) RESULT (data1c)
!     IMPLICIT NONE
!     COMPLEX(dp), INTENT(IN) :: data3c(0:mmax,-mmax:mmax,-mmax:mmax)
!     COMPLEX(dp) :: data1c(nproj)
!     INTEGER :: m, mup, mu, ip
!     DO m=0,mmax
!         DO mup=-m,m
!             DO mu=0,m,molrotsymorder
!                 ip=indp(m,mup,mu)
!                 data1c(ip)=data3c(m,mup,mu)
!             END DO
!         END DO
!     END DO
! END FUNCTION
! !
! !
! !
! FUNCTION projections_onecolumn2threecolumns (data1c) RESULT (data3c)
!     IMPLICIT NONE
!     COMPLEX(dp), INTENT(IN) :: data1c(nproj)
!     COMPLEX(dp) :: data3c(0:mmax,-mmax:mmax,-mmax:mmax)
!     INTEGER :: m, mup, mu, ip
!     DO ip=1,nproj
!         m=im(ip)
!         mup=imup(ip)
!         mu=imu(ip)
!         data3c(m,mup,mu)=data1c(ip)
!     END DO
! END FUNCTION










! SUBROUTINE read_ck_nmax( ck, normq )
!     implicit none
!     integer :: m_, n_, mu_, nu_, khi_
!     CHARACTER(3) :: mychar
!     integer, parameter :: tnalpha(0:5)=[1,4,27,79,250,549]
!     integer, parameter :: nalpha=tnalpha(mmax)
!     integer, dimension(nalpha) :: alp, m, n, mu, nu, khi
!     COMPLEX(dp), INTENT(OUT) :: ck(nalpha,nq)
!     REAL(dp), INTENT(OUT) :: normq(nq)
!     REAL(dp) :: fi(nalpha*2)
!     INTEGER :: i, iq, p, ia, ia2, ios
!     INTEGER, PARAMETER :: iounit=11
!     character(len=8), parameter :: tfilename(0:5)=["ck_nmax0","ck_nmax1","ck_nmax2","ck_nmax3","ck_nmax4","ck_nmax5"]
!     character(len=8), parameter :: filename=tfilename(mmax)
!     !
!     ! Open the file containing all the projections of the direct correlation function in the molecular frame
!     !
!     open(unit=iounit, file=filename, iostat=ios, status="old", action="read")
!     if ( ios /= 0 ) then
!         print*, "Error opening file ", filename
!         error stop "Stopping with error"
!     end if
!
!     DO i=1,10            !
!         READ(iounit,*)   ! Skip 10 lines of comments
!     END DO               !
!
!     READ(iounit,*) mychar, alp(:)
!     READ(iounit,*) mychar, m(:)
!     READ(iounit,*) mychar, n(:)
!     READ(iounit,*) mychar, mu(:)
!     READ(iounit,*) mychar, nu(:)
!     READ(iounit,*) mychar, khi(:)
!     READ(iounit,*)
!
!     print*,
!     PRINT*,"Non-zero projections of the direct correlation function in molecular frame c^{m n}_{mu nu, khi}:"
!     PRINT*,"        index         m           n           mu          nu          khi"
!     PRINT*,"        -----        ---         ---          --          --          ---"
!     DO ia=1,nalpha
!         print*, alp(ia), m(ia), n(ia), mu(ia), nu(ia), khi(ia)
!     END DO
!     PRINT*
!
!     ialpha=errorgenerator
!     DO ia=1,nalpha
!         ialpha( m(ia), n(ia), mu(ia), nu(ia), khi(ia) ) = alp(ia)
!     END DO
!
!     DO iq=1,nq
!         READ(iounit,*) normq(iq), fi
!         DO ia=1,nalpha
!             ck(ia,iq)= CMPLX(fi(2*ia-1),fi(2*ia),dp)
!         END DO
!     END DO
!
!     !
!     ! Close the file containing all the projections of the direct correlation function in the molecular frame
!     !
!     close(unit=iounit, iostat=ios, status="keep")
!     if ( ios /= 0 ) then
!         print*, "Error closing file ", filename, "in read_ck_nmax"
!         error stop "Stopping with error"
!     end if
! END SUBROUTINE read_ck_nmax











! SUBROUTINE read_ck_toutes_nmax( ck , normq )
!     IMPLICIT NONE
!     INTEGER :: m_, n_, mu_, nu_, khi_
!     INTEGER, PARAMETER :: nalpha=SUM( &
!     [([([([([(1,nu_=-n_,n_)],mu_=-m_,m_)],n_=abs(khi_),mmax)],&
!     m_=abs(khi_),mmax)],khi_=-mmax,mmax)]  )
!     !INTEGER, PARAMETER :: nalpha=(1+mmax)*(1+2*mmax)*(3+2*mmax)*(5+4*mmax*(2+mmax))/15 ! number of projections for c^{m n}_{mu nu,khi}(q)
!     CHARACTER(3) :: mychar
!     INTEGER, DIMENSION(nalpha) :: alp, m, n, mu, nu, khi
!     COMPLEX(dp), INTENT(OUT) :: ck(nalpha,nq)
!     REAL(dp), INTENT(OUT) :: normq(nq)
!     REAL(dp) :: fi(nalpha*2)
!     INTEGER :: i, iq, p, ia, ia2
!     SELECT CASE (mmax)
!     CASE (1)
!         OPEN(11, file="ck_toutes_nmax1")
!     CASE (2)
!         OPEN(11, file="ck_toutes_nmax2")
!     CASE (3)
!         OPEN(11, file="ck_toutes_nmax3")
!     CASE DEFAULT
!         ERROR STOP "I don't ck_toutes_nmax file with all projections for nmax = 0 or nmax > 3 yet"
!     END SELECT
!     ! skip 10 lignes of comments
!     DO i=1,10
!         READ(11,*)
!     END DO
!     ! then read useful stuff
!     alp=0
!     m=0
!     n=0
!     mu=0
!     nu=0
!     khi=0
!     ! PRINT*,"nalpha=",nalpha
!     READ(11,*) mychar, alp(:)
!     READ(11,*) mychar, m(:)
!     READ(11,*) mychar, n(:)
!     READ(11,*) mychar, mu(:)
!     READ(11,*) mychar, nu(:)
!     READ(11,*) mychar, khi(:)
!     READ(11,*)
!     ! PRINT*,int(alp,1)
!     ! PRINT*,int(m,1)
!     ! print*,int(n,1)
!     ! print*,int(mu,1)
!     ! print*,int(nu,1)
!     ! print*,int(khi,1)
!     ialpha=errorgenerator
!     DO ia=1,nalpha
!         ialpha( m(ia), n(ia), mu(ia), nu(ia), khi(ia) ) = alp(ia)
!     END DO
!     DO iq=1,nq
!         READ(11,*) normq(iq), fi
!         DO ia=1,nalpha
!             ck(ia,iq)= CMPLX(fi(2*ia-1),fi(2*ia),dp)
!         END DO
!     END DO
!     !
!     !
!     !
!     ! print*,"== TEST READING =="
!     ! do iq=3,3
!     !     do ia=1,nalpha
!     !         ia2=ialpha( m(ia), n(ia), -mu(ia), -nu(ia), khi(ia) )
!     !         block
!     !             real :: diff
!     !             diff = abs(  ck(ia2,iq)  -   (-1)**(m(ia) + n(ia) + mu(ia) + nu(ia))*CONJG(ck(ia,iq))    )
!     !             if( diff == 0 ) then
!     !                 print*, ia, ia2, m(ia), n(ia), mu(ia), nu(ia), khi(ia)
!     !             else
!     !                 print*, ia, ia2, m(ia), n(ia), mu(ia), nu(ia), khi(ia), " #", ck(ia2,iq), (-1)**(m(ia) + n(ia) + mu(ia) + nu(ia))*CONJG(ck(ia,iq))
!     !             end if
!     !         end block
!     !     !if( abs(  ck(ia2,iq)  -   (-1)**(m(ia) + n(ia) + mu(ia) + nu(ia))*CONJG(ck(ia,iq))      ) /= 0) print*,"ERROR FOR ia=",ia
!     ! !    if(abs(ck(ia,iq))/=0) print*,int(iq,1),ia,int([m(ia),n(ia),mu(ia),nu(ia),khi(ia)],1), ck(ia,iq)
!     !     end do
!     ! end do
!     ! print*,"== END TEST READING =="
!     ! CECI EST UN TEST DE LA ROUTINE DE GENERATION
!     do ia=1,nalpha
!         print*,m(ia), n(ia), mu(ia), nu(ia), khi(ia), ck(ia,2)
!     end do
!     ! CECI EST LA FIN DU TEST DE LA ROUTINE DE GENERATION
! END SUBROUTINE read_ck_toutes_nmax

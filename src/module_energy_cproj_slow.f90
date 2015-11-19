module module_energy_cproj_slow

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
    ! Arrays to translate index of vector to tuple of projections, and inverse
    ! alpha(m,mup,mu), m(alpha), mup(alpha), mu(alpha)
    !
    type :: p3_type
        integer :: np ! how many projections
        integer, allocatable :: p(:,:,:) ! p of m, mup, mu
        integer, allocatable :: m(:) ! m of p
        integer, allocatable :: mup(:) ! mu' of p
        integer, allocatable :: mu(:) ! mu of p
        real(dp), allocatable :: tharm_sph(:,:) ! Associated Legendre function of np and theta
        complex(dp), allocatable :: drho(:,:,:,:) ! delta rho en projections
        complex(dp), allocatable :: gam (:,:,:,:) ! gamma en projections
        complex(dp), allocatable :: drho_q(:) ! drho for a given q
        complex(dp), allocatable :: gam_q(:) ! gam for a given q
    end type p3_type
    type (p3_type) :: p3

    type :: o3_type
        complex(dp), allocatable :: gam(:,:,:,:) ! gamma(orientation,x,y,z )
        integer, allocatable :: o(:,:,:) ! index of the orientation corresponding to theta, phi, psi
        complex(dp), allocatable :: in(:,:,:,:) ! for test purpose only
    end type o3_type
    type (o3_type) :: o3

    integer, allocatable :: indp(:,:,:) ! index of the projection corresponding to m, mup, mu
    integer, allocatable :: im(:) ! m for projection 1 to np
    integer, allocatable :: imup(:) ! mup for projection 1 to np. mup corresponds to phi
    integer, allocatable :: imu(:) ! mu for projection 1 to np. mu corresponds to psi

    !
    ! Direct correlation functions from Luc
    !
    type :: cq_type
        logical :: isok = .false.
        real(dp) :: q
        integer :: na ! number of projections for the direct correlation function
        integer :: nq=200 ! number of q points in c(q)
        integer, allocatable :: inda(:,:,:,:,:) ! inda(0:mmax,0:mmax,-mmax:mmax,-mmax:mmax,-mmax:mmax) ! m n mu nu khi
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

    public :: energy_cproj_slow

contains

    subroutine energy_cproj_slow (ff,df)
        ! use ieee_arithmetic
        use iso_c_binding, only: c_ptr, dp=>c_double
        use module_grid, only: grid
        use module_thermo, only: thermo
        implicit none
        real(dp), intent(out) :: ff
        real(dp), intent(out) :: df(:,:,:,:,:)
        real(dp) :: dv, kT
        logical :: q_eq_mq
        integer :: ix, iy, iz, ix_q, iy_q, iz_q, ix_mq, iy_mq, iz_mq
        real(dp), allocatable :: gamma(:,:,:,:,:,:), gamma_o(:,:,:,:), deltarho(:,:,:,:,:,:)
        complex(dp), allocatable :: deltarho_p(:,:,:,:), deltarho_p_q(:), deltarho_p_mq(:), gamma_p(:,:,:,:)
        integer :: nx, ny, nz, np, no, ns, ntheta, nphi, npsi, mmax, na, nq, p, molrotsymorder
        integer :: m, n, mu, nu, khi, mup, ia, ip, io, ipsi, iphi, itheta, iq, is, i,j,k
        complex(dp), allocatable :: gamma_p_q(:), gamma_p_mq(:)
        complex(dp) :: gamma_m_khi_mu_q, gamma_m_khi_mu_mq
        real(dp) :: tharm_sph(grid%ntheta,grid%np)
        real(dp) :: q(3), mq(3)
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
        real(dp) :: theta(grid%ntheta), wtheta(grid%ntheta)
        real(dp) :: lx, ly, lz, rho0
        logical, allocatable :: gamma_p_isok(:,:,:)
        complex(dp), allocatable :: foo(:,:,:)

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

        if (.not.allocated(fm)) allocate (fm(0:mmax) ,source= [( sqrt(real(2*m+1,dp)) ,m=0,mmax  )])

        allocate (deltarho(nx,ny,nz,ntheta,nphi,npsi) ,source=zero)
        allocate (gamma(nx,ny,nz,ntheta,nphi,npsi) ,source=zero)
        allocate (gamma_o(nx,ny,nz,no) ,source=zero)

        allocate (gamma_p_q(np), source=zeroc)
        allocate (gamma_p_mq(np), source=zeroc)
        allocate (gamma_p(np,nx,ny,nz) ,source=zeroc)
        allocate (deltarho_p(np,nx,ny,nz) ,source=zeroc)
        allocate (deltarho_p_q(np) ,source=zeroc)
        allocate (deltarho_p_mq(np) ,source=zeroc)
        allocate (gamma_p_isok(nx,ny,nz), source=.false.)

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
        ! indp(m,mup,mu) lie le triplet (m,mup,mu) a l'unique indice de projection ip
        ! im(p) donne inversement le m correspondant à l'indice p dans le tableau des projections
        ! On a donc indp(im(p),imup(p),imu(p)) == p
        !
        if (.not.allocated(indp)) allocate ( indp(0:mmax,-mmax:mmax, 0:mmax) ,source=-huge(1))
        if (.not.allocated(im))   allocate ( im(np) ,source=-huge(1))
        if (.not.allocated(imup)) allocate ( imup(np) ,source=-huge(1))
        if (.not.allocated(imu))  allocate ( imu(np) ,source=-huge(1))


        !
        ! Count projections without hermitian symetry
        !
        if (.not.allocated(p3%p)) then
            p3%np=0
            do m=0,mmax
                do mup=-m,m
                    do mu=-m,m
                        p3%np=p3%np+1
                    end do
                end do
            end do
            if (.not.allocated(p3%p)) allocate(p3%p(0:mmax,-mmax:mmax,-mmax:mmax), source=-huge(1))
            if (.not.allocated(p3%m)) allocate(p3%m(1:p3%np), source=-huge(1))
            if (.not.allocated(p3%mup)) allocate(p3%mup(1:p3%np), source=-huge(1))
            if (.not.allocated(p3%mu)) allocate(p3%mu(1:p3%np), source=-huge(1))
            p=0
            do m=0,mmax
                do mup=-m,m
                    do mu=-m,m
                        p=p+1
                        p3%p(m,mup,mu)=p
                        p3%m(p)=m
                        p3%mup(p)=mup
                        p3%mu(p)=mu
                    end do
                end do
            end do
            if (p3%np/=p) error stop "9866751823"
        end if


        ip=0
        do m=0,mmax
            do mup=-m,m
                !
                ! We chose to have indp and imu to contain the true value of mu, not mu/molrotsymorder.
                ! Thus, we loop over mu with steps of molrotsymorder
                ! Thus, we have holes in indp, but all calls to indp and imu should be done with this true mu
                ! For instance, if molrotsymorder=2 (for water or any C2V molecule), any call to indp(m,mup,mu) with mu=1 will return -999
                ! Also, the array imu only contains even values (des valeurs paires)
                !
                do mu=0,m!,molrotsymorder
                    ip=ip+1
                    IF (ip > np) ERROR STOP "p > np at line 166"
                    indp(m,mup,mu) = ip
                    im(ip) = m
                    imup(ip) = mup
                    imu(ip) = mu
                end do
            end do
        end do
        if (ip /= np) error stop "ip /= np in energy_cproj_slow"
        if ( any(abs(im)>mmax) .or. any(abs(imup)>mmax) .or. any(abs(imu)>mmax) ) then
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
            write(11,*) ip, im(ip), imup(ip), imu(ip)
        end do
        close(11)

        open(12, file="output/nonzero-projections_nohermitiansymetry.out")
        write(12,*)"Non-zero projections:"
        write(12,*)"        index         m          mup          mu"
        write(12,*)"        -----        ---         ---          --"
        do p=1,p3%np
            write(12,*) p, p3%m(p), p3%mup(p), p3%mu(p)
        end do
        close(12)

        !
        ! Initie poids et racines de l'integration angulaire par la methode de Gauss-Legendre
        ! Les racines sont les racines d'un polynomes en cos(theta).
        ! Comme nous voulons les theta correspondant, on en prend l'arccos(cos(theta))
        !
        theta=grid%thetaofntheta
        wtheta=grid%wthetaofntheta


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
        call test_routines_calcul_de_Rm_mup_mu_q

        tharm_sph = 0._dp
        do ip=1,np
            m = im(ip)
            mup = imup(ip)
            mu = imu(ip)
            do itheta=1,ntheta
                tharm_sph(itheta,ip) = harm_sph(m,mup,mu,theta(itheta))
            end do
        end do

        if (.not. allocated(p3%tharm_sph)) then
            allocate (p3%tharm_sph(p3%np,ntheta), source=0._dp)
            do itheta=1,ntheta
                do p=1,p3%np
                    m=p3%m(p)
                    mup=p3%mup(p)
                    mu=p3%mu(p)
                    p3%tharm_sph(p,itheta) = harm_sph(m,mup,mu,theta(itheta))
                end do
            end do
        end if


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! DEFINTION DE Δρ(r,ω) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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



        if (.not. allocated(o3%in)) allocate(o3%in(no,nx,ny,nz), source=zeroc)
        do concurrent (ix=1:nx, iy=1:ny, iz=1:nz, io=1:no)
            o3%in(io,ix,iy,iz) = complex(solvent(1)%density(ix,iy,iz,io)-rho0, 0._dp)
        end do



        if (.not. allocated(foo)) allocate (foo(ntheta,-mmax:mmax,-mmax:mmax) ,source=zeroc)
        if (.not. allocated(p3%drho)) allocate(p3%drho(p3%np,nx,ny,nz), source=zeroc)
        do iz=1,nz
            do iy=1,ny
                do ix=1,nx

                    foo=zeroc
                    do itheta=1,ntheta
                        do mup=-mmax,mmax
                            do mu=-mmax,mmax

                                do iphi=1,nphi
                                    do ipsi=1,npsi
                                        j=iphi-1
                                        k=ipsi-1
                                        io=grid%indo(itheta,iphi,ipsi)
                                        foo(itheta,mup,mu)=foo(itheta,mup,mu) &
                                             +complex(solvent(1)%density(ix,iy,iz,io)-rho0, 0._dp)&
                                             *exp(twopi*ii*(real(mup*j,dp)/real(nphi,dp)+real(mu*k,dp)/real(npsi,dp)))&
                                             /real(nphi*npsi,dp)
                                    end do
                                end do

                            end do
                        end do
                    end do

                    p3%drho(:,ix,iy,iz)=zeroc
                    do mup=-mmax,mmax
                        do mu=-mmax,mmax
                            do m=max(abs(mup),abs(mu)),mmax
                                p=p3%p(m,mup,mu)
                                do i=1,ntheta
                                    p3%drho(p,ix,iy,iz) = p3%drho(p,ix,iy,iz) &
                                    +fm(m)*wtheta(i)*foo(i,mup,mu)*p3%tharm_sph(p,i)
                                end do
                            end do
                        end do
                    end do

                end do
            end do
        end do

        do iz=1,nz
            do iy=1,ny
                do ix=1,nx
                    deltarho_p(1:np,ix,iy,iz) = euler2proj( deltarho(ix,iy,iz,1:ntheta,1:nphi,1:npsi) ) /real(nphi*npsi,dp)
                end do
            end do
        end do

        block
            complex(dp) :: a, b
            do ix=1,nx
                do iy=1,ny
                    do iz=1,nz
                        do m=0,mmax
                            do mup=-m,m
                                do mu=-m,m
                                    if(mu>=0) then
                                        a= p3%drho(p3%p(m,mup,mu),ix,iy,iz)
                                        b= deltarho_p(indp(m,mup,mu),ix,iy,iz)
                                        if (abs(a-b)>epsdp) then
                                            print*, a, b, a-b
                                            error stop "ifhdyr"
                                        end if
                                    else
                                        a= p3%drho(p3%p(m,mup,mu),ix,iy,iz)
                                        b= (-1)**(mup+mu)*conjg(deltarho_p(indp(m,-mup,-mu),ix,iy,iz))
                                        if (abs(a-b)>epsdp) then
                                            print*, a, b, a-b
                                            error stop "098766512723"
                                        end if
                                    end if
                                end do
                            end do
                        end do
                    end do
                end do
            end do
        end block


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

        do p=1,p3%np
            call dfftw_execute_dft( fft3d%plan_backward, p3%drho(p,:,:,:), p3%drho(p,:,:,:) )
        end do

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
        ! See function qproj
        ! And the inverse index that gives vector -q
        ! See function qproj_inv
        !


        call for_all_q_find_indices_of_mq




        ! rotation to molecular frame without symetry
        if (.not. allocated( p3%drho_q) ) allocate (p3%drho_q(p3%np) , source=zeroc)
        block
            complex(dp) :: sum_over_mup
            complex(dp) :: R(0:mmax,-mmax:mmax,-mmax:mmax)
            ! pour chaque q
            do ix_q=1,nx
                do iy_q=1,ny
                    do iz_q=1,nz

                        q = [ qx(ix_q), qy(iy_q), qz(iz_q) ]
                        p3%drho_q(:) = zeroc
                        R = rotation_matrix_between_complex_spherical_harmonics_lu (q, mmax)
                        do m=0,mmax
                            do khi=-m,m
                                do mu=-m,m

                                    sum_over_mup = zeroc
                                    do mup=-m,m
                                        sum_over_mup = sum_over_mup &
                                            + p3%drho(p3%p(m,mup,mu),ix_q,iy_q,iz_q) * R(m,mup,khi)
                                    end do
                                    p3%drho_q(p3%p(m,khi,mu)) = sum_over_mup

                                end do
                            end do
                        end do
                        p3%drho(:,ix_q,iy_q,iz_q) = p3%drho_q(:)

                    end do
                end do
            end do
        end block






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
!
!         !OZ sans sym
        if (.not. allocated(p3%gam) ) allocate ( p3%gam(p3%np,nx,ny,nz) )
p3%gam = p3%drho ! TODO FOR TESTING ONLY
!         p3%gam = zeroc
!         do ix_q=1,nx
!             do iy_q=1,ny
!                 do iz_q=1,nz
!                     iq = int( sqrt( qx(ix_q)**2 + qy(iy_q)**2 + qz(iz_q)**2 ) /cq%dq )+1
!                     if (iq>cq%nq) iq=cq%nq
!
!                     do khi=-mmax,mmax
!                         do m=abs(khi),mmax
!                             do mu=-m,m
!
!                                 do n=abs(khi),mmax
!                                     do nu=-n,n
!                                         if (mod(mu,2)/=0 .or. mod(nu,2)/=0) cycle
!                                         ia = cq%inda(m,n,mu,nu,khi)
! p3%gam( p3%p(m,khi,mu) ,ix_q,iy_q,iz_q) = p3%gam( p3%p(m,khi,mu) ,ix_q,iy_q,iz_q) &
!     + (-1)**(khi+nu) * ck(ia,iq) * p3%drho( p3%p(n,khi,-nu), ix_q,iy_q,iz_q)
!                                     end do
!                                 end do
!                             end do
!                         end do
!                     end do
!
!                 end do
!             end do
!         end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! rotation to fixed frame without symetry
        if (.not. allocated(p3%gam_q) ) allocate ( p3%gam_q(p3%np) )
        block
            complex(dp) :: sum_over_khi
            complex(dp) :: R(0:mmax,-mmax:mmax,-mmax:mmax)
            ! pour chaque q
            do ix_q=1,nx
                do iy_q=1,ny
                    do iz_q=1,nz

                        q = [ qx(ix_q), qy(iy_q), qz(iz_q) ]
                        p3%gam_q(:) = zeroc
                        R = rotation_matrix_between_complex_spherical_harmonics_lu (q, mmax)
                        do m=0,mmax
                            do khi=-m,m
                                do mu=-m,m

                                    sum_over_khi = zeroc
                                    do mup=-m,m
                                        sum_over_khi = sum_over_khi &
                                            + p3%gam(p3%p(m,khi,mu),ix_q,iy_q,iz_q) * conjg( R(m,mup,khi) )
                                    end do
                                    p3%gam_q(p3%p(m,khi,mu)) = sum_over_khi

                                end do
                            end do
                        end do
                        p3%gam(:,ix_q,iy_q,iz_q) = p3%gam_q(:)

                    end do
                end do
            end do
        end block



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! retour à l'espace reel
        do p=1,p3%np
            call dfftw_execute_dft( fft3d%plan_forward, p3%gam(p,:,:,:), p3%gam(p,:,:,:) )
        end do
        p3%gam=p3%gam/real(nx*ny*nz,dp)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! retour en orientations
        if (.not. allocated(o3%gam)) allocate (o3%gam(no,nx,ny,nz) ,source=zeroc)
        do iz=1,nz
            do iy=1,ny
                do ix=1,nx

                    foo(:,:,:)=zeroc
                    do mup=-mmax,mmax
                        do mu=-mmax,mmax
                            do itheta=1,ntheta

                                do m= max(abs(mup),abs(mu)) ,mmax
                                    p = p3%p(m,mup,mu)
                                    foo(itheta,mup,mu) = foo(itheta,mup,mu) +p3%gam(p,ix,iy,iz) * p3%tharm_sph(p,itheta)
                                end do

                            end do
                        end do
                    end do

                    do iphi=1,nphi
                        do ipsi=1,npsi
                            do itheta=1,ntheta

                                io=grid%indo(itheta,iphi,ipsi)
                                do mup=-mmax,mmax
                                    do mu=-mmax,mmax

                                        j=iphi-1
                                        i=ipsi-1
                                        o3%gam(io,ix,iy,iz) = o3%gam(io,ix,iy,iz)&
                                            + foo(itheta,mup,mu) *exp(-twopi*ii*real(mup*j,dp)/real(nphi,dp)&
                                            + real(mu*k,dp)/real(npsi,dp)) /real(nphi*npsi,dp)

                                    end do
                                end do

                            end do
                        end do
                    end do

                end do
            end do
        end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do ix=1,nx
    do iy=1,ny
        do iz=1,nz
            do io=1,grid%no
                print*, "ix,iy,iz,itheta,iphi,ipsi",ix,iy,iz,grid%theta(io), grid%phi(io), grid%psi(io),&
                    o3%in(io,ix,iy,iz), o3%gam(io,ix,iy,iz)
            end do
        end do
    end do
end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! df=0._dp
        ! ff=0._dp
        ! do is=1,ns
        !     do ix=1,nx
        !         do iy=1,ny
        !             do iz=1,nz
        !                 do itheta=1,ntheta
        !                     do iphi=1,nphi
        !                         do ipsi=1,npsi
        !                             io=grid%indo(itheta,iphi,ipsi)
        ! ff=ff-kT/2._dp*dv*gamma(ix,iy,iz,itheta,iphi,ipsi)&
        !         *(solvent(is)%density(ix,iy,iz,io)-solvent(is)%rho0)/0.0333_dp*grid%w(io)**2&
        !             *real(nphi*npsi*ntheta,dp)
        ! df(ix,iy,iz,io,is)=-kT*dv*gamma(ix,iy,iz,itheta,iphi,ipsi)&
        !                             /0.0333_dp*grid%w(io)**2*real(nphi*npsi*ntheta,dp)
        !                         end do
        !                     end do
        !                 end do
        !             end do
        !         end do
        !     end do
        ! end do















        stop "i would gladly abord"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

                    if (gamma_p_isok(ix_q,iy_q,iz_q).eqv..true.) then
                        if (gamma_p_isok(ix_mq, iy_mq,iz_mq).eqv..true.) then
                            cycle
                        else
                            print*,"in the main loop of energy_cproj_slow, i have already done q but not mq!"
                            print*, "ix_q, iy_q, iz_q =",ix_q, iy_q, iz_q
                            print*, "ix_mq, iy_mq, iz_mq =",ix_mq, iy_mq, iz_mq
                            error stop
                        end if
                    end if

                    q = [qx(ix_q), qy(iy_q), qz(iz_q)]
                    mq = [qx(ix_mq), qy(iy_mq), qz(iz_mq)]

                    if (ix_mq==ix_q .and. iy_mq==iy_q .and. iz_mq==iz_q) then ! this should only happen for ix=1 and ix=nx/2
                        q_eq_mq=.true.
                    else
                        q_eq_mq=.false.
                    end if

                    !
                    ! Move all projections to a smaller temporary array: deltarho_p_q (for q) and deltarho_p_mq (for -q)
                    !
                    deltarho_p_q (1:np) = deltarho_p(1:np,ix_q,iy_q,iz_q)
                    deltarho_p_mq(1:np) = deltarho_p(1:np,ix_mq,iy_mq,iz_mq)

                    if(any(deltarho_p_mq/=deltarho_p_mq)) then
                        do ip=1,np
                            print*, ip, ix_q, iy_q, iz_q, deltarho_p_mq(ip)
                        end do
                        error stop "lisfwlhukwlkuhwzlef"
                    end if
                    if(any(deltarho_p_q/=deltarho_p_q)) then
                        do iq=1,np
                            print*, ip, ix_q, iy_q, iz_q, deltarho_p_q(ip)
                        end do
                        error stop "hkhwkjhselkuhfkjhl"
                    end if

                    call rotate_to_molecular_frame ! does q and mq at the same time.

                    ! NE FONCTIONNE PAS
                    !block
                    !    complex(dp) :: a, b
                    !    do m=0,mmax
                    !        do khi=-m,m
                    !            do mu=0,m
                    !                a = deltarho_p_q( indp(m,khi,mu) )
                    !                b = p3%drho( p3%p(m,khi,mu) ,ix_q,iy_q,iz_q)
                    !                print*, ix_q,iy_q,iz_q,m,khi,mu, a, b, a-b
                    !            end do
                    !        end do
                    !    end do
                    !end block




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


                    if (q_eq_mq .and. any(deltarho_p_q/=deltarho_p_mq) ) then
                        print*, "q = mq mais apres la rotation vers le molecular frame, on n'a plus deltarho_p_q=deltarho_p_mq"
                        print*, "ix_q,iy_q,iz_q=", ix_q,iy_q,iz_q
                        print*, "ix_mq,iy_mq,iz_mq=", ix_mq, iy_mq, iz_mq
                        print*, "deltarho_p_q", deltarho_p_q
                        print*, "deltarho_p_mq", deltarho_p_mq
                        error stop "oijazd"
                    end if

                    !
                    ! c^{m,n}_{mu,nu,chi}(|q|) is tabulated for 200 values of |q|.
                    ! Find the tabulated value that is closest to |q|. Its index is iq.
                    ! Note |q| = |-q| so iq is the same for both vectors.
                    !
                    iq = int(norm2(q)/cq%dq)+1
                    if (iq>nq) then
                        iq=nq
                    else if (iq<=0) then
                        print*, "in energy_cproj_slow, ik <=0"
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

                                        ia = cq%inda(m,n,mu,nu,khi) ! the index of the projection of c(q). 1<=ia<na
                                        if (ia<=0) then
                                            print*,"ia=",ia
                                            print*,"for m, n, mu, nu, khi=",m,n,mu,nu,khi
                                            error stop "ia<=0"
                                        else if (ia>cq%na) then
                                            print*,"ia=",ia
                                            print*,"for m, n, mu, nu, khi=",m,n,mu,nu,khi
                                            error stop "ia>na"
                                        end if

                                        ip = indp(n,khi,abs(nu))

                                        if (nu<0) then ! no problem with delta rho (n, khi, -nu) since -nu>0. Thus, we apply eq. 1.30 directly
                                            gamma_m_khi_mu_q  = gamma_m_khi_mu_q  + (-1)**(khi+nu) *ck(ia,iq) *deltarho_p_q(ip)
                                            gamma_m_khi_mu_mq = gamma_m_khi_mu_mq + (-1)**(khi+nu) *ck(ia,iq) *deltarho_p_mq(ip)
                                        else ! transform delta rho (n, khi, -nu)(q) into conjg( deltarho(n,khi,nu)(-q) )
                                            gamma_m_khi_mu_q  = gamma_m_khi_mu_q  + (-1)**(n) *ck(ia,iq) *conjg(deltarho_p_mq(ip))
                                            gamma_m_khi_mu_mq = gamma_m_khi_mu_mq + (-1)**(n) *ck(ia,iq) *conjg(deltarho_p_q(ip))
                                        end if

                                    end do
                                end do

                                if (gamma_m_khi_mu_q/=gamma_m_khi_mu_q) error stop "jhwdlijweflkhwke"
                                if (gamma_m_khi_mu_mq/=gamma_m_khi_mu_mq) error stop "lijsfdkuglserhk"

                                ip=indp(m,khi,mu)
                                gamma_p_q(ip)  = gamma_m_khi_mu_q
                                gamma_p_mq(ip) = gamma_m_khi_mu_mq

                                if (q_eq_mq .and. any(gamma_p_q/=gamma_p_mq) ) then
                                    print*, "we have q=mq"
                                    print*, "for ix_q,iy_q,iz_q=",  ix_q, iy_q, iz_q
                                    print*, "and ix_mq,iy_mq,iz_mq=", ix_mq, iy_mq, iz_mq
                                    print*, "q=",q
                                    print*, "mq=",mq
                                    print*, "gamma_p_q=",gamma_p_q
                                    print*, "gamma_p_mq=",gamma_p_mq
                                    print*, "deltarho_p_q=",deltarho_p_q
                                    print*, "deltarho_p_mq=",deltarho_p_mq
                                    error stop "apzokdukhsf"
                                end if

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
                            print*, "old gamma_p =",ip, im(ip), imup(ip), imu(ip), gamma_p(ip,ix_q,iy_q,iz_q)
                            print*, "new gamma_p =",ip, im(ip), imup(ip), imu(ip), gamma_p_q(ip)
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
                            print*, "old gamma_mq =",ip, im(ip), imup(ip), imu(ip), gamma_p(ip,ix_mq,iy_mq,iz_mq)
                            print*, "new gamma_mq =",ip, im(ip), imup(ip), imu(ip), gamma_p_mq(ip)
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
        gamma=0._dp
        call dfftw_plan_dft_c2r_2d( ifft2d%plan, npsi, nphi, ifft2d%in, ifft2d%out, FFTW_EXHAUSTIVE)
        do iz=1,nz
            do iy=1,ny
                do ix=1,nx
                    gamma(ix,iy,iz,1:ntheta,1:nphi,1:npsi) = proj2euler( gamma_p(1:np,ix,iy,iz) )
                end do
            end do
        end do

        df=0._dp
        ff=0._dp
        do is=1,ns
            do ix=1,nx
                do iy=1,ny
                    do iz=1,nz
                        do itheta=1,ntheta
                            do iphi=1,nphi
                                do ipsi=1,npsi
                                    io=grid%indo(itheta,iphi,ipsi)
        ff=ff-kT/2._dp*dv*gamma(ix,iy,iz,itheta,iphi,ipsi)&
                *(solvent(is)%density(ix,iy,iz,io)-solvent(is)%rho0)/0.0333_dp*grid%w(io)**2&
                    *real(nphi*npsi*ntheta,dp)
        df(ix,iy,iz,io,is)=-kT*dv*gamma(ix,iy,iz,itheta,iphi,ipsi)&
                                    /0.0333_dp*grid%w(io)**2*real(nphi*npsi*ntheta,dp)
                                end do
                            end do
                        end do
                    end do
                end do
            end do
        end do


contains


        subroutine rotate_to_fixed_frame
            implicit none
            complex(dp) :: gamma_p_q_temp(1:np)
            complex(dp) :: gamma_p_mq_temp(1:np)
            complex(dp) :: gshrot (0:grid%mmax, -grid%mmax:grid%mmax, -grid%mmax:grid%mmax)
            complex(dp) :: gshrot_mq (0:grid%mmax, -grid%mmax:grid%mmax, -grid%mmax:grid%mmax)
            integer :: m, khi, mu, mup, ip, ip2
            if (ix_q==1 .and. iy_q==1 .and. iz_q==1) then ! that is if |q|=0. Rotation is arbitrary: we chose the Identity and thus have nothing to do.
                !return
            else
                gamma_p_q_temp = zeroc
                gamma_p_mq_temp = zeroc
                gshrot = conjg(rotation_matrix_between_complex_spherical_harmonics_lu(q,mmax) )
                gshrot_mq = conjg(rotation_matrix_between_complex_spherical_harmonics_lu(mq,mmax) )
                do m=0,mmax
                    do mup=-m,m
                        do mu=0,m
                            ip=indp(m,mup,mu)
                            ! Equation 1.22
                            do khi=-m,m
                                ip2=indp(m,khi,mu)
                                gamma_p_q_temp(ip) = gamma_p_q_temp(ip)   + gamma_p_q(ip2)*gshrot(m,mup,khi)
                                ! gamma_p_mq_temp(ip) = gamma_p_mq_temp(ip) + gamma_p_mq(ip2)*(-1)**m*gshrot(m,mup,-khi)
                                gamma_p_mq_temp(ip) = gamma_p_mq_temp(ip) + gamma_p_mq(ip2)*gshrot_mq(m,mup,khi)

                                if (q_eq_mq .and. gamma_p_q_temp(ip)/=gamma_p_mq_temp(ip)) then
                                    print*, "q == mq mais pourtant deltarho_p_q_temp /= deltarho_p_mq_temp"
                                    print*, "m,khi,mu,ip,mup,ip2,gamma_p_q_temp(ip),gamma_p_mq_temp(ip)"
                                    print*,m,khi,mu,ip,mup,ip2,gamma_p_q_temp(ip),gamma_p_mq_temp(ip)
                                    error stop
                                end if

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
            complex(dp) :: gshrot_mq (0:grid%mmax, -grid%mmax:grid%mmax, -grid%mmax:grid%mmax)
            integer :: m, khi, mu, mup, ip2, ip
            !
            ! Rotate projections from laboratory frame to molecular frame
            ! If q=0 then chose the rotation you want: the identity is the simplest option.
            !
            if (ix_q==1 .and. iy_q==1 .and. iz_q==1) then ! that is if |q|=0. Rotation is arbitrary: we chose the Identity and thus have nothing to do.
                ! return
            else
                gshrot    = rotation_matrix_between_complex_spherical_harmonics_lu ( q, mmax)
                gshrot_mq = rotation_matrix_between_complex_spherical_harmonics_lu ( mq, mmax)

                deltarho_p_q_temp = zeroc
                deltarho_p_mq_temp = zeroc
                do m=0,mmax
                    do khi=-m,m
                        do mu=0,m
                            ip=indp(m,khi,mu)
                            ! Equation 1.22
                            do mup=-m,m
                                ip2=indp(m,mup,mu)
                                deltarho_p_q_temp(ip) = deltarho_p_q_temp(ip)  + deltarho_p_q(ip2)*gshrot(m,mup,khi)
                                ! Eq. 1.23 We don't need to compute gshrot for -q. We do q and -q at the same time.
                                ! deltarho_p_mq_temp(ip) = deltarho_p_mq_temp(ip) + deltarho_p_mq(ip2)*(-1)**m*gshrot(m,mup,-khi)
                                deltarho_p_mq_temp(ip) = deltarho_p_mq_temp(ip) + deltarho_p_mq(ip2)*gshrot_mq(m,mup,khi)
                                if (q_eq_mq .and. deltarho_p_q_temp(ip)/=deltarho_p_mq_temp(ip)) then
                                    print*, "q == mq mais pourtant deltarho_p_q_temp /= deltarho_p_mq_temp"
                                    print*,m,khi,mu,ip,mup,ip2,deltarho_p_q_temp(ip),deltarho_p_mq_temp(ip)
                                    error stop
                                end if
                                ! print*, "ix, iy, iz =", ix_q, iy_q, iz_q,"m,mup,mu=", m, mup, khi
                                ! print*, "q=", q
                                ! print*, "mq=", mq
                                ! print*, "gshrot(q,m,mup,khi)                 =", gshrot(m,mup,khi)
                                ! print*, "gshrot(mq,m,mup,khi) par sym 1.23/1 =", (-1)**m*gshrot(m,mup,-khi)
                                ! print*, "gshrot(mq,m,mup,khi) par sym 1.23/2 =", (-1)**(m+mup+khi)*conjg(gshrot(m,-mup,khi))
                                ! print*, "gshrot(mq,m,mup,khi) brutal         =", gshrot_mq(m,mup,khi)
                                ! print*,
                            end do
                        end do
                    end do
                end do
                deltarho_p_q = deltarho_p_q_temp
                deltarho_p_mq = deltarho_p_mq_temp
            end if
        end subroutine rotate_to_molecular_frame


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
            REAL(dp), PARAMETER :: twopi=2*acos(-1._dp)
            INTEGER, INTENT(IN) :: i, imax ! index the node and number of nodes in the direction you wish
            ! If you use REAL to complex transforms, and thus have in the first direction nx/2+1 nodes, you should still pass nx to this
            ! routine, not nx/2+1
            REAL(dp), INTENT(IN) :: length ! size of the supercell
            IF( i == 1 ) THEN
                qproj = 0._dp
            ELSE IF( i <= imax/2  ) THEN
                qproj = twopi*REAL(i-1,dp)/length
            ELSE
                qproj = twopi*REAL(i-1-imax,dp)/length
            END IF
        END FUNCTION qproj

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


            i=(1+mmax)*(1+2*mmax)*(3+2*mmax)*(5+4*mmax*(2+mmax))/15
          !  if (na/=i) then
          !      print*, "in read_ck_cproj nalpha est bizarre"
          !      print*, "na=",na
          !      print*, "it should be (1+mmax)*(1+2*mmax)*(3+2*mmax)*(5+4*mmax*(2+mmax))/15 =",i
          !      ! error stop
          !  end if

            i=sum([([([([([( 1 ,nu=-n,n)] ,mu=-m,m)], n=abs(khi),mmax)] ,m=abs(khi),mmax)] ,khi=-mmax,mmax)]  )
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

            allocate (cq%inda(0:mmax,0:mmax,-mmax:mmax,-mmax:mmax,-mmax:mmax), source=-huge(1)) ! m n mu nu khi. -huge is used to spot more easily bugs that may come after
            do ia=1,na
                m = m_ck(ia)
                n = n_ck(ia)
                mu = mu_ck(ia)
                nu = nu_ck(ia)
                khi = khi_ck(ia)
                cq%inda(m,n,mu,nu,khi) = ia
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

        FUNCTION euler2proj (eulerangles) RESULT (projections)
            use module_rotation, only: harm_sph
            IMPLICIT NONE
            REAL(dp), INTENT(IN) :: eulerangles(1:grid%ntheta,1:grid%nphi,1:grid%npsi)
            COMPLEX(dp) :: proj_theta(1:grid%ntheta,0:grid%mmax/grid%molrotsymorder,-grid%mmax:grid%mmax) ! itheta,mu,mup    note we changed mu(psi) to the 2nd position
            COMPLEX(dp) :: projections(1:grid%np)
            INTEGER :: itheta, iphi, ipsi, m, mup, mu, ip, mmax, molrotsymorder
            complex(dp), allocatable :: test_explicit(:,:,:), proj_theta_full(:,:,:)
            complex(dp), allocatable :: proj_m_mup_mu(:,:,:)
            mmax=grid%mmax
            molrotsymorder=grid%molrotsymorder
            projections = zeroc
            proj_theta = zeroc

            allocate (proj_theta_full(ntheta,-mmax:mmax,-mmax:mmax), source=zeroc)
            allocate (test_explicit(ntheta,-mmax:mmax,-mmax:mmax), source=zeroc)
            do itheta=1,ntheta
                do mu=-mmax,mmax
                    do mup=-mmax,mmax
                        do ipsi=1,grid%npsi
                            do iphi=1,grid%nphi
                                test_explicit(itheta,mu,mup) = test_explicit(itheta,mu,mup)&
                                +eulerangles(itheta,iphi,ipsi) &
                                *exp(ii*mup*grid%phiofnphi(iphi)) *exp(ii*mu*grid%psiofnpsi(ipsi))
                            end do
                        end do
                    end do
                end do
            end do


            do itheta=1,ntheta
                do iphi=1,nphi
                    do ipsi=1,npsi
                        fft2d%in(ipsi,iphi) = eulerangles(itheta,iphi,ipsi)
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
                proj_theta(itheta,0:mmax/molrotsymorder,0:mmax)   = CONJG( fft2d%out(:,1:mmax+1) )
                proj_theta(itheta,0:mmax/molrotsymorder,-mmax:-1) = CONJG( fft2d%out(:,mmax+2:) )

                proj_theta_full(itheta,0:mmax,0:mmax)     = fft2d_c%out(1:npsi/2+1, 1:nphi/2+1)
                proj_theta_full(itheta,0:mmax,-mmax:-1)   = fft2d_c%out(1:npsi/2+1, nphi/2+2:)
                proj_theta_full(itheta,-mmax:-1,0:mmax)   = fft2d_c%out(npsi/2+2:, 1:nphi/2+1)
                proj_theta_full(itheta,-mmax:-1,-mmax:-1) = fft2d_c%out(npsi/2+2:, nphi/2+2:)

            end do
            !
            !
            ! do itheta=1,ntheta
            !     do mup=-mmax,mmax
            !         do mu=-mmax,mmax
            !             print*,"itheta, mup, mu =",itheta,mup,mu
            !             print*,"explicit                      =",test_explicit(itheta,mu,mup)
            !             print*,"proj_theta_full(itheta,mu,mup)=",proj_theta_full(itheta,mu,mup)
            !             print*,"proj_theta_full reconstruit   =",(-1)**(-mup-mu)*conjg(proj_theta_full(itheta,-mu,-mup))," WWWW"
            !             if(mu>=0) print*,"proj_theta(itheta,mu,mup)     =",proj_theta(itheta,mu,mup)
            !             if(mu<0 ) print*,"proj_theta(itheta,mu,mup)     =",(-1)**(-mup-mu)*conjg(proj_theta(itheta,-mu,-mup))," reconstruit"
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
                            proj_theta_full(itheta,mu,mup)*wtheta(itheta)*harm_sph(m,mup,mu,theta(itheta))*fm(m)
                        end do
                    end do
                end do
            end do
            do m=0,mmax
                do mup=-m,m
                    do mu=-m,m
                        ! Check relation 1.7
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
                        ip = indp( m,mup,mu )
                        projections(ip)= sum(proj_theta(:,mu,mup)*tharm_sph(:,ip)*wtheta(:))*fm(m)
                        if (abs(projections(ip)-proj_m_mup_mu(m,mup,mu))>epsilon(1._dp)) then
                            print*, "projections(ip)/=proj_m_mup_mu(m,mup,mu)"
                            print*, "projections(ip)         =",projections(ip)
                            print*, "proj_m_mup_mu(m,mup,mu) =",proj_m_mup_mu(m,mup,mu)
                            error stop
                        end if
                    end do
                end do
            end do
        END FUNCTION euler2proj

        FUNCTION proj2euler (projections) RESULT (eulerangles)
            IMPLICIT NONE
            REAL(dp) :: eulerangles(1:grid%ntheta,1:grid%nphi,1:grid%npsi)
            COMPLEX(dp), INTENT(IN) :: projections(1:grid%np)
            COMPLEX(dp) :: proj_theta(1:grid%ntheta,0:grid%mmax/grid%molrotsymorder,-grid%mmax:grid%mmax)
            INTEGER :: m, mup, mu, ip, itheta, iphi, ipsi, mmax, molrotsymorder
            mmax=grid%mmax
            molrotsymorder=grid%molrotsymorder
            proj_theta = CMPLX(0,0)
            do itheta=1,ntheta
                do mup=-mmax,mmax
                    do mu=0,mmax/molrotsymorder
                        do m= MAX(ABS(mup),ABS(mu)), mmax
                            ip=indp(m,mup,mu)
                            proj_theta(itheta,mu,mup) = proj_theta(itheta,mu,mup)&
                            +projections(ip)*tharm_sph(itheta,ip)*fm(m)
                        end do
                    end do
                end do
            end do
            do itheta=1,ntheta
                ifft2d%in(:,1:mmax+1) = CONJG( proj_theta(itheta,:,0:mmax)   )
                ifft2d%in(:,mmax+2:)  = CONJG( proj_theta(itheta,:,-mmax:-1) )
                call dfftw_execute( ifft2d%plan )
                do iphi=1,nphi
                    do ipsi=1,npsi
                        eulerangles(itheta,iphi,ipsi) = ifft2d%out(ipsi,iphi) *(nphi*npsi)
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
            qx(1:nx) = [( qproj(ix,nx,lx), ix=1,nx )]
            qy(1:ny) = [( qproj(iy,ny,ly), iy=1,ny )]
            qz(1:nz) = [( qproj(iz,nz,lz), iz=1,nz )]

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
                print*, "that is not compatible with our energy_cproj_slow for now"
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
                            qx=qproj(ix_q,nx,lx)
                            qy=qproj(iy_q,ny,ly)
                            qz=qproj(iz_q,nz,lz)
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


    end subroutine energy_cproj_slow

end module module_energy_cproj_slow

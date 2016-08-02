module module_orientation_projection_transform

    use iso_c_binding, only: C_PTR, C_INT, C_INT32_T, C_INTPTR_T, C_DOUBLE_COMPLEX, C_DOUBLE, C_FUNPTR, C_SIZE_T, C_FLOAT, &
    C_FLOAT_COMPLEX, C_CHAR
    use precision_kinds, only: dp
    use module_grid, only: grid
    use module_solvent, only: solvent
    use module_wigner_d, only: wigner_small_d

    implicit none
    private

    !
    ! FFTW3 header - modern (fortran 2003) version. Expects iso_c_binding
    !
    include 'fftw3.f03'

    complex(dp), allocatable, protected :: f_theta_mu2_mup(:,:,:)
    type :: p3_type
        real(dp), allocatable :: wigner_small_d(:,:) ! tabulation des harmoniques sphériques r(m,mup,mu,theta) en un tableau r(itheta,p)
        integer, allocatable :: p(:,:,:) ! index of the projection corresponding to m, mup, mu
        integer, allocatable :: m(:) ! m for projection 1 to np
        integer, allocatable :: mup(:) ! mup for projection 1 to np. mup corresponds to phi
        integer, allocatable :: mu(:) ! mu for projection 1 to np. mu corresponds to psi
    end type p3_type
    type (p3_type) :: p3
    type :: fft_type
        type(c_ptr) :: plan2dp, plan2dm ! plans for 2D FFTs with sign +(p) or -(m) in the exponential
    end type fft_type
    type(fft_type), protected :: fft
    real(dp), allocatable :: r2d(:,:)
    complex(dp), allocatable :: c2d(:,:)
    integer :: nx, ny, nz, np, no, ntheta, nphi, npsi, mmax, mrso
    real(dp), allocatable :: wtheta(:)
    real(dp), parameter :: fm(0:5) = [ 1._dp, sqrt(3._dp), sqrt(5._dp), sqrt(7._dp), sqrt(9._dp), sqrt(11._dp) ] ! sqrt(2m+1)
    logical :: is_init = .false.

    public :: angl2proj, proj2angl

contains

    subroutine init
        use omp_lib
        use iso_c_binding, only: c_ptr
        use precision_kinds, only: dp
        use module_grid, only: grid
        implicit none

        if(is_init) return

        ntheta=grid%ntheta
        nphi=grid%nphi
        npsi=grid%npsi
        mmax=grid%mmax
        mrso=grid%molrotsymorder
        nx=grid%nx
        ny=grid%ny
        nz=grid%nz
        np=grid%np
        no=grid%no

        allocate (wtheta(grid%ntheta))
        allocate( r2d(npsi,nphi) )
        allocate( c2d(npsi/2+1,nphi))
        allocate( f_theta_mu2_mup(ntheta,0:mmax/mrso,-mmax:mmax)   ,  source=(0._dp,0._dp)   )


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
        block
            integer :: ip, m, mup, mu
            allocate ( p3%p(0:mmax,-mmax:mmax, 0:mmax/mrso) ,source=-huge(1)) ! Dans p3%p, on met mu2=2*mu, pas mu
            allocate ( p3%m(np) ,source=-huge(1)) ! mettre des -huge comme valeur initiale permet de plus tard verifier que s'il reste un -huge qqpart, il y a un problème. Si on avait mis 0, on ne peut pas vérifier.
            allocate ( p3%mup(np) ,source=-huge(1))
            allocate ( p3%mu(np) ,source=-huge(1)) ! c'est bien mu qu'on met dans p3%mu(), pas mu2
            ip=0
            do m=0,mmax
                do mup=-m,m
                    do mu=0,m,mrso
                        ip=ip+1
                        if (ip > np) ERROR STOP "p > np at line 166"
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
        end block

        !
        ! L'initialisation des racines et des poids de la quadrature angulaire a été faite beaucoup plus tot
        ! dans le module_grid.
        ! On fait attention qu'on a deux tableaux des racines et des poids:
        ! - le premier tableau, grid%theteaofntheta, est un tableau de taille  le nombre de theta (ntheta)
        ! - le deuxième, grid%theta, est de taille le nombre total d'angle. Il retourne le theta qui correspond à l'indice de chaque angle.
        ! Si on veut la liste de tous les theta, on utilisera donc le tableau grid%thetaofntheta.
        !
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
        block
            integer :: p, m, mup, mu, i
            allocate ( p3%wigner_small_d(1:ntheta, 1:np) ,source=0._dp)
            do p=1,np
                m = p3%m(p)
                mup = p3%mup(p)
                mu = p3%mu(p)
                do i=1,ntheta
                    ! Pour chaque theta, calcule la fonction de Wigner-d correspondant à toutes les projections avec la méthode de Wigner.
                    p3%wigner_small_d(i,p) =  wigner_small_d( m,mup,mu,  grid%thetaofntheta(i)  )
                end do
            end do
        end block

        !
        ! Les projections sont complexes, mais la symetrie hermitienne
        ! permet de ne garder que les mup>=0 ou les mu>=0.
        ! On choisit les mu>=0 (cf doc de Luc)
        !

        ! ON PREPARE LES FFT
        ! Pay attention to the definition of "forward" and "backward" for the FFTW library. Source: http://www.fftw.org/doc/Real_002ddata-DFTs.html
        ! FFTW_FORWARD = -1
        ! FFTW_BACKWARD = +1
        ! R2C == FORWARD = -1
        ! C2R == BACKWARD = +1
        select case(dp)
        case(c_double)
            call dfftw_plan_dft_r2c_2d(  fft%plan2dm, npsi, nphi, r2d, c2d, FFTW_EXHAUSTIVE ) ! npsi est en premier indice
            call dfftw_plan_dft_c2r_2d(  fft%plan2dp, npsi, nphi, c2d, r2d, FFTW_EXHAUSTIVE )
        case(c_float)
            call sfftw_plan_dft_r2c_2d(  fft%plan2dm, npsi, nphi, r2d, c2d, FFTW_EXHAUSTIVE ) ! npsi est en premier indice
            call sfftw_plan_dft_c2r_2d(  fft%plan2dp, npsi, nphi, c2d, r2d, FFTW_EXHAUSTIVE )
        end select

        is_init = .true.
    end subroutine init


    subroutine angl2proj(o,p)
        implicit none
        real(dp), intent(in) :: o(:) ! orientations from 1 to no
        complex(dp), intent(out) :: p(:) ! projections from 1 to no
        integer :: itheta, iphi, ipsi, m, mup, mu2, ip, io
        if( .not. is_init) call init
        p = (0._dp,0._dp)
        f_theta_mu2_mup = (0._dp,0._dp)
        io=0
        do itheta=1,ntheta
            do iphi=1,nphi
                do ipsi=1,npsi
                    io=io+1
                    r2d(ipsi,iphi) = o(io)
                end do
            end do
            select case(dp)
            case(c_double)
                call dfftw_execute (fft%plan2dm)
            case(c_float)
                call sfftw_execute (fft%plan2dm)
            end select
            c2d = conjg(c2d)/real(nphi*npsi,dp) ! we wanted plan with + in exponential but did not have choice since r2c is done with sign - in fftw.
            f_theta_mu2_mup(itheta,0:mmax/mrso,0:mmax)   = c2d(:,1:mmax+1)
            f_theta_mu2_mup(itheta,0:mmax/mrso,-mmax:-1) = c2d(:,mmax+2:)
        end do
        ip=0
        do m=0,mmax
            do mup=-m,m
                do mu2=0,m/mrso
                    ip=ip+1
                    p(ip)= sum(f_theta_mu2_mup(:,mu2,mup)*p3%wigner_small_d(:,ip)*wtheta(:))*fm(m)
                end do
            end do
        end do
    end subroutine angl2proj


    subroutine proj2angl(p,o)
        implicit none
        complex(dp), intent(in) :: p(:) ! np
        real(dp), intent(out) :: o(:) ! no
        integer :: itheta, iphi, ipsi, m, mup, mu2, ip, io
        if(.not.is_init) call init
        o = 0._dp
        f_theta_mu2_mup = (0._dp,0._dp)
        do mup=-mmax,mmax
            do mu2=0,mmax/mrso
                do m= max(abs(mup),mrso*abs(mu2)), mmax
                    ip=p3%p(m,mup,mu2)
                    do itheta=1,ntheta
                        f_theta_mu2_mup(itheta,mu2,mup) = f_theta_mu2_mup(itheta,mu2,mup) +p(ip)*p3%wigner_small_d(itheta,ip)*fm(m)
                    end do
                end do
            end do
        end do
        io=0
        do itheta=1,ntheta
            c2d(:,1:mmax+1) = conjg( f_theta_mu2_mup(itheta,0:mmax/mrso,0:mmax)   )
            c2d(:,mmax+2:)  = conjg( f_theta_mu2_mup(itheta,0:mmax/mrso,-mmax:-1) )
            select case(dp)
            case(c_double)
                call dfftw_execute( fft%plan2dp )
            case(c_float)
                call sfftw_execute( fft%plan2dp )
            end select
            do iphi=1,nphi
                do ipsi=1,npsi
                    io=io+1
                    o(io) = r2d(ipsi,iphi)
                end do
            end do
        end do
    end subroutine proj2angl


end module module_orientation_projection_transform

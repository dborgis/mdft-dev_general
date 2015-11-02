module module_energy_cproj
    use precision_kinds, only: dp
    implicit none
    private
    integer :: nx, ny, nz, no, mmax, p, m, mup, mu
    integer :: np ! number of projections => should be moved to grid%np
    integer, allocatable :: indp(:,:,:) ! index of the projection corresponding to triplet m,mup,mu
    integer, allocatable :: indm(:) ! m(alpha)
    integer, allocatable :: indmup(:) ! mup(alpha)
    integer, allocatable :: indmu(:) ! mu(alpha)
    real(dp), allocatable :: tharm_sph(:,:)
    real(dp), allocatable :: fm(:)
    public :: energy_cproj
contains
    subroutine energy_cproj (ff, df)
        use module_debug, only: debugmode
        use precision_kinds, only: dp
        use mathematica, only: fact ! factorial(n) (type=double precision)
        use module_solvent, only: solvent, print_solvent_not_allocated
        use module_grid, only: grid
        use module_fft, only: fft2d
        implicit none
        real(dp), intent(out) :: ff
        real(dp), intent(in) :: df(:,:,:,:,:)
        real(dp) :: ti, tf
        integer :: ix, iy, iz, p
        complex(dp), allocatable :: rhop(:,:,:,:) ! ijk,io
        integer :: ip, itheta, iphi, ipsi
        real(dp), allocatable :: test_density(:)
        complex(dp), allocatable :: test_rhop(:)
print*,"AAAAAAAAAA"
        allocate (fm(0:grid%mmax), source= [( sqrt(real(2*m+1,dp))/real(grid%nphi*grid%npsi,dp) ,m=0,grid%mmax)])
print*,"BBBBBBBBB"
        ff = 0._dp

        if (.not.allocated(solvent)) call print_solvent_not_allocated("In module_energy_cproj")
        call cpu_time (ti)
        !
        ! Toutes les projections sont stockées dans un meme vecteur de taille nproj
        ! iproj(m,mup,mu) lie le triplet (m,mup,mu) a l'unique indice de projection
        ! im(p) donne inversement le m correspondant à l'indice p dans le tableau des projections
        ! On a donc iproj(im(p),imup(p),imu(p)) == p
        !
        mmax=grid%mmax
        nx=grid%nx
        ny=grid%ny
        nz=grid%nz
        no=grid%no
        np=grid%np

        print*, "in energy_cproj, we have "
        print*, no,"orientations and",np,"projections"

        !
        ! We chose to have iproj and imu to contain the true value of mu, not mu/molrotsymorder.
        ! Thus, we loop over mu with steps of molrotsymorder
        ! Thus, we have holes in iproj, but all calls to iproj and imu should be done with this true mu
        ! For instance, if molrotsymorder=2 (for water or any C2V molecule), any call to iproj(m,mup,mu) with mu=1 will return -999
        ! Also, the array imu only contains even values (des valeurs paires)
        !
        allocate (indp(0:mmax,-mmax:mmax,0:mmax), source=huge(1)) ! huge is here to make it clear something's wrong while code execution p(m,mup,mu)
        allocate (indm(np) ,source=huge(1))
        allocate (indmup(np) ,source=huge(1))
        allocate (indmu(np) ,source=huge(1))
print*,"CCCCCCCCCCCCCCCC"
        p=0
        do m=0,mmax
            do mup=-m,m
                do mu=0,m,grid%molrotsymorder
                    p=p+1
                    if (p>np) error stop "p > nombre de projections in energy_cproj"
                    indp(m,mup,mu)=p
                    indm(p)=m
                    indmup(p)=mup
                    indmu(p)=mu
                end do
            end do
        end do
        if (p/=np) error stop "in energy_cproj, p/=np after loop over all projections"

        if (any(abs(indm)>mmax) .or. any(abs(indmup)>mmax) .or. any(abs(indmu)>mmax)) then
            print*, "in energy_cproj, we have some incorrect value for m, mup or mu"
            error stop
        end if

print*,"DDDDDDDDDDDDDDDDDD"

        !
        ! Print all projections that will be kept in memory
        !
        if (debugmode) then
            open(11, file="output/nonzero-projections.out")
            write(11,*)"Non-zero projections:"
            write(11,*)"        index         m          mup          mu"
            write(11,*)"        -----        ---         ---          --"
            do p=1,np
                write(11,*) p, indm(p), indmup(p), indmu(p)
            end do
            close(11)
        end if
print*,"EEEEEEEEEEEEEEEEEEEEEEE"

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
        allocate( tharm_sph(grid%ntheta,grid%np), source=0._dp)
        tharm_sph = 0._dp
        do p=1,np
            m = indm(p)
            mup = indmup(p)
            mu = indmu(p)
            do itheta=1,grid%ntheta
                tharm_sph(itheta,p) = harm_sph(m,mup,mu,grid%thetaofitheta(itheta))
            end do
        end do
print*, "FFFFFFFFFFFFFFFFFFFFFFFF"

        !
        ! On passe tout de suite en projection, dans le repère cartesien
        ! ces projections sont complexes, mais la symetrie hermitienne
        ! permet de ne garder que les mup>=0 ou les mu>=0.
        ! On choisit les mu>=0 (cf doc de Luc)
        !
        ! call dfftw_plan_dft_r2c_2d (plan2dr2c, grid%npsi, grid%nphi, fft2d)
        if (.not.fft2d%isok) call fft2d%init
print*, "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG"
        !
        ! Project the density
        !
        if (.not.allocated(rhop)) allocate (rhop(np,nx,ny,nz))
        rhop=cmplx(0,0)
        do iz=1,nz
            do iy=1,ny
                do ix=1,nx
                    print*, "let's call o2p"
                    call o2p (solvent(1)%density(ix,iy,iz,1:no), rhop(1:np,ix,iy,iz))
                    print*, ix, iy,iz, "ok"
                end do
            end do
        end do
        ! catch NaN or Inf
        if (any(rhop/=rhop)) then
            print*, "In module_energy_cproj, I caugth a NaN or Inf in rhop"
            error stop
        end if
print*, "HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH"


        !
        ! Gather projections
        !
        if (.not.allocated(test_density)) allocate(test_density(1:no))
        if (.not.allocated(test_rhop)) allocate(test_rhop(1:np))
        do iz=1,nz
            do iy=1,ny
                do ix=1,nx
                    call p2o (rhop(1:np,ix,iy,iz), test_density(1:no))
                    call o2p (test_density, test_rhop)
                    print*, "testing o2p+p2o"
                    print*, norm2(solvent(1)%density(ix,iy,iz,1:no)), norm2(test_density(1:no))
                    print*, "testing o2p+p2o+o2p"
                    print*, norm2(abs(rhop)), norm2(abs(test_rhop))
                    print*,
                end do
            end do
        end do
print*, "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII"
        stop "ok boy energy_cproj"
    end subroutine energy_cproj


    subroutine o2p (ao, ap)
        use module_fft, only: fft2d
        use module_grid, only: grid
        use precision_kinds, only: dp
        implicit none
        real(dp), intent(in) :: ao(:) ! array of orientations
        complex(dp), intent(out) :: ap(:) ! array of projections
        integer :: mmax, itheta, iphi, ipsi, m, mup, mu, molrotsymorder, p, ip, io
        complex(dp), allocatable :: proj_theta(:,:,:)
        print*, "début de o2p"
        allocate (proj_theta(grid%ntheta,0:grid%mmax/grid%molrotsymorder,0:grid%mmax))
        print*, "allocation de proj_theta ok"
        mmax=grid%mmax
        molrotsymorder=grid%molrotsymorder
        do itheta=1,grid%ntheta
            do iphi=1,grid%nphi
                do ipsi=1,grid%npsi
                    io=grid%indo(itheta,iphi,ipsi)
                    print*,"io=",io
                    fft2d%in_forward(ipsi,iphi) = ao(io)
                end do
            end do
            print*, "boucle de remplissage de fft2d%in_forward ok"
            if (fft2d%isok) then
                call dfftw_execute (fft2d%plan_forward)
            else
                print*, "fft2d%is NOT ok in o2p"
                error stop "in o2p < energy_cproj"
            end if
            print*, "fftw execute ok"
            proj_theta(itheta,0:mmax/molrotsymorder,0:mmax)   = conjg( fft2d%out_forward(:,1:mmax+1) )
            proj_theta(itheta,0:mmax/molrotsymorder,-mmax:-1) = conjg( fft2d%out_forward(:,mmax+2:) )
            print*, "remplissage de proj_theta ok"
        end do
        ap=cmplx(0,0)
        do m=0,mmax
            do mup=-m,m
                do mu=0,m,molrotsymorder
                    p=indp(m,mup,mu)
                    if (p>size(ap)) error stop "in o2p, p>size(ap)"
                    ap(p)=sum(proj_theta(:,mu,mup)*tharm_sph(:,ip)*grid%wtheta(:))*fm(m)
                end do
            end do
        end do
        print*, "remplissage de ap ok"
        stop "fin de o2p"
    end subroutine o2p

    subroutine p2o (ap, ao)
        use module_fft, only: fft2d
        use module_grid, only: grid
        use precision_kinds, only: dp
        implicit none
        real(dp), intent(out) :: ao(:)
        complex(dp),intent(in) :: ap(:)
        complex(dp), allocatable, save :: proj_theta(:,:,:)
        integer :: itheta, iphi, ipsi, mup, mu, m, p, io
        if (.not. allocated(proj_theta)) then
            allocate (proj_theta(1:grid%ntheta,0:grid%mmax/grid%molrotsymorder,-grid%mmax:grid%mmax))
        end if
        proj_theta=complex(0._dp,0._dp)
        do itheta=1,grid%ntheta
            do mup=-grid%mmax,grid%mmax
                do mu=0,grid%mmax/grid%molrotsymorder
                    do m=max(abs(mup),abs(mu)), grid%mmax
                        p=indp(m,mup,mu)
                        proj_theta(itheta,mu,mup)=proj_theta(itheta,mu,mup)+ap(p)*tharm_sph(itheta,p)*fm(m)
                    end do
                end do
            end do
        end do
        do itheta=1,grid%ntheta
            fft2d%in_backward(:,1:grid%mmax+1) = conjg( proj_theta(itheta,:,0:grid%mmax))
            fft2d%in_backward(:,grid%mmax+2:) =  conjg( proj_theta(itheta,:,-grid%mmax:-1) )
            call dfftw_execute (fft2d%plan_backward)
            do iphi=1,grid%nphi
                do ipsi=1,grid%npsi
                    io=grid%indo(itheta,iphi,ipsi)
                    ao(io)=fft2d%out_backward(ipsi,iphi)*real(grid%nphi*grid%npsi,dp)
                end do
            end do
        end do
    end subroutine p2o





    PURE FUNCTION harm_sph( m, mu, mup, beta ) ! Luc's luc72p143
        !
        ! R^m_{mu,mup}(\beta)
        ! beta is the angle in radian
        !
        use mathematica, only: fact
        implicit none
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



end module module_energy_cproj

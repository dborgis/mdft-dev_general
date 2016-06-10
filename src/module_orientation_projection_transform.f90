module module_orientation_projection_transform

    use iso_c_binding, only: C_PTR, C_INT, C_INT32_T, C_INTPTR_T, C_DOUBLE_COMPLEX, C_DOUBLE, C_FUNPTR, C_SIZE_T, C_FLOAT, &
                             C_FLOAT_COMPLEX, C_CHAR
    use precision_kinds, only: dp
    use module_grid, only: grid
    use module_solvent, only: solvent
    use module_rotation, only: rotation_matrix_between_complex_spherical_harmonics_lu
    use module_wigner_d, only: wigner_small_d

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
    complex(dp), parameter :: zeroc=(0._dp, 0._dp)

    complex(dp), allocatable, protected :: foo_theta_mu_mup(:,:,:)


    type :: p3_type
        real(dp), allocatable :: wigner_small_d(:,:) ! tabulation des harmoniques sphériques r(m,mup,mu,theta) en un tableau r(itheta,p)
        integer, allocatable :: p(:,:,:) ! index of the projection corresponding to m, mup, mu
        integer, allocatable :: m(:) ! m for projection 1 to np
        integer, allocatable :: mup(:) ! mup for projection 1 to np. mup corresponds to phi
        integer, allocatable :: mu(:) ! mu for projection 1 to np. mu corresponds to psi
!        complex(dp), allocatable :: foo_q(:) ! foo (:) is a temporary array of size np
!        complex(dp), allocatable :: foo_mq(:) ! foo (:) is a temporary array of size np
    end type p3_type
    type (p3_type), protected :: p3


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
    type(  fft2d_type ), save, protected :: fft2d
    type( ifft2d_type ), save, protected :: ifft2d


    integer :: i, p
    integer :: nx, ny, nz, np, no, ntheta, nphi, npsi, mmax, mrso
    integer :: m, mu, mup, ip, ipsi, iphi, itheta
    real(dp), allocatable :: theta(:), wtheta(:)
    real(dp), parameter :: fm(0:7) = [( sqrt(real(2*m+1,dp)) ,m=0,7 )]
    logical :: is_init = .false.


    public :: angl2proj, proj2angl

contains

    subroutine init
      use omp_lib
      use iso_c_binding, only: c_ptr
      use precision_kinds, only: dp
      use module_grid, only: grid
      implicit none


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

      if (.not.allocated(theta)) allocate (theta(grid%ntheta))
      if (.not.allocated(wtheta)) allocate (wtheta(grid%ntheta))
      if (.not.allocated(fft2d%in)) allocate (fft2d%in(npsi,nphi))
      if (.not.allocated(fft2d%out)) allocate (fft2d%out(npsi/2+1,nphi)) ! pay attention, for fft2d we have psi as the first index, while it is the second index everywhere else.
      if (.not.allocated(ifft2d%in)) allocate (ifft2d%in(npsi/2+1,nphi)) ! this is because of our choice of doing half of mu thanks to hermitian symetry
      if (.not.allocated(ifft2d%out)) allocate (ifft2d%out(npsi,nphi))
      if (.not. allocated( foo_theta_mu_mup) ) allocate ( foo_theta_mu_mup(1:ntheta,0:mmax/mrso,-mmax:mmax))


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

        if (.not. fft2d%isalreadyplanned) then
          select case(dp)
          case(c_double)
            call dfftw_plan_dft_r2c_2d(  fft2d%plan, npsi, nphi,  fft2d%in,  fft2d%out, FFTW_EXHAUSTIVE ) ! npsi est en premier indice
            call dfftw_plan_dft_c2r_2d(  ifft2d%plan, npsi, nphi, ifft2d%in, ifft2d%out, FFTW_EXHAUSTIVE )
          case(c_float)
            call sfftw_plan_dft_r2c_2d(  fft2d%plan, npsi, nphi,  fft2d%in,  fft2d%out, FFTW_EXHAUSTIVE ) ! npsi est en premier indice
            call sfftw_plan_dft_c2r_2d(  ifft2d%plan, npsi, nphi, ifft2d%in, ifft2d%out, FFTW_EXHAUSTIVE )
          end select
          ! call dfftw_plan_dft_2d (fft2d_c%plan, npsi, nphi, fft2d_c%in, fft2d_c%out, FFTW_BACKWARD, FFTW_EXHAUSTIVE)
          fft2d%isalreadyplanned =.true.
          ifft2d%isalreadyplanned =.true.
        end if

      is_init = .true.

    end subroutine init


    function angl2proj (foo_o) result (foo_p)
        implicit none
        real(dp), intent(in) :: foo_o(:) ! orientations from 1 to no
        complex(dp) :: foo_p(1:grid%np) ! orientations from 1 to no
        integer :: itheta, iphi, ipsi, m, mup, mu, ip

        if(is_init .eqv. .false.) then
          call init
        end if

        foo_theta_mu_mup = zeroc
        do itheta=1,ntheta
            do iphi=1,nphi
                do ipsi=1,npsi
                    fft2d%in(ipsi,iphi) = foo_o(grid%indo(itheta,iphi,ipsi))
                end do
            end do
            select case(dp)
            case(c_double)
              call dfftw_execute (fft2d%plan)
            case(c_float)
              call sfftw_execute (fft2d%plan)
            end select
            foo_theta_mu_mup(itheta,0:mmax/mrso,0:mmax)   = CONJG( fft2d%out(:,1:mmax+1) )/real(nphi*npsi,dp)
            foo_theta_mu_mup(itheta,0:mmax/mrso,-mmax:-1) = CONJG( fft2d%out(:,mmax+2:) )/real(nphi*npsi,dp)
        end do
        foo_p = zeroc
        do m=0,mmax
            do mup=-m,m
                do mu=0,m,mrso
                    ip = p3%p( m,mup,mu/mrso )
                    foo_p(ip)= sum(foo_theta_mu_mup(:,mu/mrso,mup)*p3%wigner_small_d(:,ip)*wtheta(:))*fm(m)
                end do
            end do
        end do
    end function angl2proj


    function proj2angl (foo_p) result (foo_o)
        implicit none
        complex(dp), intent(in) :: foo_p(:) ! np
        real(dp) :: foo_o (1:grid%no) ! np

        if(is_init .eqv. .false.) then
          call init
        end if

        foo_theta_mu_mup = zeroc

        foo_o = 0._dp
        do mup=-mmax,mmax
            do mu=0,mmax,mrso
                do m= max(abs(mup),abs(mu)), mmax
                    ip=p3%p(m,mup,mu/mrso)
                    do itheta=1,ntheta
                        foo_theta_mu_mup(itheta,mu/mrso,mup) = foo_theta_mu_mup(itheta,mu/mrso,mup)&
                        +foo_p(ip)*p3%wigner_small_d(itheta,ip)*fm(m)
                    end do
                end do
            end do
        end do
        do itheta=1,ntheta
            ifft2d%in(:,1:mmax+1) = CONJG( foo_theta_mu_mup(itheta,0:mmax/mrso,0:mmax)   )
            ifft2d%in(:,mmax+2:)  = CONJG( foo_theta_mu_mup(itheta,0:mmax/mrso,-mmax:-1) )
            select case(dp)
            case(c_double)
              call dfftw_execute( ifft2d%plan )
            case(c_float)
              call sfftw_execute( ifft2d%plan )
            end select
            do iphi=1,nphi
                do ipsi=1,npsi
                    foo_o (grid%indo(itheta,iphi,ipsi)) = ifft2d%out(ipsi,iphi)
                end do
            end do
        end do
    end function proj2angl


end module module_orientation_projection_transform

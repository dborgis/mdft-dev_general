!This routine works properly for a single solvent or a mixture of it, but is more RAM consuming thatn the single solvent one
!module_energy_cproj_mrso, thus you should use it only with more than 1 solvent. This choice is done at initialization so 
!users should not worry about it. 
!TLDR: If you are using only one solvent then you are not using that routine but the one in module_energy_cproj_mrso.f90

module module_energy_cproj_mrso_mixture

    use iso_c_binding, only: C_PTR, C_INT, C_INT32_T, C_INTPTR_T, C_DOUBLE_COMPLEX, C_DOUBLE, C_FUNPTR, C_SIZE_T, C_FLOAT, &
                             C_FLOAT_COMPLEX, C_CHAR
    use precision_kinds, only: dp
    use module_grid, only: grid
    use module_solvent, only: solvent
    use module_rotation, only: rotation_matrix_between_complex_spherical_harmonics_lu, init_module_rotation => init

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
    type(c_type),protected,allocatable :: c(:)

    complex(dp), allocatable, protected :: deltarho_p(:,:,:,:,:) ! deltarho_p(np,nx,ny,nz)
    complex(dp), allocatable:: gammatmp(:,:,:,:,:) ! deltarho_p(np,nx,ny,nz)
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
    complex(dp), allocatable :: c3d(:,:,:), buf(:,:,:)
    !complex(dp), allocatable, protected :: R(:,:,:) ! Table of generalized spherical harmonics of m, mup, mu

    public :: energy_cproj_mrso_mixture


contains

    subroutine energy_cproj_mrso_mixture (ff, df, print_timers)

        use omp_lib
        use precision_kinds, only: dp
        use module_grid, only: grid
        use module_thermo, only: thermo
        use module_orientation_projection_transform, only: angl2proj, proj2angl, init_module_orientation_projection_transform => init
        use module_wigner_d, only: wigner_small_d

        implicit none

        real(dp), intent(out) :: ff
        real(dp), contiguous, intent(inout), optional :: df(:,:,:,:,:) ! x y z o s
        logical, intent(in), optional :: print_timers
        real(dp) :: dv, kT
        logical :: q_eq_mq
        integer :: ix, iy, iz, ix_q, iy_q, iz_q, ix_mq, iy_mq, iz_mq, ip
        integer :: nx, ny, nz, np, no, ns, ntheta, nphi, npsi, mmax, mrso
        integer :: m, n, mu, nu, khi, mup, iq, mu2, nu2, p, i, ia, s,s2
        real(dp) :: q(3), lx, ly, lz, rho0
        real(dp) :: theta(grid%ntheta), wtheta(grid%ntheta)
        logical, allocatable :: gamma_p_isok(:,:,:,:,:)
        real :: time(20)
        real(dp) :: vexc(grid%no)
        real :: total_time_in_subroutine
        complex(dp) :: deltarho_p_q_loc, deltarho_p_mq_loc
        complex(dp), parameter :: zeroc=(0._dp, 0._dp)
        real(dp), allocatable :: o(:)
        integer :: ierr
        real(dp) :: effectiveiq, alpha
        complex(dp), allocatable :: ceff(:)
        real(dp) :: prefactor
        real(dp), parameter :: fourpisq = 4._dp*acos(-1._dp)**2
        call cpu_time (time(1))

        ntheta = grid%ntheta
        nphi   = grid%nphi
        npsi   = grid%npsi
        mmax   = grid%mmax
        ! In the case of C∞v or any symetry with ∞, the user enters 0 for ∞ in the input file in molrotsymorder.
        ! To ease the loops over all mu=-m/mrso,m/mrso, it is easier to make it large, like 100, so that m/mrso:=0.
        ! With this trick, we have the same loops as before, with only one value of mu, 0.
        mrso = grid%molrotsymorder
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
    
        if (.not. allocated(c)) then
          allocate(c(solvent(1)%nspec))
        end if

        if (.not. allocated (deltarho_p) ) then
            allocate (deltarho_p(np,nx,ny,nz,size(solvent)) ,source=zeroc, stat=ierr)
            if (ierr/=0) PRINT*,"Allocate deltarho_p returns error ",ierr            
        end if
        if (.not. allocated (gammatmp) ) then
            allocate (gammatmp(np,nx,ny,nz,size(solvent)) ,source=zeroc, stat=ierr)
            if (ierr/=0) PRINT*,"Allocate deltarho_p returns error ",ierr            
        end if


        if( .not. allocated (gamma_p_isok) ) then 
           allocate (gamma_p_isok(nx,ny,nz,size(solvent),size(solvent)), source=.false., stat=ierr)
           if( ierr/=0) PRINT*,"Allocate gamma_p_isok returns error ", ierr
        end if

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

            !
            ! Prepare all callings to spherical harmonics generations
            !
            call init_module_rotation
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
        theta  = grid%thetaofntheta
        wtheta = grid%wthetaofntheta

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
            allocate( p3%wigner_small_d(1:ntheta, 1:np) ,source=0._dp)
            do p = 1, np
                m = p3%m(p)
                mup = p3%mup(p)
                mu = p3%mu(p)
                do i=1,ntheta
                    ! Pour chaque theta, calcule la fonction de Wigner-d correspondant à toutes les projections avec la méthode de Wigner.
                    p3%wigner_small_d(i,p) = wigner_small_d(m,mup,mu,theta(i))
                end do
            end do
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

        !
        ! 2/ ON PROJETTE delta_rho
        !
        call init_module_orientation_projection_transform
          
        !$omp parallel private (iz, iy, ix, o )
        allocate( o(no), stat=ierr)
        if (ierr/=0) PRINT*,"Allocate o returns error ",ierr
        !$omp do
        do s=1,size(solvent)
          do iz=1,nz
             do iy=1,ny
                do ix=1,nx
                   o = solvent(s)%rho0*(solvent(s)%xi(:,ix,iy,iz)**2 -1._dp)
                   call angl2proj( o, deltarho_p(:,ix,iy,iz,s) )
                end do
             end do
          end do
        end do
        !$omp end do
        deallocate(o)
        !$omp end parallel


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
            
        !$omp parallel private ( ip, buf )
        allocate(buf(nx,ny,nz),source=zeroC, stat=ierr)
        if (ierr/=0) PRINT*,"Allocate buf returns error ",ierr

        select case(dp);
        case(c_double)
            !$omp do
            do s=1,size(solvent)
              do ip = 1, np
                  buf=deltarho_p(ip,:,:,:,s)
                  call dfftw_execute_dft(fft%plan3dp, buf, buf)
                  deltarho_p(ip,:,:,:,s)=buf
              end do
            end do
            !$omp end do
        case(c_float)
            !$omp do
            do s=1,size(solvent)
              do ip=1,np
                  buf=deltarho_p(ip,:,:,:,s)
                  call sfftw_execute_dft(fft%plan3dp, buf, buf)
                  deltarho_p(ip,:,:,:,s)=buf
              end do
            end do
           !$omp end do
        end select
        deallocate(buf)
        !$omp end parallel

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
        !Guillaume: Here we read the c's, in principle we should read the pure solvent dcf and the mixture dcf
        !the mixture dcf are not available yet thus i code with the assumpution that c_ij=c_i and i read only the 
        !pure c function. This should be of course change to do proper mixture
        do s=1,size(solvent)
          if (.not.c(s)%isok) then
              block
                  use module_read_c_luc, only: read_c_luc
                  real(dp) :: qmaxnecessary
                  qmaxnecessary = norm2([maxval(grid%kx(1:nx)), maxval(grid%ky(1:ny)), maxval(grid%kz(1:nz/2+1))])
                  call read_c_luc(s,c(s)%mnmunukhi_q,mmax,mrso,qmaxnecessary,c(s)%np,c(s)%nq,c(s)%dq,c(s)%m,c(s)%n,c(s)%mu,c(s)%nu,c(s)%khi,c(s)%ip)
                  ! The c(m,n,mu,nu,khi) that we read from Luc had 2 drawbacks:
                  ! - we read it the way Luc give it, not the optimal way for our loops
                  ! - it does not contain our use of the symetries, and here we decide to have mu>0.
                  ! We could have chosen, mu'>0 or to keep all mu and mu' but work for half q-points.
                  ! That's our choice, and we must thus change the representation of c() to make it fit our purpose.
                  ! That's what we do below.
                  block
                      integer :: np_new, i
                      integer, allocatable :: m_new(:), n_new(:), mu_new(:), nu_new(:), khi_new(:), ip_new(:,:,:,:,:)
                      complex(dp), allocatable :: cnu2nmu2khim_q_new(:,:)
                      np_new = 0
                      do m=0,mmax; do khi=-m,m; do mu2=0,m/mrso; do n=abs(khi),mmax; do nu2=-n/mrso,n/mrso
                          np_new = np_new +1
                      end do; end do; end do; end do; end do
                      !Guillaume for now it is assume that nq are identical for
                      !all solvent and that the number of line of the c files
                      !coincid this is obviously not the case for interesting
                      !stuff but it will ease the coding for water in water
                      !This should be of course change to do proper mixture
                      allocate( cnu2nmu2khim_q_new(np_new,c(1)%nq) )
                      allocate( m_new(np_new) )
                      allocate( n_new(np_new) )
                      allocate( mu_new(np_new) )
                      allocate( nu_new(np_new) )
                      allocate( khi_new(np_new) )
                      allocate( ip_new(-mmax:mmax,0:mmax,0:mmax/mrso,-mmax:mmax,-mmax:mmax), source=-huge(1) )
                      i=0
                      do m=0,mmax; do khi=-m,m; do mu2=0,m/mrso; do n=abs(khi),mmax; do nu2=-n/mrso,n/mrso
                          i=i+1
                          m_new(i)=m
                          n_new(i)=n
                          mu_new(i)=mu2*mrso
                          nu_new(i)=nu2*mrso
                          khi_new(i)=khi
                          ip_new(nu2,n,mu2,khi,m)=i
                          cnu2nmu2khim_q_new(i,:) = c(s)%mnmunukhi_q(c(s)%ip(m,n,mu2,nu2,khi),:)
                      end do; end do; end do; end do; end do
                      deallocate( c(s)%m,c(s)%n,c(s)%mu,c(s)%nu,c(s)%khi,c(s)%ip, c(s)%mnmunukhi_q)
                      c(s)%np = np_new
                      allocate( c(s)%m(c(s)%np),c(s)%n(c(s)%np),c(s)%mu(c(s)%np),c(s)%nu(c(s)%np),c(s)%khi(c(s)%np),c(s)%ip(-mmax:mmax,0:mmax,0:mmax/mrso,-mmax:mmax,-mmax:mmax) )
                      allocate( c(s)%mnmunukhi_q(c(s)%np,c(s)%nq), source=cnu2nmu2khim_q_new )
                      deallocate(cnu2nmu2khim_q_new)
                      c(s)%m = m_new
                      c(s)%n = n_new
                      c(s)%mu = mu_new
                      c(s)%nu = nu_new
                      c(s)%khi = khi_new
                      c(s)%ip = ip_new
                  end block
              end block
              c(s)%mnmunukhi_q=conjg(c(s)%mnmunukhi_q) ! this is strange, but certainly due to some error or misunderstanding with Luc. It does not appear in Luc's document.
              !
              ! Move prefactors of MOZ inside the direct correlation function so that one does not need to compute, for instance,
              !  (-1)**(khi+nu) inside the inner loop of MOZ
              !
              do ip=1,c(s)%np
                  nu = c(s)%nu(ip)
                  if( nu<0 ) then
                      khi = c(s)%khi(ip)
                      c(s)%mnmunukhi_q(ip,:) = (-1)**(khi+nu) *c(s)%mnmunukhi_q(ip,:)
                  else
                      n = c(s)%n(ip)
                      c(s)%mnmunukhi_q(ip,:) = (-1)**(n) *c(s)%mnmunukhi_q(ip,:)
                  end if
              end do
              c(s)%isok=.true.
          end if
        END DO
        call cpu_time (time(9))

        !
        ! For all vectors q and -q handled simultaneously.
        ! ix_q,iy_q,iz_q are the coordinates of vector q, while ix_mq,iy_mq_iz_mq are those of vector -q
        !
        gamma_p_isok = .false.

        block
          integer ::  m, khi, mu2, ia
          complex(dp),allocatable,dimension(:,:) :: ceff
          real(dp) :: effectiveiq, alpha
          complex(dp) :: R(0:mmax,-mmax:mmax,-mmax:mmax)
          complex(dp), allocatable :: deltarho_p_q(:)
          complex(dp), allocatable :: deltarho_p_mq(:)
          complex(dp), allocatable :: gamma_p_q(:)
          complex(dp), allocatable :: gamma_p_mq(:)
         
          !Guillaume: Again, I emphasize that this should be changed when we want to use mixture c
          allocate(ceff(size(c),maxval(c(:)%np)))
          !$omp parallel private(iz_q, iy_q, ix_q, iz_mq, iy_mq, ix_mq, q, R, q_eq_mq, gamma_p_q, gamma_p_mq, m, khi, mu2, deltarho_p_q, deltarho_p_mq, deltarho_p_q_loc, deltarho_p_mq_loc, ia, ip, effectiveiq, iq, alpha, ceff)
          allocate (deltarho_p_q(np) ,source=zeroc, stat=ierr)
          if (ierr/=0) PRINT*,"Allocate deltarho_p_q returns error ",ierr
          allocate (deltarho_p_mq(np) ,source=zeroc, stat=ierr)
          if (ierr/=0) PRINT*,"Allocate deltarho_p_mq returns error ",ierr
          allocate (gamma_p_q(np), source=zeroc, stat=ierr)
          if (ierr/=0) PRINT*,"Allocate gamma_p_q returns error ",ierr
          allocate (gamma_p_mq(np), source=zeroc, stat=ierr)
          if (ierr/=0) PRINT*,"Allocate gamma_p_mq returns error ",ierr          
          !This is a temporary array of the same size than deltarho_p, it is used to do not overwrite deltarho_p
          !it is deallocated just after the big loop
          gammatmp(np,nx,ny,nz,size(solvent))=zeroC
          
          do s=1,size(solvent)
            do s2=1,size(solvent)
          !$omp do
              do iz_q=1,nz/2+1
                 iz_mq = grid%iz_mq(iz_q)
                 q(3) = grid%kz(iz_q) ! cartesian coordinates of vector q in lab frame

                 do iy_q=1,ny
                    iy_mq = grid%iy_mq(iy_q)
                    q(2) = grid%ky(iy_q) ! cartesian coordinates of vector q in lab frame
                    do ix_q=1,nx
                        ix_mq = grid%ix_mq(ix_q)
                        q(1) = grid%kx(ix_q) ! cartesian coordinates of vector q in lab frame

                        !
                        ! gamma_p_isok is a logical array. If gamma(ix_q,iy_q,iz_q) has already been calculated, it is .true.
                        !
                        if ( gamma_p_isok(ix_q,iy_q,iz_q,s,s2) .and. gamma_p_isok(ix_mq, iy_mq, iz_mq,s,s2) ) cycle

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
                        where( abs(real(R)) < epsilon(1.) )  R = cmplx( 0., aimag(R) )
                        where( abs(aimag(R)) < epsilon(1.) ) R = cmplx( real(R), 0.)

                        ! Eq. 1.23 We don't need to compute gshrot for -q since there are symetries between R(q) and R(-q).
                        ! Thus, we do q and -q at the same time. That's the most important point in doing all q but half of mu.
                        ! Lu decided to do all mu but half of q in her code

                        !
                        ! Rotation to molecular (q) frame
                        !

                        !  on  a       deltarho_p_q(m,khi,mu2) =  sum/mup  @   gamma_p_q(m,mup,mu2) * R(m,mup,khi)
                        !=>              gamma_p_q(mup,m,mu2) * R(mup,m,khi)
                        
                        !Guillaume: Here we are computing gamma_ij with i j beeing possible 2 diffrent species
                        !But for now gamma_p_q hold DeltaRho
                        gamma_p_q  = deltarho_p(:,ix_q,iy_q,iz_q,s) ! this a temporary array
                        gamma_p_mq = deltarho_p(:,ix_mq,iy_mq,iz_mq,s) ! this a temporary array
                        ip=0
                        do m=0,mmax
                            do khi=-m,m
                                do mu2=0,m/mrso
                                    ip=ip+1
                                    ! Equation 1.22
                                    deltarho_p_q_loc  = (0._dp,0._dp)
                                    deltarho_p_mq_loc = (0._dp,0._dp)
                                    do mup=-m,m
                                        deltarho_p_q_loc  = deltarho_p_q_loc  + gamma_p_q (p3%p(m,mup,mu2)) * R(m,mup,khi)
                                        deltarho_p_mq_loc = deltarho_p_mq_loc + gamma_p_mq(p3%p(m,mup,mu2)) * R(m,mup,khi)
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
                        ! iq = int( norm2(q) /c%dq +0.5) +1
                        ! ceff(:) = c%mnmunukhi_q(:,iq)

                        !Guillaume: Since we used s to describe solvent 1 in rho, we must use s2 to descrobe solvent 2 in c, to get gamma_ss2
                        effectiveiq = norm2(q)/c(s2)%dq +1  ! norm(q)/dq is in [0,n] while our iq should be in [1,n+1]. Thus, add +1.
                        iq = int(effectiveiq) ! the lower bound. The upper bound is iq+1

                        alpha = effectiveiq - iq ! linear interpolation    y=alpha*upperbound + (1-alpha)*lowerbound
                        ceff(s2,:) =         alpha  * c(s2)%mnmunukhi_q(:,iq+1) &
                             + (1._dp-alpha) * c(s2)%mnmunukhi_q(:,iq)

                        !
                        ! Ornstein-Zernike in the molecular frame
                        ! We do OZ for q and -q at the same time
                        !

                        ! ! Prepare deltarho_p_q for MOZ
                        ! block
                        !     integer :: ip
                        !     complex(dp) :: dummy4swap
                        !     do ip=1,np
                        !         if( p3%mu(ip)>=0 ) then
                        !             dummy4swap        = conjg(deltarho_p_q(ip))
                        !             deltarho_p_q(ip)  = conjg(deltarho_p_mq(ip))
                        !             deltarho_p_mq(ip) = dummy4swap
                        !         else
                        !             stop "coucou"
                        !         end if
                        !     end do
                        ! end block

                        gamma_p_q  = zeroc
                        gamma_p_mq = zeroc !1:np
                        ia = 0
                        do ip = 1, np
                            m   = p3%m(ip)       ! m=0,mmax
                            khi = p3%mup(ip)     ! khi=-m,m
                            mu2 = p3%mu(ip)/mrso ! mu2=-m/mrso,m/mrso
                            do n = abs(khi), mmax
                                !ia = ia +1
                                !select case (n)
                                !case(0, 1) ! nu2 == 0
                                !    gamma_p_q(ip)  = gamma_p_q(ip)    + ceff(ia  ) *conjg(deltarho_p_mq(p3%p(n,khi,0)))
                                !    gamma_p_mq(ip) = gamma_p_mq(ip)   + ceff(ia  ) *conjg(deltarho_p_q (p3%p(n,khi,0)))
                                !case(2, 3) ! nu2 = -1,0,1
                                !    gamma_p_q(ip)  = gamma_p_q(ip)    + ceff(ia  ) *      deltarho_p_q (p3%p(n,khi,1))  &
                                !                                      + ceff(ia+1) *conjg(deltarho_p_mq(p3%p(n,khi,0))) &
                                !                                      + ceff(ia+2) *conjg(deltarho_p_mq(p3%p(n,khi,1)))
                                !    gamma_p_mq(ip) = gamma_p_mq(ip)   + ceff(ia  ) *      deltarho_p_mq(p3%p(n,khi,1))  &
                                !                                      + ceff(ia+1) *conjg(deltarho_p_q (p3%p(n,khi,0))) &
                                !                                      + ceff(ia+2) *conjg(deltarho_p_q (p3%p(n,khi,1)))
                                !    ia = ia +2
                                !case(4, 5) ! nu2 = -2,-1,0,1,2
                                !    gamma_p_q(ip)  = gamma_p_q(ip)    + ceff(ia  ) *      deltarho_p_q (p3%p(n,khi,2))  &
                                !                                      + ceff(ia+1) *      deltarho_p_q (p3%p(n,khi,1))  &
                                !                                      + ceff(ia+2) *conjg(deltarho_p_mq(p3%p(n,khi,0))) &
                                !                                      + ceff(ia+3) *conjg(deltarho_p_mq(p3%p(n,khi,1))) &
                                !                                      + ceff(ia+4) *conjg(deltarho_p_mq(p3%p(n,khi,2)))
                                !    gamma_p_mq(ip) = gamma_p_mq(ip)   + ceff(ia  ) *      deltarho_p_mq(p3%p(n,khi,2))  &
                                !                                      + ceff(ia+1) *      deltarho_p_mq(p3%p(n,khi,1))  &
                                !                                      + ceff(ia+2) *conjg(deltarho_p_q (p3%p(n,khi,0))) &
                                !                                      + ceff(ia+3) *conjg(deltarho_p_q (p3%p(n,khi,1))) &
                                !                                      + ceff(ia+4) *conjg(deltarho_p_q (p3%p(n,khi,2)))
                                !    ia = ia +4
                                !end select
                                !Guillaume: This does not appear explicitely here but deltarho_p_q and deltarho_p_mq are s depandant, thus 
                                !gamma_p_q and gamma_p_mq correspond indeed to gamma_ss2
                                do nu2=-n/mrso ,n/mrso 
                                  ia=ia+1
                                    if (nu2<0) then
                                      gamma_p_q(ip)  = gamma_p_q(ip)    + ceff(s2,ia  ) * deltarho_p_q (p3%p(n,khi,abs(nu2)))  
                                      gamma_p_mq(ip) = gamma_p_mq(ip)   + ceff(s2,ia  ) * deltarho_p_mq(p3%p(n,khi,abs(nu2)))  
                                    else
                                      gamma_p_q(ip)  = gamma_p_q(ip)    + ceff(s2,ia) *conjg(deltarho_p_mq(p3%p(n,khi,nu2))) 
                                      gamma_p_mq(ip) = gamma_p_mq(ip)   + ceff(s2,ia) *conjg(deltarho_p_q (p3%p(n,khi,nu2))) 
                                    end if
                                end do
                            end do
                        end do

                        !
                        ! Rotation from molecular frame to fix frame
                        !
                        R = conjg(R) ! le passage retour au repaire fixe se fait avec simplement le conjugue complexe de l'harm sph generalisee
                        ! we use deltarho_p_q and deltarho_p_mq as temp arrays since they're not used after MOZ
                        

                        !
                        ! prevent underflow in gamma_p_q/mq * R if gamma_p is very low
                        !
                        where( abs(real(gamma_p_q))   < epsilon(1._dp) .and. abs(aimag(gamma_p_q))  < epsilon(1._dp) ) gamma_p_q  = cmplx( 0._dp, 0._dp )
                        where( abs(real(gamma_p_mq))  < epsilon(1._dp) .and. abs(aimag(gamma_p_mq)) < epsilon(1._dp) ) gamma_p_mq = cmplx( 0._dp, 0._dp )
                        where( abs(real(gamma_p_q))   < epsilon(1._dp) ) gamma_p_q  = cmplx( 0._dp, aimag(gamma_p_q) )
                        where( abs(real(gamma_p_mq))  < epsilon(1._dp) ) gamma_p_mq = cmplx( 0._dp, aimag(gamma_p_mq) )
                        where( abs(aimag(gamma_p_q))  < epsilon(1._dp) ) gamma_p_q  = cmplx( real(gamma_p_q), 0._dp )
                        where( abs(aimag(gamma_p_mq)) < epsilon(1._dp) ) gamma_p_mq = cmplx( real(gamma_p_mq), 0._dp )
                        

                        deltarho_p_q = (0._dp,0._dp)
                        deltarho_p_mq = (0._dp,0._dp)
                        ip=0
                        do m=0,mmax
                            do mup=-m,m
                                do mu2=0,m/mrso
                                    ip=ip+1
                                    ! Equation 1.22
                                    do khi=-m,m
                                        deltarho_p_q(ip)  = deltarho_p_q(ip)  + gamma_p_q (p3%p(m,khi,mu2)) *R(m,mup,khi)
                                        deltarho_p_mq(ip) = deltarho_p_mq(ip) + gamma_p_mq(p3%p(m,khi,mu2)) *R(m,mup,khi)
                                    end do
                                    deltarho_p_mq(ip) = deltarho_p_mq(ip) *(-1)**m
                                end do
                            end do
                        end do


                        !
                        ! Move the result for this given vector q to the big array containing all results.
                        ! First, for q,

                        !Guillaume:Here we shoud add all the contribtuion coming froms s2 in order to have a simple integration in ff
                        gammatmp(1:np, ix_q, iy_q, iz_q,s2) =gammatmp(1:np, ix_q, iy_q, iz_q,s2) + deltarho_p_q(1:np)
                        !gammatmp(1:np, ix_q, iy_q, iz_q,s) =deltarho_p_q(1:np)
                        !deltarho_p(1:np, ix_q, iy_q, iz_q,s) =deltarho_p_q(1:np)
                        !
                        ! Then, for -q. Again, pay attention to the singular mid-k point
                        !
                        if( q_eq_mq .and. (ix_q==nx/2+1 .or. iy_q==ny/2+1 .or. iz_q==nz/2+1)) then
                          !if (any(gammatmp(1:np, ix_mq, iy_mq, iz_mq,s)/=zeroC)) then
                          !  print*, "first gammatmp(1:np, ix_mq, iy_mq, iz_mq,s) is not zero",  ix_mq, iy_mq, iz_mq, q_eq_mq
                          !else
                          !  print*, "another one that actualy count in first", ix_mq, iy_mq, iz_mq, q_eq_mq
                          !end if
                          if (.not. q_eq_mq ) gammatmp(1:np, ix_mq, iy_mq, iz_mq,s2) = gammatmp(1:np, ix_mq, iy_mq, iz_mq,s2) + conjg(deltarho_p_mq(1:np))
                          !gammatmp(1:np, ix_mq, iy_mq, iz_mq,s) = conjg(deltarho_p_mq(1:np))
                          !deltarho_p(1:np, ix_mq, iy_mq, iz_mq,s) = conjg(deltarho_p_mq(1:np))
                        else
                          !if (any(gammatmp(1:np, ix_mq, iy_mq, iz_mq,s)/=zeroC)) then
                          !  print*, "second gammatmp(1:np, ix_mq, iy_mq, iz_mq,s) is not zero",  ix_mq, iy_mq, iz_mq,s, q_eq_mq
                          !else if (q_eq_mq) then
                          !  print*, "this q=-q counts sound weird"
                          !end if
                          if (.not. q_eq_mq )  gammatmp(1:np, ix_mq, iy_mq, iz_mq,s2) = gammatmp(1:np, ix_mq, iy_mq, iz_mq,s2) + deltarho_p_mq(1:np)
                          !gammatmp(1:np, ix_mq, iy_mq, iz_mq,s) =  deltarho_p_mq(1:np)
                          !deltarho_p(1:np, ix_mq, iy_mq, iz_mq,s) = deltarho_p_mq(1:np)
                        end if
                        !
                        ! And store you have already done the job
                        !
                        gamma_p_isok(ix_q,iy_q,iz_q,s,s2)=.true.
                        gamma_p_isok(ix_mq,iy_mq,iz_mq,s,s2)=.true.

                  end do
              end do
            end do
        !$omp end do
          end do!solvent 2 loop
        end do!solvent 1 loop
        ! deallocate : (not necessary)
        deallocate (deltarho_p_q)
        deallocate (deltarho_p_mq)
        deallocate (gamma_p_q)
        deallocate (gamma_p_mq)
        !$omp end parallel
        deltarho_p=gammatmp
        deallocate (gammatmp)
    end block



        call cpu_time(time(10))

        if ( any(gamma_p_isok .eqv. .false.) ) then
            error stop "not all gamma_p(projections,ix,iy,iz) have not been computed"
        end if
        deallocate(gamma_p_isok)
        call cpu_time(time(11))


        !
        ! FFT3D from Fourier space to real space
        !
        block
          real(dp) :: cst
            
            !$omp parallel private ( ip, cst, buf )
            allocate(buf(nx,ny,nz), stat=ierr)
            if (ierr/=0) PRINT*,"Allocate buf returns error ",ierr

            cst=1._dp/real(nx*ny*nz,dp)
            
            select case(dp)
            case(c_double)
               !$omp do
               do s=1,size(solvent)
                do ip=1,np
                   buf = deltarho_p(ip,:,:,:,s) 
                   call dfftw_execute_dft( fft%plan3dm, buf, buf)
                   deltarho_p(ip,:,:,:,s) = buf*cst
                end do
              end do
               !$omp end do
            case(c_float)
               !$omp do
               do s=1,size(solvent)
                do ip=1,np
                   buf = deltarho_p(ip,:,:,:,s)
                   call sfftw_execute_dft( fft%plan3dm, buf, buf)
                   deltarho_p(ip,:,:,:,s) = buf*cst
                end do
               end do
               !$omp end do
            end select
            deallocate(buf)
            !$omp end parallel
        end block

        call cpu_time(time(12))

        !
        ! Gather projections into gamma
        ! Note that gamma==df
        !

            if(present(df)) then
                ff=0._dp
                do s2=1,size(solvent)
                  prefactor = -kT*2.0_dp*fourpisq/mrso /solvent(s2)%n0*solvent(s2)%mole_fraction!/2.0! the division by n0 comes from Luc's normalization of c
                !$omp parallel private(iz, iy, ix, vexc) reduction(+:ff)
                !$omp do
                !do s=1,size(solvent)
                  do iz=1,nz
                      do iy=1,ny
                          do ix=1,nx
                              call proj2angl( deltarho_p(:,ix,iy,iz,s2), vexc)
                              vexc = prefactor*grid%w*vexc
                              ff = ff + sum((solvent(s2)%xi(:,ix,iy,iz)**2*solvent(s2)%rho0-solvent(s2)%rho0)*vexc)
                              df(:,ix,iy,iz,s2) = df(:,ix,iy,iz,s2) + 2._dp*solvent(s2)%rho0*solvent(s2)%xi(:,ix,iy,iz)*vexc
                          end do
                      end do
                  end do
                !$omp end do
                !$omp end parallel
                end do
                              !stop
                ff = ff*0.5_dp*dv
            else
                !Guillaume:I did not check mixture code for this loop 
                ff=0._dp
                !$omp parallel private(iz, iy, ix, vexc) reduction(+:ff)
                !$omp do
                do s2=1,size(solvent)
                  prefactor = -kT*2.0_dp*fourpisq/mrso /solvent(s2)%n0*solvent(s2)%mole_fraction!/2.0! the division by n0 comes from Luc's normalization of c
                  do iz=1,nz
                      do iy=1,ny
                          do ix=1,nx
                              call proj2angl( deltarho_p(:,ix,iy,iz,s), vexc)
                              vexc = prefactor*grid%w*vexc
                              ff = ff + sum((solvent(s2)%xi(:,ix,iy,iz)**2*solvent(s2)%rho0-solvent(s2)%rho0)*vexc)
                          end do
                      end do
                  end do
                end do
                !$omp end do
                !$omp end parallel
                ff = ff*0.5_dp*dv
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

    end subroutine energy_cproj_mrso_mixture

end module module_energy_cproj_mrso_mixture

module module_fastpoissonsolver

    use precision_kinds, only: dp
    implicit none
    private

    real(dp), allocatable, dimension(:,:,:), protected :: sourcedistrib ! charge distribution that generates the electrostatic potential
    real(dp), allocatable, dimension(:,:,:), protected :: phi ! phi is the electrostatic potential
    type :: psgrid_type
        integer :: n(3)
        real(dp) :: len(3)
    end type psgrid_type
    type (psgrid_type), protected :: psgrid

    public :: init_fastpoissonsolver

contains

    !===================================================================================================================================

    subroutine init_fastpoissonsolver
        use precision_kinds, only: dp
        use module_grid, only: grid
        use module_solute, only: soluteChargeDensityFromSoluteChargeCoordinates
        implicit none
        type :: myerror_type
            integer :: i
            character(180) :: msg
        end type
        type (myerror_type) :: er
        real(dp), parameter :: epsdp=epsilon(1._dp)

        !
        ! Check allocations and initializations
        !
        if (.not. grid%isinitiated) then
            print*, "Vous voulez utiliser le type grid dans module_fastpoissonsolver mais il n'est pas init"
            error stop
        end if

        ! Pour l'instant, psgrid est exactement de la dimension de grid. L'objectif ici est de pouvoir résoudre poisson sur une meilleure grid
        psgrid%n = grid%n_nodes
        psgrid%len  = grid%length

        allocate( sourcedistrib (psgrid%n(1),psgrid%n(2),psgrid%n(3)) ,SOURCE=0._dp, stat=er%i, errmsg=er%msg)
        if (er%i/=0) then
            print*, er%msg
            error stop "This problem arises in module fastPoissonSolver"
        end if

        allocate( phi (psgrid%n(1),psgrid%n(2),psgrid%n(3)) ,SOURCE=0._dp, stat=er%i, errmsg=er%msg)
        if (er%i/=0) then
            print*, er%msg
            error stop "This problem arises in module fastPoissonSolver"
        end if

        call soluteChargeDensityFromSoluteChargeCoordinates (psgrid%n, psgrid%len, sourcedistrib)

        if (all(abs(sourcedistrib)<=epsdp)) then
            phi = 0._dp
        end if

        call poissonSolver (psgrid%n, psgrid%len, sourcedistrib, phi)


        !
        ! Deallocate unnecessary stuff
        !
        if (allocated(sourcedistrib)) deallocate (sourcedistrib)
        if (allocated(phi)) deallocate (phi)

        ! at this point, one has Vext_q (old fashioned construction of Poisson potential and extrapolation to solvent sites outside grid nodes)
        ! and of solvent(s)%vext(i,j,k,o,p)%q, the new construction with multipole expansion of the solvent molecular charge density.

        ! Since vext_q is used later in the code, if one wants the newer construction, move solvent%vext%q to vext_q
        ! if (getinput%log('better_poisson_solver')) then
        !     do concurrent (s = 1:size(solvent))
        !         vext_q(:,:,:,:,s) = solvent(s)%vextq
        !     end do
        ! ! else ! this should be removed soon once we're used to this new construction
        ! !     print*,"BETTER POISSON SOLVER [OFF]"
        ! !     call vext_q_from_v_c (psgrid%n, psgrid%len, phi) ! defines Vext_q(i,j,k,o,p) that we used in the past for Poisson solver
        ! end if


    end subroutine init_fastpoissonsolver

    !===================================================================================================================================

    ! This SUBROUTINE uses the Gauss theorem associated to the expression of the derivation in Fourier transforms to get V(r).
    ! Field E(r)=-grad(V(r))
    ! Local expression of Gauss theorem : div(E(r))=sourcedistrib(r)/eps0
    ! => laplacien(V(r))=-sourcedistrib(r)/eps0
    ! => V(k)=sourcedistrib(k)/(esp0*k^2)
    ! FFT(V(k)) = V(r)

    SUBROUTINE poissonSolver (gridnode, gridlength, sourcedistrib, phi)

        use precision_kinds, only: dp, i4b
        use constants, only: qfact
        use module_solvent, only: solvent
        use module_input, only: getinput
        use module_grid, only: grid

        ! phi = electrostatic potential from charge density and poisson equation

        implicit none

        integer, intent(in) :: gridnode(3)
        real(dp), intent(in) :: gridlength(3)
        real(dp), intent(in) :: sourcedistrib (gridnode(1),gridnode(2),gridnode(3))
        real(dp), intent(out) :: phi (gridnode(1),gridnode(2),gridnode(3))
        complex(dp), allocatable, dimension(:,:,:) :: sourcedistrib_k, phi_k
        integer  :: i,j,k,m1,m2,m3,s,io
        real(dp) :: k2
        real(dp), allocatable :: fftw3inforward(:,:,:), fftw3outbackward(:,:,:)
        complex(dp), allocatable :: fftw3outforward(:,:,:), fftw3inbackward(:,:,:)
        integer(i4b) :: fpspf, fpspb ! fast poisson solver plan forward or backward
        complex(dp), parameter :: zeroc=complex(0._dp,0._dp)
        real(dp), parameter :: twopi  =2._dp*acos(-1._dp)
        real(dp), parameter :: fourpi =2._dp*twopi
        real(dp), parameter :: espdp=epsilon(1._dp)
        integer :: nx, ny, nz

        include "fftw3.f"

        nx = gridnode(1)
        ny = gridnode(2)
        nz = gridnode(3)


        ! allocate the arrays needed as input for fft (in_forward) or output for fft (out_forward)
        ! or needed as input for inverse fft (in_backward) etc.
        allocate ( fftw3inforward   (nx    , ny,nz))
        allocate ( fftw3outforward  (nx/2+1, ny,nz))
        allocate ( fftw3outbackward (nx    , ny,nz))
        allocate ( fftw3inbackward  (nx/2+1, ny,nz))
        ! prepare plans needed by fftw3
        call dfftw_plan_dft_r2c_3d (fpspf, nx, ny, nz, fftw3inforward, fftw3outforward, fftw_estimate)
        call dfftw_plan_dft_c2r_3d (fpspb, nx, ny, nz, fftw3inbackward, fftw3outbackward, fftw_estimate )
        ! note that since the fast poisson solver implies only 1 fft in each direction, it is useless to use fftw_measure or even
        ! more rigorous planning-flags. see http://www.fftw.org/doc/planner-flags.html



        if ( all(abs(sourcedistrib)<=epsilon(1._dp)) ) then
            phi = 0._dp
            return
        end if

        allocate( sourcedistrib_k (nx/2+1,gridnode(2),gridnode(3)) ,source=zeroc)
        allocate( phi_k (nx/2+1,gridnode(2),gridnode(3)) ,source=zeroc)

        ! fourier transform of the solute charge density
        fftw3inforward = sourcedistrib
        call dfftw_execute (fpspf)
        sourcedistrib_k = fftw3OutForward ! It is verified that at this point, FFT-1(sourcedistrib_k)/ (nfft1*nfft2*nfft3) = sourcedistrib
        ! FFT(Laplacian(V(r))) = FFT( - 4Pi charge density(r) ) in elecUnits = (ik)^2 V(k) = -4pi rho(k)
        ! V(k) = 4Pi rho(k) / k^2


        DO k = 1, gridnode(3)
            DO j = 1, gridnode(2)
                DO i = 1, nx/2+1

                    IF ( i<=nx/2 ) THEN
                        m1 = i-1
                    ELSE
                        m1 = i-1-nx
                    END IF

                    IF ( j<=gridnode(2)/2 ) THEN
                        m2 = j-1
                    ELSE
                        m2 = j-1-gridnode(2)
                    END IF

                    IF ( k<=gridnode(3)/2 ) THEN
                        m3 = k-1
                    ELSE
                        m3 = k-1-gridnode(3)
                    END IF

                    k2 = sum(   (twopi*real([m1,m2,m3],dp)/gridlength)**2   )

                    IF ( abs(k2) > epsilon(1._dp) ) THEN
                        phi_k(i,j,k) = sourcedistrib_k(i,j,k) * fourpi/k2 ! in electrostatic units : V=-4pi rho
                    ELSE
                        phi_k(i,j,k) = complex(0._dp,0._dp)
                    END IF

                END DO
            END DO
        END DO


        ! old construction
        if (getinput%log('better_poisson_solver', defaultvalue=.false.)) then
            ! get electrostatic potentiel in real space, that is the true Poisson potentiel. It is not solvent dependent
            fftw3InBackward = phi_k
            call dfftw_execute (fpspb)
            phi = qfact * fftw3OutBackward / REAL(PRODUCT(gridnode),dp) ! kJ/mol
        else
            ! new construction
            ! get multipolar electrostatic potential in real space. It is already solvent dependent here
            do s=1,solvent(1)%nspec
                if (.not. allocated(solvent(s)%sigma_k) ) then
                    call solvent(s)%init_chargedensity_molecularpolarization
                end if
                do io = 1, grid%no
                    fftw3InBackward = phi_k * conjg( solvent(s)%sigma_k(:,:,:,io) )
                    call dfftw_execute (fpspb)
                    solvent(s)%vextq(io,:,:,:) = qfact * fftw3OutBackward / real( product(gridnode) ,dp) ! kJ/mol
                end do
            end do
        end if

        call dfftw_destroy_plan (fpspf)
        call dfftw_destroy_plan (fpspb)

        deallocate( sourcedistrib_k, phi_k, fftw3InForward, fftw3OutForward, fftw3OutBackward, fftw3InBackward)

    END SUBROUTINE poissonSolver

    !===================================================================================================================================

    SUBROUTINE vext_q_from_v_c (gridnode, gridlen, phi)

        use precision_kinds     ,ONLY: dp
        use module_solute              ,ONLY: solute
        use module_solvent, only: solvent
        use module_grid, only: grid
        use constants           ,ONLY: qfact, zero
        use mathematica, only: TriLinearInterpolation, UTest_TrilinearInterpolation, UTest_floorNode,&
        floorNode, ceilingNode, UTest_ceilingNode, UTest_distToFloorNode, distToFloorNode, chop

        IMPLICIT NONE

        INTEGER, INTENT(IN) :: gridnode(3)
        REAL(dp), DIMENSION(gridnode(1),gridnode(2),gridnode(3)), INTENT(IN) :: phi
        REAL(dp), INTENT(IN) :: gridlen(3)
        INTEGER :: i, j, k, m, s, nfft(3), l(3), u(3), io
        real(dp), dimension(3,solvent(1)%nspec,solvent(1)%nsite,grid%no) :: xmod
        REAL(dp) :: vext_q_of_r_and_omega ! external potential for a given i,j,k,omega & psi.
        REAL(dp) :: r(3), cube(0:1,0:1,0:1), dl(3), x(3)
        TYPE :: testtype
            LOGICAL :: pb
            CHARACTER(180) :: msg
            real(dp) :: l(3),u(3),r(3),x(3)
        END TYPE
        TYPE(testtype) :: err

        nfft = grid%n_nodes
        dl = grid%dl

        IF(.NOT. ALLOCATED(solvent(1)%vextq)) STOP "solvent(1)%vextq should already be allocated in vext_q_from_v_c.f90"
        IF( ANY(solvent(1)%vextq/=0._dp) ) STOP "vext_q should be zero everywhere in vext_q_from_v_c.f90"

        IF ( ALL(solute%site%q == zero) .OR. ALL(solvent(1)%site%q == zero) ) RETURN

        ! Tabulate the cartesian coordinates of all solvent sites, for all molecular orientations, centered on any MDFT's grid node.
        do concurrent (s=1:size(solvent))
            DO CONCURRENT( io=1:grid%no, m=1:size(solvent(s)%site) )
                xmod(1,s,m,io) = DOT_PRODUCT( [grid%Rotxx(io),grid%Rotxy(io),grid%Rotxz(io)] , solvent(s)%site(m)%r )
                xmod(2,s,m,io) = DOT_PRODUCT( [grid%Rotyx(io),grid%Rotyy(io),grid%Rotyz(io)] , solvent(s)%site(m)%r )
                xmod(3,s,m,io) = DOT_PRODUCT( [grid%Rotzx(io),grid%Rotzy(io),grid%Rotzz(io)] , solvent(s)%site(m)%r )
            end do
        END DO

        CALL UTest_floorNode
        CALL UTest_ceilingNode
        CALL UTest_TrilinearInterpolation
        CALL UTest_distToFloorNode
        ! Compute external potential for each combination of solvent side and grid node and orientation
        err%pb=.FALSE. ! becomes TRUE if a problem is detected during execution.


        DO CONCURRENT( s=1:solvent(1)%nspec)
            do concurrent ( io=1:grid%no, i=1:nfft(1), j=1:nfft(2), k=1:nfft(3) )
                vext_q_of_r_and_omega = 0.0_dp

                DO CONCURRENT (m=1:solvent(1)%nsite, abs(solvent(s)%site(m)%q)>epsilon(1.0_dp))

                    x = (REAL([i,j,k],dp)-1.0_dp)*dl + [xmod(1,s,m,io),xmod(2,s,m,io),xmod(3,s,m,io)]! cartesian coordinate x of the solvent site m. May be outside the supercell.
                    x(1)=chop(x(1))
                    x(2)=chop(x(2))
                    x(3)=chop(x(3))

                    l = floorNode(gridnode,gridlen,x,.TRUE.) ! r should be in full cartesian coordinates between -infty and +infty
                    u = ceilingNode(gridnode,gridlen,x,.TRUE.)
                    r = distToFloorNode(gridnode,gridlen,x,.TRUE.) ! 0 <= distToFloorNode < 1

                    if ( ANY(r<0._dp) .or. ANY(r>1._dp) ) THEN
                        err%pb=.TRUE.
                        err%msg="Problem with r in vext_q_from_v_c.f90"
                        err%l = l
                        err%u = u
                        err%r = r
                        err%x = x
                    end if

                    if ( ANY(l<LBOUND(phi)) .or. ANY(l>UBOUND(phi)) ) THEN
                        err%pb=.TRUE.
                        err%msg="Problem with l in vext_q_from_v_c.f90"
                        err%l = l
                        err%u = u
                        err%r = r
                        err%x = x
                    END IF

                    if ( ANY(u<LBOUND(phi)) .or. ANY(u>UBOUND(phi)) ) THEN
                        err%pb=.TRUE.
                        err%msg="Problem with u in vext_q_from_v_c.f90"
                        err%l = l
                        err%u = u
                        err%r = r
                        err%x = x
                    end if

                    cube(0,0,0) = phi (l(1),l(2),l(3))
                    cube(1,0,0) = phi (u(1),l(2),l(3))
                    cube(0,1,0) = phi (l(1),u(2),l(3))
                    cube(0,0,1) = phi (l(1),l(2),u(3))
                    cube(1,1,0) = phi (u(1),u(2),l(3))
                    cube(1,0,1) = phi (u(1),l(2),u(3))
                    cube(0,1,1) = phi (l(1),u(2),u(3))
                    cube(1,1,1) = phi (u(1),u(2),u(3))

                    vext_q_of_r_and_omega = vext_q_of_r_and_omega + solvent(s)%site(m)%q * TriLinearInterpolation(cube,r)
                end do
                solvent(s)%vextq(io,i,j,k) = vext_q_of_r_and_omega

            end do
        end do
        if (err%pb) then
            print*,"msg:",err%msg
            print*,"Lxyz:",gridlen
            print*,"x:",err%x
            print*,"l/u min:",1,1,1
            print*,"l/u max:",gridnode
            print*,"l:",err%l
            print*,"u:",err%u
            print*,"r:",err%r
            stop
        end if

    end subroutine vext_q_from_v_c

end module module_fastpoissonsolver
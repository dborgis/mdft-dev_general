module fastPoissonSolver
    ! This module contains everything related to the fast poisson solver(s) of MDFT.

    use precision_kinds     ,only: dp, i2b
    use system              ,only: solute, solvent
    use module_grid, only: grid
    use external_potential  ,only: vext_q
    use input               ,only: getinput

    implicit none

    real(dp), allocatable, dimension(:,:,:), private :: soluteChargeDensity, Vpoisson ! TODO THIS ON FINE grid
    type :: PoissonGridType
        integer(i2b) :: n(3)
        real(dp) :: len(3)
    end type PoissonGridType
    type (PoissonGridType), private :: poissonSolverGrid
    private
    public :: init

contains

    !===================================================================================================================================

    subroutine init
        implicit none
        integer(i2b) :: i,j,k,o,p,s,n1,n2,n3
        type :: ertype
            integer :: i
            character(180) :: msg
        end type
        type (ertype) :: er
        poissonSolverGrid%n = grid%n_nodes
        poissonSolverGrid%len  = grid%length
        allocate( soluteChargeDensity (poissonSolverGrid%n(1),poissonSolverGrid%n(2),poissonSolverGrid%n(3)) ,SOURCE=0._dp, &
        stat=er%i, errmsg=er%msg)
        if (er%i/=0) then
            print*, er%msg
            stop "Cannot This problem arises in module fastPoissonSolver"
        end if
        allocate( Vpoisson            (poissonSolverGrid%n(1),poissonSolverGrid%n(2),poissonSolverGrid%n(3)) ,SOURCE=0._dp, &
        stat=er%i, errmsg=er%msg)
        if (er%i/=0) then
            print*, er%msg
            stop "Cannot This problem arises in module fastPoissonSolver"
        end if
        call soluteChargeDensityFromSoluteChargeCoordinates (poissonSolverGrid%n, poissonSolverGrid%len, soluteChargeDensity)
        call poissonSolver (poissonSolverGrid%n, poissonSolverGrid%len, soluteChargeDensity, Vpoisson)

        ! at this point, one has Vext_q (old fashioned construction of Poisson potential and extrapolation to solvent sites outside grid nodes)
        ! and of solvent(s)%vext(i,j,k,o,p)%q, the new construction with multipole expansion of the solvent molecular charge density.

        ! Since vext_q is used later in the code, if one wants the newer construction, move solvent%vext%q to vext_q
        if (getinput%log('better_poisson_solver')) then
            print*,"BETTER POISSON SOLVER [ON]"
            do concurrent (s = 1:size(solvent))
                vext_q(:,:,:,:,s) = solvent(s)%vext%q
            end do
        else ! this should be removed soon once we're used to this new construction
            print*,"BETTER POISSON SOLVER [OFF]"
            call vext_q_from_v_c (poissonSolverGrid%n, poissonSolverGrid%len, Vpoisson) ! defines Vext_q(i,j,k,o,p) that we used in the past for Poisson solver
        end if

        ! cutoff seems useless right now with the new construction
        ! where ( vext_q > 100._dp )
        !   vext_q = 100._dp
        ! else where ( vext_q < -50._dp )
        !   vext_q = -50_dp
        ! end where

    end subroutine init

    !===================================================================================================================================

    ! This SUBROUTINE uses the Gauss theorem associated to the expression of the derivation in Fourier transforms to get V(r).
    ! Field E(r)=-grad(V(r))
    ! Local expression of Gauss theorem : div(E(r))=soluteChargeDensity(r)/eps0
    ! => laplacien(V(r))=-soluteChargeDensity(r)/eps0
    ! => V(k)=soluteChargeDensity(k)/(esp0*k^2)
    ! FFT(V(k)) = V(r)

    SUBROUTINE poissonSolver (gridnode, gridlength, soluteChargeDensity, Vpoisson)

        use precision_kinds,    only: dp, i2b, i4b
        use constants,          only: fourpi , zeroC, twopi, qfact
        use system,             only: solvent
        ! Vpoisson = electrostatic potential from charge density and poisson equation

        IMPLICIT NONE

        INTEGER(i2b), INTENT(IN) :: gridnode(3)
        REAL(dp), DIMENSION(gridnode(1),gridnode(2),gridnode(3)), INTENT(IN) :: soluteChargeDensity
        REAL(dp), INTENT(IN) :: gridlength(3)
        REAL(dp), DIMENSION(gridnode(1),gridnode(2),gridnode(3)), INTENT(OUT) :: Vpoisson
        COMPLEX(dp), ALLOCATABLE, DIMENSION(:,:,:) :: soluteChargeDensity_k, Vpoisson_k
        INTEGER (i2b) :: i,j,k,m1,m2,m3,s,o,p,io
        REAL(dp) :: k2
        REAL(dp), ALLOCATABLE :: fftw3InForward(:,:,:), fftw3OutBackward(:,:,:)
        COMPLEX(dp), ALLOCATABLE :: fftw3OutForward(:,:,:), fftw3InBackward(:,:,:)
        INTEGER(i4b) :: fpspf, fpspb ! fast Poisson Solver Plan Forward or Backward
        INCLUDE "fftw3.f"


        ! IF ( ALL(abs(soluteChargeDensity)<=epsilon(1._dp)) ) THEN
        !     Vpoisson = 0._dp
        !     RETURN
        ! end if

        ALLOCATE( soluteChargeDensity_k (gridnode(1)/2+1,gridnode(2),gridnode(3)) ,SOURCE=zeroC)
        ALLOCATE( Vpoisson_k (gridnode(1)/2+1,gridnode(2),gridnode(3)) ,SOURCE=zeroC)

        CALL prepare_fftw3_for_PoissonGrid

        ! Fourier transform of the solute charge density
        fftw3InForward = soluteChargeDensity
        CALL dfftw_execute (fpspf)
        soluteChargeDensity_k = fftw3OutForward ! It is verified that at this point, FFT-1(soluteChargeDensity_k)/ (nfft1*nfft2*nfft3) = soluteChargeDensity
        ! FFT(Laplacian(V(r))) = FFT( - 4Pi charge density(r) ) in elecUnits = (ik)^2 V(k) = -4pi rho(k)
        ! V(k) = 4Pi rho(k) / k^2


        DO k = 1, gridnode(3)
            DO j = 1, gridnode(2)
                DO i = 1, gridnode(1)/2+1

                    IF ( i<=gridnode(1)/2 ) THEN
                        m1 = i-1
                    ELSE
                        m1 = i-1-gridnode(1)
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
                        Vpoisson_k(i,j,k) = soluteChargeDensity_k(i,j,k) * fourpi/k2 ! in electrostatic units : V=-4pi rho
                    ELSE
                        Vpoisson_k(i,j,k) = cmplx(0._dp,0._dp)
                    END IF

                END DO
            END DO
        END DO


        ! old construction
        if (.not. getinput%log('better_poisson_solver')) then
            ! get electrostatic potentiel in real space, that is the true Poisson potentiel. It is not solvent dependent
            fftw3InBackward = Vpoisson_k
            call dfftw_execute (fpspb)
            Vpoisson = qfact * fftw3OutBackward / REAL(PRODUCT(gridnode),dp) ! kJ/mol
        else
            ! new construction
            ! get multipolar electrostatic potential in real space. It is already solvent dependent here
            do s=1,size(solvent)
                do io = 1, grid%no
                    fftw3InBackward = Vpoisson_k * conjg( solvent(s)%sigma_k(:,:,:,io) )
                    call dfftw_execute (fpspb)
                    solvent(s)%vext(:,:,:,io)%q = qfact * fftw3OutBackward / real( product(gridnode) ,dp) ! kJ/mol
                end do
            end do
        end if

    deallocate( soluteChargeDensity_k, Vpoisson_k, fftw3InForward, fftw3OutForward, fftw3OutBackward, fftw3InBackward)

CONTAINS

    !===================================================================================================================================

    SUBROUTINE prepare_fftw3_for_poissongrid
        implicit none
        ! allocate the arrays needed as input for FFT (in_forward) or output for FFT (out_forward)
        ! or needed as input for inverse FFT (in_backward) etc.
        ALLOCATE ( fftw3InForward   ( gridnode(1)      , gridnode(2) , gridnode(3) ) )
        ALLOCATE ( fftw3OutForward  ( gridnode(1)/2 +1 , gridnode(2) , gridnode(3) ) )
        ALLOCATE ( fftw3OutBackward ( gridnode(1)      , gridnode(2) , gridnode(3) ) )
        ALLOCATE ( fftw3InBackward  ( gridnode(1)/2 +1 , gridnode(2) , gridnode(3) ) )
        ! prepare plans needed by fftw3
        CALL dfftw_plan_dft_r2c_3d &
        ( fpspf, gridnode(1), gridnode(2), gridnode(3), fftw3InForward, fftw3OutForward, FFTW_ESTIMATE )
        CALL dfftw_plan_dft_c2r_3d &
        ( fpspb, gridnode(1), gridnode(2), gridnode(3), fftw3InBackward, fftw3OutBackward, FFTW_ESTIMATE )
        ! Note that since the fast Poisson solver implies only 1 FFT in each direction, it is useless to use FFTW_MEASURE or even
        ! more rigorous planning-flags. See http://www.fftw.org/doc/Planner-Flags.html
    END SUBROUTINE prepare_fftw3_for_poissongrid

    !===================================================================================================================================

END SUBROUTINE poissonSolver

!===================================================================================================================================

SUBROUTINE vext_q_from_v_c (gridnode, gridlen, Vpoisson)

    use precision_kinds     ,ONLY: dp, i2b
    use system              ,ONLY: nb_solvent_sites, solute, solvent
    use module_grid, only: grid
    use external_potential  ,ONLY: vext_q
    use constants           ,ONLY: qfact, zero
    use mathematica         ,ONLY: TriLinearInterpolation, UTest_TrilinearInterpolation, UTest_floorNode, floorNode, ceilingNode,&
    UTest_ceilingNode, UTest_distToFloorNode, distToFloorNode, chop

    IMPLICIT NONE

    INTEGER(i2b), INTENT(IN) :: gridnode(3)
    REAL(dp), DIMENSION(gridnode(1),gridnode(2),gridnode(3)), INTENT(IN) :: Vpoisson
    REAL(dp), INTENT(IN) :: gridlen(3)
    INTEGER(i2b) :: i, j, k, o, p, m, s, nfft(3), l(3), u(3), d, io, no, ns
    real(dp), dimension(3,solvent(1)%nspec,nb_solvent_sites,grid%no) :: xmod
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

    IF(.NOT. ALLOCATED(vext_q)) STOP "vext_q should already be allocated in vext_q_from_v_c.f90"
    IF( ANY(vext_q/=0._dp) ) STOP "vext_q should be zero everywhere in vext_q_from_v_c.f90"

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


    DO CONCURRENT( s=1:solvent(1)%nspec, i=1:nfft(1), j=1:nfft(2), k=1:nfft(3), io=1:grid%no )
        vext_q_of_r_and_omega = 0.0_dp

        DO CONCURRENT (m=1:nb_solvent_sites, abs(solvent(s)%site(m)%q)>epsilon(1.0_dp))

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

            if ( ANY(l<LBOUND(Vpoisson)) .or. ANY(l>UBOUND(Vpoisson)) ) THEN
                err%pb=.TRUE.
                err%msg="Problem with l in vext_q_from_v_c.f90"
                err%l = l
                err%u = u
                err%r = r
                err%x = x
            END IF

            if ( ANY(u<LBOUND(Vpoisson)) .or. ANY(u>UBOUND(Vpoisson)) ) THEN
                err%pb=.TRUE.
                err%msg="Problem with u in vext_q_from_v_c.f90"
                err%l = l
                err%u = u
                err%r = r
                err%x = x
            end if

            cube(0,0,0) = Vpoisson (l(1),l(2),l(3))
            cube(1,0,0) = Vpoisson (u(1),l(2),l(3))
            cube(0,1,0) = Vpoisson (l(1),u(2),l(3))
            cube(0,0,1) = Vpoisson (l(1),l(2),u(3))
            cube(1,1,0) = Vpoisson (u(1),u(2),l(3))
            cube(1,0,1) = Vpoisson (u(1),l(2),u(3))
            cube(0,1,1) = Vpoisson (l(1),u(2),u(3))
            cube(1,1,1) = Vpoisson (u(1),u(2),u(3))

            vext_q_of_r_and_omega = vext_q_of_r_and_omega + solvent(s)%site(m)%q * TriLinearInterpolation(cube,r)
        end do
        Vext_q(i,j,k,io,s) = vext_q_of_r_and_omega

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

end module fastPoissonSolver

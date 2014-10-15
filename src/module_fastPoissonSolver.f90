!===================================================================================================================================
MODULE fastPoissonSolver
!===================================================================================================================================
! This module contains everything related to the fast poisson solver(s) of MDFT.

    USE precision_kinds     ,ONLY: dp,i2b
    USE system              ,ONLY: spaceGrid
    
    IMPLICIT NONE
    
    INTEGER(i2b), PRIVATE :: i
    CHARACTER(180), PRIVATE :: j
    REAL(dp), ALLOCATABLE, DIMENSION(:,:,:), PRIVATE :: soluteChargeDensity, Vpoisson ! TODO THIS ON FINE GRID
    TYPE :: PoissonGrid
        INTEGER(i2b) :: nnod(3)
        REAL(dp) :: len(3)
    END TYPE
    TYPE(PoissonGrid), PRIVATE :: pgrid
    PRIVATE
    PUBLIC :: init

    CONTAINS
    
!===================================================================================================================================

    SUBROUTINE init
        IMPLICIT NONE
        pgrid%nnod = spaceGrid%n_nodes
        pgrid%len = spaceGrid%length
        ALLOCATE ( soluteChargeDensity (pgrid%nnod(1),pgrid%nnod(2),pgrid%nnod(3)) ,SOURCE=0._dp, STAT=i, ERRMSG=j)
            IF (i/=0) THEN; PRINT*,j; STOP "This problem arises in module fastPoissonSolver"; END IF
        ALLOCATE ( Vpoisson            (pgrid%nnod(1),pgrid%nnod(2),pgrid%nnod(3)) ,SOURCE=0._dp, STAT=i, ERRMSG=j)
            IF (i/=0) THEN; PRINT*,j; STOP "This problem arises in module fastPoissonSolver"; END IF
        CALL soluteChargeDensityFromSoluteChargeCoordinates (pgrid%nnod, pgrid%len, soluteChargeDensity)
        CALL poissonSolver (pgrid%nnod, pgrid%len, soluteChargeDensity, Vpoisson)
        CALL vext_q_from_v_c (pgrid%nnod, pgrid%len, Vpoisson)
    END SUBROUTINE init

!===========================================================================================================================

! This SUBROUTINE uses the Gauss theorem associated to the expression of the derivation in Fourier transforms to get V(r).
! Field E(r)=-grad(V(r))
! Local expression of Gauss theorem : div(E(r))=soluteChargeDensity(r)/eps0
! => laplacien(V(r))=-soluteChargeDensity(r)/eps0
! => V(k)=soluteChargeDensity(k)/(esp0*k^2)
! FFT(V(k)) = V(r)

    SUBROUTINE poissonSolver (gridnode, gridlen, soluteChargeDensity, Vpoisson)

        USE precision_kinds,    ONLY: dp, i2b, i4b
        USE constants,          ONLY: fourpi , zeroC, twopi
        ! Vpoisson = electrostatic potential from charge density and poisson equation

        IMPLICIT NONE
        
        INTEGER(i2b), INTENT(IN) :: gridnode(3)
        REAL(dp), DIMENSION(gridnode(1),gridnode(2),gridnode(3)), INTENT(IN) :: soluteChargeDensity
        REAL(dp), INTENT(IN) :: gridlen(3)
        REAL(dp), DIMENSION(gridnode(1),gridnode(2),gridnode(3)), INTENT(OUT) :: Vpoisson
        COMPLEX(dp), ALLOCATABLE, DIMENSION(:,:,:) :: soluteChargeDensity_k,Vpoisson_k
        INTEGER (i2b) :: i,j,k,m1,m2,m3
        REAL(dp) :: k2
        REAL(dp), ALLOCATABLE :: fftw3InForward(:,:,:), fftw3OutBackward(:,:,:)
        COMPLEX(dp), ALLOCATABLE :: fftw3OutForward(:,:,:), fftw3InBackward(:,:,:)
        INTEGER(i4b) :: fpspf, fpspb ! fast Poisson Solver Plan Forward or Backward
        INCLUDE "fftw3.f"


        IF ( ALL(abs(soluteChargeDensity)>epsilon(1._dp)) ) THEN
            Vpoisson = 0._dp
            RETURN
        end if

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
                    
                    k2 =  (twopi/gridlen(1)*REAL(m1,dp))**2 &
                        + (twopi/gridlen(2)*REAL(m2,dp))**2 &
                        + (twopi/gridlen(3)*REAL(m3,dp))**2

                    IF ( abs(k2) > epsilon(1._dp) ) THEN
                        Vpoisson_k(i,j,k) = soluteChargeDensity_k(i,j,k) * fourpi/k2 ! in electrostatic units : V=-4pi rho
                    ELSE
                        Vpoisson_k(i,j,k) = CMPLX(0._dp,0._dp)
                    END IF
                
                END DO
            END DO
        END DO

        ! get real space potential V(r)
        fftw3InBackward = Vpoisson_k
        CALL dfftw_execute (fpspb)
        Vpoisson = fftw3OutBackward / REAL(PRODUCT(gridnode),dp)
    
        DEALLOCATE (soluteChargeDensity_k, Vpoisson_k, fftw3InForward, fftw3OutForward, fftw3OutBackward, fftw3InBackward)

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
            ! Note that since the fast Poisson solver implies only 1 FFT in each direct, it is useless to use FFTW_MEASURE or even 
            ! more rigorous planning-flags. See http://www.fftw.org/doc/Planner-Flags.html
        END SUBROUTINE prepare_fftw3_for_poissongrid

!===================================================================================================================================

    END SUBROUTINE poissonSolver

!===================================================================================================================================

    SUBROUTINE vext_q_from_v_c (gridnode, gridlen, Vpoisson)

    USE precision_kinds     ,ONLY: dp, i2b
    USE system              ,ONLY: nb_solvent_sites, nb_species, spaceGrid, solute, solvent
    USE quadrature          ,ONLY: angGrid, molRotGrid, Rotxx,Rotxy,Rotxz,Rotyx,Rotyy,Rotyz,Rotzx,Rotzy,Rotzz
    USE external_potential  ,ONLY: vext_q
    USE constants           ,ONLY: qfact, zero
    USE mathematica         ,ONLY: TriLinearInterpolation, UTest_TrilinearInterpolation, UTest_floorNode, floorNode, ceilingNode,&
                                   UTest_ceilingNode, UTest_distToFloorNode, distToFloorNode, chop

    IMPLICIT NONE

    INTEGER(i2b), INTENT(IN) :: gridnode(3)
    REAL(dp), DIMENSION(gridnode(1),gridnode(2),gridnode(3)), INTENT(IN) :: Vpoisson
    REAL(dp), INTENT(IN) :: gridlen(3)
    INTEGER(i2b) :: i, j, k, o, p, m, s, nfft(3), l(3), u(3)
    REAL(dp), DIMENSION(nb_solvent_sites,molRotGrid%n_angles,angGrid%n_angles) :: xmod, ymod, zmod
    REAL(dp) :: vext_q_of_r_and_omega ! external potential for a given i,j,k,omega & psi.
    REAL(dp) :: r(3), cube(0:1,0:1,0:1), dl(3), x(3)
    TYPE :: testtype
        LOGICAL :: pb
        CHARACTER(180) :: msg
        real(dp) :: l(3),u(3),r(3),x(3)
    END TYPE
    TYPE(testtype) :: err

    nfft = spaceGrid%n_nodes
    dl = spaceGrid%dl

    IF(.NOT. ALLOCATED(vext_q)) STOP "vext_q should already be allocated in vext_q_from_v_c.f90"
    IF( ANY(vext_q/=0._dp) ) STOP "vext_q should be zero everywhere in vext_q_from_v_c.f90"

    IF ( ALL(solute%site%q == zero) .OR. ALL(solvent(1)%site%q == zero) ) RETURN

    ! Tabulate the cartesian coordinates of all solvent sites, for all molecular orientations, centered on any MDFT's grid node.
    DO CONCURRENT( o=1:angGrid%n_angles , p=1:molRotGrid%n_angles , m=1:SIZE(solvent(1)%site) )
        xmod(m,p,o) = DOT_PRODUCT( [Rotxx(o,p),Rotxy(o,p),Rotxz(o,p)] , solvent(1)%site(m)%r )
        ymod(m,p,o) = DOT_PRODUCT( [Rotyx(o,p),Rotyy(o,p),Rotyz(o,p)] , solvent(1)%site(m)%r )
        zmod(m,p,o) = DOT_PRODUCT( [Rotzx(o,p),Rotzy(o,p),Rotzz(o,p)] , solvent(1)%site(m)%r )
    END DO

    CALL UTest_floorNode
    CALL UTest_ceilingNode
    CALL UTest_TrilinearInterpolation
    CALL UTest_distToFloorNode
    ! Compute external potential for each combination of solvent side and grid node and orientation
    err%pb=.FALSE. ! becomes TRUE if a problem is detected during execution.
    
    
    DO CONCURRENT( s=1:nb_species, i=1:nfft(1), j=1:nfft(2), k=1:nfft(3), o=1:angGrid%n_angles, p=1:molRotGrid%n_angles )
        vext_q_of_r_and_omega = 0.0_dp
        
        DO CONCURRENT (m=1:nb_solvent_sites, abs(solvent(1)%site(m)%q)>epsilon(1.0_dp))
        
            x = (REAL([i,j,k],dp)-1.0_dp)*dl + [xmod(m,p,o),ymod(m,p,o),zmod(m,p,o)]! cartesian coordinate x of the solvent site m. May be outside the supercell.
            x(1)=chop(x(1))
            x(2)=chop(x(2))
            x(3)=chop(x(3))

            l = floorNode(gridnode,gridlen,x,.TRUE.) ! r should be in full cartesian coordinates between -infty and +infty
            u = ceilingNode(gridnode,gridlen,x,.TRUE.)
            r = distToFloorNode(gridnode,gridlen,x,.TRUE.) ! 0 <= distToFloorNode < 1
            
            IF( ANY(r<0._dp) .or. ANY(r>1._dp) ) THEN
                err%pb=.TRUE.
                err%msg="Problem with r in vext_q_from_v_c.f90"
                err%l = l
                err%u = u
                err%r = r
                err%x = x
            END IF
            
            IF( ANY(l<LBOUND(Vpoisson)) .or. ANY(l>UBOUND(Vpoisson)) ) THEN
                err%pb=.TRUE.
                err%msg="Problem with l in vext_q_from_v_c.f90"
                err%l = l
                err%u = u
                err%r = r
                err%x = x
            END IF
            
            IF( ANY(u<LBOUND(Vpoisson)) .or. ANY(u>UBOUND(Vpoisson)) ) THEN
                err%pb=.TRUE.
                err%msg="Problem with u in vext_q_from_v_c.f90"
                err%l = l
                err%u = u
                err%r = r
                err%x = x
            END IF
            
            cube(0,0,0) = Vpoisson (l(1),l(2),l(3))
            cube(1,0,0) = Vpoisson (u(1),l(2),l(3))
            cube(0,1,0) = Vpoisson (l(1),u(2),l(3))
            cube(0,0,1) = Vpoisson (l(1),l(2),u(3))
            cube(1,1,0) = Vpoisson (u(1),u(2),l(3))
            cube(1,0,1) = Vpoisson (u(1),l(2),u(3))
            cube(0,1,1) = Vpoisson (l(1),u(2),u(3))
            cube(1,1,1) = Vpoisson (u(1),u(2),u(3))
            
            vext_q_of_r_and_omega = vext_q_of_r_and_omega + solvent(1)%site(m)%q * TriLinearInterpolation(cube,r)
        END DO
        
        vext_q_of_r_and_omega = qfact * vext_q_of_r_and_omega
        CALL limit_potential(vext_q_of_r_and_omega)
        Vext_q(i,j,k,o,p,s) = vext_q_of_r_and_omega
        
    END DO
    
    IF (err%pb) THEN
        PRINT*,"msg:",err%msg
        print*,"Lxyz:",gridlen
        PRINT*,"x:",err%x
        PRINT*,"l/u min:",1,1,1
        print*,"l/u max:",gridnode
        PRINT*,"l:",err%l
        PRINT*,"u:",err%u
        PRINT*,"r:",err%r
        STOP
    END IF

    CONTAINS

!===================================================================================================================================

        PURE SUBROUTINE limit_potential(v)
        ! Limit the potential to -100kJ/mol<v<100kJ/mol
            IMPLICIT NONE
            REAL(dp), INTENT(INOUT) :: v
            REAL(dp), PARAMETER :: vu=100._dp, vl=-100._dp ! kJ/mol
            IF (v>vu) THEN
                v=vu
            ELSE IF (v<vl) THEN
                v=vl
            END IF
        END SUBROUTINE limit_potential

!===================================================================================================================================
        
    END SUBROUTINE vext_q_from_v_c

!===================================================================================================================================
        
END MODULE fastPoissonSolver
!===================================================================================================================================

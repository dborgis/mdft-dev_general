SUBROUTINE vext_q_from_v_c (gridnode, Vpoisson)

    USE precision_kinds     ,ONLY: dp, i2b
    USE system              ,ONLY: chg_mol , chg_solv , nb_solvent_sites , nb_species , &
                                    beta , id_solv , beta , spaceGrid, soluteSite,& 
                                    solventSite
    USE quadrature          ,ONLY: angGrid, molRotGrid,Rotxx,Rotxy,Rotxz,Rotyx,Rotyy,Rotyz,Rotzx,Rotzy,Rotzz
    USE external_potential  ,ONLY: vext_q
    USE constants, ONLY: fourpi, qfact, qunit, zero
    USE input, ONLY: verbose
    USE mathematica ,ONLY: TriLinearInterpolation, UTest_TrilinearInterpolation

    IMPLICIT NONE

    INTEGER(i2b), INTENT(IN) :: gridnode(3)
    REAL(dp), DIMENSION(gridnode(1),gridnode(2),gridnode(3)), INTENT(IN) :: Vpoisson
    INTEGER(i2b) :: i, j, k, o, p, m, s, nfft(3), l(3), u(3)
    REAL(dp), DIMENSION ( nb_solvent_sites , molRotGrid%n_angles , angGrid%n_angles ) :: xmod , ymod , zmod
    REAL(dp):: xq,yq,zq ! solvent coordinates in indices referential (ie real between 0 and nfft1+1)
    INTEGER(i2b):: im , jm , km , ip , jp , kp ! 6 indices from which all corners of cube surrounding point charge are defined
    REAL(dp):: wim , wjm , wkm , wip , wjp , wkp ! weights associated to each 'corner' index
    REAL(dp):: vpsi ! external potential for a given i,j,k,omega & psi.
    REAL(dp):: average_over_psi ! boltzmann average of vpsi over psi for a given i,j,k,omega
    REAL(dp):: charge
    REAL(dp) :: r(3), cube(0:1,0:1,0:1), dl(3)
    TYPE :: testtype
        LOGICAL :: pb
        CHARACTER(180) :: msg
    END TYPE
    TYPE(testtype) :: err

    nfft=spaceGrid%n_nodes
    dl = spaceGrid%dl

    IF(.NOT. ALLOCATED(vext_q)) STOP "vext_q should already be allocated in vext_q_from_v_c.f90"
    IF( ANY(vext_q/=0._dp) ) STOP "vext_q should be zero everywhere in vext_q_from_v_c.f90"

    IF ( ALL(soluteSite%q == zero) .OR. ALL(solventSite%q == zero) ) RETURN

    ! Tabulate the cartesian coordinates of all solvent sites, for all molecular orientations, centered on any MDFT's grid node.
    DO CONCURRENT( o=1:angGrid%n_angles , p=1:molRotGrid%n_angles , m=1:SIZE(solventSite) )
        xmod(m,p,o) = Rotxx(o,p) * solventSite(m)%r(1) + Rotxy(o,p) * solventSite(m)%r(2) + Rotxz(o,p) * solventSite(m)%r(3)
        ymod(m,p,o) = Rotyx(o,p) * solventSite(m)%r(1) + Rotyy(o,p) * solventSite(m)%r(2) + Rotyz(o,p) * solventSite(m)%r(3)
        zmod(m,p,o) = Rotzx(o,p) * solventSite(m)%r(1) + Rotzy(o,p) * solventSite(m)%r(2) + Rotzz(o,p) * solventSite(m)%r(3)
    END DO

    CALL UTest_TrilinearInterpolation
    ! Compute external potential for each combination of solvent side and grid node and orientation
    err%pb=.FALSE.
    DO CONCURRENT( s=1:nb_species, i=1:nfft(1), j=1:nfft(2), k=1:nfft(3), o=1:angGrid%n_angles, p=1:molRotGrid%n_angles )
        vpsi = 0.0_dp
        DO CONCURRENT (m=1:nb_solvent_sites, solventSite(m)%q/=0._dp)
            r = (REAL([i,j,k],dp)-1.0_dp)*dl + [xmod(m,p,o),ymod(m,p,o),zmod(m,p,o)]! cartesian coordinate x of the solvent site m. May be outside the supercell.
            r = r/dl ! cartesian coordinates in index scale
            l = FLOOR(r) ! index of node just below
            u = l+1
            r = r-REAL(l,dp)  ! cartesian coordinates between 0 and 1
            l = MODULO(l,nfft)+1 ! Note for later: here we implicitely consider periodic boundary conditions.
            u = MODULO(u,nfft)+1
            if( ANY(r<0._dp) .or. ANY(r>1._dp) ) &
                THEN; err%pb=.true.; err%msg="Problem with r in vext_q_from_v_c.f90"; END IF
            if( ANY(l<LBOUND(Vpoisson)) .or. ANY(l>UBOUND(Vpoisson)) )&
                THEN; err%pb=.true.; err%msg="Problem with l in vext_q_from_v_c.f90"; END IF
            if( ANY(u<LBOUND(Vpoisson)) .or. ANY(u>UBOUND(Vpoisson)) )&
                THEN; err%pb=.true.; err%msg="Problem with u in vext_q_from_v_c.f90"; END IF
            cube(0,0,0) = Vpoisson (l(1),l(2),l(3))
            cube(1,0,0) = Vpoisson (u(1),l(2),l(3))
            cube(0,1,0) = Vpoisson (l(1),u(2),l(3))
            cube(0,0,1) = Vpoisson (l(1),l(2),u(3))
            cube(1,1,0) = Vpoisson (u(1),u(2),l(3))
            cube(1,0,1) = Vpoisson (u(1),l(2),u(3))
            cube(0,1,1) = Vpoisson (l(1),u(2),u(3))
            cube(1,1,1) = Vpoisson (u(1),u(2),u(3))
            vpsi = vpsi + solventSite(m)%q * TriLinearInterpolation(cube,r)
        END DO
        vpsi = qfact * vpsi
        CALL limit_potential(vpsi)
        Vext_q(i,j,k,o,p,s) = vpsi
    END DO
    IF (err%pb) THEN
        PRINT*,err%msg
        STOP
    END IF

    CONTAINS

        !===========================================================================================================================
        PURE SUBROUTINE limit_potential(v)
        !===========================================================================================================================
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
        !===========================================================================================================================
        
END SUBROUTINE vext_q_from_v_c

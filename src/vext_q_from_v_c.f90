SUBROUTINE vext_q_from_v_c (gridnode, Vpoisson)

    USE precision_kinds     ,ONLY: dp, i2b
    USE system              ,ONLY: chg_mol , chg_solv , nb_solvent_sites , nb_species , &
                                    beta , id_solv , beta , spaceGrid, soluteSite,& 
                                    solventSite
    USE quadrature          ,ONLY: angGrid, molRotGrid,Rotxx,Rotxy,Rotxz,Rotyx,Rotyy,Rotyz,Rotzx,Rotzy,Rotzz
    USE external_potential  ,ONLY: vext_q
    USE constants, ONLY: fourpi, qfact, qunit, zero
    USE input, ONLY: verbose

    IMPLICIT NONE

    INTEGER(i2b), INTENT(IN) :: gridnode(3)
    REAL(dp), DIMENSION(gridnode(1),gridnode(2),gridnode(3)), INTENT(IN) :: Vpoisson
    INTEGER(i2b) :: i, j, k, o, p, m, z,s,nfft1 , nfft2 , nfft3
    REAL(dp), DIMENSION ( nb_solvent_sites , molRotGrid%n_angles , angGrid%n_angles ) :: xmod , ymod , zmod
    REAL(dp):: xq,yq,zq ! solvent coordinates in indices referential (ie real between 0 and nfft1+1)
    INTEGER(i2b):: im , jm , km , ip , jp , kp ! 6 indices from which all corners of cube surrounding point charge are defined
    REAL(dp):: wim , wjm , wkm , wip , wjp , wkp ! weights associated to each 'corner' index
    REAL(dp):: vpsi ! external potential for a given i,j,k,omega & psi.
    REAL(dp):: average_over_psi ! boltzmann average of vpsi over psi for a given i,j,k,omega
    REAL(dp):: charge

    nfft1=spaceGrid%n_nodes(1)
    nfft2=spaceGrid%n_nodes(2)
    nfft3=spaceGrid%n_nodes(3)
    z=0
    IF(.NOT. ALLOCATED(vext_q)) STOP "vext_q should already be allocated in vext_q_from_v_c.f90"
    IF( ANY(vext_q/=0._dp) ) STOP "vext_q should be zero everywhere in vext_q_from_v_c.f90"

    IF ( ALL(soluteSite%q == zero) .OR. ALL(solventSite%q == zero) ) RETURN

    ! Tabulate the cartesian coordinates of all solvent sites, for all molecular orientations, centered on any MDFT's grid node.
    DO o = 1 , angGrid%n_angles
        DO p = 1 , molRotGrid%n_angles
            DO m = 1 , SIZE(solventSite)
                xmod(m,p,o) = Rotxx(o,p) * solventSite(m)%r(1) + Rotxy(o,p) * solventSite(m)%r(2) + Rotxz(o,p) * solventSite(m)%r(3)
                ymod(m,p,o) = Rotyx(o,p) * solventSite(m)%r(1) + Rotyy(o,p) * solventSite(m)%r(2) + Rotyz(o,p) * solventSite(m)%r(3)
                zmod(m,p,o) = Rotzx(o,p) * solventSite(m)%r(1) + Rotzy(o,p) * solventSite(m)%r(2) + Rotzz(o,p) * solventSite(m)%r(3)
            END DO
        END DO
    END DO
    
    ! Compute external potential due to charge density
    DO CONCURRENT( s=1:nb_species, i=1:nfft1, j=1:nfft2, k=1:nfft3, o=1:angGrid%n_angles, p=1:molRotGrid%n_angles )
        vpsi = 0.0_dp
        DO CONCURRENT (m=1:nb_solvent_sites, solventSite(m)%q/=0._dp)
            xq = MODULO( REAL(i-1,dp) * spaceGrid%dl(1) + xmod(m,p,o) ,spaceGrid%length(1)) /spaceGrid%dl(1)
            yq = MODULO( REAL(j-1,dp) * spaceGrid%dl(2) + ymod(m,p,o) ,spaceGrid%length(2)) /spaceGrid%dl(2)
            zq = MODULO( REAL(k-1,dp) * spaceGrid%dl(3) + zmod(m,p,o) ,spaceGrid%length(3)) /spaceGrid%dl(3)
            im = INT(xq) +1
            jm = INT(yq) +1
            km = INT(zq) +1
            IF ( im == nfft1 + 1 ) im = 1
            IF ( jm == nfft2 + 1 ) jm = 1
            IF ( km == nfft3 + 1 ) km = 1
            ip = im + 1
            jp = jm + 1
            kp = km + 1
            IF ( ip == nfft1 + 1 ) ip = 1
            IF ( jp == nfft2 + 1 ) jp = 1
            IF ( kp == nfft3 + 1 ) kp = 1
            wip = xq - REAL(INT( xq ,i2b),dp)
            wjp = yq - REAL(INT( yq ,i2b),dp)
            wkp = zq - REAL(INT( zq ,i2b),dp)
            wim = 1.0_dp - wip
            wjm = 1.0_dp - wjp
            wkm = 1.0_dp - wkp
            vpsi = vpsi + solventSite(m)%q * (  Vpoisson (im,jm,km) * wim * wjm * wkm &
                                              + Vpoisson (ip,jm,km) * wip * wjm * wkm &
                                              + Vpoisson (im,jp,km) * wim * wjp * wkm &
                                              + Vpoisson (im,jm,kp) * wim * wjm * wkp &
                                              + Vpoisson (ip,jp,km) * wip * wjp * wkm &
                                              + Vpoisson (ip,jm,kp) * wip * wjm * wkp &
                                              + Vpoisson (im,jp,kp) * wim * wjp * wkp &
                                              + Vpoisson (ip,jp,kp) * wip * wjp * wkp )
        END DO
        vpsi = qfact * vpsi
        IF (vpsi>100._dp) THEN
            Vext_q(i,j,k,o,p,s)=100._dp
        ELSE IF (vpsi<-100_dp) THEN
            Vext_q(i,j,k,o,p,s)=-100._dp
        ELSE
            Vext_q(i,j,k,o,p,s)=vpsi
        END IF
    END DO

END SUBROUTINE vext_q_from_v_c

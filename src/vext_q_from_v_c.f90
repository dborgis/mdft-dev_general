SUBROUTINE vext_q_from_v_c (V_c)

    USE precision_kinds, ONLY: dp, i2b
    USE system, ONLY: chg_mol , chg_solv , x_solv , y_solv , z_solv , nb_solvent_sites , nb_species , Lx , Ly , &
                        Lz , deltax , deltay , deltaz , beta , id_solv , beta , nfft1 , nfft2 , nfft3, spaceGrid, soluteSite, &
                        solventSite
    USE quadrature, ONLY: angGrid, molRotGrid,Rotxx,Rotxy,Rotxz,Rotyx,Rotyy,Rotyz,Rotzx,Rotzy,Rotzz
    USE external_potential, ONLY: vext_q , vext_lj
    ! v_c = electrostatic potential from charge density and poisson equation
    ! vext_q = electrostatic potential energy in general and as used in the calculation of the total external potential
    USE constants, ONLY: fourpi, qfact, qunit, zero
    ! qfact = qunit ** 2 * 1.0e-3_dp * Navo / ( fourpi * eps0 * 1.0e-10_dp ) ! electrostatic potential unit so that QFACT*q*q/r is kJ/mol
    USE input, ONLY: verbose

    IMPLICIT NONE

    REAL(dp), DIMENSION (spaceGrid%n_nodes(1),spaceGrid%n_nodes(2),spaceGrid%n_nodes(3)), INTENT(IN) :: V_c
    INTEGER(i2b) :: i, j, k, o, p, m, z,s
    REAL(dp), DIMENSION ( nb_solvent_sites , molRotGrid%n_angles , angGrid%n_angles ) :: xmod , ymod , zmod
    REAL(dp):: xq,yq,zq ! solvent coordinates in indices referential (ie real between 0 and nfft1+1)
    INTEGER(i2b):: im , jm , km , ip , jp , kp ! 6 indices from which all corners of cube surrounding point charge are defined
    REAL(dp):: wim , wjm , wkm , wip , wjp , wkp ! weights associated to each 'corner' index
    REAL(dp):: vpsi ! external potential for a given i,j,k,omega & psi.
    REAL(dp):: average_over_psi ! boltzmann average of vpsi over psi for a given i,j,k,omega
    REAL(dp):: charge

    z=0
    ! be sure this subroutine comes at right moment
    ! init (again during debug) vext_q
    IF(.NOT. ALLOCATED(vext_q)) ALLOCATE( vext_q (spaceGrid%n_nodes(1), spaceGrid%n_nodes(2), spaceGrid%n_nodes(3),&
        angGrid%n_angles, molRotGrid%n_angles, nb_species ), SOURCE=0._dp )
    
    IF ( ALL(soluteSite%q == zero) .OR. ALL(solventSite%q == zero) ) RETURN

    ! Tabulate rotation matrix * solvent coordinates
    do o = 1 , angGrid%n_angles
        do p = 1 , molRotGrid%n_angles
            do m = 1 , nb_solvent_sites
                xmod ( m, p , o ) = Rotxx ( o , p ) * x_solv ( m ) + Rotxy ( o , p ) * y_solv ( m ) + Rotxz ( o , p ) * z_solv (m)
                ymod ( m, p , o ) = Rotyx ( o , p ) * x_solv ( m ) + Rotyy ( o , p ) * y_solv ( m ) + Rotyz ( o , p ) * z_solv (m)
                zmod ( m, p , o ) = Rotzx ( o , p ) * x_solv ( m ) + Rotzy ( o , p ) * y_solv ( m ) + Rotzz ( o , p ) * z_solv (m)
            end do
        end do
    end do
    
    ! Compute external potential due to charge density
    DO CONCURRENT( s=1:nb_species, i=1:nfft1, j=1:nfft2, k=1:nfft3, o=1:angGrid%n_angles, p=1:molRotGrid%n_angles )
        vpsi = 0.0_dp
        DO m = 1 , nb_solvent_sites
            charge = chg_solv(id_solv(m))
            IF ( charge==0.0_dp ) CYCLE
            xq = MODULO( REAL(i-1,dp) * deltax + xmod(m,p,o) ,Lx)/deltax
            yq = MODULO( REAL(j-1,dp) * deltay + ymod(m,p,o) ,Ly)/deltay
            zq = MODULO( REAL(k-1,dp) * deltaz + zmod(m,p,o) ,Lz)/deltaz
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
            wim = (    1.0_dp - (   xq - real( int( xq ,i2b),dp)   )    )
            wjm = (    1.0_dp - (   yq - real( int( yq ,i2b),dp)   )    )
            wkm = (    1.0_dp - (   zq - real( int( zq ,i2b),dp)   )    )
            wip = (             (   xq - real( int( xq ,i2b),dp)   )    )
            wjp = (             (   yq - real( int( yq ,i2b),dp)   )    )
            wkp = (             (   zq - real( int( zq ,i2b),dp)   )    )
            vpsi = vpsi + (chg_solv ( id_solv (m) ) * qfact * ( V_c (im,jm,km) * wim * wjm * wkm &
                                                              + V_c (ip,jm,km) * wip * wjm * wkm &
                                                              + V_c (im,jp,km) * wim * wjp * wkm &
                                                              + V_c (im,jm,kp) * wim * wjm * wkp &
                                                              + V_c (ip,jp,km) * wip * wjp * wkm &
                                                              + V_c (ip,jm,kp) * wip * wjm * wkp &
                                                              + V_c (im,jp,kp) * wim * wjp * wkp &
                                                              + V_c (ip,jp,kp) * wip * wjp * wkp ))
        END DO
        IF (vpsi>100._dp) THEN
            Vext_q(i,j,k,o,p,s)=100._dp
        ELSE IF (vpsi<-100_dp) THEN
            Vext_q(i,j,k,o,p,s)=-100._dp
        ELSE
            Vext_q(i,j,k,o,p,s)=vpsi
        END IF
    END DO

END SUBROUTINE vext_q_from_v_c

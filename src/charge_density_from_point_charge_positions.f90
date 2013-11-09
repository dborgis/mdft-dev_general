! This subroutine reads the positions of the charges of the solute sites. It extrapolates them to the grid points in order to 
! get a charge density rho_c ( i , j , k )

SUBROUTINE charge_density_from_point_charge_positions
    
    USE precision_kinds, ONLY: i2b , dp
    USE system, ONLY: rho_c, spaceGrid, soluteSite
    
    IMPLICIT NONE
    
    INTEGER(i2b):: s
    REAL(dp) :: xq , yq, zq ! coordinates of the charge in indicial coordinates
    INTEGER(i2b) :: im , jm , km , ip , jp , kp ! indices of corner in indicial coordinates
    REAL(dp) :: wim , wjm , wkm , wip , wjp , wkp ! weight associated to each index
    CHARACTER(50) :: filename

    ALLOCATE( rho_c ( spaceGrid%n_nodes(1), spaceGrid%n_nodes(2), spaceGrid%n_nodes(3) ), SOURCE=0._dp )

    ! extrapolate each solute point charge to grid nodes
    DO s = 1 , SIZE(soluteSite)

        IF ( soluteSite(s)%q == 0.0_dp ) CYCLE ! if the solute does not have charge, go to next solute
        
        ! transform cartesian coordinates 0 <= x < Lx in 'indicial coordinates' 0 <= xq < nfft1
        xq = modulo ( soluteSite(s)%r(1) , spaceGrid%length(1) ) / spaceGrid%dl(1)
        yq = modulo ( soluteSite(s)%r(2) , spaceGrid%length(2) ) / spaceGrid%dl(2)
        zq = modulo ( soluteSite(s)%r(3) , spaceGrid%length(3) ) / spaceGrid%dl(3)
        ! get coordinates of grid node juste below (corner of the cube with smallest indices)
        ! +1 is because indexation does not begin to 0 but to 1. Thus, when coordinate xq=0, it corresponds to index 1.
        im = int ( xq ) + 1
        jm = int ( yq ) + 1
        km = int ( zq ) + 1
        ! grid node juste above (corner of the cube with highest indices)
        ip = im + 1
        jp = jm + 1
        kp = km + 1
        ! this corner should never have coordinates higher than nfft
        IF ( ip == spaceGrid%n_nodes(1) + 1 ) ip = 1
        IF ( jp == spaceGrid%n_nodes(2) + 1 ) jp = 1
        IF ( kp == spaceGrid%n_nodes(3) + 1 ) kp = 1
        ! define weights associated with each corner
        wim = (    1.0_dp - (   xq - real(int(xq,i2b),dp)   )    )
        wjm = (    1.0_dp - (   yq - real(int(yq,i2b),dp)   )    )
        wkm = (    1.0_dp - (   zq - real(int(zq,i2b),dp)   )    )
        wip = (             (   xq - real(int(xq,i2b),dp)   )    )
        wjp = (             (   yq - real(int(yq,i2b),dp)   )    )
        wkp = (             (   zq - real(int(zq,i2b),dp)   )    )
        ! increase density accordingly
        rho_c ( im , jm , km ) = rho_c ( im , jm , km ) + soluteSite(s)%q * wim * wjm * wkm
        rho_c ( ip , jm , km ) = rho_c ( ip , jm , km ) + soluteSite(s)%q * wip * wjm * wkm
        rho_c ( im , jp , km ) = rho_c ( im , jp , km ) + soluteSite(s)%q * wim * wjp * wkm
        rho_c ( im , jm , kp ) = rho_c ( im , jm , kp ) + soluteSite(s)%q * wim * wjm * wkp
        rho_c ( ip , jp , km ) = rho_c ( ip , jp , km ) + soluteSite(s)%q * wip * wjp * wkm
        rho_c ( ip , jm , kp ) = rho_c ( ip , jm , kp ) + soluteSite(s)%q * wip * wjm * wkp
        rho_c ( im , jp , kp ) = rho_c ( im , jp , kp ) + soluteSite(s)%q * wim * wjp * wkp
        rho_c ( ip , jp , kp ) = rho_c ( ip , jp , kp ) + soluteSite(s)%q * wip * wjp * wkp
    END DO
    rho_c = rho_c / spaceGrid%dv ! charge density is in charge per unit volume

    CALL write_to_cube_file ( rho_c , 'output/charge_density.cube' )
    
end subroutine charge_density_from_point_charge_positions

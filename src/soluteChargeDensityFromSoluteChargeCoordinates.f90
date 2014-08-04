! This subroutine uses the positions of the charges of the solute to extrapolate them on the grid nodes in order to 
! get the solute charge density soluteChargeDensity (i,j,k).

SUBROUTINE soluteChargeDensityFromSoluteChargeCoordinates (gridnode, gridlen, soluteChargeDensity)
    
    USE precision_kinds, ONLY: i2b,dp
    USE system,          ONLY: spaceGrid,soluteSite,nb_species
    USE input,           ONLY: verbose
    USE mathematica     ,ONLY: distToFloorNode, floorNode, ceilingNode
    
    IMPLICIT NONE

    INTEGER(i2b), INTENT(IN) :: gridnode(3)
    REAL(dp), INTENT(IN) :: gridlen(3)
    REAL(dp), DIMENSION(gridnode(1),gridnode(2),gridnode(3)), INTENT(OUT) :: soluteChargeDensity
    INTEGER(i2b):: s, m(3), p(3)
    REAL(dp) :: r(3) ! coordinates of the charge in indicial coordinates
    REAL(dp) :: wm(3), wp(3) ! weight associated to each index
    CHARACTER(50) :: filename

    soluteChargeDensity = 0._dp

    ! extrapolate each solute point charge to grid nodes
    DO s = 1 , SIZE(soluteSite)

        IF ( soluteSite(s)%q == 0.0_dp ) CYCLE ! if the solute does not have charge, go to next solute

        r = distToFloorNode (spaceGrid%n_nodes,spaceGrid%length,soluteSite(s)%r,.TRUE.)
        m = floorNode       (spaceGrid%n_nodes,spaceGrid%length,soluteSite(s)%r,.TRUE.)
        p = ceilingNode     (spaceGrid%n_nodes,spaceGrid%length,soluteSite(s)%r,.TRUE.)

        wp = r ! weights
        wm = 1._dp - r

        ! increase density accordingly
        soluteChargeDensity (m(1),m(2),m(3)) = soluteChargeDensity (m(1),m(2),m(3)) + soluteSite(s)%q * wm(1) * wm(2) * wm(3)
        soluteChargeDensity (p(1),m(2),m(3)) = soluteChargeDensity (p(1),m(2),m(3)) + soluteSite(s)%q * wp(1) * wm(2) * wm(3)
        soluteChargeDensity (m(1),p(2),m(3)) = soluteChargeDensity (m(1),p(2),m(3)) + soluteSite(s)%q * wm(1) * wp(2) * wm(3)
        soluteChargeDensity (m(1),m(2),p(3)) = soluteChargeDensity (m(1),m(2),p(3)) + soluteSite(s)%q * wm(1) * wm(2) * wp(3)
        soluteChargeDensity (p(1),p(2),m(3)) = soluteChargeDensity (p(1),p(2),m(3)) + soluteSite(s)%q * wp(1) * wp(2) * wm(3)
        soluteChargeDensity (p(1),m(2),p(3)) = soluteChargeDensity (p(1),m(2),p(3)) + soluteSite(s)%q * wp(1) * wm(2) * wp(3)
        soluteChargeDensity (m(1),p(2),p(3)) = soluteChargeDensity (m(1),p(2),p(3)) + soluteSite(s)%q * wm(1) * wp(2) * wp(3)
        soluteChargeDensity (p(1),p(2),p(3)) = soluteChargeDensity (p(1),p(2),p(3)) + soluteSite(s)%q * wp(1) * wp(2) * wp(3)
    END DO
    soluteChargeDensity = soluteChargeDensity / spaceGrid%dv ! charge density is in charge per unit volume

    IF (verbose) THEN
        filename='output/soluteChargeDensity.cube'
        CALL write_to_cube_file ( soluteChargeDensity, filename  )
    END IF

END SUBROUTINE soluteChargeDensityFromSoluteChargeCoordinates

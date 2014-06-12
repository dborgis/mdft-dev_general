! Compute mean density perpendicular to plan xy

SUBROUTINE compute_z_density (array,filename)
    
    USE precision_kinds ,ONLY: dp, i2b
    USE system          ,ONLY: spaceGrid
    USE mathematica     ,ONLY: chop
    
    IMPLICIT NONE
    
    REAL(dp) :: mean_density, z
    REAL(dp), INTENT(IN) :: array(spaceGrid%n_nodes(1),spaceGrid%n_nodes(2),spaceGrid%n_nodes(3))
    INTEGER(i2b) :: k
    CHARACTER(50), INTENT(IN):: filename
    OPEN (10, FILE=filename)
        ! Compute mean density over x and y
        DO k = 1, spaceGrid%n_nodes(3)
            z = REAL(k-1,dp) * spaceGrid%dl(3)
            mean_density = chop(SUM(array(:,:,k)) / REAL(PRODUCT(spaceGrid%n_nodes(1:2)),dp))
            WRITE (10,*)z,mean_density
        END DO
    CLOSE (10)

END SUBROUTINE compute_z_density

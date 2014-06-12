! Compute mean density perpendicular to plan xy

SUBROUTINE compute_z_density (array,filename)
    
    USE precision_kinds ,ONLY: dp, i2b
    USE system          ,ONLY: nfft1, nfft2, nfft3, spaceGrid
    USE mathematica     ,ONLY: chop
    
    IMPLICIT NONE
    
    REAL(dp) :: mean_density
    REAL(dp), DIMENSION (nfft1,nfft2,nfft3), INTENT(IN) :: array
    INTEGER(i2b) :: k
    REAL(dp) :: z
    CHARACTER(50), INTENT(IN):: filename
    
    OPEN (10, FILE=filename)
        ! Compute mean density over x and y
        mean_density=0.0_dp
        DO k=1,nfft3
            z = REAL(k-1,dp) * spaceGrid%dl(3)
            mean_density = chop(SUM(array(:,:,k)) / REAL(nfft1*nfft2,dp))
            WRITE (10,*)z,mean_density
        END DO
    CLOSE (10)

END SUBROUTINE compute_z_density

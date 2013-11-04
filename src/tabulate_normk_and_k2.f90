! This subroutine tabulates the value of norm(k) for each vector k (l,m,n) ((nfft1,nfft2,nfft3)) and its square
SUBROUTINE tabulate_normk_and_k2

    USE precision_kinds , ONLY: i2b , dp
    USE system , ONLY: spaceGrid
    USE fft , ONLY: kx , ky , kz ! what we want to tabulate

    IMPLICIT NONE

    INTEGER(i2b) :: nfft1, nfft2, nfft3
    INTEGER(i2b):: l , m , n , m1 , m2 , m3
    REAL(dp) :: lx, ly, lz
    REAL(dp), PARAMETER :: twopi = ACOS(-1._dp)*2._dp

    nfft1 = spaceGrid%n_nodes(1)
    nfft2 = spaceGrid%n_nodes(2)
    nfft3 = spaceGrid%n_nodes(3)
    lx = spaceGrid%length(1)
    ly = spaceGrid%length(2)
    lz = spaceGrid%length(3)

    ALLOCATE ( kx ( nfft1/2+1 ) )
    ALLOCATE ( ky ( nfft2 ) )
    ALLOCATE ( kz ( nfft3 ) )

    DO CONCURRENT ( n=1:nfft3 )
        IF ( n <= nfft3/2 ) THEN
            m3 = n - 1
        ELSE
            m3 = n - 1 - nfft3
        END IF
        kz(n) = twopi/Lz*real(m3,dp)
    END DO
    
    DO CONCURRENT ( m=1:nfft2 )
        IF ( m <= nfft2/2 ) THEN
            m2 = m - 1
        ELSE
            m2 = m - 1 - nfft2
        END IF
        ky(m) = twopi/Ly*real(m2,dp)
    END DO
    
    DO CONCURRENT ( l=1:nfft1/2+1 )
        IF ( l <= nfft1/2 ) THEN
            m1 = l - 1
        ELSE
            m1 = l - 1 - nfft1
        END IF
        kx(l) = twopi/Lx*real(m1,dp)
    END DO
    
END SUBROUTINE tabulate_normk_and_k2

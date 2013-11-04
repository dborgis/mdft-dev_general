! This subroutine tabulates the value of norm(k) for each vector k (l,m,n) ((nfft1,nfft2,nfft3)) and its square
SUBROUTINE tabulate_normk_and_k2

    USE precision_kinds , ONLY: i2b , dp
    USE system , ONLY: spaceGrid
    USE fft , ONLY: kx , ky , kz ! what we want to tabulate

    IMPLICIT NONE

    INTEGER(i2b) :: nfft1, nfft2, nfft3
    INTEGER(i2b):: l, m1
    REAL(dp), PARAMETER :: twopi = ACOS(-1._dp)*2._dp

    nfft1 = spaceGrid%n_nodes(1)
    nfft2 = spaceGrid%n_nodes(2)
    nfft3 = spaceGrid%n_nodes(3)

    ALLOCATE ( kx ( nfft1/2+1 ) )
    ALLOCATE ( ky ( nfft2 ) )
    ALLOCATE ( kz ( nfft3 ) )

    DO CONCURRENT ( l=1:nfft1/2+1 )
        IF ( l <= nfft1/2 ) THEN
            m1 = l - 1
        ELSE
            m1 = l - 1 - nfft1
        END IF
        kx(l) = twopi*REAL(m1,dp)/spaceGrid%length(1)
    END DO

    DO CONCURRENT ( l=1:nfft2 )
        IF ( l <= nfft2/2 ) THEN
            m1 = l - 1
        ELSE
            m1 = l - 1 - nfft2
        END IF
        ky(l) = twopi*REAL(m1,dp)/spaceGrid%length(2)
    END DO
    
    DO CONCURRENT ( l=1:nfft3 )
        IF ( l <= nfft3/2 ) THEN
            m1 = l - 1
        ELSE
            m1 = l - 1 - nfft3
        END IF
        kz(l) = twopi*REAL(m1,dp)/spaceGrid%length(3)
    END DO
    
END SUBROUTINE tabulate_normk_and_k2

! This subroutine tabulates the value of norm(k) for each vector k (l,m,n) ((nfft1,nfft2,nfft3)) and its square
SUBROUTINE tabulate_normk_and_k2

    USE precision_kinds , ONLY: i2b , dp
    USE system , ONLY: spaceGrid
    USE fft , ONLY: kx , ky , kz, kproj

    IMPLICIT NONE

    INTEGER(i2b), DIMENSION(3) :: nfft
    INTEGER(i2b):: l
    REAL(dp), PARAMETER :: twopi = ACOS(-1._dp)*2._dp

    nfft = spaceGrid%n_nodes

    ALLOCATE ( kx ( nfft(1)/2+1 ) )
    ALLOCATE ( ky ( nfft(2) ) )
    ALLOCATE ( kz ( nfft(3) ) )

    DO CONCURRENT ( l=1:nfft(1)/2+1 )
        kx(l) = kproj(1,l)
    END DO

    DO CONCURRENT ( l=1:nfft(2) )
        ky(l) = kproj(2,l)
    END DO
    
    DO CONCURRENT ( l=1:nfft(3) )
        kz(l) = kproj(3,l)
    END DO
    
END SUBROUTINE tabulate_normk_and_k2

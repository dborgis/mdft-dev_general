!> Variables needed by FFTW3.
MODULE fft
    
    USE precision_kinds , ONLY : i4b, i2b, dp

    IMPLICIT NONE

    INTEGER(i4b) :: plan_forward, plan_backward ! descriptors of our FFTs
    REAL(dp), ALLOCATABLE, DIMENSION (:,:,:) :: in_forward, out_backward ! input of FFT and output (result) of inverse FFT
    COMPLEX(dp), ALLOCATABLE , DIMENSION (:,:,:) :: out_forward, in_backward ! output (result) of FFT and input of inverse FFT
    REAL(dp), ALLOCATABLE, DIMENSION (:,:,:) :: k2_nocoef ! norm squared of vector k without any coefficient (L or nfft)
    REAL(dp), PRIVATE, PARAMETER :: twopi = 2._dp*ACOS(-1._dp)

    REAL(dp), ALLOCATABLE, DIMENSION (:) :: kx, ky, kz ! projection of k
    !~ REAL(dp), ALLOCATABLE , DIMENSION (:,:,:) :: norm_k ! norm of vector k tabulated for l,m,n  (nfft1,nfft2,nfft3)
    !~ REAL(dp), ALLOCATABLE , DIMENSION (:,:,:) :: k2 ! norm squared of vector k (nfft1,nfft2,nfft3)

CONTAINS

    PURE FUNCTION k2 (l,m,n)
        INTEGER(i2b), INTENT(IN) :: l,m,n
        REAL(dp) :: k2
        k2 = kproj(1,l)**2 + kproj(2,m)**2 + kproj(3,n)**2
    END FUNCTION k2
    
    PURE FUNCTION norm_k (l,m,n)
        INTEGER(i2b), INTENT(IN) :: l,m,n
        REAL(dp) :: norm_k
        norm_k = SQRT( k2(l,m,n) )
    END FUNCTION norm_k

    PURE FUNCTION kproj (dir,l)
        USE system, ONLY: spaceGrid
        INTEGER(i2b), INTENT(IN) :: dir, l ! dir is 1 for x, 2 for y, 3 for z
        REAL(dp) :: kproj
        INTEGER(i2b) :: m1
        IF ( l <= spaceGrid%n_nodes(dir)/2 ) THEN
            m1 = l - 1
        ELSE
            m1 = l - 1 - spaceGrid%n_nodes(dir)
        END IF
        kproj = twopi/spaceGrid%length(1)*REAL(m1,dp)
    END FUNCTION
    
    PURE FUNCTION kvec (l,m,n)
        INTEGER(i2b), INTENT(IN) :: l,m,n
        REAL(dp), DIMENSION(3) :: kvec
        kvec = [ kproj(1,l), kproj(2,m), kproj(3,n) ]
    END FUNCTION kvec

    SUBROUTINE deallocate_everything_fft
        CALL dfftw_destroy_plan ( plan_forward ) ! destroy FFTW3 plans 
        CALL dfftw_destroy_plan ( plan_backward ) ! destroy FFTW3 plans 
    END SUBROUTINE deallocate_everything_fft
    
    PURE FUNCTION timesExpPrefactork2 (array3D, prefactor)
        COMPLEX(dp), DIMENSION(:,:,:), INTENT(IN) :: array3D
        COMPLEX(dp), DIMENSION(SIZE(array3D,1),SIZE(array3D,2),SIZE(array3D,3)) :: timesExpPrefactork2
        REAL(dp), INTENT(IN) :: prefactor
        INTEGER(i2b) :: i,j,k,imax,jmax,kmax
        imax = SIZE(array3D,1)
        jmax = SIZE(array3D,2)
        kmax = SIZE(array3D,3)
        DO CONCURRENT ( i=1:imax, j=1:jmax, k=1:kmax )
            timesExpPrefactork2 (i,j,k) = array3D (i,j,k) * EXP( prefactor* k2 (i,j,k) )
        END DO
    END FUNCTION

END MODULE fft

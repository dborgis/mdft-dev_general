! this SUBROUTINE initiates the plans needed by fftw3.
! it also allocates the input and output arrays of fftw3.

SUBROUTINE init_fftw3_plans
    USE fft, ONLY: fftw3
    USE system, ONLY: spaceGrid
 
    IMPLICIT NONE

    INTEGER, DIMENSION(3) :: nfft
    INCLUDE "fftw3.f"

    ! allocate the arrays needed as input for FFT (in_forward) or output for FFT (out_forward)
    ! or needed as input for inverse FFT (in_backward) etc.
    nfft = spaceGrid%n_nodes
    ALLOCATE ( fftw3%in_forward   ( nfft(1)      , nfft(2) , nfft(3) ) )
    ALLOCATE ( fftw3%out_forward  ( nfft(1)/2 +1 , nfft(2) , nfft(3) ) )
    ALLOCATE ( fftw3%out_backward ( nfft(1)      , nfft(2) , nfft(3) ) )
    ALLOCATE ( fftw3%in_backward  ( nfft(1)/2 +1 , nfft(2) , nfft(3) ) )
    ! prepare plans needed by fftw3
    CALL dfftw_plan_dft_r2c_3d ( fftw3%plan_forward, nfft(1), nfft(2), nfft(3), fftw3%in_forward, fftw3%out_forward, FFTW_MEASURE )
    CALL dfftw_plan_dft_c2r_3d ( fftw3%plan_backward, nfft(1), nfft(2), nfft(3), fftw3%in_backward, fftw3%out_backward,&
                                                                                                                FFTW_MEASURE )
    
END SUBROUTINE init_fftw3_plans

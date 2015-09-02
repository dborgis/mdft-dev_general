! this SUBROUTINE initiates the plans needed by fftw3.
! it also allocates the input and output arrays of fftw3.

SUBROUTINE init_fftw3_plans
    use fft, ONLY: fftw3
    use module_grid, only: grid

    IMPLICIT NONE

    INTEGER, DIMENSION(3) :: nfft
    INCLUDE "fftw3.f"

    ! allocate the arrays needed as input for FFT (in_forward) or output for FFT (out_forward)
    ! or needed as input for inverse FFT (in_backward) etc.
    nfft = grid%n_nodes
    ALLOCATE ( fftw3%in_forward   ( nfft(1)      , nfft(2) , nfft(3) ) )
    ALLOCATE ( fftw3%out_forward  ( nfft(1)/2 +1 , nfft(2) , nfft(3) ) )
    ALLOCATE ( fftw3%out_backward ( nfft(1)      , nfft(2) , nfft(3) ) )
    ALLOCATE ( fftw3%in_backward  ( nfft(1)/2 +1 , nfft(2) , nfft(3) ) )

    ! prepare plans needed by fftw3
    CALL dfftw_plan_dft_r2c_3d &
        ( fftw3%plan_forward, nfft(1), nfft(2), nfft(3), fftw3%in_forward, fftw3%out_forward, FFTW_MEASURE )
    CALL dfftw_plan_dft_c2r_3d &
        ( fftw3%plan_backward, nfft(1), nfft(2), nfft(3), fftw3%in_backward, fftw3%out_backward, FFTW_MEASURE )
    ! Note about fftw planning-flags:
    ! Since lots of FFT will be done in each direction with these two plans, we use the quite rigorous planner flag (FFTW_MEASURE)
    ! to find the optimal plan. This, of course, costs a substantial planning time. We do not use more rigorous planner because
    ! it then becomes the most time-consuming step of MDFT !
    ! See http://www.fftw.org/doc/Planner-Flags.html

END SUBROUTINE init_fftw3_plans

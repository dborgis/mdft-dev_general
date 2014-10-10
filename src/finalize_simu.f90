subroutine finalize_simu

    use fft, only: finalize_fftw

    implicit none
    
    call finalize_fftw

end subroutine finalize_simu

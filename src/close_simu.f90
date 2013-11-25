! This SUBROUTINE closes the program properly.
! It deallocates memory and destroys FFTW3 plans.
SUBROUTINE close_simu

    USE fft, ONLY : deallocate_everything_fft

    IMPLICIT NONE

    CALL deallocate_everything_fft

END SUBROUTINE close_simu

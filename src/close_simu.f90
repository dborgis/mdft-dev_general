! This subroutine closes the program properly.
! It deallocates memory and destroys FFTW3 plans.
SUBROUTINE close_simu
    USE cg, ONLY: finalizeMinimizer
    USE fft, ONLY : deallocate_everything_fft
    USE external_potential, ONLY : deallocate_everything_external_potential

    IMPLICIT NONE

    CALL deallocate_everything_external_potential
    CALL deallocate_everything_fft
    CALL finalizeMinimizer

END SUBROUTINE close_simu

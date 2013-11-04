!> This subroutine closes the program properly.
!! It deallocates memory and destroys FFTW3 plans.
subroutine close_simu
use cg , only : deallocate_everything_cg
use fft , only : deallocate_everything_fft
use external_potential , only : deallocate_everything_external_potential
implicit none
call deallocate_everything_external_potential
call deallocate_everything_fft
call deallocate_everything_cg
end subroutine close_simu

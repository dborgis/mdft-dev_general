!> Variables needed by FFTW3.
module fft
use precision_kinds , only : i4b , dp
implicit none
integer ( kind = i4b ) :: plan_forward !> @var indicates kind of planed FFT
integer ( kind = i4b ) :: plan_backward !> @var indicates kind of planed FFT
real(dp), allocatable , dimension ( : , : , : ) :: in_forward !> @var Array input in FFT
complex(dp), allocatable , dimension ( : , : , : ) :: out_forward !> @var Array output in FFT
complex(dp), allocatable , dimension ( : , : , : ) :: in_backward !> @var Array input of FFT-1
real(dp), allocatable , dimension ( : , : , : ) :: out_backward !> @var Array output for FFT-1
real(dp), allocatable , dimension ( : , : , : ) :: norm_k ! norm of vector k tabulated for l,m,n  (nfft1,nfft2,nfft3)
real(dp), allocatable , dimension ( : , : , : ) :: k2 ! norm squared of vector k (nfft1,nfft2,nfft3)
real(dp), allocatable , dimension ( : , : , : ) :: k2_nocoef ! norm squared of vector k without any coefficient (L or nfft)
real(dp), allocatable , dimension ( : ) :: kx, ky, kz ! projection of k
contains
  subroutine deallocate_everything_fft
    implicit none
    if ( allocated ( in_forward ) ) deallocate ( in_forward ) ! if X is allocated then deallocate it
    if ( allocated ( out_forward ) ) deallocate ( out_forward )
    if ( allocated ( in_backward ) ) deallocate ( in_backward )
    if ( allocated ( out_backward ) ) deallocate ( out_backward )
    if ( allocated ( norm_k ) ) deallocate ( norm_k )
    if ( allocated ( k2 ) ) deallocate ( k2 )
    if ( allocated ( kx ) ) deallocate ( kx )
    if ( allocated ( ky ) ) deallocate ( ky )
    if ( allocated ( kz ) ) deallocate ( kz )
    ! destroy FFTW3 plans 
    call dfftw_destroy_plan ( plan_forward )
    call dfftw_destroy_plan ( plan_backward )
  end subroutine deallocate_everything_fft
end module fft

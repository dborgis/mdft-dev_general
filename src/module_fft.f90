module module_fft
    ! This module deals with everything related to the fast fourier transforms and all Fourier space related functions
    use precision_kinds, only: i2b, dp, i4b
    use iso_c_binding, only: c_int
    implicit none
    private

    type :: fftw3_type3d_r2c2r
        integer(i4b) :: plan_forward, plan_backward ! descriptors of our FFTs
        real(dp), allocatable, dimension(:,:,:) :: in_forward, out_backward ! input of FFT and output (result) of inverse FFT
        complex(dp), allocatable, dimension(:,:,:) :: out_forward, in_backward ! output (result) of FFT and input of inverse FFT
        integer :: nthread = 1
    end type


    type :: fftw3_type2d_r2c2r
        integer(i4b) :: plan_forward, plan_backward ! descriptors of our FFTs
        real(dp), allocatable, dimension(:,:) :: in_forward, out_backward    ! input of FFT and output (result) of inverse FFT
        complex(dp), allocatable, dimension(:,:) :: in_backward, out_forward ! output (result) of FFT and input of inverse FFT
        integer :: nthread = 1
        logical :: isok=.false.
    contains
        procedure, nopass :: init => init_fft2d
    end type

    type(fftw3_type3d_r2c2r) :: fftw3
    type(fftw3_type2d_r2c2r) :: fft2d

    public :: init, fft2d, fftw3

contains

    subroutine init_fft2d
        use iso_c_binding
        use precision_kinds, only: dp
        use module_grid, only: grid
        implicit none
        integer :: npsi, nphi
        include "fftw3.f03"
        if (fft2d%isok) then
            print*, "in init_fft2d, fft2d%isok is already ok"
            print*, "but still, i have received a call to fft2d%init"
            error stop "in module_fft > init_fft2d = fft2d%init"
        end if
        ! allocate the arrays needed as input for FFT (in_forward) or output for FFT (out_forward)
        ! or needed as input for inverse FFT (in_backward) etc.
        npsi = grid%npsi
        nphi = grid%nphi
        allocate (fft2d%in_forward  (npsi,nphi), source=0._dp)
        allocate (fft2d%out_backward(npsi,nphi), source=0._dp)
        allocate (fft2d%in_backward (npsi/2+1,nphi), source=(0._dp,0._dp))
        allocate (fft2d%out_forward (npsi/2+1,nphi), source=(0._dp,0._dp))
        ! prepare plans needed by fftw3
        ! in case of using energy_cproj, we will do a HUGE number of these fft2d.
        ! It is worth doing FFTW_EXHAUSTIVE
        select case (dp)
        case(c_double)
          call dfftw_plan_dft_r2c_2d (fft2d%plan_forward, npsi, nphi, fft2d%in_forward, fft2d%out_forward, FFTW_EXHAUSTIVE)
          call dfftw_plan_dft_c2r_2d (fft2d%plan_backward, npsi, nphi, fft2d%in_backward, fft2d%out_backward, FFTW_EXHAUSTIVE)
        case(c_float)
          call sfftw_plan_dft_r2c_2d (fft2d%plan_forward, npsi, nphi, fft2d%in_forward, fft2d%out_forward, FFTW_EXHAUSTIVE)
          call sfftw_plan_dft_c2r_2d (fft2d%plan_backward, npsi, nphi, fft2d%in_backward, fft2d%out_backward, FFTW_EXHAUSTIVE)
        end select
        fft2d%isok=.true.
    end subroutine init_fft2d

    subroutine init_threads
        use module_input, only: getinput
        implicit none
        integer :: iRet
        fftw3%nthread = getinput%int('number_of_fftw3_threads',defaultvalue=1,assert=">0")
        if (fftw3%nthread/=1) then
            print*,"Number of threads for FFTW3:", fftw3%nthread
        end if
        call dfftw_init_threads(iRet)
        if (iRet/=1) then
            print*, "Problem in dfftw_init_threads(), returned value is ",iRet
            stop
        end if
        call dfftw_plan_with_nthreads (fftw3%nthread)
    end subroutine


    ! subroutine finalize_fftw
    !     call dfftw_cleanup_threads()
    ! end subroutine finalize_fftw


    subroutine init
        call init_threads()
        call init_fftw3_plans() ! initialize FFTW3 plans
    end subroutine


    ! This subroutine initiates the plans needed by fftw3.
    ! it also allocates the input and output arrays of fftw3.
    subroutine init_fftw3_plans
        use iso_c_binding
        use precision_kinds, only: dp
        use module_grid, only: grid
        IMPLICIT NONE
        integer :: nx, ny, nz
        include "fftw3.f03"
        ! allocate the arrays needed as input for FFT (in_forward) or output for FFT (out_forward)
        ! or needed as input for inverse FFT (in_backward) etc.
        nx = grid%nx
        ny = grid%ny
        nz = grid%nz
        allocate (fftw3%in_forward(nx,ny,nz))
        allocate (fftw3%out_backward(nx,ny,nz))
        allocate (fftw3%in_backward(nx/2+1,ny,nz))
        allocate (fftw3%out_forward(nx/2+1,ny,nz))

        !
        ! prepare FFTW3 plans. http://www.fftw.org/doc/Planner-Flags.html
        !
        select case (dp)
        case(c_double)
          call dfftw_plan_dft_r2c_3d( fftw3%plan_forward, nx, ny, nz, fftw3%in_forward, fftw3%out_forward, FFTW_MEASURE )
          call dfftw_plan_dft_c2r_3d( fftw3%plan_backward, nx, ny, nz, fftw3%in_backward, fftw3%out_backward, FFTW_MEASURE )
        case(c_float)
          call sfftw_plan_dft_r2c_3d( fftw3%plan_forward, nx, ny, nz, fftw3%in_forward, fftw3%out_forward, FFTW_MEASURE )
          call sfftw_plan_dft_c2r_3d( fftw3%plan_backward, nx, ny, nz, fftw3%in_backward, fftw3%out_backward, FFTW_MEASURE )
        case default
          error stop "You are working on neither single or double precision ... ?"
        end select
    end subroutine init_fftw3_plans

end module module_fft

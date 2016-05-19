module module_fft
    ! This module deals with everything related to the fast fourier transforms and all Fourier space related functions
    use precision_kinds, only: i2b, dp, i4b
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
        use precision_kinds, only: dp
        use module_grid, only: grid
        implicit none
        integer :: npsi, nphi
        include "fftw3.f"
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
        call dfftw_plan_dft_r2c_2d (fft2d%plan_forward, npsi, nphi, fft2d%in_forward, fft2d%out_forward, FFTW_EXHAUSTIVE)
        call dfftw_plan_dft_c2r_2d (fft2d%plan_backward, npsi, nphi, fft2d%in_backward, fft2d%out_backward, FFTW_EXHAUSTIVE)
        ! in case of using energy_cproj, we will do a HUGE number of these fft2d.
        ! It is worth doing FFTW_EXHAUSTIVE
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


    SUBROUTINE init
        call init_threads()
        call init_fftw3_plans() ! initialize FFTW3 plans
    END SUBROUTINE


    ! this SUBROUTINE initiates the plans needed by fftw3.
    ! it also allocates the input and output arrays of fftw3.
    SUBROUTINE init_fftw3_plans
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
        ! prepare plans needed by fftw3
        CALL dfftw_plan_dft_r2c_3d &
        ( fftw3%plan_forward, nx, ny, nz, fftw3%in_forward, fftw3%out_forward, FFTW_MEASURE )
        CALL dfftw_plan_dft_c2r_3d &
        ( fftw3%plan_backward, nx, ny, nz, fftw3%in_backward, fftw3%out_backward, FFTW_MEASURE )
        ! Note about fftw planning-flags:
        ! Since lots of FFT will be done in each direction with these two plans, we use the quite rigorous planner flag (FFTW_MEASURE)
        ! to find the optimal plan. This, of course, costs a substantial planning time. We do not use more rigorous planner because
        ! it then becomes the most time-consuming step of MDFT !
        ! See http://www.fftw.org/doc/Planner-Flags.html
    END SUBROUTINE init_fftw3_plans


    ! SUBROUTINE deallocate_everything_fft
    !     CALL dfftw_destroy_plan ( fftw3%plan_forward ) ! destroy FFTW3 plans
    !     CALL dfftw_destroy_plan ( fftw3%plan_backward ) ! destroy FFTW3 plans
    ! END SUBROUTINE deallocate_everything_fft

end module module_fft

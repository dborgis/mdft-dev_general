module fft
    ! This module deals with everything related to the fast fourier transforms and all Fourier space related functions
    use precision_kinds, only: i2b, dp, i4b
    implicit none
    private
    type :: fftw3type
        INTEGER(i4b) :: plan_forward, plan_backward ! descriptors of our FFTs
        REAL(dp), allocatable, dimension(:,:,:) :: in_forward, out_backward ! input of FFT and output (result) of inverse FFT
        COMPLEX(dp), allocatable, dimension(:,:,:) :: out_forward, in_backward ! output (result) of FFT and input of inverse FFT
        integer :: nthread
    end type
    type(fftw3type), public :: fftw3
    public :: init

contains

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


    subroutine finalize_fftw
        call dfftw_cleanup_threads()
    end subroutine finalize_fftw


    SUBROUTINE init
        call init_threads()
        call init_fftw3_plans() ! initialize FFTW3 plans
    END SUBROUTINE


    ! this SUBROUTINE initiates the plans needed by fftw3.
    ! it also allocates the input and output arrays of fftw3.
    SUBROUTINE init_fftw3_plans
        use module_grid, only: grid
        IMPLICIT NONE
        integer :: nx, ny, nz
        INCLUDE "fftw3.f"
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


    SUBROUTINE deallocate_everything_fft
        CALL dfftw_destroy_plan ( fftw3%plan_forward ) ! destroy FFTW3 plans
        CALL dfftw_destroy_plan ( fftw3%plan_backward ) ! destroy FFTW3 plans
    END SUBROUTINE deallocate_everything_fft

END MODULE fft

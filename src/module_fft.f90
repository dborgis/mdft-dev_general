!===================================================================================================================================
MODULE fft
!===================================================================================================================================
! This module deals with everything related to the fast fourier transforms and all Fourier space related functions

    USE precision_kinds ,ONLY : i2b, dp, i4b

    IMPLICIT NONE

    TYPE :: fftw3Needs
        INTEGER(i4b) :: plan_forward, plan_backward ! descriptors of our FFTs
        REAL(dp), ALLOCATABLE, DIMENSION (:,:,:) :: in_forward, out_backward ! input of FFT and output (result) of inverse FFT
        COMPLEX(dp), ALLOCATABLE, DIMENSION (:,:,:) :: out_forward, in_backward ! output (result) of FFT and input of inverse FFT
    END TYPE
    TYPE( fftw3Needs ) :: fftw3
    REAL(dp), PRIVATE, PARAMETER :: twopi = 2._dp*ACOS(-1._dp)
    REAL(dp), ALLOCATABLE, DIMENSION (:) :: kx, ky, kz ! projection of k

    CONTAINS

    !===============================================================================================================================
    subroutine init_threads
        use input, only: getinput
        implicit none
        integer :: n_threads, iRet
        n_threads = getinput%int('number_of_fftw3_threads',1)
        print*,"Number of threads for FFTW3:", n_threads
        call dfftw_init_threads(iRet)
            if (iRet/=1) then
                print*, "Problem in dfftw_init_threads(), returned value is ",iRet
                stop
            end if
        call dfftw_plan_with_nthreads (n_threads)
    end subroutine
    !===============================================================================================================================


    !===============================================================================================================================
    subroutine finalize_fftw
        implicit none
        call dfftw_cleanup_threads()
    end subroutine finalize_fftw
    !===============================================================================================================================


    !===============================================================================================================================
    SUBROUTINE init
    !===============================================================================================================================
        CALL init_fftw3_plans! initialize FFTW3 plans
        CALL tabulate_kx_ky_kz! tabulate the value of norm(\vec{k}) and k2 = norm_k**2 in order to speed up program
    END SUBROUTINE
    !===============================================================================================================================


    !===============================================================================================================================
    SUBROUTINE tabulate_kx_ky_kz
    !===============================================================================================================================
        use module_grid, only: grid
        INTEGER(i2b), DIMENSION(3) :: nfft
        INTEGER(i2b):: l
        nfft = grid%n_nodes
        ALLOCATE ( kx (nfft(1)/2+1), ky (nfft(2)), kz (nfft(3)), SOURCE=0._dp)
        DO CONCURRENT ( l=1:nfft(1)/2+1 )
            kx(l) = kproj(1,l)
        END DO
        DO CONCURRENT ( l=1:nfft(2) )
            ky(l) = kproj(2,l)
        END DO
        DO CONCURRENT ( l=1:nfft(3) )
            kz(l) = kproj(3,l)
        END DO
    END SUBROUTINE
    !===============================================================================================================================


    !===============================================================================================================================
    PURE FUNCTION k2 (l,m,n)
    !===============================================================================================================================
        INTEGER(i2b), INTENT(IN) :: l,m,n
        REAL(dp) :: k2
        k2 = kx(l)**2 + ky(m)**2 + kz(n)**2
    END FUNCTION k2
    !===============================================================================================================================


    !===============================================================================================================================
    PURE FUNCTION norm_k (l,m,n)
    !===============================================================================================================================
        INTEGER(i2b), INTENT(IN) :: l,m,n
        REAL(dp) :: norm_k
        norm_k = SQRT( k2(l,m,n) )
    END FUNCTION norm_k
    !===============================================================================================================================


    !===============================================================================================================================
    PURE FUNCTION kproj (dir,l)
    !===============================================================================================================================
    ! note the special ordering for negative values. See
    ! http://www.fftw.org/doc/Real_002ddata-DFT-Array-Format.html#Real_002ddata-DFT-Array-Format
    ! http://www.fftw.org/doc/The-1d-Discrete-Fourier-Transform-_0028DFT_0029.html#The-1d-Discrete-Fourier-Transform-_0028DFT_0029
        use module_grid, only: grid
        INTEGER(i2b), INTENT(IN) :: dir, l ! dir is 1 for x, 2 for y, 3 for z
        REAL(dp) :: kproj
        INTEGER(i2b) :: m1
        IF ( l <= grid%n_nodes(dir)/2 ) THEN
            m1 = l - 1
        ELSE
            m1 = l - 1 - grid%n_nodes(dir)
        END IF
        kproj = twopi/grid%length(dir)*REAL(m1,dp)
    END FUNCTION
    !===============================================================================================================================


    !===============================================================================================================================
    PURE FUNCTION kvec (l,m,n)
    !===============================================================================================================================
        INTEGER(i2b), INTENT(IN) :: l,m,n
        REAL(dp), DIMENSION(3) :: kvec
        kvec(1) = kproj(1,l)
        kvec(2) = kproj(2,m)
        kvec(3) = kproj(3,n)
    END FUNCTION kvec
    !===============================================================================================================================


    !===============================================================================================================================
    SUBROUTINE deallocate_everything_fft
    !===============================================================================================================================
        CALL dfftw_destroy_plan ( fftw3%plan_forward ) ! destroy FFTW3 plans
        CALL dfftw_destroy_plan ( fftw3%plan_backward ) ! destroy FFTW3 plans
    END SUBROUTINE deallocate_everything_fft
    !===============================================================================================================================


    !===============================================================================================================================
    PURE FUNCTION timesExpPrefactork2 (array3D, prefactor)
    !===============================================================================================================================
        COMPLEX(dp), DIMENSION(:,:,:), INTENT(IN) :: array3D
        COMPLEX(dp), DIMENSION(SIZE(array3D,1),SIZE(array3D,2),SIZE(array3D,3)) :: timesExpPrefactork2
        REAL(dp), INTENT(IN) :: prefactor
        INTEGER(i2b) :: i,j,k,imax,jmax,kmax
        imax = SIZE(array3D,1)
        jmax = SIZE(array3D,2)
        kmax = SIZE(array3D,3)
        DO CONCURRENT ( i=1:imax, j=1:jmax, k=1:kmax )
            timesExpPrefactork2 (i,j,k) = array3D (i,j,k) * EXP( prefactor* k2 (i,j,k) )
        END DO
    END FUNCTION
    !===============================================================================================================================


END MODULE fft

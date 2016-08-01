module module_read_c_luc
    implicit none
    private
    public :: read_c_luc
contains
    subroutine read_c_luc( cmnmunukhi, mmax, mrso, qmaxwanted, np, nq, m, n, mu, nu, khi, p )
        use precision_kinds, only: dp
        use module_input, only: n_linesInFile
        implicit none
        complex(dp), allocatable :: cmnmunukhi(:,:)
        character(len=40), parameter :: filename(0:5) = &
        ["input/dcf/water/SPCE/ck_nonzero_nmax0_ml",&
         "input/dcf/water/SPCE/ck_nonzero_nmax1_ml",&
         "input/dcf/water/SPCE/ck_nonzero_nmax2_ml",&
         "input/dcf/water/SPCE/ck_nonzero_nmax3_ml",&
         "input/dcf/water/SPCE/ck_nonzero_nmax4_ml",&
         "input/dcf/water/SPCE/ck_nonzero_nmax5_ml" ]
        integer, intent(out) :: np(0:5)
        integer, intent(in) :: mmax
        integer, intent(in) :: mrso ! Symetry of the main axis. For instance, mrso=2 for a molecule of symetry C2v (like water).
        integer, intent(out) :: nq
        real(dp), intent(in) :: qmaxwanted
        integer, intent(out), allocatable :: m(:), n(:), mu(:), nu(:), khi(:), p(:,:,:,:,:)
        ! Number of projections for mmax = 0 to 5
        np(0:5) = [1,6,75,252,877,2002]
        ! Inquire that the file exists
        block
            logical :: exist
            inquire(file=filename(mmax), exist=exist)
            if( .not. exist) then
                error stop "In module_read_c_luc/read_c_luc, filename does not exist"
            end if
        end block
        ! Verify it contains 1024 values of q
        block
            integer :: nline
            nline = n_linesInFile(filename(mmax)) -17 ! header contains 17 lines
            if( nline /= nq) then
                error stop "In module_read_c_luc/read_c_luc, filename does not contain 1024 values of q"
            end if
        end block
        ! Open the file. Its unit will be 88.
        block
            integer :: ios
            open(88, file=filename(mmax), iostat=ios, status="old", action="read")
            if( ios/=0) error stop "Cant open file in module_read_c_luc"
        end block
        ! Prepare the arrays that will contain the value of m, n, mu, nu and khi
        block
            integer :: npm
            npm = np(mmax)
            allocate( m(npm), source=-huge(1)) ! -huge(1) makes it easier to detect problems later.
            allocate( n(npm), source=-huge(1)) ! since that's a huge out of bound.
            allocate( mu(npm), source=-huge(1))
            allocate( nu(npm), source=-huge(1))
            allocate( khi(npm), source=-huge(1))
            allocate( p(0:mmax,0:mmax,-mmax/mrso:mmax/mrso,-mmax/mrso:mmax/mrso,-mmax:mmax) , source=-huge(1))
        end block
        ! Skip 10 lines of comments
        block
            integer :: i
            do i=1,10
                read(88,*)
            end do
        end block
        ! Read a comment line, then the whole arrays of m, n, mu, nu and khi, one per line.
        block
            character(3) :: somechar
            read(8,*) somechar
            read(8,*) somechar, m ! m, n, mu, nu and khi are arrays of dimension np(mmax)
            read(8,*) somechar, n
            read(8,*) somechar, mu
            read(8,*) somechar, nu
            read(8,*) somechar, khi
            read(8,*)
        end block
        ! Fill the array that hash map from a projection index to m, n, mu, nu, khi
        block
            integer :: ip, im, in, imu, inu, ikhi
            do ip=1,np(mmax)
                im = m(ip)
                in = n(ip)
                imu = mu(ip) ! mu(:) and nu(:) contains the read value of mu and nu
                inu = nu(ip)
                ikhi = khi(ip)
                p(im,in,imu/mrso,inu/mrso,ikhi) = ip ! but p(:,:,:,:,:) uses mu/mrso and nu/mrso, we often call them mu2 and nu2
            end do
        end block
        ! Read q, cmnmunukhi(q)
        block
            integer, parameter  :: nqinfile = 1024
            integer :: iq
            real(dp) :: q(nqinfile)
            allocate( cmnmunukhi(np(mmax),nq) ,source=(0._dp,0._dp)  )
            do iq=1,nq
                read(88,*) q(iq), cmnmunukhi(np(mmax),iq)
                if( q(iq)>qmaxwanted ) then
                    nq = iq
                    exit
                end if
            end do
        end block
        close(88)
    end subroutine read_c_luc

end module module_read_c_luc

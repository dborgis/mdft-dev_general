module module_read_c_luc
    implicit none
    private
    public :: read_c_luc
contains
    subroutine read_c_luc( cmnmunukhi, mmax, mrso, qmaxwanted, np, nq, dq, m, n, mu, nu, khi, p )
        use precision_kinds, only: dp
        use module_input, only: n_linesInFile
        use module_solvent, only: solvent
        implicit none
        complex(dp), allocatable, intent(out) :: cmnmunukhi(:,:)
        character(len=4) :: endOfFilename
        character(len=100) :: filename
        integer, intent(in) :: mmax
        integer, intent(in) :: mrso ! Symetry of the main axis. For instance, mrso=2 for a molecule of symetry C2v (like water).
        real(dp), intent(in) :: qmaxwanted
        integer, intent(out) :: np
        integer, intent(out) :: nq
        real(dp), intent(out) :: dq
        integer, intent(out), allocatable :: m(:), n(:), mu(:), nu(:), khi(:), p(:,:,:,:,:)
        integer :: npluc(0:6) != [1,6,75,252,877,2002] ! $ grep alpha input/dcf/water/SPCE/ck_nonzero_nmax5_ml

        npluc=solvent(1)%npluc

        write (endOfFilename, "(I1,A3)") mmax,"_ml"
        filename = "data/dcf/" // trim(solvent(1)%name) // "/" // trim(solvent(1)%name) // "-ck_nonzero_nmax" // trim(endOfFilename)
		
        ! c(mnmunukhi) should not be already allocated
        if( allocated(cmnmunukhi)) error stop "c(m,n,mu,nu,khi) is already allocated in module_read_c_luc"
        ! Number of projections for mmax = 0 to 5
        np = npluc(mmax)
        print*, 'np when read' , np
        ! Inquire that the file exists
        block
            logical :: exist
            inquire(file=filename, exist=exist)
            if( .not. exist) then
                print*, "Entering module_read_c_luc/read_c_luc. You're looking for file:", filename
                print*, "But this filename does not exists"
                print*, "For TIP3P REMEMBER YOU CAN ONLY USE MMAX5 FOR NOW. DONT LOOK FOR FILENAME WITH NMAX 3 FOR INSTANCE"
                ERROR STOP
            end if
        end block
        ! Verify it contains the expected values of q
        block
            integer :: nline
            nline = n_linesInFile(filename) -17 ! header contains 17 lines
            if( nline /= solvent(1)%n_line_cfile) then
                print*, "WARNING:!!!!! len of c file not verified here !!!!!!      module_read_c_luc.f90 at line 44"
                !error stop "In module_read_c_luc/read_c_luc, filename does not contain the expected number of  values of q"
            end if
        end block
        ! Open the file. Its unit will be 88.
        block
            integer :: ios
            open(88, file=filename, iostat=ios, status="old", action="read")
            if( ios/=0) error stop "Cant open file in module_read_c_luc"
        end block
        ! Prepare the arrays that will contain the value of m, n, mu, nu and khi
        block
            allocate( m(np), source=-huge(1)) ! -huge(1) makes it easier to detect problems later.
            allocate( n(np), source=-huge(1)) ! since that's a huge out of bound.
            allocate( mu(np), source=-huge(1))
            allocate( nu(np), source=-huge(1))
            allocate( khi(np), source=-huge(1))
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
            character(len=6) :: somechar
            read(88,*) somechar
            read(88,*) somechar, m(1:np) ! m, n, mu, nu and khi are arrays of dimension np(mmax)
            read(88,*) somechar, n(1:np)
            read(88,*) somechar, mu(1:np)
            read(88,*) somechar, nu(1:np)
            read(88,*) somechar, khi(1:np)
            read(88,*)
        end block
        ! Fill the array that hash map from a projection index to m, n, mu, nu, khi
        block
            integer :: ip, im, in, imu, inu, ikhi
            do ip=1,np
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
            integer  :: nqinfile 
            real(dp) :: q, actualdq
            integer :: iq
            nqinfile = solvent(1)%n_line_cfile
            dq = solvent(1)%dq
            nq = int(qmaxwanted/dq +0.5)+2
            if( nq>nqinfile) error stop "You want more values of q that are available in module_read_c_luc"
            allocate( cmnmunukhi(np,nq) ,source=(0._dp,0._dp)  )
            do iq=1,nq ! if you want more than available, use all that is available. If you need less, use less.
                read(88,*) q, cmnmunukhi(1:np,iq)
                if (iq==1) actualdq=q
                if (iq==2) actualdq=q-actualdq
                if( q>(qmaxwanted+2._dp*dq) ) error stop "q>qmaxwanted in module_read_c_luc" ! One needs a little bit of tolerance: We must have qmaxwanted within our bins.
            end do
            if (actualdq/=dq) stop("the dq in you cfilesdoes correspond to the one assumed by the code it might be due to single precision so please check")
            block
                open(89,file="output/arraysinmemory_module_read_c_luc")
                write(89,*)"From module_read_c_luc:"
                write(89,*)"cmnmunukhi(np,nq)",size(cmnmunukhi)*storage_size(cmnmunukhi),"Bytes"
                write(89,*)"m(np)",size(m)*storage_size(m),"Bytes"
                write(89,*)"n(np)",size(n)*storage_size(n),"Bytes"
                write(89,*)"mu(np)",size(mu)*storage_size(mu),"Bytes"
                write(89,*)"nu(np)",size(nu)*storage_size(nu),"Bytes"
                write(89,*)"khi(np)",size(khi)*storage_size(khi),"Bytes"
                write(89,*)"p(0:mmax,0:mmax,-mmax/mrso:mmax/mrso,-mmax/mrso:mmax/mrso,0:mmax)",size(p)*storage_size(p),"Bytes"
                close(89)
            end block
        end block
        close(88)
    end subroutine read_c_luc

end module module_read_c_luc

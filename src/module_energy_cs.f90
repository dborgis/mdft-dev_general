module module_energy_cs
    implicit none
    private
    public :: energy_cs
contains
    subroutine energy_cs (fexc, df)

        use precision_kinds, only: dp
        use module_thermo, only: thermo
        use module_solvent, only: solvent, print_solvent_not_allocated
        use module_fft, only: fftw3
        use module_grid, only: grid

        implicit none
        integer :: nx, ny, nz, no, ns, ix, iy, iz, is, io, ik, ikmax_incsin
        real(dp), intent(out) :: Fexc
        real(dp) :: kT, dV, ksqmax, k, dk
        real(dp) :: kxsq(grid%nx/2+1), kysq(grid%ny), kzsq(grid%nz)
        real(dp), intent(inout) :: dF(:,:,:,:,:)
        real(dp), allocatable :: delta_n(:,:,:)

        if (.not.allocated(solvent)) call print_solvent_not_allocated ("Dans module_energy_cs")
        if (solvent(1)%nspec >1) then
            print*, "in energy_cs, you want several species but that's not implemented"
            error stop
        end if

        nx = grid%nx
        ny = grid%ny
        nz = grid%nz
        no = grid%no
        ns = solvent(1)%nspec
        kT = thermo%kbT
        dV = grid%dV

        !
        ! fill array fftw3 with deltan = n(x,y,z)-n0
        ! with n(x,y,z) built quadrature over orientational
        !
        allocate (delta_n(nx,ny,nz), source=0._dp)
        do iz=1,nz
            do iy=1,ny
                do ix=1,nx
                    delta_n(ix,iy,iz) = sum((solvent(1)%xi(:,ix,iy,iz)**2)*solvent(1)%rho0*grid%w(:))  - solvent(1)%n0 ! =n(r)-n0=
                end do
            end do
        end do
        fftw3%in_forward = delta_n

        !
        ! FFT(deltan(x)) => deltan(k)
        !
        call dfftw_execute( fftw3%plan_forward )

        !
        ! Now we should take care of cs(k). Do we already loaded it?
        !
        if ( .not. solvent(1)%cs%isok ) call read_cs

        !
        ! Since the distance between two points is constant in cs
        !
        dk = solvent(1)%cs%x(2) - solvent(1)%cs%x(1)
        if (dk<=epsilon(1._dp)) error stop "In module_energy_cs, dk dans cs(k) <= eps machine"

        !
        ! multiply c_s(k) by deltan(k) and fill in_backward of fft by it
        !
        ksqmax = solvent(1)%cs%x(ubound(solvent(1)%cs%x,1))**2
        ikmax_incsin = ubound(solvent(1)%cs%x,1)
        kxsq(1:nx/2+1) = [( grid%kx(ix)**2  ,ix=1,nx/2+1)]
        kysq(1:ny)     = [( grid%ky(iy)**2  ,iy=1,ny)]
        kzsq(1:nz)     = [( grid%kz(iz)**2  ,iz=1,nz)]

        if (maxval(kxsq) + maxval(kysq) + maxval(kzsq) > ksqmax) then
            print*, "Dans module_energy_cs > energy_cs, on a tellement de points k que le plus grand |k| est plus grand"
            print*, "que celui du fichier c(|k|) en input"
            print*, "k max in code =",sqrt(maxval(kxsq) + maxval(kysq) + maxval(kzsq))
            print*, "k max in cs.in =", sqrt(ksqmax)
            error stop
        end if

        fftw3%in_backward = complex(0._dp,0._dp)
        do iz=1,nz
            do iy=1,ny
                do ix=1,nx/2+1
                    k = sqrt(kxsq(ix) + kysq(iy) + kzsq(iz))
                    ik = int(k/dk) +1
                    fftw3%in_backward(ix,iy,iz) = fftw3%out_forward(ix,iy,iz) * solvent(1)%cs%y(ik)
                end do
            end do
        end do
        ! at this point we have gamma(k) in fftw3%in_backward

        call dfftw_execute( fftw3%plan_backward ) ! this fill fftw3%out_backward with gamma(x)
        fftw3%out_backward = fftw3%out_backward  /real(nx*ny*nz,dp) ! gamma(x)  note: the normalization factor /real(nx,ny,nz) comes from the way FFTW3 does NOT normalize its FFT

        ! excess free energy
        ! in_forward still contains deltan(x)=n(x)-n0
        ! out_forward contains gamma(x) = convolution of cs(x) and deltan(x)
        ! we want the integral of deltan(x)*gamma(x)
        Fexc = -kT/2._dp*dv*sum (delta_n*fftw3%out_backward)
        ! print*, "       Fexc_cs in real space, vs Fexc_cs in k space:",Fexc, -kT/2._dp*dv*sum(fftw3%in_backward*fftw3%out_forward)
        ! stop

        do is=1,ns
          do iz=1,nz
            do iy=1,ny
              do ix=1,nx
                do io=1,no
                  df(io,ix,iy,iz,is) = df(io,ix,iy,iz,is) &
                      -kT*grid%w(io)*fftw3%out_backward(ix,iy,iz)*2._dp*solvent(is)%rho0*solvent(is)%xi(io,ix,iy,iz)
                end do
              end do
            end do
          end do
        end do

    end subroutine energy_cs



    SUBROUTINE read_cs ! c(k)
        use precision_kinds, only: dp
        use module_solvent, only: solvent
        use module_input, only: n_linesInFile
        implicit none
        integer :: ios, i, is, s, nb_k
        real(dp), parameter :: onedp=1._dp, zerodp=0._dp
        type cfile_type
            real(dp), allocatable :: x(:), y(:), y2(:)
            character(150) :: filename
        end type cfile_type
        type(cfile_type) :: c_s

        if (.not. allocated(solvent)) then
            print*, "In read_cs, solvent(:)% is not allocated"
            print*, "It should be initiated before!"
            stop
        end if

        do is=1,solvent(1)%nspec
            if (solvent(1)%cs%isok) then
                print*, "I am here in module_dcf > read_cs to read cs(k)"
                print*, "but solvent(",s,")%cs%isok is already .true."
                error stop
            end if
        end do

        ! if (any(solvent%is_initiated.eq..false.))  then
        !     print*, "In read_cs, solvent has already been initiated"
        !     print*, "This is a bug"
        !     stop
        ! end if
        do s=1,solvent(1)%nspec
            if (.not.solvent(s)%is_initiated) then
                print*, "In read_cs, solvent has already been initiated"
                print*, "this is a bug"
                stop
            end if
        end do

        if (solvent(1)%nspec/=1) then
            print*, "I am in read_cs"
            print*, "This subroutine is valid for 1 solvent species at most"
            print*, "You have ",solvent(1)%nspec,"species"
            stop
        end if

        do is=1, solvent(1)%nspec
            select case (solvent(is)%name)
            case ("spce")
                c_s%filename = 'input/direct_correlation_functions/water/SPCE/cs.in'
            case ("spc")
                c_s%filename = 'input/direct_correlation_functions/water/SPC_Lionel_Daniel/cs.in'
            case ("stockmayer")
                c_s%filename = 'input/direct_correlation_functions/stockmayer/cs.in'
            case ("perso")
                c_s%filename = 'input/cs.in'
            case default
                print*, "solvent seems to be, from solvent.in", solvent(1)%name
                stop "this is not understood by module_dcf"
            end select

            nb_k = n_linesInFile(c_s%filename)

            allocate( c_s%x(nb_k) , source=zerodp)
            allocate( c_s%y(nb_k) , source=zerodp)
            allocate( c_s%y2(nb_k), source=zerodp)


            ! read c(k) as given by user and print to output folder
            OPEN (13, FILE=c_s%filename, IOSTAT=ios)
            if (ios/=0) then
                write(*,*) 'Cant open file ',c_s%filename,' in read_cs(c_s)'
                error stop
            end if
            open (14, file='output/cs.in', iostat=ios)
            if (ios/=0) stop 'Cant open file output/cs.in in read_cs(c_s)'
            do i=1,nb_k
                READ (13,*,IOSTAT=ios) c_s%x(i), c_s%y(i)
                IF (ios/=0) THEN
                    WRITE(*,*)'Error while reading line',i,"of",c_s%filename, 'in read_cs(c_s)'
                    STOP
                END IF
                WRITE(14,*,IOSTAT=ios) c_s%x(i), c_s%y(i)
                if (ios/=0) then
                    print*,'Something is wrong while writing c_s%x and c_s%y in read_cs'
                    print*,'for i=',i
                    print*,'c_s%x(i)=',c_s%x(i)
                    print*,'and c_s%y(i)=',c_s%y(i)
                    stop
                end if
            END DO
            CLOSE(13)
            close(14)

            solvent(is)%cs%filename = c_s%filename
            allocate (solvent(is)%cs%x(nb_k) ,source=c_s%x)
            allocate (solvent(is)%cs%y(nb_k) ,source=c_s%y)
            solvent(is)%cs%isok = .true.

            deallocate(c_s%x)
            deallocate(c_s%y)
            deallocate(c_s%y2)

        end do ! the loop over all solvent species

    END SUBROUTINE read_cs
end module module_energy_cs

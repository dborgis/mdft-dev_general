module module_energy_cdeltacd
    !
    ! See Borgis et al., doi:10.1103/PhysRevE.66.031206
    !
    implicit none
    private
    public :: energy_cdeltacd
contains
    subroutine energy_cdeltacd (fexc, df)
        use precision_kinds, only: dp
        use module_thermo, only: thermo
        use module_solvent, only: solvent, print_solvent_not_allocated
        use module_fft, only: fftw3
        use module_grid, only: grid
        use module_dcf, only: init_dcf

        implicit none
        integer :: nx, ny, nz, no, ix, iy, iz, is, io, ik, ikmax
        real(dp), intent(out) :: Fexc
        real(dp) :: kT, dv, ksq, k, dk, vexc, kmax, xi, rho0
        complex(dp) :: kP
        real(dp) :: kx(grid%nx/2+1), ky(grid%ny), kz(grid%nz)
        real(dp) :: kxsq(grid%nx/2+1), kysq(grid%ny), kzsq(grid%nz)
        real(dp), intent(inout) :: dF(:,:,:,:,:)
        real(dp), allocatable :: Px(:,:,:), Py(:,:,:), Pz(:,:,:)
        complex(dp), allocatable :: Px_k(:,:,:), Py_k(:,:,:), Pz_k(:,:,:)
        complex(dp), allocatable :: Ex_k(:,:,:), Ey_k(:,:,:), Ez_k(:,:,:)
        real(dp), allocatable :: Ex(:,:,:), Ey(:,:,:), Ez(:,:,:)
        real(dp), parameter :: epsdp=epsilon(1._dp)
        real(dp), parameter :: eightpisq=8._dp*acos(-1._dp)**2

        nx=grid%nx
        ny=grid%ny
        nz=grid%nz
        no=grid%no
        dv=grid%dv
        kT=thermo%kbT

        if (.not.allocated(solvent)) call print_solvent_not_allocated("Dans module_energy_cdeltacd")
        rho0=solvent(1)%rho0

        !
        ! Now we should take care of cdelta(k) and cd(k). Do we already loaded them?
        !
        if (.not. solvent(1)%cdelta%isok .or. .not. solvent(1)%cd%isok) call read_cdeltacd

        !
        ! We build the polarization vector field Px, Py, Pz
        !
        allocate (Px(nx,ny,nz), Py(nx,ny,nz), Pz(nx,ny,nz), source=0._dp)
        if (.not.allocated(grid%OMx)) then
            print*, "OMx is not allocated in module_energy_cdeltacd"
            error stop
        end if

        do iz=1,nz
            do iy=1,ny
                do ix=1,nx
                    Px(ix,iy,iz) =sum(grid%w(:)*grid%OMx(:)*solvent(1)%xi(:,ix,iy,iz)**2*rho0)
                    Py(ix,iy,iz) =sum(grid%w(:)*grid%OMy(:)*solvent(1)%xi(:,ix,iy,iz)**2*rho0)
                    Pz(ix,iy,iz) =sum(grid%w(:)*grid%OMz(:)*solvent(1)%xi(:,ix,iy,iz)**2*rho0)
                end do
            end do
        end do


        !
        ! Fourier transform the polarization vector field
        !
        fftw3%in_forward = Px
        deallocate (Px)
        call dfftw_execute (fftw3%plan_forward)
        allocate (Px_k(nx/2+1,ny,nz) ,source=fftw3%out_forward)

        fftw3%in_forward = Py
        deallocate (Py)
        call dfftw_execute (fftw3%plan_forward)
        allocate (Py_k(nx/2+1,ny,nz) ,source=fftw3%out_forward)

        fftw3%in_forward = Pz
        deallocate (Pz)
        call dfftw_execute (fftw3%plan_forward)
        allocate (Pz_k(nx/2+1,ny,nz) ,source=fftw3%out_forward)

        !
        ! Build excess electrical field in Fourier space
        !
        allocate (Ex_k(nx/2+1,ny,nz), source=(0._dp,0._dp))
        allocate (Ey_k(nx/2+1,ny,nz), source=(0._dp,0._dp))
        allocate (Ez_k(nx/2+1,ny,nz), source=(0._dp,0._dp))

        kx(1:nx/2+1) = [( grid%kx(ix)  ,ix=1,nx/2+1)]
        ky(1:ny)     = [( grid%ky(iy)  ,iy=1,ny)]
        kz(1:nz)     = [( grid%kz(iz)  ,iz=1,nz)]

        kxsq = kx**2
        kysq = ky**2
        kzsq = kz**2

        !
        ! Since the distance between two points is constant in cs
        !
        dk = solvent(1)%cdelta%x(2) - solvent(1)%cdelta%x(1)
        if (dk<=epsilon(1._dp)) error stop "In module_energy_cs, dk dans cdelta(k) <= eps machine"
        ikmax =size(solvent(1)%cdelta%x)
        kmax =solvent(1)%cdelta%x(ikmax)

        do iz=1,nz
            do iy=1,ny
                do ix=1,nx/2+1
                    ksq=kxsq(ix)+kysq(iy)+kzsq(iz)
                    if (ksq>epsdp) then
                        kP=(kx(ix)*Px_k(ix,iy,iz)+ky(iy)*Py_k(ix,iy,iz)+kz(iz)*Pz_k(ix,iy,iz))/ksq
                    else
                        kP=complex(0._dp,0._dp)
                    end if
                    k=sqrt(ksq)
                    ik=int(k/dk)+1
                    if (ik>ikmax) ik=ikmax
                    Ex_k(ix,iy,iz)=solvent(1)%cdelta%y(ik)*Px_k(ix,iy,iz)+solvent(1)%cd%y(ik)*(3._dp*kP*kx(ix)-Px_k(ix,iy,iz))
                    Ey_k(ix,iy,iz)=solvent(1)%cdelta%y(ik)*Py_k(ix,iy,iz)+solvent(1)%cd%y(ik)*(3._dp*kP*ky(iy)-Py_k(ix,iy,iz))
                    Ez_k(ix,iy,iz)=solvent(1)%cdelta%y(ik)*Pz_k(ix,iy,iz)+solvent(1)%cd%y(ik)*(3._dp*kP*kz(iz)-Pz_k(ix,iy,iz))
                end do
            end do
        end do

        !
        ! Fourier transform back the excess electrical field in real space
        !
        fftw3%in_backward = Ex_k
        deallocate (Ex_k)
        call dfftw_execute (fftw3%plan_backward)
        allocate (Ex(nx,ny,nz), source=fftw3%out_backward/real(nx*ny*nz,dp))

        fftw3%in_backward = Ey_k
        deallocate (Ey_k)
        call dfftw_execute (fftw3%plan_backward)
        allocate (Ey(nx,ny,nz), source=fftw3%out_backward/real(nx*ny*nz,dp))

        fftw3%in_backward = Ez_k
        deallocate (Ez_k)
        call dfftw_execute (fftw3%plan_backward)
        allocate (Ez(nx,ny,nz), source=fftw3%out_backward/real(nx*ny*nz,dp))


        fexc=0._dp
        is=1
        do ix=1,nx
            do iy=1,ny
                do iz=1,nz
                    do io=1,no
                        xi=solvent(1)%xi(io,ix,iy,iz)
                        vexc=-kT*grid%w(io)*(grid%OMx(io)*Ex(ix,iy,iz)+grid%OMy(io)*Ey(ix,iy,iz)+grid%OMz(io)*Ez(ix,iy,iz))
                        fexc=fexc+(rho0*xi**2-rho0)*0.5_dp*vexc*dv
                        df(io,ix,iy,iz,is)=df(io,ix,iy,iz,is)+vexc*2._dp*rho0*xi
                    end do
                end do
            end do
        end do

        deallocate (Ex,Ey,Ez)

    end subroutine energy_cdeltacd



    subroutine read_cdeltacd
        use precision_kinds, only: dp
        use module_solvent, only: solvent
        use module_input, only: n_linesInFile
        implicit none
        integer :: nk, i, ios
        if (solvent(1)%nspec/=1) then
            print*, "In read_cdeltacd, nspec>1 not implemented"
            error stop
        end if
        if (solvent(1)%cd%isok .or. solvent(1)%cdelta%isok) then
            print*, "In read_cdeltacd, cd and cdelta seem to be already ok solvent%cd%isok and solvent%cdelta%isok"
            error stop
        end if
        select case (solvent(1)%name)
        case ("spce")
            solvent(1)%cdelta%filename='input/direct_correlation_functions/water/SPCE/cdelta.in'
            solvent(1)%cd%filename='input/direct_correlation_functions/water/SPCE/cd.in'
        case ("spc")
            solvent(1)%cdelta%filename='input/direct_correlation_functions/water/SPC_Lionel_Daniel/cdelta.in'
            solvent(1)%cd%filename='input/direct_correlation_functions/water/SPC_Lionel_Daniel/cd.in'
        case ("stockmayer")
            solvent(1)%cdelta%filename='input/direct_correlation_functions/stockmayer/cdelta.in'
            solvent(1)%cd%filename='input/direct_correlation_functions/stockmayer/cd.in'
        case ("perso")
            solvent(1)%cdelta%filename='input/cdelta.in'
            solvent(1)%cd%filename='input/cd.in'
        case default
            print*, "In read_cdeltacd, solvent(1)%name=", solvent(1)%name
            print*, "I don't know how to get cdelta(k) and cd(k) for this name"
            error stop
        end select
        nk = n_linesInFile(solvent(1)%cdelta%filename)
        if (nk /= n_linesInFile(solvent(1)%cd%filename)) then
            print*, "In read_cdeltacd, these two files don't have the same number of points"
            print*, solvent(1)%cdelta%filename, solvent(1)%cd%filename
        end if
        allocate (solvent(1)%cd%x(nk), solvent(1)%cd%y(nk), solvent(1)%cdelta%x(nk), solvent(1)%cdelta%y(nk), source=0._dp)

        open (13, file=solvent(1)%cdelta%filename, iostat=ios)
        IF (ios/=0) THEN
            WRITE(*,*)'Cant open file ',solvent(1)%cdelta%filename,' in read_cdeltacd'
            STOP
        END IF
        open( 30, file="./output/cdelta.in")
        DO i = 1, size(solvent(1)%cdelta%x)
            READ (13,*,IOSTAT=ios) solvent(1)%cdelta%x(i), solvent(1)%cdelta%y(i)
            IF (ios/=0) THEN
                WRITE(*,*)'Error while reading ',solvent(1)%cdelta%filename, 'in read_cdeltacd'
                STOP
            END IF
            write(30,*) solvent(1)%cdelta%x(i) , solvent(1)%cdelta%y(i)
        END DO
        close (13)
        close (30)

        OPEN (13, FILE=solvent(1)%cd%filename, IOSTAT=ios)
        IF (ios/=0) THEN
            WRITE(*,*)'Cant open file ',solvent(1)%cd%filename,' in readPolarizationPolarizationCorrelationFunction'
            STOP
        END IF
        open (30, file="./output/cd.in")
        DO i = 1, size(solvent(1)%cd%x)
            READ (13,*,IOSTAT=ios) solvent(1)%cd%x(i), solvent(1)%cd%y(i)
            IF (ios/=0) THEN
                WRITE(*,*)'Error while reading ',solvent(1)%cd%filename, 'in readPolarizationPolarizationCorrelationFunction (c_d)'
                STOP
            END IF
            write(30,*) solvent(1)%cd%x(i) , solvent(1)%cd%y(i)
        END DO
        close (13)
        close (30)

        if (solvent(1)%cd%x(2)-solvent(1)%cd%x(1) /= solvent(1)%cdelta%x(2)-solvent(1)%cdelta%x(1)) then
            print*, "In read_cdeltacd, the distance between two k points is not the same in files"
            print*, solvent(1)%cd%filename, solvent(1)%cdelta%filename
            print*, "It is:", solvent(1)%cd%x(2)-solvent(1)%cd%x(1), solvent(1)%cdelta%x(2)-solvent(1)%cdelta%x(1)
            error stop
        end if

        solvent(1)%cdelta%isok = .true.
        solvent(1)%cd%isok = .true.

    end subroutine read_cdeltacd

end module module_energy_cdeltacd

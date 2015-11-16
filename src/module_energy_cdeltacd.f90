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
        real(dp) :: kT, dv, ksq, k, dk, vexc, kmax, rho, rho0
        complex(dp) :: kP
        real(dp) :: kx(grid%nx/2+1), ky(grid%ny), kz(grid%nz)
        real(dp) :: kxsq(grid%nx/2+1), kysq(grid%ny), kzsq(grid%nz)
        real(dp), intent(out) :: dF(:,:,:,:,:)
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
        if (.not. solvent(1)%cdelta%isok .or. .not. solvent(1)%cd%isok) call init_dcf ("cdeltacd")

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
                    Px(ix,iy,iz) =sum(grid%w(:)*grid%OMx(:)*solvent(1)%density(ix,iy,iz,:))!/rho0
                    Py(ix,iy,iz) =sum(grid%w(:)*grid%OMy(:)*solvent(1)%density(ix,iy,iz,:))!/rho0
                    Pz(ix,iy,iz) =sum(grid%w(:)*grid%OMz(:)*solvent(1)%density(ix,iy,iz,:))!/rho0
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
        if (dk<=epsilon(1._dp)) then
            print*, "In module_energy_cs, dk dans cdelta(k) <= eps machine"
            error stop
        end if
        ikmax =size(solvent(1)%cdelta%x)
        kmax =solvent(1)%cdelta%x(ikmax)

        if (sqrt(maxval(kx)**2+maxval(ky)**2+maxval(kz)**2)>kmax) then
            print*, "In module_energy_cdeltacs"
            print*, "because of our grid",nx,ny,nz,"and deltax:",grid%dx,grid%dy,grid%dz
            print*, "|k| will be as high as:",sqrt(maxval(kx)**2+maxval(ky)**2+maxval(kz)**2)
            print*, "that is higher than what we have tabulated :",kmax
            print*, "instead of this kmax, we use the last value we have in the file :",kmax
            error stop
        end if

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


        rho0=solvent(1)%rho0
        fexc=0._dp
        df = 0._dp
        is=1
        do ix=1,nx
            do iy=1,ny
                do iz=1,nz
                    do io=1,no
                        rho=solvent(1)%density(ix,iy,iz,io)
                        vexc=-kT*grid%w(io)*(grid%OMx(io)*Ex(ix,iy,iz)+grid%OMy(io)*Ey(ix,iy,iz)+grid%OMz(io)*Ez(ix,iy,iz))
                        fexc=fexc+(rho-rho0)*0.5_dp*vexc*dv
                        df(ix,iy,iz,io,is)=df(ix,iy,iz,io,is)+vexc*dv
                    end do
                end do
            end do
        end do

        deallocate (Ex,Ey,Ez)

    end subroutine energy_cdeltacd
end module module_energy_cdeltacd

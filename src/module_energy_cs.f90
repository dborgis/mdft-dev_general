module module_energy_cs
    implicit none
    private
    public :: energy_cs
contains
    subroutine energy_cs (fexc, df)

        use precision_kinds, only: dp
        use module_thermo, only: thermo
        use module_solvent, only: solvent, print_solvent_not_allocated
        use fft, only: fftw3
        use module_grid, only: grid
        use module_dcf, only: init_dcf

        implicit none
        integer :: nx, ny, nz, no, ns, ix, iy, iz, is, io, ik
        real(dp), intent(out) :: Fexc
        real(dp) :: kT, dV, ksqmax, ksq, c_loc, a, k, dk
        real(dp) :: kxsq(grid%nx/2+1), kysq(grid%ny), kzsq(grid%nz)
        real(dp), intent(inout) :: dF(:,:,:,:,:)


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
        do iz=1,nz
            do iy=1,ny
                do ix=1,nx
                    fftw3%in_forward(ix,iy,iz) = sum (solvent(1)%density(ix,iy,iz,:) * grid%w(:)) - solvent(1)%n0
                end do
            end do
        end do

        !
        ! FFT(deltan(x)) => deltan(k)
        !
        call dfftw_execute( fftw3%plan_forward )

        !
        ! Now we should take care of cs(k). Do we already loaded it?
        !
        if ( .not. solvent(1)%cs%isok ) call init_dcf ("cs")

        !
        ! Since the distance between two points is constant in cs
        !
        dk = solvent(1)%cs%x(2) - solvent(1)%cs%x(1)
        if (dk<=epsilon(1._dp)) then
            print*, "In module_energy_cs, dk dans cs(k) <= eps machine"
            error stop
        end if

        !
        ! multiply c_s(k) by deltan(k) and fill in_backward of fft by it
        !
        ksqmax = maxval(solvent(1)%cs%x)**2
        kxsq(1:nx/2+1) = [( grid%kx(ix)**2  ,ix=1,nx/2+1)]
        kysq(1:ny)     = [( grid%ky(iy)**2  ,iy=1,ny)]
        kzsq(1:nz)     = [( grid%kz(iz)**2  ,iz=1,nz)]

        if (maxval(kxsq) + maxval(kysq) + maxval(kzsq) > ksqmax) then
            print*, "Dans module_energy_cs > energy_cs, on a tellement de points k que le plus grand |k| est plus grand"
            print*, "que celui du fichier c(|k|) en input"
            error stop
        end if

        fftw3%in_backward(:,:,:) = cmplx(0._dp,0._dp)
        do iz=1,nz
            do iy=1,ny
                do ix=1,nx/2+1
                    ksq = kxsq(ix) + kysq(iy) + kzsq(iz)
                    k = sqrt(ksq)
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
        Fexc = -kT/2._dp*dv*sum (fftw3%in_forward*fftw3%out_backward)

        do is=1,ns
            do io=1,no
                do iz=1,nz
                    do iy=1,ny
                        do ix=1,nx
                            df(ix,iy,iz,io,is) = df(ix,iy,iz,io,is) -kT*dv*grid%w(io)*fftw3%out_backward(ix,iy,iz)
                        end do
                    end do
                end do
            end do
        end do

    end subroutine energy_cs
end module module_energy_cs

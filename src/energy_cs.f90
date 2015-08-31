module energy
    implicit none
    private
    public :: energy_cs
contains
    subroutine energy_cs (cs, Fexcnn, dFexcnn, exitstatus)

        use precision_kinds , only: i2b, dp
        use system          , only: spacegrid, solvent, thermocond, nb_species
        use minimizer       , only: cg_vect_new, dF_new
        use fft             , only: fftw3, kx, ky, kz
        use dcf             , only: cfile
        use mathematica     , only: splint

        implicit none
        integer :: nx, ny, nz, no, ns, i,j,k,s, icg, io 
        real(dp), intent(out) :: Fexcnn ! what we want to compute in this routine
        real(dp) :: kT, dV, kz2, kz2_ky2, k2, c_loc, psi, fact, Vint, k2max
        integer, intent(out) :: exitstatus
        type(cfile), intent(in) :: cs
        real(dp), intent(out) :: dFexcnn(:,:,:,:,:)
        real(dp) :: n(spacegrid%nx,spacegrid%ny,spacegrid%nz), n_loc

        exitstatus = 1 ! everything is fine
        nx = spacegrid%nx
        ny = spacegrid%ny
        nz = spacegrid%nz
        no = spacegrid%no
        ns = nb_species
        ! allocate( dFexcnn(nx,ny,nz,o,p,nb_species) ,source=0._dp)
        kT = thermocond%kbT
        dV = spacegrid%dV

        if( size(solvent)/=1 ) then
            exitstatus=-1
            return
        else
            s=1
        end if

        n(:,:,:) = 0._dp
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    do io = 1, no
                        n(i,j,k) = n(i,j,k) + cg_vect_new(i,j,k,io,s)**2 * solvent(s)%rho0 * spacegrid%w(io)
                    end do
                end do
            end do
        end do

        fftw3%in_forward = n - solvent(s)%n0
        call dfftw_execute( fftw3%plan_forward )

        k2max=maxval(cs%x)
        do k=1,nz
            kz2 = kz(k)**2
            do j=1,ny
                kz2_ky2=kz2+ky(j)**2
                do i=1,nx/2+1
                    k2 = sqrt(kz2_ky2+kx(i)**2)
                    if( k2<k2max ) then
                        call splint( xa=cs%x, ya=cs%y, y2a=cs%y2, n=size(cs%y), x=k2, y=c_loc)
                    else
                        c_loc=0
                    end if
                    fftw3%in_backward(i,j,k) = fftw3%out_forward(i,j,k) * c_loc
                end do
            end do
        end do
        ! at this point we have gamma(k)

        call dfftw_execute( fftw3%plan_backward )
        fftw3%out_backward = fftw3%out_backward /real(nx*ny*nz,dp) ! gamma(x)

        ! excess free energy
        Fexcnn = -kT/2 *sum( (n-solvent(1)%n0 )*fftw3%out_backward )*dV

        ! gradient
        do concurrent( i=1:nx, j=1:ny, k=1:nz, io=1:no, s=1:nb_species )
            dfexcnn(i,j,k,io,s) = 2*cg_vect_new(i,j,k,io,s) * spacegrid%w(io) * fftw3%out_backward(i,j,k) *(-kT*dV*solvent(s)%rho0)
        end do

    end subroutine energy_cs
end module energy

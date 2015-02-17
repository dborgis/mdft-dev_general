subroutine energy_cs (cs, Fexcnn, dFexcnn, exitstatus)

  USE ISO_FORTRAN_ENV
  use precision_kinds , only: i2b, dp
  use system          , only: spacegrid, solvent, thermocond
  use quadrature      , only: anggrid, molrotgrid
  use minimizer       , only: from_cgvect_get_rho, from_rho_get_n, cg_vect, deallocate_solvent_rho
  use fft             , only: fftw3, kx, ky, kz
  use dcf             , only: cfile
  use mathematica     , only: splint

  implicit none
  integer :: nfft1, nfft2, nfft3, i,j,k,s,icg,o,p
  real(dp), intent(out) :: Fexcnn ! what we want to compute in this routine
  real(dp) :: kT, dV, kz2, kz2_ky2, k2, c_loc, psi, fact, Vint, k2max
  integer, intent(inout) :: exitstatus
  type(cfile), intent(in) :: cs
  real(dp), intent(out) :: dFexcnn(size(cg_vect)) ! gradient

  exitstatus=1 ! everything is fine
  nfft1 = spacegrid%n_nodes(1)
  nfft2 = spacegrid%n_nodes(2)
  nfft3 = spacegrid%n_nodes(3)
  kT = thermocond%kbT
  dV = spacegrid%dV

  call from_cgvect_get_rho ! returns
  call from_rho_get_n      ! returns solvent(:)%n(nfft1,nfft2,nfft3)
  call deallocate_solvent_rho

  if( size(solvent)/=1 ) then
    exitstatus=-1
  else
    s=1
  end if

  fftw3%in_forward = solvent(s)%n - solvent(s)%n0
  call dfftw_execute( fftw3%plan_forward )

  k2max=-1
  do k=1,nfft3
    kz2 = kz(k)**2
    do j=1,nfft2
      kz2_ky2=kz2+ky(j)**2
      do i=1,nfft1/2+1
        k2 = sqrt(kz2_ky2+kx(i)**2)
        call splint( xa=cs%x, ya=cs%y, y2a=cs%y2, n=size(cs%y), x=k2, y=c_loc)
        fftw3%in_backward(i,j,k) = fftw3%out_forward(i,j,k) * c_loc
      end do
    end do
  end do
  ! at this point we have gamma(k)

  call dfftw_execute( fftw3%plan_backward )
  fftw3%out_backward = fftw3%out_backward /(nfft1*nfft2*nfft3)

  ! excess free energy
  Fexcnn = -kT/2 *sum( (solvent(1)%n-solvent(1)%n0 )*fftw3%out_backward )*dV

  ! gradient
  icg=0
  do s=1,size(solvent)
    fact= -kT*dV*solvent(s)%rho0
    do i=1,nfft1
      do j=1,nfft2
        do k=1,nfft3
          Vint=fftw3%out_backward(i,j,k) *fact
          do o=1,anggrid%n_angles
            do p=1,molrotgrid%n_angles
              icg=icg +1
              psi=cg_vect(icg)
              dFexcnn(icg)= 2*psi*anggrid%weight(o)*molrotgrid%weight(p)*Vint
            end do
          end do
        end do
      end do
    end do
  end do

end subroutine energy_cs

module energy
    implicit none
contains
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine energy_cs (cs, Fexcnn, dFexcnn, exitstatus)

  use precision_kinds , only: i2b, dp
  use system          , only: spacegrid, solvent, thermocond, nb_species
  use quadrature      , only: anggrid, molrotgrid
  use minimizer       , only: from_cgvect_get_rho, from_rho_get_n, cg_vect_new, deallocate_solvent_rho, dF_new
  use fft             , only: fftw3, kx, ky, kz
  use dcf             , only: cfile
  use mathematica     , only: splint

  implicit none
  integer :: nfft1, nfft2, nfft3, i,j,k,s,icg,o,p
  real(dp), intent(out) :: Fexcnn ! what we want to compute in this routine
  real(dp) :: kT, dV, kz2, kz2_ky2, k2, c_loc, psi, fact, Vint, k2max
  integer, intent(out) :: exitstatus
  type(cfile), intent(in) :: cs
  real(dp), intent(out) :: dFexcnn(:,:,:,:,:,:)
  real(dp) :: n(spacegrid%n_nodes(1),spacegrid%n_nodes(2),spacegrid%n_nodes(3)), n_loc

  exitstatus=1 ! everything is fine
  nfft1 = spacegrid%n_nodes(1)
  nfft2 = spacegrid%n_nodes(2)
  nfft3 = spacegrid%n_nodes(3)
  o = anggrid%n_angles
  p = molrotgrid%n_angles
 ! allocate( dFexcnn(nfft1,nfft2,nfft3,o,p,nb_species) ,source=0._dp)
  kT = thermocond%kbT
  dV = spacegrid%dV

  if( size(solvent)/=1 ) then
    exitstatus=-1
    return
  else
    s=1
  end if

  do concurrent( i=1:nfft1, j=1:nfft2, k=1:nfft3 )
    n_loc = 0
    do concurrent( o=1:anggrid%n_angles, p=1:molrotgrid%n_angles )
      n_loc = n_loc + cg_vect_new(i,j,k,o,p,s)**2*solvent(s)%rho0*molrotgrid%weight(p)*anggrid%weight(o)
    end do
    n(i,j,k) = n_loc
  end do

  fftw3%in_forward = n - solvent(s)%n0
  call dfftw_execute( fftw3%plan_forward )

  k2max=maxval(cs%x)
  do k=1,nfft3
    kz2 = kz(k)**2
    do j=1,nfft2
      kz2_ky2=kz2+ky(j)**2
      do i=1,nfft1/2+1
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
  fftw3%out_backward = fftw3%out_backward /(nfft1*nfft2*nfft3) ! gamma(x)

  ! excess free energy
  Fexcnn = -kT/2 *sum( (n-solvent(1)%n0 )*fftw3%out_backward )*dV

  ! gradient
  do concurrent( i=1:nfft1, j=1:nfft2, k=1:nfft3, o=1:anggrid%n_angles, p=1:molrotgrid%n_angles, s=1:nb_species )
    dfexcnn(i,j,k,o,p,s) = 2 *cg_vect_new(i,j,k,o,p,s) *anggrid%weight(o) *molrotgrid%weight(p) *fftw3%out_backward(i,j,k) *(-kT*dV*solvent(s)%rho0)
  end do

end subroutine energy_cs
end module energy

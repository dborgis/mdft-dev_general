!This routine compute : -The solvent molecular charge density, which can be used into Vcoul_from_solvent_charge_density.f90 to
!evaluate the electrostatic potential.
!                       -The solvent molecular polarization (from Ranieriet al : J. Chem. Phys. 98 (11) 1993) that can be used
!into energy_polarization_..._.f90 to compute the (multipolar) polarization Free energy.

subroutine chargeDensityAndMolecularPolarizationOfASolventMoleculeAtOrigin (Rotxx,Rotxy,Rotxz,Rotyx,Rotyy,Rotyz,Rotzx,Rotzy,Rotzz)

  use precision_kinds     ,ONLY: i2b, dp
  use constants           ,ONLY: iC=>i_complex, zeroC
  use system              ,ONLY: solvent, spaceGrid
  use quadrature          ,ONLY: angGrid, molRotGrid
  use fft                 ,ONLY: kx, ky, kz, k2
  use mathematica         ,ONLY: factorial

  implicit none

  real(dp), dimension(angGrid%n_angles,molRotGrid%n_angles), intent(in) :: Rotxx,Rotxy,Rotxz,Rotyx,Rotyy,Rotyz,Rotzx,Rotzy,Rotzz
  integer :: nfft1, nfft2, nfft3
  integer(i2b) :: i, j, k, o, p, n, s
  real(dp)     :: r(3), kr, kvec(3)
  complex(dp)  :: fac, X
  type :: smoother_type
    real(dp) :: radius = 0.5_dp ! dramaticaly important
    real(dp) :: factor
  end type smoother_type
  type (smoother_type) :: smoother

  nfft1 = spacegrid%n_nodes(1)
  nfft2 = spacegrid%n_nodes(2)
  nfft3 = spacegrid%n_nodes(3)

  ! sigma_k is the Fourier transformed charge density of a single solvent molecule in the reference frame defined by solvent.in
  ! molec_polar_k is the Fourier transformed molecular polarization
  do concurrent( s=1:size(solvent) , sum(abs(solvent(s)%site%q))>epsilon(1._dp)) ! mask elimitates solvent molecules without point charges
    allocate( solvent(s)%sigma_k(nfft1/2+1, nfft2, nfft3, angGrid%n_angles, molRotGrid%n_angles), SOURCE=zeroC)
    allocate( solvent(s)%molec_polar_k(3,nfft1/2+1, nfft2, nfft3, angGrid%n_angles, molRotGrid%n_angles), SOURCE=zeroC)
  end do

  do concurrent ( i=1:nfft1/2+1, j=1:nfft2, k=1:nfft3, s=1:size(solvent)  , sum(abs(solvent(s)%site%q))>epsilon(1._dp)) ! mask elimitates solvent molecules without point charges)

    kvec = [ kx(i), ky(j), kz(k) ]
    smoother%factor =  exp(-smoother%radius**2 *k2(i,j,k)/2)

    do concurrent ( o=1:angGrid%n_angles, p=1:molRotGrid%n_angles, &
            n=1:SIZE(solvent(s)%site), (abs(solvent(s)%site(n)%q)>epsilon(1.0_dp)) )

      r(1) = dot_product(   [Rotxx(o,p),Rotxy(o,p),Rotxz(o,p)]  ,  solvent(s)%site(n)%r  )
      r(2) = dot_product(   [Rotyx(o,p),Rotyy(o,p),Rotyz(o,p)]  ,  solvent(s)%site(n)%r  )
      r(3) = dot_product(   [Rotzx(o,p),Rotzy(o,p),Rotzz(o,p)]  ,  solvent(s)%site(n)%r  )
      kr = dot_product( kvec, r )
      X = -iC*kr
      solvent(s)%sigma_k(i,j,k,o,p) = solvent(s)%sigma_k(i,j,k,o,p) + solvent(s)%site(n)%q *exp(X) *smoother%factor ! exact
      ! solvent(s)%sigma_k(i,j,k,o,p) = solvent(s)%sigma_k(i,j,k,o,p) +solvent(s)%site(n)%q* sum([(X**i/factorial(i), i=0,4)])&
            ! * smoother%factor ! Series expansion of exp(x) at 0 => multipole expansion of Vcoul(x). i=4 :: hexadecapole (16)
      if ( abs(kr)<=epsilon(1.0_dp) ) then
        solvent(s)%molec_polar_k(:,i,j,k,o,p) = solvent(s)%molec_polar_k(:,i,j,k,o,p) + solvent(s)%site(n)%q *r
      else
        fac = -iC*(exp(iC*kr)-1._dp)/kr *smoother%factor
        solvent(s)%molec_polar_k(:,i,j,k,o,p) = solvent(s)%molec_polar_k(:,i,j,k,o,p) + fac*solvent(s)%site(n)%q *r
      end if

    end do
  end do

  call removeTraceOfMolecularPolarization

contains

  subroutine removeTraceOfMolecularPolarization
    use precision_kinds, only: i2b
    use quadrature, only: molRotGrid, angGrid
    implicit none
    integer(i2b) :: i,j,k,s,d
    do concurrent (i=1:nfft1/2+1, j=1:nfft2, k=1:nfft3, s=1:size(solvent), d=1:3)
        solvent(s)%molec_polar_k(d,i,j,k,:,:) = solvent(s)%molec_polar_k(d,i,j,k,:,:)  &
                      -sum( solvent(s)%molec_polar_k(d,i,j,k,:,:) ) /real(angGrid%n_angles*molRotGrid%n_angles,dp)
    end do
  end subroutine removeTraceOfMolecularPolarization

END SUBROUTINE chargeDensityAndMolecularPolarizationOfASolventMoleculeAtOrigin

SUBROUTINE energy_fmt (Ffmt)

  USE precision_kinds  ,ONLY: dp , i2b
  USE system           ,ONLY: thermocond, nb_species, mole_fraction, spaceGrid, solvent
  USE quadrature       ,ONLY: molRotSymOrder, angGrid, molRotGrid
  USE minimizer        ,ONLY: cg_vect_new, dF_new
  USE constants        ,ONLY: pi, fourPi, twopi, zeroC, zerodp, onedp, epsdp
  USE fft              ,ONLY: fftw3
  USE input            ,ONLY: input_char
  USE hardspheres      ,ONLY: hs, weight_functions

  IMPLICIT NONE

  real(dp), intent(out) :: Ffmt ! Internal part of free energy
  integer(i2b) :: icg,i,j,k,o,p,s,nfft1,nfft2,nfft3,wdl,wdu
  real(dp) :: local_density, psi, dV, nb_molecules(size(solvent)), time0, time1, kT
  real(dp)   , ALLOCATABLE, DIMENSION(:,:,:,:) :: rho ! density per angle (recall : rho_0 = n_0 / 4pi ) ! x y z nb_species
  real(dp)   , ALLOCATABLE , DIMENSION(:,:,:,:) :: wd ! weighted density at node i,j,k, for index 0:3
  real(dp)   , ALLOCATABLE , DIMENSION(:,:,:,:) :: dFHS
  complex(dp), ALLOCATABLE , DIMENSION(:,:,:,:) :: dFHS_k
  real(dp)   , ALLOCATABLE , DIMENSION(:,:,:,:) :: dFex ! gradient in real space
  character(2) :: hs_functional ! hard sphere functional = PY for Percus-Yevick or CS for Carnahan-Starling
  real(dp), PARAMETER :: inv8pi  = 1.0_dp/(  8.0_dp * pi)
  real(dp), PARAMETER :: inv12pi = 1.0_dp/( 12.0_dp * pi)
  real(dp), PARAMETER :: inv18pi = 1.0_dp/( 18.0_dp * pi)
  real(dp), PARAMETER :: inv24pi = 1.0_dp/( 24.0_dp * pi)
  real(dp), PARAMETER :: inv36pi = 1.0_dp/( 36.0_dp * pi)

  IF (nb_species/=1) STOP "I stop because FMT calculations are only possible with one solvent species."

  CALL CPU_TIME ( time0 ) ! init timer

  nfft1 = spaceGrid%n_nodes(1)
  nfft2 = spaceGrid%n_nodes(2)
  nfft3 = spaceGrid%n_nodes(3)
  dV = spaceGrid%dv
  kT = thermocond%kbT

  ! in case of empty supercell
  if( all(abs(cg_vect_new)<=epsdp) ) then
    Ffmt = -sum( hs%Fexc0 * mole_fraction )
    return
  end if

  allocate( rho (nfft1,nfft2,nfft3,nb_species) ,SOURCE=0._dp)
  icg = 0
  do s=1,size(solvent)
    do i=1,nfft1
      do j=1,nfft2
        do k=1,nfft3
          local_density=0
          do o=1,angGrid%n_angles
            do p=1,molRotGrid%n_angles
              icg = icg + 1
              local_density = local_density &
                + cg_vect_new(i,j,k,o,p,s)**2 *molRotGrid%weight(p) *angGrid%weight(o)
            end do
          end do
          ! correct by *(8*pi²/n)**-1 as the integral over all orientations o and psi is 4pi and 2pi/n
          ! at the same time integrate rho in order to count the total number of implicit molecules.
          rho(i,j,k,s) = local_density*molRotSymOrder/(fourpi*twopi) ! n/n0 at this stage
        end do
      end do
    end do
  end do

  ! total number of molecules of each species
  do concurrent ( s=1:nb_species )
    nb_molecules(s) = sum(rho(:,:,:,s))*dV *solvent(s)%n0 *mole_fraction(s)
  end do

  ! compute the weight functions
  call weight_functions

  ! weighted densities
  wdl=0 ! 0:3 for Kierliek Rosinberg
  wdu=3
  if( lbound(hs(1)%w_k,4) /= 0 ) then
    stop "problem with lbound l95 of energy_hard_sphere_fmt.f90"
  else if( ubound(hs(1)%w_k,4) /= 3 ) then
    stop "problem with ubound l98 of energy_hard_sphere_fmt.f90"
  else
    allocate( wd (nfft1,nfft2,nfft3, 0:3) ,SOURCE=0._dp)
  end if
  do s=1,nb_species
    fftw3%in_forward = rho(:,:,:,s)
    call dfftw_execute ( fftw3%plan_forward ) ! fourier transform the density and put it into fftw3%out_forward
    do i=0,3
      fftw3%in_backward = fftw3%out_forward * solvent(s)%n0 * mole_fraction(s) * cmplx( hs(s)%w_k(:,:,:,i) ,0)! rho_k(s) * w_k(i)
      call dfftw_execute( fftw3%plan_backward )
      wd(:,:,:,i) = wd(:,:,:,i) + fftw3%out_backward /(nfft1*nfft2*nfft3)
    end do
  end do
  deallocate(rho)

  ! check if the hard sphere functional is Percus-Yevick or Carnahan-Starling
  ! Get the free energy functional that should be used. For now Percus Yevick and Carnahan Starling only. May be expanded.
  hs_functional=input_char('hs_functional')
  if( hs_functional(1:2)=='PY') then
    if( any(wd(:,:,:,3)>=1) ) stop "problem in log(1-w3) in energy_fmt"
    Ffmt =  -sum(wd(:,:,:,0)*log(1-wd(:,:,:,3))) &
            +sum(wd(:,:,:,1)*wd(:,:,:,2)/(1-wd(:,:,:,3))) &
            +sum(wd(:,:,:,2)**3/(1-wd(:,:,:,3))**2)/(24*pi)
  else if( hs_functional(1:2)=='CS') then
    if( any(wd(:,:,:,3)>=1) ) stop "problem in log(1-w3) in energy_fmt"
    if( any(wd(:,:,:,3)<=epsdp) ) stop "I stop before dividing by 0 in energy_fmt"
    Ffmt = sum(wd(:,:,:,1)*wd(:,:,:,2)/(1-wd(:,:,:,3))) &
          +sum( (inv36pi*wd(:,:,:,2)**3 / wd(:,:,:,3)**2 -wd(:,:,:,0)) * LOG((1-wd(:,:,:,3)))  ) &
          +sum( inv36pi*wd(:,:,:,2)**3/(wd(:,:,:,3)*(1-wd(:,:,:,3))**2)  )
  else
    stop "hs_functional not initialized in energy_hard_spheres_fmt.f90:121"
  end if

  Ffmt = Ffmt*kT*dV -sum(hs%excchempot*nb_molecules) -sum(hs%Fexc0*mole_fraction)

  ! gradients
  ! dFHS_i and weighted_density_j are arrays of dimension (nfft1,nfft2,nfft3)
  ! excess free energy per atom derived wrt weighted density 0 (eq. A5 in Sears2003)
  allocate ( dFHS(nfft1,nfft2,nfft3,wdl:wdu) ,SOURCE=0._dp)

  ! Perkus Yevick
  if( hs_functional(1:2)=='PY') then
    IF( ANY((1-wd(:,:,:,3))<=epsdp)) stop "I found some value in 1-weightedensity3 that cannot go in log"
    dFHS(:,:,:,0) = -log(1-wd(:,:,:,3))
    dFHS(:,:,:,1) = wd(:,:,:,2)/(1-wd(:,:,:,3))
    dFHS(:,:,:,2) = wd(:,:,:,1)/(1-wd(:,:,:,3)) +inv8pi*dFHS(:,:,:,1)**2
    dFHS(:,:,:,3) = (wd(:,:,:,0)+wd(:,:,:,1)*dFHS(:,:,:,1))/(1-wd(:,:,:,3)) +inv12pi*dFHS(:,:,:,1)**3
  else if( hs_functional(1:2)=='CS') then
    IF( ANY((1-wd(:,:,:,3))<=epsdp)) STOP "I found some value of 1-weightedensity3 that cannot go in log"
    IF( ANY(ABS(wd(:,:,:,3))<=epsdp)) STOP "I found some value in w3 that would lead to divide by 0"

    dFHS(:,:,:,0) = - log ( (1-wd(:,:,:,3)) )

    dFHS(:,:,:,1) = wd(:,:,:,2) / (1-wd(:,:,:,3))

    dFHS(:,:,:,2) = - inv12pi * ( wd(:,:,:,2) / wd(:,:,:,3) ) ** 2 * dFHS(:,:,:,0) &
    + wd(:,:,:,1) / (1-wd(:,:,:,3)) + inv12pi * dFHS(:,:,:,1) ** 2 / wd(:,:,:,3)

    dFHS(:,:,:,3) = inv18pi * dFHS(:,:,:,0) * ( wd(:,:,:,2) / wd(:,:,:,3) ) ** 3 &
    - ( inv36pi * wd(:,:,:,2) ** 3 / wd(:,:,:,3) ** 2 - wd(:,:,:,0) ) /&
    (1-wd(:,:,:,3)) &
    + wd(:,:,:,1) * dFHS(:,:,:,1) / (1-wd(:,:,:,3)) + inv36pi * &
    wd(:,:,:,2) ** 3 / wd(:,:,:,3) ** 2 * &
    ( 3.0_dp * wd(:,:,:,3) - 1.0_dp ) / (1-wd(:,:,:,3)) ** 3
  else
    stop "hs_functional not initialized in energy_hard_spheres_fmt.f90:121"
  end if
  deallocate(wd) ! deallocate weighted_densities

  ! compute gradients in k space
  allocate( dFHS_k(nfft1/2+1,nfft2,nfft3,0:3) ,SOURCE=zeroC)
  ! FFT dFHS for computing convolution
  do i= wdl,wdu
    fftw3%in_forward = dFHS(:,:,:,i)
    CALL dfftw_execute ( fftw3%plan_forward )
    dFHS_k(:,:,:,i) = fftw3%out_forward
  end do
  deallocate( dFHS )

  ! compute final gradient in k-space
  allocate( dFex (nfft1,nfft2,nfft3,nb_species) ,source=zerodp)
  do s=1,nb_species
    fftw3%in_backward= dFHS_k(:,:,:,0)*hs(s)%w_k(:,:,:,0) &
                      +dFHS_k(:,:,:,1)*hs(s)%w_k(:,:,:,1) &
                      +dFHS_k(:,:,:,2)*hs(s)%w_k(:,:,:,2) &
                      +dFHS_k(:,:,:,3)*hs(s)%w_k(:,:,:,3)
    call dfftw_execute ( fftw3%plan_backward )
    dFex(:,:,:,s) = fftw3%out_backward / (nfft1*nfft2*nfft3)
  end do

  deallocate( dFHS_k )
  do s=1,nb_species
    deallocate( hs(s)%w_k )
  end do

  ! transfer in rank 1 vector dF
  icg=0
  do s=1,nb_species
    do i=1,nfft1
      do j=1,nfft2
        do k=1,nfft3
          do o=1,angGrid%n_angles
            do p=1,molRotGrid%n_angles
              icg=icg+1
              psi=cg_vect_new(i,j,k,o,p,s)
              dF_new(i,j,k,o,p,s) = dF_new(i,j,k,o,p,s)&
                  + 2*psi*solvent(s)%rho0*dV*( kT*dFex(i,j,k,s)-hs(s)%excchempot )*angGrid%weight(o)*molRotGrid%weight(p)
            end do
          end do
        end do
      end do
    end do
  end do
  deallocate (dFex)

  call cpu_time( time1 )

end subroutine energy_fmt

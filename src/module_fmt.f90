module module_fmt

contains

subroutine  energy_fmt (Ffmt,df)

  use precision_kinds  ,ONLY: dp, i2b
  use module_thermo, only: thermo
  use module_solvent, only: solvent
  use module_fft              ,ONLY: fftw3
  use module_input            ,ONLY: getinput
  use hardspheres      ,ONLY: hs, weight_functions
  use module_grid, only: grid,dq, norm_k
  use constants
  use module_dcf, only :  c_s_hs

  IMPLICIT NONE

  real(dp), intent(out) :: Ffmt ! Internal part of free energy
  real(dp), intent(inout), contiguous, optional :: df(:,:,:,:,:)
  integer(i2b) :: icg,i,j,k,s,nx,ny,nz,wdl,wdu,io,iq
  real(dp) :: local_density, psi, dV, nb_molecules(size(solvent)), time0, time1, kT,q
  real(dp)   , ALLOCATABLE, DIMENSION(:,:,:,:) :: rho ! density per angle (recall : rho_0 = n_0 / 4pi ) ! x y z solvent(1)%nspec
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

  IF (solvent(1)%nspec/=1) STOP "I stop because FMT calculations are only possible with one solvent species."
  CALL CPU_TIME ( time0 ) ! init timer

  nx = grid%nx
  ny = grid%ny
  nz = grid%nz
  dV = grid%dv
  kT = thermo%kBT 

  ! in case of empty supercell
  if( all(abs(solvent(1)%xi)<=epsdp) ) then
    Ffmt = -sum( hs%Fexc0 * solvent(1)%mole_fraction )
    return
  end if

  allocate( rho (nx,ny,nz,solvent(1)%nspec) ,SOURCE=0._dp)
  do s=1,size(solvent)
    do i=1,nx
      do j=1,ny
        do k=1,nz
          local_density=0
          do io=1,grid%no
            local_density = local_density + solvent(s)%xi(io,i,j,k)**2 * grid%w(io)
          end do
          ! correct by *(8*pi²/n)**-1 as the integral over all orientations o and psi is 4pi and 2pi/n
          ! at the same time integrate rho in order to count the total number of implicit molecules.
          rho(i,j,k,s) = local_density*grid%molRotSymOrder/(fourpi*twopi) ! n/n0 at this stage
        end do
      end do
    end do
  end do

  ! total number of molecules of each species
  do concurrent ( s=1:solvent(1)%nspec )
    nb_molecules(s) = sum(rho(:,:,:,s))*dV *solvent(s)%n0 *solvent(s)%mole_fraction
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
    allocate( wd (nx,ny,nz, 0:3) ,SOURCE=0._dp)
  end if
  do s=1,solvent(1)%nspec
    fftw3%in_forward = rho(:,:,:,s)
    call dfftw_execute ( fftw3%plan_forward ) ! fourier transform the density and put it into fftw3%out_forward
    do i=0,3
      fftw3%in_backward = fftw3%out_forward * solvent(s)%n0 * solvent(s)%mole_fraction * cmplx( hs(s)%w_k(:,:,:,i) ,0)! rho_k(s) * w_k(i)
      call dfftw_execute( fftw3%plan_backward )
      wd(:,:,:,i) = wd(:,:,:,i) + fftw3%out_backward /(nx*ny*nz)
    end do
  end do
  


  ! check if the hard sphere functional is Percus-Yevick or Carnahan-Starling
  ! Get the free energy functional that should be used. For now Percus Yevick and Carnahan Starling only. May be expanded.
  hs_functional=getinput%char('hs_functional')
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


  Ffmt = Ffmt*kT*dV -sum(hs%excchempot*nb_molecules) -sum(hs%Fexc0*solvent(1)%mole_fraction)
  

  ! gradients
  ! dFHS_i and weighted_density_j are arrays of dimension (nx,ny,nz)
  ! excess free energy per atom derived wrt weighted density 0 (eq. A5 in Sears2003)
  allocate ( dFHS(nx,ny,nz,wdl:wdu) ,SOURCE=0._dp)

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
  allocate( dFHS_k(nx/2+1,ny,nz,0:3) ,SOURCE=zeroC)
  ! FFT dFHS for computing convolution
  do i= wdl,wdu
    fftw3%in_forward = dFHS(:,:,:,i)
    CALL dfftw_execute ( fftw3%plan_forward )
    dFHS_k(:,:,:,i) = fftw3%out_forward
  end do
  deallocate( dFHS )

  ! compute final gradient in k-space
  allocate( dFex (nx,ny,nz,solvent(1)%nspec) ,source=zerodp)
  do s=1,solvent(1)%nspec
    fftw3%in_backward= dFHS_k(:,:,:,0)*hs(s)%w_k(:,:,:,0) &
                      +dFHS_k(:,:,:,1)*hs(s)%w_k(:,:,:,1) &
                      +dFHS_k(:,:,:,2)*hs(s)%w_k(:,:,:,2) &
                      +dFHS_k(:,:,:,3)*hs(s)%w_k(:,:,:,3)
    call dfftw_execute ( fftw3%plan_backward )
    dFex(:,:,:,s) = fftw3%out_backward / (nx*ny*nz)
  end do

  deallocate( dFHS_k )
  do s=1,solvent(1)%nspec
    deallocate( hs(s)%w_k )
  end do

  ! transfer in rank 1 vector dF
  icg=0
  do s=1,solvent(1)%nspec
    do i=1,nx
      do j=1,ny
        do k=1,nz
          do io=1,grid%no
              icg=icg+1
              psi=solvent(1)%xi(io,i,j,k)
              df(io,i,j,k,s) = df(io,i,j,k,s) &
                  + 2*psi*solvent(s)%rho0*dV*( kT*dFex(i,j,k,s)-hs(s)%excchempot )*grid%w(io)
          end do
        end do
      end do
    end do
  end do
  deallocate (dFex)


  call cpu_time( time1 )


!!!Now remove the c_hs contribution
  !rho is bearing deltaN now
  rho=rho*solvent(1)%n0*solvent(1)%mole_fraction-solvent(1)%n0
  fftw3%in_forward = rho(:,:,:,1)
  call dfftw_execute( fftw3%plan_forward )
  fftw3%in_backward = (0._dp,0._dp)
  do k=1,nx
      do j=1,ny
          do i=1,nz/2+1
              q =  norm_k(i,j,k) 
              iq = int(q/dq) +1
              fftw3%in_backward(i,j,k) = fftw3%out_forward(i,j,k) * c_s_hs%y(iq)
          end do
      end do
  end do
  ! at this point we have gamma(k) in fftw3%in_backward

  call dfftw_execute( fftw3%plan_backward ) ! this fill fftw3%out_backward with gamma(x)
  fftw3%out_backward = fftw3%out_backward  /real(nx*ny*nz,dp) ! gamma(x)  note: the normalization factor /real(nx,ny,nz) comes from the way FFTW3 does NOT normalize its FFT
  
  !This is a HS bridge term so watch the sign
  Ffmt=Ffmt+kT/2._dp*dv*sum (rho(:,:,:,1)*fftw3%out_backward)

  do s=1,size(solvent)
    do k=1,nz
      do j=1,ny
        do i=1,nx
          do io=1,grid%no
            df(io,i,j,k,s) = df(io,i,j,k,s) &
                +kT*grid%w(io)*fftw3%out_backward(i,j,k)*2._dp*solvent(s)%rho0*solvent(s)%xi(i,j,k,s)
          end do
        end do
      end do
    end do
  end do

  deallocate(rho)
end subroutine energy_fmt



end module module_fmt

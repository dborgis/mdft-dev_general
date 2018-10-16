module hardspheres
  use precision_kinds  ,only: dp
  implicit none
  private
  type type_hs
    real(dp) :: r          ! radius
    real(dp) :: excchempot ! excess chemical potential
    real(dp) :: fexc0
    real(dp) :: pf         ! packing fraction
    real(dp), allocatable :: w_k(:,:,:,:) ! weight function in Fourier space (i,j,k,0:3)
  end type type_hs
  type(type_hs), public, allocatable :: hs(:) ! one per constituant
  public :: packing_fraction, compute_hard_spheres_parameters, weight_functions, read_hard_sphere_radius

contains

!===================================================================================================================================
! packing fraction : eta = 4/3 * pi * R^3 * number density of the constituant ( /= total solvent density)
pure function packing_fraction (numberdensity, HSradius)
  use precision_kinds ,only: dp
  use constants       ,only: pi
  implicit none
  real(dp) :: packing_fraction
  real(dp), parameter  :: fourthird=4./3.
  real(dp), intent(in) :: numberdensity ! density of the constituant (/= total solvant density)
  real(dp), intent(in) :: HSradius ! radius of the constituant
  packing_fraction = fourthird*pi*numberdensity*HSradius**3
end function packing_fraction

!===================================================================================================================================
!Compute the weight functions defined by Kierlik and Rosinberg, Phys. Rev. A 42, 3382-3387 (1990)
! The four weight functions are here defined in k-space. They are known analyticaly and only depend on the so called fundamental
! measures of the hard spheres.
! w_3i(k=0) = V_i the volume of constituant i
! w_2i(k=0) = S_i the surface area of constituant i
! w_1i(k=0) = R_i the radius of constituant i
! w_0i(k=0) = 1
! They are scalar numbers in opposition to scalar and vector weight functions by Rosenfeld in its seminal Phys. Rev. Lett.
! introducing the fundamental measure theory (FMT)
subroutine weight_functions
  use precision_kinds ,only: dp, i2b
  use constants       ,only: fourpi, zeroC, zerodp
  use module_solvent, only: solvent
  use module_grid, only: grid, norm_k
  implicit none
  real(dp)     :: kR, sinkR, coskR, norm_k_loc, R, k
  integer(i2b) :: l,m,n,s,nfft1,nfft2,nfft3
  nfft1=grid%nx
  nfft2=grid%ny
  nfft3=grid%nz
  do concurrent( s=1:solvent(1)%nspec ,.not.allocated(hs(1)%w_k) )
    allocate( hs(s)%w_k(nfft1/2+1,nfft2,nfft3,0:3) ,source=zerodp)
  end do
  print*, hs(:)%R
  ! Weight functions of a fluid of hard spheres are known analyticaly.
  do concurrent( s=1:solvent(1)%nspec, l=1:nfft1/2+1, m=1:nfft2, n=1:nfft3 )
    R=hs(s)%R
    k=norm_k(l,m,n)
    if( abs(k)>epsilon(1._dp) ) then ! k /= 0.0_dp
      kR=k*R
      sinkR=sin(kR)
      coskR=cos(kR)
      hs(s)%w_k(l,m,n,3) = fourpi*(sinkR-kR*coskR)/(k**3)
      hs(s)%w_k(l,m,n,2) = fourpi*R*sinkR/k
      hs(s)%w_k(l,m,n,1) = (sinkR+kR*coskR)/(2*k)
      hs(s)%w_k(l,m,n,0) = coskR+kR*sinkR/2
    else ! for k=0
      hs(s)%w_k(l,m,n,3) = fourpi/3*R**3 ! volume
      hs(s)%w_k(l,m,n,2) = fourpi  *R**2 ! surface area
      hs(s)%w_k(l,m,n,1) = R             ! radius
      hs(s)%w_k(l,m,n,0) = 1             ! unity
    end if
  end do
end subroutine weight_functions

!===================================================================================================================================
! Read, allocate, compute the hard sphere parameters
! read hard sphere radius
! read the mole fraction of each hard sphere
! compute their weight functions (using their fundamental measures)
! compute their packing fraction
! read the excess functional
! compute accordingly the chemical potential and the reference bulk density
subroutine compute_hard_spheres_parameters
  use module_solvent, only: solvent
  use module_input           ,only: getinput
  implicit none
  character(4) :: hs_functional
  if(.not.allocated(hs)) allocate(hs(solvent(1)%nspec))
  hs_functional=getinput%char('hs_functional')
  if(hs_functional(1:2)/='PY' .and. hs_functional(1:2)/='CS') then
    stop "functional not implemented. See module_hardspheres:97"
  end if
  call read_hard_sphere_radius ! read hard sphere radius and allocate if necessary
  call compute_packing_fractions
  ! compute excess chemical potential and grand potential at reference bulk density
  call excess_chemical_potential_and_reference_bulk_grand_potential
end subroutine compute_hard_spheres_parameters

!===============================================================================================================================
! Subroutine that deduct the packing fraction of all reference bulk fluids of the consitutuants of the mixture.
! We want the reference bulk fluids to have physical meaning, that's why we pay attention to them not to be greater than 0.74, where
! they would be unphysical. 0.74 is the maximum packing of a solid crystal.
subroutine compute_packing_fractions
  use precision_kinds ,only: dp,i2b
  use module_solvent, only: solvent
  use module_input           ,only: verbose
  implicit none
  integer(i2b) :: s
  do s=1,solvent(1)%nspec
    hs(s)%pf = packing_fraction( numberdensity=solvent(s)%n0, HSradius=hs(s)%R )
    if(verbose) write(*,'(A,i1,A,f5.3)') 'Packing fraction of species ',s,' is ',hs(s)%pf
    ! compute homogeneous fluid reference with Perkus Yevick
    ! It is important to keep in mind it is the packing fraction of the REFERENCE fluid(s), not a partial packing fraction of our mixture.
    ! although the Percus-Yevick equation shows no singularities for eta < 1 , the region beyond eta = pi / (3 sqrt(2) ) = 0.74 is unphysical, since the fluid then has a packing density greater than that of a closed packed solid.
    if( hs(s)%pf >= 0.492 .and. hs(s)%pf <0.74 ) then
      print*,"packing fraction of species",s,"is",hs(s)%pf,"which is higher than freezing packing fraction of 0.492"
    else if( hs(s)%pf >= 0.74_dp ) then
      print*,'packing fraction of species',s,"is",hs(s)%pf,"which is higher than metastable fluid-solid pack. frac. of 0.72"
      stop
    end if
  end do
end subroutine compute_packing_fractions

!===============================================================================================================================
! excess_chemical_potential_and_reference_bulk_grand_potential calculates the excess chemical potential and reference bulk grand-pot
! Here we calculate the excess chemical potential which is defined so that the difference in the bulk homogeneous systeme grand-pot.
! and the reference bulk homogeneous grand-potential is zero.
! We also calculate the reference bulk grand potential
SUBROUTINE excess_chemical_potential_and_reference_bulk_grand_potential
  use precision_kinds ,only: dp, i2b
  use constants       ,only: fourpi, pi
  use module_thermo, only: thermo
  use module_solvent, only: solvent
  use module_grid, only: grid
  use module_input           ,only: verbose, getinput
  implicit none
  real(dp) :: n0, n1, n2, n3 ! weighted densities in the case of constant density = ref bulk density
  real(dp) :: dphidn(0:3)    ! partial derivative of phi w.r.t. weighted densities
  real(dp) :: dndrho(0:3)    ! partial derivative of weighted densities w.r.t. density of constituant i
  real(dp) :: R,kT
  integer  :: s
  character(len=4) :: hs_functional
  hs_functional=getinput%char('hs_functional')
  kT=thermo%kbT
  do s=1,size(solvent) ! compute excess chemical potential, so that bulk grand potential is zero for density = constant = ref bulk density
    ! weighted densities in the case of constant density = ref bulk density
    R=hs(s)%R
    n0 = solvent(s)%n0
    n1 =            n0*R
    n2 =            n0*4*pi*R**2
    n3 =            n0*4*pi*R**3/3
    ! partial derivative of phi w.r.t. weighted densities
    if( hs_functional(1:2)=='PY' ) then
      dphidn(0) = -log(1-n3)
      dphidn(1) = n2/(1-n3)
      dphidn(2) = n1/(1-n3) + n2**2/(8*pi*(1-n3)**2)
      dphidn(3) = n0/(1-n3) + n1*n2/(1-n3)**2 + n2**3/(12*pi*(1-n3)**3)
    else if( hs_functional(1:2)=='CS' .or. hs_functional(1:4)=='MCSL' ) then
      dphidn(0) = -log(1-n3)
      dphidn(1) = n2/(1-n3)
      dphidn(2) = (n3*(n2**2-12*n1*(n3-1)*n3*Pi)+n2**2*(n3-1)**2*log(1-n3)) / (12*pi*(n3-1)**2*n3**2)
      dphidn(3) = (n3*(n2**3*(2-5*n3+n3**2) +36*n1*n2*(n3-1)*n3**2*Pi-36*n0*(n3-1)**2*n3**2*Pi) -2*n2**3*(n3-1)**3*log(1-n3)) &
                  /(36*pi*(n3-1)**3*n3**3)
    end if
    ! partial derivative of weighted densities w.r.t. density of constituant i. It may be shown it is weight function (k=0)
    dndrho(0) = 1
    dndrho(1) = R
    dndrho(2) = 4*pi*R**2
    dndrho(3) = 4*pi*R**3/3.
    ! excess chemical potential
    hs(s)%excchempot = kT*sum(dphidn*dndrho)
    if( verbose) then
      write(*,'(A,i1,A,f10.2)') 'chemical potential mu_exc0 (',s,') = ', hs(s)%excchempot
    end if
    ! grand potential of the homogeneous reference fluid
    ! compute reference bulk grand-potential Omega(rho = rho_0) !! Do not forget the solver minimizes Omega[rho]-Omega[rho_0] = Fsolvatation
    IF ( hs_functional(1:2)=='PY' ) THEN
      hs(s)%Fexc0 = kT*( -n0*log(1-n3) + n1*n2/(1-n3) + n2**3/(24*pi*(1-n3)**2) )
    ELSE IF ( hs_functional(1:2)=='CS' .or. hs_functional(1:4)=='MCSL' ) THEN
      hs(s)%Fexc0 = kT*(&
       ((1./(36*pi))*n2**3/n3**2-n0) *log(1-n3) + n1*n2/(1-n3) + (1./(36*pi))*n2**3/((1-n3)**2*n3)    )
    END IF
    ! integration factors
    hs(s)%Fexc0 = (hs(s)%Fexc0 - hs(s)%excchempot * solvent(s)%n0) * PRODUCT(grid%length)
    if( verbose) then
      write(*,'(A,i1,A,f10.2)') 'Fexc0 ( ' , s , ' ) = ' , hs(s)%Fexc0
    end if
  end do
end subroutine excess_chemical_potential_and_reference_bulk_grand_potential

!===================================================================================================================================
! Read hard sphere radius of every constituant of the fluid. If necessary.
! It begins by checking if it has already been allocated. If not it allocates it.
! Then, it reads every line of input_line which contains the inputs for the tag 'hard_sphere_radius"
! It then reads line after line the hard sphere radius of each constituant
subroutine read_hard_sphere_radius
  use precision_kinds  ,only: dp, i2b
  use module_input            ,only: input_line
  use module_thermo ,only: thermo
  use module_solvent, only: solvent
  use module_input           ,only: verbose, getinput
  implicit none
  integer(i2b) :: i,j,s
  real(dp)     :: d_wca ! optimal diameter for hard spheres in the case of lennard jones perturbation as defined by Verlet and Weis, Phys Rev A 1972
  real(dp)     :: dmin
  character(50) :: dummychar

  !allocate(hs(size(solvent)))
  dummychar= getinput%char_multiple('hs_radius', defaultvalue="0.0 0.0 0.0 0.0 0.0 0.0")
  read(dummychar,*) hs(:)%R
! Get hard sphere radius
  ! two possibilities : 1/ pure hard spheres : read radius in input/dft.in    2/ perturbated HS by a lennard jones
  ! next few lines is a test to know if we have
  ! Algo : one looks for line containing 'hard_sphere_radius'. Next solvent(1)%nspec lines contain the radius of each hard sphere species.
  ! If the radius is negative, the user means that the radius has to be calculated using the week chandler anderson model
  ! by reading the lennard jones sigma and epsilon values accordingly to integer species between 1 and solvent(1)%nspec
! read hard sphere radius and if WCA the compute the WCA radius
  !do i= 1, SIZE(input_line)
  !  do s = 1 , solvent(1)%nspec
  !      !print*, "You ask for a WCA calculation of hard sphere diameter"
  !      !CALL compute_wca_diameter (solvent(s)%n0 , thermo%T, solvent(s)%site(1)%sig , solvent(s)%site(1)%eps, d_wca)
  !      !hs(s)%R = d_wca / 2.0_dp
  !  end do
  !  exit
  !end do
end subroutine read_hard_sphere_radius

end module hardspheres

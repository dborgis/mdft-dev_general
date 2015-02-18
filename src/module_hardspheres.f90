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
  public :: packfrac, compute_hard_spheres_parameters

contains

!===================================================================================================================================
! packing fraction : eta = 4/3 * pi * R^3 * number density of the constituant ( /= total solvent density)
pure function packfrac (numberdensity, HSradius)
  use precision_kinds ,only: dp
  use constants       ,only: pi
  implicit none
  real(dp) :: packfrac
  real(dp), parameter  :: fourthird=4./3.
  real(dp), intent(in) :: numberdensity ! density of the constituant (/= total solvant density)
  real(dp), intent(in) :: HSradius ! radius of the constituant
  packfrac = fourthird*pi*numberdensity*HSradius**3
end function packfrac

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
  use precision_kinds ,ONLY: dp, i2b
  use constants       ,ONLY: fourpi, zeroC, zerodp
  use system          ,ONLY: nb_species, spaceGrid
  use fft             ,ONLY: norm_k
  implicit none
  real(dp)     :: kR, sinkR, coskR, norm_k_loc, R, k
  integer(i2b) :: l,m,n,s,nfft1,nfft2,nfft3
  nfft1=spacegrid%n_nodes(1)
  nfft2=spacegrid%n_nodes(2)
  nfft3=spacegrid%n_nodes(3)
  do concurrent( s=1:nb_species )
    allocate( hs(s)%w_k(nfft1/2+1,nfft2,nfft3,0:3) ,source=zerodp)
  end do
  ! Weight functions of a fluid of hard spheres are known analyticaly.
  do concurrent( s=1:nb_species, l=1:nfft1/2+1, m=1:nfft2, n=1:nfft3 )
    R=hs(s)%R
    k=norm_k(l,m,n)
    if( abs(k)>epsilon(1._dp) ) then !  => k /= 0.0_dp
      kR=k*R
      sinkR=sin(kR)
      coskR=cos(kR)
      hs(s)%w_k(l,m,n,3) = fourpi*(sinkR-kR*coskR)/(k**3)
      hs(s)%w_k(l,m,n,2) = fourpi*R*sinkR/k
      hs(s)%w_k(l,m,n,1) = (sinkR+kR*coskR)/(2*k)
      hs(s)%w_k(l,m,n,0) = coskR+kR*sinkR/2
    else
      hs(s)%w_k(l,m,n,3) = fourpi/3*R**3 ! volume
      hs(s)%w_k(l,m,n,2) = fourpi  *R**2 ! surface area
      hs(s)%w_k(l,m,n,1) = R             ! radius
      hs(s)%w_k(l,m,n,0) = 1             ! unity
    end if
  end do
end subroutine weight_functions

  !===================================================================================================================================
  ! Read, allocate, compute the hard sphere parameters
  !===================================================================================================================================
  ! read hard sphere radius
  ! read the mole fraction of each hard sphere
  ! compute their weight functions (using their fundamental measures)
  ! compute their packing fraction
  ! read the excess functional
  ! compute accordingly the chemical potential and the reference bulk density
  !===================================================================================================================================
  SUBROUTINE compute_hard_spheres_parameters

    use system          ,ONLY: nb_species, solvent
    use input           ,ONLY: input_char

    implicit none

    CHARACTER(4) :: hs_functional

    hs_functional = input_char('hs_functional')
    CALL read_hard_sphere_radius_and_allocate_if_necessary ! read hard sphere radius and allocate if necessary
    CALL weight_functions
    CALL compute_packing_fractions_and_check_legality
    CALL check_functional_legality ( hs_functional ) ! Get the free energy functional that should be used. For now Percus Yevick and Carnahan Starling only. May be expanded.

    ! compute excess chemical potential and grand potential at reference bulk density
    CALL excess_chemical_potential_and_reference_bulk_grand_potential (hs_functional)

  CONTAINS

    !===============================================================================================================================
    ! Subroutine that deduct the packing fraction of all reference bulk fluids of the consitutuants of the mixture.
    !===============================================================================================================================
    ! We want the reference bulk fluids to have physical meaning, that's why we pay attention to them not to be greater than 0.74, where
    ! they would be unphysical. 0.74 is the maximum packing of a solid crystal.
    !===============================================================================================================================
    SUBROUTINE compute_packing_fractions_and_check_legality
      use precision_kinds ,ONLY: dp , i2b
      use system          ,ONLY: nb_species
      use input           ,ONLY: verbose
      implicit none
      integer(i2b) :: s
      DO s = 1, nb_species
        hs(s)%pf = packfrac( numberdensity=solvent(s)%n0, HSradius=hs(s)%R )
        IF (verbose) PRINT*,'Packing fraction of species ',s,') is ',hs(s)%pf
        ! compute homogeneous fluid reference with Perkus Yevick
        ! It is important to keep in mind it is the packing fraction of the REFERENCE fluid(s), not a partial packing fraction of our mixture.
        ! although the Percus-Yevick equation shows no singularities for eta < 1 , the region beyond eta = pi / (3 sqrt(2) ) = 0.74 is unphysical, since the fluid then has a packing density greater than that of a closed packed solid.
        if (hs(s)%pf >= 0.492 .and. hs(s)%pf <0.74 ) then
          print*,"packing fraction of species",s,"is",hs(s)%pf,"which is higher than freezing packing fraction of 0.492"
        else IF ( hs(s)%pf >= 0.74_dp ) then
          print*,'packing fraction of species',s,"is",hs(s)%pf,"which is higher than metastable fluid-solid pack. frac. of 0.72"
          STOP
        END IF
      END DO
    END SUBROUTINE compute_packing_fractions_and_check_legality

    !===============================================================================================================================
    ! this SUBROUTINE checks if the functional asked in input file is legal. Else, stop execution.
    !===============================================================================================================================
    SUBROUTINE check_functional_legality ( hs_functional )
      use precision_kinds  ,ONLY: i2b
      use system           ,ONLY: nb_species
      implicit none
      CHARACTER(4), INTENT(INOUT) :: hs_functional
      integer(i2b):: i

      i = 0
      IF ( hs_functional(1:2) == 'CS' ) i = i+1
      IF ( hs_functional(1:2) == 'PY' ) i = i+1
      IF ( hs_functional(1:4) == 'MCSL' ) i = i+1
      ! if the hs functional is unknown, i is still 0
      IF ( i == 0 ) THEN
        PRINT*,'Asked functional is ' , hs_functional(1:4)
        PRINT*,'This functional is not implemeted'
        PRINT*,'Default value will be used: Carnahan-starling for pure fluids and Mansoori-Carnahan-Starling-Leland for multi.'
        IF ( nb_species == 1 ) hs_functional(1:2) = 'CS'
        IF ( nb_species  > 1 ) hs_functional(1:4) = 'MCSL'
      END IF
    END SUBROUTINE check_functional_legality

    !===============================================================================================================================
    ! excess_chemical_potential_and_reference_bulk_grand_potential calculates the excess chemical potential and reference bulk grand-pot
    !===============================================================================================================================
    ! Here we calculate the excess chemical potential which is defined so that the difference in the bulk homogeneous systeme grand-pot.
    ! and the reference bulk homogeneous grand-potential is zero.
    ! We also calculate the reference bulk grand potential
    !===============================================================================================================================
    SUBROUTINE excess_chemical_potential_and_reference_bulk_grand_potential &
      ( hs_functional )

      use precision_kinds ,ONLY: dp, i2b
      use constants       ,ONLY: fourpi, pi
      use system          ,ONLY: thermoCond, spaceGrid
      use input           ,ONLY: verbose

      implicit none

      CHARACTER(4), INTENT(IN) :: hs_functional
      real(dp) :: n0, n1, n2, n3 ! weighted densities in the case of constant density = ref bulk density
      real(dp) :: dphidn(0:3) ! partial derivative of phi w.r.t. weighted densities
      real(dp) :: dndrho(0:3) ! partial derivative of weighted densities w.r.t. density of constituant i
      integer(i2b) :: s ! dummy

      DO s=1,nb_species ! compute excess chemical potential, so that bulk grand potential is zero for density = constant = ref bulk density

        ! weighted densities in the case of constant density = ref bulk density
        n0 = solvent(s)%n0
        n1 = hs(s)%R * n0
        n2 = 4*pi * hs(s)%R**2 *n0
        n3 = 4./3.*pi *hs(s)%R**3 *n0

        ! partial derivative of phi w.r.t. weighted densities
        if ( hs_functional(1:2)=='PY' ) then
          dphidn(0) = -log(1-n3)
          dphidn(1) = n2/(1-n3)
          dphidn(2) = n1/(1-n3) + n2**2/(8*pi*(1-n3)**2)
          dphidn(3) = n0/(1-n3) + n1*n2/(1-n3)**2 + n2**3/(12*pi*(1-n3)**3)

        ELSE IF ( hs_functional(1:2)=='CS' .OR. hs_functional(1:4)=='MCSL' ) THEN
          dphidn(0) = -log(1-n3)
          dphidn(1) = n2/(1-n3)
          dphidn(2) = (n3 * ( n2 ** 2 - 12.0_dp * n1 * ( -1.0_dp + n3 ) * n3 * Pi ) &
          + n2 ** 2 * ( -1.0_dp + n3 ) ** 2 * log ( 1.0_dp - n3 ) )/( 12.0_dp * ( -1.0_dp + n3 ) ** 2 * n3 ** 2 * pi )
          dphidn(3) = ( n3 * ( n2 ** 3 * ( 2.0_dp-5.0_dp*n3 + n3 ** 2 ) + 36.0_dp * n1 * n2 * (-1.0_dp+n3)&
          * n3 ** 2 * Pi - 36.0_dp * n0 * (-1.0_dp + n3) ** 2 * n3**2 * Pi) - 2.0_dp * n2 ** 3 *&
          (-1.0_dp + n3)**3 * log ( 1.0_dp - n3 ) ) / ( 36.0_dp *(-1.0_dp + n3)**3 * n3**3 *Pi)
        END IF

        ! partial derivative of weighted densities w.r.t. density of constituant i. It may be shown it is weight function (k=0)
        dndrho(0) = 1
        dndrho(1) = hs(s)%R
        dndrho(2) = fourpi*hs(s)%R**2
        dndrho(3) = fourpi/3*hs(s)%R**3

        ! excess chemical potential
        hs(s)%excchempot = thermoCond%kbT * SUM(dphidn*dndrho)
        IF (verbose) PRINT*,'chemical potential mu_exc0 ( ' , s , ' ) = ' , hs(s)%excchempot

        ! compute reference bulk grand-potential Omega(rho = rho_0) !! Do not forget the solver minimizes Omega[rho]-Omega[rho_0] = Fsolvatation
        IF ( hs_functional(1:2)=='PY' ) THEN
          hs(s)%Fexc0 = thermoCond%kbT*( -n0*log(1-n3) + n1*n2/(1-n3) + n2**3/(24*pi*(1-n3)**2) )

        ELSE IF ( hs_functional(1:2)=='CS' .or. hs_functional(1:4)=='MCSL' ) THEN
          hs(s)%Fexc0 = thermoCond%kbT*(&
           ((1./(36*pi))*n2**3/n3**2-n0) *log(1-n3) + n1*n2/(1-n3) + (1./(36*pi))*n2**3/((1-n3)**2*n3)    )
        END IF

        ! integration factors
        hs(s)%Fexc0 = (hs(s)%Fexc0 - hs(s)%excchempot * solvent(s)%n0) * PRODUCT(spaceGrid%length)
        IF (verbose) PRINT*,'Fexc0 ( ' , s , ' ) = ' , hs(s)%Fexc0
      END DO
    END SUBROUTINE excess_chemical_potential_and_reference_bulk_grand_potential


  END SUBROUTINE compute_hard_spheres_parameters


END MODULE hardspheres

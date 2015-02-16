MODULE hardspheres

  USE precision_kinds  ,ONLY: dp
  IMPLICIT NONE

  PRIVATE

  TYPE type_hs
    REAL(dp) :: R          ! radius
    REAL(dp) :: excchempot ! excess chemical potential
    REAL(dp) :: Fexc0
    REAL(dp) :: pf         ! packing fraction
    COMPLEX(dp), ALLOCATABLE :: w_k(:,:,:,:)
  END TYPE type_hs
  TYPE(type_hs), PUBLIC, ALLOCATABLE :: hs(:) ! one per constituant

  PUBLIC :: packfrac, populate_weight_functions_in_Fourier_space, compute_hard_spheres_parameters



CONTAINS

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

  !===========================================================================================================================
  ! This SUBROUTINE computes the density independant weight functions as defined by Kierlik and Rosinberg in 1990
  !===========================================================================================================================
  ! The four weight functions are here defined in k-space. They are known analyticaly and only depend on the so called fundamental
  ! measures of the hard spheres.
  ! w_3i(k=0) = V_i the volume of constituant i
  ! w_2i(k=0) = S_i the surface area of constituant i
  ! w_1i(k=0) = R_i the radius of constituant i
  ! w_0i(k=0) = 1
  ! They are scalar numbers in opposition to scalar and vector weight functions by Rosenfeld in its seminal Phys. Rev. Lett.
  ! introducing the fundamental measure theory (FMT)
  !===========================================================================================================================
  SUBROUTINE populate_weight_functions_in_Fourier_space

    USE precision_kinds ,ONLY: dp, i2b
    USE constants       ,ONLY: FourPi, zeroC
    USE system          ,ONLY: nb_species, spaceGrid
    USE fft             ,ONLY: norm_k
    IMPLICIT NONE

    REAL(dp) :: kR , FourPiR , sinkR , coskR, norm_k_local
    INTEGER(i2b) :: l,m,n,s,nfft(3)

    nfft = spaceGrid%n_nodes

    DO CONCURRENT (s=1:nb_species)
      ALLOCATE( hs(s)%w_k(nfft(1)/2+1, nfft(2), nfft(3), 0:3 ) ,SOURCE=zeroC)
    END DO

    ! Weight functions of a fluid of hard spheres are known analyticaly. They are easily defined in Fourier space.
    ! They depend on fundamental measures of the fluid.
    ! Here is the Kierlik and Rosinberg'FMT : 4 scalar weight functions by species. Ref: Kierlik and Rosinberg, PRA 1990

    DO CONCURRENT ( s=1:nb_species, l=1:nfft(1)/2+1, m=1:nfft(2), n=1:nfft(3) )
      FourPiR = FourPi * hs(s)%R
      IF ( l+m+n/=3 ) THEN !  => norm_k_local /= 0.0_dp
        norm_k_local = norm_k (l,m,n)
        kR = norm_k_local * hs(s)%R
        sinkR = SIN(kR)
        coskR = COS(kR)
        hs(s)%w_k(l,m,n,3) = CMPLX(FourPi * ( sinkR - kR * coskR ) / ( norm_k_local ** 3 ) ,0._dp, dp)
        hs(s)%w_k(l,m,n,2) = CMPLX(FourPiR * sinkR / norm_k_local ,0._dp, dp)
        hs(s)%w_k(l,m,n,1) = CMPLX(( sinkR + kR * coskR ) / ( 2.0_dp * norm_k_local ) ,0._dp, dp)
        hs(s)%w_k(l,m,n,0) = CMPLX(coskR + 0.5_dp * kR * sinkR ,0._dp, dp)
      ELSE
        hs(s)%w_k(l,m,n,3) = CMPLX(FourPi / 3.0_dp * hs(s)%R** 3  ,0._dp, dp)! volume
        hs(s)%w_k(l,m,n,2) = CMPLX(FourPi * hs(s)%R** 2  ,0._dp, dp)! surface area
        hs(s)%w_k(l,m,n,1) = CMPLX(hs(s)%R  ,0._dp, dp)! radius
        hs(s)%w_k(l,m,n,0) = CMPLX(1.0_dp  ,0._dp, dp)! unity
      END IF
    END DO
  END SUBROUTINE populate_weight_functions_in_Fourier_space

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

    USE system          ,ONLY: nb_species, solvent
    USE input           ,ONLY: input_char

    IMPLICIT NONE

    CHARACTER(4) :: hs_functional

    hs_functional = input_char('hs_functional')
    CALL read_hard_sphere_radius_and_allocate_if_necessary ! read hard sphere radius and allocate if necessary
    CALL populate_weight_functions_in_Fourier_space
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
      USE precision_kinds ,ONLY: dp , i2b
      USE system          ,ONLY: nb_species
      USE input           ,ONLY: verbose
      IMPLICIT NONE
      INTEGER(i2b) :: s
      DO s = 1, nb_species
        hs(s)%pf = packfrac( numberdensity=solvent(s)%n0, HSradius=hs(s)%R )
        IF (verbose) PRINT*,'Packing fraction of species ',s,') is ',hs(s)%pf
        ! compute homogeneous fluid reference with Perkus Yevick
        ! It is important to keep in mind it is the packing fraction of the REFERENCE fluid(s), not a partial packing fraction of our mixture.
        ! although the Percus-Yevick equation shows no singularities for eta < 1 , the region beyond eta = pi / (3 sqrt(2) ) = 0.74 is unphysical, since the fluid then has a packing density greater than that of a closed packed solid.
        IF ( hs(s)%pf >= 0.74_dp ) then
          PRINT*,'packing fraction of species ',s,'>= 0.74 , ie closed packed solid. unphysical region explored. stop'
          STOP
        END IF
      END DO
    END SUBROUTINE compute_packing_fractions_and_check_legality

    !===============================================================================================================================
    ! this SUBROUTINE checks if the functional asked in input file is legal. Else, stop execution.
    !===============================================================================================================================
    SUBROUTINE check_functional_legality ( hs_functional )
      USE precision_kinds  ,ONLY: i2b
      USE system           ,ONLY: nb_species
      IMPLICIT NONE
      CHARACTER(4), INTENT(INOUT) :: hs_functional
      INTEGER(i2b):: i

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

      USE precision_kinds ,ONLY: dp, i2b
      USE constants       ,ONLY: fourpi, pi
      USE system          ,ONLY: thermoCond, spaceGrid
      USE input           ,ONLY: verbose

      IMPLICIT NONE

      CHARACTER(4), INTENT(IN) :: hs_functional
      REAL(dp) :: n0, n1, n2, n3 ! weighted densities in the case of constant density = ref bulk density
      REAL(dp) :: dphidn(0:3) ! partial derivative of phi w.r.t. weighted densities
      REAL(dp) :: dndrho(0:3) ! partial derivative of weighted densities w.r.t. density of constituant i
      INTEGER(i2b) :: s ! dummy

      DO s=1,nb_species ! compute excess chemical potential, so that bulk grand potential is zero for density = constant = ref bulk density

        ! weighted densities in the case of constant density = ref bulk density
        n0 = 1.0_dp * solvent(s)%n0
        n1 = hs(s)%R * solvent(s)%n0
        n2 = 4.0_dp * pi * hs(s)%R ** 2 * solvent(s)%n0
        n3 = 4.0_dp / 3.0_dp * pi * hs(s)%R ** 3 * solvent(s)%n0

        ! partial derivative of phi w.r.t. weighted densities
        IF ( hs_functional(1:2)=='PY' ) THEN
          dphidn(0) = -log ( 1.0_dp - n3 )
          dphidn(1) = n2 / ( 1.0_dp - n3 )
          dphidn(2) = n1 / ( 1.0_dp - n3 ) + n2 ** 2 / ( 8.0_dp * pi * ( 1.0_dp - n3 ) ** 2 )
          dphidn(3) = n0 / ( 1.0_dp - n3 ) + n1 * n2 / ( 1.0_dp - n3 ) ** 2 &
          - n2 ** 3 / ( 12.0_dp * pi * ( n3 - 1.0_dp ) ** 3 )
        ELSE IF ( hs_functional(1:2)=='CS' .OR. hs_functional(1:4)=='MCSL' ) THEN
          dphidn(0) = - log ( 1.0_dp - n3 )
          dphidn(1) = n2 / ( 1.0_dp - n3 )
          dphidn(2) = ( n3 * ( n2 ** 2 - 12.0_dp * n1 * ( -1.0_dp + n3 ) * n3 * Pi ) &
          + n2 ** 2 * ( -1.0_dp + n3 ) ** 2 * log ( 1.0_dp - n3 ) )/( 12.0_dp * ( -1.0_dp + n3 ) ** 2 * n3 ** 2 * pi )
          dphidn(3) = ( n3 * ( n2 ** 3 * ( 2.0_dp-5.0_dp*n3 + n3 ** 2 ) + 36.0_dp * n1 * n2 * (-1.0_dp+n3)&
          * n3 ** 2 * Pi - 36.0_dp * n0 * (-1.0_dp + n3) ** 2 * n3**2 * Pi) - 2.0_dp * n2 ** 3 *&
          (-1.0_dp + n3)**3 * log ( 1.0_dp - n3 ) ) / ( 36.0_dp *(-1.0_dp + n3)**3 * n3**3 *Pi)
        END IF

        ! partial derivative of weighted densities w.r.t. density of constituant i. It may be shown it is weight function (k=0)
        dndrho(0) = 1.0_dp
        dndrho(1) = hs(s)%R
        dndrho(2) = fourpi * hs(s)%R ** 2
        dndrho(3) = fourpi / 3.0_dp * hs(s)%R ** 3

        ! excess chemical potential
        hs(s)%excchempot = thermoCond%kbT * SUM(dphidn*dndrho)
        IF (verbose) PRINT*,'chemical potential mu_exc0 ( ' , s , ' ) = ' , hs(s)%excchempot

        ! compute reference bulk grand-potential Omega(rho = rho_0) !! Do not forget the solver minimizes Omega[rho]-Omega[rho_0] = Fsolvatation
        IF ( hs_functional(1:2)=='PY' ) THEN
          hs(s)%Fexc0 = thermoCond%kbT * ( - n0 * log ( 1.0_dp - n3 )                            &
          + n1 * n2 / ( 1.0_dp - n3 )                           &
          + n2 ** 3 / ( 24.0_dp * pi * ( 1.0_dp - n3 ) ** 2 ) )
        ELSE IF ( hs_functional(1:2)=='CS' .or. hs_functional(1:4)=='MCSL' ) THEN
          hs(s)%Fexc0 = thermoCond%kbT * ( ( ( 1.0_dp / ( 36.0_dp * pi ) ) * n2 ** 3 / n3 ** 2 - n0 ) * log ( 1.0_dp - n3 ) &
          + n1 * n2 / ( 1.0_dp - n3 )                                                &
          + ( 1.0_dp / ( 36.0_dp * pi ) ) * n2 ** 3 / ( ( 1.0_dp - n3 ) ** 2 * n3 )   )
        END IF
        ! integration factors
        hs(s)%Fexc0 = (hs(s)%Fexc0 - hs(s)%excchempot * solvent(s)%n0) * PRODUCT(spaceGrid%length)
        IF (verbose) PRINT*,'Fexc0 ( ' , s , ' ) = ' , hs(s)%Fexc0
      END DO
    END SUBROUTINE excess_chemical_potential_and_reference_bulk_grand_potential


  END SUBROUTINE compute_hard_spheres_parameters


END MODULE hardspheres

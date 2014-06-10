!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Read, allocate, compute the hard sphere parameters
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! read hard sphere radius
! read the mole fraction of each hard sphere
! compute their weight functions (using their fundamental measures)
! compute their packing fraction
! read the excess functional
! compute accordingly the chemical potential and the reference bulk density
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE compute_hard_spheres_parameters
    
    USE precision_kinds ,ONLY: dp
    USE system          ,ONLY: nb_species, muexc_0_multispec, Fexc_0_multispec, n_0_multispec
    USE hardspheres     ,ONLY: populate_weight_functions_in_Fourier_space, hs
    IMPLICIT NONE
    
    CHARACTER(4) :: hs_functional
    
    CALL read_hard_sphere_radius_and_allocate_if_necessary ! read hard sphere radius and allocate if necessary
    CALL populate_weight_functions_in_Fourier_space
    CALL compute_packing_fractions_and_check_legality ! compute packing fraction of each constituant and the total packing fraction in order to check if the calculation is physical
    CALL check_functional_legality ( hs_functional ) ! Get the free energy functional that should be used. For now Percus Yevick and Carnahan Starling only. May be expanded.

    ! compute excess chemical potential and grand potential at reference bulk density
    ALLOCATE ( muexc_0_multispec(nb_species) ,SOURCE=0._dp )
    ALLOCATE ( Fexc_0_multispec (nb_species) ,SOURCE=0._dp )
    CALL excess_chemical_potential_and_reference_bulk_grand_potential &
                ( nb_species , n_0_multispec , muexc_0_multispec , Fexc_0_multispec , hs_functional )

    CONTAINS
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine that deduct the packing fraction of all reference bulk fluids of the consitutuants of the mixture.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! We want the reference bulk fluids to have physical meaning, that's why we pay attention to them not to be greater than 0.74, where
    ! they would be unphysical. 0.74 is the maximum packing of a solid crystal.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE compute_packing_fractions_and_check_legality
        USE precision_kinds ,ONLY: dp , i2b
        USE constants       ,ONLY: fourpi
        USE system          ,ONLY: nb_species , n=>n_0_multispec
        USE input           ,ONLY: verbose
        USE hardspheres     ,ONLY: packfrac, hs
        IMPLICIT NONE
        INTEGER(i2b) :: s
        DO s = 1, nb_species
            IF (verbose) PRINT*,'Packing fraction of species ',s,') is ',packfrac(n(s),hs(s)%R)
            ! compute homogeneous fluid reference with Perkus Yevick
            ! It is important to keep in mind it is the packing fraction of the REFERENCE fluid(s), not a partial packing fraction of our mixture.
            ! although the Percus-Yevick equation shows no singularities for eta < 1 , the region beyond eta = pi / (3 sqrt(2) ) = 0.74 is unphysical, since the fluid then has a packing density greater than that of a closed packed solid.
            IF ( packfrac(n(s),hs(s)%R) >= 0.74_dp ) then
                PRINT*,'packing fraction of species ',s, '>= 0.74 , ie closed packed solid. unphysical region explored. stop'
                STOP
            END IF
        END DO
    END SUBROUTINE compute_packing_fractions_and_check_legality
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! this SUBROUTINE checks if the functional asked in input file is legal. Else, stop execution.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! excess_chemical_potential_and_reference_bulk_grand_potential calculates the excess chemical potential and reference bulk grand-pot
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Here we calculate the excess chemical potential which is defined so that the difference in the bulk homogeneous systeme grand-pot.
    ! and the reference bulk homogeneous grand-potential is zero.
    ! We also calculate the reference bulk grand potential
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE excess_chemical_potential_and_reference_bulk_grand_potential &
                ( nb_species , n_0_multispec , muexc_0_multispec , Fexc_0_multispec , hs_functional )
        
        USE precision_kinds ,ONLY: dp, i2b
        USE constants       ,ONLY: fourpi, pi
        USE system          ,ONLY: kbT, spaceGrid
        USE input           ,ONLY: verbose
        USE hardspheres     ,ONLY: hs
        
        IMPLICIT NONE
    
        INTEGER(i2b), INTENT(in) :: nb_species
        REAL(dp), dimension (nb_species) , intent(in) :: n_0_multispec ! ref bulk densities
        CHARACTER(4), intent(in) :: hs_functional
        REAL(dp), dimension (nb_species) , intent(out) :: muexc_0_multispec ! excess chemical potential defined so that grand potential is zero at ref bulk density
        REAL(dp), dimension (nb_species) , intent(out) :: Fexc_0_multispec ! excess helmotz free energy of reference bulk system
        INTEGER(i2b):: s ! dummy
        REAL(dp):: n0 , n1 , n2 , n3 ! weighted densities in the case of constant density = ref bulk density
        REAL(dp):: partial_phi_over_partial_n0 , partial_phi_over_partial_n1 ! partial derivative of phi w.r.t. weighted densities
        REAL(dp):: partial_phi_over_partial_n2 , partial_phi_over_partial_n3
        REAL(dp):: partial_n0_over_partial_rho , partial_n1_over_partial_rho ! partial derivative of weighted densities w.r.t. density of constituant i
        REAL(dp):: partial_n2_over_partial_rho , partial_n3_over_partial_rho
        REAL(dp) :: lx,ly,lz
        lx = spaceGrid%length(1)
        ly = spaceGrid%length(2)
        lz = spaceGrid%length(3)

        DO s=1,nb_species ! compute excess chemical potential, so that bulk grand potential is zero for density = constant = ref bulk density

            ! weighted densities in the case of constant density = ref bulk density
            n0 = 1.0_dp * n_0_multispec(s)
            n1 = hs(s)%R * n_0_multispec(s)
            n2 = 4.0_dp * pi * hs(s)%R ** 2 * n_0_multispec(s)
            n3 = 4.0_dp / 3.0_dp * pi * hs(s)%R ** 3 * n_0_multispec(s)
            ! partial derivative of phi w.r.t. weighted densities
            IF ( hs_functional ( 1 : 2 ) == 'PY' ) THEN
                partial_phi_over_partial_n0 = -log ( 1.0_dp - n3 )
                partial_phi_over_partial_n1 = n2 / ( 1.0_dp - n3 )
                partial_phi_over_partial_n2 = n1 / ( 1.0_dp - n3 ) + n2 ** 2 / ( 8.0_dp * pi * ( 1.0_dp - n3 ) ** 2 )
                partial_phi_over_partial_n3 = n0 / ( 1.0_dp - n3 ) + n1 * n2 / ( 1.0_dp - n3 ) ** 2 &
                                            - n2 ** 3 / ( 12.0_dp * pi * ( n3 - 1.0_dp ) ** 3 )
            ELSE IF ( hs_functional ( 1 : 2 ) == 'CS' .or. hs_functional ( 1 : 4 ) == 'MCSL' ) THEN
                partial_phi_over_partial_n0 = - log ( 1.0_dp - n3 )
                partial_phi_over_partial_n1 = n2 / ( 1.0_dp - n3 )
                partial_phi_over_partial_n2 = ( n3 * ( n2 ** 2 - 12.0_dp * n1 * ( -1.0_dp + n3 ) * n3 * Pi ) &
                        + n2 ** 2 * ( -1.0_dp + n3 ) ** 2 * log ( 1.0_dp - n3 ) )/( 12.0_dp * ( -1.0_dp + n3 ) ** 2 * n3 ** 2 * pi )
                partial_phi_over_partial_n3 = ( n3 * ( n2 ** 3 * ( 2.0_dp-5.0_dp*n3 + n3 ** 2 ) + 36.0_dp * n1 * n2 * (-1.0_dp+n3)&
                        * n3 ** 2 * Pi - 36.0_dp * n0 * (-1.0_dp + n3) ** 2 * n3**2 * Pi) - 2.0_dp * n2 ** 3 *&
                        (-1.0_dp + n3)**3 * log ( 1.0_dp - n3 ) ) / ( 36.0_dp *(-1.0_dp + n3)**3 * n3**3 *Pi)
            END IF
            ! partial derivative of weighted densities w.r.t. density of constituant i. It may be shown it is weight function (k=0)
            partial_n0_over_partial_rho = 1.0_dp
            partial_n1_over_partial_rho = hs(s)%R
            partial_n2_over_partial_rho = fourpi * hs(s)%R ** 2
            partial_n3_over_partial_rho = fourpi / 3.0_dp * hs(s)%R ** 3
            ! excess chemical potential
            muexc_0_multispec(s) = kBT * ( partial_phi_over_partial_n0 * partial_n0_over_partial_rho &
                                                    + partial_phi_over_partial_n1 * partial_n1_over_partial_rho &
                                                    + partial_phi_over_partial_n2 * partial_n2_over_partial_rho &
                                                    + partial_phi_over_partial_n3 * partial_n3_over_partial_rho )
            IF (verbose) PRINT*,'chemical potential mu_exc0 ( ' , s , ' ) = ' , muexc_0_multispec(s)
            ! compute reference bulk grand-potential Omega(rho = rho_0) !! Do not forget the solver minimizes Omega[rho]-Omega[rho_0] = Fsolvatation
            IF ( hs_functional ( 1 : 2 ) == 'PY' ) THEN
                Fexc_0_multispec(s) = kBT * ( - n0 * log ( 1.0_dp - n3 )                            &
                                                    + n1 * n2 / ( 1.0_dp - n3 )                           &
                                                    + n2 ** 3 / ( 24.0_dp * pi * ( 1.0_dp - n3 ) ** 2 ) )
            ELSE IF ( hs_functional ( 1 : 2 ) == 'CS' .or. hs_functional ( 1 : 4 ) == 'MCSL' ) THEN
                Fexc_0_multispec(s) = kBT * ( ( ( 1.0_dp / ( 36.0_dp * pi ) ) * n2 ** 3 / n3 ** 2 - n0 ) * log ( 1.0_dp - n3 ) &
                                                    + n1 * n2 / ( 1.0_dp - n3 )                                                &
                                                    + ( 1.0_dp / ( 36.0_dp * pi ) ) * n2 ** 3 / ( ( 1.0_dp - n3 ) ** 2 * n3 )   )
            END IF
            ! integration factors
            Fexc_0_multispec(s) = Fexc_0_multispec(s) * Lx * Ly * Lz &
                                        - muexc_0_multispec ( s) * Lx * Ly * Lz * n_0_multispec(s)  ! integration factor
            IF (verbose) PRINT*,'Fexc0 ( ' , s , ' ) = ' , Fexc_0_multispec(s)
        END DO
        END SUBROUTINE excess_chemical_potential_and_reference_bulk_grand_potential
    

END SUBROUTINE compute_hard_spheres_parameters

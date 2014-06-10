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
    USE system ,ONLY : radius, nb_species, muexc_0_multispec, Fexc_0_multispec, n_0_multispec
    IMPLICIT NONE
    
    CHARACTER(4) :: hs_functional
    
    CALL read_hard_sphere_radius_and_allocate_if_necessary ! read hard sphere radius and allocate if necessary
    CALL compute_hard_sphere_weight_functions_k ! Compute hard sphere weight functions in Fourier space
    CALL compute_packing_fractions_and_check_legality ! compute packing fraction of each constituant and the total packing fraction in order to check if the calculation is physical
    CALL read_hs_functional ( hs_functional ) ! Get the free energy functional that should be used. For now Percus Yevick and Carnahan Starling only. May be expanded.

    ! compute excess chemical potential and grand potential at reference bulk density
    ALLOCATE ( muexc_0_multispec(nb_species) ,SOURCE=0._dp )
    ALLOCATE ( Fexc_0_multispec (nb_species) ,SOURCE=0._dp )
    CALL excess_chemical_potential_and_reference_bulk_grand_potential &
                ( nb_species , n_0_multispec , radius , muexc_0_multispec , Fexc_0_multispec , hs_functional )

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
        USE system          ,ONLY: nb_species , n_0_multispec , radius
        USE input           ,ONLY: verbose
        IMPLICIT NONE
        REAL(dp) :: eta(nb_species)
        INTEGER(i2b) :: s
        ! packing fraction : eta = 4/3 * pi * R^3 * solvent density of the constituant ( /= total solvent density)
        DO s = 1 , nb_species
            eta ( s ) = fourpi/3.0_dp * n_0_multispec(s) * radius(s)**3
            IF (verbose) PRINT*,'Packing fraction of species ',s,') is ',eta(s)
            ! compute homogeneous fluid reference with Perkus Yevick
            ! It is important to keep in mind it is the packing fraction of the REFERENCE fluid(s), not a partial packing fraction of our mixture.
            ! although the Percus-Yevick equation shows no singularities for eta < 1 , the region beyond eta = pi / (3 sqrt(2) ) = 0.74 is unphysical, since the fluid then has a packing density greater than that of a closed packed solid.
            IF ( eta ( s ) >= 0.74_dp ) then
                write (*,*) 'packing fraction of species ' ,s, ' >= 0.74 , ie closed packed solid. unphysical region explored. stop'
                stop
            END IF
        END DO
    END SUBROUTINE compute_packing_fractions_and_check_legality
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Reads hard sphere excess functional
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE read_hs_functional ( hs_functional )
        USE precision_kinds,only : i2b
        use input,only : input_line, input_char
        IMPLICIT NONE
        integer(i2b):: i , j ! dummy
        character ( 4 ) , intent(out) :: hs_functional
    !    do i = 1 , size ( input_line )
    !      j = len ( 'hs_functional' )
    !      if ( input_line (i) (1:j) == 'hs_functional' ) read ( input_line (i) (j+4:j+7) , * ) hs_functional
    !    END DO
        hs_functional=input_char('hs_functional')
        ! Check legality of hs_functional
        call check_functional_legality ( hs_functional )
    END SUBROUTINE read_hs_functional
        
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! this SUBROUTINE checks if the functional asked in input file is legal. Else, stop execution.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE check_functional_legality ( hs_functional )
        USE precision_kinds,only : i2b
        use system,only : nb_species
        IMPLICIT NONE
        integer(i2b):: i ! dummy
        character ( 4 ) , intent ( inout ) :: hs_functional
        i = 0
        if ( hs_functional(1:2) == 'CS' ) i = i+1
        if ( hs_functional(1:2) == 'PY' ) i = i+1
        if ( hs_functional(1:4) == 'MCSL' ) i = i+1
        ! if the hs functional is unknown, i is still 0
        if ( i == 0 ) then
        write (*,*) 'Asked functional is ' , hs_functional(1:4)
        write (*,*) 'This functional is not implemeted'
        write (*,*) 'Default value will be used: Carnahan-starling for pure fluids and Mansoori-Carnahan-Starling-Leland for multi.'
        if ( nb_species == 1 ) hs_functional(1:2) = 'CS'
        if ( nb_species > 1 ) hs_functional(1:4) = 'MCSL'
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
                ( nb_species , n_0_multispec , radius , muexc_0_multispec , Fexc_0_multispec , hs_functional )
        USE precision_kinds,only : dp , i2b
        use constants,only : fourpi , pi
        use system,only : kbT , Lx , Ly , Lz
        USE input, ONLY: verbose
        IMPLICIT NONE
    
        integer(i2b), intent(in) :: nb_species
        real(dp), dimension ( nb_species ) , intent(in) :: n_0_multispec ! ref bulk densities
        real(dp), dimension ( nb_species ) , intent(in) :: radius ! hard sphere radius
        character ( 4 ) , intent(in) :: hs_functional
        real(dp), dimension ( nb_species ) , intent(out) :: muexc_0_multispec ! excess chemical potential defined so that grand potential is zero at ref bulk density
        real(dp), dimension ( nb_species ) , intent(out) :: Fexc_0_multispec ! excess helmotz free energy of reference bulk system
        integer(i2b):: s ! dummy
        real(dp):: n0 , n1 , n2 , n3 ! weighted densities in the case of constant density = ref bulk density
        real(dp):: partial_phi_over_partial_n0 , partial_phi_over_partial_n1 ! partial derivative of phi w.r.t. weighted densities
        real(dp):: partial_phi_over_partial_n2 , partial_phi_over_partial_n3
        real(dp):: partial_n0_over_partial_rho , partial_n1_over_partial_rho ! partial derivative of weighted densities w.r.t. density of constituant i
        real(dp):: partial_n2_over_partial_rho , partial_n3_over_partial_rho
        ! compute excess chemical potential, so that bulk grand potential is zero for density = constant = ref bulk density
        do s = 1 , nb_species
        
        ! weighted densities in the case of constant density = ref bulk density
        n0 = 1.0_dp * n_0_multispec ( s )
        n1 = radius ( s ) * n_0_multispec ( s )
        n2 = 4.0_dp * pi * radius ( s ) ** 2 * n_0_multispec ( s )
        n3 = 4.0_dp / 3.0_dp * pi * radius ( s ) ** 3 * n_0_multispec ( s )
        ! partial derivative of phi w.r.t. weighted densities
        if ( hs_functional ( 1 : 2 ) == 'PY' ) then
            partial_phi_over_partial_n0 = - log ( 1.0_dp - n3 )
            partial_phi_over_partial_n1 = n2 / ( 1.0_dp - n3 )
            partial_phi_over_partial_n2 = n1 / ( 1.0_dp - n3 ) + n2 ** 2 / ( 8.0_dp * pi * ( 1.0_dp - n3 ) ** 2 )
            partial_phi_over_partial_n3 = n0 / ( 1.0_dp - n3 ) + n1 * n2 / ( 1.0_dp - n3 ) ** 2 &
                                        - n2 ** 3 / ( 12.0_dp * pi * ( n3 - 1.0_dp ) ** 3 )
        ELSE IF ( hs_functional ( 1 : 2 ) == 'CS' .or. hs_functional ( 1 : 4 ) == 'MCSL' ) then
            partial_phi_over_partial_n0 = - log ( 1.0_dp - n3 )
            partial_phi_over_partial_n1 = n2 / ( 1.0_dp - n3 )
            partial_phi_over_partial_n2 = ( n3 * ( n2 ** 2 - 12.0_dp * n1 * ( -1.0_dp + n3 ) * n3 * Pi ) &
                    + n2 ** 2 * ( -1.0_dp + n3 ) ** 2 * log ( 1.0_dp - n3 ) ) / ( 12.0_dp * ( -1.0_dp + n3 ) ** 2 * n3 ** 2 * pi )
            partial_phi_over_partial_n3 = ( n3 * ( n2 ** 3 * ( 2.0_dp - 5.0_dp * n3 + n3 ** 2 ) + 36.0_dp * n1 * n2 * (-1.0_dp+n3)&
                    * n3 ** 2 * Pi - 36.0_dp * n0 * (-1.0_dp + n3) ** 2 * n3**2 * Pi) - 2.0_dp * n2 ** 3 *&
                    (-1.0_dp + n3)**3 * log ( 1.0_dp - n3 ) ) / ( 36.0_dp *(-1.0_dp + n3)**3 * n3**3 *Pi)
        END IF
        ! partial derivative of weighted densities w.r.t. density of constituant i. It may be shown it is weight function (k=0)
        partial_n0_over_partial_rho = 1.0_dp
        partial_n1_over_partial_rho = radius ( s )
        partial_n2_over_partial_rho = fourpi * radius ( s ) ** 2
        partial_n3_over_partial_rho = fourpi / 3.0_dp * radius ( s ) ** 3
        ! excess chemical potential
        muexc_0_multispec ( s ) = kBT * ( partial_phi_over_partial_n0 * partial_n0_over_partial_rho &
                                                + partial_phi_over_partial_n1 * partial_n1_over_partial_rho &
                                                + partial_phi_over_partial_n2 * partial_n2_over_partial_rho &
                                                + partial_phi_over_partial_n3 * partial_n3_over_partial_rho )
        IF (verbose) write ( * , * ) 'chemical potential mu_exc0 ( ' , s , ' ) = ' , muexc_0_multispec ( s )
        ! compute reference bulk grand-potential Omega(rho = rho_0) !! Do not forget the solver minimizes Omega[rho]-Omega[rho_0] = Fsolvatation
        if ( hs_functional ( 1 : 2 ) == 'PY' ) then
            Fexc_0_multispec ( s ) = kBT * ( - n0 * log ( 1.0_dp - n3 )                            &
                                                + n1 * n2 / ( 1.0_dp - n3 )                           &
                                                + n2 ** 3 / ( 24.0_dp * pi * ( 1.0_dp - n3 ) ** 2 ) )
        ELSE IF ( hs_functional ( 1 : 2 ) == 'CS' .or. hs_functional ( 1 : 4 ) == 'MCSL' ) then
            Fexc_0_multispec ( s ) = kBT * ( ( ( 1.0_dp / ( 36.0_dp * pi ) ) * n2 ** 3 / n3 ** 2 - n0 ) * log ( 1.0_dp - n3 ) &
                                                + n1 * n2 / ( 1.0_dp - n3 )                                                        &
                                                + ( 1.0_dp / ( 36.0_dp * pi ) ) * n2 ** 3 / ( ( 1.0_dp - n3 ) ** 2 * n3 )          )
        END IF
        ! integration factors
        Fexc_0_multispec ( s ) = Fexc_0_multispec ( s ) * Lx * Ly * Lz &
                                    - muexc_0_multispec ( s) * Lx * Ly * Lz * n_0_multispec ( s )  ! integration factor
        IF (verbose) write ( * , * ) 'Fexc0 ( ' , s , ' ) = ' , Fexc_0_multispec ( s )
        END DO
        END SUBROUTINE excess_chemical_potential_and_reference_bulk_grand_potential
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! This SUBROUTINE computes the density independant weight functions as defined by Kierlik and Rosinberg in 1990
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! The four weight functions are here defined in k-space. They are known analyticaly and only depend on the so called fundamental
    ! measures of the hard spheres.
    ! w_3i(k=0) = V_i the volume of constituant i
    ! w_2i(k=0) = S_i the surface area of constituant i
    ! w_1i(k=0) = R_i the radius of constituant i
    ! w_0i(k=0) = 1
    ! They are scalar numbers in opposition to scalar and vector weight functions by Rosenfeld in its seminal Phys. Rev. Lett.
    ! introducing the fundamental measure theory (FMT)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE compute_hard_sphere_weight_functions_k
    
            USE precision_kinds,only : dp , i2b
            USE constants,only : FourPi
            USE system,only : nb_species, radius, spaceGrid,&
                        weight_function_3_k , weight_function_2_k , weight_function_1_k , weight_function_0_k
            USE fft,only : norm_k
    
            IMPLICIT NONE
    
            REAL(dp) :: kR , FourPiR , sinkR , coskR, norm_k_local ! dummy for speeding up
            INTEGER(i2b):: l , m , n, s ! dummy for loops
            INTEGER(i2b), DIMENSION(3) :: nfft
            nfft = spaceGrid%n_nodes
            ALLOCATE ( weight_function_3_k ( nfft(1)/2+1, nfft(2), nfft(3), nb_species ) )
            ALLOCATE ( weight_function_2_k ( nfft(1)/2+1, nfft(2), nfft(3), nb_species ) )
            ALLOCATE ( weight_function_1_k ( nfft(1)/2+1, nfft(2), nfft(3), nb_species ) )
            ALLOCATE ( weight_function_0_k ( nfft(1)/2+1, nfft(2), nfft(3), nb_species ) )
    
            ! density weights for hard spheres are known analyticaly
            ! they only depends on fundamental measures of hard spheres
            ! here is the Kierlik and Rosinberg FMT : 4 scalar weight function by species
            ! Compute hard sphere weight functions as defined by Kierlik and Rosinberg, PRA1990
            DO CONCURRENT ( s=1:nb_species, l=1:nfft(1)/2+1, m=1:nfft(2), n=1:nfft(3) )
                FourPiR = FourPi * radius (s)
                norm_k_local = norm_k (l,m,n)
                IF ( norm_k_local /= 0.0_dp ) THEN
                    kR = norm_k_local * radius (s)
                    sinkR = sin ( kR )
                    coskR = cos ( kR )
                    weight_function_3_k (l,m,n,s) = FourPi * ( sinkR - kR * coskR ) / ( norm_k_local ** 3 )
                    weight_function_2_k (l,m,n,s) = FourPiR * sinkR / norm_k_local
                    weight_function_1_k (l,m,n,s) = ( sinkR + kR * coskR ) / ( 2.0_dp * norm_k_local )
                    weight_function_0_k (l,m,n,s) = coskR + 0.5_dp * kR * sinkR
                ELSE
                    weight_function_3_k (l,m,n,s) = FourPi / 3.0_dp * radius(s)** 3 ! volume
                    weight_function_2_k (l,m,n,s) = FourPi * radius(s)** 2 ! surface area
                    weight_function_1_k (l,m,n,s) = radius(s) ! radius
                    weight_function_0_k (l,m,n,s) = 1.0_dp ! unity
                END IF
            END DO
        END SUBROUTINE compute_hard_sphere_weight_functions_k

END SUBROUTINE compute_hard_spheres_parameters

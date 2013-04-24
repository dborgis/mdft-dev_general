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

subroutine compute_hard_spheres_parameters


use system , only : radius , nb_species , muexc_0_multispec , Fexc_0_multispec , n_0_multispec


implicit none

character ( 4 ) :: hs_functional



!> Warn user

write(*,*)'>>> Compute hard sphere parameters in compute_hs_parameters.f90'




! read hard sphere radius and allocate if necessary

call read_hard_sphere_radius_and_allocate_if_necessary



!> Compute hard sphere weight functions in k-space. There are needed for later convolutions for calculation of weighted densities n_{alpha}^i

call compute_hard_sphere_weight_functions_k




! compute packing fraction of each constituant and the total packing fraction in order to check if the calculation is physical

call compute_packing_fractions_and_check_legality




! Get the free energy functional that should be used. For now Percus Yevick and Carnahan Starling only. May be expanded.

call read_hs_functional ( hs_functional )



! compute excess chemical potential and grand potential at reference bulk density

allocate ( muexc_0_multispec ( nb_species ) )

allocate ( Fexc_0_multispec ( nb_species ) )

call excess_chemical_potential_and_reference_bulk_grand_potential &
               ( nb_species , n_0_multispec , radius , muexc_0_multispec , Fexc_0_multispec , hs_functional )


















contains



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subroutine that deduct the packing fraction of all reference bulk fluids of the consitutuants of the mixture.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! We want the reference bulk fluids to have physical meaning, that's why we pay attention to them not to be greater than 0.74, where
! they would be unphysical. 0.74 is the maximum packing of a solid crystal.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine compute_packing_fractions_and_check_legality

    use precision_kinds , only : dp , i2b

    use constants , only : fourpi

    use system , only : nb_species , n_0_multispec , radius

    implicit none

    real ( kind = dp ) , dimension ( nb_species ) :: eta

    integer ( kind = i2b ) :: species ! dummy



    ! packing fraction : eta = 4/3 * pi * R^3 * solvent density of the constituant ( /= total solvent density)

    do species = 1 , nb_species

      eta ( species ) = fourpi / 3.0_dp * n_0_multispec ( species ) * radius ( species ) ** 3

      write ( * , * ) 'Packing fraction (eta) (' , species , ') = ' , eta ( species )



    ! compute homogeneous fluid reference with Perkus Yevick

    ! It is important to keep in mind it is the packing fraction of the REFERENCE fluid(s), not a partial packing fraction of our mixture.

    ! although the Percus-Yevick equation shows no singularities for eta < 1 , the region beyond eta = pi / (3 sqrt(2) ) = 0.74 is unphysical, since the fluid then has a packing density greater than that of a closed packed solid.

      if ( eta ( species ) >= 0.74_dp ) then

        write (*,*) 'packing fraction of species ' , species , ' >= 0.74 , ie closed packed solid. unphysical region explored. stop'

        stop

      end if

    end do

  end subroutine compute_packing_fractions_and_check_legality












!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Reads hard sphere excess functional
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine read_hs_functional ( hs_functional )

    use precision_kinds , only : i2b

    use input , only : input_line, input_char

    implicit none

    integer ( kind = i2b ) :: i , j ! dummy

    character ( 4 ) , intent ( out ) :: hs_functional

!    do i = 1 , size ( input_line )

!      j = len ( 'hs_functional' )

!      if ( input_line (i) (1:j) == 'hs_functional' ) read ( input_line (i) (j+4:j+7) , * ) hs_functional

!    end do
    hs_functional=trim(adjustl(input_char('hs_functional')))

    ! Check legality of hs_functional

    call check_functional_legality ( hs_functional )

  end subroutine read_hs_functional
    








!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! this subroutine checks if the functional asked in input file is legal. Else, stop execution.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine check_functional_legality ( hs_functional )

    use precision_kinds , only : i2b

    use system , only : nb_species

    implicit none

    integer ( kind = i2b ) :: i ! dummy

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

    end if

  end subroutine check_functional_legality










!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! excess_chemical_potential_and_reference_bulk_grand_potential calculates the excess chemical potential and reference bulk grand-pot
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Here we calculate the excess chemical potential which is defined so that the difference in the bulk homogeneous systeme grand-pot.
! and the reference bulk homogeneous grand-potential is zero.
! We also calculate the reference bulk grand potential
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine excess_chemical_potential_and_reference_bulk_grand_potential &
               ( nb_species , n_0_multispec , radius , muexc_0_multispec , Fexc_0_multispec , hs_functional )

    use precision_kinds , only : dp , i2b

    use constants , only : fourpi , pi

    use system , only : kbT , Lx , Ly , Lz

    implicit none
 
    integer ( kind = i2b ) , intent ( in ) :: nb_species

    real ( kind = dp ) , dimension ( nb_species ) , intent ( in ) :: n_0_multispec ! ref bulk densities

    real ( kind = dp ) , dimension ( nb_species ) , intent ( in ) :: radius ! hard sphere radius

    character ( 4 ) , intent ( in ) :: hs_functional

    real ( kind = dp ) , dimension ( nb_species ) , intent ( out ) :: muexc_0_multispec ! excess chemical potential defined so that grand potential is zero at ref bulk density

    real ( kind = dp ) , dimension ( nb_species ) , intent ( out ) :: Fexc_0_multispec ! excess helmotz free energy of reference bulk system

    integer ( kind = i2b ) :: species ! dummy

    real ( kind = dp ) :: n0 , n1 , n2 , n3 ! weighted densities in the case of constant density = ref bulk density

    real ( kind = dp ) :: partial_phi_over_partial_n0 , partial_phi_over_partial_n1 ! partial derivative of phi w.r.t. weighted densities

    real ( kind = dp ) :: partial_phi_over_partial_n2 , partial_phi_over_partial_n3

    real ( kind = dp ) :: partial_n0_over_partial_rho , partial_n1_over_partial_rho ! partial derivative of weighted densities w.r.t. density of constituant i

    real ( kind = dp ) :: partial_n2_over_partial_rho , partial_n3_over_partial_rho




    ! compute excess chemical potential, so that bulk grand potential is zero for density = constant = ref bulk density

    do species = 1 , nb_species

    
      ! weighted densities in the case of constant density = ref bulk density

      n0 = 1.0_dp * n_0_multispec ( species )

      n1 = radius ( species ) * n_0_multispec ( species )

      n2 = 4.0_dp * pi * radius ( species ) ** 2 * n_0_multispec ( species )

      n3 = 4.0_dp / 3.0_dp * pi * radius ( species ) ** 3 * n_0_multispec ( species )



      ! partial derivative of phi w.r.t. weighted densities

      if ( hs_functional ( 1 : 2 ) == 'PY' ) then

        partial_phi_over_partial_n0 = - log ( 1.0_dp - n3 )

        partial_phi_over_partial_n1 = n2 / ( 1.0_dp - n3 )

        partial_phi_over_partial_n2 = n1 / ( 1.0_dp - n3 ) + n2 ** 2 / ( 8.0_dp * pi * ( 1.0_dp - n3 ) ** 2 )

        partial_phi_over_partial_n3 = n0 / ( 1.0_dp - n3 ) + n1 * n2 / ( 1.0_dp - n3 ) ** 2 &
                                    - n2 ** 3 / ( 12.0_dp * pi * ( n3 - 1.0_dp ) ** 3 )

      else if ( hs_functional ( 1 : 2 ) == 'CS' .or. hs_functional ( 1 : 4 ) == 'MCSL' ) then

        partial_phi_over_partial_n0 = - log ( 1.0_dp - n3 )

        partial_phi_over_partial_n1 = n2 / ( 1.0_dp - n3 )

        partial_phi_over_partial_n2 = ( n3 * ( n2 ** 2 - 12.0_dp * n1 * ( -1.0_dp + n3 ) * n3 * Pi ) &
                   + n2 ** 2 * ( -1.0_dp + n3 ) ** 2 * log ( 1.0_dp - n3 ) ) / ( 12.0_dp * ( -1.0_dp + n3 ) ** 2 * n3 ** 2 * pi )

        partial_phi_over_partial_n3 = ( n3 * ( n2 ** 3 * ( 2.0_dp - 5.0_dp * n3 + n3 ** 2 ) + 36.0_dp * n1 * n2 * ( -1.0_dp + n3 )&
                  * n3 ** 2 * Pi - 36.0_dp * n0 * (-1.0_dp + n3) ** 2 * n3**2 * Pi) - 2.0_dp * n2 ** 3 *&
                  (-1.0_dp + n3)**3 * log ( 1.0_dp - n3 ) ) / ( 36.0_dp *(-1.0_dp + n3)**3 * n3**3 *Pi)

      end if



      ! partial derivative of weighted densities w.r.t. density of constituant i. It may be shown it is weight function (k=0)

      partial_n0_over_partial_rho = 1.0_dp

      partial_n1_over_partial_rho = radius ( species )

      partial_n2_over_partial_rho = fourpi * radius ( species ) ** 2

      partial_n3_over_partial_rho = fourpi / 3.0_dp * radius ( species ) ** 3



      ! excess chemical potential

      muexc_0_multispec ( species ) = kBT * ( partial_phi_over_partial_n0 * partial_n0_over_partial_rho &
                                            + partial_phi_over_partial_n1 * partial_n1_over_partial_rho &
                                            + partial_phi_over_partial_n2 * partial_n2_over_partial_rho &
                                            + partial_phi_over_partial_n3 * partial_n3_over_partial_rho )

      write ( * , * ) 'chemical potential mu_exc0 ( ' , species , ' ) = ' , muexc_0_multispec ( species )



      ! compute reference bulk grand-potential Omega(rho = rho_0) !! Do not forget the solver minimizes Omega[rho]-Omega[rho_0] = Fsolvatation

      if ( hs_functional ( 1 : 2 ) == 'PY' ) then

        Fexc_0_multispec ( species ) = kBT * ( - n0 * log ( 1.0_dp - n3 )                            &
                                               + n1 * n2 / ( 1.0_dp - n3 )                           &
                                               + n2 ** 3 / ( 24.0_dp * pi * ( 1.0_dp - n3 ) ** 2 ) )

      else if ( hs_functional ( 1 : 2 ) == 'CS' .or. hs_functional ( 1 : 4 ) == 'MCSL' ) then

        Fexc_0_multispec ( species ) = kBT * ( ( ( 1.0_dp / ( 36.0_dp * pi ) ) * n2 ** 3 / n3 ** 2 - n0 ) * log ( 1.0_dp - n3 ) &
                                             + n1 * n2 / ( 1.0_dp - n3 )                                                        &
                                             + ( 1.0_dp / ( 36.0_dp * pi ) ) * n2 ** 3 / ( ( 1.0_dp - n3 ) ** 2 * n3 )          )

      end if

      ! integration factors

      Fexc_0_multispec ( species ) = Fexc_0_multispec ( species ) * Lx * Ly * Lz &
                                   - muexc_0_multispec ( species) * Lx * Ly * Lz * n_0_multispec ( species )  ! integration factor

      write ( * , * ) 'Fexc0 ( ' , species , ' ) = ' , Fexc_0_multispec ( species )

    end do


  end subroutine excess_chemical_potential_and_reference_bulk_grand_potential
























!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine computes the density independant weight functions as defined by Kierlik and Rosinberg in 1990
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

  subroutine compute_hard_sphere_weight_functions_k
  
  use precision_kinds , only : dp , i2b
  
  use constants , only : FourPi
  
  use system , only : nfft1 , nfft2 , nfft3 , nb_species , radius , &
                      weight_function_3_k , weight_function_2_k , weight_function_1_k , weight_function_0_k
  
  use fft , only : norm_k
  
  
  implicit none
  
  
  real ( kind = dp ) :: kR , FourPiR , sinkR , coskR ! dummy for speeding up
  
  integer ( kind = i2b ) :: species ! dummy between 1 and nb_species
  
  integer ( kind = i2b ) :: l , m , n ! dummy for loops
  
  real ( kind = dp ) :: norm_k_local ! dummy local variable
  
  
  
  ! Allocate the weight functions.
  
  allocate ( weight_function_3_k ( nfft1 / 2 + 1 , nfft2 , nfft3 , nb_species ) )
  
  allocate ( weight_function_2_k ( nfft1 / 2 + 1 , nfft2 , nfft3 , nb_species ) )
  
  allocate ( weight_function_1_k ( nfft1 / 2 + 1 , nfft2 , nfft3 , nb_species ) )
  
  allocate ( weight_function_0_k ( nfft1 / 2 + 1 , nfft2 , nfft3 , nb_species ) )
  
  ! density weights for hard spheres are known analyticaly
  ! they only depends on fundamental measures of hard spheres
  ! here is the Kierlik and Rosinberg FMT : 4 scalar weight function by species
  
  !> Warn user
  
  write (*,*)  '>>> Compute hard sphere weight functions as defined by Kierlik and Rosinberg, PRA1990'
  
  ! For each species (ie each radius), compute weight functions in k space
  
  do species = 1 , nb_species
  
    FourPiR = FourPi * radius ( species )
  
    do n = 1 , nfft3
  
      do m = 1 , nfft2
  
        do l = 1 , nfft1 / 2 + 1
  
          norm_k_local = norm_k ( l , m , n )
  
          if ( norm_k_local /= 0.0_dp ) then
  
            kR = norm_k_local * radius ( species )
  
            sinkR = sin ( kR )
  
            coskR = cos ( kR )
  
            weight_function_3_k ( l , m , n , species ) = FourPi * ( sinkR - kR * coskR ) / ( norm_k_local ** 3 )
  
            weight_function_2_k ( l , m , n , species ) = FourPiR * sinkR / norm_k_local
  
            weight_function_1_k ( l , m , n , species ) = ( sinkR + kR * coskR ) / ( 2.0_dp * norm_k_local )
  
            weight_function_0_k ( l , m , n , species ) = coskR + 0.5_dp * kR * sinkR
  
          else if ( norm_k_local == 0.0_dp ) then
  
            weight_function_3_k ( l , m , n , species ) = FourPi / 3.0_dp * radius ( species ) ** 3 ! volume
  
            weight_function_2_k ( l , m , n , species ) = FourPi * radius ( species ) ** 2 ! surface area
  
            weight_function_1_k ( l , m , n , species ) = radius ( species ) ! radius
  
            weight_function_0_k ( l , m , n , species ) = 1.0_dp ! unity
  
          end if
  
        end do ! n nfft3 end loop over k-vectors
  
      end do ! m nfft2
  
    end do ! l nfft1
  
  end do ! species
  
  
  
  end subroutine compute_hard_sphere_weight_functions_k


















end subroutine compute_hard_spheres_parameters



















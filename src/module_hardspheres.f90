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
    
    PUBLIC :: packfrac, populate_weight_functions_in_Fourier_space
    
    
    
    CONTAINS

!===================================================================================================================================    
    
        ! packing fraction : eta = 4/3 * pi * R^3 * solvent density of the constituant ( /= total solvent density)
        PURE REAL(dp) FUNCTION packfrac (n,R)
            REAL(dp), INTENT(IN) :: n ! density of the constituant (/= total solvant density)
            REAL(dp), INTENT(IN) :: R ! radius of the constituant
            packfrac = 4._dp/3._dp * ACOS(-1._dp) * R**3 * n
        END FUNCTION packfrac

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



END MODULE hardspheres

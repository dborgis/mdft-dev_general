MODULE hardspheres

    USE precision_kinds  ,ONLY: dp
    IMPLICIT NONE
    
    COMPLEX(dp), ALLOCATABLE , DIMENSION (:,:,:,:,:) :: weightfun_k
    
    TYPE type_hs
        REAL(dp) :: R          ! radius
        REAL(dp) :: excchempot ! excess chemical potential
        REAL(dp) :: Fexc0
        REAL(dp) :: pf         ! packing fraction
    END TYPE type_hs
    TYPE (type_hs), ALLOCATABLE :: hs(:) ! one per constituant
    
    
    
    
    CONTAINS
    
    
        ! packing fraction : eta = 4/3 * pi * R^3 * solvent density of the constituant ( /= total solvent density)
        PURE REAL(dp) FUNCTION packfrac (n,R)
            REAL(dp), INTENT(IN) :: n ! density of the constituant (/= total solvant density)
            REAL(dp), INTENT(IN) :: R ! radius of the constituant
            packfrac = 4._dp/3._dp * ACOS(-1._dp) * R**3 * n
        END FUNCTION packfrac










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
    
            USE precision_kinds ,ONLY: dp , i2b
            USE constants       ,ONLY: FourPi, zeroC
            USE system          ,ONLY: nb_species, spaceGrid
            USE fft             ,ONLY: norm_k
            IMPLICIT NONE
    
            REAL(dp) :: kR , FourPiR , sinkR , coskR, norm_k_local
            INTEGER(i2b) :: l,m,n,s,nfft(3)
            nfft = spaceGrid%n_nodes
            
            ! Weight functions of a fluid of hard spheres are known analyticaly. They are easily defined in Fourier space.
            ! They depend on fundamental measures of the fluid.
            ! Here is the Kierlik and Rosinberg'FMT : 4 scalar weight functions by species. Ref: Kierlik and Rosinberg, PRA 1990
            ALLOCATE ( weightfun_k (nfft(1)/2+1, nfft(2), nfft(3), nb_species, 0:3 ) ,SOURCE=zeroC) !0:3 for Kierlik-Rosinberg
            DO CONCURRENT ( s=1:nb_species, l=1:nfft(1)/2+1, m=1:nfft(2), n=1:nfft(3) )
                FourPiR = FourPi * hs(s)%R
                norm_k_local = norm_k (l,m,n)
                IF ( norm_k_local /= 0.0_dp ) THEN
                    kR = norm_k_local * hs(s)%R
                    sinkR = SIN(kR)
                    coskR = COS(kR)
                    weightfun_k(l,m,n,s,3) = FourPi * ( sinkR - kR * coskR ) / ( norm_k_local ** 3 )
                    weightfun_k(l,m,n,s,2) = FourPiR * sinkR / norm_k_local
                    weightfun_k(l,m,n,s,1) = ( sinkR + kR * coskR ) / ( 2.0_dp * norm_k_local )
                    weightfun_k(l,m,n,s,0) = coskR + 0.5_dp * kR * sinkR
                ELSE
                    weightfun_k(l,m,n,s,3) = FourPi / 3.0_dp * hs(s)%R** 3 ! volume
                    weightfun_k(l,m,n,s,2) = FourPi * hs(s)%R** 2 ! surface area
                    weightfun_k(l,m,n,s,1) = hs(s)%R ! radius
                    weightfun_k(l,m,n,s,0) = 1.0_dp ! unity
                END IF
            END DO
        END SUBROUTINE populate_weight_functions_in_Fourier_space



END MODULE hardspheres

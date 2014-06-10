SUBROUTINE energy_hard_sphere_fmt (Fint)

    USE precision_kinds  ,ONLY: dp , i2b
    USE system           ,ONLY: kBT , nb_species , n_0_multispec , mole_fraction , rho_0_multispec, spaceGrid
    USE quadrature       ,ONLY: molRotSymOrder , angGrid, molRotGrid
    USE minimizer        ,ONLY: cg_vect , FF , dF
    USE constants        ,ONLY: pi , FourPi , twopi, zeroC
    USE fft              ,ONLY: fftw3
    USE input            ,ONLY: verbose, input_char, input_dp
    USE hardspheres      ,ONLY: weightfun_k, hs
    
    IMPLICIT NONE
    
    REAL(dp), INTENT(OUT) :: Fint ! Internal part of free energy
    INTEGER(i2b) :: icg , i , j , k , o , p,s,nfft1,nfft2,nfft3
    REAL(dp) :: Nk ! Total number of k points = nfft1*nfft2*nfft3
    REAL(dp) :: local_density, psi, deltav
    REAL(dp), ALLOCATABLE, DIMENSION(:) :: nb_molecules ! temp before full implementation of multispecies
    REAL(dp), ALLOCATABLE, DIMENSION(:,:,:,:) :: rho ! density per angle (recall : rho_0 = n_0 / 4pi ) ! x y z nb_species
    COMPLEX(dp), ALLOCATABLE , DIMENSION(:,:,:,:) :: rho_k !> @var Density in k space
    REAL(dp) :: time0 , time1 ! time stamps
    REAL(dp) :: w0,w1,w2,w3,omw3
    REAL(dp) :: F_HS ! Perkus Yevick (or Carnahan Starling) expression for excess energy
    REAL(dp)   , ALLOCATABLE , DIMENSION(:,:,:,:) :: wd ! weighted density at node i,j,k, for index 0:3
    REAL(dp)   , ALLOCATABLE , DIMENSION(:,:,:) :: one_min_wd_3 ! dummy for 1 - weighted density 3
    REAL(dp)   , ALLOCATABLE , DIMENSION(:,:,:) :: dFHS_0 , dFHS_1 , dFHS_2 , dFHS_3
    COMPLEX(dp), ALLOCATABLE , DIMENSION(:,:,:) :: dFHS_0_k , dFHS_1_k , dFHS_2_k , dFHS_3_k
    REAL(dp)   , ALLOCATABLE , DIMENSION(:,:,:,:) :: dFHS
    COMPLEX(dp), ALLOCATABLE , DIMENSION(:,:,:,:) :: dFHS_k
    COMPLEX(dp), ALLOCATABLE , DIMENSION(:,:,:,:) :: dFex_k ! gradient in k space
    REAL(dp)   , ALLOCATABLE , DIMENSION(:,:,:,:) :: dFex ! gradient in real space
    CHARACTER(2) :: hs_functional ! hard sphere functional = PY for Percus-Yevick or CS for Carnahan-Starling
    REAL(dp), PARAMETER :: inv8pi  = 1.0_dp/(  8.0_dp * pi)
    REAL(dp), PARAMETER :: inv12pi = 1.0_dp/( 12.0_dp * pi)
    REAL(dp), PARAMETER :: inv18pi = 1.0_dp/( 18.0_dp * pi)
    REAL(dp), PARAMETER :: inv24pi = 1.0_dp/( 24.0_dp * pi)
    REAL(dp), PARAMETER :: inv36pi = 1.0_dp/( 36.0_dp * pi)
    REAL(dp) :: GCdmu
    
    CALL CPU_TIME ( time0 ) ! init timer

    nfft1 = spaceGrid%n_nodes(1)
    nfft2 = spaceGrid%n_nodes(2)
    nfft3 = spaceGrid%n_nodes(3)
    Nk = REAL(PRODUCT(spaceGrid%n_nodes),dp) ! total number of k points needed for inverse fft normalization
    deltav = spaceGrid%dv

    GCdmu = input_dp('impose_some_deltamu')

    ALLOCATE ( rho (nfft1,nfft2,nfft3,nb_species) ,SOURCE=0._dp)
    icg = 0
    DO s = 1 , nb_species
        DO i = 1 , nfft1
            DO j = 1 , nfft2
                DO k = 1 , nfft3
                    local_density = 0.0_dp
                    DO o = 1, angGrid%n_angles
                        DO p=1, molRotGrid%n_angles
                            icg = icg + 1
                            local_density = local_density + angGrid%weight (o) * cg_vect (icg) ** 2*molRotGrid%weight(p)
                        END DO          
                    END DO
                    ! correct by 8*pi²/n as the integral over all orientations o and psi is 4pi and 2pi/n
                    local_density = local_density*REAL(molRotSymOrder, dp) / (fourpi*twopi)
                    ! at the same time integrate rho in order to count the total number of implicit molecules. here we forget the integration factor = n_0 * deltav
                    rho(i,j,k,s) = local_density
                END DO
            END DO
        END DO
    END DO

    ! total number of molecules of each species
    ALLOCATE ( nb_molecules ( nb_species ) ,SOURCE=0._dp)
    DO CONCURRENT ( s=1:nb_species )
        nb_molecules(s) = SUM( rho(:,:,:,s ) ) * n_0_multispec(s) * mole_fraction(s) * deltav
    END DO
!~     IF (verbose) THEN     ! tell user about the number of molecule of each species in the supercell
!~         DO s = 1 , nb_species
!~             PRINT*,'There are ',nb_molecules(s),' molecules of type',s
!~         END DO
!~     END IF

    ! fourier transform the density rho => rho_k
    ALLOCATE ( rho_k (nfft1/2+1,nfft2,nfft3,nb_species) ,SOURCE=zeroC)
    DO s = 1 , nb_species
        fftw3%in_forward = rho(:,:,:,s)
        CALL dfftw_execute ( fftw3%plan_forward )
        rho_k(:,:,:,s) = fftw3%out_forward * n_0_multispec(s) * mole_fraction(s)
    END DO
    DEALLOCATE ( rho )

    ! inverse fourier transform the weighted densities
    ! allocate the arrays for receiving the FFT-1
    ALLOCATE ( wd (nfft1,nfft2,nfft3,0:3) ,SOURCE=0._dp) ! 0:3 for Kierliek Rosinberg
    DO s= 1, nb_species
        DO i= LBOUND(wd,4) , UBOUND(wd,4)
            fftw3%in_backward = weightfun_k(:,:,:,s,i) * rho_k(:,:,:,s)
            CALL dfftw_execute ( fftw3%plan_backward )
            wd(:,:,:,i) = wd(:,:,:,i) + fftw3%out_backward / Nk
        END DO
    END DO
    
    ! check if the hard sphere functional is Percus-Yevick or Carnahan-Starling
    ! Get the free energy functional that should be used. For now Percus Yevick and Carnahan Starling only. May be expanded.
    hs_functional = input_char('hs_functional')
    
    ! compute free intrinsic energy
    Fint = 0.0_dp
    DO k = 1 , nfft3 ! please pay attention to inner / outer loop.
        DO j = 1 , nfft2
            DO i = 1 , nfft1

                w0 = wd(i,j,k,0)
                w1 = wd(i,j,k,1)
                w2 = wd(i,j,k,2)
                w3 = wd(i,j,k,3)

                omw3 = 1.0_dp - w3 ! nowhere should w3 be lower or equal than 1
                IF ( omw3 <= 0.0_dp ) THEN
                    CALL error_message_energy_hard_sphere_fmt ( i , j , k , w0 , w1 , w2 , w3 )
                    CALL process_output ! process output so that we know where we are !
                END IF

                IF ( hs_functional == 'PY' .OR. hs_functional == 'py' ) THEN
                    F_HS = -w0*LOG(omw3) &
                           +w1*w2/omw3 &
                           +inv24pi*w2**3/omw3**2
                ELSE IF ( hs_functional == 'CS' .OR. hs_functional == 'cs' ) THEN
                    F_HS = ( inv36pi * w2 ** 3 / w3 ** 2 - w0 ) * LOG(omw3) &
                           + w1 * w2 / omw3 &
                           + inv36pi * w2 ** 3 / ( omw3 ** 2 * w3 )
                END IF
  
                Fint = Fint + F_HS
                
            END DO
        END DO
    END DO

    IF( SIZE(nb_molecules)/=1 ) STOP "HERE GCdmu is for one component only"
    Fint = Fint * kBT * DeltaV -SUM(hs%excchempot * nb_molecules) -SUM( hs%Fexc0 * mole_fraction ) -GCdmu*nb_molecules(1)
    FF = FF + Fint

    ! gradients
    ! dFHS_i and weighted_density_j are arrays of dimension (nfft1,nfft2,nfft3)
    ! one_min_weighted_density_3 is dummy for speeding up and simplifying while using unnecessary memory.
    ALLOCATE ( one_min_wd_3 ( nfft1 , nfft2 , nfft3 ) ,SOURCE=1.0_dp-wd(:,:,:,3))
    
    ! excess free energy per atom derived wrt weighted density 0 (eq. A5 in Sears2003)
    ALLOCATE ( dFHS(nfft1,nfft2,nfft3,0:3) ,SOURCE=0._dp)
    
    ! Perkus Yevick
    IF ( hs_functional == 'PY' .or. hs_functional == 'py' ) THEN
        dFHS(:,:,:,0) = - log ( one_min_wd_3 )
        dFHS(:,:,:,1) = wd(:,:,:,2) / one_min_wd_3
        dFHS(:,:,:,2) = wd(:,:,:,1) / one_min_wd_3 + inv8pi * dFHS(:,:,:,1) ** 2
        dFHS(:,:,:,3) = ( wd(:,:,:,0) + wd(:,:,:,1) * dFHS(:,:,:,1) )/one_min_wd_3  + inv12pi*dFHS(:,:,:,1)**3
    
    ! Carnahan Starling
    ELSE IF ( hs_functional == 'CS' .or. hs_functional == 'cs' ) THEN
    
        dFHS(:,:,:,0) = - log ( one_min_wd_3 )
    
        dFHS(:,:,:,1) = wd(:,:,:,2) / one_min_wd_3
    
        dFHS(:,:,:,2) = - inv12pi * ( wd(:,:,:,2) / wd(:,:,:,3) ) ** 2 * dFHS(:,:,:,0) &
                + wd(:,:,:,1) / one_min_wd_3 + inv12pi * dFHS(:,:,:,1) ** 2 / wd(:,:,:,3)
    
        dFHS(:,:,:,3) = inv18pi * dFHS(:,:,:,0) * ( wd(:,:,:,2) / wd(:,:,:,3) ) ** 3 &
                - ( inv36pi * wd(:,:,:,2) ** 3 / wd(:,:,:,3) ** 2 - wd(:,:,:,0) ) /&
                                    one_min_wd_3 &
                + wd(:,:,:,1) * dFHS(:,:,:,1) / one_min_wd_3 + inv36pi * &
                wd(:,:,:,2) ** 3 / wd(:,:,:,3) ** 2 * &
                ( 3.0_dp * wd(:,:,:,3) - 1.0_dp ) / one_min_wd_3 ** 3
    END IF
    
    ! deallocate weighted_densities
    DEALLOCATE ( wd )
    DEALLOCATE ( one_min_wd_3 )
    
    ! compute gradients in k space
    ALLOCATE ( dFHS_k(nfft1/2+1,nfft2,nfft3,0:3) ,SOURCE=zeroC)
    
    ! FFT dFHS for computing convolution
    DO i= LBOUND(dFHS,4) , UBOUND(dFHS,4)
        fftw3%in_forward = dFHS(:,:,:,i)
        CALL dfftw_execute ( fftw3%plan_forward )
        dFHS_k(:,:,:,i) = fftw3%out_forward
    END DO
    
    ! deallocate useless
    DEALLOCATE ( dFHS )
    
    ! compute final gradient in k-space
    ALLOCATE ( dFex_k (nfft1/2+1,nfft2,nfft3,nb_species) ,SOURCE=zeroC) ! FFT of dFex (in fact dFex_k is known from which FFT-1 gives dFex)
    DO CONCURRENT (s=1:nb_species)
        dFex_k(:,:,:,s) = SUM( dFHS_k*weightfun_k(:,:,:,s,:) ,4)
    END DO
    DEALLOCATE ( dFHS_k )
    
    ! inverse fourier transform gradient
    ALLOCATE ( dFex (nfft1,nfft2,nfft3,nb_species) ,SOURCE=0._dp)
    DO s = 1, nb_species
        fftw3%in_backward = dFex_k(:,:,:,s)
        CALL dfftw_execute ( fftw3%plan_backward )
        dFex(:,:,:,s) = fftw3%out_backward / Nk
    END DO
    DEALLOCATE ( dFex_k )
    
    ! transfer in rank 1 vector dF
    icg = 0
    DO s = 1 , nb_species
        DO i = 1, nfft1
            DO j = 1, nfft2
                DO k = 1, nfft3
                    DO o = 1, angGrid%n_angles
                        DO p=1, molRotGrid%n_angles
                            icg = icg + 1
                            psi = cg_vect ( icg )
    ! AHAAAAAAAAAAAAATTENTION ICI LE 12 SEPTEMBRE 2011 JE ME RENDS COMPTE QU iL Y A PEUT ETRE UN FACTEUR WEIGHT(o) QUI MANQUE, CA DEPEND DE LA DEF DE N0 par RAPPORT à rho_0
                            dF ( icg ) = dF ( icg ) &
                            + 2.0_dp * psi * rho_0_multispec ( s ) * deltav * ( kBT * dFex ( i , j , k , s ) - hs(s)%excchempot )&
                            *angGrid%weight(o)*molRotGrid%weight(p)
    ! ATTENTION J'AI MULTIPLIE PAR WEIGHT(O) QD PASSAGE A angGrid%n_angles /=1 LE 18 JUILLET 2011
                        END DO !p
                    END DO !o
                END DO !i
            END DO !j
        END DO !k
    END DO !species
    DEALLOCATE (dFex) ! deallocate useless gradient
    
    CALL CPU_TIME ( time1 ) !stop timer
    

    CONTAINS

    !===============================================================================================================================
    ! this SUBROUTINE prints error message related to SUBROUTINE excess_cs_hard_sphere
    ! it may stop program execution depending on the error.
    !===============================================================================================================================
    SUBROUTINE error_message_energy_hard_sphere_fmt ( i , j , k , w0 , w1 , w2 , w3 )
        USE precision_kinds ,ONLY: i2b, dp
        IMPLICIT NONE
        INTEGER(i2b), INTENT(IN) :: i,j,k
        REAL(dp), INTENT(IN) :: w0,w1,w2,w3
        PRINT*,'i , j , k = ',i,j,k
        PRINT*,'log (1-w3<=0) in energy_hard_sphere_fmt.f90. Critical stop'
        PRINT*,'w0, w1, w2, w3 = ' ,w0,w1,w2,w3
        STOP
    END SUBROUTINE error_message_energy_hard_sphere_fmt

END SUBROUTINE energy_hard_sphere_fmt

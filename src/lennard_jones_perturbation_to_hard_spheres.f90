! this SUBROUTINE computes the lennard jones perturbative contribution to the hard sphere Helmotz energy
! 201109121332 creation by Maximilien Levesque
! 201109151545 added the calculation of the perturbation potential
SUBROUTINE lennard_jones_perturbation_to_hard_spheres

    USE precision_kinds ,ONLY: dp,i2b
    USE system          ,ONLY: nfft1,nfft2,nfft3,Lx,Ly,Lz,n_0,sig_solv,eps_solv,v_perturbation_k,spaceGrid,nb_species
    USE quadrature      ,ONLY: angGrid
    USE minimizer       ,ONLY: cg_vect,dF,FF
    USE constants       ,ONLY: fourpi,twopi,zeroC
    USE fft             ,ONLY: fftw3, norm_k

    IMPLICIT NONE
    
    REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) :: rho_n ! local density
    COMPLEX(dp), ALLOCATABLE, DIMENSION(:,:,:) :: rho_k ! fourier transformed rho_n
    REAL(dp) :: local_density ! dummy for rho_n
    INTEGER(i2b) :: icg , i , j , k , o , l , m , n , m1 , m2 , m3 ! dummy
    INTEGER(i2b) :: nf1 , nf2 , nf3 ! dummy for nfft1 /2 , nfft2 /2, nfft3 /2
    REAL(dp) :: Nk ! total number of k points
    REAL(dp) :: kx2 , ky2 , kz2 
    COMPLEX(dp), ALLOCATABLE, DIMENSION(:,:,:) :: vk ! fourier transform of the lennard jones perturbation (WCA)
    REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) :: v_perturbation_r ! vk in REAL space
    REAL(dp) :: Fperturbation ! what we want : the perturbative contribution of the lennard jones attractive potential tail
    REAL(dp) :: time0 , time1 ! timers
    REAL(dp) :: twopiolx , twopioly , twopiolz ! dummy for speeding up loops
    REAL(dp) :: potential ! dFp / drho at rho=rho_0 in order the grand potential to be zero at rho = rho_0
    REAL(dp) :: nb_molecule ! total number of hard spheres ie integral of density over all space

    IF (nb_species/=1) STOP "When dealing with LJ as perturbation in lennard_jones_..._spheres.f90, only nb_species=1 is ok."
    
    ! init
    CALL cpu_time ( time0 )
    nf1 = nfft1 / 2
    nf2 = nfft2 / 2
    nf3 = nfft3 / 2

    ! put result from last minimization step as density from which one computes energy and gradients
    ALLOCATE( rho_n (nfft1,nfft2,nfft3 ) ,SOURCE=0._dp)
    icg = 0
    nb_molecule = 0.0_dp
    DO i = 1 , nfft1
        DO j = 1 , nfft2
            DO k = 1 , nfft3
                local_density = 0.0_dp
                DO o = 1, angGrid%n_angles ! angGrid%n_angles=1
                    icg = icg + 1
                    local_density = local_density + angGrid%weight (o) * cg_vect (icg) ** 2
                END DO
                local_density = local_density / fourpi ! correct by fourpi as the integral over all orientations o is 4pi
                ! at the same time integrate rho_n in order to count the total number of implicit molecules. here we forget the integration factor = n_0 * deltav
                rho_n ( i , j , k ) = local_density
                nb_molecule = nb_molecule + local_density * n_0 * spaceGrid%dv
            END DO
        END DO
    END DO

    ! total number of k points needed for inverse fft normalization
    Nk = REAL ( nfft1 * nfft2 * nfft3 , dp )
    ! fourier transform the density rho_n => rho_k
    fftw3%in_forward = rho_n
    CALL dfftw_execute ( fftw3%plan_forward )
    ALLOCATE( rho_k ( nfft1 / 2 + 1 , nfft2 , nfft3 ) ,SOURCE=zeroC)
    rho_k = fftw3%out_forward
    
    ! compute lennard jones perturbation in k space
    ! if v_perturbation_k doesn't exist, then compute it and put it in Vk. 
    IF ( .NOT. ALLOCATED( v_perturbation_k ) ) THEN
        ALLOCATE ( v_perturbation_k (nfft1/2+1,nfft2,nfft3) ,SOURCE=zeroC )
        DO n=1,nfft3
            DO m=1,nfft2
                DO l=1,nfft1/2+1
                    v_perturbation_k(l,m,n) = vlj_wca_k ( norm_k(l,m,n), sig_solv(1), eps_solv(1) )
                END DO
            END DO
        END DO
    END IF

    ! put the backup in what we use (perhaps redondant)
    ALLOCATE( Vk ( nfft1 / 2 + 1 , nfft2 , nfft3 ) ,SOURCE=v_perturbation_k) ! Vk is the one we use in this routine. It may be saved in v_perturbation_k in order not to compute it each time.
    ! once equation written, dFp / drho at rho = rho_0 is shown to be equal to rho_0 * Vk(k=0)
    potential = n_0 * REAL ( Vk ( 1 , 1 , 1 ) )

    ! FFT-1 of perturbation
    fftw3%in_backward = Vk * rho_k

    DEALLOCATE( Vk )
    DEALLOCATE( V_perturbation_k )
    DEALLOCATE( rho_k )

    CALL dfftw_execute (fftw3%plan_backward)
    ALLOCATE( v_perturbation_r ( nfft1 , nfft2 , nfft3 ) ,SOURCE=(fftw3%out_backward/Nk))
    
    ! Compute the free energy due to the perturbation
    Fperturbation = 0.0_dp
    DO i = 1 , nfft1
        DO j = 1 , nfft2
            DO k = 1 , nfft3
                DO o = 1, angGrid%n_angles ! angGrid%n_angles=1
                    Fperturbation = Fperturbation + rho_n ( i , j , k ) * v_perturbation_r ( i , j , k )
                END DO
            END DO
        END DO
    END DO
    Fperturbation = Fperturbation * 0.5_dp * spaceGrid%dv ! normalization
    DEALLOCATE( rho_n )

    ! add perturbation energy to total energy
    PRINT*,'nb_molecule = ',nb_molecule
    PRINT*,'mu_p = ', potential
    PRINT*,'- nb_molecule * potential = ', - nb_molecule * potential
    FF = FF + Fperturbation - nb_molecule * potential

    !===============================================================================================================================
    ! gradient
    !===============================================================================================================================
    icg = 0
    DO i = 1 , nfft1
        DO j = 1 , nfft2
            DO k = 1 , nfft3
                DO o = 1 , angGrid%n_angles
                    icg = icg + 1
                    dF(icg) = dF(icg) + 2.0_dp * cg_vect(icg) * spaceGrid%dv * v_perturbation_r(i,j,k)
                END DO
            END DO
        END DO
    END DO
    DEALLOCATE( v_perturbation_r )

    CALL cpu_time ( time1 )
    !IF(verbose) PRINT*,'Fperturbati = ' , Fperturbation , 'computed in (sec)' , time1 - time0



    CONTAINS

    
        !===========================================================================================================================
        ! Here we compute Uperturbation in kspace. it's an integration we DO numericaly
        FUNCTION vlj_wca_k ( k , sigma_lj , epsilon_lj )
        !===========================================================================================================================
        
            USE precision_kinds ,ONLY: dp, i2b
            USE constants       ,ONLY: fourpi, zeroC
            
            IMPLICIT NONE
            
            COMPLEX(dp):: vlj_wca_k ! which computes the reciprocal value of the potential 'vk'
            REAL(dp), INTENT(IN) :: k ! one gives the k point 'k' to eat to the routine
            REAL(dp), INTENT(IN) :: sigma_lj , epsilon_lj ! lennard jones parameters in Angstroms and KJ/mol
            REAL(dp):: cutoff ! cut off under which U(r) = constant = -epsilon_lj and after which U(r)=Vlj(r)
            REAL(dp):: dx ! width of the integration step
            REAL(dp):: borne_sup , borne_inf !sup and inf limits of the integration
            INTEGER(i2b), PARAMETER :: nstep = 1000 ! nb of integration steps ! TODO CHECK EFFECT
            REAL(dp):: ri !integrants
            REAL(dp):: sigmaori6 !dummy for (sigma_lj/ri)**6
            INTEGER(i2b):: i !dummy for loop
            
            !the integration is theoreticaly between 2^(1/6)*sigma and infinity
            !we do it between sigma_lj and N*sigma for a beginning
            !perhaps 10*sigma is more than enough
            !we use 10**3 steps for the integration.
            
            borne_sup = 10.0_dp * sigma_lj ! TODO 40 is better and converged. check later and do analyticaly
            borne_inf = 0.0_dp
            
            ! TODO here one integrates over the whole range of r. one could calculate the first part (U=-epsilon_lj) analyticaly and thus speed up everything by a factor of d_wca/borne_sup
            dx = ( borne_sup - borne_inf ) / REAL(nstep,dp)
            
            cutoff = 2.0_dp**(1.0_dp/6.0_dp)*sigma_lj ! compute the value after which U = Vlj, which is the value of x for which Vlj is minimum, thus 2^(1/6)sigma

            ! init vk
            vlj_wca_k = zeroC
            
            ! integrate v(r) over all r in R^3
            DO i = 1 , nstep
                ri = borne_inf + REAL ( i - 1 , dp ) * dx
                IF ( ri < cutoff ) then ! IF under cutoff
                    IF ( k > 0.0_dp ) then
                        vlj_wca_k = vlj_wca_k + CMPLX(-epsilon_lj*ri*SIN(k*ri)/k , 0.0_dp )
                    ELSE ! if k = 0 then lim sin(k*r)/k = r
                        vlj_wca_k = vlj_wca_k + CMPLX(-epsilon_lj*ri**2 , 0.0_dp )
                    END IF
                ELSE
                    sigmaori6 = (sigma_lj/ri)**6
                    IF ( k > 0.0_dp ) then ! general case
                        vlj_wca_k = vlj_wca_k + CMPLX(4.0_dp*epsilon_lj*ri*(sigmaori6**2-sigmaori6)*SIN(k*ri)/k  ,0.0_dp)
                    ELSE ! if k = 0 then lim sin(k*r)/k = r 
                        vlj_wca_k = vlj_wca_k + CMPLX(4.0_dp*epsilon_lj*ri*(sigmaori6**2-sigmaori6)*ri , 0.0_dp )
                    END IF
                END IF
            END DO
            vlj_wca_k = vlj_wca_k * CMPLX( dx*fourpi , 0.0_dp ) ! integration factors

        END FUNCTION vlj_wca_k
        !===========================================================================================================================

END SUBROUTINE lennard_jones_perturbation_to_hard_spheres

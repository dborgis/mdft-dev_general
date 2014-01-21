!> Gets the final density from the last minimizer step.
SUBROUTINE get_final_density ( neq )

    USE precision_kinds, ONLY: dp , i2b
    USE system,          ONLY: nb_species, spaceGrid
    USE constants,       ONLY: twopi
    USE minimizer,       ONLY: CG_vect
    USE quadrature,      ONLY: molRotSymOrder, angGrid, molRotGrid
    USE fft,             ONLY: fftw3 , timesExpPrefactork2
    USE input,           ONLY: verbose
    
    IMPLICIT NONE

    INTEGER(i2b) :: i,j,k,omega,icg,s,p,nfft1,nfft2,nfft3
    REAL(dp), INTENT(OUT) :: neq (spaceGrid%n_nodes(1),spaceGrid%n_nodes(2),spaceGrid%n_nodes(3),nb_species) ! equilibrium density(position)
    REAL(dp) :: rho,local_density,Nk
    REAL(dp), DIMENSION(3) :: dl
    REAL(dp), PARAMETER :: gaussianWidth = 0._dp
    COMPLEX(dp), ALLOCATABLE, DIMENSION(:,:,:,:) :: rho_k 
    
    nfft1 = spaceGrid%n_nodes(1); nfft2 = spaceGrid%n_nodes(2); nfft3 = spaceGrid%n_nodes(3)
    dl = spaceGrid%dl
    Nk = real( sum( spaceGrid%n_nodes ), dp)

    icg = 0
    DO s =1,nb_species
        DO i =1,nfft1
            DO j =1,nfft2
                DO k =1,nfft3
                    local_density = 0.0_dp
                    DO omega =1,angGrid%n_angles
                        DO p =1,molRotGrid%n_angles
                            icg = icg + 1
                            rho = cg_vect ( icg ) ** 2 / (twopi**2*2.0_dp/molRotSymOrder)
                            local_density = local_density + angGrid%weight ( omega ) *molRotGrid%weight(p)* rho ! integral of rho over all orientations ie 'n'
                        END DO
                    END DO
                    neq (i,j,k,s) = local_density
                END DO
            END DO
        END DO
    END DO

!~     IF ( gaussianWidth /= 0._dp ) THEN !convolute with a gaussian
!~         ALLOCATE ( rho_k (nfft1/2+1, nfft2, nfft3, nb_species) )
!~         DO s = 1, nb_species
!~             fftw3%in_forward = neq ( : , : , : , s )
!~             CALL dfftw_execute ( fftw3%plan_forward )
!~             rho_k (:,:,:,s) = timesExpPrefactork2 ( fftw3%out_forward, gaussianWidth**2/2.0_dp )
!~             fftw3%in_backward = rho_k (:,:,:,s)
!~             DEALLOCATE( rho_k )
!~             CALL dfftw_execute ( fftw3%plan_backward )
!~             neq (:,:,:,s ) = fftw3%out_backward/Nk 
!~         END DO
!~     END IF

END SUBROUTINE get_final_density

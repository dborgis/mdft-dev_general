!> Gets the final density from the last minimizer step.
SUBROUTINE get_final_density ( neq )

    USE precision_kinds, ONLY: dp , i2b
    USE system, ONLY : nb_species, spaceGrid
    USE constants, ONLY : twopi
    USE minimizer, ONLY: CG_vect
    USE quadrature, ONLY : sym_order, angGrid, molRotGrid
    USE fft, ONLY: fftw3 , timesExpPrefactork2
    USE input, ONLY: verbose
    
    IMPLICIT NONE

    REAL(dp), INTENT(OUT) :: neq (spaceGrid%n_nodes(1),spaceGrid%n_nodes(2),spaceGrid%n_nodes(3),nb_species) ! equilibrium density(position)
    INTEGER(i2b) :: i, j, k, omega, icg, species, p
    REAL(dp) :: rho_over_fourpi, local_density
    COMPLEX(dp), ALLOCATABLE, DIMENSION(:,:,:,:) :: rho_k 
    REAL(dp) :: Nk ! Total number of k points = nfft1*nfft2*nfft3
    INTEGER(i2b) :: nfft1, nfft2, nfft3
    REAL(dp), DIMENSION(3) :: dl
    REAL(dp), PARAMETER :: gaussianWidth = 0._dp
    
    nfft1 = spaceGrid%n_nodes(1); nfft2 = spaceGrid%n_nodes(2); nfft3 = spaceGrid%n_nodes(3)
    dl = spaceGrid%dl
    Nk = real( sum( spaceGrid%n_nodes ), dp)
        
    icg = 0
    do species = 1 , nb_species
        do i = 1 , nfft1 ; do j = 1 , nfft2 ; do k = 1 , nfft3
            local_density = 0.0_dp
            do omega = 1 , angGrid%n_angles
                do p=1, molRotGrid%n_angles
                    icg = icg + 1
                    rho_over_fourpi = cg_vect ( icg ) ** 2 / (twopi**2*2.0_dp/sym_order)
                    local_density = local_density + angGrid%weight ( omega ) *molRotGrid%weight(p)* rho_over_fourpi ! integral of rho over all orientations ie 'n'
                END DO
            END DO
            neq (i,j,k,species) = local_density
        END DO; END DO ; END DO
    END DO

    IF ( gaussianWidth /= 0._dp ) THEN !convolute with a gaussian
        ALLOCATE ( rho_k (nfft1/2+1, nfft2, nfft3, nb_species) )
        DO species = 1, nb_species
            fftw3%in_forward = neq ( : , : , : , species )
            CALL dfftw_execute ( fftw3%plan_forward )
            rho_k (:,:,:,species) = timesExpPrefactork2 ( fftw3%out_forward, gaussianWidth**2/2.0_dp )
            fftw3%in_backward = rho_k (:,:,:,species)
            DEALLOCATE( rho_k )
            CALL dfftw_execute ( fftw3%plan_backward )
            neq (:,:,:,species ) = fftw3%out_backward/Nk 
        END DO
    END IF

    IF (verbose) THEN
        open(11,file='output/density_along_x.dat')
            do i=1,nfft1
                write(11,*) (i-1)*dl(1), neq(i,nfft2/2+1,nfft3/2+1,1)
            END DO
        close(11)
    
        open(12,file='output/density_along_y.dat')
            do i=1,nfft2
                write(12,*) (i-1)*dl(2), neq(nfft1/2+1,i,nfft3/2+1,1)
            END DO
        close(12)
    
        open(13,file='output/density_along_z.dat')
            do i=1,nfft3
                write(13,*) (i-1)*dl(3), neq(nfft1/2+1,nfft2/2+1,i,1)
            END DO
        close(13)
    END IF

END SUBROUTINE get_final_density

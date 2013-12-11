! Compute total energy and gradients using direct correlation functions c_s_hs of a hard sphere fluid
SUBROUTINE energy_cs_hard_sphere (Fint)

    USE precision_kinds, only: i2b,dp
    use system, only: nfft1 , nfft2 , nfft3 , Lx , Ly , Lz , c_s_hs , kBT , deltaV , rho_0_multispec , nb_species
    use quadrature, only: sym_order, angGrid, molRotGrid
    USE minimizer, ONLY: cg_vect , FF , dF
    use constants, only: fourpi , pi , twopi
    use fft, only: fftw3, norm_k
    USE dcf, ONLY: nb_k , delta_k
    USE input, ONLY: verbose
    
    IMPLICIT NONE
    integer(i2b) :: i, j, k, l, m, n, o, icg, species,p !> Dummy
    integer(i2b) :: k_index, Nk
    real(dp), INTENT(OUT) :: Fint !> Internal part of the free energy
    real(dp) :: Vint !> Dummy for calculation of Vint
    real(dp) :: fact !> facteur d'integration
    real(dp) :: psi ! Dummy
    real(dp), allocatable, dimension(:,:,:) :: Delta_rho
    complex(dp), allocatable, dimension(:,:,:) :: rho_k , Vpolarization_k
    real(dp) :: time1, time0
    real(dp) , allocatable , dimension(:,:,:) :: Vpolarization
    
    
    call cpu_time ( time0 )
    Nk = nfft1*nfft2*nfft3 ! nombre de points k
    allocate ( Delta_rho ( nfft1 , nfft2 , nfft3 ) )
    Delta_rho = 0.0_dp ! density
    !> Put density of last minimization step in delta_rho
    icg=0
    do i=1,nfft1
    do j=1,nfft2
        do k=1,nfft3
        do o = 1, angGrid%n_angles
            do p=1, molRotGrid%n_angles
            icg=icg+1
            Delta_rho(i,j,k) = Delta_rho(i,j,k) + angGrid%weight(o) * cg_vect(icg)**2*molRotGrid%weight(p)
        END DO
        END DO
        END DO
    END DO
    END DO
    Delta_rho = Delta_rho-(twopi*fourpi)/real(sym_order,dp)
    !> Next FFT sequences can be done on multiple threads
    !> Compute rho in k-space
    fftw3%in_forward = Delta_rho
    call dfftw_execute ( fftw3%plan_forward )
    allocate ( rho_k ( nfft1 / 2 + 1 , nfft2 , nfft3 ) )
    rho_k = fftw3%out_forward
    ! Compute polarisation in k-space
    allocate ( Vpolarization_k ( nfft1 / 2 + 1 , nfft2 , nfft3 ) )
    ! V(k)=cs(k)*rho(k)
    do n = 1 , nfft3
    do m = 1 , nfft2
        do l = 1 , nfft1 / 2 + 1
        k_index = int ( norm_k ( l , m , n ) / delta_k ) + 1
        ! Here it happens that k_index gets higher than the highest c_k index.
        ! In this case one imposes k_index = k_index_max
        if ( k_index > nb_k ) k_index = nb_k
        ! V(k)=cs(k)*rho(k)
        Vpolarization_k ( l , m , n ) = rho_k ( l , m , n ) * c_s_hs ( k_index )
        END DO
    END DO
    END DO
    ! since rho(k) is now useless, deallocate associated array
    deallocate ( rho_k )
    ! FFT-1
    fftw3%in_backward = Vpolarization_k
    deallocate (Vpolarization_k)
    call dfftw_execute (fftw3%plan_backward)
    allocate ( Vpolarization ( nfft1 , nfft2 , nfft3 ) )
    Vpolarization = fftw3%out_backward / real(Nk,dp)
    ! compute excess energy and its gradient
    Fint = 0.0_dp ! excess energy
    icg = 0 ! index of cg_vect
    do species = 1 , nb_species
    fact = DeltaV * rho_0_multispec ( species ) !> facteur d'integration
    do i = 1 , nfft1
        do j = 1 , nfft2
        do k = 1 , nfft3
            Vint   = -kBT * rho_0_multispec ( species ) * Vpolarization(i,j,k)
            do o = 1 , angGrid%n_angles
            do p=1, molRotGrid%n_angles
            icg = icg + 1
            psi = CG_vect ( icg )
            Fint   = Fint   + angGrid%weight(o)*molRotGrid%weight(p) * fact * 0.5_dp * ( psi ** 2 - 1.0_dp) * Vint
    !         dF (icg) = dF ( icg ) + 2.0_dp * psi * angGrid%weight(o) * fact * Vint ! in case of bridge calculation, one deduces the pair contribution of hard spheres + => - and FF=FF-Fint
            dF (icg) = dF ( icg ) - 2.0_dp * psi * angGrid%weight(o) *molRotGrid%weight(p)* fact * Vint
            END DO
            END DO
        END DO
        END DO
    END DO
    END DO ! species
    deallocate(Vpolarization)
    ! conclude
    FF = FF - Fint
    call cpu_time(time1)
    IF (verbose) write(*,*) 'Fexc c_hs   = ' , Fint , 'computed in (sec)' , time1 - time0
 
END SUBROUTINE energy_cs_hard_sphere

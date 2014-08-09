!===================================================================================================================================
SUBROUTINE energy_nn_cs (Fint)
!===================================================================================================================================
! This subroutine computes the excess free energy in the HNC approximation of the homogeneous reference fluid approximation of MDFT.
! It uses the spherical projection c_s(k) that must be known at this point of the program.

    use precision_kinds ,only: i2b, dp
    use system          ,only: spaceGrid, nb_species, solvent, thermocond
    use dcf             ,only: c_s, delta_k
    use quadrature      ,only: molRotSymOrder, angGrid, molRotGrid
    use minimizer       ,only: cg_vect, FF, dF, from_cgvect_get_rho, from_rho_get_n, deallocate_solvent_rho
    use fft             ,only: fftw3, norm_k

    IMPLICIT NONE

    real(dp), intent(out)    :: Fint
    integer(i2b)             :: i, j, k, l, m, n, o, p, icg, s
    integer(i2b), pointer    :: nfft1, nfft2, nfft3
    real(dp)                 :: Vint, fact, psi, time1, time0
    logical, save            :: c_s_isOK =.FALSE.

    nfft1 => spacegrid%n_nodes(1)
    nfft2 => spacegrid%n_nodes(2)
    nfft3 => spacegrid%n_nodes(3)

    call cpu_time (time0)

    call from_cgvect_get_rho                        ! get solvent%rho from previous MDFT's minimization step
    call from_rho_get_n                             ! transform solvent%rho into solvent%n
    call deallocate_solvent_rho                     ! solvent%rho becomes useless
    fftw3%in_forward = solvent(1)%n - solvent(1)%n0 ! solvent%Dn
    call dfftw_execute (fftw3%plan_forward)         ! build solvent%Dn(k) in fftw3%outward
    call check_c_s
    do concurrent (i=1:nfft1/2+1, j=1:nfft2, k=1:nfft3)
        fftw3%in_backward(i,j,k) = fftw3%out_forward(i,j,k) * c_s(k_index(i,j,k))     ! gamma(k) = cs(k) * Dn(k)
    end do
    call dfftw_execute (fftw3%plan_backward)        ! Inverse Fourier transform gamma(k) => gamma(r)
    fftw3%out_backward = fftw3%out_backward / real(product(spacegrid%n_nodes),dp) ! Normalize the iFFT with FFTW3 conventions.

    ! excess free energy and its gradient
    Fint = -0.5_dp*thermocond%kbT * sum((solvent(1)%n-solvent(1)%n0) * fftw3%out_backward) * spacegrid%dv ! Gloria OOP !
    
    icg = 0
    DO s =1, nb_species
        fact = -thermocond%kbT * spaceGrid%dV * solvent(s)%rho0
        DO i =1,nfft1
            DO j =1,nfft2
                DO k =1,nfft3
                    Vint = fftw3%out_backward(i,j,k) * fact
                    DO o =1,angGrid%n_angles
                        DO p =1,molRotGrid%n_angles
                            icg = icg + 1
                            psi = CG_vect(icg)
                            dF(icg) = dF(icg) + 2.0_dp * psi * angGrid%weight(o) * molRotGrid%weight(p) * Vint
                        END DO
                    END DO
                END DO
            END DO
        END DO
    END DO

    FF = FF + Fint
    call cpu_time (time1)

 
    contains

        !===========================================================================================================================
        subroutine check_c_s
        !===========================================================================================================================   
            IF (.NOT. c_s_isOK) THEN
                IF (.NOT. ALLOCATED(c_s)) STOP "c_s(k) not allocated in energy_nn_cs.f90"
                IF (ALL(c_s==0._dp)) STOP "I find c_s(k) to be 0 everywhere in energy_nn_cs.f90"
                c_s_isOK = .TRUE.
            END IF
        end subroutine check_c_s
        !===========================================================================================================================
        
        
        !===========================================================================================================================
        pure function k_index (i,j,k)
        !===========================================================================================================================
            integer(i2b) :: k_index
            integer(i2b), intent(in) :: i, j, k
            k_index = MIN(  INT(norm_k(i,j,k)/delta_k)+1  ,  size(c_s)  ) ! Min prevents looking for undefined bin
        end function k_index
        !===========================================================================================================================
 
end subroutine energy_nn_cs

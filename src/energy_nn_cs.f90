! This subroutine computes the excess free energy in the HNC approximation of the homogeneous reference fluid approximation of MDFT.
! It uses the rotational invariant c_s(r), which is an input of the code.

SUBROUTINE energy_nn_cs (Fint)

    USE precision_kinds ,ONLY: i2b, dp
    USE system          ,ONLY: spaceGrid, nb_species, solvent, thermocond
    USE dcf             ,ONLY: c_s, nb_k, delta_k
    USE quadrature      ,ONLY: molRotSymOrder, angGrid, molRotGrid
    USE minimizer       ,ONLY: cg_vect, FF, dF, from_cgvect_get_rho, from_rho_get_n, deallocate_solvent_rho
    USE constants       ,ONLY: twopi, zeroC
    USE fft             ,ONLY: fftw3, norm_k

    IMPLICIT NONE

    REAL(dp), INTENT(OUT)    :: Fint
    INTEGER(i2b)             :: i, j, k, l, m, n, o, p, icg, s
    INTEGER(i2b), POINTER    :: nfft1, nfft2, nfft3
    REAL(dp)                 :: Vint, fact, psi, time1, time0
    LOGICAL    , SAVE        :: c_s_isOK =.FALSE.

    nfft1 => spacegrid%n_nodes(1)
    nfft2 => spacegrid%n_nodes(2)
    nfft3 => spacegrid%n_nodes(3)

    call cpu_time ( time0 )

    call from_cgvect_get_rho        ! get solvent%rho from previous MDFT's minimization step
    call from_rho_get_n             ! transform solvent%rho into solvent%n
    call deallocate_solvent_rho     ! solvent%rho becomes useless
    fftw3%in_forward = solvent(1)%n - solvent(1)%n0 ! solvent%Dn
    call dfftw_execute (fftw3%plan_forward)
    call check_c_s
    do concurrent (i=1:nfft1/2+1, j=1:nfft2, k=1:nfft3)
        fftw3%in_backward(i,j,k) = fftw3%out_forward(i,j,k) * c_s(k_index(i,j,k))     ! gamma(k)=cs(k)*Dn(k)
    end do
    call dfftw_execute (fftw3%plan_backward)
    fftw3%out_backward = fftw3%out_backward / real(product(spacegrid%n_nodes),dp)

    ! excess free energy and its gradient
    Fint = 0.0_dp
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
                            Fint = Fint + 0.5_dp * ( psi**2 - 1.0_dp) * angGrid%weight(o) * molRotGrid%weight(p) * Vint
                            dF(icg) = dF(icg) + 2.0_dp * psi          * angGrid%weight(o) * molRotGrid%weight(p) * Vint
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
            k_index = MIN(  INT(norm_k(i,j,k)/delta_k)+1  ,  nb_k  )
        end function k_index
        !===========================================================================================================================
 
end subroutine energy_nn_cs

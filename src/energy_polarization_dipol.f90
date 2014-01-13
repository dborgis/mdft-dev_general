SUBROUTINE energy_polarization_dipol (Fint)

    USE precision_kinds, ONLY : i2b, dp
    USE system,          ONLY : kBT, rho_0, spaceGrid
    USE quadrature,      ONLY : Omx, Omy, Omz, angGrid, molRotGrid
    USE minimizer,       ONLY : cg_vect , FF , dF
    USE constants,       ONLY : twopi
    USE fft,             ONLY : fftw3, k2, kproj, norm_k
    USE input,           ONLY : input_log, input_char, verbose
    USE dcf,             ONLY : c_delta , c_d, delta_k, nb_k
    
    IMPLICIT NONE
    
    INTEGER(i2b):: icg, i, j, k, l, m, n, m1, m2, m3, o, p
    INTEGER(i2b):: nfft1, nfft2, nfft3, nf1 , nf2 , nf3
    INTEGER(i2b):: k_index
    REAL(dp):: Lx, Ly, Lz, dx, dy, dz
    REAL(dp), INTENT(OUT) :: Fint ! Internal part of the free energy due to polarization
    REAL(dp):: Vint ! Dummy for calculation of Vint
    REAL(dp):: fact ! facteur d'integration
    REAL(dp):: rho, psi
    REAL(dp), ALLOCATABLE, DIMENSION (:,:,:) :: Px , Py , Pz , Ex , Ey , Ez
    REAL(dp), ALLOCATABLE, DIMENSION (:,:,:,:) :: polascal
    COMPLEX(dp), ALLOCATABLE, DIMENSION (:,:,:) :: Pkx , Pky , Pkz , Ekx , Eky , Ekz
    COMPLEX(dp) :: k_dot_P, pxt_k , pyt_k , pzt_k
    REAL(dp) :: c_deltat, c_dt ! dummy local values of c_delta and c_d in loops
    REAL(dp) :: time1, time0, time2, time3! timestamps
    REAL(dp) :: twopioLx, twopioLy, twopioLz ! dummy for 2pi/Lx, 2pi/Ly, 2pi/Lz
    REAL(dp) :: pxt, pyt, pzt, r, Ex_tmp, Ey_tmp, Ez_tmp
    REAL(dp) :: Nk
    REAL(dp), ALLOCATABLE, DIMENSION(:) :: weight_omx , weight_omy , weight_omz
    CHARACTER(50) :: filename
    
    nfft1 = spaceGrid%n_nodes(1); nfft2 = spaceGrid%n_nodes(2); nfft3 = spaceGrid%n_nodes(3)
    nf1 = nfft1/2; nf2 = nfft2/2; nf3 = nfft3/2
    Lx = spaceGrid%length(1); Ly = spaceGrid%length(2); Lz = spaceGrid%length(3)
    twopioLx = twopi/Lx; twopioLy = twopi/Ly; twopioLz = twopi/Lz
    dx = spaceGrid%dl(1); dy = spaceGrid%dl(2); dz = spaceGrid%dl(3)
    Nk = REAL ( nfft1*nfft2*nfft3 , dp ) ! total number of k grid points
    ! allocate and init polarization vector
    ALLOCATE ( Px (nfft1,nfft2,nfft3), SOURCE=0.0_dp )
    ALLOCATE ( Py (nfft1,nfft2,nfft3), SOURCE=0.0_dp )
    ALLOCATE ( Pz (nfft1,nfft2,nfft3), SOURCE=0.0_dp )

    IF (verbose) ALLOCATE ( Polascal ( nfft1 , nfft2 , nfft3 , 1), SOURCE=0.0_dp )
    
    ! put density of last minimization step in delta_rho and P
    ! but first prepare the product angGrid%weight(omega)*Omx in order not to repeat it
    ALLOCATE ( weight_omx ( angGrid%n_angles ), SOURCE=angGrid%weight*Omx )
    ALLOCATE ( weight_omy ( angGrid%n_angles ), SOURCE=angGrid%weight*Omy )
    ALLOCATE ( weight_omz ( angGrid%n_angles ), SOURCE=angGrid%weight*Omz )

    icg = 0
    DO i = 1 , nfft1
        DO j = 1 , nfft2
            DO k = 1 , nfft3    
                pxt = 0.0_dp
                pyt = 0.0_dp
                pzt = 0.0_dp
                DO o = 1 , angGrid%n_angles
                    DO p=1, molRotGrid%n_angles
                        icg = icg + 1
                        rho = cg_vect(icg) ** 2
                        pxt = pxt + weight_Omx(o) * molRotGrid%weight(p) * rho
                        pyt = pyt + weight_Omy(o) * molRotGrid%weight(p) * rho
                        pzt = pzt + weight_Omz(o) * molRotGrid%weight(p) * rho
                    END DO
                END DO
                Px (i,j,k) = pxt
                Py (i,j,k) = pyt
                Pz (i,j,k) = pzt
                IF (verbose) THEN
                    IF (i>nfft1/2) THEN; m1=i-1-nfft1; ELSE; m1=i-1; END IF
                    IF (j>nfft2/2) THEN; m2=j-1-nfft2; ELSE; m2=j-1; END IF
                    IF (k>nfft3/2) THEN; m3=k-1-nfft3; ELSE; m3=k-1; END IF
                    r = SQRT((m1*dx)**2+(m2*dy)**2+(m3*dz)**2)
                    polascal(i,j,k,1) = polascal(i,j,k,1) +pxt*m1*dx/r +pyt*m2*dy/r +pzt*m3*dz/r
                END IF
            END DO
        END DO
    END DO

    DEALLOCATE ( weight_omx )
    DEALLOCATE ( weight_omy )
    DEALLOCATE ( weight_omz )

    IF (verbose) THEN
        BLOCK
            REAL(dp), DIMENSION (nfft1,nfft2,nfft3,1) ::  polatot
            OPEN(11,FILE='output/polatotxmax')
                DO i=1,nfft1
                    WRITE(11,*) i*dx, 0.4894_dp*Px(i,nfft2/2+1,nfft3/2+1)
                END DO
            CLOSE(11)
            !Compute Radial Polarization
            filename='output/radial_polarization_dipolar'
            polatot(:,:,:,1)=sqrt(0.4894_dp**2*(Px(:,:,:)**2+Py(:,:,:)**2+Pz(:,:,:)**2))*rho_0
            CALL compute_rdf(polatot, filename)
            filename='output/radial_polarization_scalar'
            CALL compute_rdf(polascal*rho_0*0.4894_dp, filename)
            DEALLOCATE (polascal)
        END BLOCK
    END IF
  
    ! fourier transform px , py and pz
    ALLOCATE ( Pkx (nfft1/2+1, nfft2, nfft3) )
    ALLOCATE ( Pky (nfft1/2+1, nfft2, nfft3) )
    ALLOCATE ( Pkz (nfft1/2+1, nfft2, nfft3) )
    
    fftw3%in_forward = Px
    CALL dfftw_execute ( fftw3%plan_forward )
    Pkx = fftw3%out_forward
    
    fftw3%in_forward = Py
    CALL dfftw_execute ( fftw3%plan_forward )
    Pky = fftw3%out_forward
    
    fftw3%in_forward = Pz
    CALL dfftw_execute ( fftw3%plan_forward )
    Pkz = fftw3%out_forward 

    DEALLOCATE ( Px )
    DEALLOCATE ( Py )
    DEALLOCATE ( Pz )
    
    ! compute polarisation in k-space
    ALLOCATE ( Ekx (nfft1/2+1, nfft2, nfft3) )
    ALLOCATE ( Eky (nfft1/2+1, nfft2, nfft3) )
    ALLOCATE ( Ekz (nfft1/2+1, nfft2, nfft3) )
    
    ! get maximum number of k points as inputed in c_delta and c_d
    DO CONCURRENT (l=1:nfft1/2+1, m=1:nfft2, n=1:nfft3)
        pxt_k = Pkx(l,m,n)
        pyt_k = Pky(l,m,n)
        pzt_k = Pkz(l,m,n)
        IF ( k2(l,m,n) /= 0.0_dp) THEN
            k_dot_P = ( kproj(1,l)*pxt_k + kproj(2,m)*pyt_k + kproj(3,n)*pzt_k )/ k2(l,m,n)
        ELSE
            k_dot_P = CMPLX( tiny(1.0_dp),tiny(1.0_dp) ,dp )
        END IF
        k_index = MIN( INT( norm_k(l,m,n)/delta_k ) +1, nb_k)
        c_deltat = c_delta ( k_index )
        c_dt = c_d ( k_index )
        Ekx ( l,m,n ) = c_deltat * pxt_k + c_dt * ( 3.0_dp * k_dot_P * kproj(1,l) - pxt_k )
        Eky ( l,m,n ) = c_deltat * pyt_k + c_dt * ( 3.0_dp * k_dot_P * kproj(2,m) - pyt_k )
        Ekz ( l,m,n ) = c_deltat * pzt_k + c_dt * ( 3.0_dp * k_dot_P * kproj(3,n) - pzt_k )
    END DO
    
    DEALLOCATE ( Pkx )
    DEALLOCATE ( Pky )
    DEALLOCATE ( Pkz )
    
    ! inverse fourier transform the polarization field
    ! next inverse fourier transform sequence could be done in parallel
    ALLOCATE ( Ex ( nfft1 , nfft2 , nfft3 ) )
    ALLOCATE ( Ey ( nfft1 , nfft2 , nfft3 ) )
    ALLOCATE ( Ez ( nfft1 , nfft2 , nfft3 ) )
    
    fftw3%in_backward = Ekx
    CALL dfftw_execute ( fftw3%plan_backward )
    Ex = fftw3%out_backward / Nk
    
    fftw3%in_backward = Eky
    CALL dfftw_execute ( fftw3%plan_backward )
    Ey = fftw3%out_backward / Nk 
    
    fftw3%in_backward = Ekz
    CALL dfftw_execute ( fftw3%plan_backward )
    Ez = fftw3%out_backward / Nk 
    
    DEALLOCATE (Ekx)
    DEALLOCATE (Eky)
    DEALLOCATE (Ekz)
    
    fact = spaceGrid%dv * rho_0! integration factor
    
    Fint = 0.0_dp
    icg = 0
    DO i = 1, nfft1
    DO j = 1, nfft2
    DO k = 1, nfft3
        Ex_tmp = Ex(i,j,k)
        Ey_tmp = Ey(i,j,k)
        Ez_tmp = Ez(i,j,k)
        DO o = 1 , angGrid%n_angles
            Vint = -kBT*rho_0*( Omx(o)*Ex_tmp + Omy(o)*Ey_tmp + Omz(o)*Ez_tmp ) *angGrid%weight(o)
            DO p=1 , molRotGrid%n_angles
                icg = icg + 1
                psi = cg_vect ( icg )
                rho = psi ** 2
                Fint = Fint + (rho - 1.0_dp) * molRotGrid%weight(p) * Vint
                dF(icg) = dF(icg) + 2.0_dp * psi * molRotGrid%weight(p) * Vint * fact
            END DO  
        END DO
    END DO
    END DO
    END DO
    Fint = Fint * 0.5_dp * spaceGrid%dv * rho_0

    FF = FF + Fint

    DEALLOCATE ( Ex, Ey, Ez )

END SUBROUTINE energy_polarization_dipol

SUBROUTINE energy_polarization_dipol (Fint)

    USE precision_kinds, ONLY : i2b, dp
    USE system,          ONLY : thermocond, spaceGrid, solvent
    USE quadrature,      ONLY : Omx, Omy, Omz, angGrid, molRotGrid
    USE minimizer,       ONLY : cg_vect , FF , dF
    USE constants,       ONLY : twopi, zeroC
    USE fft,             ONLY : fftw3, k2, kproj, norm_k
    USE input,           ONLY : verbose
    use mathematica,     only : splint

    IMPLICIT NONE

    REAL(dp), INTENT(OUT) :: Fint ! Internal part of the free energy due to polarization
    REAL(dp), ALLOCATABLE, DIMENSION (:,:,:)    :: Px , Py , Pz , Ex , Ey , Ez
    COMPLEX(dp), ALLOCATABLE, DIMENSION (:,:,:) :: Pkx, Pky, Pkz, Ekx, Eky, Ekz
    REAL(dp) :: rho, psi, fact, Vint, Lx, Ly, Lz, Nk
    REAL(dp) :: pxt, pyt, pzt, r, Ex_tmp, Ey_tmp, Ez_tmp
    INTEGER(i2b):: icg, i, j, k, l, m, n, o, p, nfft1, nfft2, nfft3

    nfft1 = spaceGrid%n_nodes(1); nfft2 = spaceGrid%n_nodes(2); nfft3 = spaceGrid%n_nodes(3)
    Lx = spaceGrid%length(1); Ly = spaceGrid%length(2); Lz = spaceGrid%length(3)
    Nk = REAL ( nfft1*nfft2*nfft3 , dp ) ! total number of k grid points


    CALL build_polarization_vector_field_Px_Py_Pz ! Px(r) = \int_\Omega \Omega * \rho(r,\Omega) d\Omega  see Borgis et al., DOI:10.1103/PhysRevE.66.031206
    CALL build_Fourier_transformed_polarization_vector_field_Pkx_Pky_Pkz
    CALL build_excess_electric_field_in_Fourier_space
    CALL build_excess_electric_field_from_inverse_Fourier_transform


    fact = spaceGrid%dv * solvent(1)%rho0 ! integration factor
    Fint = 0.0_dp
    icg = 0
    DO i = 1, nfft1
    DO j = 1, nfft2
    DO k = 1, nfft3
        Ex_tmp = Ex(i,j,k)
        Ey_tmp = Ey(i,j,k)
        Ez_tmp = Ez(i,j,k)
        DO o = 1 , angGrid%n_angles
            Vint = -thermocond%kbT*solvent(1)%rho0 *( Omx(o)*Ex_tmp + Omy(o)*Ey_tmp + Omz(o)*Ez_tmp ) *angGrid%weight(o)
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
    Fint = Fint * 0.5_dp * spaceGrid%dv * solvent(1)%rho0
    FF = FF + Fint

    DEALLOCATE ( Ex, Ey, Ez )



    CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        SUBROUTINE build_polarization_vector_field_Px_Py_Pz
            REAL(dp), ALLOCATABLE, DIMENSION(:) :: weight_omx , weight_omy , weight_omz
            CHARACTER(50) :: filename
            ALLOCATE ( Px (nfft1,nfft2,nfft3), SOURCE=0.0_dp )
            ALLOCATE ( Py (nfft1,nfft2,nfft3), SOURCE=0.0_dp )
            ALLOCATE ( Pz (nfft1,nfft2,nfft3), SOURCE=0.0_dp )
            ! put density of last minimization step in delta_rho and P
            ! but first prepare the product angGrid%weight(omega)*Omx in order not to repeat it
            ALLOCATE ( weight_omx ( angGrid%n_angles ), SOURCE=angGrid%weight*Omx )
            ALLOCATE ( weight_omy ( angGrid%n_angles ), SOURCE=angGrid%weight*Omy )
            ALLOCATE ( weight_omz ( angGrid%n_angles ), SOURCE=angGrid%weight*Omz )
!            print*,anggrid%weight
!            print*,OMx
!            print*,OMy
!            print*,OMz
!            error stop "totoalaa"
            icg = 0
            DO i =1,nfft1
                DO j =1,nfft2
                    DO k = 1,nfft3
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
                    END DO
                END DO
            END DO
            DEALLOCATE ( weight_omx, weight_omy, weight_omz )
            IF (verbose) THEN
                BLOCK
                    REAL(dp), DIMENSION (nfft1,nfft2,nfft3,1) ::  polatot
        !~             OPEN(11,FILE='output/polatotxmax')
        !~                 DO i=1,nfft1
        !~                     WRITE(11,*) i*dx, 0.4894_dp*Px(i,nfft2/2+1,nfft3/2+1)
        !~                 END DO
        !~             CLOSE(11)
                    !Compute Radial Polarization
                    filename='output/radial_polarization_dipolar'
                    polatot(:,:,:,1)=sqrt(0.4894_dp**2*(Px(:,:,:)**2+Py(:,:,:)**2+Pz(:,:,:)**2))*solvent(1)%rho0
                    CALL output_rdf(polatot(:,:,:,1), filename)
!~                     filename='output/radial_polarization_scalar'
                END BLOCK
            END IF
        END SUBROUTINE build_polarization_vector_field_Px_Py_Pz

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        SUBROUTINE build_Fourier_transformed_polarization_vector_field_Pkx_Pky_Pkz
            ! fourier transform px , py and pz
            fftw3%in_forward = Px
            DEALLOCATE ( Px )
            CALL dfftw_execute ( fftw3%plan_forward )
            ALLOCATE ( Pkx (nfft1/2+1, nfft2, nfft3) ,SOURCE=fftw3%out_forward)

            fftw3%in_forward = Py
            DEALLOCATE ( Py )
            CALL dfftw_execute ( fftw3%plan_forward )
            ALLOCATE ( Pky (nfft1/2+1, nfft2, nfft3) ,SOURCE=fftw3%out_forward)

            fftw3%in_forward = Pz
            DEALLOCATE ( Pz )
            CALL dfftw_execute ( fftw3%plan_forward )
            ALLOCATE ( Pkz (nfft1/2+1, nfft2, nfft3) ,SOURCE=fftw3%out_forward)
        END SUBROUTINE build_Fourier_transformed_polarization_vector_field_Pkx_Pky_Pkz

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        SUBROUTINE build_excess_electric_field_in_Fourier_space
            USE dcf, ONLY : c_delta , c_d
            COMPLEX(dp) :: k_dot_P, pxt_k , pyt_k , pzt_k
            REAL(dp) :: c_delta_loc, c_d_loc, norm_k_loc
            ALLOCATE ( Ekx (nfft1/2+1, nfft2, nfft3) ,SOURCE=zeroC)
            ALLOCATE ( Eky (nfft1/2+1, nfft2, nfft3) ,SOURCE=zeroC)
            ALLOCATE ( Ekz (nfft1/2+1, nfft2, nfft3) ,SOURCE=zeroC)
            ! get maximum number of k points as inputed in c_delta and c_d
            DO CONCURRENT (l=1:nfft1/2+1, m=1:nfft2, n=1:nfft3)
                pxt_k = Pkx(l,m,n)
                pyt_k = Pky(l,m,n)
                pzt_k = Pkz(l,m,n)
                IF ( k2(l,m,n) > epsilon(1._dp) ) then
                    k_dot_P = ( kproj(1,l)*pxt_k + kproj(2,m)*pyt_k + kproj(3,n)*pzt_k )/ k2(l,m,n)
                ELSE
                    k_dot_P = CMPLX( 0._dp,0._dp ,dp )
                END IF
                norm_k_loc = norm_k(l,m,n)
                call splint( xa=c_delta%x, ya=c_delta%y, y2a=c_delta%y2, n=size(c_delta%y), x=norm_k_loc, y=c_delta_loc)
                call splint( xa=c_d%x, ya=c_d%y, y2a=c_d%y2, n=size(c_d%y), x=norm_k_loc, y=c_d_loc)
                Ekx ( l,m,n ) = c_delta_loc * pxt_k + c_d_loc * ( 3.0_dp * k_dot_P * kproj(1,l) - pxt_k )
                Eky ( l,m,n ) = c_delta_loc * pyt_k + c_d_loc * ( 3.0_dp * k_dot_P * kproj(2,m) - pyt_k )
                Ekz ( l,m,n ) = c_delta_loc * pzt_k + c_d_loc * ( 3.0_dp * k_dot_P * kproj(3,n) - pzt_k )
            END DO
            DEALLOCATE ( Pkx, Pky, Pkz )
        END SUBROUTINE build_excess_electric_field_in_Fourier_space

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        SUBROUTINE build_excess_electric_field_from_inverse_Fourier_transform
            ! inverse fourier transform the electric field
            fftw3%in_backward = Ekx
            DEALLOCATE (Ekx)
            CALL dfftw_execute ( fftw3%plan_backward )
            ALLOCATE ( Ex ( nfft1 , nfft2 , nfft3 ) ,SOURCE=fftw3%out_backward / Nk)
            fftw3%in_backward = Eky
            DEALLOCATE (Eky)
            CALL dfftw_execute ( fftw3%plan_backward )
            ALLOCATE ( Ey ( nfft1 , nfft2 , nfft3 ) ,SOURCE=fftw3%out_backward / Nk)
            fftw3%in_backward = Ekz
            DEALLOCATE (Ekz)
            CALL dfftw_execute ( fftw3%plan_backward )
            ALLOCATE ( Ez ( nfft1 , nfft2 , nfft3 ) ,SOURCE=fftw3%out_backward / Nk)
        END SUBROUTINE build_excess_electric_field_from_inverse_Fourier_transform

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE energy_polarization_dipol

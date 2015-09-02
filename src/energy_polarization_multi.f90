SUBROUTINE energy_polarization_multi (F_pol)
!
    use precision_kinds, ONLY: i2b, dp
    use module_grid, only: grid
!     use system, ONLY: thermocond, solvent(1)%nspec, grid, solvent
!     use dcf, ONLY: chi_l, chi_t, nb_k, delta_k
!     use quadrature,ONLY : angGrid, molRotGrid
!     use module_minimizer, ONLY: cg_vect_new, FF, dF_new
!     use constants, ONLY : twopi, fourpi, qfact, zeroC
!     use fft, ONLY: fftw3, kx, ky, kz, k2, norm_k
!     use module_input, ONLY: verbose
!
!     IMPLICIT NONE
!
    REAL(dp), INTENT(OUT) :: F_pol
!     REAL(dp) :: F_pol_long, F_pol_trans , F_pol_tot  ! Longitudinal, transverse and total Polarization free energy
!     REAL(dp) :: facsym, deltaVkn, Lweight, kvec(3)
!     REAL(dp),    ALLOCATABLE, DIMENSION (:,:,:,:,:,:) :: rho, dF_pol_tot
!     COMPLEX(dp), ALLOCATABLE, DIMENSION (:,:,:,:,:,:) :: dF_pol_tot_k
!     COMPLEX(dp), ALLOCATABLE, DIMENSION (:,:,:,:,:) :: Ptrans_k, Plong_k, Ptot_k
!     complex(dp), allocatable :: rho_k(:,:,:)
!     INTEGER(i2b) :: icg, i, j, k, o, p, n, s, k_index, d
!     COMPLEX(dp) :: k_tens_k_P(3), toto, temp, Ptrans_k_loc(3), Plong_k_loc(3), Ptot_k_loc(3), molec_polar_k_loc(3),&
!                    dF_pol_long_k, dF_pol_trans_k
!     integer(i2b), pointer :: nfft1 => grid%n_nodes(1), nfft2 => grid%n_nodes(2), nfft3 => grid%n_nodes(3)
!
!     if (size(solvent)/=1) STOP 'energy_polarization_multi.f90 not working for multisolvent species'
!
!     IF ( (SIZE(chi_l)/=SIZE(chi_t))) THEN
!         WRITE(*,*)"chi_l and chi_t should have the same number of points, at least for now"
!         PRINT*, SIZE(chi_l), SIZE(chi_t)
!         STOP
!     END IF
!
! !===================================================================================================================================
! !                	Initialization
! !
! !===================================================================================================================================
!     deltaVkn = 1._dp/PRODUCT(grid%length)
!     F_pol = 0.0_dp
!     F_pol_long = 0.0_dp
!     F_pol_trans = 0.0_dp
!     F_pol_tot = 0.0_dp
! !======================================================================================================
! !            ====================================================
! !                 	Compute density in real and 		        !
! !       		            Fourier space		            	!
! !							                                    !
! !            ====================================================
!     ALLOCATE ( rho (grid%n_nodes(1), grid%n_nodes(2), grid%n_nodes(3),&
!                         angGrid%n_angles, molRotGrid%n_angles, solvent(1)%nspec), SOURCE=0._dp )
!     icg=0
!     DO s=1, solvent(1)%nspec
!         DO i=1, nfft1
!             DO j=1, nfft2
!                 DO k=1, nfft3
!                     DO o=1, angGrid%n_angles
!                         DO p=1, molRotGrid%n_angles
!                             icg = icg + 1
!                             rho(i,j,k,o,p,s) = cg_vect_new(i,j,k,o,p,s)**2 ! this is rho(r,o)/rho0
!                         END DO
!                     END DO
!                 END DO
!             END DO
!         END DO
!     END DO
!
! !===================================================================================================================================
! !            !    		get Density 			!
! !            !			in Fourier Space		and total polarization!
! !===================================================================================================================================
!     allocate (rho_k (nfft1/2+1, nfft2, nfft3), SOURCE=zeroC)
!     allocate (Ptot_k (3,nfft1/2+1, nfft2, nfft3, solvent(1)%nspec), SOURCE=zeroC )
!     DO s = 1, solvent(1)%nspec
!         DO p = 1, molRotGrid%n_angles
!             DO o = 1, angGrid%n_angles
!                 fftw3%in_forward = rho(:,:,:,o,p,s)
!                 CALL dfftw_execute (fftw3%plan_forward)
!                 rho_k = fftw3%out_forward *grid%dv
!                 do concurrent (d=1:3)
!                   Ptot_k(d,:,:,:,s) =  Ptot_k(d,:,:,:,s)+&
!                     rho_k * angGrid%weight(o) * molRotGrid%weight(p) * solvent(s)%molec_polar_k(d,:,:,:,o,p)
!                 end do
!             END DO
!         END DO
!     END DO
!     DEALLOCATE (rho, rho_k)
!
! !... CONSTRUCT THE ELECTRONIC POLARIZATION VECTOR FIELD IN K-SPACE, M(k)
! !...
! !...
! ! Ptot_k = Ptot_k + M_k
!
!
! !            ====================================================
! !            !    	Compute 				!
! !            !	Transverse and longitudinal Polarization	!
! !            ====================================================
!     allocate (Plong_k (3,nfft1/2+1, nfft2, nfft3, solvent(1)%nspec), SOURCE=zeroC )
!     do concurrent ( i=1:nfft1/2+1, j=1:nfft2, k=1:nfft3, s=1:solvent(1)%nspec )
!         if (abs(k2(i,j,k))<=epsilon(1._dp)) THEN
!             Plong_k(:,i,j,k,s) = zeroC
!         else
!             kvec = [kx(i),ky(j),kz(k)]
!             Plong_k(:,i,j,k,s) = kvec * sum(kvec * Ptot_k(:,i,j,k,s)) / k2(i,j,k)
!         end if
!     end do
!
!     allocate (Ptrans_k (3,nfft1/2+1, nfft2, nfft3, solvent(1)%nspec), SOURCE=(Ptot_k - Plong_k) )
!
! !            ====================================================
! !            !    	Compute Free energy due to		!
! !            !	Transverse and longitudinal Polarization	!
! !            ====================================================
!     ALLOCATE( dF_pol_tot_k (nfft1/2+1, nfft2, nfft3, angGrid%n_angles, molRotGrid%n_angles, solvent(1)%nspec), source=zeroC)
!
!     DO CONCURRENT ( i=1:nfft1/2+1, j=1:nfft2, k=1:nfft3, s=1:solvent(1)%nspec )
!
!         IF (i>1 .AND. i<nfft1/2+1) THEN
!             facsym=2.0_dp
!         ELSE
!             facsym=1.0_dp
!         END IF
!
!         Ptrans_k_loc = Ptrans_k(:,i,j,k,s)
!         Plong_k_loc  = Plong_k(:,i,j,k,s)
!         Ptot_k_loc   = Ptot_k(:,i,j,k,s)
!
!         k_index = INT ( norm_k(i,j,k) / delta_k ) + 1
!         IF ( k_index > nb_k ) k_index = nb_k ! Here it happens that k_index gets higher than the highest c_k index. In this case one imposes k_index = k_index_max
!
!         toto = deltaVkn*0.5_dp*qfact*solvent(s)%rho0**2*facsym!/(twopi)**3
!
!         F_pol_long = F_pol_long + fourpi * toto * dot_product( Plong_k_loc, Plong_k_loc)    /chi_l(k_index)
!
!         F_pol_trans = F_pol_trans +        toto * dot_product( Ptrans_k_loc, Ptrans_k_loc)  /chi_t(k_index)
!
!         F_pol_tot = F_pol_tot &
!                    - thermocond%kbT*3/(2*norm2(solvent(s)%dipole)**2*solvent(s)%n0) *deltaVkn*facsym*solvent(s)%rho0**2&
!                         * dot_product( Ptot_k_loc , Ptot_k_loc )
!
!         DO CONCURRENT ( o=1:angGrid%n_angles, p=1:molRotGrid%n_angles )
!             Lweight = angGrid%weight(o) * molRotGrid%weight(p)
!             molec_polar_k_loc = solvent(s)%molec_polar_k(:,i,j,k,o,p)
! !========================================================================================================================
! !Evaluate gradient
! !========================================================================================================================
!             IF ( k2(i,j,k)<=epsilon(1.0_dp) ) THEN
!                 dF_pol_long_k = zeroC
!                 dF_pol_trans_k = solvent(s)%rho0*0.5_dp*qfact/chi_t(k_index)*&
!                     dot_product( molec_polar_k_loc , Ptrans_k_loc) *Lweight*2.0_dp
!             ELSE
!                 kvec = [kx(i),ky(j),kz(k)]
!                 k_tens_k_P = kvec * dot_product(    molec_polar_k_loc   ,kvec)
!
!                 dF_pol_long_k = solvent(s)%rho0*0.5_dp* qfact/chi_l(k_index)*fourpi*Lweight*&
!                     sum( k_tens_k_P * Plong_k_loc ) /k2(i,j,k)*2.0_dp
!
!                 dF_pol_trans_k = 0.5_dp*solvent(s)%rho0 *qfact /chi_t(k_index) *Lweight *2.0_dp *&
!                     sum(   Ptrans_k_loc * (conjg(molec_polar_k_loc) - k_tens_k_P/k2(i,j,k))      )
!             END IF
!
!             dF_pol_tot_k(i,j,k,o,p,s) = dF_pol_long_k + dF_pol_trans_k &
!                 -thermocond%kbT*3._dp*solvent(s)%rho0/(2._dp*norm2(solvent(s)%dipole)**2&
!                 *solvent(s)%n0)*Lweight*2.0_dp*&
!                 dot_product(   molec_polar_k_loc     ,     Ptot_k_loc   )
!
!         END DO
!     END DO
!
! !========================================================================================================================
! !						Get gradient in real space
! !========================================================================================================================
!     ALLOCATE( dF_pol_tot (nfft1, nfft2, nfft3, angGrid%n_angles, molRotGrid%n_angles, solvent(1)%nspec), source=0._dp)
!     DO s=1,solvent(1)%nspec
!         DO p=1, molRotGrid%n_angles
!             DO o=1,angGrid%n_angles
!                 fftw3%in_backward= dF_pol_tot_k (:,:,:,o,p,s)
!                 call dfftw_execute (fftw3%plan_backward)
!                 dF_pol_tot (:,:,:,o,p,s) = fftw3%out_backward *deltaVkn
!             END DO
!         END DO
!     END DO
! !========================================================================================================================
! !						Allocate it for minimizing
! !========================================================================================================================
!     icg = 0
!     DO s=1, solvent(1)%nspec
!         DO i=1, nfft1
!             DO j=1, nfft2
!                 DO k=1,nfft3
!                     DO o=1,angGrid%n_angles
!                         DO p=1, molRotGrid%n_angles
!                             icg = icg + 1
!                             dF_new(i,j,k,o,p,s) = dF_new(i,j,k,o,p,s) + dF_pol_tot (i,j,k,o,p,s) * cg_vect_new(i,j,k,o,p,s) * solvent(s)%rho0 * 2.0_dp * grid%dv
!                         END DO
!                     END DO
!                 END DO
!             END DO
!         END DO
!     END DO
!
!     F_pol = F_pol_tot + F_pol_long + F_pol_trans
!     FF = FF + F_pol
!
!     DEALLOCATE (dF_pol_tot, dF_pol_tot_k)
stop "energy_polarization_multi desactivated since o,p => io"
END SUBROUTINE energy_polarization_multi

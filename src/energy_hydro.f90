!This routine evaluate the excess free energy plus an hydrophobic part
SUBROUTINE energy_hydro (Fint)
    !
    use precision_kinds     ,ONLY: dp, i2b
    use module_grid, only: grid
    ! use system              ,ONLY: thermocond, solvent(1)%nspec, grid, solvent, solvent(1)%nspec
    ! use dcf                 ,ONLY: c_s, nb_k, delta_k, c_s_hs
    ! use constants           ,ONLY: fourpi,i_complex,twopi,zerodp,zeroC
    ! use minimizer           ,ONLY: cg_vect_new, FF, dF_new
    ! use quadrature          ,ONLY: molRotSymOrder, angGrid, molRotGrid
    ! use fft                 ,ONLY: fftw3,norm_k,kx,ky,kz,k2,timesExpPrefactork2
    ! use module_input               ,ONLY: getinput%log, verbose, getinput%dp
    ! use mathematica         ,only: splint
    !
    IMPLICIT NONE
    !
    real(dp), intent(out) :: Fint
    ! REAL(dp) :: mu_0 ! phenomenological potential
    ! REAL(dp) :: Nk ! total number of k points. Real because really often used as a divisor of a real number
    ! REAL(dp) :: deltaVk ! volume unit per k point
    ! REAL(dp) :: fact_n ! integration factor = deltaV*n_0 = deltaV*rho_0/fourpi
    ! REAL(dp) :: prefactor ! dummy for a prefactor calculated a large number of times  =-R_cg**2/2.0_dp
    ! COMPLEX(dp), ALLOCATABLE, DIMENSION(:,:,:) :: dS_cgk, dF_cgk ! coarse grained part of dF in k space
    ! COMPLEX(dp), ALLOCATABLE, DIMENSION(:,:,:) :: delta_nk, delta_nbark, V_nk, V_nbark ! all names in k are Fourier transformed versions of non-k
    ! REAL(dp)   , ALLOCATABLE, DIMENSION(:,:,:) :: delta_n, delta_nbar, V_n, sum_grad_delta_nbar_sq ! all names in " "bar are coarse grained version of non-bar
    ! REAL(dp)   , ALLOCATABLE, DIMENSION(:,:,:) :: dF_cg ! coarse grained part of dF
    ! REAL(dp) :: Vint ! excess energy one wants to compute here
    ! REAL(dp) :: S_cg, F_cg ! surface part of the energy (Cahn Hilliard)
    ! REAL(dp) :: time0 , time1 ! timesteps
    ! REAL(dp) :: R_cg ! radius of the coarse graining gaussian function
    ! REAL(dp) :: a, mm, b1, b2, factn, c_s_loc, c_s_hs_loc
    ! INTEGER(i2b) :: icg,i,j,k,o,l,m,p,nfft1,nfft2,nfft3,s
    !
    ! CALL CPU_TIME (time0)
    ! !print*, c_s(1) , c_s_hs(1)
    ! !a=c_s(1)-c_s_hs(1)
    ! !print*, 'a = ' , a
    ! !factn=4/3*r_HS^3*pi
    ! ! for r_HS=1.27 angstrom
    ! factn=8.575896827
    ! a=5.228*2.0_dp*factn
    ! b1=5.80697224_dp*factn**3
    ! b2=9.609322_dp*factn**4
    ! mm=100.0_dp!/(n_0**2)/30.0_dp
    ! IF (verbose) THEN
    !     print*, a,12.3*2.0_dp*8.1812242666880284_dp
    !     print*, 'Hydrophobicity Van der Waals way'
    ! END IF
    !
    !
    ! nfft1 = grid%n_nodes(1)
    ! nfft2 = grid%n_nodes(2)
    ! nfft3 = grid%n_nodes(3)
    !
    ! !In this routine we use the pair correlation function of a hard sphere fluid.
    ! IF (.NOT. ALLOCATED(c_s_hs%x)) THEN
    !     PRINT*,"I'm going through cs_of_k_hard_sphere"
    !     CALL cs_of_k_hard_sphere
    !     IF (ALL(c_s_hs%y==zerodp)) STOP "problem in the cs hard sphere generated from call from energy_hydro.f90"
    ! END IF
    ! IF (.not. getinput%log('hard_sphere_fluid')) THEN
    !     PRINT*, 'energy_hydro MUST be used with hard sphere solvent correction  TAG of hard_sphere_solute in dft.in must be T'
    !     STOP
    ! END IF
    !
    ! ! following count does not take into account all cases where angGrid%n_angles /=1
    ! IF ( angGrid%n_angles /= 1 ) THEN
    !     PRINT*, 'angGrid%n_angles /= 1 and hydrophobic part is calculated while not still implemented'
    !     PRINT*, 'not possible for now. error in cs_plus_hydro.f90'
    !     STOP
    ! END IF
    !
    ! ! this SUBROUTINE only works for solvent(1)%nspec = 1
    ! IF ( solvent(1)%nspec /= 1 ) THEN
    !     PRINT*, 'solvent(1)%nspec /= 1 and this is not implemented yet in cs_plus_hydro.f90'
    !     STOP
    ! END IF
    !
    ! ! macroscopic water parameters, from Chandler
    ! mu_0 = 7.16d-4*thermocond%kbT ! phenomenological potential
    ! R_cg = getinput%dp("hydro_coarsegrainingradius") ! Gaussian radius for coarse graining
    !
    ! ! total number of kpoints and volume per kpoint
    ! Nk = REAL(nfft1*nfft2*nfft3,dp)
    ! DeltaVk = grid%dV/Nk !twopi**3/(Lx*ly*Lz*Nk)!
    ! fact_n = grid%dV*solvent(1)%n0 ! integration factor
    ! prefactor = -R_cg**2 / 2.0_dp ! dummy for later in the gaussian in k space g(k)=exp(prefactor*k2)
    !
    ! BLOCK
    !     REAL(dp) :: delta_n_ijk ! dummy for local delta_n
    !     ! get density from last minimization step
    !     ALLOCATE(delta_n(nfft1,nfft2,nfft3) ,SOURCE=zerodp)
    !     ! HERE WE WORK WITH delta_n = delta_n normalized wrt n0 = delta_n/n0 = (n-n0)/n0
    !     icg=0
    !     s=1
    !     DO i=1,nfft1
    !         DO j=1,nfft2
    !             DO k=1,nfft3
    !                 delta_n_ijk =0.0_dp ! init delta_n
    !                 DO o=1,angGrid%n_angles
    !                     DO p=1, molRotGrid%n_angles
    !                         icg = icg +1
    !                         delta_n_ijk = delta_n_ijk + angGrid%weight(o)*molRotGrid%weight(p) * cg_vect_new(i,j,k,o,p,s) ** 2 ! sum over all orientations
    !                     END DO
    !                 END DO
    !                 delta_n(i,j,k) = delta_n_ijk
    !             END DO
    !         END DO
    !     END DO
    !     delta_n = delta_n *molRotSymOrder/(twopi*fourpi) -1.0_dp ! normalize (n=1/fourpi int_o rho(r,o))
    ! END BLOCK
    !
    ! ! TODO Next FFT sequences can be done on multiple threads
    ! ! compute delta_nk which is FFT(delta_n)
    ! fftw3%in_forward = delta_n
    ! CALL dfftw_execute ( fftw3%plan_forward )
    ! ALLOCATE ( delta_nk (nfft1/2+1,nfft2,nfft3) ,SOURCE=fftw3%out_forward)
    !
    ! ! the coarse-grained delta_n is simply delta_n multiplied by a gaussian distribution of empirical
    ! ALLOCATE ( delta_nbark ( nfft1/2+1,nfft2,nfft3) ,SOURCE=zeroC)
    ! delta_nbark = timesExpPrefactork2 (delta_nk,prefactor)
    !
    ! ! surface energy (cahn hilliard part of the intrinsic excess energy)
    ! ALLOCATE ( dS_cgk (nfft1/2+1,nfft2,nfft3) ,SOURCE=delta_nbark)
    ! BLOCK
    !     REAL(dp) :: facsym ! 1 or 2 for taking into account time reversal symetry of FFT(real)
    !     DO CONCURRENT (l=1:nfft1/2+1, m=1:nfft2, p=1:nfft3)
    !         IF ( l >= 2 .and. l <= nfft1/2 ) THEN ! take time reversal symetry into account small part
    !             facsym = 2.0_dp
    !         ELSE
    !             facsym = 1.0_dp
    !         END IF
    !         dS_cgk(l,m,p) = dS_cgk(l,m,p) * facsym * k2 (l,m,p)
    !     END DO
    ! dS_cgk = dS_cgk * solvent(1)%n0**2 * thermocond%kbT * mm * DeltaVk
    ! END BLOCK
    !
    ! ! gradient of delta_nbar in kspace
    ! ! remember that the fourier transform of the gradient of f(r) is i*vector_k*f(k)
    ! DO CONCURRENT (l=1:nfft1/2+1)
    !     fftw3%in_backward(l,:,:) = delta_nbark(l,:,:) * kx(l) * i_complex
    ! END DO
    ! CALL dfftw_execute(fftw3%plan_backward)
    ! ALLOCATE ( sum_grad_delta_nbar_sq (nfft1,nfft2,nfft3) ,SOURCE=(fftw3%out_backward/Nk)**2)
    ! DO CONCURRENT (l=1:nfft2)
    !     fftw3%in_backward(:,l,:) = delta_nbark(:,l,:) * ky(l) * i_complex
    ! END DO
    ! CALL dfftw_execute(fftw3%plan_backward)
    ! sum_grad_delta_nbar_sq = sum_grad_delta_nbar_sq + (fftw3%out_backward/Nk)**2
    ! DO l=1,nfft3
    !     fftw3%in_backward(:,:,l) = delta_nbark(:,:,l) * kz(l) * i_complex
    ! END DO
    ! CALL dfftw_execute(fftw3%plan_backward)
    ! sum_grad_delta_nbar_sq = sum_grad_delta_nbar_sq + (fftw3%out_backward/Nk)**2
    !
    ! ! large gradient part
    ! S_cg = 0.5_dp*mm*grid%dv*SUM(sum_grad_delta_nbar_sq)*solvent(1)%n0**2*thermocond%kBT
    ! DEALLOCATE (sum_grad_delta_nbar_sq)
    !
    ! fftw3%in_backward = delta_nbark
    ! CALL dfftw_execute(fftw3%plan_backward)
    ! ALLOCATE ( delta_nbar(nfft1,nfft2,nfft3) ,SOURCE=fftw3%out_backward/Nk)
    !
    ! F_cg = -0.5_dp*a*grid%dv*SUM(delta_nbar**2)*solvent(1)%n0**2*thermocond%kBT&
    !              +b1*grid%dv*SUM(delta_nbar**3)*solvent(1)%n0**3*thermocond%kBT&
    !              +b2*grid%dv*SUM(delta_nbar**4)*solvent(1)%n0**4*thermocond%kBT
    ! IF (verbose) THEN
    !     print *, 'F_CG =  ',F_cg
    !     print *, 'S_CG =  ',S_cg
    ! END IF
    !
    ! ! intrinsic local excess energy (second order term in Fexc) in k spae
    ! ALLOCATE ( V_nk (nfft1/2+1,nfft2,nfft3) ,SOURCE=(delta_nk-delta_nbark))
    ! DEALLOCATE (delta_nk, delta_nbark)
    ! DO CONCURRENT (l=1:nfft1/2+1, m=1:nfft2, p=1:nfft3)
    !   call splint( xa=c_s_hs%x, ya=c_s_hs%y, y2a=c_s_hs%y2, n=size(c_s_hs%y), x=norm_k(l,m,p), y=c_s_hs_loc)
    !   call splint( xa=c_s%x, ya=c_s%y, y2a=c_s%y2, n=size(c_s%y), x=norm_k(l,m,p), y=c_s_loc)
    !     V_nk (l,m,p) = V_nk(l,m,p)*(c_s_loc-c_s_hs_loc) ! we suppose c_s(nbar)=cst=c_s(n_0) which is a crude approximation ! compute radial excess potential in kspace
    ! END DO
    !
    ! ! coarse grain V_nk
    ! ALLOCATE ( V_nbark (nfft1/2+1,nfft2,nfft3) ,SOURCE=zeroC)
    ! V_nbark = timesExpPrefactork2 (V_nk,prefactor)
    !
    ! ! compute V_n(r) by FFT-1
    ! fftw3%in_backward = V_nk
    ! DEALLOCATE (V_nk)
    ! CALL dfftw_execute(fftw3%plan_backward)
    ! ALLOCATE ( V_n (nfft1,nfft2,nfft3) ,SOURCE=fftw3%out_backward/Nk)
    !
    ! ! compute coarse grained objects in real space by FFT-1
    ! fftw3%in_backward = V_nbark
    ! DEALLOCATE(V_nbark)
    ! CALL dfftw_execute(fftw3%plan_backward)
    ! V_n = V_n - (fftw3%out_backward/Nk) !V_n = V_n - V_nbar
    !
    ! ! finaly, the total energy
    ! Fint = 0.5_dp*SUM( (delta_n-delta_nbar)*V_n )*(-thermocond%kbT*solvent(1)%n0*fact_n)
    ! DEALLOCATE(delta_n)
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!
    ! !       gradients
    ! !!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! ! compute coarse grained free energy and gradient small part
    ! ALLOCATE ( dF_cg (nfft1,nfft2,nfft3) ,SOURCE=zerodp)
    ! dF_cg = (  delta_nbar    *(-a*solvent(1)%n0**2) &
    !          + delta_nbar**2 *(3._dp*b1*solvent(1)%n0**3) &
    !          + delta_nbar**3 *(4._dp*b2*solvent(1)%n0**4) ) *thermocond%kbT*grid%dv
    ! DEALLOCATE (delta_nbar)
    !
    ! ! FFT (dF_cg) in order to work in kspace
    ! fftw3%in_forward = dF_cg
    ! CALL dfftw_execute(fftw3%plan_forward)
    ! ALLOCATE ( dF_cgk (nfft1/2+1,nfft2,nfft3) ,SOURCE=zeroC)
    !
    ! ! define coarse grained gradients in k space
    ! dF_cgk = timesExpPrefactork2 ( fftw3%out_forward + dS_cgk, prefactor)
    ! DEALLOCATE(dS_cgk)
    !
    ! ! get back coarse grained gradients in r space
    ! fftw3%in_backward = dF_cgk
    ! DEALLOCATE ( dF_cgk )
    ! CALL dfftw_execute ( fftw3%plan_backward )
    ! dF_cg = fftw3%out_backward/Nk
    !
    !
    ! ! compute final FF and dF_new
    ! BLOCK
    !     REAL(dp) :: aa, dF_cg_ijk, psi
    !     icg = 0
    !     s=1
    !     DO i = 1, nfft1
    !         DO j = 1, nfft2
    !             DO k = 1, nfft3
    !                 Vint = -thermocond%kBT * solvent(1)%n0 * fact_n * V_n(i,j,k)
    !                 dF_cg_ijk = dF_cg(i,j,k)
    !                 aa = 1._dp/(fourpi*twopi/molRotSymOrder)*(Vint+dF_cg_ijk)
    !                 DO o = 1, angGrid%n_angles
    !                     DO p=1,molRotGrid%n_angles
    !                         icg = icg + 1
    !                         psi = cg_vect_new(i,j,k,o,p,s)
    !                         dF_new(i,j,k,o,p,s) = dF_new(i,j,k,o,p,s) +2.0_dp*psi*molRotGrid%weight(p)*angGrid%weight(o)*aa
    !                     END DO
    !                 END DO
    !             END DO
    !         END DO
    !     END DO
    ! END BLOCK
    !
    !
    ! Fint = Fint + F_cg + S_cg
    ! FF = FF + Fint
    !
    !
    ! CALL CPU_TIME ( time1 )
stop "energy hydro is desactivated since o,p => io"
END SUBROUTINE

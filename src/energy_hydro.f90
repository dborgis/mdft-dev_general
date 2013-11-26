!This routine evaluate the excess free energy plus an hydrophobic part
SUBROUTINE energy_hydro (Fint)
 USE precision_kinds,only : dp , i2b
  use system,only : nfft1 , nfft2 , nfft3 , deltaV, c_s_hs ,  kBT , nb_species,n_0,   Lx,Ly,Lz
    USE dcf, ONLY: c_s, nb_k, delta_k
  use constants,only : fourpi , i_complex,twopi
  USE minimizer, ONLY: cg_vect , FF , dF
  use quadrature, only: sym_order, angGrid, molRotGrid
  use fft,only : fftw3 , norm_k,kx,ky,kz,k2,&
                timesExpPrefactork2
  use input, only : input_log, verbose
  
  IMPLICIT NONE
  real(dp) :: mu_0 ! phenomenological potential
  real(dp) :: Nk ! total number of k points. Real because really often used as a divisor of a real number
  real(dp) :: deltaVk ! volume unit per k point
  real(dp) :: fact_n ! integration factor = deltaV*n_0 = deltaV*rho_0/fourpi
  integer(i2b) :: icg,i,j,k,o,l,m,p ! dummy for loops
  real(dp), allocatable, dimension(:,:,:) :: delta_n ! =n-1
  real(dp), allocatable, dimension(:,:,:) :: delta_nbar ! coarse grained delta_n
  complex(dp):: delta_nbark_loc ! dummy local value for delta_nbark(i,j,k)
  complex(dp), allocatable, dimension(:,:,:) :: delta_nk ! delta_n in kspace
  complex(dp), allocatable, dimension(:,:,:) :: delta_nbark ! delta_nbar in kspace
  real(dp) :: facsym ! 1 or 2 for taking into account time reversal symetry of FFT(real)
  real(dp) :: k2_loc ! dummy local k2
  integer(i2b) :: k_index ! index of the k point for reading c_s(k_index)
  !real(dp) :: norm_k ! norm of k vector == sqrt(kx**2+ky**2+kz**2)
  real(dp) :: prefactor ! dummy for a prefactor calculated a large number of times  =-R_cg**2/2.0_dp
  real(dp), allocatable, dimension(:,:,:) :: V_n ! radial excess potential
  complex(dp), allocatable, dimension(:,:,:) :: V_nk ! FFT(V_n)
  real(dp), allocatable, dimension(:,:,:) :: V_nbar ! coarse grained V_n
  complex(dp), allocatable, dimension(:,:,:) :: V_nbark ! FFT(V_n)
  real(dp), allocatable, dimension(:,:,:) :: gradx_delta_nbar , grady_delta_nbar , gradz_delta_nbar ! grad delta_nbar
  complex(dp), allocatable, dimension(:,:,:) :: gradx_delta_nbark , grady_delta_nbark , gradz_delta_nbark ! gradient of delta_nbark
  real(dp) :: Fint ! excess energy one wants to compute here
  real(dp) :: Vint ! dummy part of Fint
  real(dp) :: S_cg, F_cg ! surface part of the energy (Cahn Hilliard)
  complex(dp), allocatable, dimension(:,:,:) :: dS_cgk ! coarse grained part of dF in k space
  complex(dp), allocatable, dimension(:,:,:) :: dF_cgk ! coarse grained part of dF
  real(dp), allocatable, dimension(:,:,:) :: dF_cg ! coarse grained part of dF
  real(dp) :: time0 , time1 ! timesteps
  real(dp) :: R_cg ! radius of the coarse graining gaussian function  
  real(dp) :: dF_cg_ijk ! dummy for local dF_cg
  real(dp) :: delta_n_ijk ! dummy for local delta_n 
  real(dp):: psi ! local cg_vect
  real (dp) :: a, mm
  
  call cpu_time(time0)
  !print*, c_s(1) , c_s_hs(1)
  !a=c_s(1)-c_s_hs(1)
  print*, 'a = ' , a
  a=12.3*2.0_dp*8.1812242666880284_dp!/1.03_dp
  print*, a,12.3*2.0_dp*8.1812242666880284_dp
  mm=300.0_dp!/(n_0**2)/30.0_dp
  print*, 'Hydrophobicity Van der Waals way'
  !In this routine we use c_s for HS fluid and not the HS bridge that is why c_s_hs has to be calculated
  if (.not. allocated(c_s_hs)) then
    call cs_of_k_hard_sphere
  END IF
  if (.not. input_log('hard_sphere_fluid')) then
    print*, 'energy_hydro MUST be used with hard sphere solvent correction  TAG of hard_sphere_solute in dft.in must be T'
    stop
  END IF
  ! macroscopic water parameters, from Chandler
  mu_0 = 7.16d-4*kBT ! phenomenological potential
  R_cg=4.0_dp ! ARBITRARY gaussian radius for coarse graining
  
  !TEST ZONE START
  R_cg = R_cg
  !TEST ZONE STOP
  ! total number of kpoints and volume per kpoint
  Nk = real( nfft1 * nfft2 * nfft3 ,dp)
  DeltaVk = DeltaV / Nk!twopi**3/(Lx*ly*Lz*Nk)!
  fact_n = DeltaV * n_0 ! integration factor
  prefactor = - R_cg ** 2 / 2.0_dp ! dummy for later in the gaussian in k space g(k)=exp(prefactor*k2)
  
  ! get density from last minimization step
  allocate(delta_n(nfft1,nfft2,nfft3))
  ! HERE WE WORK WITH delta_n = delta_n normalized wrt n0 = delta_n/n0 = (n-n0)/n0
  icg=0
  do i=1,nfft1
    do j=1,nfft2
      do k=1,nfft3
        delta_n_ijk=0.0_dp ! init delta_n
        do o=1,angGrid%n_angles
          do p=1, molRotGrid%n_angles
            icg=icg+1
            delta_n_ijk = delta_n_ijk + angGrid%weight(o)*molRotGrid%weight(p) * cg_vect(icg) ** 2 ! sum over all orientations
          END DO  
        END DO
          delta_n(i,j,k)=delta_n_ijk*sym_order/(twopi*fourpi) - 1.0_dp ! normalize (n=1/fourpi int_o rho(r,o))
      END DO
    END DO
  END DO
  
  
  ! TODO Next FFT sequences can be done on multiple threads
  ! compute delta_nk which is FFT(delta_n)
  fftw3%in_forward = delta_n
  call dfftw_execute ( fftw3%plan_forward )
  allocate ( delta_nk ( nfft1 / 2 + 1 , nfft2 , nfft3 ) )
  delta_nk = fftw3%out_forward
  
  ! the coarse-grained delta_n is simply delta_n multiplied by a gaussian distribution of empirical
  allocate ( delta_nbark ( nfft1/2+1 , nfft2 , nfft3 ) )
  delta_nbark = timesExpPrefactork2 (delta_nk,prefactor)
  ! surface energy (cahn hilliard part of the intrinsic excess energy)
  S_cg = 0.0_dp
  allocate ( dS_cgk ( nfft1 / 2 + 1 , nfft2 , nfft3 ) )
  do l = 1 , nfft1 / 2 + 1
  
    ! take time reversal symetry into account small part
    if ( l >= 2 .and. l <= nfft1/2 ) then
      facsym = 2.0_dp
    ELSE
      facsym = 1.0_dp
    END IF
  
    do m = 1 , nfft2
      do p = 1 , nfft3
        k2_loc = k2 ( l , m , p )
        delta_nbark_loc = Delta_nbark ( l , m , p )
        dS_cgk ( l , m , p ) = facsym * mm*DeltaVk* k2_loc * Delta_nbark_loc*n_0**2*kBT
      END DO
    END DO
  END DO
  print *, 'S_cg in k space = ' , S_cg
  
  ! gradient of delta_nbar in kspace
  ! remember that the fourier transform of the gradient of f(r) is i*vector_k*f(k)
  allocate(gradx_delta_nbark(nfft1/2+1,nfft2,nfft3))
  allocate(grady_delta_nbark(nfft1/2+1,nfft2,nfft3))
  allocate(gradz_delta_nbark(nfft1/2+1,nfft2,nfft3))
  
  do l=1,nfft1/2+1
    gradx_delta_nbark(l,:,:) = i_complex * kx(l) * delta_nbark(l,:,:)
  END DO
  
  do l=1,nfft2
    grady_delta_nbark(:,l,:) = i_complex * ky(l) * delta_nbark(:,l,:)
  END DO
  
  do l=1,nfft3
    gradz_delta_nbark(:,:,l) = i_complex * kz(l) * delta_nbark(:,:,l)
  END DO
  
  
  
  ! intrinsic local excess energy (second order term in Fexc) in k spae
  allocate(V_nk(nfft1/2+1,nfft2,nfft3))
  do p = 1 , nfft3
    do m = 1 , nfft2
      do l = 1 , nfft1 / 2 + 1
  
        ! get index of k point wrt input/cs.in
        k_index = int ( norm_k ( l , m , p ) / delta_k ) + 1
  
        ! if index is superior to biggest kpoint in cs.in, stay with largest in cs.in
        if ( k_index > nb_k ) k_index = nb_k
  
        ! compute radial excess potential in kspace
        ! we suppose c_s(nbar)=cst=c_s(n_0) which is a crude approximation
        V_nk ( l , m , p ) = (c_s ( k_index )-c_s_hs(k_index)) * ( delta_nk ( l , m , p ) - delta_nbark ( l , m , p ) )
  
      END DO
    END DO
  END DO
  
  
  ! coarse grain V_nk
  allocate( V_nbark ( nfft1 / 2 + 1 , nfft2 , nfft3 ) )
  V_nbark = timesExpPrefactork2 (V_nk,prefactor)
  
  ! compute V_n(r) by FFT-1
  fftw3%in_backward = V_nk
  deallocate(V_nk)
  call dfftw_execute(fftw3%plan_backward)
  allocate(V_n(nfft1,nfft2,nfft3))
  V_n = fftw3%out_backward/Nk
  
  ! compute coarse grained objects in real space by FFT-1
  fftw3%in_backward = V_nbark
  deallocate(V_nbark)
  call dfftw_execute(fftw3%plan_backward)
  allocate(V_nbar(nfft1,nfft2,nfft3))
  V_nbar = fftw3%out_backward/Nk
  
  fftw3%in_backward = delta_nbark
  deallocate(delta_nbark)
  call dfftw_execute(fftw3%plan_backward)
  allocate(delta_nbar(nfft1,nfft2,nfft3))
  delta_nbar = fftw3%out_backward/Nk
  
  fftw3%in_backward = gradx_delta_nbark
  deallocate(gradx_delta_nbark)
  call dfftw_execute(fftw3%plan_backward)
  allocate(gradx_delta_nbar(nfft1,nfft2,nfft3))
  gradx_delta_nbar = fftw3%out_backward/Nk
  
  fftw3%in_backward = grady_delta_nbark
  deallocate(grady_delta_nbark)
  call dfftw_execute(fftw3%plan_backward)
  allocate(grady_delta_nbar(nfft1,nfft2,nfft3))
  grady_delta_nbar = fftw3%out_backward/Nk
  
  fftw3%in_backward = gradz_delta_nbark
  deallocate(gradz_delta_nbark)
  call dfftw_execute(fftw3%plan_backward)
  allocate(gradz_delta_nbar(nfft1,nfft2,nfft3))
  gradz_delta_nbar = fftw3%out_backward/Nk
  ! compute coarse grained free energy and gradient small part
  allocate ( dF_cg ( nfft1 , nfft2 , nfft3 ) )
  F_cg = 0.0_dp
  do i=1,nfft1
    do j=1,nfft2
      do k=1,nfft3
         dF_CG(i,j,k)=-a*delta_nbar(i,j,k)*deltaV*n_0**2*kBT
      END DO
    END DO
  END DO
!large gradient part
  S_cg = 0.5_dp*mm*DeltaV *sum( gradx_delta_nbar**2 + grady_delta_nbar**2+ gradz_delta_nbar**2 )*n_0**2*kBT
  F_cg= - 0.5_dp*a*DeltaV*sum(delta_nbar**2)*n_0**2*kBT
  deallocate(gradx_delta_nbar)
  deallocate(grady_delta_nbar)
  deallocate(gradz_delta_nbar)
  print *, 'F_CG =  ',F_CG
  print *, 'S_CG =  ',S_CG
  
  ! FFT (dF_cg) in order to work in kspace
  fftw3%in_forward = dF_cg
  call dfftw_execute(fftw3%plan_forward)
  allocate(dF_cgk(nfft1/2+1,nfft2,nfft3))
  ! define coarse grained gradients in k space
  dF_cgk = timesExpPrefactork2 ( fftw3%out_forward + dS_cgk, prefactor)
  deallocate(dS_cgk)
  
  ! get back coarse grained gradients in r space
  fftw3%in_backward = dF_cgk
  deallocate ( dF_cgk )
  call dfftw_execute ( fftw3%plan_backward ) 
  dF_cg = fftw3%out_backward/Nk 
  
  ! following count does not take into account all cases where angGrid%n_angles /=1
  if ( angGrid%n_angles /= 1 ) then
    print *, 'angGrid%n_angles /= 1 and hydrophobic part is calculated while not still implemented'
    print *, 'not possible for now. error in cs_plus_hydro.f90'
    stop
  END IF
  
  ! this SUBROUTINE only works for nb_species = 1
  if ( nb_species /= 1 ) then
    print *, 'nb_species /= 1 and this is not implemented yet in cs_plus_hydro.f90'
    stop
  END IF
  ! compute final FF and dF
  Fint = 0.0_dp
  icg = 0
  do i = 1, nfft1
    do j = 1, nfft2
      do k = 1, nfft3
        Vint = -kBT * n_0 * fact_n * ( V_n(i,j,k) - V_nbar(i,j,k) ) 
        Fint = Fint + 0.5_dp* ( Delta_n (i,j,k) - Delta_nbar (i,j,k) ) * Vint ! TODO doesnt work in case of angGrid%n_angles /= 1)
        dF_cg_ijk = dF_cg(i,j,k)
        do o = 1, angGrid%n_angles
          do p=1,molRotGrid%n_angles
            icg = icg + 1
            psi = cg_vect ( icg )
            dF (icg) = dF ( icg )+2.0_dp*psi*molRotGrid%weight(p)*angGrid%weight(o)/(fourpi*twopi/sym_order)*( Vint + dF_cg_ijk )
          END DO
        END DO
      END DO
    END DO
  END DO
  
  ! conclude
    FF = FF + Fint + F_cg + S_cg
    Fint = Fint + F_cg + S_cg
  ! warn user
    call cpu_time ( time1 )
    IF (verbose) THEN
        WRITE(*,'(''    Exces / radial     = '',f11.3,'' in '',I5,'' sec'')') Fint , NINT(time1-time0)
        WRITE(*,'(''    Exces / radial+cg  = '',f11.3,'' in '',I5,'' sec'')') F_cg + S_cg , NINT(time1-time0)
    END IF

END SUBROUTINE

! this subroutine computes the radial excess free energy + its associated hydrophobic part.
! TODO This subroutine should be merged in one way or another with cs_from_dcf
subroutine cs_plus_hydro
  use precision_kinds , only : dp , i2b
  use system , only : nfft1 , nfft2 , nfft3 , deltaV, rho_0 , nb_k , c_s , kBT , delta_k , nb_species,n_0, Lx,Ly,Lz
  use constants , only : fourpi , i_complex,twopi
  use cg , only : cg_vect , FF , dF
  use quadrature, only: sym_order,angGrid, molRotGrid
  use fft , only : in_forward , in_backward , out_forward , out_backward , plan_forward , plan_backward , norm_k,kx,ky,kz,k2,&
                    timesExpPrefactork2
  use input, only : input_log
  
  implicit none
  
  real(dp) :: d_0 ! Angstroms, from Chandler
  real(dp) :: gamma_0 ! surface tension of water (0.174kBT/A² or 0.431 KJ/A²)
  real(dp) :: mu_0 ! phenomenological potential
  
  real(dp) :: Nk ! total number of k points. Real because really often used as a divisor of a real number
  real(dp) :: deltaVk ! volume unit per k point
  real(dp) :: fact_n ! integration factor = deltaV*n_0 = deltaV*rho_0/fourpi
  
  integer(i2b) :: icg,i,j,k,o,l,m,p ! dummy for loops
  
  real(dp) :: nbar ! temp for local coarse grained n
  
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
  real(dp) :: lognbar ! dummy for log(nbar)
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
  
  call cpu_time(time0)
  if (input_log('bridge_hard_sphere')) then
    print*,'Chandler Version of Hydrophobicity is not consistent with the use of hard sphere bridge, be aware of what you are doing'
  end if
  
  ! macroscopic water parameters, from Chandler
  d_0 = 1.27_dp ! Angstroms
  gamma_0 = 0.174_dp*kBT ! surface tension of water (0.174kBT/A² or 0.431 KJ/A²)
  mu_0 = 7.16d-4*kBT ! phenomenological potential
  R_cg=4.0_dp ! ARBITRARY gaussian radius for coarse graining
  
  !TEST ZONE START
  d_0 = d_0
  gamma_0 = gamma_0/10.0_dp  !10. !10  22 is good for R=10A
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
          end do  
        end do
          delta_n(i,j,k)=delta_n_ijk*sym_order/(twopi*fourpi) - 1.0_dp ! normalize (n=1/fourpi int_o rho(r,o))
      end do
    end do
  end do
  
  
  ! TODO Next FFT sequences can be done on multiple threads
  ! compute delta_nk which is FFT(delta_n)
  in_forward = delta_n
  call dfftw_execute ( plan_forward )
  allocate ( delta_nk ( nfft1 / 2 + 1 , nfft2 , nfft3 ) )
  delta_nk = out_forward
  
  ! the coarse-grained delta_n is simply delta_n multiplied by a gaussian distribution of empirical
    allocate ( delta_nbark ( nfft1/2+1 , nfft2 , nfft3 ) )
    delta_nbark = timesExpPrefactork2 ( delta_nk,prefactor ) ! ~n=n*exp(prefactor*k²)
  ! surface energy (cahn hilliard part of the intrinsic excess energy)
  S_cg = 0.0_dp
  allocate ( dS_cgk ( nfft1 / 2 + 1 , nfft2 , nfft3 ) )
  do l = 1 , nfft1 / 2 + 1
  
    ! take time reversal symetry into account small part
    if ( l >= 2 .and. l <= nfft1/2 ) then
      facsym = 2.0_dp
    else
      facsym = 1.0_dp
    end if
  
    do m = 1 , nfft2
      do p = 1 , nfft3
        k2_loc = k2 (l,m,p)
        delta_nbark_loc = Delta_nbark (l,m,p)
        S_cg = S_cg + facsym * DeltaVk * 1.5_dp * d_0 * gamma_0 * k2_loc * real( Delta_nbark_loc , dp ) ** 2
        dS_cgk ( l , m , p ) = facsym * DeltaVk * 3.0 * d_0 * gamma_0 * k2_loc * Delta_nbark_loc
      end do
    end do
  end do
  print *, 'S_cg in k space = ' , S_cg
  
  ! gradient of delta_nbar in kspace
  ! remember that the fourier transform of the gradient of f(r) is i*vector_k*f(k)
  allocate(gradx_delta_nbark(nfft1/2+1,nfft2,nfft3))
  allocate(grady_delta_nbark(nfft1/2+1,nfft2,nfft3))
  allocate(gradz_delta_nbark(nfft1/2+1,nfft2,nfft3))
  
  do l=1,nfft1/2+1
    gradx_delta_nbark(l,:,:) = i_complex * kx(l) * delta_nbark(l,:,:)
  end do
  
  do l=1,nfft2
    grady_delta_nbark(:,l,:) = i_complex * ky(l) * delta_nbark(:,l,:)
  end do
  
  do l=1,nfft3
    gradz_delta_nbark(:,:,l) = i_complex * kz(l) * delta_nbark(:,:,l)
  end do
  
  
  
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
        V_nk ( l , m , p ) = c_s ( k_index ) * ( delta_nk ( l , m , p ) - delta_nbark ( l , m , p ) )
  
      end do
    end do
  end do
  
  
  ! coarse grain V_nk
    ALLOCATE( V_nbark ( nfft1 / 2 + 1 , nfft2 , nfft3 ) )
    V_nbark = timesExpPrefactork2 (V_nk,prefactor)
  
  ! compute V_n(r) by FFT-1
  in_backward = V_nk
  deallocate(V_nk)
  call dfftw_execute(plan_backward)
  allocate(V_n(nfft1,nfft2,nfft3))
  V_n = out_backward/Nk
  
  ! compute coarse grained objects in real space by FFT-1
  in_backward = V_nbark
  deallocate(V_nbark)
  call dfftw_execute(plan_backward)
  allocate(V_nbar(nfft1,nfft2,nfft3))
  V_nbar = out_backward/Nk
  
  in_backward = delta_nbark
  deallocate(delta_nbark)
  call dfftw_execute(plan_backward)
  allocate(delta_nbar(nfft1,nfft2,nfft3))
  delta_nbar = out_backward/Nk
  
  in_backward = gradx_delta_nbark
  deallocate(gradx_delta_nbark)
  call dfftw_execute(plan_backward)
  allocate(gradx_delta_nbar(nfft1,nfft2,nfft3))
  gradx_delta_nbar = out_backward/Nk
  
  in_backward = grady_delta_nbark
  deallocate(grady_delta_nbark)
  call dfftw_execute(plan_backward)
  allocate(grady_delta_nbar(nfft1,nfft2,nfft3))
  grady_delta_nbar = out_backward/Nk
  
  in_backward = gradz_delta_nbark
  deallocate(gradz_delta_nbark)
  call dfftw_execute(plan_backward)
  allocate(gradz_delta_nbar(nfft1,nfft2,nfft3))
  gradz_delta_nbar = out_backward/Nk
  ! compute coarse grained free energy and gradient small part
  allocate ( dF_cg ( nfft1 , nfft2 , nfft3 ) )
  F_cg = 0.0_dp
  do i=1,nfft1
    do j=1,nfft2
      do k=1,nfft3
  
        nbar = delta_nbar(i,j,k) + 1.0_dp
        if (nbar==0.0) stop 'problem in cs_plus_hydro.f90 : nbar should not be =0.0. STOP'
        lognbar = log(nbar)
  
        F_cg = F_cg + DeltaV*(6.0*gamma_0/d_0)*(nbar-1.0_dp)**2*nbar**2 - mu_0*n_0*(nbar-1.0_dp)*deltaV
  
        ! + ideal part associated to nbar
        F_cg = F_cg - kBT*fact_n*(nbar*lognbar- nbar + 1.0_dp)
  
        dF_cg(i,j,k) = (DeltaV*(6.0*gamma_0/d_0)*2.0_dp*(nbar-1.0_dp)*nbar*(2.0_dp*nbar-1.0_dp) - mu_0*n_0*deltaV)&
                        - kBT*fact_n*lognbar
      end do
    end do
  end do
!large gradient part
  S_cg = 3.0_dp/2.0_dp*d_0*gamma_0*DeltaV *sum( gradx_delta_nbar**2 + grady_delta_nbar**2&
                  + gradz_delta_nbar**2 )
  deallocate(gradx_delta_nbar)
  deallocate(grady_delta_nbar)
  deallocate(gradz_delta_nbar)
  print *, 'F_CG =  ',F_CG
  print *, 'S_CG =  ',S_CG
  
  ! FFT (dF_cg) in order to work in kspace
  in_forward = dF_cg
  call dfftw_execute(plan_forward)
  allocate(dF_cgk(nfft1/2+1,nfft2,nfft3))
  dF_cgk = out_forward
  
  ! define coarse grained gradients in k space
  dF_cgk = timesExpPrefactork2 ( dF_cgk + dS_cgk, prefactor )
  deallocate(dS_cgk)
  
  ! get back coarse grained gradients in r space
  in_backward = dF_cgk
  deallocate ( dF_cgk )
  call dfftw_execute ( plan_backward ) 
  dF_cg = out_backward/Nk 
  
  ! following count does not take into account all cases where angGrid%n_angles /=1
  if ( angGrid%n_angles /= 1 ) then
    print *, 'angGrid%n_angles /= 1 and hydrophobic part is calculated while not still implemented'
    print *, 'not possible for now. error in cs_plus_hydro.f90'
    stop
  end if
  
  ! this subroutine only works for nb_species = 1
  if ( nb_species /= 1 ) then
    print *, 'nb_species /= 1 and this is not implemented yet in cs_plus_hydro.f90'
    stop
  end if
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
          end do
        end do
      end do
    end do
  end do
  
  ! conclude
  FF = FF + Fint + F_cg + S_cg
  
  ! warn user
  call cpu_time ( time1 )
  write(*,*) 'Fexc(rad)   = ' , Fint , 'computed in (sec)'
  write(*,*) 'Fexc(rad+cg)= ' , F_cg + S_cg , 'computed in (sec)' , time1 - time0
end subroutine cs_plus_hydro

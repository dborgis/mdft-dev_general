!===================================================================================================================================
module energy_ck_en_harm_sph
!-----------------------------------------------------------------------------------------------------------------------------------
  use precision_kinds, only : i2b,dp
  use minimizer,       only : cg_vect,FF,dF
  implicit none
  complex(dp), dimension(:,:,:,:), allocatable :: proj_tot,proj_test    ! i'm not the only one

  contains


subroutine compute_Fexc_ck_en_harm_sph_asymm(Fexc_ck_proj) ! test usage

    use rotational_invariants, only : num_proj,num_proj2,projInd2,gsh,rotM_gen,fgsht_fwd,fgsht_bwd,nmup0
    use quadrature,            only : mmax,mrso=>molRotSymOrder,intScheme,angGrid,molRotGrid
    use dcf,                   only : nb_k,delta_k,ck_proj_chi
    use fft,                   only : fft3c,kx,ky,kz
    use system,                only : kBT,spaceGrid,rho_0,n_0
    use constants,             only : twopi,eightpiSQ
    implicit none

    real(dp),  intent(out)  :: Fexc_ck_proj ! excess part of the free energy
    integer(i2b)            :: nfft1,nfft2,nfft3,nk,num_cos,num_phi,num_psi
    integer(i2b)            :: i,j,k,l,m,n,o1,o2,p,icg,ik,a1,a2,ic,nmm,nmup,nmu,nchi,nnn,nnu,counter,ntheta,nphi,npsi
    real(dp)                :: den,local_weight, time1, time2
    real(dp), dimension(3)  :: q
    complex(dp),  parameter :: imag=(0._dp,1._dp)

    real(dp),    dimension(mmax/mrso*2+1,mmax*2+1,mmax+1)          :: foo_temp              ! delta_rho_temp
    complex(dp), dimension(mmax/mrso*2+1,mmax*2+1,mmax+1)          :: gamma_temp            ! gamma_temp
    complex(dp), dimension(num_proj2)                              :: foo_proj2
    complex(dp), dimension(-mmax/mrso:mmax/mrso,-mmax:mmax,0:mmax) :: foo_proj              ! delta_rho_proj & gamma_proj
    complex(dp), dimension(num_proj)                               :: proj_temp,proj_temp_p ! temporary full proj array
    complex(dp), dimension(0:mmax,-mmax:mmax,-mmax:mmax)           :: rotM
    real(dp)                :: theta(mmax+1),phi(mmax*2+1),psi(mmax/mrso*2+1)
    integer(i2b)            :: projInd(-mmax/mrso:mmax/mrso,-mmax:mmax,0:mmax),projIndChi(-mmax/mrso:mmax/mrso,0:mmax,-mmax:mmax)
    integer(i2b), parameter :: ltest=2,mtest=3,ntest=5 ! to be modified
    character(33)           :: filetemp

  ! initialisation .................................................................................................................

  call cpu_time(time1)

    nfft1 = spaceGrid%n_nodes(1); nfft2 = spaceGrid%n_nodes(2); nfft3 = spaceGrid%n_nodes(3); nk = nfft1*nfft2*nfft3
    num_cos = mmax+1; num_phi = mmax*2+1; num_psi = mmax/mrso*2+1
    Fexc_ck_proj = 0._dp
    if (.not. allocated(proj_tot)) allocate(proj_tot(num_proj,nfft3,nfft2,nfft1),source=(0._dp,0._dp))
  !  allocate(proj_test(num_proj,nfft3,nfft2,nfft1),source=(0._dp,0._dp))

  ! tabulate theta,phi,psi
    do ic=1,mmax+1;        theta(ic) = acos(intScheme%root(ic));                   end do
    do ic=1,mmax*2+1;      phi(ic)   = real(ic-1,dp)*twopi/real(mmax*2+1,dp);      end do
    do ic=1,mmax/mrso*2+1; psi(ic)   = real(ic-1,dp)*twopi/real(mmax/mrso*2+1,dp); end do ! twopi to adapt the nmu index

  ! tabulate indices
    counter = 0
    do nmm=0,mmax; do nmup=-nmm,nmm; do nmu=-nmm/mrso,nmm/mrso
      counter = counter + 1
      projInd(nmu,nmup,nmm) = counter
    end do; end do; end do

    counter = 0
    do nchi=-mmax,mmax; do nmm=abs(nchi),mmax; do nmu=-nmm/mrso,nmm/mrso
        counter = counter + 1
        projIndChi(nmu,nmm,nchi) = counter
    end do; end do; end do

  ! expended rho onto generalized spherical harmonics ..............................................................................
    icg = 0
    do i=1,nfft1; do j=1,nfft2; do k=1,nfft3

    ! copy rho
      do o1=1,num_cos; do o2=1,num_phi; do p=1,num_psi
        icg = icg + 1
        foo_temp(p,o2,o1) = cg_vect(icg) ** 2 - 1._dp
      end do; end do; end do


    ! fgsht-forward
      call fgsht_fwd(foo_temp,foo_proj2)
      foo_proj2 = foo_proj2


    ! copy foo_proj2 to proj_tot
      do nmm=0,mmax; do nmup=0,2*nmm
        do nmu=-nmm/mrso,-1
          a1 = projInd(nmu,nmup0(nmup,nmm),nmm)
          a2 = projInd2(-nmu,mod(2*nmm+1-nmup,2*nmm+1),nmm)
          proj_tot(a1,k,j,i) = (-1)**(nmup0(nmup,nmm)+nmu*mrso)*conjg(foo_proj2(a2))
        end do
        do nmu=0,nmm/mrso
          a1 = projInd(nmu,nmup0(nmup,nmm),nmm)
          a2 = projInd2(nmu,nmup,nmm)
          proj_tot(a1,k,j,i) = foo_proj2(a2)
        end do
      end do; end do

    end do; end do; end do

  ! fourier transform of proj_tot ..................................................................................................
    do counter=1,num_proj
      fft3c%in_forward = proj_tot(counter,:,:,:) ! (k,j,i)
      call dfftw_execute(fft3c%plan_forward)
      proj_tot(counter,:,:,:) = fft3c%out_forward
    end do


  ! open loop k ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    do l=1,nfft1; do m=1,nfft2; do n=1,nfft3
      q(1) = kx(l); q(2) = ky(m); q(3) = kz(n)
      ik = int(norm2(q) / delta_k + 0.5_dp) + 1 ! index k


  ! generate rotM(nmm,nmup,nchi) for (l,m,n) .......................................................................................
      rotM = rotM_gen(mmax,q)

  ! rotated into local coordinates system ..........................................................................................
      proj_temp = (0._dp,0._dp)
      do nchi=-mmax,mmax; do nmm=abs(nchi),mmax; do nmup=-nmm,nmm; do nmu=-nmm/mrso,nmm/mrso
        a1 = projInd(nmu,nmup,nmm)
        a2 = projIndChi(nmu,nmm,nchi)
        proj_temp(a2) = proj_temp(a2) + proj_tot(a1,n,m,l) * rotM(nmm,nmup,nchi)
      end do; end do; end do; end do


      ic = 0
      proj_temp_p = (0._dp,0._dp)
      do nchi=-mmax,mmax; do nnn=abs(nchi),mmax; do nnu=-nnn/mrso,nnn/mrso; do nmm=abs(nchi),mmax; do nmu=-nmm/mrso,nmm/mrso
        ic = ic + 1
        a1 = projIndChi(nmu,nmm,nchi)
        a2 = projIndChi(-nnu,nnn,nchi)
        proj_temp_p(a1) = proj_temp_p(a1) + ck_proj_chi(ic,ik) * proj_temp(a2)
      end do; end do; end do; end do; end do


  ! rotate gamma back to global coordinates system .................................................................................
      proj_tot(:,n,m,l) = (0._dp,0._dp)
      do nchi=-mmax,mmax; do nmm=abs(nchi),mmax; do nmup=-nmm,nmm; do nmu=-nmm/mrso,nmm/mrso
        a1 = projIndChi(nmu,nmm,nchi)
        a2 = projInd(nmu,nmup,nmm)
        proj_tot(a2,n,m,l) = proj_tot(a2,n,m,l) + proj_temp_p(a1) * conjg(rotM(nmm,nmup,nchi))
      end do; end do; end do; end do

    end do; end do; end do !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


  ! inverse fourier transform of gamma .............................................................................................
    do counter=1,num_proj
      fft3c%in_backward = proj_tot(counter,:,:,:) ! (k,j,i)
      call dfftw_execute(fft3c%plan_backward)
      proj_tot(counter,:,:,:) = fft3c%out_backward
    end do
    proj_tot = proj_tot / nk

  ! rebuild gamma in angular frame / calculation of excess part of the free energy functional ......................................
    icg = 0
    do i=1,nfft1; do j=1,nfft2; do k=1,nfft3


    ! copy proj_tot to foo_proj2
      do nmm=0,mmax; do nmup=0,2*nmm
        do nmu=0,nmm/mrso
          a1 = projInd(nmu,nmup0(nmup,nmm),nmm)
          a2 = projInd2(nmu,nmup,nmm)
          foo_proj2(a2) = proj_tot(a1,k,j,i)
        end do
      end do; end do

    ! fgsht-backward
      call fgsht_bwd(foo_proj2,foo_temp)

    ! calculate FF and dF
      foo_temp = foo_temp * (-kBT) * rho_0 ! gamma
      do o1=1,num_cos; do o2=1,num_phi; do p=1,num_psi
        icg = icg + 1
        den = cg_vect(icg)
        local_weight = angGrid%weight((o1-1)*num_phi+o2) * molRotGrid%weight(p)
        Fexc_ck_proj = Fexc_ck_proj + foo_temp(p,o2,o1) * (den**2-1) * local_weight * eightpiSQ / mrso ! Fexc
        dF(icg) = dF(icg) + 2._dp * den * foo_temp(p,o2,o1) * spaceGrid%dv * rho_0 * local_weight * eightpiSQ / mrso ! dF
      end do; end do; end do
    end do; end do; end do

    Fexc_ck_proj = Fexc_ck_proj / 2._dp * spaceGrid%dv * rho_0
    FF = FF + Fexc_ck_proj

    call cpu_time(time2)
  !  write(*,*) 'tous les q, time= ',time2 - time1
  !.................................................................................................................................
  end subroutine compute_Fexc_ck_en_harm_sph_asymm




  integer function mq(id,iq) ! find the index of -q for each dimension id

    use precision_kinds, only : i2b
    use system,          only : spaceGrid
    implicit none
    integer(i2b), intent(in) :: id,iq

    if (iq==1) then
      mq = 1
    else
      mq = spaceGrid%n_nodes(id) + 2 - iq
    end if

  end function mq
!-----------------------------------------------------------------------------------------------------------------------------------
end module energy_ck_en_harm_sph
!===================================================================================================================================

subroutine energy_polarization
use precision_kinds , only : i2b , dp
use system , only : nfft1 , nfft2 , nfft3 , nb_omega , Lx , Ly , Lz , c_delta , c_d , kBT , rho_0 , delta_k , nb_k ,&
                   deltav, nb_psi, deltax,deltay,deltaz
use quadrature , only : weight , Omx , Omy , Omz, sym_order , weight_psi
use cg , only : cg_vect , FF , dF
use constants , only : twopi
use fft , only : in_forward , in_backward , out_forward , out_backward , plan_forward , plan_backward
use input , only : input_line, input_log,input_char
implicit none
integer(i2b):: icg , i , j , k , l , m , n , m1 , m2 , m3 , o , p!> Dummy
integer(i2b):: nf1 , nf2 , nf3 ! dummy nfft1/2 , nfft2/2 , nfft3/2
integer(i2b):: k_index
real(dp):: kx , ky , kz , kx2 , ky2 , kz2 , k2 , norm_k
real(dp):: Fint !> Internal part of the free energy due to polarization
real(dp):: Vint !> Dummy for calculation of Vint
real(dp):: fact !> facteur d'integration
real(dp):: rho , psi !> Dummy
real(dp), allocatable , dimension ( : , : , : ) :: Px , Py , Pz , Ex , Ey , Ez
real(dp), allocatable , dimension ( : , : , : , : ) :: polatot, polascal
complex(dp), allocatable , dimension ( : , : , : ) :: Pkx , Pky , Pkz , Ekx , Eky , Ekz
complex(dp):: k_dot_P
real(dp):: c_deltat , c_dt ! dummy local values of c_delta and c_d in loops
real(dp):: time1 , time0 , time2, time3! timestamps
real(dp):: twopioLx , twopioLy , twopioLz ! dummy for 2pi/Lx, 2pi/Ly, 2pi/Lz
real(dp):: pxt , pyt , pzt , r! dummy
complex(dp):: pxt_k , pyt_k , pzt_k ! dummy
real(dp):: Nk ! total number of k grid points
real(dp), allocatable , dimension ( : ) :: weight_omx , weight_omy , weight_omz ! dummy
character (50) :: filename
! look for tag polarization in input
if (.not. input_log('polarization')) then
return
end if
!Check if you want to compute Polarization from a macroscopic point of view
!do i = 1 , size ( input_line )
!  j = len ( 'evaluate_polarization' )
!  if ( input_line (i) (1:j) == 'evaluate_polarization' .and. input_line (i) (j+4:j+8) == 'Micro' ) return ! exit this routine to get back to energy calculation skeleton.  
!end do
if (trim(adjustl(input_char('evaluate_polarization')))== 'multi') return
! init timer
call cpu_time(time0)
! init variables
nf1 = nfft1 / 2
nf2 = nfft2 / 2
nf3 = nfft3 / 2
twopioLx = twopi / Lx
twopioLy = twopi / Ly
twopioLz = twopi / Lz
! total number of k grid points
Nk = real ( nfft1 * nfft2 * nfft3 , dp )
! allocate and init polarization vector
allocate ( Px ( nfft1 , nfft2 , nfft3 ) )
allocate ( Py ( nfft1 , nfft2 , nfft3 ) )
allocate ( Pz ( nfft1 , nfft2 , nfft3 ) )
allocate ( Polatot ( nfft1 , nfft2 , nfft3 , 1) )
Polatot=0.0_dp
allocate ( Polascal ( nfft1 , nfft2 , nfft3 , 1) )
Polascal=0.0_dp
Px = 0.0_dp
Py = 0.0_dp
Pz = 0.0_dp
! put density of last minimization step in delta_rho and P
! but first prepare the product weight(omega)*Omx in order not to repeat it
allocate ( weight_omx ( nb_omega ) )
allocate ( weight_omy ( nb_omega ) )
allocate ( weight_omz ( nb_omega ) )
weight_omx = weight * Omx
weight_omy = weight * Omy
weight_omz = weight * Omz
icg = 0
do i = 1 , nfft1
    
    m1=i-1
    if (i> nfft1/2) m1=i-1-nfft1
  do j = 1 , nfft2
    m2=j-1
    if (i>nfft2/2) m2=j-1-nfft2
    do k = 1 , nfft3    
    m3=k-1
    if (k>nfft3/2) m3=k-1-nfft3
      ! init dummy variables tpx , tpy and tpz in order not to loop directly over big arrays
      pxt = 0.0_dp
      pyt = 0.0_dp
      pzt = 0.0_dp
      do o = 1 , nb_omega
        do p=1, nb_psi
          icg = icg + 1
          rho = cg_vect (icg) ** 2
          pxt = pxt + weight_Omx ( o ) * weight_psi(p) * rho
          pyt = pyt + weight_Omy ( o ) * weight_psi(p) * rho
          pzt = pzt + weight_Omz ( o ) * weight_psi(p) * rho
          r=sqrt((m1*deltax)**2+(m2*deltay)**2+(m3*deltaz)**2)
        end do
      end do
      Px ( i , j , k ) = pxt
      Py ( i , j , k ) = pyt
      Pz ( i , j , k ) = pzt
      polascal(i,j,k,1)=polascal(i,j,k,1)+pxt*m1*deltax/r+pyt*m2*deltay/r+pzt*m3*deltaz/r
    end do
  end do
end do
open(11,file='output/polatotxmax')
do i=1,nfft1
write(11,*) i*deltax, 0.4894_dp*Px(i,nfft2/2+1,nfft3/2+1)
end do
close(11)
!Compute Radial Polarization
filename='output/radial_polarization_dipolar'
polatot(:,:,:,1)=sqrt(0.4894_dp**2*(Px(:,:,:)**2+Py(:,:,:)**2+Pz(:,:,:)**2))*rho_0
call compute_rdf(polatot, filename)
filename='output/radial_polarization_scalar'
call compute_rdf(polascal*rho_0*0.4894_dp, filename)
! deallocate now useless dummy arrays
deallocate ( weight_omx )
deallocate ( weight_omy )
deallocate ( weight_omz )
! fourier transform px , py and pz
allocate ( Pkx ( nfft1 / 2 + 1 , nfft2 , nfft3 ) )
allocate ( Pky ( nfft1 / 2 + 1 , nfft2 , nfft3 ) )
allocate ( Pkz ( nfft1 / 2 + 1 , nfft2 , nfft3 ) )
call cpu_time(time2)
in_forward = Px
call dfftw_execute ( plan_forward )
Pkx = out_forward
call cpu_time(time3)
print*, 'Duree FFT=', time3-time2
in_forward = Py
call dfftw_execute ( plan_forward )
Pky = out_forward
in_forward = Pz
call dfftw_execute ( plan_forward )
! deallocate useless arrays
deallocate ( Px )
deallocate ( Py )
deallocate ( Pz )
Pkz = out_forward 
! compute polarisation in k-space
allocate ( Ekx ( nfft1 / 2 + 1 , nfft2 , nfft3 ) )
allocate ( Eky ( nfft1 / 2 + 1 , nfft2 , nfft3 ) )
allocate ( Ekz ( nfft1 / 2 + 1 , nfft2 , nfft3 ) )
! get maximum number of k points as inputed in c_delta and c_d
do l = 1 , nf1 + 1
  m1 = l - 1
  if ( l > nf1 ) m1 = l - 1 - nfft1
  kx = twopioLx * real ( m1 , dp )
  kx2 = kx ** 2
  
  do m = 1 , nfft2
    m2 = m - 1
    if ( m > nf2 ) m2 = m - 1 - nfft2
    ky = twopioLy * real ( m2 , dp )
    ky2 = ky ** 2
    do n = 1 , nfft3
      m3 = n - 1
      if ( n > nf3 ) m3 = n - 1 - nfft3
      kz = twopioLz * real ( m3 , dp )
      kz2 = kz ** 2
      ! squared norm of k vector
      k2 = kx2 + ky2 + kz2
  
      norm_k = sqrt ( k2 )
      ! index of k number
      k_index = int ( norm_k / delta_k ) + 1
      ! Here it happens that k_index gets higher than the highest c_k index.
      ! In this case one imposes k_index = k_index_max
      if ( k_index > nb_k ) k_index = nb_k
      ! use dummy pxt_k instead of Pkx to reduce acces to big arrays
      pxt_k = Pkx ( l , m , n )
      pyt_k = Pky ( l , m , n )
      pzt_k = Pkz ( l , m , n )
      ! pay attention to first moment value in order not to divise by zero
      if ( k2 /= 0.0_dp) then
        k_dot_P = ( kx * pxt_k + ky * pyt_k + kz * pzt_k ) / k2
      else
        k_dot_P = cmplx ( tiny ( 1.0_dp ) , tiny ( 1.0_dp ) , dp )
      end if
      c_deltat = c_delta ( k_index )
      c_dt = c_d ( k_index )
      Ekx ( l , m , n ) = c_deltat * pxt_k + c_dt * ( 3.0_dp * k_dot_P * kx - pxt_k )
      Eky ( l , m , n ) = c_deltat * pyt_k + c_dt * ( 3.0_dp * k_dot_P * ky - pyt_k )
      Ekz ( l , m , n ) = c_deltat * pzt_k + c_dt * ( 3.0_dp * k_dot_P * kz - pzt_k )
    end do
  end do
end do
! deallocate useless Pkx , Pky and Pkz
deallocate ( Pkx )
deallocate ( Pky )
deallocate ( Pkz )
! inverse fourier transform the polarization field
! next inverse fourier transform sequence could be done in parallel
allocate ( Ex ( nfft1 , nfft2 , nfft3 ) )
allocate ( Ey ( nfft1 , nfft2 , nfft3 ) )
allocate ( Ez ( nfft1 , nfft2 , nfft3 ) )
in_backward = Ekx
call dfftw_execute ( plan_backward )
Ex = out_backward / Nk
in_backward = Eky
call dfftw_execute ( plan_backward )
Ey = out_backward / Nk 
in_backward = Ekz
call dfftw_execute ( plan_backward )
Ez = out_backward / Nk 
deallocate(Ekx)
deallocate(Eky)
deallocate(Ekz)
! integration factor
fact = deltav * rho_0
! init energy and gradient
Fint = 0.0_dp
icg = 0
do i = 1 , nfft1
  do j = 1 , nfft2
    do k = 1 , nfft3
      do o = 1 , nb_omega
        do p=1 , nb_psi
          icg = icg + 1
          psi = cg_vect ( icg )
          rho = psi ** 2
          Vint = - kBT * rho_0 * ( Omx ( o ) * Ex ( i ,  j , k ) + Omy ( o ) * Ey ( i , j , k ) + Omz ( o ) * Ez ( i , j , k ) )
          Fint = Fint + (rho - 1.0_dp) * weight(o) * weight_psi(p) * Vint
          dF (icg) = dF ( icg ) + 2.0_dp * psi * weight(o) * weight_psi(p) * fact * Vint
        end do  
      end do
    end do
  end do
end do
Fint = Fint * 0.5_dp * deltav * rho_0
FF = FF + Fint
! deallocate useless arrays
deallocate ( Ex )
deallocate ( Ey )
deallocate ( Ez )
! stop timer
call cpu_time ( time1 )
! inform user
write (*,*) 'Fexc polar  = ' , Fint , 'computed in (sec)' , time1 - time0
end subroutine energy_polarization

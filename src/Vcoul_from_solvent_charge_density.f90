 !In this routine we get the electrostatic potential generated by the solvent charge density
subroutine Vcoul_from_solvent_charge_density
use precision_kinds , only : i2b , dp
use input , only : input_line
use cg , only : FF , dF , cg_vect
use fft , only : in_forward , in_backward , out_forward , out_backward , plan_forward , plan_backward
use system, only : sigma_k, nfft1, nfft2, nfft3, Lx, Ly, Lz, deltax, rho_c, deltaV, nb_species,&
                   rho_0_multispec , rho_c_k_myway
use constants, only : twopi, fourpi, qfact
use quadrature, only : sym_order, angGrid, molRotGrid
implicit none
integer(i2b):: i , j , k , o , p, icg , m1, m2, m3, nf1, nf2, nf3, species, i1, i2
! dummy
real (dp) :: kx, ky, kz, kx2, ky2, kz2, k2
real (dp) :: twopioLx , twopioLy , twopioLz, deltaVk, Vcoul
complex(dp), dimension ( nfft1 / 2 + 1 , nfft2 , nfft3 ) :: V_c_k
real (dp) , allocatable , dimension ( : , : , :  ) ::  v_c2, rho_c_solv
real (dp) , allocatable , dimension ( : , : , :  , :,:,: ) :: rho
complex (dp) , allocatable , dimension ( : , : , : , :, : , :) :: rho_k
real(dp):: rho_temp
! init total energy and gradient to 0
! FF is the TOTAL ENERGY of the system, it is thus the functional of the density that is minimized by solver
! dF is the gradient of FF with respect to all coordinates. Remember it is of the kind dF ( number of variables over density (ie angles etc))
FF = 0.0_dp ! scalar
dF = 0.0_dp ! array
if (.not. allocated (rho_c_k_myway) ) allocate (rho_c_k_myway(nfft1/2+1, nfft2, nfft3))
rho_c_k_myway=(0.0_dp, 0.0_dp)
icg=0
twopioLx=twopi/Lx
twopioLy=twopi/Ly
twopioLz=twopi/Lz
nf1=nfft1/2
nf2=nfft2/2
nf3=nfft3/2
deltaVk=twopi**3/(Lx*Ly*Lz)
Vcoul=0.0_dp
allocate ( rho_k (nf1+1, nfft2 , nfft3, angGrid%n_angles , molRotGrid%n_angles, nb_species ))
rho_k = 0.0_dp
allocate(v_c2(nfft1, nfft2, nfft3))
allocate(rho_c_solv(nfft1, nfft2, nfft3))
rho_c_solv=0.0_dp
allocate(rho (nfft1, nfft2, nfft3,angGrid%n_angles,molRotGrid%n_angles, nb_species ))
! compute part of the energy due to the external potential
icg = 0
!check if we want to compute electrostatic using this way
do i1 = 1 , size ( input_line )
  i2 = len ( 'evaluate_polarization' )
  if ( input_line (i1) (1:i2) == 'evaluate_polarization' .and. input_line (i1) (i2+4:i2+8) == 'Macro' ) then
    return ! it is useless to continue the loop over i if you've found the tag
  end if
end do
!Compute rho_k
do species =1 , nb_species
   DO i=1,nfft1
     do j=1, nfft2
       do k=1, nfft3
        do o = 1 , angGrid%n_angles
         
          do p=1, molRotGrid%n_angles
            icg = icg + 1
            rho_temp = cg_vect ( icg ) ** 2 
            rho ( i , j , k , o , p , species ) = rho_temp*rho_0_multispec ( species ) 
          end do
        end do
      end do
    end do
  end do
end do
!Fourier transform backward solvent density
do species = 1 , nb_species
   do o =1, angGrid%n_angles
     do p=1, molRotGrid%n_angles
      in_forward = rho ( : , : , : , o , p , species )
      call dfftw_execute ( plan_forward )
      rho_k ( : , : , : , o, p ,species ) = out_forward*deltaV
     end do
   end do
end do
!Compute rho_c(k)
do species = 1, nb_species
do i = 1 , nf1 + 1
  m1 = i - 1
  if ( i > nf1 ) m1 = i - 1 - nfft1
  kx = twopioLx * real ( m1 , dp )
  kx2 = kx ** 2
  
  do j = 1 , nfft2
    m2 = j - 1
    if ( j > nf2 ) m2 = j - 1 - nfft2
    ky = twopioLy * real ( m2 , dp )
    ky2 = ky ** 2
    do k = 1 , nfft3
      m3 = k - 1
      if ( k > nf3 ) m3 = k - 1 - nfft3
      kz = twopioLz * real ( m3 , dp )
      kz2 = kz ** 2
      ! squared norm of k vector
         do o = 1 , angGrid%n_angles
            do p=1 , molRotGrid%n_angles
              rho_c_k_myway(i,j,k)=rho_c_k_myway(i,j,k)+angGrid%weight(o)*molRotGrid%weight(p)*sigma_k(i,j,k,o,p,species)*&
              rho_k (i , j , k , o , p , species)!*deltavK
         
          end do  !p
        end do   !o
     end do  !k
   end do  !j
end do  !i
end do
!compute V_c(k)
do i = 1 , nf1 + 1
  m1 = i - 1
  if ( i > nf1 ) m1 = i - 1 - nfft1
  kx = twopioLx * real ( m1 , dp )
  kx2 = kx ** 2
  
  do j = 1 , nfft2
    m2 = j - 1
    if ( j > nf2 ) m2 = j - 1 - nfft2
    ky = twopioLy * real ( m2 , dp )
    ky2 = ky ** 2
    do k = 1 , nfft3
      m3 = k - 1
      if ( k > nf3 ) m3 = k - 1 - nfft3
      kz = twopioLz * real ( m3 , dp )
      kz2 = kz ** 2
      ! squared norm of k vector
      k2 = kx2 + ky2 + kz2
    
          if (k2/=0.0_dp) then
          V_c_k(i,j,k) = rho_c_k_myway(i,j,k) * fourpi / k2 
          else
          V_c_k(i,j,k) = 0.0_dp
          end if
        
     end do  !k
   end do  !j
end do  !i
! get real space potential V(r)
if ( .not. allocated ( V_c2 ) ) allocate ( V_c2 ( nfft1 , nfft2 , nfft3 ) )
in_backward = V_c_k
call dfftw_execute ( plan_backward )
V_c2 = out_backward*deltaVk/(twopi)**3 
!Compute electrostatic energy
do i =1, nfft1
   do j=1, nfft2
     do k=1, nfft3
       Vcoul=Vcoul+rho_c( i , j , k )*V_c2( i , j , k  )*qfact*deltaV
     end do
   end do
end do
Print*, Vcoul , ' Vcoul calculé depuis l espace des k'
!Trace V_c(r) 
open(11, file='output/v_c2.dat')
do i=1, nfft1
m1=i-1
  if (m1>nf1) then
  m1= i-1-nfft1
  end if
write(11,*) i*deltax, V_c2( i , j , k  )*qfact*deltaV
end do
print*, 'v_c2 written'
!Trace V_c(r) selon x
open(17, file='output/V_c_along_z.dat')
do i=1, nfft1
write(17,*), i*deltax, V_c2(i,nf2,nf3)
end do
close(17)
end subroutine

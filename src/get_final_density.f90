!> Gets the final density from the last minimizer step.
subroutine get_final_density ( neq )
use precision_kinds , only: dp , i2b
use system , only : nfft1 , nfft2 , nfft3 , nb_species , Lx , Ly , Lz , mole_fraction, n_0_multispec, nb_psi,deltax,&
deltay,deltaz
use constants , only : fourpi, pi , twopi
use cg , only : CG_vect
use quadrature , only : weight, sym_order, weight_psi, angGrid
use fft , only : in_forward , out_forward , in_backward , out_backward , plan_forward , plan_backward , k2
implicit none
real(dp), intent(out) , dimension ( nfft1 , nfft2 , nfft3 , nb_species ) :: neq !> equilibrium density(position)
integer(i2b):: i , j , k , omega , icg , species , p ! dummy
real(dp):: rho_over_fourpi !> = CG_vect(i)**2/fourpi
real(dp):: local_density ! dummy for speeding loops
complex(dp), allocatable , dimension ( : , : , : , : ) :: rho_k 
!integer(i2b):: m1 , m2 , m3 , nf1 , l , m , n , nf2 , nf3
!real  ( dp )  :: kx , ky , kz , k2 , kx2 , ky2 , kz2 , twopioLx , twopioLz , twopioLy , Rc
real(dp):: Nk , Rc ! Total number of k points = nfft1*nfft2*nfft3
!nf1=nfft1/2
!nf2=nfft2/2
!nf3=nfft3/2
Nk=real (nfft1*nfft2*nfft3 , dp)
!twopioLx=twopi/Lx
!twopioLy=twopi/Ly
!twopioLz=twopi/Lz
Rc=0.0_dp
! init outputs
neq = 0.0_dp
! read cg_vect and get density and polarization from it
! note again that rho is the density per angle so that 
icg = 0
do species = 1 , nb_species
  do i = 1 , nfft1
    do j = 1 , nfft2
      do k = 1 , nfft3
        local_density = 0.0_dp
        do omega = 1 , angGrid%n_angles
         
          do p=1, nb_psi
            icg = icg + 1
            rho_over_fourpi = cg_vect ( icg ) ** 2 / (twopi**2*2.0_dp/sym_order)
            local_density = local_density + weight ( omega ) *weight_psi(p)* rho_over_fourpi ! integral of rho over all orientations ie 'n'
          end do
        end do
        neq ( i , j , k , species ) = local_density
      end do
    end do
  end do
end do
  !convolute with a gaussian
allocate ( rho_k ( nfft1 / 2 + 1 , nfft2 , nfft3 , nb_species ) )
do species = 1 , nb_species
  in_forward = neq ( : , : , : , species )
  call dfftw_execute ( plan_forward )
  rho_k ( : , : , : , species ) = out_forward 
end do
do species = 1 , nb_species
  rho_k ( : , : , : , nb_species ) = rho_k ( : , : , : , nb_species ) * exp ( - k2 ( : , : , : ) * Rc **2 / 2.0_dp )
end do
do species = 1 , nb_species
  in_backward = rho_k ( : , : , : , species )
  call dfftw_execute ( plan_backward )
  neq ( : , : , : , species ) = out_backward/Nk 
end do 
open(11,file='output/density_along_x.dat')
do i=1,nfft1
write(11,*) (i-1)*deltax, neq(i,nfft2/2+1,nfft3/2+1,1)
end do
close(11)
open(12,file='output/density_along_y.dat')
do i=1,nfft2
write(12,*) (i-1)*deltay, neq(nfft1/2+1,i,nfft3/2+1,1)
end do
close(12)
open(13,file='output/density_along_z.dat')
do i=1,nfft3
write(13,*) (i-1)*deltaz, neq(nfft1/2+1,nfft2/2+1,i,1)
end do
close(13)
print*, 'blablabla' , Rc
end subroutine get_final_density

! this subroutine init the density as a function of position and orientation
! it should be initiated to exp(-beta*Vext_total) but
! Vext_q is the electrostatic part of Vext_total, and is pathologic (it sometimes diverges)
! we thus init the density not using vext, but Vext_total - Vext_q
subroutine init_density
use precision_kinds , only : dp , i2b
use system , only : nfft1 , nfft2 , nfft3 , beta , nb_species, nb_psi
use quadrature, only: angGrid
use cg , only : cg_vect
use external_potential , only : Vext_total , Vext_q
use input , only : input_line, input_log
implicit none
real(dp):: local_density0 !> @var local_density0 is the density at a space and angular grid point
integer(i2b):: i , j , k , o , p , icg ! dummy
real(dp):: Vext_total_local , Vext_total_local_min_Vext_q ! dummy
integer(i2b):: species ! dummy between 1 and nb_species
! If reuse_density then just read density in binary format
if (input_log('reuse_density')) then
  open( 10 , file = 'input/density.bin' , form = 'unformatted' , iostat = k )
  if ( k /= 0 ) then
    print *, 'density.bin not found. stop. bug at init_density.f90'
    stop
  end if
  read ( 10 ) cg_vect
  print *, 'cg_vect has been read'
  return ! get out of subroutine
end if
  
! We minimize with respect to sqrt(density) and not directly with respect to the density.
! This allows the variable (sqrt(density)) to be strictly positive.
! icg is the index of the cg_vect 1dimensional array. it obliges us to loop over x then y then z then omega
icg = 0
do species = 1 , nb_species
  do i = 1 , nfft1
    do j = 1 , nfft2
      do k = 1 , nfft3
        do o = 1 , angGrid%n_angles
          do p= 1, nb_psi
          icg = icg + 1
          Vext_total_local = Vext_total ( i , j , k , o ,p , species )
          if ( allocated ( Vext_q ) ) then
            Vext_total_local_min_Vext_q = Vext_total_local - Vext_q ( i , j , k , o , p , species )
          else
            Vext_total_local_min_Vext_q = Vext_total_local
          end if
          if ( Vext_total_local >= 100.0_dp .or. Vext_total_local_min_Vext_q >= 100.0_dp ) then
            local_density0 = tiny ( 0.0_dp ) ! Don't put 0 as it induces problems in the calculations of log(density) in ideal part of F.
          else
            local_density0 = exp ( - beta * Vext_total_local_min_Vext_q )
          end if
          ! put result in cg_vect. note that we do prefer to minimize sqrt(density) than the density in order to avoid sign problems
          cg_vect ( icg ) = sqrt ( local_density0 )
          end do !psi
        end do ! omega
      end do ! nfft3
    end do ! nfft2
  end do ! nfft1
end do ! species
! only Vext_total is used in the functional. We may deallocate it now.
!if ( allocated ( Vext_q ) ) deallocate ( Vext_q )
end subroutine init_density

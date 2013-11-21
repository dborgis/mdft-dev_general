subroutine get_charge_density
use precision_kinds, only : i2b, dp
use system, only : nfft1, nfft2, nfft3, deltaV , nb_species, rho_0_multispec
use quadrature, only : molRotGrid, angGrid, molRotGrid
use external_potential, only : x_charge ,y_charge ,z_charge ,q_charge ,nb_of_interpolation,Fcoul, v_c
USE cg, ONLY: CG_vect
use constants, only : qfact
implicit none
real (dp ) :: rho, psi, time1,time0
real ( dp) , allocatable , dimension ( : , : ,  :) :: charge_density
integer ( i2b) :: i , j , k , o , p , species , n,convert_coordinate_into_icg, xtemp, ytemp, ztemp
  allocate ( charge_density ( nfft1 , nfft2 , nfft3 ) )
  charge_density=0.0_dp
Fcoul=0.0_dp
call cpu_time(time0)
print*, rho_0_multispec (1)
do species = 1 , nb_species
  do i = 1 , nfft1
    do j = 1 , nfft2
      do k = 1 , nfft3
        do o = 1 , angGrid%n_angles
          do p=1 , molRotGrid%n_angles            
   
            do n=1, nb_of_interpolation
            xtemp= modulo(-(x_charge(n , p , o)-i),nfft1)+1
            ytemp= modulo(-(y_charge(n , p , o)-j),nfft2)+1
            ztemp= modulo(-(z_charge(n , p , o)-k),nfft3)+1
!print*, i,j,k,xtemp, ytemp, ztemp 
!print*, xtemp, ytemp, ztemp , o , p, convert_coordinate_into_icg ( species , xtemp, ytemp, ztemp, o , p )
            psi=cg_vect( convert_coordinate_into_icg ( species , xtemp, ytemp, ztemp, o , p ) )
            rho = psi ** 2
            
            charge_density (i , j , k ) =  charge_density (i , j , k ) + deltaV * angGrid%weight ( o ) * molRotGrid%weight ( p )&
                                              * rho *q_charge (n , p , o )*rho_0_multispec (species)
            end do   !n
          end do    !psi
        end do  !omega
      end do    !i
    end do      !j
  end do        !k
end do        !species
print*, 'charge_density', minval(charge_density), maxval(charge_density)
!call write_to_cube_file ( charge_density(:,:,:), 'output/charge_density.cube')
do i =1, nfft1
  do j=1, nfft2
    do k=1, nfft3
     Fcoul=Fcoul + charge_density (i ,j , k)*V_c(i,j,k)*qfact
 
    end do
  end do
end do
call cpu_time(time1)
end subroutine

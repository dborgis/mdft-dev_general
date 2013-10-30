!> This subroutine compute the external purely repulsive potential as defined by Dzubiella and Hansen in J Chem Phys 121 (2004)
! This potential has the form V(r)=kbT*(r-R0)^-(12)
! Written by Maximilien Levesque, June 2011 @ Ecole Normale Superieure
! For now This potential is named Vcoul which means as soon as we will have to use both electrostatic and (NOT TRUE ANYMORE : VREP)
! purely repulsive potential, one will HAVE TO re-write this subroutine with unique variable names
subroutine compute_purely_repulsive_potential ( Rotxx , Rotxy , Rotxz , Rotyx , Rotyy , Rotyz , Rotzx , Rotzy , Rotzz )
use precision_kinds , only : dp , i2b
use input , only : input_line,input_dp 
use system, only: nfft1 , nfft2 , nfft3 , nb_psi , x_solv , y_solv , z_solv , x_mol , y_mol , z_mol , &
                        beta , nb_solute_sites , nb_solvent_sites , Lx , Ly , Lz
use constants , only : fourpi
use external_potential , only : Vext_total
use quadrature, only: angGrid
implicit none
real(dp), dimension ( angGrid%n_angles , nb_psi ) , intent (in) :: Rotxx , Rotxy , Rotxz , Rotyx , Rotyy , Rotyz , Rotzx , &
                                                                      Rotzy , Rotzz ! rotationnal matrix for sites
integer(i2b):: i , j , k , o , p , m , n ! dummy
real(dp), dimension ( nb_solvent_sites , nb_psi , angGrid%n_angles ) :: xmod , ymod , zmod
real(dp):: x_grid , y_grid , z_grid ! coordinates of grid nodes
real(dp):: x_m , y_m , z_m ! solvent sites coordinates
real(dp):: x_nm , y_nm , z_nm ! coordinate of vecteur solute-solvent
real(dp):: r_nm2 ! norm**2 of vector x_nm;y_nm;z_nm
real(dp):: V_psi ! dummy
real(dp):: time0 , time1 ! timer start and end
real ( dp ), dimension (:,:,:) , allocatable :: temparray ! dummy
character (50) :: filename
real(dp):: rc, rc2 ! radius of the hard wall in the purely repulsive potential of dzubiella
real(dp):: tempVrep ! temporary Vrep(i,j,k,o)
real(dp):: kbT ! kbT of the purely repulsive potential V such as V(R0+1)=kBT => kbT = kBT
real(dp), dimension ( nfft1 , nfft2 , nfft3 , angGrid%n_angles , nb_psi) :: Vrep
real(dp):: deltax , deltay , deltaz ! == Lx / nfft1 , Ly / nfft2 , Lz / nfft3
! init timer
call cpu_time(time0)
! get the radius of the purely repulsive solute
! the radius is defined such as in Dzubiella and Hansen, J Chem Phys 121 , 2011
! look for tag 'purely_repulsive_solute_radius' in dft.in for hard wall thickness
Rc=input_dp('radius_of_purely_repulsive_solute')
write (*,*) 'the radius of the purely repulsive solute is ' , Rc
print*, size(Vrep)
! init variables
Vrep = 0.0_dp
Rc2 = Rc**2
kbT = 1.0_dp/beta
! tabulate coordinates of solvent sites for each omega and psi angles
do o = 1 , angGrid%n_angles
  do p = 1 , nb_psi
    do m = 1 , nb_solvent_sites
      xmod (m,p,o) = Rotxx (o,p) * x_solv (m) + Rotxy (o,p) * y_solv (m) + Rotxz (o,p) * z_solv (m)
      ymod (m,p,o) = Rotyx (o,p) * x_solv (m) + Rotyy (o,p) * y_solv (m) + Rotyz (o,p) * z_solv (m)
      zmod (m,p,o) = Rotzx (o,p) * x_solv (m) + Rotzy (o,p) * y_solv (m) + Rotzz (o,p) * z_solv (m)
    end do
  end do
end do
! distance between two grid points in x y z directions
deltax = Lx / real ( nfft1 , dp )
deltay = Ly / real ( nfft2 , dp )
deltaz = Lz / real ( nfft3 , dp )
! TODO CARE HERE ONLY ONE SPECIES
do k = 1 , nfft3
  ! print screen how many grid points in z have still to be calculated
  write (*,*) nfft3 - k
  z_grid = real(k-1,dp) * deltaz
  do j = 1 , nfft2
    y_grid = real(j-1,dp) * deltay
    do i = 1 , nfft1
      x_grid = real(i-1,dp) * deltax
      do o = 1 , angGrid%n_angles
        tempVrep = 0.0_dp
        do p = 1 , nb_psi
          V_psi = 0.0_dp
          do m = 1 , nb_solvent_sites
            x_m = x_grid + xmod (m,p,o)
            y_m = y_grid + ymod (m,p,o)
            z_m = z_grid + zmod (m,p,o)
            do n = 1,nb_solute_sites
              x_nm = x_m - x_mol (n)
              y_nm = y_m - y_mol (n)
              z_nm = z_m - z_mol (n)
              r_nm2 = x_nm**2+y_nm**2+z_nm**2
              ! if we're in the hard repulsive zone v is huge, else its purely repulsive
              if ( r_nm2 <= Rc2 ) then
                V_psi = huge(1.0_dp)
              else
                V_psi = V_psi + kbT * ( sqrt ( r_nm2 ) - Rc )**(-12)
              end if
            end do ! solute
          end do ! solvent
!          tempVrep = tempVrep + exp ( - beta * V_psi )
      
        ! care about not having log(0)        
      !  if ( tempVrep < 1.e-10_dp ) then
      !    Vrep (i,j,k,o) = 25.0_dp * kbT
      !  else
          Vrep (i,j,k,o, p) = V_psi
        if ( Vrep (i,j,k,o, p) > 100.0_dp ) Vrep (i,j,k,o, p) = 100.0_dp
      !  end if
        end do !ploop ! psi
        ! limit maximum value of Vrep to 100.0
      end do
    end do
  end do
end do
print*, 'BETA' , beta
! vrep is only part of vext so put vrep in vext
Vext_total (:,:,:,:, :,1) = Vext_total (:,:,:, :,:,1) + vrep (:,:,:,:, :)
! warn user about vrep extrema for debugging
write (*,*) 'minval(Vrep) = ' , minval ( Vrep )
write (*,*) 'maxval(Vrep) = ' , maxval ( Vrep )
! stop timer
call cpu_time(time1)
write (*,*) 'time for compute_purely_repulsive_potential ' , time1 - time0
! mean over orientations and print
allocate ( temparray ( nfft1 , nfft2 , nfft3 ) )
call mean_over_orientations ( Vrep , temparray )
temparray = temparray / fourpi
filename = 'output/Vrep.cube'
call write_to_cube_file ( temparray , filename )
filename = 'output/Vrep_radial.dat'
call compute_rdf ( temparray , filename )
deallocate ( temparray )
end subroutine compute_purely_repulsive_potential

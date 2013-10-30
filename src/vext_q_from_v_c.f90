subroutine vext_q_from_v_c (Rotxx,Rotxy,Rotxz,Rotyx,Rotyy,Rotyz,Rotzx,Rotzy,Rotzz)
use precision_kinds , only : dp , i2b
use system , only : chg_mol , chg_solv , nb_psi, x_solv , y_solv , z_solv , nb_solvent_sites , nb_species , Lx , Ly , &
                    Lz , deltax , deltay , deltaz , beta , id_solv , beta , nfft1 , nfft2 , nfft3
! chg_mol (nb_id_mol ) = charge of each solute site
! chg_solv ( nb_id_solv ) = charge of each solvent site
! nb_psi = total number of psi angles
use quadrature, only: angGrid
use external_potential , only : v_c , vext_q , vext_lj
! v_c = electrostatic potential from charge density and poisson equation
! vext_q = electrostatic potential energy in general and as used in the calculation of the total external potential
use constants , only : fourpi , qfact , qunit
! fourpi = 4pi
! qfact = qunit ** 2 * 1.0e-3_dp * Navo / ( fourpi * eps0 * 1.0e-10_dp ) ! electrostatic potential unit so that QFACT*q*q/r is kJ/mol
implicit none
! Declaration zone
real(dp), dimension ( angGrid%n_angles , nb_psi ) , intent(in) :: Rotxx , Rotxy , Rotxz , Rotyx , Rotyy , Rotyz , Rotzx ,&
                                                                        Rotzy , Rotzz
integer(i2b):: i , j , k , o , p , m,z ! dummy
real(dp), dimension ( nb_solvent_sites , nb_psi , angGrid%n_angles ) :: xmod , ymod , zmod
integer(i2b):: species ! dummy
real(dp):: xq , yq , zq ! solvent coordinates in indices referential (ie real between 0 and nfft1+1)
integer(i2b):: im , jm , km , ip , jp , kp ! 6 indices from which all corners of cube surrounding point charge are defined
real(dp):: wim , wjm , wkm , wip , wjp , wkp ! weights associated to each 'corner' index
real(dp):: vpsi ! external potential for a given i,j,k,omega & psi.
real(dp):: average_over_psi ! boltzmann average of vpsi over psi for a given i,j,k,omega
real(dp):: time0 , time1 ! timers
write(*,*) '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
write(*,*) 'qunit' , qunit , beta ,qfact
write(*,*) '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
! time 0
call cpu_time ( time0 )
z=0
! be sure this subroutine comes at right moment
! init (again during debug) vext_q
if ( .not. allocated ( vext_q ) ) allocate ( vext_q ( nfft1 , nfft2 , nfft3 , angGrid%n_angles , nb_psi , nb_species ) )
vext_q = 0.0_dp
print *, 'vext_q_from_v_c : min = ' , minval ( Vext_q ) , ' ; max = ' , maxval ( Vext_q ) 
! test if all solutes or solvent sites have zero charge then don't waste your time : go to end of subroutine
if ( minval ( chg_mol ) == 0.0_dp .and. maxval ( chg_mol ) == 0.0_dp ) then   ! solute is not charged
 print *, 'solute has no charge'
  return
else if ( minval ( chg_solv ) == 0.0_dp .and. maxval ( chg_solv ) == 0.0_dp ) then   ! solvent is not charged
  print *, 'solvent has no charge'
  return
end if
open(11, file='output/v_c.dat')
 do i=1, nfft1
   do j=1, nfft2
     do k=1, nfft3
write(11,*) v_c(i , j , k)
 end do
   end do
     end do
     
close(11)
! Tabulate rotation matrix * solvent coordinates
do o = 1 , angGrid%n_angles
  do p = 1 , nb_psi
    do m = 1 , nb_solvent_sites
      xmod ( m , p , o ) = Rotxx ( o , p ) * x_solv ( m ) + Rotxy ( o , p ) * y_solv ( m ) + Rotxz ( o , p ) * z_solv ( m )
      ymod ( m , p , o ) = Rotyx ( o , p ) * x_solv ( m ) + Rotyy ( o , p ) * y_solv ( m ) + Rotyz ( o , p ) * z_solv ( m )
      zmod ( m , p , o ) = Rotzx ( o , p ) * x_solv ( m ) + Rotzy ( o , p ) * y_solv ( m ) + Rotzz ( o , p ) * z_solv ( m )
!print*, species, m , o , p ,  xmod ( m , p , o ) ,  ymod ( m , p , o )  , zmod ( m , p , o )
    end do
  end do
end do

! Compute external potential due to charge density
do species = 1 , nb_species
  do i = 1 , nfft1
    do j = 1 , nfft2
      do k = 1 , nfft3
        do o = 1 , angGrid%n_angles
          average_over_psi = 0.0_dp ! init average over psi which is for all psi for a given i,j,k,omega
          do p = 1 , nb_psi
            vpsi = 0.0_dp ! vpsi is initiated. We make a boltzmann average over psi angles
            do m = 1 , nb_solvent_sites
              if ( chg_solv ( id_solv (m) ) == 0.0_dp ) cycle
              xq = modulo ( real( i-1 ,dp) * deltax + xmod (m,p,o) , Lx ) / deltax
              yq = modulo ( real( j-1 ,dp) * deltay + ymod (m,p,o) , Ly ) / deltay
              zq = modulo ( real( k-1 ,dp) * deltaz + zmod (m,p,o) , Lz ) / deltaz
              im = int ( xq ) + 1
              jm = int ( yq ) + 1
              km = int ( zq ) + 1
              if ( im == nfft1 + 1 ) im = 1
              if ( jm == nfft2 + 1 ) jm = 1
              if ( km == nfft3 + 1 ) km = 1
              ip = im + 1
              jp = jm + 1
              kp = km + 1
              if ( ip == nfft1 + 1 ) ip = 1
              if ( jp == nfft2 + 1 ) jp = 1
              if ( kp == nfft3 + 1 ) kp = 1
              wim = (    1.0_dp - (   xq - real( int( xq ,i2b),dp)   )    )
              wjm = (    1.0_dp - (   yq - real( int( yq ,i2b),dp)   )    )
              wkm = (    1.0_dp - (   zq - real( int( zq ,i2b),dp)   )    )
              wip = (             (   xq - real( int( xq ,i2b),dp)   )    )
              wjp = (             (   yq - real( int( yq ,i2b),dp)   )    )
              wkp = (             (   zq - real( int( zq ,i2b),dp)   )    )
              vpsi = vpsi + (chg_solv ( id_solv (m) ) * qfact * ( V_c (im,jm,km) * wim * wjm * wkm &
                                                               + V_c (ip,jm,km) * wip * wjm * wkm &
                                                               + V_c (im,jp,km) * wim * wjp * wkm &
                                                               + V_c (im,jm,kp) * wim * wjm * wkp &
                                                               + V_c (ip,jp,km) * wip * wjp * wkm &
                                                               + V_c (ip,jm,kp) * wip * wjm * wkp &
                                                               + V_c (im,jp,kp) * wim * wjp * wkp &
                                                               + V_c (ip,jp,kp) * wip * wjp * wkp ))
!print*, species, p , o ,k ,j ,i, vpsi
       !     print*, vpsi
            end do ! solvent sites
!Potential depends of psi now
!            average_over_psi = average_over_psi + exp ( - beta * vpsi ) ! this is a boltzmann average over psi
       !     write(*,*)  average_over_psi
 
          ! Things are clear but surely useless here.
!          if ( average_over_psi < 1.e-10_dp ) then
!            Vext_q ( i , j , k , o , species ) = 25.0_dp / beta ! TODO magic numbers are not welcome
!          else
!            Vext_q ( i , j , k , o , species ) = - log ( average_over_psi / nb_psi ) / beta
!          end if
!          if ( average_over_psi > 0.0_dp ) then
            Vext_q ( i , j , k , o , p , species ) = vpsi !- log ( average_over_psi / nb_psi ) / (beta)
!          else
!            Vext_q ( i , j , k , o , species ) = huge(1.0_dp)
!          end if
          if ( Vext_q ( i , j , k , o , p , species )        > 100.0_dp ) then
             Vext_q ( i , j , k , o , p  , species )         = 100.0_dp
          end if
          
          if ( Vext_q ( i , j , k , o , p , species )        < -100.0_dp ) then
            Vext_q ( i , j , k , o , p , species )        = -100.0_dp
          end if
       
          end do ! psi 
  
        end do ! omega
      end do ! k
    end do ! j
  end do ! i
end do ! species
! warn user about min and max values and time
call cpu_time(time1)
print *, 'vext_q_from_v_c : min = ' , minval ( Vext_q ) , ' ; max = ' , maxval ( Vext_q ) , ' ; in (sec) ' , time1 - time0
! Get the external potential over orientations and print it
!allocate ( temparray ( nfft1 , nfft2 , nfft3 ) )
!call mean_over_orientations ( Vext_q (:,:,:,:,1) , temparray )
!temparray = temparray / fourpi
!!filename = 'output/vext_q_from_v_c.cube'
!!call write_to_cube_file ( temparray , filename )
!filename = 'output/vext_q_from_v_c_along-z.dat'
!call compute_z_density ( temparray , filename )
!deallocate ( temparray )
! DEBUG ONLY
open ( unit = 10 , file = 'output/vext_q_from_v_c_along-z_no_planar_average.dat' )
do k = 1 , nfft3
  write ( 10 , * ) (real(k-1,dp)*deltaz-Lz/2.0_dp) , (Vext_q(nfft1/2,nfft2/2,k,1,1,1) / QFACT)
end do
close(10)
end subroutine vext_q_from_v_c

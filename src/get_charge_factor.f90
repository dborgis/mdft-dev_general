!TO DO use a smaller array to store sigma by finding the larger indexes (on x, y, z) on which some charges are projected
subroutine get_charge_factor ( Rotxx,Rotxy,Rotxz,Rotyx,Rotyy,Rotyz,Rotzx,Rotzy,Rotzz ) 
use precision_kinds , only : i2b , dp
use system , only : nfft1 , nfft2 , nfft3 , x_solv  ,y_solv , z_solv , chg_solv , id_solv , &
 nb_solvent_sites, &
 Lx , Ly , Lz , deltax , deltay , deltaz, wigma
use external_potential, only : x_charge, y_charge, z_charge, q_charge, nb_of_interpolation
use quadrature, only: angGrid, molRotGrid
implicit none
integer(i2b):: i , j , o , p , m
real(dp), dimension(angGrid%n_angles,molRotGrid%n_angles), intent(in) :: Rotxx,Rotxy,Rotxz,Rotyx,Rotyy,Rotyz,Rotzx,Rotzy,Rotzz
real (dp ) :: xq , yq , zq , wim , wjm , wip , wjp , wkm , wkp, Dx, Dy, Dz
integer (kind = i2b ) :: im , jm , km , ip , jp , kp 
real (kind= dp ) , allocatable , dimension ( : , : , : ) :: xmod , ymod , zmod
 character (len=1) , allocatable , dimension ( : ) :: coordinate
integer ( i2b) , allocatable , dimension ( : ) :: nombre_image
!allocate ( wigma ( nfft1, nfft2, nfft3 , angGrid%n_angles , molRotGrid%n_angles ) ) 
allocate ( nombre_image ( nb_solvent_sites ) )
allocate ( coordinate ( nb_solvent_sites ) )
!check how many coordinates it is neccessary to store if there is one site of the solvent in ( 0 , 0 , 0 ) it will not move by rotation
do m = 1 , nb_solvent_sites
   xq = modulo ( x_solv ( m )  , Lx ) / deltax
   yq = modulo ( y_solv ( m ) , Ly ) / deltay
   zq = modulo ( z_solv ( m ) , Lz ) / deltaz
   Dx = xq - real(int(xq,i2b),dp) 
   Dy = yq - real(int(yq,i2b),dp) 
   Dz = zq - real(int(zq,i2b),dp) 
   if ( Dx==0.0_dp .and. Dy==0.0_dp .and. Dz==0.0_dp ) then
     nombre_image ( m ) = 1
!   else if ( Dx==0.0_dp .and. Dy==0.0_dp) then
!     nombre_image( m ) =  2
!     coordinate (m ) = 'z'
!   else if ( Dx==0.0_dp .and. Dz==0.0_dp) then
!     nombre_image( m ) =  2
!     coordinate (m ) = 'y'
!   else if ( Dy==0.0_dp .and. Dz==0.0_dp) then
!     nombre_image( m ) =  2
!     coordinate (m ) = 'x'
!   else if ( Dx==0.0_dp ) then
!     nombre_image( m ) =  4
!     coordinate (m ) = 'x'  
!   else if ( Dy==0.0_dp ) then
!     nombre_image( m ) =  4
!     coordinate (m ) = 'y'  
!   else if ( Dz==0.0_dp ) then
!     nombre_image( m ) =  4
!     coordinate (m ) = 'z'     
   else 
     nombre_image( m ) =  8
   end if
end do
print*, nombre_image
!evaluate the total number of images which is necessary to store
nb_of_interpolation = sum ( nombre_image )
print*, nb_of_interpolation
allocate (x_charge ( nb_of_interpolation , molRotGrid%n_angles , angGrid%n_angles ) )
allocate (y_charge (nb_of_interpolation , molRotGrid%n_angles , angGrid%n_angles ) )
allocate (z_charge (nb_of_interpolation , molRotGrid%n_angles , angGrid%n_angles ) )
allocate (q_charge (nb_of_interpolation , molRotGrid%n_angles , angGrid%n_angles ) )
x_charge=0
y_charge=0
z_charge=0
q_charge=0.0_dp
allocate ( xmod ( nb_solvent_sites  , molRotGrid%n_angles , angGrid%n_angles) )
allocate ( ymod ( nb_solvent_sites  , molRotGrid%n_angles , angGrid%n_angles) )
allocate ( zmod ( nb_solvent_sites  , molRotGrid%n_angles , angGrid%n_angles) )
xmod=0.0_dp
ymod=0.0_dp
zmod=0.0_dp
do m=1, nb_solvent_sites                                                                         !loop on the different sites of the solvent
  
   if (nombre_image(m) == 8 ) then                                                              !check how many images it is necessary to calculate
print*, 'boucle8'      
       do o =1 , angGrid%n_angles
       
         do p=1 , molRotGrid%n_angles
              i=0
              do j = 1 , m-1
              i= i+nombre_image(j)                                                              !Find the value of index i for this solute site
              end do
              xmod(m,p,o)= Rotxx(o,p)*x_solv(m) + Rotxy(o,p)*y_solv(m) + Rotxz(o,p)*z_solv(m)   !new  coordinates of this site after rotation
              ymod(m,p,o)= Rotyx(o,p)*x_solv(m) + Rotyy(o,p)*y_solv(m) + Rotyz(o,p)*z_solv(m)
              zmod(m,p,o)= Rotzx(o,p)*x_solv(m) + Rotzy(o,p)*y_solv(m) + Rotzz(o,p)*z_solv(m)  
  
              if ( chg_solv ( id_solv ( m ) ) == 0.0_dp ) cycle                                 !if charge of this solvent site is 0 it is not necessary to continue
                                                                                                ! transform cartesian coordinates 0 <= x < Lx in 'indicial coordinates' 0 <= xq < nfft1
              xq = modulo ( xmod(m,p,o) , Lx ) / deltax                                         ! define weights associated with each corner
              yq = modulo ( ymod(m,p,o) , Ly ) / deltay
              zq = modulo ( zmod(m,p,o) , Lz ) / deltaz
       
              if (xq==80.0_dp) xq=0.0_dp
               if (yq==80.0_dp) yq=0.0_dp
                if (zq==80.0_dp) zq=0.0_dp
              im = int ( xq ) + 1                                                                ! get coordinates of grid node juste below (corner of the cube with smallest indices) 
              jm = int ( yq ) + 1                                                                ! +1 is because indexation does not begin to 0 but to 1. Thus, when coordinate xq=0, it corresponds to index 1.
              km = int ( zq ) + 1
              wim = (    1.0_dp - (   xq - real(int(xq,i2b),dp)   )    )
              wjm = (    1.0_dp - (   yq - real(int(yq,i2b),dp)   )    )
              wkm = (    1.0_dp - (   zq - real(int(zq,i2b),dp)   )    )
 
              ip = im + 1
              jp = jm + 1                                                                         ! grid node juste above (corner of the cube with highest indices)
              kp = km + 1
              wip = (             (   xq - real(int(xq,i2b),dp)   )    )
              wjp = (             (   yq - real(int(yq,i2b),dp)   )    )
              wkp = (             (   zq - real(int(zq,i2b),dp)   )    )
                                                                                                  ! this corner should never have coordinates higher than nfft
              if ( ip == nfft1 + 1 ) ip = 1
              if ( jp == nfft2 + 1 ) jp = 1
              if ( kp == nfft3 + 1 ) kp = 1
  
             i=i+1                                                                                !Store charges and coordiantes of those images
             
             x_charge ( i,  p, o ) =   im
             y_charge ( i,  p, o ) =   jm
             z_charge ( i,  p, o ) =   km
             q_charge ( i,  p, o ) = wim*wjm*wkm* chg_solv ( id_solv ( m ) ) 
             i=i+1
             x_charge ( i,  p, o ) =   ip
             y_charge ( i,  p, o ) =   jm
             z_charge ( i,  p, o ) =   km
             q_charge ( i,  p, o ) = wip*wjm*wkm* chg_solv ( id_solv ( m ) ) 
             i=i+1
             x_charge ( i,  p, o ) =   im
             y_charge ( i,  p, o ) =   jp
             z_charge ( i,  p, o ) =   km
             q_charge ( i,  p, o ) = wim*wjp*wkm* chg_solv ( id_solv ( m ) ) 
             i=i+1
             x_charge ( i,  p, o ) =   im
             y_charge ( i,  p, o ) =   jm
             z_charge ( i,  p, o ) =   kp
             q_charge ( i,  p, o ) = wim*wjm*wkp* chg_solv ( id_solv ( m ) ) 
             i=i+1
             x_charge ( i,  p, o ) =   ip
             y_charge ( i,  p, o ) =   jp
             z_charge ( i,  p, o ) =   km
             q_charge ( i,  p, o ) = wip*wjp*wkm* chg_solv ( id_solv ( m ) ) 
             i=i+1
             x_charge ( i,  p, o ) =   im
             y_charge ( i,  p, o ) =   jp
             z_charge ( i,  p, o ) =   kp
             q_charge ( i,  p, o ) = wim*wjp*wkp* chg_solv ( id_solv ( m ) ) 
             i=i+1
             x_charge ( i,  p, o ) =   ip
             y_charge ( i,  p, o ) =   jm
             z_charge ( i,  p, o ) =   kp
             q_charge ( i,  p, o ) = wip*wjm*wkp* chg_solv ( id_solv ( m ) ) 
             i=i+1
             x_charge ( i,  p, o ) =   ip
             y_charge ( i,  p, o ) =   jp
             z_charge ( i,  p, o ) =   kp
             q_charge ( i,  p, o ) = wip*wjp*wkp* chg_solv ( id_solv ( m ) ) 
         end do
       end do
   else if (nombre_image(m) == 1 ) then
     !print*, 'boucle1'         
       do o =1 , angGrid%n_angles
       
          do p=1 , molRotGrid%n_angles
             i=0
             do j = 1 , m-1
               i= i+nombre_image(j)
             end do
             xmod(m,p,o)= Rotxx(o,p)*x_solv(m) + Rotxy(o,p)*y_solv(m) + Rotxz(o,p)*z_solv(m)
             ymod(m,p,o)= Rotyx(o,p)*x_solv(m) + Rotyy(o,p)*y_solv(m) + Rotyz(o,p)*z_solv(m)   
             zmod(m,p,o)= Rotzx(o,p)*x_solv(m) + Rotzy(o,p)*y_solv(m) + Rotzz(o,p)*z_solv(m)  
  
             if ( chg_solv ( id_solv ( m ) ) == 0.0_dp ) cycle
             xq = modulo ( xmod(m,p,o) , Lx ) / deltax
             yq = modulo ( ymod(m,p,o) , Ly ) / deltay
             zq = modulo ( zmod(m,p,o) , Lz ) / deltaz
             im = int ( xq ) + 1
             jm = int ( yq ) + 1
             km = int ( zq ) + 1
             wim = (    1.0_dp - (   xq - real(int(xq,i2b),dp)   )    )
             wjm = (    1.0_dp - (   yq - real(int(yq,i2b),dp)   )    )
             wkm = (    1.0_dp - (   zq - real(int(zq,i2b),dp)   )    )
             ip = im + 1
             jp = jm + 1
             kp = km + 1
             wip = (             (   xq - real(int(xq,i2b),dp)   )    )
             wjp = (             (   yq - real(int(yq,i2b),dp)   )    )
             wkp = (             (   zq - real(int(zq,i2b),dp)   )    )
            if ( ip == nfft1 + 1 ) ip = 1
            if ( jp == nfft2 + 1 ) jp = 1
            if ( kp == nfft3 + 1 ) kp = 1
            i=i+1
             
            x_charge ( i,  p, o ) =   im
            y_charge ( i,  p, o ) =   jm
            z_charge ( i,  p, o ) =   km
            q_charge ( i,  p, o ) = wim*wjm*wkm* chg_solv ( id_solv ( m ) ) 
          end do  !end loop on psi
       end do     !end loop on omega
   end if         !end condition on number of images
end do    !end loop on solvent site                                                     
print*,  'HERE IS I :::: ' ,i
print*, 'Somme charge', Sum(q_charge)
open (11, file='output/x_charge')
write(11,*) , x_charge
 close(11)
open (11, file='output/z_charge')
write(11,*) , z_charge
 close(11)
open (12, file='output/q_charge')
write(12,*) , q_charge
 close(12)
!check if both ways to evaluate charge factor are similar.
end subroutine

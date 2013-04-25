subroutine compute_angular_grid_leb( Rotxx, Rotxy, Rotxz, Rotyx, Rotyy, Rotyz, Rotzx, Rotzy, Rotzz)

use precision_kinds , only: i2b,dp

use quadrature , only :  Omx , Omy , Omz , weight , x_psi , weight_psi , sym_order,x_leb , y_leb, z_leb , weight_leb 

use system , only :  nb_psi , nb_omega

use constants , only : pi , twopi , fourpi





implicit none

integer(i2b) ::  n !dummy

integer(i2b) :: n_psi

real(dp) :: psii

real(dp) :: phi , theta

real(dp) :: cos_theta

real(dp) :: sin_theta

real(dp) :: cos_phi

real(dp) :: sin_phi

real(dp) :: cos_psi

real(dp) :: sin_psi

real(dp), dimension ( nb_omega, nb_psi ), intent(out) :: Rotxx

real(dp), dimension ( nb_omega, nb_psi ), intent(out) :: Rotxy

real(dp), dimension ( nb_omega, nb_psi ), intent(out) :: Rotxz

real(dp), dimension ( nb_omega, nb_psi ), intent(out) :: Rotyx

real(dp), dimension ( nb_omega, nb_psi ), intent(out) :: Rotyy

real(dp), dimension ( nb_omega, nb_psi ), intent(out) :: Rotyz

real(dp), dimension ( nb_omega, nb_psi ), intent(out) :: Rotzx

real(dp), dimension ( nb_omega, nb_psi ), intent(out) :: Rotzy

real(dp), dimension ( nb_omega, nb_psi ), intent(out) :: Rotzz


allocate (weight (nb_omega ) )

if ( .not. allocated ( Omx ) ) allocate ( Omx ( nb_omega ) ) ! orientational vector along x
if ( .not. allocated ( Omy ) ) allocate ( Omy ( nb_omega ) ) ! orientational vector along x
if ( .not. allocated ( Omz ) ) allocate ( Omz ( nb_omega ) ) ! orientational vector along x


do n=1 , nb_omega

  weight(n) = weight_leb(n)

  theta=acos(z_leb(n))
 
    if    ( x_leb(n)  == 0.0_dp .and. y_leb(n)  == 0.0_dp ) then 

          phi=0.0_dp

    else if  (y_leb(n)>=0.0_dp) then 

          phi=acos(x_leb(n)/(sqrt ( x_leb(n)**2 + y_leb(n)**2 ) ) )
   
    else  

          phi=twopi - acos(x_leb(n)/(sqrt ( x_leb(n)**2 + y_leb(n)**2 ) ) )
  
    end if

     cos_theta=cos(theta)

     sin_theta=sin(theta) 

     cos_phi=cos(phi)

     sin_phi=sin(phi)    

  !   print*, cos_theta , sin_theta ,  cos_phi ,  sin_phi , n

     OMx ( n ) = sin_theta * cos_phi

     OMy ( n ) = sin_theta * sin_phi

     OMz ( n ) = cos_theta

       do n_psi = 1 , nb_psi

         psii = x_psi(n_psi)
       
         cos_psi = cos(psii) !1

         sin_psi = sin(psii) !0

        
    
        Rotxx(n,n_psi) =  cos_theta*cos_phi*cos_psi-sin_phi*sin_psi
       !Rotxx(1,1)           =  0
       !Rotxx(2,1)           =  0
        Rotxy(n,n_psi) = -cos_theta*cos_phi*sin_psi-sin_phi*cos_psi
       !Rotxy(1,1)           =  0
       !Rotxy(2,1)           =  0
        Rotxz(n,n_psi) =  sin_theta*cos_phi
       !Rotxz(1,1)           =  1
       !Rotxz(2,1)           =  -1
        Rotyx(n,n_psi) =  cos_theta*sin_phi*cos_psi+cos_phi*sin_psi
        Rotyy(n,n_psi) = -cos_theta*sin_phi*sin_psi+cos_phi*cos_psi
        Rotyz(n,n_psi) =  sin_theta*sin_phi
        Rotzx(n,n_psi) = -sin_theta*cos_psi
        Rotzy(n,n_psi) =  sin_theta*sin_psi
        Rotzz(n,n_psi) =  cos_theta

   end do



end do



! check if sum over all omega of weight(omega) is fourpi

if ( abs ( sum ( weight ( : ) ) - fourpi )  > 1.0e-10_dp ) then

  print *, 'problem detected in compute_angular_grid.f90 :'

  print *, 'sum over omegas of weight(omega) is not 4pi. stop'

  stop

end if



if ( abs ( sum ( weight_psi ( : ) ) - twopi/sym_order )  > 1.0e-10_dp ) then

  print *, 'problem detected in compute_angular_grid.f90 :'

  print *, 'sum over omegas of weight_psi(omega) is not 2pi/sym_order. stop'

  stop

end if




end subroutine

!> Compute angular grid properties : Omx, Omy, Omz, weight
subroutine compute_angular_grid ( Rotxx, Rotxy, Rotxz, Rotyx, Rotyy, Rotyz, Rotzx, Rotzy, Rotzz)
use precision_kinds , only: i2b,dp
use quadrature , only : w_legendre , x_legendre , Omx , Omy , Omz , weight , x_psi , weight_psi , sym_order
use system , only : nb_omega , nb_psi , nb_legendre 
use constants , only : pi , twopi , fourpi
implicit none
integer(i2b) :: n_psi
integer(i2b) :: n_theta
integer(i2b) :: n_phi
integer(i2b) :: n_omega
real(dp) :: psii
real(dp) :: phi 
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
write(*,*)'>> computing angular grid'
!> allocate what we want to compute
if ( .not. allocated ( Omx ) ) allocate ( Omx ( nb_omega ) ) ! orientatioal vector along x
if ( .not. allocated ( Omy ) ) allocate ( Omy ( nb_omega ) ) ! orientational vector along x
if ( .not. allocated ( Omz ) ) allocate ( Omz ( nb_omega ) ) ! orientational vector along x
if ( .not. allocated ( weight ) ) allocate ( weight ( nb_omega ) ) ! weight of each orientation
!> initiate
!> Calcul de Daniel pour backup : fastend up
! do n_theta = 1, nb_legendre
!  cos_theta = x_legendre(n_theta,nb_legendre)
!  sin_theta = sqrt(one - cos_theta**2)
!  OMz(n_omega) = cos_theta
!  weight(n_omega) = w_legendre(n_theta,nb_legendre)*pi/nb_legendre
!  do n_phi = 1, 2*nb_legendre
!   n_omega = (n_theta-1)*2*nb_legendre + n_phi
!   angle_phi = (n_phi-1)*twopi/(2*nb_legendre)
!   cos_phi = cos(angle_phi)
!   sin_phi = sin(angle_phi)
!   OMx(n_omega) = sin_theta*cos_phi
!   OMy(n_omega) = sin_theta*sin_phi
!   totalweight = totalweight + weight(n_omega)
!   do n_psi = 1, nb_psi
!    angle_psi = (n_psi-1)*pi/nb_psi
!    cos_psi = cos(angle_psi)
!    sin_psi = sin(angle_psi)
!    Rotxx(n_omega,n_psi) =  cos_theta*cos_phi*cos_psi-sin_phi*sin_psi
!    Rotxy(n_omega,n_psi) = -cos_theta*cos_phi*sin_psi-sin_phi*cos_psi
!    Rotxz(n_omega,n_psi) =  sin_theta*cos_phi
!    Rotyx(n_omega,n_psi) =  cos_theta*sin_phi*cos_psi+cos_phi*sin_psi
!    Rotyy(n_omega,n_psi) = -cos_theta*sin_phi*sin_psi+cos_phi*cos_psi
!    Rotyz(n_omega,n_psi) =  sin_theta*sin_phi
!    Rotzx(n_omega,n_psi) = -sin_theta*cos_psi
!    Rotzy(n_omega,n_psi) =  sin_theta*sin_psi
!    Rotzz(n_omega,n_psi) =  cos_theta
!   end do ! n_psi
!  end do ! n_phi
! end do ! n_theta
!test if we use GL quadrature
if ( nb_legendre == 1 ) then
  weight ( 1 )  = fourpi
  Rotxx = 1.0_dp ; Rotxy = 0.0_dp ; Rotxz = 0.0_dp
  Rotyx = 0.0_dp ; Rotyy = 1.0_dp ; Rotyz = 0.0_dp
  Rotzx = 0.0_dp ; Rotzy = 0.0_dp ; Rotzz = 1.0_dp
else if ( nb_legendre >= 2 ) then 
  do  n_theta = 1, nb_legendre !1 <= si nb_legendre=1
    cos_theta = x_legendre ( n_theta , nb_legendre ) !0
    sin_theta = sqrt ( 1.0_dp - cos_theta ** 2 ) !1
    do  n_phi = 1, 2*nb_legendre !1,2
      n_omega = ( n_theta - 1 ) * 2 * nb_legendre + n_phi !1,2
      phi = real ( n_phi - 1 , dp ) * twopi / real ( 2 * nb_legendre , dp ) !0,pi
      cos_phi = cos ( phi ) !1,-1
      sin_phi = sin ( phi ) !0,0
      OMx ( n_omega ) = sin_theta * cos_phi
      OMy ( n_omega ) = sin_theta * sin_phi
      OMz ( n_omega ) = cos_theta
      weight ( n_omega ) = w_legendre ( n_theta , nb_legendre ) * pi / real ( nb_legendre , dp ) ! 2pi,2pi
!Commented since we know use psi as a variable and use symmetries to reduce the number of angle
  !    do n_psi = 1, nb_psi !1
  !      psii = real(n_psi-1,dp) * pi / real(nb_psi,dp) !0
  
      do n_psi = 1 , nb_psi
         psii = x_psi(n_psi)
        
         cos_psi = cos(psii) !1
         sin_psi = sin(psii) !0
         
     
        Rotxx(n_omega,n_psi) =  cos_theta*cos_phi*cos_psi-sin_phi*sin_psi
        !Rotxx(1,1)           =  0
        !Rotxx(2,1)           =  0
        Rotxy(n_omega,n_psi) = -cos_theta*cos_phi*sin_psi-sin_phi*cos_psi
        !Rotxy(1,1)           =  0
        !Rotxy(2,1)           =  0
        Rotxz(n_omega,n_psi) =  sin_theta*cos_phi
        !Rotxz(1,1)           =  1
        !Rotxz(2,1)           =  -1
        Rotyx(n_omega,n_psi) =  cos_theta*sin_phi*cos_psi+cos_phi*sin_psi
        Rotyy(n_omega,n_psi) = -cos_theta*sin_phi*sin_psi+cos_phi*cos_psi
        Rotyz(n_omega,n_psi) =  sin_theta*sin_phi
        Rotzx(n_omega,n_psi) = -sin_theta*cos_psi
        Rotzy(n_omega,n_psi) =  sin_theta*sin_psi
        Rotzz(n_omega,n_psi) =  cos_theta
      end do
    end do
  end do
else
  write(*,*)'Error detected in compute_angular_grid.f90'
  stop
end if
! check if sum over all omega of weight(omega) is fourpi
if ( abs ( sum ( weight ( : ) ) - fourpi ) / fourpi > 1.0e-10_dp ) then
  print *, 'problem detected in compute_angular_grid.f90 :'
  print *, 'sum over omegas of weight(omega) is not 4pi. stop'
  stop
end if
if ( abs ( sum ( weight_psi ( : ) ) - twopi/sym_order )  > 1.0e-10_dp ) then
  print *, 'problem detected in compute_angular_grid.f90 :'
  print *, 'sum over omegas of weight(omega) is not 4pi. stop'
  stop
end if
end subroutine compute_angular_grid

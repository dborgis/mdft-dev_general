! Module for angular grid and Gauss-Legendre integration.
module quadrature

    use precision_kinds
    use system , only : nb_omega , nb_psi , nb_legendre 
    use constants , only : pi , twopi , fourpi

    implicit none

    real(dp), allocatable, dimension(:,:) :: w_legendre , x_legendre ! w(i,L) (weights) and x(i,L) (roots) for order L integration
    real(dp), allocatable, dimension(:) :: Omx , Omy , Omz , weight  ! unit vector for orientation OMEGA and associated weight
    real(dp), allocatable, dimension(:) :: weight_psi , x_psi
    real(dp), allocatable, dimension(:) :: x_leb, y_leb , z_leb , weight_leb 
    integer(i2b) :: sym_order
    
    contains
    
        subroutine deallocate_everything_gauss_legendre
            if ( allocated (weight_psi ) ) deallocate ( weight_psi)
            if ( allocated (x_psi ) ) deallocate ( x_psi)
            if ( allocated ( w_legendre ) ) deallocate ( w_legendre )
            if ( allocated ( x_legendre ) ) deallocate ( x_legendre ) 
            if ( allocated ( Omx ) ) deallocate ( Omx )
            if ( allocated ( Omy ) ) deallocate ( Omy )
            if ( allocated ( Omz ) ) deallocate ( Omz )
            if ( allocated ( x_leb ) ) deallocate ( y_leb)
            if ( allocated ( y_leb ) ) deallocate ( x_leb)
            if ( allocated ( z_leb ) ) deallocate ( z_leb)
            if ( allocated ( weight ) ) deallocate ( weight)
            if ( allocated ( weight_leb ) ) deallocate ( weight_leb)
        end subroutine
        
        ! Compute angular grid properties : Omx, Omy, Omz, weight
        subroutine gauss ( Rotxx, Rotxy, Rotxz, Rotyx, Rotyy, Rotyz, Rotzx, Rotzy, Rotzz)

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
            if ( abs( sum( weight(:) ) - fourpi ) / fourpi > 1.0e-10_dp &
                    .or.  abs ( sum ( weight_psi ( : ) ) - twopi/sym_order )  > 1.0e-10_dp ) then
                print *, 'problem detected in compute_angular_grid.f90 :'
                print *, 'sum over omegas of weight(omega) is not 4pi. stop'
                stop
            end if
        end subroutine

        subroutine lebedev( Rotxx, Rotxy, Rotxz, Rotyx, Rotyy, Rotyz, Rotzx, Rotzy, Rotzz)

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
        
            do n= 1, nb_omega
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

end module quadrature

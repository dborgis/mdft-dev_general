! Module for numerical integration
module quadrature

    use precision_kinds, only: dp, i2b
    use constants , only : pi, twopi, fourpi
    use input, only: input_log, input_char, input_int

    implicit none

    real(dp), allocatable, dimension(:), private :: x_leb, y_leb , z_leb
    real(dp), allocatable, dimension(:), public :: Omx , Omy , Omz  ! unit vector for orientation OMEGA and associated weight
    real(dp), allocatable, dimension(:,:), public :: Rotxx, Rotxy, Rotxz, Rotyx, Rotyy, Rotyz, Rotzx, Rotzy, Rotzz
    type angularGrid
        integer(i2b) :: n_angles
        real(dp), allocatable, dimension(:) :: weight, root
    end type
    type (angularGrid), public :: angGrid ! angular grid
    type (angularGrid), public :: molRotGrid ! rotation of molecule around its main axis, e.g., around C2v axis for H2O.
    type integrationScheme
        character(80) :: name
        integer(i2b) :: order
        real(dp), allocatable, dimension(:) :: weight, root
    end type
    type (integrationScheme), public :: intScheme
    integer(i2b), public :: sym_order

    contains
    
    
        subroutine init
            implicit none

            sym_order = input_int('sym_order')
            if ( sym_order <= 0 ) then
                print*,"in module_quadrature, sym_order, readen from input/dft.in, is unphysical. critical stop"
                stop
            end if

            ! molecular rotation grid
            call get_psi_integration_roots_and_weights (molRotGrid, sym_order)
            call check_weights_psi(molRotGrid%weight)

            ! integration scheme
            intScheme%name = trim(adjustl(input_char('quadrature')))
            if( intScheme%name /= 'L' .and. intScheme%name /= 'GL' ) then
                print*,'You ask for a integration scheme called ',intScheme%name
                print*,'it is not implemeted for now. Check readme for more information.'
                stop
            end if
            
            call read_order_of_quadrature (intScheme%order)
            if( intScheme%order <= 0 ) then
                print*,'You ask for a quadrature of order less than 1',intScheme%order
                STOP 'CRITICAL STOP. UNPHYSICAL QUADRATURE ORDER'
            end if
            
            allocate (intScheme%weight(intScheme%order), source=0._dp)
            allocate (intScheme%root(intScheme%order), source=0._dp)

            ! LEBEDEV
            if (intScheme%name=='L') then
                angGrid%n_angles = intScheme%order
                call allocate_Rotij (angGrid%n_angles,molRotGrid%n_angles,Rotxx,Rotxy,Rotxz,Rotyx,Rotyy,Rotyz,Rotzx,Rotzy,Rotzz)
                allocate( x_leb(intScheme%order), y_leb(intScheme%order), z_leb(intScheme%order) )
                allocate (angGrid%weight(angGrid%n_angles), source=0._dp)
                call lebedev_integration_roots_and_weights (intScheme%order, x_leb ,y_leb , z_leb, intScheme%weight)
                call lebedev (angGrid, intScheme, molRotGrid, Rotxx, Rotxy, Rotxz, Rotyx, Rotyy, Rotyz, Rotzx, Rotzy, Rotzz)
                deallocate( x_leb, y_leb, z_leb)

            ! GAUSS-LEGENDRE
            else if (intScheme%name=='GL') then
                select case (intScheme%order)
                case (1)
                    angGrid%n_angles = 1
                case default
                    angGrid%n_angles = 2*(intScheme%order**2)
                end select
                call allocate_Rotij (angGrid%n_angles,molRotGrid%n_angles,Rotxx,Rotxy,Rotxz,Rotyx,Rotyy,Rotyz,Rotzx,Rotzy,Rotzz)
                call gauss_legendre_integration_roots_and_weights (intScheme%order, intScheme%weight , intScheme%root)
                allocate (angGrid%weight(angGrid%n_angles), source=0._dp)
                call gauss (angGrid, intScheme, molRotGrid, Rotxx, Rotxy, Rotxz, Rotyx, Rotyy, Rotyz, Rotzx, Rotzy, Rotzz)
            end if
            
        end subroutine init

        ! GAUSS
        subroutine gauss (angGrid, intScheme, molRotGrid, Rotxx, Rotxy, Rotxz, Rotyx, Rotyy, Rotyz, Rotzx, Rotzy, Rotzz)
            type (angularGrid), intent(inout) :: angGrid
            type (angularGrid), intent(in) :: molRotGrid
            type (integrationScheme), intent(in) :: intScheme
            real(dp), dimension(:,:), intent(out) :: Rotxx,Rotxy,Rotxz,Rotyx,Rotyy,Rotyz,Rotzx,Rotzy,Rotzz
            integer(i2b) :: omega, n_psi, theta, n_phi
            real(dp) :: phi, cos_theta, sin_theta, cos_phi, sin_phi, cos_psi, sin_psi

            allocate(OMx(angGrid%n_angles))
            allocate(OMy(angGrid%n_angles))
            allocate(OMz(angGrid%n_angles))

            select case (intScheme%order)
            case (1)
                angGrid%weight (1)  = fourpi
                Rotxx = 1.0_dp ; Rotxy = 0.0_dp ; Rotxz = 0.0_dp
                Rotyx = 0.0_dp ; Rotyy = 1.0_dp ; Rotyz = 0.0_dp
                Rotzx = 0.0_dp ; Rotzy = 0.0_dp ; Rotzz = 1.0_dp
            case default
                do  theta = 1, intScheme%order
                    cos_theta = intScheme%root(theta)
                    sin_theta = sqrt ( 1.0_dp - cos_theta ** 2 )
                    do  n_phi = 1, 2*intScheme%order
                        omega = (theta-1)*2*intScheme%order + n_phi
                        phi = real(n_phi-1,dp) * twopi / real ( 2 * intScheme%order , dp ) !0,pi
                        cos_phi = cos ( phi ) ; sin_phi = sin ( phi )
                        OMx ( omega ) = sin_theta * cos_phi
                        OMy ( omega ) = sin_theta * sin_phi
                        OMz ( omega ) = cos_theta
                        angGrid%weight(omega) = intScheme%weight(theta) *pi /real(intScheme%order,dp)
                        do n_psi = 1 , molRotGrid%n_angles
                            cos_psi = cos(  molRotGrid%root(n_psi)  )
                            sin_psi = sin(  molRotGrid%root(n_psi)  )
                            Rotxx(omega,n_psi) =  cos_theta*cos_phi*cos_psi-sin_phi*sin_psi
                            Rotxy(omega,n_psi) = -cos_theta*cos_phi*sin_psi-sin_phi*cos_psi
                            Rotxz(omega,n_psi) =  sin_theta*cos_phi
                            Rotyx(omega,n_psi) =  cos_theta*sin_phi*cos_psi+cos_phi*sin_psi
                            Rotyy(omega,n_psi) = -cos_theta*sin_phi*sin_psi+cos_phi*cos_psi
                            Rotyz(omega,n_psi) =  sin_theta*sin_phi
                            Rotzx(omega,n_psi) = -sin_theta*cos_psi
                            Rotzy(omega,n_psi) =  sin_theta*sin_psi
                            Rotzz(omega,n_psi) =  cos_theta
                        end do
                    end do
                end do
            end select

            call check_weights(angGrid%weight)

        end subroutine gauss


        ! LEBEDEV
        subroutine lebedev (angGrid, intScheme, molRotGrid, Rotxx, Rotxy, Rotxz, Rotyx, Rotyy, Rotyz, Rotzx, Rotzy, Rotzz)
            
            type (angularGrid), intent(inout) :: angGrid
            type (angularGrid), intent(in) :: molRotGrid
            type (integrationScheme), intent(in) :: intScheme
            real(dp), dimension(:,:), intent(out) :: Rotxx,Rotxy,Rotxz,Rotyx,Rotyy,Rotyz,Rotzx,Rotzy,Rotzz
            integer(i2b) ::  n, psi
            real(dp) :: phi, theta, cos_theta, sin_theta, cos_phi, sin_phi, cos_psi, sin_psi

            allocate(OMx(angGrid%n_angles))
            allocate(OMy(angGrid%n_angles))
            allocate(OMz(angGrid%n_angles))
            
            angGrid%weight(1:angGrid%n_angles) = intScheme%weight(1:angGrid%n_angles)
            call check_weights(angGrid%weight)

            do concurrent (n=1:angGrid%n_angles)
                if ( x_leb(n)==0.0_dp .and. y_leb(n)==0.0_dp ) then 
                    phi = 0.0_dp
                else if (y_leb(n)>=0.0_dp) then 
                    phi = acos(x_leb(n)/(sqrt ( x_leb(n)**2 + y_leb(n)**2 ) ) )
                else
                    phi = twopi - acos(x_leb(n)/(sqrt ( x_leb(n)**2 + y_leb(n)**2 ) ) )
                end if
                theta = acos(z_leb(n))
                OMx (n) = sin(theta) * cos(phi)
                OMy (n) = sin(theta) * sin(phi)
                OMz (n) = cos(theta)
            end do
            
            do concurrent (n=1:angGrid%n_angles)
                theta=acos(z_leb(n))
                cos_theta=cos(theta)
                sin_theta=sin(theta)
                if    ( x_leb(n)  == 0.0_dp .and. y_leb(n)  == 0.0_dp ) then 
                    phi=0.0_dp
                else if  (y_leb(n)>=0.0_dp) then 
                    phi=acos(x_leb(n)/(sqrt ( x_leb(n)**2 + y_leb(n)**2 ) ) )
                else  
                    phi=twopi - acos(x_leb(n)/(sqrt ( x_leb(n)**2 + y_leb(n)**2 ) ) )
                end if
                cos_phi=cos(phi)
                sin_phi=sin(phi)    
                do concurrent (psi=1:molRotGrid%n_angles)
                    cos_psi = cos(  molRotGrid%root(psi)  )
                    sin_psi = sin(  molRotGrid%root(psi)  )
                    Rotxx(n,psi) =  cos_theta*cos_phi*cos_psi-sin_phi*sin_psi
                    Rotxy(n,psi) = -cos_theta*cos_phi*sin_psi-sin_phi*cos_psi
                    Rotxz(n,psi) =  sin_theta*cos_phi
                    Rotyx(n,psi) =  cos_theta*sin_phi*cos_psi+cos_phi*sin_psi
                    Rotyy(n,psi) = -cos_theta*sin_phi*sin_psi+cos_phi*cos_psi
                    Rotyz(n,psi) =  sin_theta*sin_phi
                    Rotzx(n,psi) = -sin_theta*cos_psi
                    Rotzy(n,psi) =  sin_theta*sin_psi
                    Rotzz(n,psi) =  cos_theta
                end do
            end do
        end subroutine lebedev


        subroutine check_weights_psi(weight_over_molecular_rotations)
            implicit none
            real(dp), dimension(:), intent(in) :: weight_over_molecular_rotations
            if ( abs ( sum ( weight_over_molecular_rotations( : ) ) - twopi/sym_order )  > 1.0e-10_dp ) then
                print*, 'problem detected in module_quadrature.f90 :'
                print*, 'sum over omegas of molRotGrid%weight(omega) is not 2pi/sym_order, it is ',&
                    sum ( weight_over_molecular_rotations )
                stop 'CRITICAL'
            end if
        end subroutine check_weights_psi

        
        subroutine check_weights(weight)
            implicit none
            real(dp), dimension(:), intent(in) :: weight
            ! check if sum over all omega of weight(omega) is fourpi
            if ( abs ( sum ( weight ( : ) ) - fourpi )  > 1.0e-10_dp ) then
                print *, 'problem detected in compute_angular_grid.f90 :'
                print *, 'sum over omegas of weight(omega) is not 4pi. it is ',sum ( weight )
                stop 'CRITICAL'
            end if            
        end subroutine check_weights
        
        

        subroutine get_psi_integration_roots_and_weights (molRotGrid, sym_order)
            type (angularGrid), intent(inout) :: molRotGrid
            integer(i2b), intent(out) :: sym_order
            integer(i2b) :: i
            molRotGrid%n_angles = input_int('nb_psi')
            if ( molRotGrid%n_angles < 1 ) then
                print*,"in module_quadrature, nb_psi, readen from input/dft.in is unphysical. critical stop."
                stop
            end if
            allocate( molRotGrid%weight(molRotGrid%n_angles), source=0._dp)
            molRotGrid%weight = twopi/real(molRotGrid%n_angles*sym_order,dp) ! equidistant repartition between 0 and 2Pi => all weights are equal
            allocate( molRotGrid%root(molRotGrid%n_angles), source=0._dp)
            do concurrent (i=1:molRotGrid%n_angles) ! equidistant repartition between 0 and 2Pi
                molRotGrid%root(i) = real(i-1,dp)*twopi/real(molRotGrid%n_angles*sym_order,dp) ! roots
            end do
        end subroutine get_psi_integration_roots_and_weights





        pure subroutine  allocate_Rotij (n_angles, n_psi,Rotxx,Rotxy,Rotxz,Rotyx,Rotyy,Rotyz,Rotzx,Rotzy,Rotzz)
            integer(i2b), intent(in) :: n_angles, n_psi
            real(dp), allocatable, dimension(:,:), intent(out) :: Rotxx,Rotxy,Rotxz,Rotyx,Rotyy,Rotyz,Rotzx,Rotzy,Rotzz
            allocate( Rotxx(n_angles,n_psi), Rotxy(n_angles,n_psi), Rotxz(n_angles,n_psi) )
            allocate( Rotyx(n_angles,n_psi), Rotyy(n_angles,n_psi), Rotyz(n_angles,n_psi) )
            allocate( Rotzx(n_angles,n_psi), Rotzy(n_angles,n_psi), Rotzz(n_angles,n_psi) )            
        end subroutine


        subroutine read_order_of_quadrature (order_of_quadrature)
            integer(i2b), intent(out) :: order_of_quadrature
            order_of_quadrature = input_int('order_of_quadrature')
            if (order_of_quadrature < 1) then
                print*,"In input/dft.in, I read order_of_quadrature is < 1. This is not physical as you must have"
                print*,"at least one point, even without angular grid. CRITICAL"
                stop
            end if
        end subroutine read_order_of_quadrature



end module quadrature

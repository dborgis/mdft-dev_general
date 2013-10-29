subroutine lebedev_integration_roots_and_weights (order, roots_x, roots_y, roots_z, weights)

!~     use quadrature , only :   weights, roots_x, roots_y, roots_z
    use precision_kinds , only : dp, i2b
!~     use input , only : input_int
!~     use constants , only : fourpi
!~     use system , only : nb_omega
    
    implicit none
    
    integer(i2b), intent(in) :: order
    real(dp), dimension(order), intent(out)  :: roots_x, roots_y, roots_z, weights
    
    real(dp), parameter :: fourpi = acos(-1._dp)*4._dp
    real(dp) :: a, b, v
    integer(i2b) ::  counter !dummy
!~ 
!~     nb_omega = input_int('order_of_quadrature')
!~     allocate (roots_x(nb_omega),roots_y(nb_omega),roots_z(nb_omega) )
!~     allocate (weights ( nb_omega ) )

    counter=0

    select case( order )
    case (6)
        V=1.0_dp/6.0_dp
        call GEN_OH( 1, roots_x, roots_y, roots_z, weights, A, B, V, counter)
    case (14)
        V=2.0_dp/30.0_dp
        call GEN_OH( 1, roots_x, roots_y, roots_z, weights, A, B, V, counter)
        V=3.0_dp/40.0_dp
        call GEN_OH( 3, roots_x, roots_y, roots_z, weights, A, B, V, counter)
    case (26)
        V=1.0_dp/21.0_dp!0.4761904761904762D-1
        call GEN_OH( 1, roots_x, roots_y, roots_z, weights, A, B, V, counter)
        V=4.0_dp/105.0_dp!0.3809523809523810D-1
        call GEN_OH( 2, roots_x, roots_y, roots_z, weights, A, B, V, counter)
        V=9.0_dp/280.0_dp!0.3214285714285714D-1
        call GEN_OH( 3, roots_x, roots_y, roots_z, weights, A, B, V, counter)
    case (38)
        V=1.0_dp/105.0_dp!0.9523809523809524D-2
        call GEN_OH( 1, roots_x, roots_y, roots_z, weights, A, B, V, counter)
        V=9.0_dp/280.0_dp!0.3214285714285714D-1
        call GEN_OH( 3, roots_x, roots_y, roots_z, weights, A, B, V, counter)
        A=0.4597008433809831D+0
        V=1.0_dp/35.0_dp!0.2857142857142857D-1
        call GEN_OH( 5, roots_x, roots_y, roots_z, weights, A, B, V, counter)
    case default
        stop "order selected for quadrature not implemented. Only 6, 14, 26 and 38 for now. STOP"
    end select
    
    weights = fourpi*weights


    contains
        
        pure subroutine gen_oh (code, roots_x, roots_y, roots_z, weights, a, b, v, counter)
            implicit none
            real(dp), dimension(:), intent(out) :: weights , roots_x , roots_y , roots_z
            real(dp), intent(inout) :: a, b, v
            integer(i2b), intent(inout) :: counter
            integer(i2b), intent(in) :: code
            select case (code)
            case (1)
                a=1.0_dp
                roots_x(counter+1) =  a
                roots_y(counter+1) =  0.0_dp
                roots_z(counter+1) =  0.0_dp
                weights(counter+1) =  v
                roots_x(counter+2) = -a
                roots_y(counter+2) =  0.0_dp
                roots_z(counter+2) =  0.0_dp
                weights(counter+2) =  v
                roots_x(counter+3) =  0.0_dp
                roots_y(counter+3) =  a
                roots_z(counter+3) =  0.0_dp
                weights(counter+3) =  v
                roots_x(counter+4) =  0.0_dp
                roots_y(counter+4) = -a
                roots_z(counter+4) =  0.0_dp
                weights(counter+4) =  v
                roots_x(counter+5) =  0.0_dp
                roots_y(counter+5) =  0.0_dp
                roots_z(counter+5) =  a
                weights(counter+5) =  v
                roots_x(counter+6) =  0.0_dp
                roots_y(counter+6) =  0.0_dp
                roots_z(counter+6) = -a
                weights(counter+6) =  v
                counter=counter+6
            case (2)
                a=sqrt(0.5_dp)
                roots_x(counter+1) =  0.0_dp
                roots_y(counter+1) =  a
                roots_z(counter+1) =  a
                weights(counter+1) =  v
                roots_x(counter+2) =  0.0_dp
                roots_y(counter+2) = -a
                roots_z(counter+2) =  a
                weights(counter+2) =  v
                roots_x(counter+3) =  0.0_dp
                roots_y(counter+3) =  a
                roots_z(counter+3) = -a
                weights(counter+3) =  v
                roots_x(counter+4) =  0.0_dp
                roots_y(counter+4) = -a
                roots_z(counter+4) = -a
                weights(counter+4) =  v
                roots_x(counter+5) =  a
                roots_y(counter+5) =  0.0_dp
                roots_z(counter+5) =  a
                weights(counter+5) =  v
                roots_x(counter+6) = -a
                roots_y(counter+6) =  0.0_dp
                roots_z(counter+6) =  a
                weights(counter+6) =  v
                roots_x(counter+7) =  a
                roots_y(counter+7) =  0.0_dp
                roots_z(counter+7) = -a
                weights(counter+7) =  v
                roots_x(counter+8) = -a
                roots_y(counter+8) =  0.0_dp
                roots_z(counter+8) = -a
                weights(counter+8) =  v
                roots_x(counter+9) =  a
                roots_y(counter+9) =  a
                roots_z(counter+9) =  0.0_dp
                weights(counter+9) =  v
                roots_x(counter+10) = -a
                roots_y(counter+10) =  a
                roots_z(counter+10) =  0.0_dp
                weights(counter+10) =  v
                roots_x(counter+11) =  a
                roots_y(counter+11) = -a
                roots_z(counter+11) =  0.0_dp
                weights(counter+11) =  v
                roots_x(counter+12) = -a
                roots_y(counter+12) = -a
                roots_z(counter+12) =  0.0_dp
                weights(counter+12) =  v
                counter=counter+12
            case (3)
                a = sqrt(1d0/3d0)
                roots_x(counter+1) =  a
                roots_y(counter+1) =  a
                roots_z(counter+1) =  a
                weights(counter+1) =  v
                roots_x(counter+2) = -a
                roots_y(counter+2) =  a
                roots_z(counter+2) =  a
                weights(counter+2) =  v
                roots_x(counter+3) =  a
                roots_y(counter+3) = -a
                roots_z(counter+3) =  a
                weights(counter+3) =  v
                roots_x(counter+4) = -a
                roots_y(counter+4) = -a
                roots_z(counter+4) =  a
                weights(counter+4) =  v
                roots_x(counter+5) =  a
                roots_y(counter+5) =  a
                roots_z(counter+5) = -a
                weights(counter+5) =  v
                roots_x(counter+6) = -a
                roots_y(counter+6) =  a
                roots_z(counter+6) = -a
                weights(counter+6) =  v
                roots_x(counter+7) =  a
                roots_y(counter+7) = -a
                roots_z(counter+7) = -a
                weights(counter+7) =  v
                roots_x(counter+8) = -a
                roots_y(counter+8) = -a
                roots_z(counter+8) = -a
                weights(counter+8) =  v
                counter=counter+8
            case (4)
                b = sqrt(1d0 - 2d0*a*a)
                roots_x(counter+1) =  a
                roots_y(counter+1) =  a
                roots_z(counter+1) =  b
                weights(counter+1) =  v
                roots_x(counter+2) = -a
                roots_y(counter+2) =  a
                roots_z(counter+2) =  b
                weights(counter+2) =  v
                roots_x(counter+3) =  a
                roots_y(counter+3) = -a
                roots_z(counter+3) =  b
                weights(counter+3) =  v
                roots_x(counter+4) = -a
                roots_y(counter+4) = -a
                roots_z(counter+4) =  b
                weights(counter+4) =  v
                roots_x(counter+5) =  a
                roots_y(counter+5) =  a
                roots_z(counter+5) = -b
                weights(counter+5) =  v
                roots_x(counter+6) = -a
                roots_y(counter+6) =  a
                roots_z(counter+6) = -b
                weights(counter+6) =  v
                roots_x(counter+7) =  a
                roots_y(counter+7) = -a
                roots_z(counter+7) = -b
                weights(counter+7) =  v
                roots_x(counter+8) = -a
                roots_y(counter+8) = -a
                roots_z(counter+8) = -b
                weights(counter+8) =  v
                roots_x(counter+9) =  a
                roots_y(counter+9) =  b
                roots_z(counter+9) =  a
                weights(counter+9) =  v
                roots_x(counter+10) = -a
                roots_y(counter+10) =  b
                roots_z(counter+10) =  a
                weights(counter+10) =  v
                roots_x(counter+11) =  a
                roots_y(counter+11) = -b
                roots_z(counter+11) =  a
                weights(counter+11) =  v
                roots_x(counter+12) = -a
                roots_y(counter+12) = -b
                roots_z(counter+12) =  a
                weights(counter+12) =  v
                roots_x(counter+13) =  a
                roots_y(counter+13) =  b
                roots_z(counter+13) = -a
                weights(counter+13) =  v
                roots_x(counter+14) = -a
                roots_y(counter+14) =  b
                roots_z(counter+14) = -a
                weights(counter+14) =  v
                roots_x(counter+15) =  a
                roots_y(counter+15) = -b
                roots_z(counter+15) = -a
                weights(counter+15) =  v
                roots_x(counter+16) = -a
                roots_y(counter+16) = -b
                roots_z(counter+16) = -a
                weights(counter+16) =  v
                roots_x(counter+17) =  b
                roots_y(counter+17) =  a
                roots_z(counter+17) =  a
                weights(counter+17) =  v
                roots_x(counter+18) = -b
                roots_y(counter+18) =  a
                roots_z(counter+18) =  a
                weights(counter+18) =  v
                roots_x(counter+19) =  b
                roots_y(counter+19) = -a
                roots_z(counter+19) =  a
                weights(counter+19) =  v
                roots_x(counter+20) = -b
                roots_y(counter+20) = -a
                roots_z(counter+20) =  a
                weights(counter+20) =  v
                roots_x(counter+21) =  b
                roots_y(counter+21) =  a
                roots_z(counter+21) = -a
                weights(counter+21) =  v
                roots_x(counter+22) = -b
                roots_y(counter+22) =  a
                roots_z(counter+22) = -a
                weights(counter+22) =  v
                roots_x(counter+23) =  b
                roots_y(counter+23) = -a
                roots_z(counter+23) = -a
                weights(counter+23) =  v
                roots_x(counter+24) = -b
                roots_y(counter+24) = -a
                roots_z(counter+24) = -a
                weights(counter+24) =  v
                counter=counter+24
            case (5)
                b=sqrt(1d0-a*a)
                roots_x(counter+1) =  a
                roots_y(counter+1) =  b
                roots_z(counter+1) =  0d0
                weights(counter+1) =  v
                roots_x(counter+2) = -a
                roots_y(counter+2) =  b
                roots_z(counter+2) =  0d0
                weights(counter+2) =  v
                roots_x(counter+3) =  a
                roots_y(counter+3) = -b
                roots_z(counter+3) =  0d0
                weights(counter+3) =  v
                roots_x(counter+4) = -a
                roots_y(counter+4) = -b
                roots_z(counter+4) =  0d0
                weights(counter+4) =  v
                roots_x(counter+5) =  b
                roots_y(counter+5) =  a
                roots_z(counter+5) =  0d0
                weights(counter+5) =  v
                roots_x(counter+6) = -b
                roots_y(counter+6) =  a
                roots_z(counter+6) =  0d0
                weights(counter+6) =  v
                roots_x(counter+7) =  b
                roots_y(counter+7) = -a
                roots_z(counter+7) =  0d0
                weights(counter+7) =  v
                roots_x(counter+8) = -b
                roots_y(counter+8) = -a
                roots_z(counter+8) =  0d0
                weights(counter+8) =  v
                roots_x(counter+9) =  a
                roots_y(counter+9) =  0d0
                roots_z(counter+9) =  b
                weights(counter+9) =  v
                roots_x(counter+10) = -a
                roots_y(counter+10) =  0d0
                roots_z(counter+10) =  b
                weights(counter+10) =  v
                roots_x(counter+11) =  a
                roots_y(counter+11) =  0d0
                roots_z(counter+11) = -b
                weights(counter+11) =  v
                roots_x(counter+12) = -a
                roots_y(counter+12) =  0d0
                roots_z(counter+12) = -b
                weights(counter+12) =  v
                roots_x(counter+13) =  b
                roots_y(counter+13) =  0d0
                roots_z(counter+13) =  a
                weights(counter+13) =  v
                roots_x(counter+14) = -b
                roots_y(counter+14) =  0d0
                roots_z(counter+14) =  a
                weights(counter+14) =  v
                roots_x(counter+15) =  b
                roots_y(counter+15) =  0d0
                roots_z(counter+15) = -a
                weights(counter+15) =  v
                roots_x(counter+16) = -b
                roots_y(counter+16) =  0d0
                roots_z(counter+16) = -a
                weights(counter+16) =  v
                roots_x(counter+17) =  0d0
                roots_y(counter+17) =  a
                roots_z(counter+17) =  b
                weights(counter+17) =  v
                roots_x(counter+18) =  0d0
                roots_y(counter+18) = -a
                roots_z(counter+18) =  b
                weights(counter+18) =  v
                roots_x(counter+19) =  0d0
                roots_y(counter+19) =  a
                roots_z(counter+19) = -b
                weights(counter+19) =  v
                roots_x(counter+20) =  0d0
                roots_y(counter+20) = -a
                roots_z(counter+20) = -b
                weights(counter+20) =  v
                roots_x(counter+21) =  0d0
                roots_y(counter+21) =  b
                roots_z(counter+21) =  a
                weights(counter+21) =  v
                roots_x(counter+22) =  0d0
                roots_y(counter+22) = -b
                roots_z(counter+22) =  a
                weights(counter+22) =  v
                roots_x(counter+23) =  0d0
                roots_y(counter+23) =  b
                roots_z(counter+23) = -a
                weights(counter+23) =  v
                roots_x(counter+24) =  0d0
                roots_y(counter+24) = -b
                roots_z(counter+24) = -a
                weights(counter+24) =  v
                counter=counter+24
            end select
        end subroutine gen_oh

end subroutine lebedev_integration_roots_and_weights 

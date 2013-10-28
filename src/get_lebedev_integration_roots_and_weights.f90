subroutine get_lebedev_integration_roots_and_weights 

    use quadrature , only :   weight_leb, x_leb, y_leb, z_leb
    use precision_kinds , only : dp, i2b
    use input , only : input_int
    use constants , only : fourpi
    use system , only : nb_omega
    
    implicit none
    
    real (dp) :: a, b, v
    integer(i2b) ::  counter !dummy

    nb_omega = input_int('order_of_quadrature')
    allocate (x_leb(nb_omega),y_leb(nb_omega),z_leb(nb_omega) )
    allocate (weight_leb ( nb_omega ) )

    counter=0

    select case( nb_omega )
        case (6)
            V=1.0_dp/6.0_dp
            call GEN_OH( 1, x_leb, y_leb, z_leb, weight_leb, A, B, V)
        case (14)
            V=2.0_dp/30.0_dp
            Call GEN_OH( 1, x_leb, y_leb, z_leb, weight_leb, A, B, V)
            V=3.0_dp/40.0_dp
            Call GEN_OH( 3 , x_leb, y_leb, z_leb, weight_leb, A, B, V)
        case (26)
            V=1.0_dp/21.0_dp!0.4761904761904762D-1
            Call GEN_OH( 1, x_leb, y_leb, z_leb, weight_leb, A, B, V)
            V=4.0_dp/105.0_dp!0.3809523809523810D-1
            Call GEN_OH( 2, x_leb, y_leb, z_leb, weight_leb, A, B, V)
            V=9.0_dp/280.0_dp!0.3214285714285714D-1
            Call GEN_OH( 3, x_leb, y_leb, z_leb, weight_leb, A, B, V)
        case (38)
            V=1.0_dp/105.0_dp!0.9523809523809524D-2
            Call GEN_OH( 1, x_leb, y_leb, z_leb, weight_leb, A, B, V)
            V=9.0_dp/280.0_dp!0.3214285714285714D-1
            Call GEN_OH( 3, x_leb, y_leb, z_leb, weight_leb, A, B, V)
            A=0.4597008433809831D+0
            V=1.0_dp/35.0_dp!0.2857142857142857D-1
            Call GEN_OH( 5, x_leb, y_leb, z_leb, weight_leb, A, B, V)
    end select
    weight_leb = fourpi*weight_leb  

    open(11, file='output/lebedev_quadrature.dat')
        write(11,*) "x_leb = ",x_leb
        write(11,*) "y_leb = ",y_leb
        write(11,*) "z_leb = ",z_leb
        write(11,*) "weight_leb = ",weight_leb 
    close(11)

    contains
        
        subroutine gen_oh(code, x_leb, y_leb, z_leb, weight_leb, a, b, v)
            implicit none
            real(dp), dimension(:), intent(out) :: weight_leb , x_leb , y_leb , z_leb
            real(dp), intent(inout) :: a, b, v
            integer(i2b), intent(in) :: code
            select case (code)
            case (1)
                a=1.0_dp
                x_leb(counter+1) =  a
                y_leb(counter+1) =  0.0_dp
                z_leb(counter+1) =  0.0_dp
                weight_leb(counter+1) =  v
                x_leb(counter+2) = -a
                y_leb(counter+2) =  0.0_dp
                z_leb(counter+2) =  0.0_dp
                weight_leb(counter+2) =  v
                x_leb(counter+3) =  0.0_dp
                y_leb(counter+3) =  a
                z_leb(counter+3) =  0.0_dp
                weight_leb(counter+3) =  v
                x_leb(counter+4) =  0.0_dp
                y_leb(counter+4) = -a
                z_leb(counter+4) =  0.0_dp
                weight_leb(counter+4) =  v
                x_leb(counter+5) =  0.0_dp
                y_leb(counter+5) =  0.0_dp
                z_leb(counter+5) =  a
                weight_leb(counter+5) =  v
                x_leb(counter+6) =  0.0_dp
                y_leb(counter+6) =  0.0_dp
                z_leb(counter+6) = -a
                weight_leb(counter+6) =  v
                counter=counter+6
            case (2)
                a=sqrt(0.5_dp)
                x_leb(counter+1) =  0.0_dp
                y_leb(counter+1) =  a
                z_leb(counter+1) =  a
                weight_leb(counter+1) =  v
                x_leb(counter+2) =  0.0_dp
                y_leb(counter+2) = -a
                z_leb(counter+2) =  a
                weight_leb(counter+2) =  v
                x_leb(counter+3) =  0.0_dp
                y_leb(counter+3) =  a
                z_leb(counter+3) = -a
                weight_leb(counter+3) =  v
                x_leb(counter+4) =  0.0_dp
                y_leb(counter+4) = -a
                z_leb(counter+4) = -a
                weight_leb(counter+4) =  v
                x_leb(counter+5) =  a
                y_leb(counter+5) =  0.0_dp
                z_leb(counter+5) =  a
                weight_leb(counter+5) =  v
                x_leb(counter+6) = -a
                y_leb(counter+6) =  0.0_dp
                z_leb(counter+6) =  a
                weight_leb(counter+6) =  v
                x_leb(counter+7) =  a
                y_leb(counter+7) =  0.0_dp
                z_leb(counter+7) = -a
                weight_leb(counter+7) =  v
                x_leb(counter+8) = -a
                y_leb(counter+8) =  0.0_dp
                z_leb(counter+8) = -a
                weight_leb(counter+8) =  v
                x_leb(counter+9) =  a
                y_leb(counter+9) =  a
                z_leb(counter+9) =  0.0_dp
                weight_leb(counter+9) =  v
                x_leb(counter+10) = -a
                y_leb(counter+10) =  a
                z_leb(counter+10) =  0.0_dp
                weight_leb(counter+10) =  v
                x_leb(counter+11) =  a
                y_leb(counter+11) = -a
                z_leb(counter+11) =  0.0_dp
                weight_leb(counter+11) =  v
                x_leb(counter+12) = -a
                y_leb(counter+12) = -a
                z_leb(counter+12) =  0.0_dp
                weight_leb(counter+12) =  v
                counter=counter+12
            case (3)
                a = sqrt(1d0/3d0)
                x_leb(counter+1) =  a
                y_leb(counter+1) =  a
                z_leb(counter+1) =  a
                weight_leb(counter+1) =  v
                x_leb(counter+2) = -a
                y_leb(counter+2) =  a
                z_leb(counter+2) =  a
                weight_leb(counter+2) =  v
                x_leb(counter+3) =  a
                y_leb(counter+3) = -a
                z_leb(counter+3) =  a
                weight_leb(counter+3) =  v
                x_leb(counter+4) = -a
                y_leb(counter+4) = -a
                z_leb(counter+4) =  a
                weight_leb(counter+4) =  v
                x_leb(counter+5) =  a
                y_leb(counter+5) =  a
                z_leb(counter+5) = -a
                weight_leb(counter+5) =  v
                x_leb(counter+6) = -a
                y_leb(counter+6) =  a
                z_leb(counter+6) = -a
                weight_leb(counter+6) =  v
                x_leb(counter+7) =  a
                y_leb(counter+7) = -a
                z_leb(counter+7) = -a
                weight_leb(counter+7) =  v
                x_leb(counter+8) = -a
                y_leb(counter+8) = -a
                z_leb(counter+8) = -a
                weight_leb(counter+8) =  v
                counter=counter+8
            case (4)
                b = sqrt(1d0 - 2d0*a*a)
                x_leb(counter+1) =  a
                y_leb(counter+1) =  a
                z_leb(counter+1) =  b
                weight_leb(counter+1) =  v
                x_leb(counter+2) = -a
                y_leb(counter+2) =  a
                z_leb(counter+2) =  b
                weight_leb(counter+2) =  v
                x_leb(counter+3) =  a
                y_leb(counter+3) = -a
                z_leb(counter+3) =  b
                weight_leb(counter+3) =  v
                x_leb(counter+4) = -a
                y_leb(counter+4) = -a
                z_leb(counter+4) =  b
                weight_leb(counter+4) =  v
                x_leb(counter+5) =  a
                y_leb(counter+5) =  a
                z_leb(counter+5) = -b
                weight_leb(counter+5) =  v
                x_leb(counter+6) = -a
                y_leb(counter+6) =  a
                z_leb(counter+6) = -b
                weight_leb(counter+6) =  v
                x_leb(counter+7) =  a
                y_leb(counter+7) = -a
                z_leb(counter+7) = -b
                weight_leb(counter+7) =  v
                x_leb(counter+8) = -a
                y_leb(counter+8) = -a
                z_leb(counter+8) = -b
                weight_leb(counter+8) =  v
                x_leb(counter+9) =  a
                y_leb(counter+9) =  b
                z_leb(counter+9) =  a
                weight_leb(counter+9) =  v
                x_leb(counter+10) = -a
                y_leb(counter+10) =  b
                z_leb(counter+10) =  a
                weight_leb(counter+10) =  v
                x_leb(counter+11) =  a
                y_leb(counter+11) = -b
                z_leb(counter+11) =  a
                weight_leb(counter+11) =  v
                x_leb(counter+12) = -a
                y_leb(counter+12) = -b
                z_leb(counter+12) =  a
                weight_leb(counter+12) =  v
                x_leb(counter+13) =  a
                y_leb(counter+13) =  b
                z_leb(counter+13) = -a
                weight_leb(counter+13) =  v
                x_leb(counter+14) = -a
                y_leb(counter+14) =  b
                z_leb(counter+14) = -a
                weight_leb(counter+14) =  v
                x_leb(counter+15) =  a
                y_leb(counter+15) = -b
                z_leb(counter+15) = -a
                weight_leb(counter+15) =  v
                x_leb(counter+16) = -a
                y_leb(counter+16) = -b
                z_leb(counter+16) = -a
                weight_leb(counter+16) =  v
                x_leb(counter+17) =  b
                y_leb(counter+17) =  a
                z_leb(counter+17) =  a
                weight_leb(counter+17) =  v
                x_leb(counter+18) = -b
                y_leb(counter+18) =  a
                z_leb(counter+18) =  a
                weight_leb(counter+18) =  v
                x_leb(counter+19) =  b
                y_leb(counter+19) = -a
                z_leb(counter+19) =  a
                weight_leb(counter+19) =  v
                x_leb(counter+20) = -b
                y_leb(counter+20) = -a
                z_leb(counter+20) =  a
                weight_leb(counter+20) =  v
                x_leb(counter+21) =  b
                y_leb(counter+21) =  a
                z_leb(counter+21) = -a
                weight_leb(counter+21) =  v
                x_leb(counter+22) = -b
                y_leb(counter+22) =  a
                z_leb(counter+22) = -a
                weight_leb(counter+22) =  v
                x_leb(counter+23) =  b
                y_leb(counter+23) = -a
                z_leb(counter+23) = -a
                weight_leb(counter+23) =  v
                x_leb(counter+24) = -b
                y_leb(counter+24) = -a
                z_leb(counter+24) = -a
                weight_leb(counter+24) =  v
                counter=counter+24
            case (5)
                b=sqrt(1d0-a*a)
                x_leb(counter+1) =  a
                y_leb(counter+1) =  b
                z_leb(counter+1) =  0d0
                weight_leb(counter+1) =  v
                x_leb(counter+2) = -a
                y_leb(counter+2) =  b
                z_leb(counter+2) =  0d0
                weight_leb(counter+2) =  v
                x_leb(counter+3) =  a
                y_leb(counter+3) = -b
                z_leb(counter+3) =  0d0
                weight_leb(counter+3) =  v
                x_leb(counter+4) = -a
                y_leb(counter+4) = -b
                z_leb(counter+4) =  0d0
                weight_leb(counter+4) =  v
                x_leb(counter+5) =  b
                y_leb(counter+5) =  a
                z_leb(counter+5) =  0d0
                weight_leb(counter+5) =  v
                x_leb(counter+6) = -b
                y_leb(counter+6) =  a
                z_leb(counter+6) =  0d0
                weight_leb(counter+6) =  v
                x_leb(counter+7) =  b
                y_leb(counter+7) = -a
                z_leb(counter+7) =  0d0
                weight_leb(counter+7) =  v
                x_leb(counter+8) = -b
                y_leb(counter+8) = -a
                z_leb(counter+8) =  0d0
                weight_leb(counter+8) =  v
                x_leb(counter+9) =  a
                y_leb(counter+9) =  0d0
                z_leb(counter+9) =  b
                weight_leb(counter+9) =  v
                x_leb(counter+10) = -a
                y_leb(counter+10) =  0d0
                z_leb(counter+10) =  b
                weight_leb(counter+10) =  v
                x_leb(counter+11) =  a
                y_leb(counter+11) =  0d0
                z_leb(counter+11) = -b
                weight_leb(counter+11) =  v
                x_leb(counter+12) = -a
                y_leb(counter+12) =  0d0
                z_leb(counter+12) = -b
                weight_leb(counter+12) =  v
                x_leb(counter+13) =  b
                y_leb(counter+13) =  0d0
                z_leb(counter+13) =  a
                weight_leb(counter+13) =  v
                x_leb(counter+14) = -b
                y_leb(counter+14) =  0d0
                z_leb(counter+14) =  a
                weight_leb(counter+14) =  v
                x_leb(counter+15) =  b
                y_leb(counter+15) =  0d0
                z_leb(counter+15) = -a
                weight_leb(counter+15) =  v
                x_leb(counter+16) = -b
                y_leb(counter+16) =  0d0
                z_leb(counter+16) = -a
                weight_leb(counter+16) =  v
                x_leb(counter+17) =  0d0
                y_leb(counter+17) =  a
                z_leb(counter+17) =  b
                weight_leb(counter+17) =  v
                x_leb(counter+18) =  0d0
                y_leb(counter+18) = -a
                z_leb(counter+18) =  b
                weight_leb(counter+18) =  v
                x_leb(counter+19) =  0d0
                y_leb(counter+19) =  a
                z_leb(counter+19) = -b
                weight_leb(counter+19) =  v
                x_leb(counter+20) =  0d0
                y_leb(counter+20) = -a
                z_leb(counter+20) = -b
                weight_leb(counter+20) =  v
                x_leb(counter+21) =  0d0
                y_leb(counter+21) =  b
                z_leb(counter+21) =  a
                weight_leb(counter+21) =  v
                x_leb(counter+22) =  0d0
                y_leb(counter+22) = -b
                z_leb(counter+22) =  a
                weight_leb(counter+22) =  v
                x_leb(counter+23) =  0d0
                y_leb(counter+23) =  b
                z_leb(counter+23) = -a
                weight_leb(counter+23) =  v
                x_leb(counter+24) =  0d0
                y_leb(counter+24) = -b
                z_leb(counter+24) = -a
                weight_leb(counter+24) =  v
                counter=counter+24
            end select
        end subroutine gen_oh

end subroutine get_lebedev_integration_roots_and_weights 

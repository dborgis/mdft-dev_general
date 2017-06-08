subroutine output_gOfRAndCosThetaAndPsi

    use module_grid, only: grid
    use precision_kinds, only: dp
    use module_solute, only: solute
    use module_solvent, only: solvent
    use mathematica,     only: deduce_optimal_histogram_properties
   
    implicit none

    real(dp), parameter :: zerodp = 0.0_dp, onedp = 1.0_dp
    real(dp), parameter :: pi = acos(-1.0_dp)
    integer :: ix, iy ,iz, io, ncostheta, icostheta, ir, nr
    real(dp) :: cosTheta, OM(3), normOM, H1(3), H2(3), Hbar(3), normHbar
    real(dp), parameter :: Oz(3) = [ 0., 0., 1. ]
    real(dp) :: costheta_low, costheta_max, dcostheta, dr, r_low, r_max
    real(dp), allocatable :: g(:,:), bin(:,:)


! Check that it is meaningfull to enter this routine, which is the case only if the solute is made of one site.
block
    integer :: countSoluteSites
    countSoluteSites = size( solute%site(:) )
    if( countSoluteSites > 1 ) then
        return ! because the number of solute site is not 1, it has no sense to make spherically symmetric approximation
    end if
end block




    nr = 100
    dr = min(grid%lz,grid%ly, grid%lz)/2. / nr
    ncostheta = 20
    dcostheta = 2. / ncostheta
    allocate(g(nr, ncostheta))
    allocate(bin(nr,ncostheta))
    g = 0.
    bin = 0

    do io = 1, grid%no
        H1(1) = dot_product( [grid%rotxx(io), grid%rotxy(io), grid%rotxz(io) ], solvent(1)%site(2)%r ) ! first H is solvent site 2
        H1(2) = dot_product( [grid%rotyx(io), grid%rotyy(io), grid%rotyz(io) ], solvent(1)%site(2)%r )
        H1(3) = dot_product( [grid%rotzx(io), grid%rotzy(io), grid%rotzz(io) ], solvent(1)%site(2)%r )
        H2(1) = dot_product( [grid%rotxx(io), grid%rotxy(io), grid%rotxz(io) ], solvent(1)%site(3)%r ) ! second H is solvent site 3
        H2(2) = dot_product( [grid%rotyx(io), grid%rotyy(io), grid%rotyz(io) ], solvent(1)%site(3)%r )
        H2(3) = dot_product( [grid%rotzx(io), grid%rotzy(io), grid%rotzz(io) ], solvent(1)%site(3)%r )
        Hbar(1) = (H1(1) + H2(1)) /2. ! Hbar is the coordinates of the barycenter
        Hbar(2) = (H1(2) + H2(2)) /2.
        Hbar(3) = (H1(3) + H2(3)) /2.
        normHbar = norm2(Hbar)
    
        do concurrent( ix = 1: grid%nx, iy = 1: grid%ny, iz = 1: grid%nz)
            OM = [ ix - 1, iy - 1, iz - 1 ] * [ grid%dx, grid%dy, grid%dz ] - solute%site(1)%r(:)
            normOM = norm2( OM )
            if( normOM < 0.0001 ) cycle
            cosTheta = dot_product( OM, Hbar ) / (normOM * normHbar)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!>>>>>>>>>>>>>>>>>>>><
block
real(dp) :: sinTheta, phi, cosPhi, sinPhi, u_kz, v_kz, psi
real(dp), parameter :: twopi = 2._dp * acos(-1._dp)
            sinTheta = sqrt( 1._dp - cosTheta **2)
            phi = angle( OM(2), OM(3) )
            cosPhi = cos(phi)
            sinPhi = sin(phi)
            u_kz = grid%rotxx(io)*sinTheta*cosPhi + grid%rotyx(io)*sinTheta*sinPhi + grid%rotzx(io)*cosTheta
            v_kz = grid%rotxy(io)*sinTheta*cosPhi + grid%rotyy(io)*sinTheta*sinPhi + grid%rotzy(io)*cosTheta
            psi = modulo( angle(-u_kz,v_kz),twopi/grid%molRotSymOrder )
end block
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!<<<<<<<<<<<<<<<<<<
            do icostheta = 1, ncostheta
                costheta_low = (icostheta-1) * dcostheta - 1.
                costheta_max = costheta_low + dcostheta
            
                do ir = 1, nr
                    r_low = (ir - 1) * dr
                    r_max = r_low + dr

                    if(       costheta >= costheta_low .and. costheta < costheta_max &
                        .and. normOM >= r_low .and. normOM <  r_max ) then

                        g(ir,icostheta) = g(ir,icostheta) + solvent(1)%xi(io,ix,iy,iz)**2 * grid%w(io)
                        bin(ir,icostheta) = bin(ir,icostheta) + 1
                    end if

                end do
            end do
        end do
    end do

    open( 77, file="output/rdf-theta.dat")
    where (bin /= 0) g = g / bin
    do icostheta = 1, ncostheta
        costheta_low = (icostheta-1) * dcostheta - 1.
        costheta_max = costheta_low + dcostheta
        do ir = 1, nr
            r_low = (ir - 1) * dr
            r_max = r_low + dr
            write(77,*) (r_low+r_max)/2., ( costheta_low + costheta_max )/2., g(ir,icostheta)
        end do
        write(77,*)
    end do
    close(77)

! gnuplot> set pm3d; set border 4095; set xlabel "r (Ang)"; set ylabel "cos theta"; set zlabel "g"; unset key; splot "output/rdf-theta.dat" w l

    stop "VICTORYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYY"


contains
    PURE FUNCTION angle(x,y)
        use precision_kinds,    ONLY: dp
        use constants,          ONLY: twopi
        IMPLICIT NONE
        REAL(dp), INTENT(IN)    :: x,y
        REAL(dp)                :: angle
        REAL(dp)                :: xx
        IF (x==0._dp .AND. y==0._dp) THEN
            angle = 0._dp
        ELSE
            xx = ACOS( x/SQRT(x**2 + y**2) )
            IF (y>=0._dp) THEN
                angle = xx
            ELSE
                angle = twopi - xx
            END IF
        END IF
    END FUNCTION angle

end subroutine output_gOfRAndCosThetaAndPsi
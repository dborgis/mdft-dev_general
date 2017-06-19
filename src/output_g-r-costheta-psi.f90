subroutine output_gOfRAndCosThetaAndPsi

    use module_grid, only: grid
    use precision_kinds, only: dp
    use module_solute, only: solute
    use module_solvent, only: solvent
    use mathematica,     only: deduce_optimal_histogram_properties
    use module_rotation, only: angle, omega_prime_fonction_de_omega
   
    implicit none

    real(dp), parameter :: zerodp = 0.0_dp, onedp = 1.0_dp
    real(dp), parameter :: pi = acos(-1.0_dp)
    integer :: ix, iy ,iz, io, ncostheta, icostheta, ir, nr, npsi, ipsi
    real(dp) :: cosTheta, OM(3), normOM, H1(3), H2(3), Hbar(3), normHbar
    real(dp), parameter :: Oz(3) = [ 0., 0., 1. ], twopi = acos(-1._dp)
    real(dp) :: costheta_low, costheta_max, dcostheta, dr, r_low, r_max, r_mean, omega(3), omegaPrime(3)
    real(dp), allocatable :: g(:,:,:), bin(:,:,:)
    real(dp) :: dpsi, psi, psi_low, psi_max


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
    ncostheta = 25
    dcostheta = 2. / ncostheta
    npsi = 25
    dpsi = twopi / npsi
    allocate(g(nr, ncostheta, npsi))
    allocate(bin(nr,ncostheta,npsi))
    g = 0.
    bin = 0

    do io = 1, grid%no
        omega = [ grid%theta(io), grid%phi(io), grid%psi(io) ]
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

            omegaPrime = omega_prime_fonction_de_omega( OM, omega ) ! omega' est le vecteur de theta',phi',psi'   ie les nouveaux theta, phi, psi dans le repère local soluté-solvant
            cosTheta = cos(omegaPrime(1))
            psi = omegaPrime(3)

            do ipsi = 1, npsi
                psi_low = (ipsi-1) * dpsi
                psi_max = psi_low + dpsi

                do icostheta = 1, ncostheta
                    costheta_low = (icostheta-1) * dcostheta - 1.
                    costheta_max = costheta_low + dcostheta
                
                    do ir = 1, nr
                        r_low = (ir - 1) * dr
                        r_max = r_low + dr

                        if(       costheta >= costheta_low .and. costheta < costheta_max &
                            .and. normOM   >= r_low        .and. normOM   < r_max &
                            .and. psi      >= psi_low      .and. psi      < psi_max ) then

                            g(ir,icostheta,ipsi) = g(ir,icostheta,ipsi) + solvent(1)%xi(io,ix,iy,iz)**2 * grid%w(io)
                            bin(ir,icostheta,ipsi) = bin(ir,icostheta,ipsi) + 1
                        end if
                    end do

                end do ! cosTheta'

            end do ! psi'

        end do ! xyz
    end do ! io

    where (bin /= 0) g = g / bin

    open( 77, file="output/g-of-r-costheta.dat")
    do icostheta = 1, ncostheta
        costheta_low = (icostheta-1) * dcostheta - 1.
        costheta_max = costheta_low + dcostheta
        do ir = 1, nr
            r_low = (ir - 1) * dr
            r_max = r_low + dr
            write(77,*) (r_low+r_max)/2., ( costheta_low + costheta_max )/2., sum(g(ir,icostheta,:)/size(g,3))
        end do
        write(77,*)
    end do
    close(77)
! gnuplot> set pm3d; set border 4095; set xlabel "r (\305)"; set ylabel "cos {/Symbol q}"; set zlabel "g"; unset key; splot "output/g-of-r-costheta.dat" w l

    open( 77, file="output/g-of-costheta-psi-at-3ang.dat")
    ir = int( 3.1/dr) +1
    block
        integer :: toto(3)
        toto = maxloc(g)
        ir = toto(1)
    end block
    do ipsi = 1, npsi
        psi_low = (ipsi-1) * dpsi
        psi_max = psi_low + dpsi
        do icostheta = 1, ncostheta
            costheta_low = (icostheta-1) * dcostheta - 1.
            costheta_max = costheta_low + dcostheta
            write(77,*) (psi_low + psi_max)/2._dp, ( costheta_low + costheta_max )/2._dp, g(ir,icostheta,ipsi)
        end do
        write(77,*)
    end do
    close(77)
! gnuplot> set pm3d; set border 4095; set xlabel "{/Symbol Y}"; set ylabel "cos {/Symbol q}"; set zlabel "g"; unset key; splot "output/g-of-costheta-psi-at-3ang.dat" w l

    stop "OOOOOOOOOOOOYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYY"

end subroutine output_gOfRAndCosThetaAndPsi
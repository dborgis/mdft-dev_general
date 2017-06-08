subroutine output_g-r-costheta-psi

    use module_grid, only: grid
    use precision_kinds, only: dp
    use module_solute, only: solute
    use module_solvent, only: solvent
    use mathematica,     only: deduce_optimal_histogram_properties
   
    implicit none

    real(dp), parameter :: zerodp = 0.0_dp, onedp = 1.0_dp
    real(dp), parameter :: pi = acos(-1.0_dp)
    integer :: ix, iy ,iz, io, ncostheta, icostheta, ir, nr
    double precision :: cosTheta, OM(3), normOM, H1(3), H2(3), Hbar(3), normHbar
    double precision, parameter :: Oz(3) = [ 0., 0., 1. ]
    double precision :: costheta_low, costheta_max, dcostheta, dr, r_low, r_max
    double precision, allocatable :: g(:,:), bin(:,:)


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

    ! The maximum of the g(r,\theta) is for r = rmaxg et theta = thetamaxg
block
    integer :: irmaxg(2), rmaxg, costhetamaxg
    irmaxg(1:2) = maxloc( g )
    print*, "irmaxg = ", irmaxg
    rmaxg = (irmaxg(1) - 0.5_dp) * dr
    costhetamaxg = (irmaxg(2)  * dcostheta
    print*, "rmaxg = ", rmaxg
    print*, "costhetamaxg = ", costhetamaxg
!    do iz = 1, grid%nz
!        z =  (iz - 1) * grid%dz
!        zsup = iz * grid%dz
!        if( z - rmaxg <= 0 .and. zsup - rmaxg > 0 ) then ! on est dans le bon interval
!            do io = 1, grid%no
!                if( grid%theta(io) - )
!            end do
!        end if
!    end do
end block

    stop "VICTORYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYY"


  ! integer  :: n,i,j,k,o,p,s,ntheta,npsi,nx,ny,nz,binthetamax,bintheta,binrmax,binr
  ! real(dp) :: dx, dy, dz, r(3), normr, costheta, theta, dr, dtheta, maxrange
  ! real(dp), allocatable :: rho(:,:,:,:)
  ! real(dp), parameter :: epsdp=epsilon(1._dp)
  ! real(dp), allocatable :: g(:,:), gcount(:,:) ! g(r,theta)
  !
  ! nx=grid%nx
  ! ny=grid%ny
  ! nz=grid%nz
  ! ntheta=grid%ntheta
  ! allocate( rho(nx,ny,nz,ntheta) )
  ! dx=grid%dx
  ! dy=grid%dy
  ! dz=grid%dz
  !
  ! ! Get the density per spherical angle Omega (integrate over psi)
  ! npsi = grid%npsi
  ! do s=1,size(solvent)
  !   do i=1,nx
  !     do j=1,ny
  !       do k=1,nz
  !         do o=1,ntheta
  !           rho(i,j,k,o) = solvent(s)%rho0 * sum( molRotGrid%weight * cg_vect_new(i,j,k,o,:,s)**2 )
  !         end do
  !       end do
  !     end do
  !   end do
  ! end do
  !
  ! ! Test that all vectors Omega are normalized (i.e., it is a unitary vector)
  ! do o=1,ntheta
  !   if( abs( norm2([OMx(o),OMy(o),OMz(o)])-1 ) > epsdp ) then
  !     print*,"ERROR: l.41 of compute_gOfRandCosTheta.f90"
  !     print*,"-----  Vector Omega is not normalized for omega index",int(o,1)
  !     stop
  !   end if
  ! end do
  !
  ! ! for the distance r between solute site and grid node:
  ! maxrange   = minval(grid%length)/2._dp
  ! call deduce_optimal_histogram_properties( nx*ny*nz, maxrange, binrmax, dr)
  ! binthetamax = 10 ! this has to be decided, but we have a very large number of angles, even if the quadrature over omega is by 1 angle only in fixed frame
  ! ! theta is between 0 and pi
  ! dtheta = pi/real(binthetamax,dp)
  ! allocate( g(binrmax,binthetamax) )
  ! allocate( gcount(binrmax,binthetamax)  )
  !
  ! open(124,file="./output/g-of-r-theta.out")
  ! do n=1,size(solute%site)
  !   g = zerodp
  !   gcount = zerodp
  !   do i=1,nx
  !     do j=1,ny
  !       do k=1,nz
  !         r = ([i,j,k]-1)*[dx,dy,dz] - solute%site(n)%r  ! vector that joins the grid node with index [i,j,k] to solute site n
  !         normr = norm2(r)
  !         if( abs(normr) <= epsdp ) cycle
  !         binr = int(normr/dr)+1
  !         if( binr > binrmax ) cycle ! don't try to bin something outside the range
  !         if( binr <= 0 ) then
  !           print*,"ERROR: l.73 of compute_gofrandcostheta.f90 binr has forbidden value",binr
  !           print*,"-----  binr = int(normr/dr)+1 with normr, dr, int(normr/dr)+1 =",normr,dr,int(normr/dr)+1
  !           stop
  !         end if
  !         do o=1,ntheta
  !           costheta = dot_product( r, [OMx(o), OMy(o), OMz(o)] )/normr
  !           if( costheta < -1 ) then
  !             costheta = -1 ! sometimes costheta=-1.000000000000001 (...round off errors)
  !           else if( costheta > 1 ) then
  !             costheta = 1
  !           end if
  !           theta = acos(costheta)
  !           bintheta = int(theta/(dtheta+2*epsdp)) +1! if dtheta is exactly theta/binthetamax, we end up with int(theta/dtheta)+1 = binthetamax+1. Numerical pb.
  !           if( bintheta <= 0 .or. bintheta > binthetamax ) then
  !             print*,"ERROR: l.82 of compute_gofrandcostheta.f90: bintheta has forbidden value",bintheta
  !             print*,"-----  bintheta = int(theta/dtheta) +1, with &
  !             & costheta, theta, dtheta, int(theta/dtheta)+1 =", costheta, theta, dtheta, int(theta/dtheta)+1
  !             stop
  !           end if
  !           g(binr,bintheta) = g(binr,bintheta) + rho(i,j,k,o)*anggrid%weight(o)
  !           gcount(binr,bintheta) = gcount(binr,bintheta) +1
  !         end do
  !       end do
  !     end do
  !   end do
  !
  !   ! test anormalities in g
  !   if( any(g/=g) ) then
  !     print*,"STOP: Nan or Infty somewhere in the histogram of g(r,theta), l.96 of compute_gofrtheta.f90"
  !     STOP
  !   end if
  !
  !   ! normalize
  !   if( size(solvent)/=1 ) stop "l.105 of compute_g-r-theta.f90. Implemented for only 1 solvent species"
  !   where( gcount /= 0 )
  !     g = g/gcount *anggrid%n_angles/(binthetamax*solvent(1)%n0) ! The last two terms here I can convince myself of but I have doubts
  !   else where
  !     g = zerodp
  !   end where
  !
  !   ! write g(r,theta) to an output file
  !   write(124,*)"#solute site",n
  !   do binr=1,binrmax
  !     do bintheta=1,binthetamax
  !       if( abs(g(binr,bintheta)) > epsdp ) then ! write to file only non-zero values. Consequences: File 10-20% smaller, but in-homogeneous printed grid
  !         write(124,*) (binr-0.5)*dr, (bintheta-0.5)*dtheta, g(binr,bintheta)
  !       end if
  !     end do
  !   end do
  !   write(124,*) ! empty line
  !
  ! end do ! loop over solute sites
  ! close(124)
  ! deallocate( g, gcount )
end subroutine output_g-r-costheta-psi

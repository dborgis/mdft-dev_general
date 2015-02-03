subroutine output_gOfRandCosTheta
  use precision_kinds, only: dp
  use system,          only: spacegrid, solute, solvent
  use quadrature,      only: anggrid, molrotgrid, molrotsymorder, OMx, OMy, OMz
  use constants,       only: zerodp, pi
  use minimizer,       only: cg_vect
  use mathematica,     only: deduce_optimal_histogram_properties

  implicit none
  integer  :: icg,n,i,j,k,o,p,s,omax,pmax,imax,jmax,kmax,binthetamax,bintheta,binrmax,binr
  real(dp) :: dx, dy, dz, r(3), normr, costheta, theta, dr, dtheta, maxrange
  real(dp), allocatable :: rho(:,:,:,:)
  real(dp), parameter :: epsdp=epsilon(1._dp)
  real(dp), allocatable :: g(:,:), gcount(:,:) ! g(r,theta)

  imax=spacegrid%n_nodes(1)
  jmax=spacegrid%n_nodes(2)
  kmax=spacegrid%n_nodes(3)
  omax=anggrid%n_angles
  allocate( rho(imax,jmax,kmax,omax) )
  dx=spacegrid%dl(1)
  dy=spacegrid%dl(2)
  dz=spacegrid%dl(3)

  ! Get the density per spherical angle Omega (integrate over psi)
  pmax=molrotgrid%n_angles
  icg=0
  do s=1,size(solvent)
    do i=1,imax
      do j=1,jmax
        do k=1,kmax
          do o=1,omax
            rho(i,j,k,o) = solvent(s)%rho0 * sum( molRotGrid%weight * cg_vect(icg+1:icg+pmax)**2 )
            icg = icg +pmax
          end do
        end do
      end do
    end do
  end do

  ! Test that all vectors Omega are normalized (i.e., it is a unitary vector)
  do o=1,omax
    if( abs( norm2([OMx(o),OMy(o),OMz(o)])-1 ) > epsdp ) then
      print*,"ERROR: l.41 of compute_gOfRandCosTheta.f90"
      print*,"-----  Vector Omega is not normalized for omega index",int(o,1)
      stop
    end if
  end do

  ! for the distance r between solute site and grid node:
  maxrange   = minval(spaceGrid%length)/2._dp
  call deduce_optimal_histogram_properties( imax*jmax*kmax, maxrange, binrmax, dr)
  binthetamax = 10 ! this has to be decided, but we have a very large number of angles, even if the quadrature over omega is by 1 angle only in fixed frame
  ! theta is between 0 and pi
  dtheta = pi/real(binthetamax,dp)
  allocate( g(binrmax,binthetamax) )
  allocate( gcount(binrmax,binthetamax)  )

  open(124,file="./output/g-of-r-theta.out")
  do n=1,size(solute%site)
    g = zerodp
    gcount = zerodp
    do i=1,imax
      do j=1,jmax
        do k=1,kmax
          r = ([i,j,k]-1)*[dx,dy,dz] - solute%site(n)%r  ! vector that joins the grid node with index [i,j,k] to solute site n
          normr = norm2(r)
          if( abs(normr) <= epsdp ) cycle
          binr = int(normr/dr)+1
          if( binr > binrmax ) cycle ! don't try to bin something outside the range
          if( binr <= 0 ) then
            print*,"ERROR: l.73 of compute_gofrandcostheta.f90 binr has forbidden value",binr
            print*,"-----  binr = int(normr/dr)+1 with normr, dr, int(normr/dr)+1 =",normr,dr,int(normr/dr)+1
            stop
          end if
          do o=1,omax
            costheta = dot_product( r, [OMx(o), OMy(o), OMz(o)] )/normr
            if( costheta < -1 ) then
              costheta = -1 ! sometimes costheta=-1.000000000000001 (...round off errors)
            else if( costheta > 1 ) then
              costheta = 1
            end if
            theta = acos(costheta)
            bintheta = int(theta/(dtheta+2*epsdp)) +1! if dtheta is exactly theta/binthetamax, we end up with int(theta/dtheta)+1 = binthetamax+1. Numerical pb.
            if( bintheta <= 0 .or. bintheta > binthetamax ) then
              print*,"ERROR: l.82 of compute_gofrandcostheta.f90: bintheta has forbidden value",bintheta
              print*,"-----  bintheta = int(theta/dtheta) +1, with &
              & costheta, theta, dtheta, int(theta/dtheta)+1 =", costheta, theta, dtheta, int(theta/dtheta)+1
              stop
            end if
            g(binr,bintheta) = g(binr,bintheta) + rho(i,j,k,o)*anggrid%weight(o)
            gcount(binr,bintheta) = gcount(binr,bintheta) +1
          end do
        end do
      end do
    end do

    ! test anormalities in g
    if( any(g/=g) ) then
      print*,"STOP: Nan or Infty somewhere in the histogram of g(r,theta), l.96 of compute_gofrtheta.f90"
      STOP
    end if

    ! normalize
    if( size(solvent)/=1 ) stop "l.105 of compute_g-r-theta.f90. Implemented for only 1 solvent species"
    where( gcount /= 0 )
      g = g/gcount *anggrid%n_angles/(binthetamax*solvent(1)%n0) ! The last two terms here I can convince myself of but I have doubts
    else where
      g = zerodp
    end where

    ! write g(r,theta) to an output file
    write(124,*)"#solute site",n
    do binr=1,binrmax
      do bintheta=1,binthetamax
        if( abs(g(binr,bintheta)) > epsdp ) then ! write to file only non-zero values. Consequences: File 10-20% smaller, but in-homogeneous printed grid
          write(124,*) (binr-0.5)*dr, (bintheta-0.5)*dtheta, g(binr,bintheta)
        end if
      end do
    end do
    write(124,*) ! empty line

  end do ! loop over solute sites
  close(124)
  deallocate( g, gcount )
end subroutine output_gOfRandCosTheta

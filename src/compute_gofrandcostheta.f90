subroutine compute_gOfRandCosTheta
  use precision_kinds, only: dp
  use system,          only: spacegrid, solute
  use quadrature,      only: anggrid, molrotgrid, molrotsymorder, OMx, OMy, OMz
  use constants,       only: zerodp, pi
  use minimizer,       only: cg_vect
  use mathematica,     only: deduce_optimal_histogram_properties

  implicit none
  integer  :: icg,n,i,j,k,o,p,omax,pmax,imax,jmax,kmax,binthetamax,bintheta,binrmax,binr
  real(dp) :: dx, dy, dz, r(3), normr, costheta, costhetamin, costhetamax, theta, dr, dtheta, maxrange
  real(dp), allocatable :: rho(:,:,:,:)
  real(dp), parameter :: epsdp=epsilon(1._dp)
  real(dp), allocatable :: g(:,:), gcount(:,:) ! g(r,theta)

  imax=spacegrid%n_nodes(1)
  jmax=spacegrid%n_nodes(2)
  kmax=spacegrid%n_nodes(3)
  omax=anggrid%n_angles
  allocate( rho(imax,jmax,kmax,omax) )
  pmax=molrotgrid%n_angles
  dx=spacegrid%dl(1)
  dy=spacegrid%dl(2)
  dz=spacegrid%dl(3)

  ! Get the density per spherical angle Omega (integrate over psi)
  icg=0
  do i=1,imax
    do j=1,jmax
      do k=1,kmax
        do o=1,omax
          rho(i,j,k,o) = zerodp ! init
          do p=1,pmax
            icg = icg +1
            rho(i,j,k,o) = rho(i,j,k,o) + molRotGrid%weight(p)*cg_vect(icg)**2/(8*pi**2/molrotsymorder) ! integrate over psi, which is the rotation around vector omega
          end do
        end do
      end do
    end do
  end do

  ! for each grid node [i,j,k], compute the distance r
  print*,"WARNING: l.35 of compute_gOfRandCosTheta.f90. Works only for 1 solute site"
  print*,"-------"

  ! Test that all vectors Omega are normalized
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
  allocate( g(binrmax,binthetamax) ,source=zerodp )
  allocate( gcount(binrmax,binthetamax) ,source=zerodp )

  ! Compute the
  open(124,file="./output/gro.out")
  costhetamin=  2
  costhetamax= -2
  do i=1,imax
    do j=1,jmax
      do k=1,kmax
        if( all([i,j,k]==1) ) cycle
        r = ([i,j,k]-1)*[dx,dy,dz] - solute%site(1)%r
        normr = norm2(r) ! vector that joins the grid node with index [i,j,k] to solute site [0,0,0] ! TODO
        if( normr > maxrange ) cycle
        binr = int(normr/dr)+1; if(binr==0) STOP "hoooooooooooooooooalllkjsdlkjazd"
        if( abs(normr) <= epsdp ) cycle
        do o=1,omax
          costheta = dot_product( r, [OMx(o), OMy(o), OMz(o)] )/normr
          if( costheta < costhetamin ) costhetamin=costheta
          if( costheta > costhetamax ) costhetamax=costheta
          theta = acos(costheta)
          bintheta = int(theta/dtheta) +1
          g(binr,bintheta) = g(binr,bintheta) + rho(i,j,k,o)
          gcount(binr,bintheta) = gcount(binr,bintheta) +1
          ! write(124,*) real([normr, theta, rho(i,j,k,o)], 4)
        end do
      end do
    end do
  end do
  do concurrent( binr=1:binrmax, bintheta=1:binthetamax )
    if( gcount(binr,bintheta) /= 0 ) then
      g(binr,bintheta) = g(binr,bintheta) / gcount(binr,bintheta)
    else
      g(binr,bintheta) = zerodp
    end if
  end do

  do binr=1,binrmax
    do bintheta=1,binthetamax
      write(124,*) (binr-0.5)*dr, (bintheta-0.5)*dtheta, real(g(binr,bintheta),4)
    end do
  end do
  close(124)
  print*,"cos theta min and max:", costhetamin, costhetamax
  stop

end subroutine compute_gOfRandCosTheta

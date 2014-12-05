subroutine output_gsitesite
  use precision_kinds, only: dp
  use system,          only: spacegrid, solute, solvent
  use quadrature,      only: anggrid, molrotgrid, molrotsymorder, OMx, OMy, OMz,&
                              Rotxx, Rotxy, Rotxz, Rotyx, Rotyy, Rotyz, Rotzx, Rotzy, Rotzz
  use constants,       only: zerodp, pi
  use minimizer,       only: cg_vect
  use mathematica,     only: deduce_optimal_histogram_properties

  implicit none
  integer  :: icg,n,i,j,k,o,p,s,omax,pmax,imax,jmax,kmax,binthetamax,bintheta,binrmax,binr,allocatestatus,ss
  real(dp) :: dx, dy, dz, r(3), normr, costheta, theta, dr, dtheta, rmax, xss, yss, zss
  real(dp), allocatable :: rho(:,:,:,:,:)
  real(dp), parameter :: epsdp=epsilon(1._dp)
  real(dp), allocatable :: g(:), gcount(:) ! g(r,theta)

  imax=spacegrid%n_nodes(1)
  jmax=spacegrid%n_nodes(2)
  kmax=spacegrid%n_nodes(3)
  omax=anggrid%n_angles
  pmax=molrotgrid%n_angles
  allocate( rho(imax,jmax,kmax,omax,pmax) )
  dx=spacegrid%dl(1)
  dy=spacegrid%dl(2)
  dz=spacegrid%dl(3)

  ! Get the density per spherical angle Omega (integrate over psi)
  icg=0
  do s=1,size(solvent)
    do i=1,imax
      do j=1,jmax
        do k=1,kmax
          do o=1,omax
            do p=1,pmax
              icg = icg+1
              rho(i,j,k,o,p) = cg_vect(icg)**2 ! rho(r,o,p)/rho0
            end do
          end do
        end do
      end do
    end do
  end do

  rmax = minval(spacegrid%length)/2._dp ! maximum range in angstroms of the radial distribution function g(r)
  call deduce_optimal_histogram_properties( imax*jmax*kmax, rmax, binrmax, dr) ! deduce optimal number of bins in r (binrmax) and, thus, dr, the binwidth

  allocate( g(binrmax) ,stat=allocatestatus) ! allocate the histogram g(r) accordingly
  if(allocatestatus /= 0) stop "ERROR: l.46 of output-gsitesite.f90. Allocation of g failed."
  allocate( gcount(binrmax) ,stat=allocatestatus)
  if(allocatestatus /= 0) stop "ERROR: l.48 of output-gsitesite.f90. Allocation of gcount failed."

  open(124,file="./output/g-sitesite.out")

  ! ! for each solute site, and each solvent site
  do n=1,size(solute%site)
    do s=1,size(solvent) ! loop over solvents
      do ss=1,size(solvent(s)%site) ! loop over all sites of the selected solvent
        write(124,*)"# g(r) between solute site",n," and solvent site",ss,"of solvent number",s
        g = zerodp
        gcount = zerodp
        do i=1,imax
          do j=1,jmax
            do k=1,kmax
              do o=1,omax
                do p=1,pmax
                  xss = dot_product( [Rotxx(o,p),Rotxy(o,p),Rotxz(o,p)] , solvent(s)%site(ss)%r )
                  yss = dot_product( [Rotyx(o,p),Rotyy(o,p),Rotyz(o,p)] , solvent(s)%site(ss)%r )
                  zss = dot_product( [Rotzx(o,p),Rotzy(o,p),Rotzz(o,p)] , solvent(s)%site(ss)%r )
                  r = ([i,j,k]-1)*[dx,dy,dz] +[xss,yss,zss] - solute%site(n)%r  ! vector that joins the grid node with index [i,j,k] to solute site n
                  normr = norm2(r)
                  if( abs(normr) <= epsdp ) cycle ! vector r is ill defined
                  binr = int(normr/dr) +1
                  if( binr > binrmax) cycle
                  if( binr <= 0) then
                    print*,"ERROR: l.79 of output_gsitesite.f90. binr is forbidden value:", binr
                    print*,"-----  binr = int(normr/dr)+1 with normr, dr, int(normr/dr)+1 =",normr,dr,int(normr/dr)+1
                    stop
                  end if
                  g(binr) = g(binr) +rho(i,j,k,o,p)*anggrid%weight(o)*molrotgrid%weight(p) !TODO shouldnt I use weights over psi and omega here?
                  gcount(binr) = gcount(binr) +1
                end do
              end do
            end do
          end do
        end do
        where( gcount /= 0 ) g = g/gcount /sum(anggrid%weight) /sum(molrotgrid%weight)
        if( any(g/=g) ) stop "STOP: Nan or Infty somewhere in the histogram of g(r), l.93 of compute_gofrtheta.f90"
        do binr=1,binrmax
          if( abs(g(binr)) > epsdp ) write(124,*) (binr-0.5)*dr, g(binr)
        end do
        write(124,*)
      end do
    end do
  end do
  close(124)
  deallocate( g, gcount )
end subroutine output_gsitesite

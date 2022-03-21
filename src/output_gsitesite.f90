subroutine output_gsitesite
     use precision_kinds, only: dp
     use module_solute, only: solute
     use module_solvent, only: solvent
     use module_grid, only: grid
     use mathematica, only: deduce_optimal_histogram_properties
    
     implicit none
     integer  :: icg,n,i,j,k,io,p,s,omax,pmax,imax,jmax,kmax,binthetamax,bintheta,binrmax,binr,allocatestatus,ss
     real(dp) :: dx, dy, dz, r(3), normr, costheta, theta, dr, dtheta, rmax, xss, yss, zss

     real(dp), parameter :: epsdp=epsilon(1._dp), zerodp = 0._dp, one =1._dp, pi = acos(-1._dp)
     real(dp), allocatable :: g(:), gcount(:) ! g(r,theta)
     integer :: nx, ny, nz, no, ns
    
     write(*,*) 'On ecrit les gr site-site'
     write(*,*) 'sumw =',sum(grid%w),8*pi**2/solvent(1)%molrotsymorder
     nx = grid%nx
     ny = grid%ny
     nz = grid%nz
     no = grid%no
     ns = solvent(1)%nspec

    
     imax=grid%n_nodes(1)
     jmax=grid%n_nodes(2)
     kmax=grid%n_nodes(3)

     dx=grid%dl(1)
     dy=grid%dl(2)
     dz=grid%dl(3)
    
         
     rmax = minval(grid%length)/2._dp ! maximum range in angstroms of the radial distribution function g(r)
     call deduce_optimal_histogram_properties( imax*jmax*kmax, rmax, binrmax, dr) ! deduce optimal number of bins in r (binrmax) and, thus, dr, the binwidth
    
     allocate( g(binrmax) ,stat=allocatestatus) ! allocate the histogram g(r) accordingly
     if(allocatestatus /= 0) stop "ERROR: l.46 of output-gsitesite.f90. Allocation of g failed."
     allocate( gcount(binrmax) ,stat=allocatestatus)
     if(allocatestatus /= 0) stop "ERROR: l.48 of output-gsitesite.f90. Allocation of gcount failed."
    
     open(124,file="output/g-sitesite.out")
    
     ! ! for each solute site, and each solvent site
     do n=1,size(solute%site)
       do s=1,size(solvent) ! loop over solvents
         do ss=1,size(solvent(s)%site) ! loop over all sites of the selected solvent
           write(124,*)"# g(r) between solute site",n," and solvent site",ss,"of solvent number",s
           g = zerodp
           gcount = zerodp
           do k=1,kmax
             do j=1,jmax
               do i=1,imax
                 do io=1,no
                  
                     xss = dot_product( [grid%rotxx(io),grid%rotxy(io),grid%rotxz(io)] , solvent(s)%site(ss)%r )
                     yss = dot_product( [grid%rotyx(io),grid%rotyy(io),grid%rotyz(io)] , solvent(s)%site(ss)%r )
                     zss = dot_product( [grid%rotzx(io),grid%rotzy(io),grid%rotzz(io)] , solvent(s)%site(ss)%r )
                     r = ([i,j,k]-1)*[dx,dy,dz] +[xss,yss,zss] - solute%site(n)%r  ! vector that joins the grid node with index [i,j,k] to solute site n
                     normr = norm2(r)
                     if( abs(normr) <= epsdp ) cycle ! vector r is ill defined
                     binr = int(normr/dr) + 1
                     if( binr > binrmax) cycle
                     if( binr <= 0) then
                       print*,"ERROR: l.79 of output_gsitesite.f90. binr is forbidden value:", binr
                       print*,"-----  binr = int(normr/dr)+1 with normr, dr, int(normr/dr)+1 =",normr,dr,int(normr/dr)+1
                       stop
                     end if
                     g(binr) = g(binr) + grid%w(io)*solvent(s)%xi(io,i,j,k)**2*solvent(s)%rho0/solvent(s)%n0
                     gcount(binr) = gcount(binr) + one
                     
                   end do  
                 end do
               end do
           end do
           where( gcount /= 0 ) g = g/gcount
           if( any(g/=g) ) stop "STOP: Nan or Infty somewhere in the histogram of g(r), l.93 of compute_gofrtheta.f90"
           write(124,*) zerodp, zerodp
           do binr=1,binrmax
             if( abs(g(binr)) > epsdp ) write(124,*) (binr-0.5)*dr, g(binr)/g(binrmax)  ! TROUVER LA BONNE NORMALISATION !!!!!
           end do
           write(124,*)
         end do
       end do
     end do
     close(124)
     deallocate( g, gcount )
end subroutine output_gsitesite
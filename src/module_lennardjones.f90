module module_lennardjones
    implicit none
    private
    public :: calcul_lennardjones
contains

    subroutine calcul_lennardjones
        use precision_kinds, only: dp
        use module_solute, only: solute
        use module_solvent, only: solvent
        use module_grid, only: grid
        use module_input, only: getinput
        implicit none
        real(dp) :: cutoff, cutoffsq
        integer :: nx, ny, nz, no, ns
        real(dp) :: lx, ly, lz
        real(dp), parameter :: fourpi=4._dp*acos(-1._dp)
        real(dp), parameter :: epsdp = epsilon(1._dp)
        integer :: i,j,k,s,u, io, ix, iy, iz, ss
        real(dp), allocatable :: x(:), y(:), z(:)
        real(dp) :: siguv6, epsuv
        real(dp) :: dx, dy, dz, xss, yss, zss
        real(dp), allocatable :: xmod(:,:,:), ymod(:,:,:), zmod(:,:,:)
        integer, allocatable :: xtab(:), ytab(:), ztab(:)
        integer :: minx, maxx, miny, maxy, minz, maxz
        integer :: xtabsize, ytabsize, ztabsize
        integer :: indextabx, indextaby, indextabz
        real(dp) :: xgrid, ygrid, zgrid

        if (.not. allocated(solvent)) error stop "solvent should be allocated in vext_lennardjones"
        if (.not. grid%isinitiated) error stop "grid is not initiated in vext_lennardjones"



        cutoff = getinput%dp('rvdw', defaultvalue=10.0_dp, assert=">=0")
        cutoffsq = cutoff**2
        nx = grid%nx
        ny = grid%ny
        nz = grid%nz
        lx = grid%lx
        ly = grid%ly
        lz = grid%lz
        no = grid%no
        ns = size(solvent)

        xtabsize = floor(2*cutoff/grid%dx)+1
        ytabsize = floor(2*cutoff/grid%dy)+1
        ztabsize = floor(2* cutoff/grid%dz)+1

        allocate( x(nx) ,source= [(real(i-1,dp)*grid%dx ,i=1,nx)] )
        allocate( y(ny) ,source= [(real(j-1,dp)*grid%dy ,j=1,ny)] )
        allocate( z(nz) ,source= [(real(k-1,dp)*grid%dz, k=1,nz)] )
        allocate( xtab(xtabsize) )
        allocate( ytab(ytabsize) )
        allocate( ztab(ztabsize) )


        block
          integer :: ssmax
          ssmax = maxval(  [ (size(solvent(i)%site) ,i=1,ns)]    )
          allocate( xmod(no,ssmax,ns))
          allocate( ymod(no,ssmax,ns))
          allocate( zmod(no,ssmax,ns))
        end block


        do s=1,ns
          do ss=1,solvent(s)%nsite
            do io=1,no
              xmod (io,ss,s) = DOT_PRODUCT( [grid%rotxx(io), grid%rotxy(io), grid%rotxz(io)] , solvent(s)%site(ss)%r )
              ymod (io,ss,s) = DOT_PRODUCT( [grid%rotyx(io), grid%rotyy(io), grid%rotyz(io)] , solvent(s)%site(ss)%r )
              zmod (io,ss,s) = DOT_PRODUCT( [grid%rotzx(io), grid%rotzy(io), grid%rotzz(io)] , solvent(s)%site(ss)%r )
            end do
          end do
        end do





        do u=1,size(solute%site)
          if( solute%site(u)%eps <= epsdp ) cycle ! if the solute site does not wear a Lennard-Jones contribution

          !prepare table to loop over
          minx = floor(mod((solute%site(u)%r(1)-cutoff+lx), lx)/grid%dx)
          maxx = floor(mod((solute%site(u)%r(1)+cutoff   ), lx)/grid%dx)
          miny = floor(mod((solute%site(u)%r(2)-cutoff+ly), ly)/grid%dy)
          maxy = floor(mod((solute%site(u)%r(2)+cutoff   ), ly)/grid%dy)
          minz = floor(mod((solute%site(u)%r(3)-cutoff+lz), lz)/grid%dz)
          maxz = floor(mod((solute%site(u)%r(3)+cutoff   ), lz)/grid%dz)

          if (minx<maxx) then
            xtab = (/ (I, I = minx, maxx) /)
          else
            xtab = [(/ (I, I = 0, maxx) /), (/ (I, I = minx, nx) /)]
          end if

          if (miny<maxy) then
            ytab = (/ (I, I = miny, maxy) /)
          else
            ytab = [(/ (I, I = 0, maxy) /), (/ (I, I = miny, ny) /)]
          end if

          if (minz<maxz) then
            ztab = (/ (I, I = minz, maxz) /)
          else
            ztab = [(/ (I, I = 0, maxz) /), (/ (I, I = minz, nz) /)]
          end if



          !end of prepare table to loop over


          do s=1,ns ! Ces boucles peuvent etre remise à l'interieur si on créé un tableau au début qui contient les indices des sites sur lesquels on peut boucler
            do ss=1,size(solvent(s)%site)
              if( solvent(s)%site(ss)%eps<=epsdp ) cycle
              epsuv=sqrt(solute%site(u)%eps * solvent(s)%site(ss)%eps)
              siguv6=(  (solute%site(u)%sig + solvent(s)%site(ss)%sig)/2._dp)**6

              do indextabz=1,ztabsize
                iz = ztab(indextabz)
                zgrid=z(iz)

                do indextaby=1,ytabsize
                  iy = ytab(indextaby)
                  ygrid=y(iy)

                  do indextabx=1,xtabsize
                    ix = xtab(indextabx)
                    xgrid=x(ix)

                    do io=1,no !! sortir io des xyz !!! a voir avec les acces memoire vext!!!
                      if( solvent(s)%vext(io,ix,iy,iz) > 1.e5 ) cycle

                      xss=xgrid+xmod(io,ss,s)
                      yss=ygrid+ymod(io,ss,s)
                      zss=zgrid+zmod(io,ss,s)

                      dx =abs(xss-solute%site(u)%r(1)); do while(dx>lx/2._dp); dx=abs(dx-lx); end do
                      dy =abs(yss-solute%site(u)%r(2)); do while(dy>ly/2._dp); dy=abs(dy-ly); end do
                      dz =abs(zss-solute%site(u)%r(3)); do while(dz>lz/2._dp); dz=abs(dz-lz); end do

                      solvent(s)%vext(io,ix,iy,iz) = solvent(s)%vext(io,ix,iy,iz) + vlj( epsuv, siguv6, dx**2+dy**2+dz**2, cutoffsq)
                    end do
                  end do
                end do
              end do
            end do
          end do
        end do

    end subroutine calcul_lennardjones




    pure function vlj(eps,sig6,rsq, cutoffsq) ! I hope this function is inlined by the compiler
        use precision_kinds, only: dp
        implicit none
        real(dp) :: vlj
        real(dp), intent(in) :: cutoffsq
        real(dp), intent(in) :: eps, sig6, rsq
        real(dp) :: div
        real(dp), parameter :: epsdp=1._dp
        if (rsq<=epsdp) then
            vlj = huge(1._dp)
        elseif (rsq>cutoffsq) then
            vlj = 0._dp
        else
            div = sig6/rsq**3 ! rsq is a distance²
            vlj = 4._dp*eps*div*(div-1._dp)
        end if
    end function vlj

end module module_lennardjones




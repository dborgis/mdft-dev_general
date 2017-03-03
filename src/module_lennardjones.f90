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
        real(dp) :: div, rsq
        real(dp) :: vlj

        if (.not. allocated(solvent)) error stop "solvent should be allocated in vext_lennardjones"
        if (.not. grid%isinitiated) error stop "grid is not initiated in vext_lennardjones"


        !
        ! by default, the Lennard-Jones potential has a cutoff of 10 Angstroms.
        ! This is also the default in Gromacs 6.
        ! The input tag rvdw of dft.in can change this behavior.
        !
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

        if (cutoff .gt. lx*0.5_dp .or. cutoff .gt. ly*0.5_dp .or. cutoff .gt. lz*0.5_dp) then
            print*, "ERROR: the cut off for the LJ potential is larger than half the box size"
            print*, "lx, ly, lz=",lx,ly,lz,"and cutoff=",cutoff
            error stop "see dft.in and src/module_lennardjones"
        end if

        xtabsize = floor(2*cutoff/grid%dx)+1
        ytabsize = floor(2*cutoff/grid%dy)+1
        ztabsize = floor(2*cutoff/grid%dz)+1

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

          ! prepare table to loop over
          ! These tables allow to loop only on a cube with l=cutoff
          minx = floor(mod((solute%site(u)%r(1)-cutoff+lx), lx)/grid%dx) + 1
          maxx = floor(mod((solute%site(u)%r(1)+cutoff   ), lx)/grid%dx) + 1
          miny = floor(mod((solute%site(u)%r(2)-cutoff+ly), ly)/grid%dy) + 1
          maxy = floor(mod((solute%site(u)%r(2)+cutoff   ), ly)/grid%dy) + 1
          minz = floor(mod((solute%site(u)%r(3)-cutoff+lz), lz)/grid%dz) + 1
          maxz = floor(mod((solute%site(u)%r(3)+cutoff   ), lz)/grid%dz) + 1
          if(minx<1) error stop "minx<1 in module_lennardjones"
          if(miny<1) error stop "miny<1 in module_lennardjones"
          if(minz<1) error stop "minz<1 in module_lennardjones"
          if(maxx>nx) error stop "maxx>nx in module_lennardjones"
          if(maxy>ny) error stop "maxy>ny in module_lennardjones"
          if(maxz>nz) error stop "maxz>nz in module_lennardjones"

          if (minx<maxx) then
            xtab = (/ (I, I = minx, maxx) /)
          else ! take into account PBC
            xtab = [(/ (I, I = 1, maxx) /), (/ (I, I = minx, nx) /)]
          end if

          if (miny<maxy) then
            ytab = (/ (I, I = miny, maxy) /)
          else ! take into account PBC
            ytab = [(/ (I, I = 1, maxy) /), (/ (I, I = miny, ny) /)]
          end if

          if (minz<maxz) then
            ztab = (/ (I, I = minz, maxz) /)
          else ! take into account PBC
            ztab = [(/ (I, I = 1, maxz) /), (/ (I, I = minz, nz) /)]
          end if
          !end of prepare table to loop over


          do s=1,ns ! Ces boucles peuvent etre remise à l'interieur si on créé un tableau au début qui contient les indices des sites sur lesquels on peut boucler
            do ss=1,size(solvent(s)%site)
              if( solvent(s)%site(ss)%eps<=epsdp ) cycle
              epsuv=sqrt(solute%site(u)%eps * solvent(s)%site(ss)%eps)
              siguv6=(  (solute%site(u)%sig + solvent(s)%site(ss)%sig)*0.5_dp)**6

              !$omp parallel private(indextabz, iz, zgrid, indextaby, iy, ygrid, indextabx, ix, xgrid, io, xss, yss, zss, dx, dy, dz, rsq, vlj, div)
              !$omp do
              do indextabz=1,ztabsize ! indextabz is the index in ztabsize
                iz = ztab(indextabz)  ! iz is the index of the point in the grid
                zgrid=z(iz)           ! zgrid is the z position of the point

                do indextaby=1,ytabsize
                  iy = ytab(indextaby)
                  ygrid=y(iy)

                  do indextabx=1,xtabsize
                    ix = xtab(indextabx)
                    xgrid=x(ix)

                    do io=1,no !! sortir io des xyz !!! a voir avec les acces memoire vext!!!
                      if( solvent(s)%vext(io,ix,iy,iz) > 1.e5 ) cycle ! TODO reflechir a un critere plus malin

                      xss=xgrid+xmod(io,ss,s)
                      yss=ygrid+ymod(io,ss,s)
                      zss=zgrid+zmod(io,ss,s)

                      dx =abs(xss-solute%site(u)%r(1)); do while(dx>lx*0.5_dp); dx=abs(dx-lx); end do
                      dy =abs(yss-solute%site(u)%r(2)); do while(dy>ly*0.5_dp); dy=abs(dy-ly); end do
                      dz =abs(zss-solute%site(u)%r(3)); do while(dz>lz*0.5_dp); dz=abs(dz-lz); end do


                      rsq = dx**2+dy**2+dz**2

                      if (rsq<=epsdp) then
                        vlj = huge(1._dp)
                      elseif (rsq>cutoffsq) then
                        vlj = 0._dp
                      else
                        div = siguv6/rsq**3 ! rsq is a distance²
                        vlj = 4._dp*epsuv*div*(div-1._dp)
                      end if

                      solvent(s)%vext(io,ix,iy,iz) = solvent(s)%vext(io,ix,iy,iz) + vlj
                    end do
                  end do
                end do
             end do
             !$omp end do
             !$omp end parallel
           end do
          end do
        end do

        deallocate(xmod, ymod, zmod, x, y, z, xtab, ytab, ztab)

    end subroutine calcul_lennardjones


end module module_lennardjones

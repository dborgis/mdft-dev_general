module module_density

    use precision_kinds, only: dp
    implicit none
    private
    public :: init_density

contains

    !
    ! Guess an initial density for the solvent, or read it from a restart file.
    !
    subroutine init_density
        use module_solvent, only: solvent
        use module_grid, only: grid
        use module_input, only: getinput
        implicit none
        integer :: s, ios

        !
        ! Be sure the object containing all information about the solvent is initiated
        ! Also, it should be the same for the grid.
        !
        if (.not. allocated(solvent)) then
            print*, "Dans module_density, solvent n'est pas allocated"
            print*, "buggy. bisous"
            error stop
        end if
        if (.not. grid%isinitiated) then
            print*, "Dans module_density, grid is not allocated"
            error stop
        end if

        ! allocate the solvent density field for each solvent species. Remember : xi**2=rho/rho0
        do s=1,solvent(1)%nspec
            if (.not. allocated( solvent(s)%xi) ) then
                allocate(solvent(s)%xi(grid%no, grid%nx, grid%ny, grid%nz),   source=1._dp, stat=ios)
                if (ios /= 0) then
                    print*, "solvent(s)%xi(grid%no,grid%nx,grid%ny,grid%nz), source=0._dp: Allocation request denied"
                    print*, "for s =", s
                end if
            end if
        end do

        !
        ! Should we restart from a restart file or guess an initial density ?
        !
        select case( getinput%log('restart', defaultvalue=.false.))
        case(.true.)
          call read_restart_file
        case(.false.)
          call guess_density
        end select

    end subroutine init_density


    !
    ! Read the restart file
    !
    subroutine read_restart_file
      use module_solvent, only: solvent
      use module_solute, only: solute
      use module_grid, only: grid
      use system, only: site_type
      use module_orientation_projection_transform, only: proj2angl, init_module_rotation => init
      implicit none
      logical :: exists
      integer :: ios

      integer :: ix, iy, iz, ix_prev, iy_prev, iz_prev
      integer :: ns_prev, no_prev, nsite_prev, isite, is, no
      integer :: mmax_prev, np_prev, io, mmax, np
      integer :: nx_prev, ny_prev, nz_prev, nx, ny, nz
      type(site_type) :: tmp_site
      real(dp) :: dx_prev, dy_prev, dz_prev, dx, dy, dz
      real(dp) :: x, y, z
      real(dp) :: x0, y0, z0, x1, y1, z1, xd, yd, zd
      real(dp) :: xi00, xi01, xi10, xi11, xi0, xi1
      real(dp) :: offset_x, offset_y, offset_z
      real(dp) :: rho0
      complex(dp), allocatable :: deltarho_p_prev(:,:,:,:)
      real(dp), allocatable :: xi_prev(:,:,:,:)
      real(dp), parameter :: zero=0._dp
      complex(dp), parameter :: zeroc=(0._dp, 0._dp)
      logical :: read_full_density
      character(len("input/density.bin")), parameter :: filename="input/density.bin"
      inquire (file=filename, EXIST=exists)
      if (.not.exists) stop "You want to restart from a density file, but input/density.bin is not found"
      OPEN (10, file = filename , form = 'unformatted' , iostat=ios, status='OLD' )
      if ( ios /= 0 ) then
        print *, 'problem while opening input/density.bin. bug at init_density.f90'
        stop
      end if

      dx=grid%dx
      dy=grid%dy
      dz=grid%dz
      nx=grid%nx
      ny=grid%ny
      nz=grid%nz
      np=grid%np
      no=grid%no
      mmax=grid%mmax

      read ( 10, iostat=ios ) ns_prev
      if (ns_prev<0) then
        read_full_density=.true.
        print*, "The denisity.bin file you are rereading contains xi(Omega,r)"
        ns_prev=-ns_prev
      else
        read_full_density=.false.
        print*,  "The denisity.bin file you are rereading is stored in projections"
      end if
      read ( 10, iostat=ios ) mmax_prev
      read ( 10, iostat=ios ) no_prev
      read ( 10, iostat=ios ) np_prev
      read ( 10, iostat=ios ) nx_prev, ny_prev, nz_prev
      read ( 10, iostat=ios ) dx_prev, dy_prev, dz_prev
      read ( 10, iostat=ios ) nsite_prev

      print*
      print*, '*** RESTARTING from input/density.bin ***'
      print*
      print *,'previous state'
      print *,'nx=', nx_prev,'  ny=', ny_prev,'  nz=', nz_prev
      print *, "mmax=", mmax
      call init_module_rotation
      offset_x = (grid%lx - dx_prev * nx_prev) / 2
      offset_y = (grid%ly - dy_prev * ny_prev) / 2
      offset_z = (grid%lz - dz_prev * nz_prev) / 2

      if(nsite_prev /= size(solute%site) ) then
        print*, "Error in restart. Solutes are differents"
        stop
      endif

      do isite=1,nsite_prev
        read ( 10, iostat=ios ) tmp_site
        if( (solute%site(isite)%r(1) /= tmp_site%r(1)+offset_x) .or. &
            (solute%site(isite)%r(2) /= tmp_site%r(2)+offset_y) .or. &
            (solute%site(isite)%r(3) /= tmp_site%r(3)+offset_z) .or. &
            (solute%site(isite)%sig /= tmp_site%sig) .or. &
            (solute%site(isite)%eps /= tmp_site%eps) ) then
          print*, "Error in restart. Solutes are differents"
          stop
        endif
      enddo

      ! TODO: test si les ns sont les mêmes

      if (.not. read_full_density) then 
      if (.not. allocated (xi_prev) ) allocate (xi_prev(no,nx_prev,ny_prev,nz_prev) ,source=0._dp)
      if (.not. allocated (deltarho_p_prev) ) allocate (deltarho_p_prev(np,nx_prev,ny_prev,nz_prev) ,source=zeroc)
      !load previous projections
        do is=1,size(solvent)
          ! read file
          do iz_prev=1,nz_prev
            do iy_prev=1,ny_prev
              do ix_prev=1,nx_prev
                read ( 10, iostat=ios ) deltarho_p_prev(1:np_prev,ix_prev,iy_prev,iz_prev)
              end do
            end do
          end do
          
          !proj2angl
          do iz_prev=1,nz_prev
            do iy_prev=1,ny_prev
              do ix_prev=1,nx_prev
                  call proj2angl( &
                      deltarho_p_prev(1:np,ix_prev,iy_prev,iz_prev) ,&
                      xi_prev(1:no,ix_prev,iy_prev,iz_prev) )
              end do
            end do
          end do
          !interpolation
          do io=1,no

            do iz_prev=1,nz_prev-1
              z0 = (iz_prev-1)*dz_prev+offset_z
              z1 = (iz_prev)*dz_prev+offset_z
              do iy_prev=1,ny_prev-1
                y0 = (iy_prev-1)*dy_prev+offset_y
                y1 = (iy_prev)*dy_prev+offset_y
                do ix_prev=1,nx_prev-1
                  x0 = (ix_prev-1)*dx_prev+offset_x
                  x1 = (ix_prev)*dx_prev+offset_x

                  iz=max(1, ceiling(z0/dz+1))
                  do while(iz<=min(nz, floor(z1/dz+1)))
                    iy=max(1, ceiling(y0/dy+1))
                    do while(iy<=min(ny, floor(y1/dy+1)))
                      ix=max(1, ceiling(x0/dx+1))
                      do while(ix<=min(nx, floor(x1/dx+1)))

                        z = (iz-1)*dz
                        y = (iy-1)*dy
                        x = (ix-1)*dx

                        xd = (x-x0)/(x1-x0)
                        yd = (y-y0)/(y1-y0)
                        zd = (z-z0)/(z1-z0)

                        xi00 = xi_prev(io, ix_prev, iy_prev, iz_prev) * (1-xd) + xi_prev(io, ix_prev+1, iy_prev, iz_prev) * xd
                        xi01 = xi_prev(io, ix_prev, iy_prev, iz_prev+1) * (1-xd) + xi_prev(io, ix_prev+1, iy_prev, iz_prev+1) * xd
                        xi10 = xi_prev(io, ix_prev, iy_prev+1, iz_prev) * (1-xd) + xi_prev(io, ix_prev+1, iy_prev+1, iz_prev) * xd
                        xi11 = xi_prev(io, ix_prev, iy_prev+1, iz_prev+1) * (1-xd) + xi_prev(io, ix_prev+1, iy_prev+1, iz_prev+1) * xd

                        xi0 = xi00 * (1-yd) + xi10 * yd
                        xi1 = xi01 * (1-yd) + xi11 * yd

                        solvent(is)%xi(io,ix,iy,iz) = xi0 * (1-zd) + xi1 * zd

                        ix = ix + 1
                      end do
                      iy = iy + 1
                    end do
                    iz = iz + 1
                  end do

                end do
              end do
            end do

          end do


        end do
      deallocate (deltarho_p_prev)
      deallocate (xi_prev)
      else
        do is=1,size(solvent)
          ! read file
          do iz_prev=1,nz_prev
            do iy_prev=1,ny_prev
              do ix_prev=1,nx_prev
                do io=1,no_prev
                  read ( 10, iostat=ios ) solvent(is)%xi(io,ix_prev,iy_prev,iz_prev)
                end do
              end do
            end do
          end do
      end do
    end if
      close (10)

    end subroutine read_restart_file


    !
    ! Guess the init density since we don't have a restart file
    !
    subroutine guess_density
      use module_solvent, only: solvent
      use module_thermo, only: thermo
      implicit none
      real(dp) :: vextmax
      integer :: s
      !
      ! If vext is high, the guessed starting density is 0. If vext is something else, the guessed density is the bulk density (xi==1).
      !
      vextmax = (-log(epsilon(1.0_dp))-0.05) * thermo%kbT
      do s=1,size(solvent)
        where (solvent(s)%vext >= vextmax)
          solvent(s)%xi = 0._dp
        else where
          solvent(s)%xi = 1._dp
        end where
      end do
    end subroutine guess_density


end module module_density

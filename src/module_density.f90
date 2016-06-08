module module_density
    use iso_c_binding, only: c_double, c_float
    use precision_kinds, only: dp
    use module_solvent, only: solvent
    implicit none
    private
    public :: init_density
contains

    ! Init the density of each species as a function of position and orientation
    ! it should be initiated to exp(-beta*vext) but
    ! vextq is the electrostatic part of vext, and is pathologic (it sometimes diverges)
    ! we thus init the density not using vext, but vext - vextq
    subroutine init_density

        use module_thermo, only: thermo
        use module_solvent, only: solvent
        use module_grid, only: grid
        use module_input, only: getinput

        implicit none

        integer :: i, j, k, s, ios
        logical :: exists, vextq_is_allocated
        real(dp) :: v, threeshold_in_betav, betav
        real(dp), allocatable :: xi_loc(:)

        real(dp), allocatable :: xi_tmp(:, :, :)
        integer :: ix, iy, iz, ix_prev, iy_prev, iz_prev
        integer :: ns_prev, no_prev, io, is, no, ns
        integer :: mmax_prev, np_prev
        integer :: nx_prev, ny_prev, nz_prev, nx, ny, nz
        real(dp) :: dx_prev, dy_prev, dz_prev, dx, dy, dz
        real(dp) :: x, y, z, x_prev, y_prev, z_prev
        real(dp) :: x0, y0, z0, x1, y1, z1, xd, yd, zd, xi00, xi01, xi10, xi11, xi0, xi1

        ! Be sure the solvent is already initiated
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

        ! Read the density from a previous run
        if (getinput%log('restart', defaultvalue=.false.)) then
            inquire (file='input/density.bin', EXIST=exists)
            if (.not.exists) stop "You want to restart from a density file, but input/density.bin is not found"
            OPEN (10, file = 'input/density.bin' , form = 'unformatted' , iostat=ios, status='OLD' )
            IF ( ios /= 0 ) then
                print *, 'problem while opening input/density.bin. bug at init_density.f90'
                stop
            END IF

            READ ( 10, iostat=ios ) ns_prev
            READ ( 10, iostat=ios ) mmax_prev
            READ ( 10, iostat=ios ) no_prev
            READ ( 10, iostat=ios ) np_prev
            READ ( 10, iostat=ios ) nx_prev, ny_prev, nz_prev
            READ ( 10, iostat=ios ) dx_prev, dy_prev, dz_prev

            ! TODO: test si les ns sont les mêmes

            allocate(xi_tmp(nx_prev, ny_prev, nz_prev))

            dx=grid%dx
            dy=grid%dy
            dz=grid%dz
            nx=grid%nx
            ny=grid%ny
            nz=grid%nz

            ! https://en.wikipedia.org/wiki/Trilinear_interpolation
            do is=1,ns_prev
              do io=1,no_prev
                READ ( 10, iostat=ios ) xi_tmp

                do iz_prev=1,nz_prev-1
                  z0 = (iz_prev-1)*dz_prev
                  z1 = (iz_prev)*dz_prev
                  do iy_prev=1,ny_prev-1
                    y0 = (iy_prev-1)*dy_prev
                    y1 = (iy_prev)*dy_prev
                    do ix_prev=1,nx_prev-1
                      x0 = (ix_prev-1)*dx_prev
                      x1 = (ix_prev)*dx_prev

                      iz=ceiling(z0/dz+1)
                      do while(iz<=min(nz, floor(z1/dz+1)))
                        iy=ceiling(y0/dy+1)
                        do while(iy<=min(ny, floor(y1/dy+1)))
                          ix=ceiling(x0/dx+1)
                          do while(ix<=min(nx, floor(x1/dx+1)))


                            z = (iz-1)*dz
                            y = (iy-1)*dy
                            x = (ix-1)*dx

                            xd = (x-x0)/(x1-x0)
                            yd = (y-y0)/(y1-y0)
                            zd = (z-z0)/(z1-z0)

                            xi00 = xi_tmp(ix_prev, iy_prev, iz_prev) * (1-xd) + xi_tmp(ix_prev+1, iy_prev, iz_prev) * xd
                            xi01 = xi_tmp(ix_prev, iy_prev, iz_prev+1) * (1-xd) + xi_tmp(ix_prev+1, iy_prev, iz_prev+1) * xd
                            xi10 = xi_tmp(ix_prev, iy_prev+1, iz_prev) * (1-xd) + xi_tmp(ix_prev+1, iy_prev+1, iz_prev) * xd
                            xi11 = xi_tmp(ix_prev, iy_prev+1, iz_prev+1) * (1-xd) + xi_tmp(ix_prev+1, iy_prev+1, iz_prev+1) * xd

                            xi0 = xi00 * (1-yd) + xi10 * yd
                            xi1 = xi01 * (1-yd) + xi11 * yd
!                            print*,  io, iz, iy, ix, xi0 * (1-zd) + xi1 * zd
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

            close (10)
            return
        end if


        select case (dp)
        case (c_double)
            threeshold_in_betav = 36.04_dp
        case (c_float)
            threeshold_in_betav = 15.9_dp
        case default
            stop "In module_density the threeshold in the exponential before underflow is not given for the KIND of real you use"
        end select
!        threeshold_in_betav = vmax_before_underflow_in_exp_minus_vmax() ! 15.9 en real, 36.04 en double precision

        ! allocate (xi_loc(1:grid%no))
        ! s=1
        ! do k = 1, grid%nz
        !   do j = 1, grid%ny
        !     do i = 1, grid%nx
        !
        !       do io = 1, grid%no
        !         if ( thermo%beta*solvent(1)%vext(io,i,j,k) >= threeshold_in_betav ) then
        !           xi_loc(io) = 0._dp
        !         else
        !           xi_loc(io) = -sqrt(exp(-thermo%beta*solvent(1)%vext(io,i,j,k))) ! xi**2=rho/rho0=exp(-beta*v)
        !         end if
        !       end do
        !       solvent(1)%xi(:,i,j,k) = xi_loc
        !
        !     end do
        !   end do
        ! end do
        ! deallocate (xi_loc)


        ! where (solvent(1)%vext > 100._dp)
        !   solvent(1)%xi = 0._dp
        ! else where
        !   solvent(1)%xi = 1._dp
        ! end where

    end subroutine init_density


    pure function vmax_before_underflow_in_exp_minus_vmax ()
        implicit none
        integer :: i
        real(dp) :: v, vmax_before_underflow_in_exp_minus_vmax
        do i=1,1000
            v = real(i,dp)
            if (exp(-v)<=epsilon(v)) then
                vmax_before_underflow_in_exp_minus_vmax = v
                exit
            end if
        end do
    end function vmax_before_underflow_in_exp_minus_vmax

end module module_density

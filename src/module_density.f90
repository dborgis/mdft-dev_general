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

        integer :: i, j, k, io, s, ios
        logical :: exists, vextq_is_allocated
        real(dp) :: v, vextmax
        real(dp), allocatable :: xi_loc(:)

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
            READ ( 10, iostat=ios ) solvent(1)%xi
            IF ( ios<0 ) THEN
                error stop "input/density.bin is null"
            ELSE IF ( ios>0 ) THEN
                print*, "problem while trying to read solvent(1)%xi in input/density.bin"
                print*, "maybe you are reading a density with different nx,ny,nz,mmax?"
                error stop
            ELSE ! fine
              print*
              print*, '*** RESTARTING from density.bin ***'
              print*
            end if
            close (10)
            return
        end if

        !
        ! If vext is high, the guessed starting density is 0. If vext is something else, the guessed density is the bulk density (xi==1).
        !
        solvent(1)%xi=1._dp
        vextmax = vmax_before_underflow_in_exp_minus_vmax() * thermo%kbT
        where (solvent(1)%vext >= vextmax)
          solvent(1)%xi = 0._dp
        else where
          solvent(1)%xi = 1._dp
        end where

    end subroutine init_density

    !
    ! Determines the maximum value of x for which exp(-x) is numericaly computable, that is for which we don't have an IEEE-underflow
    !
    pure function vmax_before_underflow_in_exp_minus_vmax()
        implicit none
        integer :: i
        real(dp) :: v, vmax_before_underflow_in_exp_minus_vmax
        do i=1,10000
            v = real(i,dp)*0.1_dp
            if (exp(-v)<=epsilon(v)) then
                vmax_before_underflow_in_exp_minus_vmax = v-0.1
                exit
            end if
        end do
    end function vmax_before_underflow_in_exp_minus_vmax

end module module_density

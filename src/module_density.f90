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
        real(dp) :: v, threeshold_in_betav, betav

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

        ! allocate the solvent density field for each solvent species
        do s=1,solvent(1)%nspec
            if (.not. allocated( solvent(s)%density) ) then
                allocate(solvent(s)%density(grid%nx,grid%ny,grid%nz,grid%no), &
                    source=solvent(s)%rho0, stat=ios)
                if (ios /= 0) then
                    print*, "solvent(s)%density(grid%nx,grid%ny,grid%nz,grid%no), source=0._dp: Allocation request denied"
                    print*, "for s =", s
                end if
            end if
        end do

        ! Read the density from a previous run
        if (getinput%log('restart', defaultvalue=.false.)) then
            stop "LOOK AT MODULE_DENSITY"
            INQUIRE (file='input/density.bin.in', EXIST=exists)
            IF ( .NOT. exists) STOP "input/density.bin.in not found"
            OPEN (10, file = 'input/density.bin.in' , form = 'unformatted' , iostat=ios, status='OLD' )
            IF ( ios /= 0 ) then
                print *, 'problem while opening input/density.bin.in. bug at init_density.f90'
                stop
            END IF
            ! READ ( 10, iostat=ios ) cg_vect_new
            IF ( ios<0 ) THEN
                STOP "input/density.bin.in is empty"
            ELSE IF ( ios>0 ) THEN
                STOP "problem while trying to read cg_vect_new in input/density.bin.in"
            END IF
            PRINT*, '*** RESTART ***'
            CLOSE (10)
            OPEN (10, FILE = 'output/density.bin.in.out', FORM = 'unformatted')
            ! WRITE ( 10 ) cg_vect_new
            CLOSE (10)
            RETURN
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


        do concurrent (s=1: solvent(1)%nspec)
            vextq_is_allocated = allocated (solvent(s)%vextq)
            do io = 1, grid%no
                do k = 1, grid%nz
                    do j = 1, grid%ny
                        do i = 1, grid%nx

                            if (vextq_is_allocated) then
                                v = max( solvent(s)%vext(i,j,k,io), solvent(s)%vext(i,j,k,io) - solvent(s)%vextq(i,j,k,io) ) ! A VERIFIER
                            else
                                v = solvent(s)%vext(i,j,k,io)
                            end if

                            betav = thermo%beta * v

                            if ( betav >= threeshold_in_betav ) then
                                solvent(s)%density(i,j,k,io) = 0.0_dp
                            else
                                solvent(s)%density(i,j,k,io) = exp(-betav)*solvent(s)%density(i,j,k,io)
                            end if

                        end do
                    end do
                end do
            end do
        end do
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

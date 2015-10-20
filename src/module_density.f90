module module_density
    use precision_kinds, only: dp
    use module_solvent, only: solvent
    implicit none
    private
    public :: init_density
contains

    ! Init the density of each species as a function of position and orientation
    ! it should be initiated to exp(-beta*Vext_total) but
    ! Vext_q is the electrostatic part of Vext_total, and is pathologic (it sometimes diverges)
    ! we thus init the density not using vext, but Vext_total - Vext_q
    subroutine init_density

        use system, only: thermocond
        use module_solvent, only: solvent
        use module_grid, only: grid
        use external_potential, only: Vext_total, Vext_q
        use module_input, only: getinput

        implicit none

        integer :: i, j, k, io, s, ios
        logical :: exists, is_vext_q_allocated
        real(dp) :: v, threeshold_in_betav, betav

        ! allocate the solvent density field for each solvent species
        do s=1,solvent(1)%nspec
            if (.not. allocated( solvent(s)%density) ) then
                allocate(solvent(s)%density(grid%nx,grid%ny,grid%nz,grid%no), source=0._dp, stat=ios)
                if (ios /= 0) then
                    print*, "solvent(s)%density(grid%nx,grid%ny,grid%nz,grid%no), source=0._dp: Allocation request denied"
                    print*, "for s =", s
                end if
            end if
        end do

        ! Read the density from a previous run
        if (getinput%log('reuse_density', defaultvalue=.false.)) then
            INQUIRE (file='input/density.bin.in', EXIST=exists)
            IF ( .NOT. exists) STOP "input/density.bin.in not found"
            OPEN (10, file = 'input/density.bin.in' , form = 'unformatted' , iostat=ios, status='OLD' )
            IF ( ios /= 0 ) then
                print *, 'problem while opening input/density.bin.in. bug at init_density.f90'
                stop
            END IF
            READ ( 10, iostat=ios ) cg_vect_new
            IF ( ios<0 ) THEN
                STOP "input/density.bin.in is empty"
            ELSE IF ( ios>0 ) THEN
                STOP "problem while trying to read cg_vect_new in input/density.bin.in"
            END IF
            PRINT*, '*** RESTART ***'
            CLOSE (10)
            OPEN (10, FILE = 'output/density.bin.in.out', FORM = 'unformatted')
            WRITE ( 10 ) cg_vect_new
            CLOSE (10)
            RETURN
        end if


        ! only Vext_total is used in the functional. We may deallocate it now.
        !if ( allocated ( Vext_q ) ) deallocate ( Vext_q )

        is_vext_q_allocated = allocated (vext_q)
        threeshold_in_betav = vmax_before_underflow_in_exp_minus_vmax()

        do s = 1, solvent(1)%nspec
            do io = 1, grid%no
                do k = 1, grid%nz
                    do j = 1, grid%ny
                        do i = 1, grid%nx

                            if ( is_vext_q_allocated ) then
                                v = max( Vext_total(i,j,k,io,s) - Vext_q (i,j,k,io,s) )
                            else
                                v = Vext_total(i,j,k,io,s)
                            end if

                            betav = thermocond%beta * v

                            if ( betav >= threeshold_in_betav ) then
                                solvent(s)%density(i,j,k,io) = 0.0_dp
                            else
                                solvent(s)%density(i,j,k,io) = exp(-betav)
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

end module density_mod

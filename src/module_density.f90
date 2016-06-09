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


    !
    ! Read the restart file
    !
    subroutine read_restart_file
      use module_solvent, only: solvent
      implicit none
      logical :: exists
      integer :: ios
      character(len("input/density.bin")), parameter :: filename="input/density.bin"
      inquire (FILE=filename, EXIST=exists)
      if(.not.exists) error stop "You want to restart from a density file, but input/density.bin is not found"
      open(10, file = 'input/density.bin' , form = 'unformatted' , iostat=ios, status='OLD' )
      if(ios /= 0) error stop 'problem while opening input/density.bin. bug at init_density.f90'
      read(10, iostat=ios ) solvent(1)%xi
      if( ios<0 ) then
          error stop "input/density.bin is null"
      else if( ios>0 ) then
          print*, "problem while trying to read solvent(1)%xi in input/density.bin"
          print*, "maybe you are reading a density with different nx,ny,nz,mmax?"
          error stop
      else
        print*
        print*, '*** RESTARTING from input/density.bin ***'
        print*
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
      vextmax = vmax_before_underflow_in_exp_minus_vmax() * thermo%kbT
      do s=1,size(solvent)
        where (solvent(s)%vext >= vextmax)
          solvent(s)%xi = 0._dp
        else where
          solvent(s)%xi = 1._dp
        end where
      end do
    end subroutine guess_density
    

end module module_density

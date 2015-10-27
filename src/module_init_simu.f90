module module_init_simu

    implicit none

    private
    public :: init_simu

contains

    !
    !   Init the simulation
    !
    subroutine init_simu
        use module_grid, only: grid
        use module_solute, only: read_solute
        use module_solvent, only: init_solvent => read_solvent
        use module_density, only: init_density
        use hardspheres, only: compute_hard_spheres_parameters
        use module_input, only: getinput
        use dcf, only: init_dcf => init
        use fft, only: init_fft => init
        use module_vext, only: init_vext
        use module_debug, only: init_debug
        use module_thermo, only: init_thermo
        implicit none
        call print_header
        call init_debug
        call grid%init
        call read_solute
        call init_thermo
        call init_solvent
        call init_dcf
        if (getinput%log('hard_sphere_fluid', defaultvalue=.false.)) call compute_hard_spheres_parameters ! If calculation based on Fundamental Measure Theory read FMT parameters and compute weight functions etc
        call init_fft
        call init_vext
        call init_density
    end subroutine init_simu

    !
    !   Print a header to stdout to welcome the user
    !.
    subroutine print_header
        character(8)  :: date
        character(10) :: time
        call date_and_time ( DATE=date,TIME=time)
        print*,' ******************************************************************'
        print*,' ******************************************************************'
        print*,' **             *************************   ',date(1:4),'/',date(5:6),'/',date(7:8),'   **********'
        print*,' **   M D F T   *************************                **********'
        print*,' **             *************************    ',time(1:2),'h',time(3:4),'m',time(5:6),'    **********'
        print*,' ******************************************************************'
        print*,' ******************************************************************'
    end subroutine print_header


end module module_init_simu

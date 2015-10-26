module module_init_simu

    use fft, only: init_fft => init
    use dcf, only: init_dcf => init
    use module_input, only: getinput
    use hardspheres, only: compute_hard_spheres_parameters
    use module_density, only: init_density
    use module_solvent, only: init_solvent => read_solvent
    use module_solute, only: read_solute

    implicit none

    private
    public :: init_simu

contains

    !
    !   Init the simulation
    !
    subroutine init_simu
        use module_grid, only: grid
        use module_quadrature, only: init_quadrature => init
        implicit none

        call print_header ! package name, day & time
        call grid%init
        call init_quadrature
        call read_solute
        call init_thermo
        call init_solvent
        call init_dcf
        if (getinput%log('hard_sphere_fluid', defaultvalue=.false.)) call compute_hard_spheres_parameters ! If calculation based on Fundamental Measure Theory read FMT parameters and compute weight functions etc
        call init_fft
        call init_external_potential
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

    ! Here are allocated variables declared in modules
    subroutine init_thermo
        use precision_kinds     ,only: dp
        use module_input        ,only: getinput
        use system              ,only: thermoCond
        use constants           ,only: boltz, navo
        implicit none
        print*,
        print*, "[Thermodynamics]====="
        thermoCond%T = getinput%dp('temperature', defaultvalue=300._dp, assert=">0") ! look for temperature in input
        thermoCond%kbT = Boltz * Navo * thermoCond%T * 1.0e-3_dp
        thermoCond%beta = 1.0_dp / thermocond%kbT
        print*, "Temperature (K)  :", thermocond%t
        print*, "kT               :", thermocond%kbt
        print*, "\beta = 1/(kT)   :", thermocond%beta
        print*, "[/Thermodynamics]====="
        print*,
    end subroutine init_thermo


end module module_init_simu

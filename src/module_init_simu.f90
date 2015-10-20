module module_init_simu

    use quadrature, only: prepare_quadrature => init
    use fft, only: init_fft => init
    use dcf, only: init_dcf => init
    use module_input, only: getinput
    use hardspheres, only: compute_hard_spheres_parameters
    use module_density, only: init_density
    use module_solvent, only: init_solvent => read_solvent

    implicit none

    private
    public :: init_simu

contains

    !
    !   Init the simulation
    !
    subroutine init_simu

        implicit none

        call print_header ! package name, day & time
        call read_solute ! AWFUL TODO trick to default lx
        call allocate_from_input ! TODO this should be removed and done only when needed ! read dft.in, solute.in and solvent.in
        call prepare_quadrature
        call init_fft
        call init_solvent ! Read solvent atomic positions, charge and Lennard-Jones param
        call init_dcf
        if (getinput%log('hard_sphere_fluid', defaultvalue=.false.)) call compute_hard_spheres_parameters ! If calculation based on Fundamental Measure Theory read FMT parameters and compute weight functions etc
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


end module module_init_simu

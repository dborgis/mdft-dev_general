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
        use module_fft, only: init_fft => init
        use module_vext, only: init_vext
        use module_thermo, only: init_thermo
        implicit none
        call print_header
        call grid%init
        call read_solute
        call init_thermo
        call init_solvent
        if (getinput%log('hard_sphere_fluid', defaultvalue=.false.)) call compute_hard_spheres_parameters ! If calculation based on Fundamental Measure Theory read FMT parameters and compute weight functions etc
        call init_fft
        call init_vext
        call init_vexc
        call init_density
    end subroutine init_simu

    !
    !   Print a header to stdout to welcome the user
    !.
    subroutine print_header

        use git, only: commit

        character(8)  :: date
        character(10) :: time
        call date_and_time ( DATE=date,TIME=time)
        print*,' ******************************************************************'
        print*,' ******************************************************************'
        print*,' **             ******    ',date(1:4),'/',date(5:6),'/',date(7:8),'   ',time(1:2),'h',time(3:4),'m',time(5:6),'                  **'
        print*,' **   M D F T   ******                                           **'
        print*,' **             ******    commit     ', commit(1:10),'                  **'
        print*,' ******************************************************************'
        print*,' ******************************************************************'
        print*

    end subroutine print_header

    subroutine init_vexc
    
      use proc_ptrs, only : excess_energy
      use module_solvent, only: solvent
      use module_energy_cproj_mrso, only: energy_cproj_mrso
      use module_energy_cproj_mrso_mixture, only: energy_cproj_mrso_mixture

      if (size(solvent)==1) then
       excess_energy=>energy_cproj_mrso
      else
        excess_energy=>energy_cproj_mrso_mixture
      end if
    !
    end subroutine

end module module_init_simu

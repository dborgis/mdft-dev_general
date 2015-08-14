module freeEnergyFunctional_mod

    use iso_c_binding, only: dp => c_double
    implicit none
    private

    type, public :: ffideal_type
        real(dp) :: energy
        real(dp), allocatable :: grad(:,:,:,:,:)
        contains
        procedure :: is_built => ffideal_test_is_built
        procedure :: get => ffideal_compute

    end type

contains

    pure function ffideal_test_is_built (this) result (is_built)
        implicit none
        logical :: is_built
        class (ffideal_type), intent(in) :: this
        if (allocated (this%grad)) then
            is_built = .true.
        else
            is_built = .false.
        end if
    end function

    function ffideal_compute (this) result (out)
        class (ffideal_type), intent(inout) :: this
        type (ffideal_type) :: out
        if (this%is_built) then
            this%
        else
            stop "this is not built"
            ! this%build
        end if
    end function

end module


faut definir un    type qui va contenir     FF et dF et pourquoi pas d2F

type :: free_energy_functional_type
    double precision :: energy
    double precision, allocatable :: grad(:)
end type

type (free_energy_functional_type) :: free_energy_functional

free_energy_functional = free_energy_functional + free_energy_functional_ideal (density,thermodynamic_conditions)

type :: energy_and_gradient
    double precision :: energy
    double precision, allocatable :: grad(:)
end type
type (value_and_gradient) :: ideal_free_energy

ideal_free_energy%get


pure function free_energy_functional_ideal (density, thermodynamic_conditions)
    use free_energy_functional_mod, only: free_energy_functional_type
    use thermodynamic_conditions_mod, only: thermodynamic_conditions_type
    implicit none
    type (free_energy_functional) :: free_energy_functional_ideal
    double precision, intent(in) :: density(:,:,:,:,:)
    type (thermodynamic_conditions_type), intent(in) :: thermodynamic_conditions
    ! local
    double precision :: temperature, kb, kbT
    temperature = thermodynamic_conditions%temperature
    kb = thermodynamic_conditions%boltzmann_constant
    kbT = kb*temperature
    free_energy_functional%energy = kbT*sum( density*log(density)-density+density_ref ,mask=density>0)

end function



free_energy_functional = opt_results

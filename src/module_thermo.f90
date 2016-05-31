module module_thermo

    use precision_kinds, only: dp
    implicit none
    private

    type :: thermo_type
        real(dp) :: t ! temperature
        real(dp) :: p ! pressure
        real(dp) :: kbt ! temperature energy unit
        real(dp) :: beta ! 1/kbt
    end type
    type (thermo_type), protected :: thermo ! everything related to the thermodynamic conditions
    public :: thermo, init_thermo

contains
        !
        ! Thermodynamics : temperature etc.
        !
        subroutine init_thermo
            use precision_kinds     ,only: dp
            use module_input        ,only: getinput
            use constants           ,only: boltz, navo
            implicit none
            print*
            print*, "===== Thermodynamics ====="
            thermo%T = getinput%dp('temperature', defaultvalue=300._dp, assert=">0")
            thermo%kbT = Boltz * Navo * thermo%T * 1.0e-3_dp
            thermo%beta = 1.0_dp / thermo%kbT
            thermo%P = getinput%dp('pressure', defaultvalue=1._dp, assert=">=0")
            print*, "   T  =", real(thermo%t), "K"
            print*, "   kT =", real(thermo%kbt),"kJ/mol"
            print*, "   β  = (kT)⁻¹ =", real(thermo%beta),"mol/kJ"
            print*, "===== Thermodynamics ====="
            print*
        end subroutine init_thermo

end module module_thermo

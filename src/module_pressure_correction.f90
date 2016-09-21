module module_pressure_correction
    implicit none
    private
    public :: pressure_correction
contains

subroutine pressure_correction
! ... See volodymyr et al., JPCL etc

    use precision_kinds, only: dp
    use module_energy_and_gradient, only: energy_and_gradient
    use module_thermo, only: thermo
    use module_solute, only: solute
    use module_solvent, only: solvent
    use module_energy_and_gradient, only: ff
    use module_input, only: getinput
    use module_grid, only: grid
    implicit none

    !... Volodymyr's partial molar volume correction. See J. Phys. Chem. Lett. 5, 1935-1942 (2014)
    real(dp) :: deltaGtotMDFT, Pbulk
    real(dp), allocatable :: density(:,:,:)
    integer :: nx, ny, nz, s, ios
    real(dp) :: nmolecules_with_solute
    real(dp) :: nmolecules_without_solute
    real(dp) :: gamma ! trace of the tensor of the quadrupole moment
    real(dp) :: numberdensity ! molecular number density, for instance 0.0332891 molecule per angstrom^3
    real(dp) :: solutecharge ! net charge of the solute, for instance -1 for Cl- ion
    real(dp), parameter :: kJpermolperang3_to_Pa = 1.66113*10**9
    real(dp), parameter :: Pa_to_atm = 9.8692327e-06
    real(dp) :: deltaG_corrected, deltaN, Pscheme_correction, correction, correction2
    real(dp), parameter :: zerodp = 0._dp
    real(dp) :: deltaG_emptybox


    numberdensity = solvent(1)%n0
    deltaGtotMDFT = ff%tot

    !... Volodymyr's partial molar volume correction. See J. Phys. Chem. Lett. 5, 1935-1942 (2014)
    if (size(solvent)/=1) stop "pressure_correction implemented for 1 solvent species only"
    allocate(  density(grid%nx,grid%ny,grid%nz),  source=0._dp)
    call grid%integrate_over_orientations( solvent(1)%xi**2 * solvent(1)%rho0, density)

    ! number of solvent molecules inside the supercell:
    ! when the supercell has the solute inside
    nmolecules_with_solute = sum(density)*grid%dv
    ! when the supercell is empty of any perturbation, ie pure solvent
    nmolecules_without_solute = solvent(1)%n0*product(grid%length)

    write(*,'(A,F12.2)') "Solvent molecules in the supercell containing the solute   ", nmolecules_with_solute
    write(*,'(A,F12.2)') "Solvent molecules in a supercell without solute", nmolecules_without_solute
    deltaN = nmolecules_without_solute - nmolecules_with_solute
    write(*,'(A,F12.2)') "ΔN solvent", deltaN
    write(*,'(A,F12.7,A)') "Bulk solvent density", solvent(1)%n0," molecule.Ang⁻³"
    write(*,'(A,F12.5,A,F9.2,A)') "Supercell volume", product(grid%length)," Ang³ = ",product(grid%length)/1000.," nm³"

    ! pressure of the bulk solvent?  GrandPotential(homogeneous system) = -PV
    ! grand potential[rho_bulk] == PV
    solvent(1)%xi = 0 ! set Density to 0
    call energy_and_gradient(deltaG_emptybox) ! this step is not a minimization step so we give a negative integeration number to avoid the printing of the not relevant obtained energies

    Pbulk = deltaG_emptybox / (grid%lx * grid%ly * grid%lz) ! Omega[rho=rho_0]=PV ! Pbulk in kJ/mol/Ang^3
    write(*,'(A,F12.2,A)') "Bulk pressure       ", Pbulk*kJpermolperang3_to_Pa*Pa_to_atm," atm"
    open(81,file="output/pressure")
    write(81,*) Pbulk
    close(81)

    correction  = -deltaN/solvent(1)%n0*Pbulk  !correction is -PV where V is excluded Volume
    correction2 =  deltaN*thermo%kbT  !correction is -PV where V is excluded Volume
    deltaG_corrected = deltaGtotMDFT + correction !+ correction2
    write(*,'(A,F12.2,A)') "PMV correction      ", correction," kJ/mol"
    write(*,'(A,F12.2,A)') "Pid correction      ", correction2," kJ/mol"
    open(79,file="output/PMV_correction")
    write(79,*) correction
    close(79)
    open(80,file="output/Pideal_PMV_correction")
    write(80,*) correction2
    close(80)

    write(*,'(A,F12.2,A)') "SFE ISc /PC         ", deltaG_corrected," kJ/mol"
    write(*,'(A,F12.2,A)') "SFE ISc*/PC+        ", deltaG_corrected + correction2," kJ/mol"

    PRINT*, "ESTIMATED SOLVATION FREE ENERGY: ", deltaG_corrected + correction2

end subroutine pressure_correction

end module module_pressure_correction

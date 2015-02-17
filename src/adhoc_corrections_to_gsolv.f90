subroutine adhoc_corrections_to_gsolv
! ... Here, we print all the adhoc corrections one should take into account before comparing MDFT results to MD and/or experiments.

    use precision_kinds, only: dp, sp, i2b
    use system, only: solute, solvent, spacegrid, thermocond
    use minimizer, only: FF , cg_vect
    use constants, only: zerodp
    use mathematica, only: chop
    use input, only: input_log
    implicit none

    real(dp) :: correction,correction2, Pbulk
    real(dp), allocatable :: neq(:,:,:,:) ! equilibrium density
    integer(i2b), pointer :: nfft1 => spacegrid%n_nodes(1), nfft2 => spacegrid%n_nodes(2), nfft3 => spacegrid%n_nodes(3)
    integer :: s, ios
    type :: nmoleculetype
        real(dp) :: withsolute
        real(dp) :: bulk
    end type nmoleculetype
    type (nmoleculetype), allocatable :: nmolecule(:)
    logical :: file_exists
    real(dp) :: gamma ! quadrupole moment trace
    real(dp) :: numberdensity ! molecular number density, for instance 0.0332891 molecule per angstrom^3
    real(dp) :: solutecharge ! net charge of the solute, for instance -1 for Cl- ion
    real(dp), parameter :: kJpermolperang3_to_Pa = 1.66113*10**9
    real(dp), parameter :: Pa_to_atm = 9.8692327e-06
    real(dp) :: FFcorrected_final, deltaN, Pscheme_correction

    FFcorrected_final = FF

    open(79,file="output/FF"); write(79,*) FF; close(79)



    !... We use P-scheme instead of M-scheme for the electrostatics in MDFT.
    ! See Kastenholz and Hunenberger, JCP 124, 124106 (2006), page 224501-8, equations 35, 35 and 37 with Ri=0
    ! "To be applied if the solvent molecule is rigid and involves a single van der Waals interaction site M,
    ! and that any scheme relying on molecular-cutoff truncation refers to this specific site for applying the truncation."
    Pscheme_correction = 0._dp
    if( input_log("poisson_solver") ) then
      gamma = solvent(1)%quadrupole(1,1)+solvent(1)%quadrupole(2,2)+solvent(1)%quadrupole(3,3) ! quadrupole moment trace
      numberdensity = solvent(1)%n0
      solutecharge = sum(solute%site%q)
      Pscheme_correction = -gamma*numberdensity*2.909857E3*solutecharge ! in kJ/mol
    end if
    write(*,'(A,F12.2,A)') "P-scheme correction ", Pscheme_correction," kJ/mol"
    open(79,file="output/Pscheme_correction"); write(79,*) Pscheme_correction; close(79)
    FFcorrected_final = FFcorrected_final + Pscheme_correction




  !... Volodymyr's partial molar volume correction. See J. Phys. Chem. Lett. 5, 1935-1942 (2014)
  correction = zerodp
  if (size(solvent)/=1) stop "CRITICAL in adhoc_corrections_to_gsolv. only 1 solvent species is implemented."
  do s=1,size(solvent)
      if (.not. allocated( solvent(s)%n )) allocate ( solvent(s)%n(nfft1,nfft2,nfft3) ,source=0._dp)
      call get_final_density ( solvent(s)%n , s) ! Get the final density(r) from the last minimizer step.
  end do
  allocate (nmolecule(size(solvent)))
  do concurrent (s=1:size(solvent))
      nmolecule%withsolute = sum(solvent(s)%n * solvent(s)%n0)  *spacegrid%dv ! number of solvent molecules inside the supercell containing the solute
  end do
  nmolecule%bulk = solvent%n0*product(spacegrid%length) ! number of solvent molecules inside the same supercell (same volume) without solute.
  write(*,'(A,F12.2)') "Solvent molecules with solute   ", nmolecule%withsolute
  write(*,'(A,F12.2)') "Solvent molecules without solute", nmolecule%bulk
  write(*,'(A,F12.2)') "ΔN solvent", nmolecule(s)%bulk - nmolecule(s)%withsolute
  write(*,'(A,F12.7,A)') "Solvent density", solvent(1)%n0," molecule.Ang⁻³"
  write(*,'(A,F12.5,A)') "Supercell volume", spacegrid%dv," Ang⁻³"



  ! pressure of the bulk solvent?  GrandPotential(homogeneous system) = -PV
  ! grand potential[rho_bulk] == PV
  cg_vect(:) = zerodp ! set Density to 0
  FF = zerodp         ! set energy to 0
  call energy_and_gradient(-10) ! this step is not a minimization step so we give a negative integeration number to avoid the printing of the not relevant obtained energies
  Pbulk = FF/product(spaceGrid%length) ! Omega[rho=rho_0]=PV ! Pbulk in kJ/mol/Ang^3
  write(*,'(A,F12.2,A)') "Bulk pressure       ", Pbulk*kJpermolperang3_to_Pa*Pa_to_atm," atm"
  open(81,file="output/bulk-pressure"); write(81,*) Pbulk; close(81)




  s = 1
  if( s /= 1 ) stop "line 61 of adhoc_corr... we have not thought of multi species case"
  deltaN = nmolecule(s)%bulk - nmolecule(s)%withsolute
  correction  = -(nmolecule(s)%bulk - nmolecule(s)%withsolute)/solvent(s)%n0*Pbulk  !correction is -PV where V is excluded Volume
  correction2 =  (nmolecule(s)%bulk - nmolecule(s)%withsolute)*thermoCond%kbT  !correction is -PV where V is excluded Volume
  FFcorrected_final = FFcorrected_final + correction !+ correction2
  write(*,'(A,F12.2,A)') "PMV correction      ", correction," kJ/mol"
  write(*,'(A,F12.2,A)') "Pid correction      ", correction2," kJ/mol"
  open(79,file="output/PMV_correction")
    write(79,*) correction
  close(79)
  open(80,file="output/Pideal_PMV_correction")
    write(80,*) correction2
  close(80)

  write(*,'(A,F12.2,A)') "SFE ISc             ", FFcorrected_final," kJ/mol"
  write(*,'(A,F12.2,A)') "SFE ISc*            ", FFcorrected_final + correction2," kJ/mol"

end subroutine adhoc_corrections_to_gsolv

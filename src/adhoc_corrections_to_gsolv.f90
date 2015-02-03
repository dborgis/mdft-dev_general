subroutine adhoc_corrections_to_gsolv
! ... Here, we print all the adhoc corrections one should take into account before comparing MDFT results to MD and/or experiments.

    use precision_kinds, only: dp, sp, i2b
    use system, only: solute, solvent, spacegrid, thermocond
    use minimizer, only: FF , cg_vect, finalizeMinimizer
    use constants, only: zerodp
    use mathematica, only: chop
    use input, only: input_log
    implicit none

    real(dp) :: correction,correction2, Pressure_bulk
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

    open(79,file="output/FF")
      write(79,*) FF
    close(79)

    !... We use P-scheme instead of M-scheme for the electrostatics in MDFT.
    ! See Kastenholz and Hunenberger, JCP 124, 124106 (2006), page 224501-8, equations 35, 35 and 37 with Ri=0
    ! "To be applied if the solvent molecule is rigid and involves a single van der Waals interaction site M,
    ! and that any scheme relying on molecular-cutoff truncation refers to this specific site for applying the truncation."
    correction = 0._dp
    if( input_log("poisson_solver") ) then
      gamma = solvent(1)%quadrupole(1,1)+solvent(1)%quadrupole(2,2)+solvent(1)%quadrupole(3,3) ! quadrupole moment trace
      numberdensity = solvent(1)%n0
      solutecharge = sum(solute%site%q)
      correction = -gamma*numberdensity*2.909857E3*solutecharge ! in kJ/mol
    end if
    print*,"You should add",real(correction,sp),"kJ/mol to FREE ENERGY because we use the P-scheme electrostatics"
    open(79,file="output/Pscheme_correction")
      write(79,*) correction
    close(79)

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

! The value of the grand potential is equal to PV, with P pressure and V volume when the system is the bulk fluid.
  cg_vect(:) = zerodp  !Set Density to 0.0_dp
  FF = zerodp
  Call energy_and_gradient(-10)  !this step is not a minimization step so we give a negative integeration number to avoid the printing of the not relevant obtained energies

  Pressure_bulk=FF/PRODUCT(spaceGrid%length) ! Omega[rho=rho_0]=PV
  print*, 'Pressurebulk=' , pressure_bulk*kJpermolperang3_to_Pa , "Pa"

  s = 1
  if( s /= 1 ) stop "line 61 of adhoc_corr... we have not thought of multi species case"
  correction  = -(nmolecule(s)%bulk - nmolecule(s)%withsolute)/solvent(s)%n0*Pressure_bulk  !correction is -PV where V is excluded Volume
  correction2 =  (nmolecule(s)%bulk - nmolecule(s)%withsolute)*thermoCond%kbT  !correction is -PV where V is excluded Volume
  print*,"You should add",correction,"kJ/mol to FREE ENERGY as partial molar volume correction" !
  print*,"You should add",correction2,"kJ/mol to FREE ENERGY as ideal partial molar volume correction" !
  open(79,file="output/PMV_correction")
  write(79,*) correction
  close(79)
  open(80,file="output/Pideal_PMV_correction")
  write(80,*) correction2
  close(80)

end subroutine adhoc_corrections_to_gsolv

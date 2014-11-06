subroutine adhoc_corrections_to_gsolv
! ... Here, we print all the adhoc corrections one should take into account before comparing MDFT results to MD and/or experiments.

    use precision_kinds, only: dp, sp, i2b
    use system, only: solute, solvent, spacegrid, thermocond
    use minimizer, only: FF , cg_vect, finalizeMinimizer

    use mathematica, only: chop
    use input, only: input_log
    implicit none

    real(dp) :: correction,  Pressure_bulk
    real(dp), allocatable :: neq(:,:,:,:) ! equilibrium density
    integer(i2b), pointer :: nfft1 => spacegrid%n_nodes(1), nfft2 => spacegrid%n_nodes(2), nfft3 => spacegrid%n_nodes(3)
    integer :: s, ios
    type :: nmoleculetype
        real(dp) :: withsolute
        real(dp) :: bulk
    end type nmoleculetype
    type (nmoleculetype), allocatable :: nmolecule(:)
    logical :: file_exists

    open(79,file="output/FF")
      write(79,*) FF
    close(79)
    !... We use P-scheme instead of M-scheme for the electrostatics in MDFT. See Kastenholz and Hunenberger, JCP 124, 124106 (2006)
    if (input_log("poisson_solver") .eqv. .true.) then
      correction = -79.8_dp*sum(solute%site%q) ! in kJ/mol
      correction = chop(correction)
    else
      correction = 0._dp
    end if
    print*,"You should add",real(correction,sp),"kJ/mol to FREE ENERGY because we use the P-scheme electrostatics"
    open(79,file="output/Pscheme_correction")
      write(79,*) correction
    close(79)

    !... Volodymyr's partial molar volume correction. See J. Phys. Chem. Lett. 5, 1935-1942 (2014)
    correction = 0._dp
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
  cg_vect(:)=0.0_dp  !Set Density to 0.0_dp
  FF=0.0_dp
  Call energy_and_gradient(-10)  !this step is not a minimization step so we give a negative integeration number to avoid the printing of the not relevant obtained energies

  Pressure_bulk=FF/PRODUCT(spaceGrid%length) ! Omega[rho=rho_0]=PV

              correction=-(nmolecule(1)%bulk - nmolecule(1)%withsolute)/solvent(1)%n0*Pressure_bulk  !correction is -PV where V is excluded Volume
              print*,"You should add",correction,"kJ/mol to FREE ENERGY as partial molar volume correction" !
              open(79,file="output/PMV_correction")
              write(79,*) correction
              close(79)



end subroutine adhoc_corrections_to_gsolv

subroutine adhoc_corrections_to_gsolv
! ... Here, we print all the adhoc corrections one should take into account before comparing MDFT results to MD and/or experiments.
  
    use precision_kinds, only: dp, sp, i2b
    use system, only: solute, solvent, spacegrid, thermocond
  
    implicit none

    real(dp) :: correction
    real(dp), allocatable :: neq(:,:,:,:) ! equilibrium density
    integer(i2b), pointer :: nfft1 => spacegrid%n_nodes(1), nfft2 => spacegrid%n_nodes(2), nfft3 => spacegrid%n_nodes(3)
    integer :: s, ios
    type :: nmoleculetype
        real(dp) :: withsolute
        real(dp) :: bulk
    end type nmoleculetype
    type (nmoleculetype), allocatable :: nmolecule(:)
    logical :: file_exists

    !... We use P-scheme instead of M-scheme for the electrostatics in MDFT. See Kastenholz and Hunenberger, JCP 124, 124106 (2006)
    correction = 0._dp
    if (  abs(sum(solute%site%q)) >= epsilon(1.0_dp)  ) then
        correction = -79.8_dp*sum(solute%site%q) ! in kJ/mol
        print*,"You should add",real(correction,sp),"kJ/mol to FREE ENERGY because we use the P-scheme electrostatics"
    end if


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
    
    inquire(file='output/cs.in',exist=file_exists);
    if (file_exists) then
        open (14, file='output/cs.in', iostat=ios)
        if (ios/=0) stop 'Cant open file output/cs.in in adhoc_corrections_to_gsolv. But it is present!'
        if (ios==0) then
            block
                real(dp) :: knorm, ck0
                read(14,*) knorm, ck0
                correction = (nmolecule(1)%bulk - nmolecule(1)%withsolute) *thermoCond%kbT * (-1._dp + 0.5_dp*solvent(1)%n0**2* ck0)
                print*,"You should add",real(correction,sp),"kJ/mol to FREE ENERGY as partial molar volume correction"
            end block
            close(14)
        end if
    else
        print*,"I could not find c(k=0) by myself. Please do the math yourself:"
        print*,"You should add",(nmolecule(1)%bulk - nmolecule(1)%withsolute) *thermoCond%kbT,&
        "*(-1+", 0.5_dp*solvent(1)%n0**2 ,"*c(k=0) ) kJ/mol to FREE ENERGY as partial molar volume correction"
    end if

end subroutine adhoc_corrections_to_gsolv
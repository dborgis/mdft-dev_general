subroutine adhoc_corrections_to_gsolv
! ... Here, we print all the adhoc corrections one should take into account before comparing MDFT results to MD and/or experiments.

    use precision_kinds, only: dp, sp, i2b
    use system, only: solute, solvent, spacegrid, thermocond
    use minimizer, only: FF
    use mathematica, only: chop

    implicit none

    real(dp) :: correction, correction2
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
    correction = -79.8_dp*sum(solute%site%q) ! in kJ/mol
    correction = chop(correction)
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

    inquire(file='output/cs.in',exist=file_exists);
    if (file_exists) then
        open (14, file='output/cs.in', iostat=ios)
        if (ios/=0) stop 'Cant open file output/cs.in in adhoc_corrections_to_gsolv. But it is present!'
        if (ios==0) then
            block
                real(dp) :: knorm, ck0
                read(14,*) knorm, ck0
                correction = (nmolecule(1)%bulk - nmolecule(1)%withsolute) *thermoCond%kbT * (-1._dp + 0.5_dp*solvent(1)%n0* ck0)
                correction2 = (nmolecule(1)%bulk - nmolecule(1)%withsolute) *thermoCond%kbT * (-1._dp + 0.5_dp*solvent(1)%n0* ck0&
              +9.2000475144588751) !MAGIC NUMBER FOR 3BODY PRESSURE TODO EXPLAIN WHAT IT IS
              print*,"You should add",correction,"kJ/mol to FREE ENERGY as partial molar volume correction" !
              open(79,file="output/PMV_correction")
                write(79,*) correction
              close(79)
              print*,"You should add",correction2,"kJ/mol to FREE ENERGY as partial molar volume correction if you work with 3Body"!
              open(79,file="output/PMV3b_correction")
                write(79,*) correction2
              close(79)
              !  print*,"You should add",correction2,"kJ/mol to FREE ENERGY as partial molar volume correction"
            end block
            close(14)
        end if
    else ! if file does not exist
        print*,"I could not find c(k=0) by myself. Please do the math yourself:"
        print*,"You should add",(nmolecule(1)%bulk - nmolecule(1)%withsolute) *thermoCond%kbT,&
        "*(-1+", 0.5_dp*solvent(1)%n0 ,"*c(k=0) ) kJ/mol to FREE ENERGY as partial molar volume correction"
    end if

  open(79,file="output/FF")
    write(79,*) FF
  close(79)

end subroutine adhoc_corrections_to_gsolv

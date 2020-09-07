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
    use module_solvent, only: solvent
    use module_energy_and_gradient, only: ff
    use module_grid, only: grid
    use module_fft, only: sum_byfft_3d
    use module_input, only: getinput
    use module_lennardjones, only: LJ_energy_correction
    use module_energy_ideal_and_external, only: energy_ideal_and_external
    use module_coarse_grained_bridge, only: coarse_grained_bridge, cgb
    use module_energy_cproj_mrso, only: c000
    use hardspheres      ,ONLY: hs

    implicit none

    !... Volodymyr's partial molar volume correction. See J. Phys. Chem. Lett. 5, 1935-1942 (2014)
    real(dp) :: deltaGtotMDFT, Pbulk, deltaG_cavity, fext, fid, fb_wda, fb_sla
    real(dp), allocatable :: density(:,:,:)
    real(dp) :: nmolecules_with_solute
    real(dp) :: nmolecules_without_solute
    real(dp) :: numberdensity ! molecular number density, for instance 0.0332891 molecule per angstrom^3
    real(dp), parameter :: kJpermolperang3_to_Pa = 1.66113*10**9
    real(dp), parameter :: Pa_to_atm = 9.8692327e-06
    real(dp), parameter :: kJpermolperang3_to_Pa_x_Pa_to_atm = 16394.078514951_dp
    real(dp) :: deltaN, PMV_correction, Volodymyr_empirical_correction, PMV, cavity_volume, PMV_correction_with_cavity_volume
    real(dp), parameter :: zerodp = 0._dp, epsilon = 1.0e-3_dp
    real(dp) :: deltaG_emptybox, deltaG
    integer :: s, nx, ny, nz, no

    numberdensity = solvent(1)%n0
    deltaGtotMDFT = ff%tot


! TRY weighted density version of cgb
!    if( getinput%char("bridge", defaultvalue="no") == "cgb" ) then
!       write(*,*) 'cgb_version=', cgb_version
!       call energy_and_gradient() ! compute the grandpotential of such system
!       call coarse_grained_bridge(fb_wda)
!       write(*,*) 'FB_wda =', fb_wda
!    end if


!... Volodymyr's partial molar volume correction. See J. Phys. Chem. Lett. 5, 1935-1942 (2014)
    if (size(solvent)/=1) print*, "pressure_correction implemented for 1 solvent species only WARNNING this does not make sense"
    allocate(  density(grid%nx,grid%ny,grid%nz),  source=0._dp)
    call grid%integrate_over_orientations( solvent(1)%xi**2 * solvent(1)%rho0, density)

! COMPUTE CAVITY Volume where rho = 0 --also FREE ENERGY; desactivated, does not work with HSB !
    cavity_volume = zerodp
    do nx = 1, grid%nx
      do ny = 1, grid%ny
        do nz = 1, grid%nz
          if( density(nx,ny,nz)/solvent(1)%n0 < epsilon ) then
             cavity_volume = cavity_volume + grid%dv
!             solvent(1)%xi(:, nx, ny, nz) = zerodp
!           else
!             solvent(1)%xi(:, nx, ny, nz) = 1.0_dp
          end if
         end do
       end do
    end do
!    call energy_and_gradient(deltaG_cavity) ! compute the grandpotential of such system
!    call energy_ideal_and_external (fid, fext)
!    deltaG_cavity = deltaG_cavity - fext
!    write(*,'(A,F12.2,A)') "Cavity free energy", deltaG_cavity," kJ/mol"


    ! number of solvent molecules inside the supercell:
    ! when the supercell has the solute inside:wq


    nmolecules_with_solute = sum_byfft_3d(density)*grid%dv  ! replaces: nmolecules_with_solute = sum(density)*grid%dv    without the rounding error of SUM() but for the price of a FFT.

    deallocate( density )
    ! when the supercell is empty of any perturbation, ie pure solvent
    nmolecules_without_solute = solvent(1)%n0*product(grid%length)
    write(*,*)
    write(*,*) "Nb iterations = ",ff%ieval
    write(*,*)
    write(*,'(A,F12.5,A,F9.2,A)') "Supercell volume", product(grid%length)," Ang³ = ",product(grid%length)/1000.," nm³"
    write(*,'(A,F12.7,A)') "Bulk solvent density", solvent(1)%n0," molecule.Ang⁻³"
    write(*,'(A,F12.2)') "Solvent molecules in the supercell containing the solute   ", nmolecules_with_solute
    write(*,'(A,F12.2)') "Solvent molecules in a supercell without solute", nmolecules_without_solute
    deltaN = nmolecules_without_solute - nmolecules_with_solute
    PMV = deltaN/solvent(1)%n0
    write(*,'(A,F8.2,A,F8.2,A)') "ΔN solvent", deltaN,"    Partial molar volume = ",PMV," Ang³"
    write(*,'(A,E6.1,A,F8.2,A)') "Cavity volume with rho < ",epsilon," = ", cavity_volume," Ang³"


block
    !
    ! Compute the pressure of the bulk solvent.
    ! We use : GrandPotential(homogeneous system) = -PV
    !
    use module_energy_and_gradient, only: ff
    logical :: oldvalue
    character(30) :: bridge_option

    bridge_option = getinput%char("bridge", defaultvalue="no")

    write(*,*) "bridge_option =  ", bridge_option
    if(bridge_option == 'hard_sphere' ) write(*,'(A,F7.3,A)') 'hard sphere radius = ',hs(1)%r,' A'
    if(bridge_option == 'cgb' ) write(*,'(A,A,A,F6.2,A,F8.3,A)') 'cgb_version = ',cgb%version,'    B6 = ',cgb%B6,'  sigma = ',cgb%sigma,' A'
    write(*,'(A,F8.2,A)') "LJ long range correction = ", LJ_energy_correction," kJ/mol"
    write(*,'(A,F12.2,A)') "Functional at min   ", deltaGtotMDFT + LJ_energy_correction," kJ/mol"
    write(*,*)

    if( getinput%char("bridge", defaultvalue="no") == "no" .or. getinput%char("bridge", defaultvalue="no") == "none" ) then

    oldvalue = ff%apply_energy_corrections_due_to_charged_solute
    ff%apply_energy_corrections_due_to_charged_solute = .false. ! Don't include corrections due to the external potential since that's a pure homogeneous solvent problem.
    do s=1,size(solvent) 
      solvent(s)%xi = 0._dp ! set density to 0 == empty the simulation box
    end do
    call energy_and_gradient(deltaG_emptybox) ! compute the grandpotential of such system
    ff%apply_energy_corrections_due_to_charged_solute = oldvalue ! put back the old value
    Pbulk = deltaG_emptybox / (grid%lx * grid%ly * grid%lz) ! Omega[rho=rho_0]=PV ! Pbulk in kJ/mol/Ang^3
    write(*,'(A,F12.2,A)') "Bulk pressure       ", Pbulk*kJpermolperang3_to_Pa_x_Pa_to_atm," atm"
    open(81,file="output/pressure")
    write(81,*) Pbulk
    close(81)
    PMV_correction  = -deltaN/solvent(1)%n0*Pbulk  !correction is -PV where V is excluded Volume
    Volodymyr_empirical_correction =  deltaN*thermo%kbT
    PMV_correction_with_cavity_volume = Pbulk*cavity_volume

    write(*,'(A,F12.2,A)') "PMV correction      ", PMV_correction," kJ/mol"
    write(*,'(A,F12.2,A)') "Pid correction      ", Volodymyr_empirical_correction," kJ/mol"
!    write(*,'(A,F12.2,A)') "PMV correction with cavity volume     ", PMV_correction_with_cavity_volume," kJ/mol"
    open(79,file="output/PMV_correction")
    write(79,*) PMV_correction
    close(79)
    open(80,file="output/Pideal_PMV_correction")
    write(80,*) Volodymyr_empirical_correction
    close(80)

    write(*,'(A,F12.2,A)') "SFE ISc /PC         ", deltaGtotMDFT + PMV_correction," kJ/mol"
!    write(*,'(A,F12.2,A)') "SFE ISc*/PC+        ", deltaGtotMDFT + PMV_correction + Volodymyr_empirical_correction," kJ/mol"
!    write(*,'(A,F12.2,A)') "SFE with cavity correction ", deltaGtotMDFT -deltaG_cavity," kJ/mol"
    end if
end block



    write(*,'(A,F12.2,A)') "ESTIMATED SOLVATION FREE ENERGY: ", deltaGtotMDFT + PMV_correction + LJ_energy_correction,"  kJ/mol"

end subroutine pressure_correction

end module module_pressure_correction

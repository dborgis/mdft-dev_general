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
    
!    use module_energy_cproj_mrso, only: c000
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
    real(dp), parameter :: atm_to_bar = 1.015
    real(dp) :: deltaN, PMV_correction, Volodymyr_empirical_correction, PMV, cavity_volume, PMV_correction_with_cavity_volume
    real(dp), parameter :: zerodp = 0._dp, epsilon = 1.0e-3_dp
    real(dp) :: deltaG_emptybox, deltaG, DeltaP, fB, Deltax, Phib
    integer :: s, nx, ny, nz, no
    logical :: apply_pressure_correction

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

    
    ! when the supercell is empty of any perturbation, ie pure solvent
    nmolecules_without_solute = solvent(1)%n0*product(grid%length)

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
    
    oldvalue = ff%apply_energy_corrections_due_to_charged_solute
    ff%apply_energy_corrections_due_to_charged_solute = .false. ! Don't include corrections due to the external potential since that's a pure homogeneous solvent problem.
    
    bridge_option = getinput%char("bridge", defaultvalue="no")
    if(bridge_option == 'hard_sphere' ) write(*,'(A,F7.3,A)') 'hard sphere radius = ',hs(1)%r,' A'
    if(bridge_option == 'cgb' ) then
    write(*,'(A,A,A,A,F8.3,A,F8.3,A,F8.3,A)') 'cgb_version = ',cgb%version,'  at 4th order','    A3 = ',cgb%A3,'    A4 = ',cgb%A4,'  sigma = ',cgb%sigma,' A'
    
    !cgb%version = getinput%char("cgb_version", defaultvalue="sla")
    if(cgb%version == 'lda') then
    
     fB = 0.0 
     do nz=1,grid%nz
       do ny=1,grid%ny
         do nx=1,grid%nx

            deltax = density( nx, ny, nz)/solvent(1)%n0 - 1.0_dp
              Phib = cgb%A3 * deltax**3 + cgb%A4 * deltax**4
            fB = fB +   Phib 
         end do !nx
      end do !ny
    end do !nz
    fB = fB * thermo%kbT * solvent(1)%n0  * grid%dv
    write(*,'(A,F8.2,A)') "LDA bridge correction = ", fB," kJ/mol  (to be subracted to result below)"  
    end if !end if lda
    end if ! end if cgb
    deallocate( density )
    
    write(*,'(A,F8.2,A)') "LJ long range correction = ", LJ_energy_correction," kJ/mol  (included in result below)"
    write(*,'(A,F12.2,A)') "Functional at min   ", deltaGtotMDFT + LJ_energy_correction," kJ/mol"
    write(*,*)

    write(*,*) 'On vide la boîte pour calculer la pression:'
    do s=1,size(solvent) 
      solvent(s)%xi = 0._dp ! set density to 0 == empty the simulation box
    end do
    call energy_and_gradient(deltaG_emptybox) ! compute the grandpotential of such system
    ff%apply_energy_corrections_due_to_charged_solute = oldvalue ! put back the old value
    Pbulk = deltaG_emptybox / (grid%lx * grid%ly * grid%lz) ! Omega[rho=rho_0]=PV ! Pbulk in kJ/mol/Ang^3
    write(*,'(A,F12.2,A)') "Pressure = ", Pbulk*kJpermolperang3_to_Pa_x_Pa_to_atm," atm"
    
    if( getinput%char("bridge", defaultvalue="no") == "no" .or. getinput%char("bridge", defaultvalue="no") == "none" ) then
     
    DeltaP = Pbulk - solvent(1)%pressure/kJpermolperang3_to_Pa_x_Pa_to_atm/atm_to_bar  
    write(*,'(A,F12.2,A)') "HNC pressure - actual pressure = ", DeltaP*kJpermolperang3_to_Pa_x_Pa_to_atm," atm"
    PMV_correction  = -deltaN/solvent(1)%n0*DeltaP  !correction is -DeltaP*V where V is excluded Volume
    write(*,'(A,F12.2,A)') "PMV correction      ", PMV_correction," kJ/mol"
    write(*,'(A,F12.2,A)') "Cavity-volume correction      ", -DeltaP*cavity_volume," kJ/mol"
   
    apply_pressure_correction = getinput%log("apply_pressure_correction", defaultvalue= .false.)
    if( apply_pressure_correction ) then
     
    write(*,'(A,F12.2,A)') "ESTIMATED SOLVATION FREE ENERGY including PMV correction: ", deltaGtotMDFT + PMV_correction + LJ_energy_correction,"  kJ/mol"
    
!    PMV_correction_with_cavity_volume = DeltaP*cavity_volume
!    write(*,'(A,F12.2,A)') "PMV correction with cavity volume     ", PMV_correction_with_cavity_volume," kJ/mol"
     
    else 
      PMV_correction = 0.0_dp
      write(*,'(A,F12.2,A)') "ESTIMATED SOLVATION FREE ENERGY without PV correction: ", deltaGtotMDFT +  LJ_energy_correction,"  kJ/mol"
    end if
    
    end if ! end bridge = none
end block


    

end subroutine pressure_correction

end module module_pressure_correction

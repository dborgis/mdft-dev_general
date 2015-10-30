module module_postprocessing
    implicit none
    private
    public :: init_postprocessing
contains
!
!     !
!     !   Post-processing of MDFT, following user's requirements
!     !
    subroutine init_postprocessing
        use precision_kinds, only: dp
        use module_solvent, only: solvent
        use module_grid, only: grid
        implicit none
        character(len=80) :: filename
        real(dp), allocatable :: density(:,:,:)
        integer :: nx, ny, nz, ix, iy, iz

        nx=grid%nx
        ny=grid%ny
        nz=grid%nz


        !
        ! print density
        !
        allocate (density(nx,ny,nz))
        call grid%integrate_over_orientations( solvent(1)%density, density)
        filename = "output/density.cube"
        call write_to_cube_file (density ,filename)

!         use system,             ONLY: thermocond
!         use module_solvent, only: solvent
!         use module_input,              ONLY: verbose, getinput
!         ! use solute_geometry,    ONLY: soluteIsPlanar => isPlanar, soluteIsLinear => isLinear
!         use constants,          ONLY: zerodp
!         use hardspheres,        only: hs
!         use module_grid, only: grid
!
!         IMPLICIT NONE
!
!         CHARACTER(50):: filename
!         REAL(dp), ALLOCATABLE , DIMENSION (:,:,:,:) :: neq, Px, Py, Pz ! equilibrium density, ie rho(r), and Pi polarization(r)
!         INTEGER(i2b) :: nfft1, nfft2, nfft3
!         INTEGER(i2b) :: s
!
!         nfft1 = grid%nx
!         nfft2 = grid%ny
!         nfft3 = grid%nz
!
!         CALL print_cg_vect_new ! print output/density.bin that contains cg_vect_new
!
!         allocate ( neq (nfft1,nfft2,nfft3,solvent(1)%nspec) ,SOURCE=zerodp)
!         do s=1,size(solvent)
!             call get_final_density (neq,s)
!         end do
!
!         DO s=1,solvent(1)%nspec
!           write(*,'(A,F12.2)') "Solvent molecules in supercell", SUM(neq)*grid%dv *solvent(s)%n0
!         END DO
!
!
!         IF (verbose) THEN
!             filename = 'output/density.cube'
!             CALL write_to_cube_file (neq(:,:,:,1), filename) ! TODO for now only write for the first species
!             IF ( getinput%char("polarization", defaultvalue="no") /= "no" ) THEN
!                 ALLOCATE ( Px (nfft1,nfft2,nfft3,solvent(1)%nspec) ,SOURCE=zerodp)
!                 ALLOCATE ( Py (nfft1,nfft2,nfft3,solvent(1)%nspec) ,SOURCE=zerodp)
!                 ALLOCATE ( Pz (nfft1,nfft2,nfft3,solvent(1)%nspec) ,SOURCE=zerodp)
!                 CALL get_final_polarization (Px,Py,Pz)
!                 filename='output/normP.cube' ; CALL write_to_cube_file ( (SQRT(Px(:,:,:,1)**2+Py(:,:,:,1)**2+Pz(:,:,:,1)**2)), filename) ! TODO for now only write for the first species
!                 filename='output/Px.cube' ; CALL write_to_cube_file(Px(:,:,:,1), filename)
!                 filename='output/z_Px.dat'; CALL compute_z_density(Px(:,:,:,1) , filename)
!                 DEALLOCATE (Px)
!                 filename='output/Py.cube' ; CALL write_to_cube_file(Py(:,:,:,1), filename)
!                 filename='output/z_Py.dat'; CALL compute_z_density(Py(:,:,:,1) , filename)
!                 DEALLOCATE (Py)
!                 filename='output/Pz.cube' ; CALL write_to_cube_file(Pz(:,:,:,1), filename)
!                 filename='output/z_Pz.dat'; CALL compute_z_density(Pz(:,:,:,1) , filename)
!                 DEALLOCATE (Pz)
!             END IF
!
!
!             ! If calculation is for hard sphere fluid in presence of a hard wall compute profile perp wall
!             ! TODO: DONT HAVE TIME TO WRITE THE TEST TODAY
!             filename = 'output/z_density.out'
!             CALL compute_z_density ( neq (:,:,:,1) , filename ) ! TODO for now only write for the first species
!
!             ! IF (soluteIsLinear()) THEN
!             !     ! nothing for now
!             ! END IF
!             !
!             ! IF( soluteIsPlanar() ) THEN
!             !     PRINT*,'This solute has planar symetry'
!             !     filename = 'output/planardensity.out'
!             !     CALL compute_planar_density ( neq (:,:,:,1) , filename ) ! TODO for now only write for the first species
!             ! END IF
!
!             IF ( getinput%char('other_predefined_vext')=='vextdef0' ) THEN
!                 filename = 'output/molecular_density_in_xy_plane.out'
!                 PRINT*,"I am writing file ",filename
!                 OPEN(378,FILE=filename)
!                     BLOCK
!                         INTEGER(i2b) :: i,j
!                         DO i=1,SIZE(neq,1)
!                             DO j=1,SIZE(neq,2)
!                                 WRITE(378,*)[i,j]*grid%length(1:2)/grid%n_nodes(1:2),neq(i,j,1,1)
!                             END DO
!                             WRITE(378,*)
!                         END DO
!                     END BLOCK
!                 CLOSE(378)
!             END IF
!         END IF
!
!
        filename = 'output/rdf.out'
        call output_rdf ( density/solvent(1)%n0 , filename ) ! Get radial distribution functions
        call output_gsitesite
        call output_gOfRandCosTheta
        deallocate (density)
!
!         CALL adhoc_corrections_to_gsolv
!
!
!         write(*,'(A,F7.2,A)') "T       ", thermocond%T,    " K"
!         write(*,'(A,F7.2,A)') "kT      ", thermocond%kbT,  " kJ/mol"
!         write(*,'(A,F7.2,A)') "β=(kT)⁻¹", thermocond%beta, " (kJ/mol)⁻¹"
!         if( allocated(hs) ) then
!           block
!             real(dp)::x
!             x=hs(1)%pf
!             write(*,'(A,F7.2)') "packing fraction η =",x
!             write(*,'(A,F7.2)') "βP/n PY by pressure route        = (1+2η+3η²)           /(1-η)² =",(1+2*x+3*x**2)/(1-x)**2
!             write(*,'(A,F7.2)') "βP/n PY by compressibility route = (1+ η+ η²)           /(1-η)³ =",(1+x+x**2)/(1-x)**3
!             write(*,'(A,F7.2)') "βP/n CS                          = (1+ η+ η²-η³)        /(1-η)³ =",(1+x+x**2-x**3)/(1-x)**3
!             write(*,'(A,F7.2)') "βP/n CSK                         = (1+ η+ η²-2(1+η)η³/3)/(1-η)³ =",&
!               (1+x+x**2-(2./3.)*(1+x)*x**3)/((1.-x)**3)
!           end block
!         end if
!
    end subroutine init_postprocessing
!
!     SUBROUTINE print_cg_vect_new
!         use module_minimizer, ONLY: cg_vect_new
!         if ( .not. allocated ( cg_vect_new ) ) then
!             print *, 'cg_vect_new is not allocated in SUBROUTINE print_cg_vect_new in process_output.f90. STOP.'
!             stop
!         END IF
!         OPEN (10, file = 'output/density.bin.out' , form = 'unformatted' )
!             write ( 10 ) cg_vect_new
!         CLOSE (10)
!     END SUBROUTINE print_cg_vect_new
!
!
!
!     subroutine adhoc_corrections_to_gsolv
!     ! ... Here, we print all the adhoc corrections one should take into account before comparing MDFT results to MD and/or experiments.
!
!         use precision_kinds, only: dp, sp, i2b
!         use system, only: thermocond
!         use module_solute, only: solute
!         use module_solvent, only: solvent
!         use module_minimizer, only: FF , cg_vect_new
!         use constants, only: zerodp
!         use mathematica, only: chop
!         use module_input, only: getinput
!         use module_grid, only: grid
!         implicit none
!
!         real(dp) :: correction,correction2, Pbulk
!         real(dp), allocatable :: neq(:,:,:,:) ! equilibrium density
!         integer :: nfft1, nfft2, nfft3
!         integer :: s, ios
!         type :: nmoleculetype
!             real(dp) :: withsolute
!             real(dp) :: bulk
!         end type nmoleculetype
!         type (nmoleculetype), allocatable :: nmolecule(:)
!         logical :: file_exists
!         real(dp) :: gamma ! quadrupole moment trace
!         real(dp) :: numberdensity ! molecular number density, for instance 0.0332891 molecule per angstrom^3
!         real(dp) :: solutecharge ! net charge of the solute, for instance -1 for Cl- ion
!         real(dp), parameter :: kJpermolperang3_to_Pa = 1.66113*10**9
!         real(dp), parameter :: Pa_to_atm = 9.8692327e-06
!         real(dp) :: FFcorrected_final, deltaN, Pscheme_correction
!
!         FFcorrected_final = FF
!
!         nfft1 = grid%n_nodes(1)
!         nfft2 = grid%n_nodes(2)
!         nfft3 = grid%n_nodes(3)
!
!         open(79,file="output/FF"); write(79,*) FF; close(79)
!
!
!
!         !... We use P-scheme instead of M-scheme for the electrostatics in MDFT.
!         ! See Kastenholz and Hunenberger, JCP 124, 124106 (2006), page 224501-8, equations 35, 35 and 37 with Ri=0
!         ! "To be applied if the solvent molecule is rigid and involves a single van der Waals interaction site M,
!         ! and that any scheme relying on molecular-cutoff truncation refers to this specific site for applying the truncation."
!         Pscheme_correction = 0._dp
!         if( getinput%log("poisson_solver") ) then
!           gamma = solvent(1)%quadrupole(1,1)+solvent(1)%quadrupole(2,2)+solvent(1)%quadrupole(3,3) ! quadrupole moment trace
!           numberdensity = solvent(1)%n0
!           solutecharge = sum(solute%site%q)
!           Pscheme_correction = -gamma*numberdensity*2.909857E3*solutecharge ! in kJ/mol
!         end if
!         write(*,'(A,F12.2,A)') "P-scheme correction ", Pscheme_correction," kJ/mol"
!         open(79,file="output/Pscheme_correction"); write(79,*) Pscheme_correction; close(79)
!         FFcorrected_final = FFcorrected_final + Pscheme_correction
!
!
!
!
!       !... Volodymyr's partial molar volume correction. See J. Phys. Chem. Lett. 5, 1935-1942 (2014)
!       correction = zerodp
!       if (size(solvent)/=1) stop "CRITICAL in adhoc_corrections_to_gsolv. only 1 solvent species is implemented."
!       do s=1,size(solvent)
!           if (.not. allocated( solvent(s)%n )) allocate ( solvent(s)%n(nfft1,nfft2,nfft3) ,source=0._dp)
!           call get_final_density ( solvent(s)%n , s) ! Get the final density(r) from the last minimizer step.
!       end do
!       allocate (nmolecule(size(solvent)))
!       do concurrent (s=1:size(solvent))
!           nmolecule%withsolute = sum(solvent(s)%n * solvent(s)%n0)  *grid%dv ! number of solvent molecules inside the supercell containing the solute
!       end do
!       nmolecule%bulk = solvent%n0*product(grid%length) ! number of solvent molecules inside the same supercell (same volume) without solute.
!       write(*,'(A,F12.2)') "Solvent molecules with solute   ", nmolecule%withsolute
!       write(*,'(A,F12.2)') "Solvent molecules without solute", nmolecule%bulk
!       write(*,'(A,F12.2)') "ΔN solvent", nmolecule(1)%bulk - nmolecule(1)%withsolute
!       write(*,'(A,F12.7,A)') "Solvent density", solvent(1)%n0," molecule.Ang⁻³"
!       write(*,'(A,F12.5,A)') "Supercell volume", product(grid%length)," Ang³"
!
!
!
!       ! pressure of the bulk solvent?  GrandPotential(homogeneous system) = -PV
!       ! grand potential[rho_bulk] == PV
!       cg_vect_new = zerodp ! set Density to 0
!       FF = zerodp         ! set energy to 0
!       call energy_and_gradient(-10) ! this step is not a minimization step so we give a negative integeration number to avoid the printing of the not relevant obtained energies
!       Pbulk = FF/product(grid%length) ! Omega[rho=rho_0]=PV ! Pbulk in kJ/mol/Ang^3
!       write(*,'(A,F12.2,A)') "Bulk pressure       ", Pbulk*kJpermolperang3_to_Pa*Pa_to_atm," atm"
!       open(81,file="output/bulk-pressure"); write(81,*) Pbulk; close(81)
!
!
!
!
!       s = 1
!       if( s /= 1 ) stop "line 61 of adhoc_corr... we have not thought of multi species case"
!       deltaN = nmolecule(s)%bulk - nmolecule(s)%withsolute
!       correction  = -(nmolecule(s)%bulk - nmolecule(s)%withsolute)/solvent(s)%n0*Pbulk  !correction is -PV where V is excluded Volume
!       correction2 =  (nmolecule(s)%bulk - nmolecule(s)%withsolute)*thermoCond%kbT  !correction is -PV where V is excluded Volume
!       FFcorrected_final = FFcorrected_final + correction !+ correction2
!       write(*,'(A,F12.2,A)') "PMV correction      ", correction," kJ/mol"
!       write(*,'(A,F12.2,A)') "Pid correction      ", correction2," kJ/mol"
!       open(79,file="output/PMV_correction")
!         write(79,*) correction
!       close(79)
!       open(80,file="output/Pideal_PMV_correction")
!         write(80,*) correction2
!       close(80)
!
!       write(*,'(A,F12.2,A)') "SFE ISc             ", FFcorrected_final," kJ/mol"
!       write(*,'(A,F12.2,A)') "SFE ISc*            ", FFcorrected_final + correction2," kJ/mol"
!
!     end subroutine adhoc_corrections_to_gsolv
!
!
!
!     pure subroutine get_final_density ( neq , s)
!         use precision_kinds, only: dp
!         use module_solvent, only: solvent
!         use module_quadrature, only: mean_over_orientations
!         implicit none
!         integer, intent(in) :: s ! the solvent species of which we want the number density
!         real(dp), intent(out) :: neq(:,:,:)
!         call mean_over_orientations( solvent(solventspecies)%density , neq)
!     end subroutine get_final_density
!
!
!
!
!
!     SUBROUTINE get_final_polarization ( Px , Py , Pz )
!
!         use precision_kinds, only: dp, i2b
!         use module_solvent, only: solvent
!         use module_grid, only: grid
!
!         IMPLICIT NONE
!         INTEGER(i2b) :: i, j, k, io, s
!         REAL(dp) :: x, local_Px, local_Py, local_Pz
!         REAL(dp), dimension(:,:,:,:), intent(out) :: Px, Py, Pz ! equilibrium polarization(r)
!         integer :: nx, ny, nz, no, ns
!         real(dp), parameter :: zerodp = 0._dp
!
!
!         Px = zerodp
!         Py = zerodp
!         Pz = zerodp
!
!         DO s =1,solvent(1)%nspec
!             DO i =1,grid%nx
!                 DO j =1,grid%ny
!                     DO k =1,grid%nz
!                         local_Px = 0.0_dp
!                         local_Py = 0.0_dp
!                         local_Pz = 0.0_dp
!                         DO io =1,grid%no
!                             x = solvent(s)%density(i,j,k,io)
!                             local_Px = local_Px + grid%omx(io) * grid%w(io) * x
!                             local_Py = local_Py + grid%omy(io) * grid%w(io) * x
!                             local_Pz = local_Pz + grid%omz(io) * grid%w(io) * x
!                         END DO
!                         Px(i,j,k,s) = local_Px
!                         Py(i,j,k,s) = local_Py
!                         Pz(i,j,k,s) = local_Pz
!                     END DO
!                 END DO
!             END DO
!         END DO
!
!     END SUBROUTINE get_final_polarization
!
!
end module

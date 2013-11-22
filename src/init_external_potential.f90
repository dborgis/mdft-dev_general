! This SUBROUTINE computes the external potential. It is one of the most time consuming routine.
!Warning there are two ways of calculating the electrostatic potential (Poisson solver and point charge) you should always have one tag
!T and one tag F for the electrostatic potential, if thera are 2 tags T, it is the last evaluation which counts, i.e Poisson solver. 

SUBROUTINE init_external_potential

    USE precision_kinds, ONLY: dp , i2b
    USE input, ONLY: input_log, input_char
    USE system , ONLY: chg_solv, soluteSite, spaceGrid
    USE external_potential, ONLY: Vext_total, Vext_q
    USE mod_lj, ONLY: initLJ => init
    USE quadrature, ONLY: Rotxx, Rotxy, Rotxz, Rotyx, Rotyy, Rotyz, Rotzx, Rotzy, Rotzz, angGrid, molRotGrid
    USE constants, ONLY: zero
    
    IMPLICIT NONE
    
    integer(i2b):: nb_id_mol , nb_id_solv ! nb of types of sites of solute and solvent
    INTEGER(i2b), DIMENSION(3) :: nfft

    nfft = spaceGrid%n_nodes

    IF( .NOT. ALLOCATED( Vext_total )) THEN
        allocate( Vext_total(spaceGrid%n_nodes(1),spaceGrid%n_nodes(2),spaceGrid%n_nodes(3),angGrid%n_angles,&
                            molRotGrid%n_angles,SIZE(soluteSite)), source=0._dp )
    ELSE
        STOP "see init_external_potential.f90 vext_total is already allocated."
    END IF

    nb_id_mol  = size ( soluteSite  ) ! total number of solute types
    nb_id_solv = size ( chg_solv ) ! total number of solvent types


    ! Hard walls
    call external_potential_hard_walls

    IF (input_log('point_charge_electrostatic')) THEN ! Charges : treatment as point charges
        IF (.NOT. ALLOCATED(Vext_q)) THEN
            BLOCK
                INTEGER(i2b), DIMENSION(6) :: al
                al(1) = nfft(1)
                al(2) = nfft(2)
                al(3) = nfft(3)
                al(4) = angGrid%n_angles
                al(5) = molRotGrid%n_angles
                al(6) = SIZE(soluteSite)
                ALLOCATE ( vext_q (al(1),al(2),al(3),al(4),al(5),al(6)), SOURCE=zero )
            END BLOCK
        END IF
        CALL compute_vcoul_as_sum_of_pointcharges( Rotxx, Rotxy, Rotxz, Rotyx, Rotyy, Rotyz, Rotzx, Rotzy, Rotzz )
    END IF

    ! Charges : Poisson solver
    if (input_log('poisson_solver')) then
        BLOCK
            REAL(dp), DIMENSION (nfft(1),nfft(2),nfft(3)) :: soluteChargeDensity
            if (.not. allocated(Vext_q) ) &
                allocate ( Vext_q ( nfft(1),nfft(2),nfft(3),angGrid%n_angles,molRotGrid%n_angles,SIZE(soluteSite)), SOURCE=zero)
            CALL soluteChargeDensityFromSoluteChargeCoordinates (soluteChargeDensity)
            call poissonSolver (soluteChargeDensity)
            call vext_q_from_v_c (Rotxx,Rotxy,Rotxz,Rotyx,Rotyy,Rotyz,Rotzx,Rotzy,Rotzz)
        END BLOCK
    END IF

    ! Lennard-Jones
    call initLJ

    ! r^-12 only
    if (input_log('purely_repulsive_solute')) then
        call compute_purely_repulsive_potential ( Rotxx , Rotxy , Rotxz , Rotyx , Rotyy , Rotyz , Rotzx , Rotzy , Rotzz )
    END IF

    ! hard spherical solute
    if (input_log('hard_sphere_solute')) then
        call compute_vext_hard_sphere
    END IF
    
    ! potential created by hard wall(s) ! CALLED TWICE BUG
    call external_potential_hard_walls

    ! hard cylinder
    if (input_log('hard_cylinder_solute')) then
        call compute_vext_hard_cylinder
    END IF
    
    ! personnal vext as implemented in personnal_vext.f90
    if (input_log('personnal_vext')) then
        call compute_vext_perso
    END IF
    
    ! compute total Vext(i,j,k,omega), the one used in the free energy functional
    call vext_total_sum
    
END SUBROUTINE init_external_potential

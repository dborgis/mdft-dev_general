! In this SUBROUTINE one calls the different parts of the total energy
! first is computed the radial_part, then ... blabla.
! this SUBROUTINE is the one called by the minimization stuff
! for computing the total energy and associated gradient
! FF is the TOTAL ENERGY of the system, it is thus the functional of the density that is minimized by solver
! dF is the gradient of FF with respect to all coordinates. Remember it is of the kind dF ( number of variables over density (ie angles etc))

SUBROUTINE energy_and_gradient (iter)

    USE precision_kinds, ONLY: i2b , dp
    USE input, ONLY: input_log, input_char, verbose
    USE minimizer, ONLY: FF , dF
    
    IMPLICIT NONE
    
    INTEGER(i2b), INTENT(INOUT) :: iter
    REAL(dp) :: Fext,Fid,FexcRad,FexcPol,F3B1,F3B2,Ffmt
    Fext = 0._dp
    Fid = 0._dp
    Ffmt = 0._dp
    FexcRad = 0._dp
    FexcPol = 0._dp
    F3B1 = 0._dp
    F3B2 = 0._dp
    Ffmt = 0._dp
    

    CALL init_functional_and_gradient_to_zero( FF, dF )

    CALL energy_external (Fext)
    CALL energy_ideal (Fid)

    ! compute radial part of the excess free energy
    IF (input_log('hard_sphere_fluid') ) CALL energy_hard_sphere_fmt ! => pure hard sphere contribution

    ! test if there is a LJ perturbation to hard spheres ! WCA model etc. to implement more intelligently
    IF (input_log('lennard_jones_perturbation_to_hard_spheres') ) CALL lennard_jones_perturbation_to_hard_spheres ! lennard jones perturbative contribution => Weeks-Chandler-Anderson

    IF (input_log('read_ck_or_chi')) THEN
        IF (input_log('hydrophobicity')) THEN
            IF (TRIM(ADJUSTL(input_char('treatment_of_hydro')))=='C')  THEN
                CALL cs_plus_hydro (FexcRad)
            ELSE IF (TRIM(ADJUSTL(input_char('treatment_of_hydro')))=='VdW')  THEN
                CALL energy_hydro (FexcRad)
            END IF
        ELSE
            CALL cs_from_dcf (FexcRad)
        END IF
    END IF

    ! bridge calculation: F(FMT)-F(c2hs)+F(c2H2O)
    IF (input_log('bridge_hard_sphere')) CALL energy_cs_hard_sphere ! better name should be given

    ! Dipolar polarization. What user wants (use it or not) is checked in SUBROUTINE for clearer code.
    CALL energy_polarization (FexcPol)
    CALL energy_polarization_myway (FexcPol)

    ! Threebody term that is needed to empiricaly force H-bonding in water. What user wants (use it or not) is checked in SUBROUTINE for clearer code
    CALL energy_threebody
    CALL energy_threebody_faster

    IF (verbose) THEN
        WRITE(*,*)'SOLVATION FREE ENERGY AT THIS STEP = ',FF
        WRITE(*,*)'-----------------------'
    END IF

    WRITE(*,'(  i3,f11.3,f11.3,f11.3,f11.3,f11.3,f11.3,f11.3,f11.3)')iter,FF,norm2(dF),Fext,Fid,FexcRad,FexcPol,F3B1,F3B2

    CONTAINS
    
    PURE SUBROUTINE init_functional_and_gradient_to_zero (FF,dF)
        REAL(dp), INTENT(OUT) :: FF, dF(:)
        REAL(dp), PARAMETER :: zero=0._dp
        FF = zero ; dF = zero ! functional and gradient
    END SUBROUTINE init_functional_and_gradient_to_zero
    
END SUBROUTINE energy_and_gradient

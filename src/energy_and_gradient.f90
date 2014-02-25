! In this SUBROUTINE one calls the different parts of the total energy
! first is computed the radial_part, then ... blabla.
! this SUBROUTINE is the one called by the minimization stuff
! for computing the total energy and associated gradient
! FF is the TOTAL ENERGY of the system, it is thus the functional of the density that is minimized by solver
! dF is the gradient of FF with respect to all coordinates. Remember it is of the kind dF ( number of variables over density (ie angles etc))
SUBROUTINE energy_and_gradient (iter)

    USE precision_kinds, ONLY: i2b , dp
    USE input, ONLY: input_log, input_char, input_dp
    USE minimizer, ONLY: FF , dF
    USE system, ONLY : lambda1_mol, lambda2_mol
    IMPLICIT NONE
    
    INTEGER(i2b), INTENT(INOUT) :: iter
    REAL(dp) :: Fext,Fid,Fexcnn,FexcPol,F3B1,F3B2,Ffmt,Ffmtcs
    LOGICAL :: opn
    Fext = 0._dp
    Fid = 0._dp
    Ffmt = 0._dp
    Ffmtcs = 0._dp
    Fexcnn = 0._dp
    FexcPol = 0._dp
    F3B1 = 0._dp
    F3B2 = 0._dp
    

    CALL init_functional_and_gradient_to_zero( FF, dF )

    CALL energy_external (Fext)
    CALL energy_ideal (Fid)

    ! compute radial part of the excess free energy
    IF (input_log('hard_sphere_fluid') ) CALL energy_hard_sphere_fmt (Ffmt) ! => pure hard sphere contribution

    ! test if there is a LJ perturbation to hard spheres ! WCA model etc. to implement more intelligently
    IF (input_log('lennard_jones_perturbation_to_hard_spheres') ) CALL lennard_jones_perturbation_to_hard_spheres ! lennard jones perturbative contribution => Weeks-Chandler-Anderson

    IF (input_log('readDensityDensityCorrelationFunction')) THEN
        IF (input_log('hydrophobicity')) THEN
            IF (input_char('treatment_of_hydro')=='C')  THEN
                CALL energy_nn_cs_plus_nbar (Fexcnn)
            ELSE IF (input_char('treatment_of_hydro')=='VdW')  THEN
                CALL energy_hydro (Fexcnn)
            END IF
        ELSE
            IF   (.NOT. input_log('include_nc_coupling')) THEN
                 CALL energy_nn_cs (Fexcnn)
            END IF
        END IF
    END IF

    ! bridge calculation: F(FMT)-F(c2hs)+F(c2H2O)
    IF (input_log('bridge_hard_sphere') .AND. .NOT. input_log('hard_sphere_fluid')) THEN
        STOP 'bridge_hard_sphere and hard_sphere_fluid should be both turned TRUE for a calculation with bridge'
    END IF
    IF (input_log('bridge_hard_sphere')) CALL energy_cs_hard_sphere (Ffmtcs) ! better name should be given


    IF ( input_log('polarization') ) THEN
        IF ( input_char('polarization_order')=='dipol' ) THEN ! cs cdelta cd ( polarization_order = dipol )
            CALL energy_polarization_dipol (FexcPol)
        ELSE IF ( input_char('polarization_order')=='multi' .AND. (.NOT. input_log('include_nc_coupling')) ) THEN ! ( polarization_order = multi )
            CALL energy_polarization_multi (FexcPol)
        ELSE IF ( input_char('polarization_order')=='multi' .AND. ( input_log('include_nc_coupling')) ) THEN ! ( polarization_order = multi )
            CALL energy_polarization_multi_with_nccoupling(FexcPol)
        ELSE
        STOP "You want to include polarization but the tag for polarization order is neither dipol nor multi"
        END IF
    END IF


    IF ( input_log('threebody') ) THEN
        IF (SUM(ABS(lambda1_mol(:))+ABS(lambda2_mol(:)))/=0.0_dp .OR. input_dp('lambda_solvent')/=0.0_dp) THEN
            CALL energy_threebody_faster (F3B1, F3B2)
        END IF
    END IF

    WRITE(*,'(  i3,f11.3,f11.3,f11.3,f11.3,f11.3,f11.3,f11.3,f11.3,f11.3,f11.3)')&
        iter,FF,norm2(dF),Fext,Fid,Fexcnn,FexcPol,F3B1,F3B2,Ffmt,Ffmtcs

    INQUIRE(FILE='output/iterate.dat', OPENED = opn)
    IF (.not. opn) THEN
        OPEN(111, file='output/iterate.dat')
        WRITE(111,*) "# ",'FF  ','|dF   |','  Fext','   Fid','    Fexc/rad','    Fexc/pol','    F3B1','      F3B2','    Ffmt',&
        '    Ffmtcs'
        WRITE(111,*) '#########################################################################'
        WRITE(111, *) iter,FF,norm2(dF),Fext,Fid,Fexcnn,FexcPol,F3B1,F3B2,Ffmt,Ffmtcs
   ELSE
        WRITE(111, *) iter,FF,norm2(dF),Fext,Fid,Fexcnn,FexcPol,F3B1,F3B2,Ffmt,Ffmtcs  
   END IF  
    
   
    CONTAINS
    
    PURE SUBROUTINE init_functional_and_gradient_to_zero (FF,dF)
        REAL(dp), INTENT(OUT) :: FF, dF(:)
        REAL(dp), PARAMETER :: zero=0._dp
        FF = zero ; dF = zero ! functional and gradient
    END SUBROUTINE init_functional_and_gradient_to_zero
    
END SUBROUTINE energy_and_gradient

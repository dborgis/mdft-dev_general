!===================================================================================================================================
subroutine energy_and_gradient (iter)
!===================================================================================================================================
! In this SUBROUTINE one calls the different parts of the total energy
! first is computed the radial_part, then ... blabla.
! this SUBROUTINE is the one called by the minimization stuff
! for computing the total energy and associated gradient
! FF is the TOTAL ENERGY of the system, it is thus the functional of the density that is minimized by solver
! dF is the gradient of FF with respect to all coordinates. Remember it is of the kind dF ( number of variables over density (ie angles etc))

use precision_kinds ,only: i2b, dp
use input           ,only: input_log, input_char, input_dp
use minimizer       ,only: FF, dF
use system          ,only: solute
use dcf             ,only: c_s, c_s_hs

implicit none

    INTEGER(i2b), INTENT(INOUT) :: iter
    REAL(dp) :: Fext,Fid,Fexcnn,FexcPol,F3B1,F3B2,F3B_ww,Ffmt,Ffmtcs,Fexc_ck_angular, dF_loc(size(dF))
    LOGICAL :: opn
    INTEGER(i2b) :: karim, exitstatus

    F3B_ww=0.0_dp
    Fext = 0._dp
    Fid = 0._dp
    Ffmt = 0._dp
    Ffmtcs = 0._dp
    Fexcnn = 0._dp
    FexcPol = 0._dp
    F3B1 = 0._dp
    F3B2 = 0._dp
    Fexc_ck_angular = 0._dp

    FF=0
    dF=0

    IF (iter==1) PRINT*,

    CALL energy_external (Fext)
    CALL energy_ideal (Fid)

    ! compute Fexc_ck_angular
    karim = 0
    IF ( input_log("ck_angular") ) karim = karim + 1
    IF ( input_log("ck_debug") ) karim = karim + 1
    IF ( input_log("ck_debug_extended") ) karim = karim + 1
    IF (karim == 1) THEN
        IF ( input_log("ck_en_harm_sph") ) STOP 'The ck_angular methods and ck_en_harm_sph are exclusive.'
        CALL energy_ck_angular (Fexc_ck_angular)
    ELSE IF (karim /= 0) THEN
        STOP 'The ck methods - ck_angular, ck_debug, ck_debug_extended are exclusive.'
    END IF

    ! compute Fexc_ck_angular en harmoniques spheres
    IF (input_log("ck_en_harm_sph") ) CALL energy_ck_en_harm_sph (Fexc_ck_angular)

    ! compute radial part of the excess free energy
    if (input_log('hard_sphere_fluid') ) then
      call energy_hard_sphere_fmt (Ffmt) ! => pure hard sphere contribution
      FF=FF+Ffmt
    end if

    ! test if there is a LJ perturbation to hard spheres ! WCA model etc. to implement more intelligently
    IF (input_log('lennard_jones_perturbation_to_hard_spheres') ) CALL lennard_jones_perturbation_to_hard_spheres ! lennard jones perturbative contribution => Weeks-Chandler-Anderson
    ! NOTE THAT THIS SHOULD USE THE NEW energy_cs.f90  instead of rewritting again and again the same things.

    IF (input_log('readDensityDensityCorrelationFunction')) THEN
        IF (input_log('hydrophobicity')) THEN

            IF (input_char('treatment_of_hydro')=='C')  THEN
                CALL energy_nn_cs_plus_nbar (Fexcnn)
            ELSE IF (input_char('treatment_of_hydro')=='VdW')  THEN
                IF (input_log('bridge_hard_sphere') ) THEN
                    PRINT*, 'You are using HSB and VdW so you are withdrawing twice the HS second order term'
                    STOP
                END IF
                CALL energy_hydro (Fexcnn)
            ELSE
                STOP "Hydrophobicity TRUE can only be associated to treatment_of_hydro == C or VdW"
            END IF

        ELSE

            if( .not. input_log('include_nc_coupling') ) then
              call energy_cs( c_s, Fexcnn, dF_loc, exitstatus)
              if( exitstatus /= 1 ) stop "problem in subroutine energy_cs( c_s, Fexcnn, dF_loc, exitstatus)"
              FF=FF+Fexcnn
              dF=dF+dF_loc
            end if

        END IF
    END IF

    ! bridge calculation: F(FMT)-F(c2hs)+F(c2H2O)
    IF (input_log('bridge_hard_sphere') .AND. .NOT. input_log('hard_sphere_fluid')) THEN
        STOP 'bridge_hard_sphere and hard_sphere_fluid should be both turned TRUE for a calculation with bridge'
    END IF
    ! if (input_log('bridge_hard_sphere')) then
    !   call energy_cs( c_s_hs, Ffmtcs, dF_loc, exitstatus)
    !   if( exitstatus /= 1 ) stop "problem in subroutine energy_cs( c_s_hs, Ffmtcs, dF_loc, exitstatus)"
    !   FF=FF-Ffmtcs
    !   dF=dF-dF_loc
    ! end if


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
        IF (SUM(ABS(solute%site%lambda1)+ABS(solute%site%lambda1))/=0.0_dp .OR. input_dp('lambda_solvent')/=0.0_dp) THEN
            CALL energy_threebody_faster (F3B1, F3B2, F3B_ww)
        END IF
    END IF

!   WRITE(*,'(  i3,f11.3,f11.3,f11.3,f11.3,f11.3,f11.3,f11.3,f11.3,f11.3,f11.3)')&
!   iter,FF,norm2(dF),Fext,Fid,Fexcnn,FexcPol,F3B1,F3B2,Ffmt,Ffmtcs
IF (iter /= -10) THEN  !!!Do not write it for the step that is used to compute ad_hoc correction, i.e rho=rho_O
    write(*,'(A,I3)') "iteration : ",iter
    write(*,'(A)')"======================================"
    write(*,'(A,F16.8)') "FF              =",FF
    write(*,'(A,F16.8)') "norm2(dF)       =",norm2(dF)
    write(*,'(A,F16.8)') "Fext            =",Fext
    write(*,'(A,F16.8)') "Fid             =",Fid
    write(*,'(A,F16.8)') "Fexcnn          =",Fexcnn
    write(*,'(A,F16.8)') "FexcPol         =",FexcPol
    write(*,'(A,F16.8)') "F3B1            =",F3B1
    write(*,'(A,F16.8)') "F3B2            =",F3B2
    write(*,'(A,F16.8)') "F3B_solvent     =",F3B_ww
    write(*,'(A,F16.8)') "Ffmt            =",Ffmt
    ! write(*,'(A,F16.8)') "Ffmtcs          =",Ffmtcs
    write(*,'(A,F16.8)') "Fexc_ck_angular =",Fexc_ck_angular
    write(*,'(A)')
END IF

    INQUIRE(FILE='output/iterate.dat', OPENED = opn)
    IF (.not. opn) THEN
        OPEN(111, file='output/iterate.dat')
        WRITE(111,*) "# ",'FF  ','|dF   |','  Fext','   Fid','    Fexc/rad','    Fexc/pol','    F3B1','      F3B2','    Ffmt',&
        '    Ffmtcs', '  Fexc_ck_angular'
        WRITE(111,*) '#########################################################################'
    END IF
    WRITE(111, *) iter,FF,norm2(dF),Fext,Fid,Fexcnn,FexcPol,F3B1,F3B2,Ffmt,Ffmtcs,Fexc_ck_angular

END SUBROUTINE energy_and_gradient
!===================================================================================================================================

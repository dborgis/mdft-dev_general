!===================================================================================================================================
subroutine energy_and_gradient (iter)

! In this SUBROUTINE one calls the different parts of the total energy
! first is computed the radial_part, then ... blabla.
! this SUBROUTINE is the one called by the minimization stuff
! for computing the total energy and associated gradient
! FF is the TOTAL ENERGY of the system, it is thus the functional of the density that is minimized by solver
! dF_new is the gradient of FF with respect to all coordinates. Remember it is of the kind dF_new ( number of variables over density (ie angles etc))

    use precision_kinds ,only: dp
    use module_input           ,only: getinput
    use module_minimizer       ,only: FF, dF_new
    use system          ,only: solute, solvent
    use dcf             ,only: c_s!, c_s_hs
    use module_grid     ,only: grid

    use ENERGY, ONLY : energy_cs

    implicit none

    INTEGER, INTENT(INOUT) :: iter
    REAL(dp) :: Fext,Fid,Fexcnn,FexcPol,F3B1,F3B2,F3B_ww,Ffmt,Ffmtcs,Fexc_ck_angular
    real(dp) :: df_loc (grid%nx, grid%ny, grid%nz, grid%no, solvent(1)%nspec)
    LOGICAL :: opn
    INTEGER :: i, exitstatus

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
    dF_new=0

    IF (iter==1) PRINT*,

    CALL energy_ideal (Fid)
    CALL energy_external (Fext)


    !
    ! ! compute Fexc_ck_angular
    ! i=0
    ! IF( getinput%log( "ck_angular"        , defaultvalue=.false.) ) i=i+1
    ! IF( getinput%log( "ck_debug"          , defaultvalue=.false.) ) i=i+1
    ! IF( getinput%log( "ck_debug_extended" , defaultvalue=.false.) ) i=i+1
    ! IF( i==1 ) THEN
    !     IF( getinput%log("ck_en_harm_sph", defaultvalue=.false.) ) STOP 'The ck_angular methods and ck_en_harm_sph are exclusive.'
    !     CALL energy_ck_angular (Fexc_ck_angular)
    ! ELSE IF (i>2) THEN
    !     STOP 'The ck methods - ck_angular, ck_debug, ck_debug_extended are exclusive.'
    ! END IF
    !
    ! ! compute Fexc_ck_angular en harmoniques spheres
    ! IF (getinput%log("ck_en_harm_sph", defaultvalue=.false.) ) CALL energy_ck_en_harm_sph (Fexc_ck_angular)

    ! compute radial part of the excess free energy
    if (getinput%log('hard_sphere_fluid', defaultvalue=.false.) ) then
      call energy_fmt (Ffmt) ! => pure hard sphere contribution
      FF=FF+Ffmt
    end if

    ! test if there is a LJ perturbation to hard spheres ! WCA model etc. to implement more intelligently
    IF (getinput%log('lennard_jones_perturbation_to_hard_spheres', defaultvalue=.false.) ) CALL lennard_jones_perturbation_to_hard_spheres ! lennard jones perturbative contribution => Weeks-Chandler-Anderson
    ! NOTE THAT THIS SHOULD use THE NEW energy_cs.f90  instead of rewritting again and again the same things.

    IF (getinput%log('readDensityDensityCorrelationFunction', defaultvalue=.true.)) THEN
        IF (getinput%log('hydrophobicity', defaultvalue=.false.)) THEN

            IF (getinput%char('treatment_of_hydro')=='C')  THEN
                CALL energy_nn_cs_plus_nbar (Fexcnn)
            ELSE IF (getinput%char('treatment_of_hydro')=='VdW')  THEN
                IF (getinput%log('bridge_hard_sphere', defaultvalue=.false.) ) THEN
                    PRINT*, 'You are using HSB and VdW so you are withdrawing twice the HS second order term'
                    STOP
                END IF
                CALL energy_hydro (Fexcnn)
            ELSE
                STOP "Hydrophobicity TRUE can only be associated to treatment_of_hydro == C or VdW"
            END IF

        ELSE

            if( .not. getinput%log('include_nc_coupling', defaultvalue=.false.) ) then
              call energy_cs( c_s, Fexcnn, dF_loc, exitstatus)
              if( exitstatus /= 1 ) stop "problem in subroutine energy_cs( c_s, Fexcnn, dF_loc, exitstatus)"
              FF=FF+Fexcnn
              dF_new = dF_new +dF_loc
            end if

        END IF
    END IF

    ! bridge calculation: F(FMT)-F(c2hs)+F(c2H2O)
    IF (getinput%log('bridge_hard_sphere', defaultvalue=.false.) .AND. .NOT. getinput%log('hard_sphere_fluid',defaultvalue=.false.)) THEN
        STOP 'bridge_hard_sphere and hard_sphere_fluid should be both turned TRUE for a calculation with bridge'
    END IF
    ! if (getinput%log('bridge_hard_sphere')) then
    !   call energy_cs( c_s_hs, Ffmtcs, dF_loc, exitstatus)
    !   if( exitstatus /= 1 ) stop "problem in subroutine energy_cs( c_s_hs, Ffmtcs, dF_loc, exitstatus)"
    !   FF=FF-Ffmtcs
    !   dF_new = dF_new -dF_loc
    ! end if


    block
        character(180) :: polarization
        polarization = getinput%char("polarization", defaultvalue="no")
        select case (polarization)
        case("no","none")
        case("dipolar")
            call energy_polarization_dipol (FexcPol)
        case("multipolar_without_coupling_to_density")
            call energy_polarization_multi (FexcPol)
        case("multipolar_with_coupling_to_density")
            call energy_polarization_multi_with_nccoupling (FexcPol)
        case default
            print*, "The tag 'polarization' in input reads ", polarization,". This is not correct"
            stop "in energy_and_gradient"
        end select
    end block





    IF ( getinput%log( 'threebody', defaultvalue=.false. ) ) THEN
        IF (SUM(ABS(solute%site%lambda1)+ABS(solute%site%lambda1))/=0.0_dp .OR. getinput%dp('lambda_solvent')/=0.0_dp) THEN
            CALL energy_threebody_faster (F3B1, F3B2, F3B_ww)
        END IF
    END IF

!   WRITE(*,'(  i3,f11.3,f11.3,f11.3,f11.3,f11.3,f11.3,f11.3,f11.3,f11.3,f11.3)')&
!   iter,FF,norm2(dF_new),Fext,Fid,Fexcnn,FexcPol,F3B1,F3B2,Ffmt,Ffmtcs
IF (iter /= -10) THEN  !!!Do not write it for the step that is used to compute ad_hoc correction, i.e rho=rho_O
    write(*,'(A)')"======================================"
    write(*,'(A,I3)')    "iteration       :",iter
    write(*,'(A,F16.8)') "FF              =",FF
    write(*,'(A,F16.8)') "norm2(dF_new)   =",norm2(dF_new)
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
        WRITE(111,*) "# ",'FF  ','|dF_new   |','  Fext','   Fid','    Fexc/rad','    Fexc/pol','    F3B1','      F3B2','    Ffmt',&
        '    Ffmtcs', '  Fexc_ck_angular'
        WRITE(111,*) '#########################################################################'
    END IF
    WRITE(111, *) iter,FF,norm2(dF_new),Fext,Fid,Fexcnn,FexcPol,F3B1,F3B2,Ffmt,Ffmtcs,Fexc_ck_angular

END SUBROUTINE energy_and_gradient

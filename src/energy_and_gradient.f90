!===================================================================================================================================
subroutine energy_and_gradient (iter, f, df)

! In this SUBROUTINE one calls the different parts of the total energy
! first is computed the radial_part, then ... blabla.
! this SUBROUTINE is the one called by the minimization stuff
! for computing the total energy and associated gradient
! FF is the TOTAL ENERGY of the system, it is thus the functional of the density that is minimized by solver
! dF_new is the gradient of FF with respect to all coordinates. Remember it is of the kind dF_new ( number of variables over density (ie angles etc))

    use precision_kinds, only: dp
    use module_input, only: getinput
    use module_solute, only: solute
    use module_solvent, only: solvent
    use dcf, only: c_s!, c_s_hs
    use module_grid, only: grid
    use energy, only : energy_cs

    implicit none

    real(dp), intent(out) :: f
    real(dp), intent(out) :: df (grid%nx, grid%ny, grid%nz, grid%no, solvent(1)%nspec)
    real(dp) :: df_bis (grid%nx, grid%ny, grid%nz, grid%no, solvent(1)%nspec)
    integer, intent(in) :: iter
    logical :: opn
    integer :: i, exitstatus
    type exc_type
        real(dp) :: cs, pol, 3b1, 3b2, 3bww, ck_angular, wca, fmt
        logical :: do_cs, do_pol, do_3b1, do_3b2, do_3bww, do_ck_angular, do_wca, do_fmt
    end type
    type ff_type
        real(dp) :: id, ext
        logical :: do_id, do_ext
        type(exc_type) :: exc
    end type
    type (ff_type) :: ff





STOP "I AM NOW IN ENERGY AND GRADIENT. I STOP HERE BECAUSE WE ARE TESTING THE INITIALIZATION"



    !   First I prepare the tree of calculations to be done
    ff%do_id = .true.
    ff%do_ext = .true.
    ff%do_fmt = .false.
    ff%do_fmtcs = .false.
    ff%do_excnn = .false.
    ff%do_excpol = .false.
    ff%do_3b1 = .false.
    ff%do_3b2 = .false.


    f=0._dp
    df=0

    IF (iter==1) PRINT*,


    !   We always compute ideal and external contributions
    ff%do_id = .true.
    ff%do_ext = .true.
    ff%exc%do_fmt = getinput%log ('hard_sphere_fluid', defaultvalue=.false.)
    ff%exc%do_wca = getinput%log ('wca', defaultvalue=.false.)
    ff%exc%do_3b = getinput%log ('threebody', defaultvalue=.false. )



    if (ff%do_id) call energy_ideal (ff%id)
    if (ff%do_ext) call energy_external (ff%ext)
    if (ff%do_fmt) call energy_fmt (ff%exc%fmt, df)
    if (ff%exc%do_wca) call lennard_jones_perturbation_to_hard_spheres (ff%exc%wca)



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
        call energy_fmt (f, df) !   pure hard sphere contribution
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
                CALL energy_hydro (fexcnn)
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
    write(*,'(A,F16.8)') "norm2(dF_new)   =",norm2(dF)
    write(*,'(A,F16.8)') "Fext            =",ff%ext
    write(*,'(A,F16.8)') "Fid             =",ff%id
    write(*,'(A,F16.8)') "Fexcnn          =",ff%excnn
    write(*,'(A,F16.8)') "FexcPol         =",ff%excPol
    write(*,'(A,F16.8)') "F3B1            =",ff%3B1
    write(*,'(A,F16.8)') "F3B2            =",ff%3B2
    write(*,'(A,F16.8)') "F3B_solvent     =",ff%3B_ww
    write(*,'(A,F16.8)') "Ffmt            =",ff%fmt
    write(*,'(A,F16.8)') "Ffmtcs          =",ff%fmtcs
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

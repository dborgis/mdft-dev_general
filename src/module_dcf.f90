MODULE dcf

    USE precision_kinds, ONLY: dp, i2b
    USE input, ONLY: input_log, input_char, n_linesInFile, deltaAbscissa

    IMPLICIT NONE

    REAL(dp) :: delta_k , delta_k_in_C! distance between two k points in cs.in, cdelta.in, cd.in
    INTEGER(i2b) :: nb_k, nb_k_in_c ! nb of k points in cs.in, cdelta.in, cd.in
    REAL(dp), ALLOCATABLE, DIMENSION (:) :: c_s ! density density correlation function
    REAL(dp), ALLOCATABLE, DIMENSION (:) :: c_delta, c_d ! polarization polarization correlation function
    REAL(dp), ALLOCATABLE, DIMENSION (:) :: chi_l, chi_t ! longitudinal and transverse dielectric susceptibilities
    REAL(dp), ALLOCATABLE, DIMENSION (:) :: Cnn, Cnc, Ccc ! longitudinal and transverse dielectric susceptibilities
    
    REAL(dp) :: delta_k_cs, delta_k_cd, delta_k_cdelta, delta_k_chi_l, delta_k_chi_t, delta_k_Cnn, delta_k_Cnc, delta_k_Ccc

! parameters for branch ck_angular ONLY
    INTEGER, PARAMETER :: i1b = SELECTED_INT_KIND(2)
    INTEGER, PARAMETER :: spc = KIND((1.0,1.0))
    COMPLEX(spc), ALLOCATABLE, DIMENSION(:,:,:,:,:,:) :: ck_angular ! angular dcf in the molecular frame
    REAL(dp), ALLOCATABLE, DIMENSION (:) :: c_q ! dipole-charge correlation function
    REAL(dp) :: delta_k_ck_angular
    TYPE TYP_angleInd; INTEGER(i1b) :: costheta,phi,psi; END TYPE
    TYPE(TYP_angleInd), ALLOCATABLE, DIMENSION(:,:,:,:,:) :: angleInd ! integer table of correspondence omega(k,Omega)
    TYPE TYP_angleVal; REAL(dp) :: costheta,phi,psi; END TYPE
    TYPE(TYP_angleVal), ALLOCATABLE, DIMENSION(:,:,:,:,:) :: angleVal ! real omega values for interpolation in energy_ck_angular
    INTEGER(i2b) :: num_phi, num_cos, num_psi ! global variables for ck_angular, used in energy_ck_angular for reproducing angles

    CONTAINS
!-----------------------------------------------------------------------------------------------------------------------------------
    
    SUBROUTINE init

        IF ( input_log('readDensityDensityCorrelationFunction') ) THEN
            CALL readDensityDensityCorrelationFunction
        END IF

        IF ( input_log('polarization') ) THEN
            IF ( input_char('polarization_order')=='dipol' ) THEN
                CALL readPolarizationPolarizationCorrelationFunction
            ELSE IF (( input_char('polarization_order')=='multi') .AND. (.NOT. input_log('include_nc_coupling'))) THEN
                CALL readDielectricSusceptibilities   
            ELSE IF (( input_char('polarization_order')=='multi') .AND. (input_log('include_nc_coupling'))) THEN
                CALL readDielectricSusceptibilities 
                CALL readTotalPolarizationCorrelationFunction
            END IF
        END IF

        IF ( input_log("ck_angular") ) THEN
            CALL read_ck_angular
            CALL c_local_to_global_coordinates
        END IF
        IF ( input_log("ck_debug") ) THEN
            CALL readDensityDensityCorrelationFunction
            CALL readPolarizationPolarizationCorrelationFunction
        END IF
        IF ( input_log("ck_debug_extended") ) THEN
            CALL read_ck_projection
        END IF

        IF ( input_log('bridge_hard_sphere')) THEN
            IF( input_log("ck_angular") .OR. input_log('ck_debug') .OR. input_log('ck_debug_extended') ) STOP 'not implemented yet'
            CALL cs_of_k_hard_sphere! in case of bridge calculation, one also need the direct correlation function c2 of the hard sphere.
        END IF

    END SUBROUTINE init

!-----------------------------------------------------------------------------------------------------------------------------------

    SUBROUTINE read_ck_angular

        USE quadrature, ONLY: molRotSymOrder
        IMPLICIT NONE

        INTEGER(i2b) :: num_symm ! psi is from 0 to pi for water, no other symetry is taken into account
        LOGICAL :: exists

        INQUIRE(FILE="input/ck_angular.in", EXIST=exists)
        IF (.NOT. exists) STOP "FATAL ERROR: File input/ck_angular.in is not found."
        OPEN(39,FILE="input/ck_angular.in",FORM='UNFORMATTED')
            READ(39) nb_k, delta_k_ck_angular, num_cos, num_phi, num_psi, num_symm ! Note that psi is from 0 to pi for water, while no other symetry is taken into account
            ALLOCATE(ck_angular(num_psi,num_psi,num_phi,num_cos,num_cos,nb_k))
            REWIND(39)
            READ(39) nb_k, delta_k_ck_angular, num_cos, num_phi, num_psi, num_symm, ck_angular
            IF (num_symm /= molRotSymOrder) THEN
                PRINT*, 'FATAL ERROR: molRotSymOrder defined as ', molRotSymOrder, ' being different with num_symm ',num_symm
                STOP
            END IF
        CLOSE(39)

    END SUBROUTINE read_ck_angular

!-----------------------------------------------------------------------------------------------------------------------------------

    SUBROUTINE read_ck_projection

        USE quadrature, ONLY: molRotSymOrder
        IMPLICIT NONE

        INTEGER(i2b) :: num_symm ! psi is from 0 to pi for water, no other symetry is taken into account
        LOGICAL :: exists

        INQUIRE(FILE="input/ck_projection.in", EXIST=exists)
        IF (.NOT. exists) STOP "FATAL ERROR: File input/ck_projection.in is not found."
        OPEN(39,FILE="input/ck_projection.in",FORM='UNFORMATTED')
            READ(39) nb_k, delta_k_ck_angular
            ALLOCATE(c_s(nb_k)); ALLOCATE(c_q(nb_k)); ALLOCATE(c_delta(nb_k)); ALLOCATE(c_d(nb_k));
            REWIND(39)
            READ(39) nb_k, delta_k_ck_angular, num_symm, c_s, c_q, c_delta, c_d
            IF (num_symm /= molRotSymOrder) THEN
                PRINT*, 'FATAL ERROR: molRotSymOrder defined as ', molRotSymOrder, ' being different with num_symm ',num_symm
                STOP
            END IF
        CLOSE(39)

    END SUBROUTINE read_ck_projection

!-----------------------------------------------------------------------------------------------------------------------------------

    SUBROUTINE c_local_to_global_coordinates

        USE fft,             ONLY: kx,ky,kz
        USE constants,       ONLY: twopi
        USE quadrature,      ONLY: AngGrid,MolRotGrid,Rotxx,Rotxy,Rotxz,Rotyx,Rotyy,Rotyz,Rotzx,Rotzy,Rotzz,molRotSymOrder
        USE system,          ONLY: spaceGrid
        USE precision_kinds, ONLY: dp

        IMPLICIT NONE

        INTEGER(i2b) :: l,m,n,o,p,nfft1,nfft2,nfft3
        REAL(dp) :: costheta_k,sintheta_k,phi_k,cosphi_k,sinphi_k,cos_value,phi_value,psi_value
        REAL(dp) :: u_kx,u_ky,u_kz,v_kx,v_ky,v_kz,w_kx,w_ky,w_kz ! Projections of solvent axes u,v,w on k-frame kx,ky,kz
        LOGICAL :: karim

        nfft1 = spaceGrid%n_nodes(1)
        nfft2 = spaceGrid%n_nodes(2)
        nfft3 = spaceGrid%n_nodes(3)

        karim = input_log("ck_angular_interpolation")
        IF (karim) THEN
            ALLOCATE(angleVal(nfft1/2+1,nfft2,nfft3,angGrid%n_angles,molRotGrid%n_angles))
            angleVal%costheta = 0._dp
            angleVal%phi      = 0._dp
            angleVal%psi      = 0._dp
        ELSE
            ALLOCATE(angleInd(nfft1/2+1,nfft2,nfft3,angGrid%n_angles,molRotGrid%n_angles))
            angleInd%costheta = 0
            angleInd%phi      = 0
            angleInd%psi      = 0
        END IF

        DO CONCURRENT (l=1:nfft1/2+1, m=1:nfft2, n=1:nfft3)

        ! Definition of theta/phi for extreme cases
            IF (kx(l)**2+ky(m)**2+kz(n)**2/=0._dp) THEN
                costheta_k = kz(n)/(kx(l)**2+ky(m)**2+kz(n)**2)**0.5_dp
                sintheta_k = (1._dp - costheta_k**2)**0.5_dp
            ELSE
                costheta_k = 1._dp
                sintheta_k = 0._dp
            END IF

            IF(kx(l)**2+ky(m)**2/=0._dp) THEN
                phi_k =  angle(kx(l),ky(m))
                cosphi_k = COS(phi_k)
                sinphi_k = SIN(phi_k)
            ELSE
                cosphi_k = 1._dp
                sinphi_k = 0._dp
            END IF

        ! Open the loop Omega,Psi
            DO CONCURRENT (o=1:angGrid%n_angles, p=1:molRotGrid%n_angles)

            ! Calculate projections of solvent axes u,v,w on k-frame kx,ky,kz
                w_kz =   Rotxz(o,p)*sintheta_k*cosphi_k &
                       + Rotyz(o,p)*sintheta_k*sinphi_k &
                       + Rotzz(o,p)*costheta_k
                w_kx =   Rotxz(o,p)*costheta_k*cosphi_k &
                       + Rotyz(o,p)*costheta_k*sinphi_k &
                       - Rotzz(o,p)*sintheta_k
                w_ky = - Rotxz(o,p)*sinphi_k &
                       + Rotyz(o,p)*cosphi_k
                u_kz =   Rotxx(o,p)*sintheta_k*cosphi_k &
                       + Rotyx(o,p)*sintheta_k*sinphi_k &
                       + Rotzx(o,p)*costheta_k
            !   u_kx =   Rotxx(o,p)*costheta_k*cosphi_k &
            !          + Rotyx(o,p)*costheta_k*sinphi_k &
            !          - Rotzx(o,p)*sintheta_k
            !   u_ky = - Rotxx(o,p)*sinphi_k &
            !          + Rotyx(o,p)*cosphi_k
                v_kz =   Rotxy(o,p)*sintheta_k*cosphi_k &
                       + Rotyy(o,p)*sintheta_k*sinphi_k &
                       + Rotzy(o,p)*costheta_k
            !   v_kx =   Rotxy(o,p)*costheta_k*cosphi_k &
            !          + Rotyy(o,p)*costheta_k*sinphi_k &
            !          - Rotzy(o,p)*sintheta_k
            !   v_ky = - Rotxy(o,p)*sinphi_k &
            !          + Rotyy(o,p)*cosphi_k

            ! Calculate angles omega
                cos_value = w_kz
                phi_value = angle(w_kx,w_ky)
                psi_value = MODULO(angle(-u_kz,v_kz), twopi/molRotSymOrder)

            ! Real omega values for interpolation
                IF (karim) THEN
                    angleVal(l,m,n,o,p)%costheta = cos_value
                    angleVal(l,m,n,o,p)%phi      = phi_value
                    angleVal(l,m,n,o,p)%psi      = psi_value

            ! Index omega(k,Omega)
                ELSE
                    angleInd(l,m,n,o,p)%costheta = MIN(INT((1._dp + cos_value)*num_cos/2._dp) + 1, num_cos)
                    angleInd(l,m,n,o,p)%phi      = MOD(INT(phi_value*num_phi/twopi), num_phi) + 1
                    angleInd(l,m,n,o,p)%psi      = MOD(INT(psi_value*num_psi*molRotSymOrder/twopi), num_psi) + 1
                END IF
            END DO
        END DO

    ! Check final nomega
        IF ( .NOT. karim ) THEN
            IF ( ANY(angleInd%costheta<=0) .OR. ANY(angleInd%phi<=0) .OR. ANY(angleInd%psi<=0) ) THEN
                STOP "Some AngleInd is negative or zero"
            END IF
            IF ( ANY(angleInd%costheta>num_cos) .OR. ANY(angleInd%phi>num_phi) .OR. ANY(angleInd%psi>num_psi) ) THEN
                STOP "Some AngleInd is > to its max authorized value"
            END IF
        END IF

    END SUBROUTINE c_local_to_global_coordinates

!-----------------------------------------------------------------------------------------------------------------------------------

    ! Compute the angle between (0,x) and (x,y).
    PURE FUNCTION angle(x,y)
    

        USE precision_kinds,    ONLY: dp
        USE constants,          ONLY: twopi
        IMPLICIT NONE

        REAL(dp), INTENT(IN)    :: x,y
        REAL(dp)                :: angle
        REAL(dp)                :: xx,r

        IF (x==0._dp .AND. y==0._dp) THEN
            angle = 0._dp
        ELSE
            xx = ACOS( x/SQRT(x**2 + y**2) )
            IF (y>=0._dp) THEN
                angle = xx
            ELSE
                angle = twopi - xx
            END IF
        END IF

    END FUNCTION angle

!-----------------------------------------------------------------------------------------------------------------------------------

    SUBROUTINE readDielectricSusceptibilities ! chi_l, chi_t
        
        IMPLICIT NONE
        
        REAL(dp) :: norm_k
        LOGICAL :: exists
        CHARACTER(80) :: file_l, file_t
        INTEGER(i2b) :: ios, n_k, i

        
        file_l = 'input/direct_correlation_functions/water/chi_SPCE_for_multi/chi_l.in'
        file_t = 'input/direct_correlation_functions/water/chi_SPCE_for_multi/chi_t.in'
        INQUIRE (FILE=file_l, EXIST=exists )
            IF (.NOT. exists) THEN
                WRITE(*,*) "chi_l not found in ", file_l
                STOP
            END IF
        INQUIRE (FILE=file_t, EXIST=exists )
            IF (.NOT. exists) THEN
                WRITE(*,*) "chi_t not found in ", file_t
                STOP
            END IF
            
        n_k = MIN( n_linesInFile(file_l), n_linesInFile(file_t) )
        
        ALLOCATE ( chi_l (n_k), SOURCE=0._dp)
        ALLOCATE ( chi_t (n_k), SOURCE=0._dp)

        delta_k_chi_l = deltaAbscissa(file_l)
        delta_k_chi_t = deltaAbscissa(file_t)
        IF ( (delta_k_chi_t-delta_k_chi_l)/delta_k_chi_t >= 1.E-10 ) THEN
            WRITE(*,*)"chi_t and chi_l should have same delta k, i.e., same step in abscissa"
            STOP
        END IF
        delta_k= delta_k_chi_t

    If (.NOT. input_log('include_nc_coupling')) THEN
        IF ( (delta_k_chi_t-delta_k_cs)/delta_k_cs >=1.E-10 ) THEN
            WRITE(*,*)"chi_l, chi_t and c_s shoud have same delta k. THIS COULD BE IMPLEMENTED. ASK GUILLAUME"
            WRITE(*,*)'delta( chi_l )=',delta_k_chi_l
            WRITE(*,*)'delta( chi_t )=',delta_k_chi_t
            WRITE(*,*)'delta( cs )=',delta_k_cs
            STOP
        END IF
   END IF

        OPEN (10, FILE=file_l, iostat=ios)
            DO i = 1, SIZE(chi_l)
                READ (10,*,IOSTAT=ios) norm_k, chi_l(i)
                    IF (ios>0 .OR. ios<0) THEN
                        WRITE(*,*)'Error while reading ',file_l, 'in readDensityDensityCorrelationFunction (c_d)'
                        STOP
                    END IF
            END DO
        CLOSE (10)
        OPEN (10, FILE=file_t, iostat=ios)
            DO i = 1, SIZE(chi_t)
                READ (10,*,IOSTAT=ios) norm_k, chi_t(i)
                    IF (ios>0 .OR. ios<0) THEN
                        WRITE(*,*)'Error while reading ',file_t, 'in readDensityDensityCorrelationFunction (c_d)'
                        STOP
                    END IF
            END DO
        CLOSE (10)
        
        nb_k=n_k

    END SUBROUTINE readDielectricSusceptibilities


    SUBROUTINE readTotalPolarizationCorrelationFunction !Cnn, Cnc, Ccc
        
        IMPLICIT NONE
        
        REAL(dp) :: norm_k
        LOGICAL :: exists
        CHARACTER(80) :: file_nn, file_nc, file_cc
        INTEGER(i2b) :: ios, nb_k_in_Cnn, i,nb_k_in_Cnc,nb_k_in_Ccc

        
        file_nn = 'input/direct_correlation_functions/water/Cnn.dat'
        file_nc = 'input/direct_correlation_functions/water/Cnc.dat'
        file_cc = 'input/direct_correlation_functions/water/Ccc.dat'
        
        INQUIRE (FILE=file_nn, EXIST=exists )
            IF (.NOT. exists) THEN
                WRITE(*,*) "Cnn not found in ", file_nn
                STOP
            END IF
        INQUIRE (FILE=file_nc, EXIST=exists )
            IF (.NOT. exists) THEN
                WRITE(*,*) "Cnc not found in ", file_nc
                STOP
            END IF
        INQUIRE (FILE=file_cc, EXIST=exists )
            IF (.NOT. exists) THEN
                WRITE(*,*) "Ccc not found in ", file_cc
                STOP
            END IF
                        
        nb_k_in_c = MIN( n_linesInFile(file_nn), n_linesInFile(file_nc), n_linesInFile(file_cc) )
        
        ALLOCATE ( Cnn (nb_k_in_c), SOURCE=0._dp)
        ALLOCATE ( Cnc (nb_k_in_c), SOURCE=0._dp)
        ALLOCATE ( Ccc (nb_k_in_c), SOURCE=0._dp)

        delta_k_Cnn = deltaAbscissa(file_nn)
        delta_k_Cnc = deltaAbscissa(file_nc)
        delta_k_Ccc = deltaAbscissa(file_cc)
        
        IF ( (delta_k_Cnn-delta_k_Cnc)/delta_k_Cnn >= 1.E-10 ) THEN
            WRITE(*,*)"Cnn and Cnc should have same delta k, i.e., same step in abscissa"
            STOP
        END IF
        IF ( (delta_k_Cnn-delta_k_Ccc)/delta_k_Cnn >= 1.E-10 ) THEN
            WRITE(*,*)"Cnn and Ccc should have same delta k, i.e., same step in abscissa"
            STOP
        END IF
        delta_k_in_C = delta_k_Cnn

        OPEN (10, FILE=file_nn, iostat=ios)
            DO i = 1, SIZE(Cnn)
                READ (10,*,IOSTAT=ios) norm_k, Cnn(i)
                    IF (ios>0 .OR. ios<0) THEN
                        WRITE(*,*)'Error while reading ',Cnn, 'in readTotalPolarizationCorrelationFunction (Cnn)'
                        STOP
                    END IF
            END DO
        CLOSE (10)
        OPEN (10, FILE=file_nc, iostat=ios)
            DO i = 1, SIZE(Cnc)
                READ (10,*,IOSTAT=ios) norm_k, Cnc(i)
                    IF (ios>0 .OR. ios<0) THEN
                        WRITE(*,*)'Error while reading ',Cnc, 'in readTotalPolarizationCorrelationFunction (Cnc)'
                        STOP
                    END IF
            END DO
        CLOSE (10)
        OPEN (10, FILE=file_cc, iostat=ios)
            DO i = 1, SIZE(Ccc)
                READ (10,*,IOSTAT=ios) norm_k, Ccc(i)
                    IF (ios>0 .OR. ios<0) THEN
                        WRITE(*,*)'Error while reading ',Ccc, 'in readTotalPolarizationCorrelationFunction (Ccc)'
                        STOP
                    END IF
            END DO
        CLOSE (10)

    END SUBROUTINE readTotalPolarizationCorrelationFunction



















    SUBROUTINE readPolarizationPolarizationCorrelationFunction ! c_delta, c_d
        INTEGER(i2b) :: ios, nb_kdelta, i, nb_kd
        REAL(dp) :: norm_k
        CHARACTER(80) :: filename, ck_species
    
        ck_species = input_char('ck_species')
        
        IF ( ck_species == 'spc  ' ) THEN
            filename='input/direct_correlation_functions/water/SPC_Lionel_Daniel/cdelta.in'
        ELSE IF ( ck_species == 'stock' ) THEN
            filename='input/direct_correlation_functions/stockmayer/cdelta.in'
        ELSE IF ( ck_species == 'perso' ) THEN
            filename='input/cdelta.in'
        ELSE IF ( ck_species == 'spce' ) then
            filename='input/direct_correlation_functions/water/SPCE/cdelta.in'
        END IF

        nb_kdelta = n_linesInFile(filename)
        ALLOCATE(c_delta(nb_kdelta), SOURCE=0._dp)
        delta_k = deltaAbscissa(filename)
        OPEN (13, FILE=filename, IOSTAT=ios)
            IF (ios/=0) THEN
                WRITE(*,*)'Cant open file ',filename,' in readDensityDensityCorrelationFunction (c_delta)'
                STOP
            END IF
            DO i = 1, SIZE(c_delta)
                READ (13,*,IOSTAT=ios) norm_k, c_delta(i)
                    IF (ios>0 .OR. ios<0) THEN
                        WRITE(*,*)'Error while reading ',filename, 'in readDensityDensityCorrelationFunction (c_delta)'
                        STOP
                    END IF
            END DO
        CLOSE (13)


        IF ( ck_species == 'spc  ' ) THEN
            filename='input/direct_correlation_functions/water/SPC_Lionel_Daniel/cd.in'
        ELSE IF ( ck_species == 'stock' ) THEN
            filename='input/direct_correlation_functions/stockmayer/cd.in'
        ELSE IF ( ck_species == 'perso' ) THEN
            filename='input/cd.in'
        ELSE IF ( ck_species == 'spce' ) then
            filename='input/direct_correlation_functions/water/SPCE/cd.in'
        END IF
        nb_kd = n_linesInFile(filename)
        ALLOCATE(c_d(nb_kd), SOURCE=0._dp)
        
        
        delta_k = deltaAbscissa(filename)
        OPEN (13, FILE=filename, IOSTAT=ios)
            IF (ios/=0) THEN
                WRITE(*,*)'Cant open file ',filename,' in readDensityDensityCorrelationFunction (c_d)'
                STOP
            END IF
            DO i = 1, SIZE(c_d)
                READ (13,*,IOSTAT=ios) norm_k, c_d(i)
                    IF (ios>0 .OR. ios<0) THEN
                        WRITE(*,*)'Error while reading ',filename, 'in readDensityDensityCorrelationFunction (c_d)'
                        STOP
                    END IF
            END DO
        CLOSE (13)

    END SUBROUTINE readPolarizationPolarizationCorrelationFunction
    




    SUBROUTINE readDensityDensityCorrelationFunction ! c_s

        CHARACTER(80) :: filename, ck_species
        INTEGER(i2b) :: ios, i
        REAL(dp) :: norm_k
        
        ck_species = input_char('ck_species')
        IF ( ck_species == 'spc  ' ) THEN
            filename = 'input/direct_correlation_functions/water/SPC_Lionel_Daniel/cs.in'
        ELSE IF ( ck_species == 'stock' ) THEN
            filename = 'input/direct_correlation_functions/stockmayer/cs.in'
        ELSE IF ( ck_species == 'perso' ) THEN
            filename = 'input/cs.in'
        ELSE IF ( ck_species == 'spce' ) then
            filename = 'input/direct_correlation_functions/water/SPCE/cs.in'
        END IF

        delta_k = deltaAbscissa(filename)
        delta_k_cs = delta_k
        nb_k = n_linesInFile(filename)
        
        ALLOCATE ( c_s(nb_k), SOURCE=0._dp )
        
        OPEN (13, FILE=filename, IOSTAT=ios)
            IF (ios/=0) THEN
                WRITE(*,*)'Cant open file ',filename,' in readDensityDensityCorrelationFunction (c_s)'
                STOP
            END IF
            DO i = 1, SIZE(c_s)
                READ (13,*,IOSTAT=ios) norm_k, c_s(i)
                    IF (ios>0 .OR. ios<0) THEN
                        WRITE(*,*)'Error while reading ',filename, 'in readDensityDensityCorrelationFunction (c_s)'
                        STOP
                    END IF
            END DO
        CLOSE (13)

    END SUBROUTINE readDensityDensityCorrelationFunction



END MODULE dcf

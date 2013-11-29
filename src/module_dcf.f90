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


    CONTAINS
    
    
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

        IF ( input_log('bridge_hard_sphere')) THEN
            CALL cs_of_k_hard_sphere! in case of bridge calculation, one also need the direct correlation function c2 of the hard sphere.
        END IF
        
 
       
    END SUBROUTINE init


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
        delta_k_chi_l = delta_k

        IF ( (delta_k_chi_t-delta_k_cs)/delta_k_cs >=1.E-10 ) THEN
            WRITE(*,*)"chi_l, chi_t and c_s shoud have same delta k. THIS COULD BE IMPLEMENTED. ASK GUILLAUME"
            WRITE(*,*)'delta( chi_l )=',delta_k_chi_l
            WRITE(*,*)'delta( chi_t )=',delta_k_chi_t
            WRITE(*,*)'delta( cs )=',delta_k_cs
            STOP
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
    
        ck_species = TRIM(ADJUSTL(input_char('ck_species')))
        
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
        
        ck_species = TRIM(ADJUSTL(input_char('ck_species')))
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

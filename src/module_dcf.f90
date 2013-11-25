MODULE dcf

    USE precision_kinds, ONLY: dp, i2b
    USE input, ONLY: input_log, input_char, n_linesInFile, deltaAbscissa
    USE system, ONLY: delta_kSYS => delta_k, nb_kSYS => nb_k

    IMPLICIT NONE

    REAL(dp) :: delta_k ! distance between two k points in cs.in, cdelta.in, cd.in
    INTEGER(i2b) :: nb_k ! nb of k points in cs.in, cdelta.in, cd.in
    REAL(dp), ALLOCATABLE, DIMENSION (:) :: c_s ! density density correlation function
    REAL(dp), ALLOCATABLE, DIMENSION (:) :: c_delta, c_d ! polarization polarization correlation function


    CONTAINS
    
    
    SUBROUTINE init

        IF ( input_log('readDensityDensityCorrelationFunction') ) THEN
            CALL readDensityDensityCorrelationFunction
        END IF

        IF ( input_log('polarization') ) THEN
            IF ( input_char('polarization_order')=='dipol' ) THEN
                CALL readPolarizationPolarizationCorrelationFunction
            ELSE IF ( input_char('polarization_order')=='multi' ) THEN
                CALL readDielectricSusceptibilities
            END IF
        END IF

        IF ( input_log('bridge_hard_sphere')) THEN
            CALL cs_of_k_hard_sphere! in case of bridge calculation, one also need the direct correlation function c2 of the hard sphere.
        END IF
        
    END SUBROUTINE init
    






    SUBROUTINE readPolarizationPolarizationCorrelationFunction ! c_delta, c_d
        INTEGER(i2b) :: ios, nb_kdelta, nk, nb_kd
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
            DO nk = 1, SIZE(c_delta)
                READ (13,*,IOSTAT=ios) norm_k, c_delta(nk)
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
            DO nk = 1, SIZE(c_d)
                READ (13,*,IOSTAT=ios) norm_k, c_d(nk)
                    IF (ios>0 .OR. ios<0) THEN
                        WRITE(*,*)'Error while reading ',filename, 'in readDensityDensityCorrelationFunction (c_d)'
                        STOP
                    END IF
            END DO
        CLOSE (13)

    END SUBROUTINE readPolarizationPolarizationCorrelationFunction
    





    SUBROUTINE readDielectricSusceptibilities ! chi_l, chi_t
    
    END SUBROUTINE readDielectricSusceptibilities




    SUBROUTINE readDensityDensityCorrelationFunction ! c_s

        CHARACTER(80) :: filename, ck_species
        INTEGER(i2b) :: ios, nk
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
        delta_kSYS = delta_k !should be removed when delta_k will be used from module_dcf everywhere
        nb_k = n_linesInFile(filename)
        nb_kSYS = nb_k !should be removed when nb_k will be used from module_dcf everywhere
        
        ALLOCATE ( c_s(nb_k), SOURCE=0._dp )
        
        OPEN (13, FILE=filename, IOSTAT=ios)
            IF (ios/=0) THEN
                WRITE(*,*)'Cant open file ',filename,' in readDensityDensityCorrelationFunction (c_s)'
                STOP
            END IF
            DO nk = 1, SIZE(c_s)
                READ (13,*,IOSTAT=ios) norm_k, c_s(nk)
                    IF (ios>0 .OR. ios<0) THEN
                        WRITE(*,*)'Error while reading ',filename, 'in readDensityDensityCorrelationFunction (c_s)'
                        STOP
                    END IF
            END DO
        CLOSE (13)

    END SUBROUTINE readDensityDensityCorrelationFunction



END MODULE dcf

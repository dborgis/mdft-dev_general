PROGRAM lancer

    IMPLICIT NONE

    INTEGER             :: i,j,k,n,njob_min,njob_max,n_lines,lng,lng_min,lng_max,lng_sp,leb,psi,psi_min,psi_max,sol,nsol,pd
    INTEGER             :: no,noption,met,nmethod,total_job_nb
    REAL                :: reso
    LOGICAL             :: karim
    CHARACTER(LEN=3)    :: foldername,solutenumber
    CHARACTER(LEN=1)    :: temp
    CHARACTER(LEN=100)  :: char_temp,text,char_option
    CHARACTER(LEN=100), ALLOCATABLE, DIMENSION(:) :: input_line,arraytemp
    CHARACTER(LEN=25),               DIMENSION(6) :: method_name = ['dipole                   ',&
                                                                    'ck_debug                 ',&
                                                                    'ck_debug_extended        ',&
                                                                    'ck_angular               ',&
                                                                    'ck_angular_interpolation ',&
                                                                    'default                  ']
    INTEGER,            ALLOCATABLE, DIMENSION(:) :: option_temp,option,method

! read options
    WRITE(*,*) 'POINCARE JOB LANCER VER_0.2'
    WRITE(*,*) '--'
    WRITE(*,'(A)',ADVANCE="NO") 'Job id range: '
    READ(*,*) njob_min, njob_max
    WRITE(*,*) 'Options: '
    WRITE(*,*) '1) Lx/Ly/Lz'
    WRITE(*,*) '2) Order_of_quadrature (Lebedev) (6/14 both)'
    WRITE(*,*) '3) Psi'
    WRITE(*,*) '4) Solute'
    WRITE(*,*) '5) Method (dipole, ck_debug, ck_debug_extended, ck_angular, ck_angular_interpolation)'
    WRITE(*,*) '6) Poisson_solver/direct_sum (both)'
    WRITE(*,'(A)',ADVANCE="NO") 'Choose number of variable options (outer to inner): '
    READ(*,'(A100)') char_option
    ALLOCATE(option_temp(6))
    karim = .FALSE.
    noption = 0
    nsol = 1; sol = 0
    DO i = 1, 100
        DO j = 1, 6
            WRITE(temp,'(i1)') j
            IF (char_option(i:i) == temp) THEN
                noption = noption + 1
                option_temp(noption) = j
            END IF
        END DO
    END DO
    ALLOCATE(option(noption))
    option = option_temp(1:noption)
    DEALLOCATE(option_temp)

! read input/dft.in
    n_lines = n_linesInFile('input/dft.in')
    ALLOCATE(input_line(n_lines))
    OPEN(11,FILE='input/dft.in')
        DO i = 1, n_lines
            READ(11,'(a)') text
            input_line(i) = trim(adjustl(text))
        END DO
    CLOSE (11)
    DO i = 1, n_lines
        DO j = 1 , LEN(text)
            IF (input_line(i)(j:j) == '#') THEN
                DO CONCURRENT (k=j:LEN(text))
                    input_line(i)(k:k) = ' '
                END DO
                EXIT
            END IF
        END DO
        input_line(i) = TRIM(ADJUSTL(input_line(i)))
    END DO
    n = 0
    DO i = 1, n_lines
        IF (input_line(i)(1:1) /= ' ') THEN
            input_line(n+1) = input_line(i)
            n = n + 1
        END IF
    END DO
    ALLOCATE (arraytemp(n),SOURCE=input_line(1:n))
    DEALLOCATE (input_line)
    ALLOCATE (input_line(n), SOURCE=arraytemp)
    DEALLOCATE (arraytemp)
    DO j = 1 , SIZE(input_line)
        IF (input_line(j)(1:2)=='Lx')                   READ(input_line(j)( 6: 9),'(i4)') lng
        IF (input_line(j)(1:4)=='nfft')                 READ(input_line(j)( 9:12),*     ) reso
        IF (input_line(j)(1:19)=='order_of_quadrature') READ(input_line(j)(23:26),'(i4)') leb
        IF (input_line(j)(1:6)=='nb_psi')               READ(input_line(j)(10:13),'(i4)') psi
    END DO
    reso = reso/lng

! read parameters
    total_job_nb = 1
    met = 1; ALLOCATE(method(1)); method(1) = 6
    DO n = 1, noption
        SELECT CASE(option(n))
        CASE(1)
            WRITE(*,'(A)',ADVANCE="NO") 'Box length (cube) range, spacing: '
            READ(*,*) lng_min, lng_max, lng_sp
            WRITE(*,'(A)',ADVANCE="NO") 'Resolution (nfft/L): '
            READ(*,*) reso
            total_job_nb = total_job_nb * ((lng_max - lng_min) / lng_sp + 1)
        CASE(2)
            total_job_nb = total_job_nb * 2
        CASE(3)
            WRITE(*,'(A)',ADVANCE="NO") 'Psi range: '
            READ(*,*) psi_min, psi_max
            total_job_nb = total_job_nb * (psi_max - psi_min + 1)
        CASE(4)
            WRITE(*,'(A)',ADVANCE="NO") 'Number of solute: '
            READ(*,*) nsol
            total_job_nb = total_job_nb * nsol
        CASE(5)
            WRITE(*,'(A)',ADVANCE="NO") 'Choose number of method: 1) dipole, 2) ck_debug, 3) ck_debug_extended, 4) ck_angular, 5) &
                                        &ck_angular_interpolation (But DCFs must be copied by hand): '
            READ(*,'(A100)') char_option
            ALLOCATE(option_temp(5))
            nmethod = 0
            DO i = 1, 100
                DO j = 1, 6
                    WRITE(temp,'(i1)') j
                    IF (char_option(i:i) == temp) THEN
                        nmethod = nmethod + 1
                        option_temp(nmethod) = j
                    END IF
                END DO
            END DO
            DEALLOCATE(method)
            ALLOCATE(method(nmethod))
            method = option_temp(1:nmethod)
            DEALLOCATE(option_temp)
            total_job_nb = total_job_nb * nmethod
        CASE(6)
            total_job_nb = total_job_nb * 2
        END SELECT
    END DO
    IF (njob_max - njob_min + 1 /= total_job_nb) STOP 'FATAL ERROR: Total job_nb error.'

! generate job folders
    CALL EXECUTE_COMMAND_LINE('mkdir -p taskgen/input')
    CALL EXECUTE_COMMAND_LINE('cp input/solvent.in taskgen/input')
    PRINT*, ' ID | SOL |  L  | nfft | Leb | psi | method'
    i = njob_min - 1
    no = noption
    CALL generator(option,noption,i,no)
    IF (no /= noption) STOP 'no /= noption'
    IF (i /= njob_max) STOP 'FATAL ERROR: Loossing jobs.'
    char_temp = 'rm -rf taskgen/input'
    CALL EXECUTE_COMMAND_LINE(char_temp)

! generate operating files
    OPEN(21,FILE='lancer.sh')
        DO i = njob_min, njob_max
            WRITE(foldername,'(i3.3)') i
            WRITE(21,*) 'echo ',foldername
            WRITE(21,*) 'cd ',foldername
            WRITE(21,*) 'llsubmit x.sh'
            WRITE(21,*) 'cd ..'
        END DO
        WRITE(21,*) 'rm lancer.sh'
    CLOSE(21)

    PRINT*, '--'
    PRINT*, 'MISSION COMPLETE POINCARE_JOB_LANCER'

    CONTAINS
!-----------------------------------------------------------------------------------------------------------------------------------
    RECURSIVE SUBROUTINE generator(option,noption,i,no)

        IMPLICIT NONE
        INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(IN)  :: option
        INTEGER, INTENT(IN)                             :: noption
        INTEGER, INTENT(INOUT)                          :: i,no
        INTEGER                                         :: j

        IF (no == 0) THEN
            i = i + 1
            WRITE(foldername,'(i3.3)') i
            char_temp = 'rm -rf ' // foldername
            CALL EXECUTE_COMMAND_LINE(char_temp)
            char_temp = 'cp -r taskgen ' // foldername
            CALL EXECUTE_COMMAND_LINE(char_temp)
            char_temp = foldername // '/input/dft.in'
            OPEN(21,FILE=char_temp)
                WRITE(21,'(a,a)') '# JOBID = ',foldername
                DO j = 1 , SIZE(input_line)
                    WRITE (21,*) input_line(j)
                END DO
            CLOSE(21)
            IF (nsol /= 1) THEN
                WRITE(solutenumber,'(i3.3)') sol
                char_temp = 'cp input/solute.' // solutenumber // ' ' // foldername // '/input/solute.in'
                CALL EXECUTE_COMMAND_LINE(char_temp)
            ELSE
                WRITE(solutenumber,'(i3.3)') sol
                char_temp = 'cp input/solute.in ' // foldername // '/input/solute.in'
                CALL EXECUTE_COMMAND_LINE(char_temp)
            END IF
            WRITE(*,'(" ",a3,"   ",a3,"   ",i3,"   ",i4,"   ",i3,"   ",i3,"   ",a25)') &
            foldername,solutenumber,lng,INT(lng*reso),leb,psi,method_name(method(met))

        ELSE IF (no > 0) THEN
            SELECT CASE(option(noption-no+1))
            CASE(1)
                no = no - 1
                DO lng = lng_min, lng_max, lng_sp
                    DO j = 1 , SIZE(input_line)
                        IF (input_line(j)(1:2)=='Lx') THEN
                            WRITE(input_line(j),   '("Lx = ",i3,"             ")') lng
                            WRITE(input_line(j+1), '("Ly = ",i3,"             ")') lng
                            WRITE(input_line(j+2), '("Lz = ",i3,"             ")') lng
                        END IF
                        IF (input_line(j)(1:5)=='nfft1') THEN
                            WRITE(input_line(j),   '("nfft1 = ",i3,"             ")') INT(lng*reso)
                            WRITE(input_line(j+1), '("nfft2 = ",i3,"             ")') INT(lng*reso)
                            WRITE(input_line(j+2), '("nfft3 = ",i3,"             ")') INT(lng*reso)
                        END IF
                    END DO
                    CALL generator(option,noption,i,no)
                END DO
                no = no + 1
            CASE(2)
                no = no - 1
                DO leb = 6, 14, 8
                    DO j = 1 , SIZE(input_line)
                        IF (input_line(j)(1:19)=='order_of_quadrature') THEN
                            WRITE(input_line(j),'("order_of_quadrature = ",i3,"             ")') leb
                        END IF
                    END DO
                    CALL generator(option,noption,i,no)
                END DO
                no = no + 1
            CASE(3)
                no = no - 1
                DO psi = psi_min, psi_max
                    DO j = 1 , SIZE(input_line)
                        IF (input_line(j)(1:6)=='nb_psi') THEN
                            WRITE(input_line(j),'("nb_psi = ",i3,"             ")') psi
                        END IF
                    END DO
                    CALL generator(option,noption,i,no)
                END DO
                no = no + 1
            CASE(4)
                no = no - 1
                DO sol = 1, nsol
                    CALL generator(option,noption,i,no)
                END DO
                no = no + 1
            CASE(5)
                no = no - 1
                DO met = 1, nmethod
                    DO j = 1 , SIZE(input_line)
                        SELECT CASE(method(met))
                        CASE(1)
                            IF (input_line(j)(1:39)=='readDensityDensityCorrelationFunction =') THEN
                                WRITE(input_line(j),'(a)') 'readDensityDensityCorrelationFunction = T'
                            END IF
                            IF (input_line(j)(1:14)=='polarization =') THEN
                                WRITE(input_line(j),'(a)') 'polarization = T'
                            END IF
                            IF (input_line(j)(1:12)=='ck_angular =') THEN
                                WRITE(input_line(j),'(a)') 'ck_angular = F'
                            END IF
                            IF (input_line(j)(1:26)=='ck_angular_interpolation =') THEN
                                WRITE(input_line(j),'(a)') 'ck_angular_interpolation = F'
                            END IF
                            IF (input_line(j)(1:10)=='ck_debug =') THEN
                                WRITE(input_line(j),'(a)') 'ck_debug = F'
                            END IF
                            IF (input_line(j)(1:19)=='ck_debug_extended =') THEN
                                WRITE(input_line(j),'(a)') 'ck_debug_extended = F'
                            END IF
                        CASE(2)
                            IF (input_line(j)(1:39)=='readDensityDensityCorrelationFunction =') THEN
                                WRITE(input_line(j),'(a)') 'readDensityDensityCorrelationFunction = F'
                            END IF
                            IF (input_line(j)(1:14)=='polarization =') THEN
                                WRITE(input_line(j),'(a)') 'polarization = F'
                            END IF
                            IF (input_line(j)(1:12)=='ck_angular =') THEN
                                WRITE(input_line(j),'(a)') 'ck_angular = F'
                            END IF
                            IF (input_line(j)(1:26)=='ck_angular_interpolation =') THEN
                                WRITE(input_line(j),'(a)') 'ck_angular_interpolation = F'
                            END IF
                            IF (input_line(j)(1:10)=='ck_debug =') THEN
                                WRITE(input_line(j),'(a)') 'ck_debug = T'
                            END IF
                            IF (input_line(j)(1:19)=='ck_debug_extended =') THEN
                                WRITE(input_line(j),'(a)') 'ck_debug_extended = F'
                            END IF
                        CASE(3)
                            IF (input_line(j)(1:39)=='readDensityDensityCorrelationFunction =') THEN
                                WRITE(input_line(j),'(a)') 'readDensityDensityCorrelationFunction = F'
                            END IF
                            IF (input_line(j)(1:14)=='polarization =') THEN
                                WRITE(input_line(j),'(a)') 'polarization = F'
                            END IF
                            IF (input_line(j)(1:12)=='ck_angular =') THEN
                                WRITE(input_line(j),'(a)') 'ck_angular = F'
                            END IF
                            IF (input_line(j)(1:26)=='ck_angular_interpolation =') THEN
                                WRITE(input_line(j),'(a)') 'ck_angular_interpolation = F'
                            END IF
                            IF (input_line(j)(1:10)=='ck_debug =') THEN
                                WRITE(input_line(j),'(a)') 'ck_debug = F'
                            END IF
                            IF (input_line(j)(1:19)=='ck_debug_extended =') THEN
                                WRITE(input_line(j),'(a)') 'ck_debug_extended = T'
                            END IF
                        CASE(4)
                            IF (input_line(j)(1:39)=='readDensityDensityCorrelationFunction =') THEN
                                WRITE(input_line(j),'(a)') 'readDensityDensityCorrelationFunction = F'
                            END IF
                            IF (input_line(j)(1:14)=='polarization =') THEN
                                WRITE(input_line(j),'(a)') 'polarization = F'
                            END IF
                            IF (input_line(j)(1:12)=='ck_angular =') THEN
                                WRITE(input_line(j),'(a)') 'ck_angular = T'
                            END IF
                            IF (input_line(j)(1:26)=='ck_angular_interpolation =') THEN
                                WRITE(input_line(j),'(a)') 'ck_angular_interpolation = F'
                            END IF
                            IF (input_line(j)(1:10)=='ck_debug =') THEN
                                WRITE(input_line(j),'(a)') 'ck_debug = F'
                            END IF
                            IF (input_line(j)(1:19)=='ck_debug_extended =') THEN
                                WRITE(input_line(j),'(a)') 'ck_debug_extended = F'
                            END IF
                        CASE(5)
                            IF (input_line(j)(1:39)=='readDensityDensityCorrelationFunction =') THEN
                                WRITE(input_line(j),'(a)') 'readDensityDensityCorrelationFunction = F'
                            END IF
                            IF (input_line(j)(1:14)=='polarization =') THEN
                                WRITE(input_line(j),'(a)') 'polarization = F'
                            END IF
                            IF (input_line(j)(1:12)=='ck_angular =') THEN
                                WRITE(input_line(j),'(a)') 'ck_angular = T'
                            END IF
                            IF (input_line(j)(1:26)=='ck_angular_interpolation =') THEN
                                WRITE(input_line(j),'(a)') 'ck_angular_interpolation = T'
                            END IF
                            IF (input_line(j)(1:10)=='ck_debug =') THEN
                                WRITE(input_line(j),'(a)') 'ck_debug = F'
                            END IF
                            IF (input_line(j)(1:19)=='ck_debug_extended =') THEN
                                WRITE(input_line(j),'(a)') 'ck_debug_extended = F'
                            END IF
                        END SELECT
                    END DO
                    CALL generator(option,noption,i,no)
                END DO
                no = no + 1
            CASE(6)
                no = no - 1
                DO pd = 1, 2
                    DO j = 1 , SIZE(input_line)
                        SELECT CASE(pd)
                        CASE(1)
                            IF (input_line(j)(1:16)=='poisson_solver =') THEN
                                WRITE(input_line(j),'(a)') 'poisson_solver = F'
                                karim = .true.
                            END IF
                            IF (input_line(j)(1:12)=='direct_sum =') THEN
                                WRITE(input_line(j),'(a)') 'direct_sum = T'
                                karim = karim .AND. .true.
                            END IF
                        CASE(2)
                            IF (input_line(j)(1:16)=='poisson_solver =') THEN
                                WRITE(input_line(j),'(a)') 'poisson_solver = T'
                                karim = .true.
                            END IF
                            IF (input_line(j)(1:12)=='direct_sum =') THEN
                                WRITE(input_line(j),'(a)') 'direct_sum = F'
                                karim = karim .AND. .true.
                            END IF
                        END SELECT
                    END DO
                    IF (.NOT. karim) STOP 'FATAL ERROR: Cannot find poisson_solver or direct_sum.'
                    CALL generator(option,noption,i,no)
                END DO
                no = no + 1
            CASE DEFAULT
                STOP 'FATAL ERROR: Cannot find option.'
            END SELECT
        END IF

    END SUBROUTINE generator
!-----------------------------------------------------------------------------------------------------------------------------------
    FUNCTION n_linesInFile (filename)

        IMPLICIT NONE
        INTEGER :: n_linesInFile
        CHARACTER(*), INTENT(IN) :: filename
        INTEGER :: ios

        OPEN (77, FILE=filename)
        n_linesInFile = 0
        DO WHILE (.true.)
            READ (77,*,IOSTAT=ios)
            IF (ios>0) THEN
                WRITE(*,*)'Error in file:',filename
                STOP
            ELSE IF (ios<0) THEN ! end of file reached
                EXIT
            ELSE
                n_linesInFile = n_linesInFile + 1
            END IF
        END DO
        CLOSE (77)

    END FUNCTION
!-----------------------------------------------------------------------------------------------------------------------------------
END PROGRAM lancer
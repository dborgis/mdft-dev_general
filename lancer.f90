!-----------------------------------------------------------------------------------------------------------------------------------
MODULE precision_kinds

    IMPLICIT NONE
! symbolic names for kind types of 4-, 2-, and 1-byte integers:
    INTEGER, PARAMETER :: i4b = SELECTED_INT_KIND(9)
    INTEGER, PARAMETER :: i2b = SELECTED_INT_KIND(4)
    INTEGER, PARAMETER :: i1b = SELECTED_INT_KIND(2)
! symbolic names for kind types of single- and double-precision reals:
    INTEGER, PARAMETER :: sp = KIND(1.0)
    INTEGER, PARAMETER :: dp = KIND(1.0d0)
! symbolic names for kind types of single- and double-precision complex:
    INTEGER, PARAMETER :: spc = KIND((1.0,1.0))
    INTEGER, PARAMETER :: dpc = KIND((1.0d0,1.0d0))
! symbolic name for kind type of default logical:
    INTEGER, PARAMETER :: lgt = KIND(.true.)
! parameters:
    REAL(dp), PARAMETER :: pi=4._dp*ATAN(1._dp), twopi=8._dp*ATAN(1._dp)
    COMPLEX(dpc), PARAMETER :: imag = (0._dp,1.0_dp)

END MODULE precision_kinds
!-----------------------------------------------------------------------------------------------------------------------------------
PROGRAM lancer

    USE precision_kinds, ONLY: i2b,dp
    IMPLICIT NONE

    INTEGER(i2b)        :: i,j,k,l,m,n,o,p,njob_min,njob_max,lng_min,lng_max,leb,psi_min,psi_max,nsol,n_lines,lng,psi,sol,lng_sp
    REAL(dp)            :: reso
    CHARACTER(LEN=3)    :: foldername,solutenumber
    CHARACTER(LEN=100)  :: char_temp,text
    CHARACTER(LEN=100), ALLOCATABLE, DIMENSION(:) :: input_line,arraytemp

! read parameters
    WRITE(*,*) 'POINCARE JOB LANCER VER_0.1'
    WRITE(*,*) '--'
    WRITE(*,'(A)',ADVANCE="NO") 'Job id range: '
    READ(*,*) njob_min, njob_max
    WRITE(*,'(A)',ADVANCE="NO") 'Box length (cube) range, spacing: '
    READ(*,*) lng_min, lng_max, lng_sp
    WRITE(*,'(A)',ADVANCE="NO") 'Resolution (nfft/L): '
    READ(*,*) reso
    WRITE(*,'(A)',ADVANCE="NO") 'Lebedev: '
    READ(*,*) leb
    WRITE(*,'(A)',ADVANCE="NO") 'Psi range: '
    READ(*,*) psi_min, psi_max
    WRITE(*,'(A)',ADVANCE="NO") 'Number of solute: '
    READ(*,*) nsol

    IF (njob_max - njob_min + 1 /= ((lng_max - lng_min) / lng_sp + 1)*(psi_max - psi_min + 1)*nsol) THEN
        STOP 'FATAL ERROR: Total job_nb error.'
    END IF

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

! generate job folders
    CALL SYSTEM('mkdir -p taskgen/input')
    CALL SYSTEM('cp input/*.in taskgen/input')
    CALL SYSTEM('cp mdft taskgen')
    CALL SYSTEM('cp x.sh taskgen')
    CALL SYSTEM('rm taskgen/input/dft.in taskgen/input/solute.in')

    i = njob_min - 1
    DO sol = 1, nsol; DO lng = lng_min, lng_max, lng_sp; DO psi = psi_min, psi_max
        i = i + 1
        WRITE(foldername,'(i3.3)') i
        char_temp = 'rm -rf ' // foldername
        CALL SYSTEM(char_temp)
        char_temp = 'cp -r taskgen ' // foldername
        CALL SYSTEM(char_temp)
        char_temp = foldername // '/input/dft.in'
        OPEN(21,FILE=char_temp)
            WRITE(21,'(a,a)') '# JOBID = ',foldername
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
                IF (input_line(j)(1:19)=='order_of_quadrature') THEN
                    WRITE(input_line(j),   '("order_of_quadrature = ",i3,"             ")') leb
                END IF
                IF (input_line(j)(1:6)=='nb_psi') THEN
                    WRITE(input_line(j),   '("nb_psi = ",i3,"             ")') psi
                END IF
                WRITE (21,*) input_line(j)
            END DO
        CLOSE(21)
        IF (nsol /= 1) THEN
            WRITE(solutenumber,'(i3.3)') sol
            char_temp = 'cp input/solute.' // solutenumber // ' ' // foldername // '/input/solute.in'
            CALL SYSTEM(char_temp)
        ELSE
            char_temp = 'cp input/solute.in ' // foldername // '/input/solute.in'
            CALL SYSTEM(char_temp)
        END IF
    END DO; END DO; END DO
    IF (i /= njob_max) STOP 'FATAL ERROR: Loossing jobs.'

    OPEN(21,FILE='lancer.sh')
        DO i = njob_min, njob_max
            WRITE(foldername,'(i3.3)') i
            WRITE(21,*) 'echo ',foldername
            WRITE(21,*) 'cd ',foldername
            WRITE(21,*) 'llsubmit x.sh'
            WRITE(21,*) 'cd ..'
        END DO
    CLOSE(21)

    OPEN(21,FILE='rm.sh')
        DO i = njob_min, njob_max
            WRITE(foldername,'(i3.3)') i
            WRITE(21,*) 'rm -rf ',foldername,'/output ',foldername,'/input/ck* ',foldername,'/mdft'
        END DO
    CLOSE(21)

    PRINT*, 'MISSION COMPLETE POINCARE_JOB_LANCER'

CONTAINS
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
MODULE input

    USE precision_kinds ,ONLY: i2b,dp
    IMPLICIT NONE
    CHARACTER(len=100), ALLOCATABLE, DIMENSION(:) :: input_line ! array containing all input lines
    LOGICAL :: verbose
    PRIVATE
    PUBLIC :: verbose, input_line, input_dp, input_int, input_log, input_char, n_linesInFile, deltaAbscissa

    CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        REAL(DP) PURE FUNCTION input_dp (That)
            IMPLICIT NONE
            CHARACTER(*), INTENT(IN) :: That
            INTEGER(i2b) :: i, j
            j=LEN(That)
            DO i =1,SIZE(input_line) 
                IF( input_line( i)( 1:j) == That  .AND. input_line(i)(j+1:j+1)==' ' ) READ(input_line(i)(j+4:j+50),*) input_dp
            END DO
        END FUNCTION input_dp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
    
        INTEGER(I2B) PURE FUNCTION input_int (that)
            IMPLICIT NONE
            CHARACTER(*), INTENT(IN) :: That
            INTEGER(i2b) :: i, j
            j=LEN(That)
            DO i = 1, SIZE( input_line) 
                IF( input_line( i)( 1:j) == That  .AND. input_line(i)(j+1:j+1)==' ' ) READ(input_line(i)(j+4:j+50),*)input_int
            END DO
        END FUNCTION input_int
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
        LOGICAL FUNCTION input_log (that)
            IMPLICIT NONE
            CHARACTER(*), INTENT(IN) :: that
            CHARACTER :: text
            INTEGER(i2b) :: i, j
            j=LEN(That)
            DO i =1,SIZE( input_line) 
                IF( input_line(i)(1:j)==that .AND. input_line(i)(j+1:j+1)==' ' ) READ( input_line (i) (j+4:j+50) , * ) text
            END DO
            j = 999 ! means error in reading
            IF( text(1:1) == 'T' ) j = 1 ! means true, 2 means false
            IF( text(1:1) == 't' ) j = 1
            IF( text(1:1) == 'F' ) j = 2
            IF( text(1:1) == 'f' ) j = 2
            IF( j == 999 ) THEN
                PRINT*, 'problem in reading logical ', that
                STOP
            END IF
            IF( j == 1 ) input_log = .TRUE.
            IF( j == 2 ) input_log = .FALSE.
        END FUNCTION input_log
            
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        CHARACTER(50) FUNCTION input_char (that)
            IMPLICIT NONE
            CHARACTER(*), INTENT(IN) :: that
            INTEGER(i2b) :: i,j,imax,iostatint
            j=LEN(That)
            i=0
            imax=SIZE(input_line)
            DO i=1,imax+1
                IF (i==imax+1) THEN
                    PRINT*,"I didnt find keyword '",That,"' in dft.in"
                    STOP
                END IF
                IF (input_line(i)(1:j)==that .AND. input_line(i)(j+1:j+1)==' ') THEN
                    READ(input_line(i)(j+4:j+50),*,IOSTAT=iostatint) input_char
                        IF (iostatint/=0) THEN
                            PRINT*,"I have a problem in reading input line:"
                            PRINT*,TRIM(ADJUSTL(input_line(i)))
                            IF (TRIM(ADJUSTL(input_line(i)(j+4:j+50)))=='') PRINT*,"I found nothing after sign ="
                            STOP
                        END IF
                    EXIT
                END IF
            END DO
            IF (input_char(1:1)==' ') THEN
                PRINT*,"First character of ",that," is a whitespace"
                STOP
            END IF
            IF (LEN(TRIM(ADJUSTL(input_char)))==0) THEN
                PRINT*,"Tag after ",That," is only whitespaces."
                STOP
            END IF
        END FUNCTION input_char
            
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
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
                        n_linesInFile = n_linesInFile +1
                    END IF
                END DO
            CLOSE (77)
        END FUNCTION

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        FUNCTION deltaAbscissa (filename)
            REAL(dp) :: abscissa, previousAbscissa, ordonates, deltaAbscissa
            CHARACTER(*), INTENT(IN) :: filename
            INTEGER(i2b) :: i, ios, n_lines
            OPEN (10, FILE=filename, IOSTAT=ios)
                IF (ios /= 0) THEN
                    WRITE(*,*)"Cant open file ",filename," in FUNCTION deltaAbscissa"
                END IF
            
                DO i= 1, 2
                    READ(10,*,IOSTAT=ios) previousAbscissa, ordonates
                        IF (ios/=0) then
                            PRINT*, 'Something went wrong in reading', filename, 'in function deltaAbscissa'
                            STOP
                        END IF
                    READ(10,*, IOSTAT=ios) abscissa, ordonates
                        IF (ios/=0) then
                            PRINT*, 'Something went wrong in reading', filename, 'in function deltaAbscissa'
                            STOP
                        END IF
                END DO
            deltaAbscissa = abscissa-previousAbscissa
            CLOSE(10)
            n_lines = n_linesInFile(filename)
            OPEN (10, FILE=filename, IOSTAT=ios)
                IF (ios /= 0) THEN
                    WRITE(*,*)"Cant open file ",filename," in FUNCTION deltaAbscissa"
                END IF
            
                READ(10,*,IOSTAT=ios) abscissa, ordonates
                    IF (ios/=0) then
                        PRINT*, 'Something went wrong in reading', filename, 'in function deltaAbscissa'
                        STOP
                    END IF
                DO i=1, n_lines-1
                    previousAbscissa = abscissa
                    READ(10,*,IOSTAT=ios) abscissa, ordonates
                    IF (ios>0) then
                        PRINT*, 'Something went wrong in reading', filename, 'in function deltaAbscissa'
                        STOP
                    ELSE IF (ios<0) THEN
                        EXIT
                    ELSE
                        IF ((abscissa-previousAbscissa-deltaAbscissa)/deltaAbscissa > 1E-5) THEN
                            PRINT*, abscissa, previousAbscissa,abscissa-previousAbscissa ,deltaAbscissa
                            PRINT*, 'STOP. Non uniform absissa in ', filename
                            STOP
                        END IF
                    END IF
                END DO
            CLOSE(10)
        END FUNCTION deltaAbscissa

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE input

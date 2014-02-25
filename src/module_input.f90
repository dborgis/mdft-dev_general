MODULE input

    USE precision_kinds ,ONLY: i2b,dp
    IMPLICIT NONE
    CHARACTER(len=100), ALLOCATABLE, DIMENSION(:) :: input_line ! array containing all input lines
    LOGICAL :: verbose
    PRIVATE
    PUBLIC :: verbose, input_line, input_dp, input_int, input_log, input_char, n_linesInFile, deltaAbscissa

    CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        REAL(DP) PURE FUNCTION INPUT_DP( That)
            IMPLICIT NONE
            CHARACTER(*), INTENT(IN) :: That
            INTEGER(i2b) :: i, j
            j=LEN(That)
            DO i =1,SIZE(input_line) 
                IF( input_line( i)( 1:j) == That  .AND. input_line(i)(j+1:j+1)==' ' ) READ(input_line(i)(j+4:j+50),*) input_dp
            END DO
        END FUNCTION INPUT_DP

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
    
        INTEGER(I2B) PURE FUNCTION INPUT_INT( That)
            IMPLICIT NONE
            CHARACTER(*), INTENT(IN) :: That
            INTEGER(i2b) :: i, j
            j=LEN(That)
            DO i = 1, SIZE( input_line) 
                IF( input_line( i)( 1:j) == That  .AND. input_line(i)(j+1:j+1)==' ' ) READ(input_line(i)(j+4:j+50),*)input_int
            END DO
        END FUNCTION INPUT_INT
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
        LOGICAL FUNCTION INPUT_LOG( That)
            IMPLICIT NONE
            CHARACTER(*), INTENT(IN) :: That
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
        END FUNCTION INPUT_LOG
            
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
        CHARACTER(50) PURE FUNCTION INPUT_CHAR( That)
            IMPLICIT NONE
            CHARACTER(*), INTENT(IN) :: That
            INTEGER(i2b) :: i, j
            j=LEN(That)
            DO i = 1,SIZE( input_line) 
                IF( input_line( i)( 1:j) == That  .AND. input_line(i)(j+1:j+1)==' ' ) READ ( input_line (i) (j+4:j+50),*)input_char
            END DO
        END FUNCTION INPUT_CHAR
            
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

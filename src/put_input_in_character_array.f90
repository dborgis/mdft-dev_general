SUBROUTINE put_input_in_character_array
    
    U precision_kinds , only : i2b
    use input , only : input_line , n_linesInFile
    implicit none
    integer(i2b):: i , j , k , n, n_lines ! dummy
    integer(i2b):: ios ! input output status of readen file
    character ( len = 100 ) :: text ! temporary input line
    character ( len = 100 ) , allocatable , dimension ( : ):: arraytemp  ! Temporary array to stock data for resizing input_line

    n_lines = n_linesInFile('input/dft.in')
    ALLOCATE ( input_line (n_lines) )
    OPEN (11, FILE='input/dft.in' )
        DO i=1, n_lines
            READ(11,'(a)') text
            input_line (i) = trim(adjustl(text))
        END DO
    CLOSE (11)
    
    !  clean up comments in the lines (expl: option = 3 # blabla)
    DO i = 1, n_lines
        DO j = 1 , len(text)
            IF ( input_line (i) (j:j) == '#' ) THEN
                DO CONCURRENT ( k=j:LEN(text) )
                    input_line(i)(k:k) = ' '
                END DO
                EXIT
            END IF
        END DO
        input_line (i) = TRIM( ADJUSTL( input_line(i) ))
    END DO
    
    !Delete blank lines and count the size of the smallest array containing initial data  
    n = 0
    do i = 1 , n_lines
        if ( input_line (i) (1:1) /= ' ' )  then
            input_line (n+1) = input_line(i)
            n = n + 1   
        endif
    end do
    
    !Resize input_line to the smallest size by using a temporary array
    allocate ( arraytemp(n) , SOURCE= input_line(1:n) )
    deallocate ( input_line )
    allocate ( input_line (n), SOURCE= arraytemp )
    deallocate (arraytemp)

end subroutine put_input_in_character_array

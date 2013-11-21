! This SUBROUTINE prints all the input parameters in output/input.out
! so that all files needed to understand the outputs are in the output folder.

SUBROUTINE print_input_to_output_folder

    USE input, ONLY: input_line
    USE precision_kinds, ONLY: i2b
    
    IMPLICIT NONE
    
    INTEGER(i2b):: i ! dummy for loop
    
    OPEN (10, file='output/dft.in.out' )
        do i = 1 , size( input_line ) ! print each line of input_line()
            write (10,*) input_line (i)
        END DO
    CLOSE (10)

END SUBROUTINE print_input_to_output_folder

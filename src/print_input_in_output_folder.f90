! This SUBROUTINE prints all the input parameters in output/input.out
! so that all files needed to understand the outputs are in the output folder.

SUBROUTINE print_input_to_output_folder

    USE input, ONLY: input_line
    USE precision_kinds, ONLY: i2b
    
    IMPLICIT NONE
    
    INTEGER(i2b):: i ! dummy for loop

    CALL test_if_output_folder_exists_and_create_one_if_necessary

    OPEN (10, FILE='output/dft.in.out' )
        DO i = 1 , SIZE( input_line ) ! print each line of input_line()
            WRITE (10,*) input_line (i)
        END DO
    CLOSE (10)


    CONTAINS
    
    SUBROUTINE  test_if_output_folder_exists_and_create_one_if_necessary
        CALL system('mkdir -p output') ! just create folder. If it already exists, nothing happens.
    END SUBROUTINE test_if_output_folder_exists_and_create_one_if_necessary

END SUBROUTINE print_input_to_output_folder

!> this SUBROUTINE prints all the input parameters in output/input.out
! so that all files needed to understand the outputs are in the output folder.
! written by Maximilien Levesque, 2011, @ Ecole Normale Superieure
SUBROUTINE print_input_to_output_folder
  use input , only : input_line
  USE precision_kinds , only : i2b
  IMPLICIT NONE
  integer(i2b):: i ! dummy for loop
 
  ! open file to write in
  open( unit = 10 , file='output/input.out' )
  ! print each line of input_line()
  do i = 1 , size( input_line )
    write (10,*) input_line (i)
  END DO
END SUBROUTINE print_input_to_output_folder

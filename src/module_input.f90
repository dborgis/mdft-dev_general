module input
  use precision_kinds , only : i2b, dp
  implicit none
  character (len = 100) , allocatable , dimension (:) :: input_line ! array containing all input lines
  integer(i2b):: TotalNumberOfInputLines
  private
  public :: input_line, input_dp, input_int, TotalNumberOfInputLines,input_log, input_char
  contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL(DP) PURE FUNCTION INPUT_DP( That)
  implicit none
  character(*), intent(in) :: That
  integer(i2b) :: i, j
  j=len(That)
  do i = 1, size( input_line) 
    if( input_line( i)( 1:j) == That ) read ( input_line (i) (j+4:j+50) , * ) input_dp
  end do
END FUNCTION INPUT_DP
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
INTEGER(I2B) PURE FUNCTION INPUT_INT( That)
  implicit none
  character(*), intent(in) :: That
  integer(i2b) :: i, j
  j=len(That)
  do i = 1, size( input_line) 
    if( input_line( i)( 1:j) == That ) read ( input_line (i) (j+4:j+50) , * ) input_int
  end do
END FUNCTION INPUT_INT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
LOGICAL FUNCTION INPUT_LOG( That)
  implicit none
  character(*), intent(in) :: That
  character :: text ! dummy
  integer(i2b) :: i, j
  j=len(That)
  do  i = 1, size( input_line) 
    if( input_line( i)( 1:j) == That ) read ( input_line (i) (j+4:j+50) , * ) text
  end do
  j = 999 ! means error in reading
  if( text(1:1) == 'T' ) j = 1 ! means true, 2 means false
  if( text(1:1) == 't' ) j = 1
  if( text(1:1) == 'F' ) j = 2
  if( text(1:1) == 'f' ) j = 2
  if( j == 999 ) then
    print*, 'problem in reading logical ', that
    stop
  end if
  if( j == 1 ) input_log = .true.
  if( j == 2 ) input_log = .false.
END FUNCTION INPUT_LOG
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CHARACTER(50) PURE FUNCTION INPUT_CHAR( That)
  implicit none
  character(*), intent(in) :: That
  integer(i2b) :: i, j
  j=len(That)
  do i = 1, size( input_line) 
    if( input_line( i)( 1:j) == That ) read ( input_line (i) (j+4:j+50) , * ) input_char
  end do
END FUNCTION INPUT_CHAR
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module input

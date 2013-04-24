subroutine put_input_in_character_array


use precision_kinds , only : i2b

use input , only : input_line , totalnumberofinputlines




implicit none

integer ( kind = i2b ) :: i , j , k , n ! dummy

integer ( kind = i2b ) :: ios ! input output status of readen file

character ( len = 100 ) :: text ! temporary input line

character ( len = 100 ) , allocatable , dimension ( : ):: arraytemp  ! Temporary array to stock data for resizing input_line





! open the file in which you have all inputs

open ( unit = 11 , file = 'input/dft.in' )




! count the total number of lines (n) in input.in
! a blank line counts for 0 if it is only a return to next line (\n)
! a line with only blanks counts for 1

totalnumberofinputlines = 0

do while (.true.)

  read ( 11 , '(a)' , iostat = ios ) text

  if ( ios > 0 ) then

    write (*,*) 'error in input. critical'

    stop

  else if ( ios < 0 ) then

    exit ! exit do loop and stop reading file

  else

    if ( text /= '\n' ) totalnumberofinputlines = totalnumberofinputlines + 1

  end if

end do




! close file in order to open it properly from start

close ( unit = 11 )
  


! allocate input_line so that it has the right number of lines

allocate ( input_line (totalnumberofinputlines) )

input_line = ' ' ! init


! open again the input file

open ( unit = 11 , file = 'input/dft.in' )



!read it and put each line in input_line

100 format(a)

do i = 1 , totalnumberofinputlines

  read ( 11 , '(a)'  ) text

  input_line (i) = trim ( adjustl (text) ) ! trim () removes trailing blanks while adjustl removes left blanks and put white blanks at the end

end do 


n=0 !init

!  clean up comments in the lines (expl: option = 3 # blabla)
  do i = 1 , totalnumberofinputlines
    do j = 1 , len(text)
      if ( input_line (i) (j:j) == '#' ) then
        forall ( k = j : len(text) )
           input_line (i) (k:k) = ' '
        end forall
        exit
      end if
    end do
    input_line (i) = trim ( adjustl (input_line (i) ) )
  end do
  

!Delete blank lines and count the size of the smallest array containing initial data  
  do i = 1 , totalnumberofinputlines
    if ( input_line (i) (1:1) /= ' ' )  then
      input_line (n+1) = input_line(i)
      n = n + 1   
    endif
  end do

!Resize input_line to the smallest size by using a temporary array

arraytemp = input_line

deallocate ( input_line )

allocate ( input_line ( n ) )

input_line = arraytemp ( 1 : n  )





end subroutine put_input_in_character_array

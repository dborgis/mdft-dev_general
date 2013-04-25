! In this subroutine, one prints the first lines of codes to stdout.
subroutine print_header
character(8)  :: date
character(10) :: time
call date_and_time ( DATE=date,TIME=time)
print*,' ******************************************************************'
print*,' ******************************************************************'
print*,' **             *************************   ',date(1:4),'/',date(5:6),'/',date(7:8),'   **********'
print*,' **   M D F T   *************************                **********'
print*,' **             *************************    ',time(1:2),'h',time(3:4),'m',time(5:6),'    **********'
print*,' ******************************************************************'
print*,' ******************************************************************'
end subroutine print_header

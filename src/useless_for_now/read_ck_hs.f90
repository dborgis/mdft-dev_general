! This SUBROUTINEs reads files containing the direct correlation functions of the first three rotational invariants
! in k-space. This files are input/cs.in input/cdelta.in and input/cd.in
SUBROUTINE read_ck_hs
!TODO
USE precision_kinds, only: i2b , dp
use system
use constants
use quadrature
IMPLICIT NONE
integer(i2b) :: nk
real(dp) :: norm_k
integer(i2b) :: ios ! iostat of the read statement: 
!    If the value of IOstatus is zero, the previous READ was executed flawlessly and all variables have received their input values. This is the normal case.
!    If the value of IOstatus is positive, the previous READ has encountered some problem. In general, without knowing the system dependent information, it is impossible to determine what the problem was. However, if hardware and I/O devices are working, a commonly seen problem would be illegal data. For example, supplying a real number to an integer variable. If IOstatus is positive, you cannot trust the values of the variables in the READ statement; they could all contain garbage values, or some of them are fine while the others are garbage.
!    If the value of IOstatus is negative, it means the end of the input has reached. Under this circumstance, some or all of the variables in the READ may not receive input values. 
real(dp):: value1 , value2 ! dummy
!> warn user
write(*,*)'>> read direct correlation functions from input/cs_hs.in'
!> read the total number of lines in input/cs.in (which is the same as in input/cd.in and input/cdelta.in
open(11,file='input/cs_hs.in')
nb_k=0
do while(.true.)
  read(11,*,iostat=ios)
  if (ios>0) then
    write(*,*)'Error in compute_ck_dipolar.f90'
    write(*,*)'something went wrong during the computation of the total number of lines in cs.in. stop'
    stop
  ELSE IF (ios<0) then
    ! end of file reached
    exit
  ELSE
    nb_k=nb_k+1
  END IF
END DO
close(11)
write(*,*)'nb_k = ',nb_k
!> read the distance between two k points in input/cs.in  (which is the same as in input/cd.in and input/cdelta.in
open(11,file='input/cs_hs.in')
read(11,*)value1, norm_k
read(11,*)value2, norm_k
close(11)
delta_k=value2-value1
write(*,*)'delta_k in cs_hs.in = ',delta_k
!> read and open the direct correlation functions previously calculated
open(11,file='input/cs_hs.in')
! Now that we know the total number of k points, allocate arrays accordingly
allocate(c_s_hs(nb_k))
do nk=1,nb_k
  read(11,*,iostat=ios) norm_k, c_s_hs(nk) !partie spherique
  if (ios>0 .or. ios<0) then
    write(*,*)'Error in read_ck_hs.f90'
    write(*,*)'something went wrong during reading c_s_hs. stop'
    stop
  END IF
END DO
! close files
close(11)
END SUBROUTINE read_ck_hs
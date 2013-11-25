! This SUBROUTINEs reads files containing the direct correlation functions of the first three rotational invariants
! in k-space. This files are input/cs.in input/cdelta.in and input/cd.in
SUBROUTINE read_ck
!~ 
!~   USE precision_kinds,only: i2b,dp
!~   use system
!~   use constants
!~   use quadrature
!~   use input,only : input_line , input_char, input_log, n_linesInFile
!~ 
!~   IMPLICIT NONE
!~ 
!~   integer(i2b) :: nk, i, j
!~   character(len=5) :: ck_species
!~   real(dp) :: norm_k
!~   integer(i2b) :: ios ! iostat of the read statement: 
!~   real(dp):: value1 , value2 ! dummy
!~   ! If the value of IOstatus is zero, the previous READ was executed flawlessly and all variables have received their input values. This is the normal case.
!~   ! If the value of IOstatus is positive, the previous READ has encountered some problem. In general, without knowing the system dependent information, it is impossible to determine what the problem was. However, if hardware and I/O devices are working, a commonly seen problem would be illegal data. For example, supplying a real number to an integer variable. If IOstatus is positive, you cannot trust the values of the variables in the READ statement; they could all contain garbage values, or some of them are fine while the others are garbage.
!~   ! If the value of IOstatus is negative, it means the end of the input has reached. Under this circumstance, some or all of the variables in the READ may not receive input values. 
!~   ! Get the correct input files depending on the solvent.
!~ 
!~    ck_species=trim(adjustl(input_char('ck_species')))
!~   ! read the total number of lines in input/cs.in (which is the same as in input/cd.in and input/cdelta.in
!~   if ( ck_species == 'spc  ' ) then
!~     open(11,file='input/direct_correlation_functions/water/SPC_Lionel_Daniel/cs.in')
!~   ELSE IF ( ck_species == 'stock' ) then
!~     open(11,file='input/direct_correlation_functions/stockmayer/cs.in')
!~   ELSE IF ( ck_species == 'perso' ) then
!~     open(11,file='input/cs.in')
!~   ELSE IF ( ck_species == 'spce' ) then
!~     open(11,file='input/direct_correlation_functions/water/SPCE/cs.in')
!~   ELSE
!~     stop 'in read_ck.f90 ck_species should be spc or sto or per'
!~   END IF
!~   nb_k=0
!~   do while(.true.)
!~     read(11,*,iostat=ios)
!~     if (ios>0) then
!~       write(*,*)'Error in compute_ck_dipolar.f90'
!~       write(*,*)'something went wrong during the computation of the total number of lines in cs.in. stop'
!~       stop
!~     ELSE IF (ios<0) then
!~       ! end of file reached
!~       exit
!~     ELSE
!~       nb_k=nb_k+1
!~     END IF
!~   END DO
!~   close(11)
!~   
!~   
!~   
!~   ! read the distance between two k points in input/cs.in  (which is the same as in input/cd.in and input/cdelta.in
!~   if ( ck_species == 'spc  ' ) then
!~     open(11,file='input/direct_correlation_functions/water/SPC_Lionel_Daniel/cs.in')
!~   ELSE IF ( ck_species == 'stock' ) then
!~     open(11,file='input/direct_correlation_functions/stockmayer/cs.in')
!~   ELSE IF ( ck_species == 'perso' ) then
!~     open(11,file='input/cs.in')
!~   ELSE IF ( ck_species == 'spce' ) then
!~     open(11,file='input/direct_correlation_functions/water/SPCE/cs.in')
!~   END IF
!~   read(11,*)value1, norm_k
!~   read(11,*)value2, norm_k
!~   close(11)
!~   delta_k=value2-value1
!~ 
!~   ! read and open the direct correlation functions previously calculated
!~   if ( ck_species == 'spc  ' ) then
!~     open(11,file='input/direct_correlation_functions/water/SPC_Lionel_Daniel/cs.in')
!~     open(12,file='input/direct_correlation_functions/water/SPC_Lionel_Daniel/cdelta.in')
!~     open(13,file='input/direct_correlation_functions/water/SPC_Lionel_Daniel/cd.in')
!~   ELSE IF ( ck_species == 'stock' ) then
!~     open(11,file='input/direct_correlation_functions/stockmayer/cs.in')
!~     open(12,file='input/direct_correlation_functions/stockmayer/cdelta.in')
!~     open(13,file='input/direct_correlation_functions/stockmayer/cd.in')
!~   ELSE IF ( ck_species == 'perso' ) then
!~     open(11,file='input/cs.in')
!~     open(12,file='input/cdelta.in')
!~     open(13,file='input/cd.in')
!~   ELSE IF ( ck_species == 'spce' ) then
!~     open(11,file='input/direct_correlation_functions/water/SPCE/cs.in')
!~     open(12,file='input/direct_correlation_functions/water/SPCE/cdelta.in')
!~     open(13,file='input/direct_correlation_functions/water/SPCE/cd.in')
!~   END IF
!~   ! Now that we know the total number of k points, allocate arrays accordingly
!~   allocate(c_s(nb_k))
!~   allocate(c_delta(nb_k))
!~   allocate(c_d(nb_k))
!~   do nk=1,nb_k
!~     read(11,*,iostat=ios) norm_k, c_s(nk) !partie spherique
!~     if (ios>0 .or. ios<0) then
!~       write(*,*)'Error in compute_ck_dipolar.f90'
!~       write(*,*)'something went wrong during reading c_s. stop'
!~       stop
!~     END IF
!~     read(12,*,iostat=ios) value1, c_delta(nk) !110
!~     if (ios>0) then
!~       write(*,*)'Error in compute_ck_dipolar.f90'
!~       write(*,*)'something went wrong during reading c_delta. stop'
!~       stop
!~     ELSE IF (ios<0) then     ! end of file reached
!~       c_delta(nk)=c_delta(nk-1)*0.5_dp ! TODO MAGIC NUMBER IN ORDER TO PUSH IT TO DECREASE TOWARD ZERO
!~     END IF
!~     read(13,*,iostat=ios) value1, c_d(nk) !211
!~     if (ios>0) then
!~       write(*,*)'Error in compute_ck_dipolar.f90'
!~       write(*,*)'something went wrong during reading c_d. stop'
!~       stop
!~     ELSE IF (ios<0) then     ! end of file reached
!~       c_d(nk)=c_d(nk-1)*0.5_dp ! TODO MAGIC NUMBER IN ORDER TO PUSH IT TO DECREASE TOWARD ZERO
!~     END IF
!~   
!~   END DO
!~   
!~   ! close files
!~   close(11)
!~   close(12)
!~   close(13)
!~ 
!~ !> read the distance between two k points in input/cs.in  (which is the same as in input/cd.in and input/cdelta.in
!~ if (input_log('read_chi')) then
!~ 
!~     open(11,file='input/direct_correlation_functions/water/chi_SPCE_for_multi/chi_t_multi.in')
!~         read(11,*)value1, norm_k
!~         read(11,*)value2, norm_k
!~     close(11)
!~     delta_k=value2-value1
!~     
!~     !> read and open the direct correlation functions previously calculated
!~     OPEN (14, file='input/direct_correlation_functions/water/chi_SPCE_for_multi/chi_L_SPCE_splinetronque.dat')        !input/chi_L.in chi_L_SPCE_splinetronque.dat
!~     OPEN (15, file='input/direct_correlation_functions/water/chi_SPCE_for_multi/chi_t_multi.in')        !chi_t.in
!~         ALLOCATE (chi_l(nb_k))
!~         ALLOCATE (chi_t(nb_k))
!~         DO nk=1,nb_k
!~             READ (14,*,IOSTAT=ios) norm_k, chi_L(nk)
!~                 IF (ios>0) then
!~                     WRITE(*,*)'Error in input/direct_correlation_functions/water/chi_SPCE_for_multi/chi_L_SPCE_splinetronque.dat'
!~                     STOP
!~                 END IF
!~             READ (15,*,IOSTAT=ios) norm_k, chi_t(nk)
!~                 IF (ios>0) then
!~                     WRITE(*,*)'Error in input/direct_correlation_functions/water/chi_SPCE_for_multi/chi_t_multi.in'
!~                     STOP
!~                 END IF
!~         END DO
!~     CLOSE (14)
!~     CLOSE (15)
!~ END IF


END SUBROUTINE read_ck

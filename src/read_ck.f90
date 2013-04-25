! This subroutines reads files containing the direct correlation functions of the first three rotational invariants
! in k-space. This files are input/cs.in input/cdelta.in and input/cd.in
subroutine read_ck
  use precision_kinds , only: i2b,dp
  use system
  use constants
  use quadrature
  use input , only : input_line , input_char, input_log! contains all input file dft.in in a character array
  implicit none
  integer(i2b) :: nk, i, j
  character(len=5) :: ck_species
  real(dp) :: norm_k
  integer(i2b) :: ios ! iostat of the read statement: 
  real(dp):: value1 , value2 ! dummy
  ! If the value of IOstatus is zero, the previous READ was executed flawlessly and all variables have received their input values. This is the normal case.
  ! If the value of IOstatus is positive, the previous READ has encountered some problem. In general, without knowing the system dependent information, it is impossible to determine what the problem was. However, if hardware and I/O devices are working, a commonly seen problem would be illegal data. For example, supplying a real number to an integer variable. If IOstatus is positive, you cannot trust the values of the variables in the READ statement; they could all contain garbage values, or some of them are fine while the others are garbage.
  ! If the value of IOstatus is negative, it means the end of the input has reached. Under this circumstance, some or all of the variables in the READ may not receive input values. 
  ! Get the correct input files depending on the solvent.
   ck_species=trim(adjustl(input_char('ck_species')))
  ! read the total number of lines in input/cs.in (which is the same as in input/cd.in and input/cdelta.in
  if ( ck_species == 'spc  ' ) then
    open(11,file='input/direct_correlation_functions/water/SPC_Lionel_Daniel/cs.in')
  else if ( ck_species == 'stock' ) then
    open(11,file='input/direct_correlation_functions/stockmayer/cs.in')
  else if ( ck_species == 'perso' ) then
    open(11,file='input/cs.in')
  else if ( ck_species == 'spce' ) then
    open(11,file='input/direct_correlation_functions/water/SPCE/cs.in')
  else
    stop 'in read_ck.f90 ck_species should be spc or sto or per'
  end if
  nb_k=0
  do while(.true.)
    read(11,*,iostat=ios)
    if (ios>0) then
      write(*,*)'Error in compute_ck_dipolar.f90'
      write(*,*)'something went wrong during the computation of the total number of lines in cs.in. stop'
      stop
    else if (ios<0) then
      ! end of file reached
      exit
    else
      nb_k=nb_k+1
    end if
  end do
  close(11)
  ! read the distance between two k points in input/cs.in  (which is the same as in input/cd.in and input/cdelta.in
  if ( ck_species == 'spc  ' ) then
    open(11,file='input/direct_correlation_functions/water/SPC_Lionel_Daniel/cs.in')
  else if ( ck_species == 'stock' ) then
    open(11,file='input/direct_correlation_functions/stockmayer/cs.in')
  else if ( ck_species == 'perso' ) then
    open(11,file='input/cs.in')
  else if ( ck_species == 'spce' ) then
    open(11,file='input/direct_correlation_functions/water/SPCE/cs.in')
  end if
  read(11,*)value1, norm_k
  read(11,*)value2, norm_k
  close(11)
  delta_k=value2-value1
  ! read and open the direct correlation functions previously calculated
  if ( ck_species == 'spc  ' ) then
    open(11,file='input/direct_correlation_functions/water/SPC_Lionel_Daniel/cs.in')
    open(12,file='input/direct_correlation_functions/water/SPC_Lionel_Daniel/cdelta.in')
    open(13,file='input/direct_correlation_functions/water/SPC_Lionel_Daniel/cd.in')
  else if ( ck_species == 'stock' ) then
    open(11,file='input/direct_correlation_functions/stockmayer/cs.in')
    open(12,file='input/direct_correlation_functions/stockmayer/cdelta.in')
    open(13,file='input/direct_correlation_functions/stockmayer/cd.in')
  else if ( ck_species == 'perso' ) then
    open(11,file='input/cs.in')
    open(12,file='input/cdelta.in')
    open(13,file='input/cd.in')
  else if ( ck_species == 'spce' ) then
    open(11,file='input/direct_correlation_functions/water/SPCE/cs.in')
    open(12,file='input/direct_correlation_functions/water/SPCE/cdelta.in')
    open(13,file='input/direct_correlation_functions/water/SPCE/cd.in')
  end if
  ! Now that we know the total number of k points, allocate arrays accordingly
  allocate(c_s(nb_k))
  allocate(c_delta(nb_k))
  allocate(c_d(nb_k))
  do nk=1,nb_k
    read(11,*,iostat=ios) norm_k, c_s(nk) !partie spherique
    if (ios>0 .or. ios<0) then
      write(*,*)'Error in compute_ck_dipolar.f90'
      write(*,*)'something went wrong during reading c_s. stop'
      stop
    end if
    read(12,*,iostat=ios) value1, c_delta(nk) !110
    if (ios>0) then
      write(*,*)'Error in compute_ck_dipolar.f90'
      write(*,*)'something went wrong during reading c_delta. stop'
      stop
    else if (ios<0) then     ! end of file reached
      c_delta(nk)=c_delta(nk-1)*0.5_dp ! TODO MAGIC NUMBER IN ORDER TO PUSH IT TO DECREASE TOWARD ZERO
    end if
    read(13,*,iostat=ios) value1, c_d(nk) !211
    if (ios>0) then
      write(*,*)'Error in compute_ck_dipolar.f90'
      write(*,*)'something went wrong during reading c_d. stop'
      stop
    else if (ios<0) then     ! end of file reached
      c_d(nk)=c_d(nk-1)*0.5_dp ! TODO MAGIC NUMBER IN ORDER TO PUSH IT TO DECREASE TOWARD ZERO
    end if
  
  end do
  
  ! close files
  close(11)
  close(12)
  close(13)
!> read the distance between two k points in input/cs.in  (which is the same as in input/cd.in and input/cdelta.in
if (input_log('read_chi')) then
open(11,file='input/direct_correlation_functions/water/chi_SPCE_for_multi/chi_t_multi.in')
read(11,*)value1, norm_k
read(11,*)value2, norm_k
close(11)
delta_k=value2-value1
!> read and open the direct correlation functions previously calculated
open(14,file='input/direct_correlation_functions/water/chi_SPCE_for_multi/chi_L_SPCE_splinetronque.dat ')        !input/chi_L.in chi_L_SPCE_splinetronque.dat
open(15,file='input/direct_correlation_functions/water/chi_SPCE_for_multi/chi_t_multi.in')        !chi_t.in
! Now that we know the total number of k points, allocate arrays accordingly
allocate(chi_l(nb_k))
allocate(chi_t(nb_k))
do nk=1,nb_k
  read(14,*,iostat=ios) norm_k, chi_L(nk) 
  read(15,*,iostat=ios) norm_k, chi_t(nk) 
end do
end if
! close files
close(11)
close(12)
close(13)
close(14)
close(15)
end subroutine read_ck

program format

  implicit none
  integer, parameter :: mmax=3, np=252
  integer, dimension(np) :: p, m, n, mu, nu, khi
  integer :: i, iq, ip ! interators
  character(10) :: sometext
  character(80) :: outfilename
  integer, parameter :: someunit=11, anotherunit=12
  integer, parameter :: nq=1024
  real, dimension(nq) :: q
  complex, dimension(np,nq) :: c
  real, dimension(np,nq) :: l2_abs_c


  open(unit=someunit,file="ck_nonzero_nmax5_ml_lu")

  do i=1,10
    read(someunit,*)
  end do
  read(someunit,*) sometext, p(:)
  read(someunit,*) sometext, m(:)
  read(someunit,*) sometext, n(:)
  read(someunit,*) sometext, mu(:)
  read(someunit,*) sometext, nu(:)
  read(someunit,*) sometext, khi(:)
  read(someunit,*)

  print*, "p  ", p(:)
  print*, "m  ", m(:)
  print*, "n  ", n(:)
  print*, "mu ", mu(:)
  print*, "nu ", nu(:)
  print*, "khi", khi(:)


  do iq=1,nq
    read(someunit,*) q(iq), c(:,iq)
  end do
  close(someunit)


  do ip=1,np

    ! if (any(abs(c(ip,:))>=1)) then
      write(outfilename,'(a,i0,a,i0,a,i0,a,i0,a,i0,a)') "cmnmunukhi_",m(ip),"_",n(ip),"_",mu(ip),"_",nu(ip),"_",khi(ip),".dat"
      open(unit=someunit,file=trim(adjustl(outfilename)))
      do iq=1,nq
        write(someunit,*) q(iq), abs(c(ip,iq))
      end do
      close(someunit)
    ! end if

  end do



end program format

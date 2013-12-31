!> Compute planar density or planar polarization
!! For now only works with plans in x, y or z direction
SUBROUTINE compute_planar_density(array,filename)
USE precision_kinds, only: dp
use system, only: nfft1,nfft2,nfft3,nb_solute_sites,x_mol,y_mol,z_mol,Lx,Ly,Lz
IMPLICIT NONE
 character(50),intent(in):: filename
real(dp),intent(in),dimension(nfft1,nfft2,nfft3):: array
integer:: plandir !> 1=yz 2=xz 3=xy
integer:: id !> grid index
integer:: i,j,k
real(dp),dimension(nfft1):: x_com !> Cartesian coordinates of grid points in i direction
real(dp),dimension(nfft2):: y_com !> Cartesian coordinates of grid points in j direction
real(dp),dimension(nfft3):: z_com !> Cartesian coordinates of grid points in k direction
!> Identify the plan coordinate which is 0 (see restrictions to using this program for now)
if (x_mol(1)==x_mol(2) .and. x_mol(1)==x_mol(3)) then
  plandir=1
ELSE IF (y_mol(1)==y_mol(2) .and. y_mol(1)==y_mol(3)) then
  plandir=2
ELSE IF (z_mol(1)==z_mol(2) .and. z_mol(1)==z_mol(3)) then
  plandir=3
ELSE
  goto 777
END IF
print*,'lalalalal'
!> Get its grid index called id
if (plandir==1) then
  id= nint(x_mol(1)*real(nfft1,dp)/Lx) +1
ELSE IF (plandir==2) then
  id= nint(y_mol(1)*real(nfft2,dp)/Ly) +1
ELSE IF (plandir==3) then
  id= nint(z_mol(1)*real(nfft3,dp)/Lz) +1
ELSE
  write(*,*)'error in compute_planar_density'
END IF
!> Compute grid points cartesian coordinates
forall(i=1:nfft1) x_com(i)=real(i-1,dp)*Lx/real(nfft1,dp)
forall(j=1:nfft2) y_com(j)=real(j-1,dp)*Ly/real(nfft2,dp)
forall(k=1:nfft3) z_com(k)=real(k-1,dp)*Lz/real(nfft3,dp)
!> Print density in this plan
open(10,file=filename,form='formatted')
100 format (3(xF10.5))
write(10,*)'# xn yn density'
if (plandir==1) then
   do j=1,nfft2 ; do k=1,nfft3
        write(10,100)y_com(j),z_com(k),array(52,j,k)
   END DO ; END DO
ELSE IF (plandir==2) then
   do i=1,nfft1 ; do k=1,nfft3
        write(10,100)x_com(i),z_com(k),array(i,id,k)
   END DO ; END DO
ELSE IF (plandir==3) then
   do i=1,nfft1 ; do j=1,nfft2
        write(10,100)x_com(i),y_com(j),array(i,j,id)
   END DO ; END DO
END IF
 close(10)
777 continue
END SUBROUTINE

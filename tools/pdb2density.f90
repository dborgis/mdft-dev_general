
! BEFORE RUNNING THIS PROGRAM, YOU SHOULD RUN IN YOUR SHELL:
! grep OW prod-short.pdb | awk '{print $6" "$7" "$8}' > xyz.dat
! this parses the pdb file: looks for lines containing the information
! about the Water Oxygen atom (named OW), then prints x y and z only.
! nothing else, nothing more.
! Each line of file xyz.dat contains the x, y and z of one Oxygen.
!
! Then you should do :
! wc -l xyz.dat
! that will count for the number of lines (and thus oxygen atoms) to
! consider. You should report this number in the integer parameter nl
! just below.

program pdb2density

implicit none
integer, parameter :: nl=18783756
integer :: ix, iy, iz, nx, ny, nz, i
real, parameter :: dx=1.,dy=dx,dz=dx
real :: x(nl), y(nl), z(nl), xmin, xmax, ymin, ymax, zmin, zmax, lx, ly, lz
real :: bxmin, bxmax, bymin, bymax, bzmin, bzmax, bx, by, bz

open(99,file="xyz.dat",action="read")
do i=1,nl
  read(99,*)x(i),y(i),z(i)
end do
close(99)
xmin=minval(x)
ymin=minval(y)
zmin=minval(z)
xmax=maxval(x)
ymax=maxval(y)
zmax=maxval(z)
lx=xmax-xmin
ly=ymax-ymin
lz=zmax-zmin
if(lx<dx) stop "lx<dx"
if(ly<dy) stop "ly<dy"
if(lz<dz) stop "lz<dz"
nx=lx/dx
ny=ly/dy
nz=lz/dz
do ix=1,nx
  bxmin=(ix-1)*dx
  bxmax=ix*dx
  bx=(bxmin+bxmax)/2.
  if( count( x>bxmin .and. x<bxmax )==0 ) cycle
  do iy=1,ny
    bymin=(iy-1)*dy
    bymax=iy*dy
    by=(bymin+bymax)/2.
    if( count( y>bymin .and. y<bymax )==0 ) cycle
    do iz=1,nz
      bzmin=(iz-1)*dz
      bzmax=iz*dz
      bz=(bzmin+bzmax)/2.
      if( count( z>bzmin .and. z<bzmax )==0 ) cycle
      !if( norm2([bx,by,bz]-[lx,ly,lz]/2.) >15 .or. bz>lz/2 ) then 
      if( bx>lx/2 .or. bz>lz/2 ) cycle
      write(*,'(F5.1,F5.1,F5.1,I5)') bx, by, bz,&
count( x>bxmin .and. x<bxmax .and. y>bymin .and. y<bymax .and. z>bzmin .and. z<bzmax )
    end do
  end do
end do
end program pdb2density

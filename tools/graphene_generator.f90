program ggen
implicit none
integer :: i, j, n
real :: y
integer, parameter :: jmax=5
integer, parameter :: imax=5
real, parameter :: q=0.
real, parameter :: sigma=3.166
real, parameter:: epsilon=0.65




real, parameter :: a=1.42, c=2.46
print*,
print*, 60
print*,
n=1
do i=0,8,2 ! 1.42*8=11.36

print*, n, q, sigma, epsilon, a*0., i*c/2., 0., 12; n=n+1
print*, n, q, sigma, epsilon, a*1., i*c/2., 0., 12; n=n+1
print*, n, q, sigma, epsilon, a*3., i*c/2., 0., 12; n=n+1
print*, n, q, sigma, epsilon, a*4., i*c/2., 0., 12; n=n+1
print*, n, q, sigma, epsilon, a*6., i*c/2., 0., 12; n=n+1
print*, n, q, sigma, epsilon, a*7., i*c/2., 0., 12; n=n+1

print*, n, q, sigma, epsilon, a*(1. + .5 ), (i+1)*c/2., 0., 12 ; n=n+1
print*, n, q, sigma, epsilon, a*(2. + .5 ), (i+1)*c/2., 0., 12 ; n=n+1
print*, n, q, sigma, epsilon, a*(4. + .5 ), (i+1)*c/2., 0., 12 ; n=n+1
print*, n, q, sigma, epsilon, a*(5. + .5 ), (i+1)*c/2., 0., 12 ; n=n+1
print*, n, q, sigma, epsilon, a*(7. + .5 ), (i+1)*c/2., 0., 12 ; n=n+1
print*, n, q, sigma, epsilon, a*(8. + .5 ), (i+1)*c/2., 0., 12 ; n=n+1

end do



stop
open(12,file="graphene.solute.in")


write(12,*) "graphene max"
n=0
do j=1,jmax
y=0.0






do i=0,imax
if(modulo(i,3)==0 .or. modulo(i,3)==1) then
  n=n+1
  write(12,*) n, q, sigma, epsilon, i*1.42, y, 0., 12
  n=n+1
  write(12,*) n, q, sigma, epsilon, (i+1.5)*1.42, y+1.42, 0., 12
  print*, i*1.42, y
  print*, (i+1.5)*1.42, y+1.42
end if
end do
y=y+2.46
end do




close(12)
end program ggen

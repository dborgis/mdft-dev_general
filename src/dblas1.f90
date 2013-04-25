function dasum ( n, x, incx )
!
!*******************************************************************************
!
!! DASUM takes the sum of the absolute values of a vector.
!
!
!  Modified:
!
!    15 February 2001
!
!  Reference:
!
!    Lawson, Hanson, Kincaid, Krogh,
!    Basic Linear Algebra Subprograms for Fortran Usage,
!    Algorithm 539, 
!    ACM Transactions on Mathematical Software, 
!    Volume 5, Number 3, September 1979, pages 308-323.
!
!  Parameters:
!
!    Input, integer N, the number of entries in the vector.
!
!    Input, real X(*), the vector to be examined.
!
!    Input, integer INCX, the increment between successive entries of X.
!    INCX must not be negative.
!
!    Output, double precision DASUM, the sum of the absolute values of X.
!
  double precision dasum
  integer incx
  integer n
  real x(*)
!
  dasum = sum ( abs ( x(1:1+(n-1)*incx:incx) ) )
  return
end
subroutine daxpy ( n, da, dx, incx, dy, incy )
!
!*******************************************************************************
!
!! DAXPY computes constant times a vector plus a vector.
!
!
!     uses unrolled loops for increments equal to one.
!     jack dongarra, linpack, 3/11/78.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
  double precision dx(*),dy(*),da
  integer i,incx,incy,ix,iy,m,n
!
  if ( n <= 0)return
  if (da  ==  0.0D+00 ) return
  if ( incx == 1.and.incy == 1)go to 20
!
!        code for unequal increments or equal increments
!          not equal to 1
!
  ix = 1
  iy = 1
  if ( incx < 0)ix = (-n+1)*incx + 1
  if ( incy < 0)iy = (-n+1)*incy + 1
  do 10 i = 1,n
    dy(iy) = dy(iy) + da*dx(ix)
    ix = ix + incx
    iy = iy + incy
   10 continue
  return
!
!        code for both increments equal to 1
!
!
!        clean-up loop
!
   20 m = mod(n,4)
  if (  m  ==  0 ) go to 40
  do 30 i = 1,m
    dy(i) = dy(i) + da*dx(i)
   30 continue
  if (  n  <  4 ) return
   40 continue
  do i = m+1, n, 4
    dy(i) = dy(i) + da*dx(i)
    dy(i + 1) = dy(i + 1) + da*dx(i + 1)
    dy(i + 2) = dy(i + 2) + da*dx(i + 2)
    dy(i + 3) = dy(i + 3) + da*dx(i + 3)
  end do
  return
end
subroutine dcopy ( n, dx, incx, dy, incy )
!
!*******************************************************************************
!
!! DCOPY copies a vector, x, to a vector, y.
!
!
!     uses unrolled loops for increments equal to one.
!     jack dongarra, linpack, 3/11/78.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
  double precision dx(*),dy(*)
  integer i,incx,incy,ix,iy,m,n
!
  if ( n <= 0)return
  if ( incx == 1.and.incy == 1)go to 20
!
!        code for unequal increments or equal increments
!          not equal to 1
!
  ix = 1
  iy = 1
  if ( incx < 0)ix = (-n+1)*incx + 1
  if ( incy < 0)iy = (-n+1)*incy + 1
  do i = 1,n
    dy(iy) = dx(ix)
    ix = ix + incx
    iy = iy + incy
  end do
  return
!
!        code for both increments equal to 1
!
!
!        clean-up loop
!
   20 m = mod(n,7)
  if (  m  ==  0 ) go to 40
  dy(1:m) = dx(1:m)
  if (  n  <  7 ) return
   40 continue
  do i = m+1, n ,7
    dy(i) = dx(i)
    dy(i + 1) = dx(i + 1)
    dy(i + 2) = dx(i + 2)
    dy(i + 3) = dx(i + 3)
    dy(i + 4) = dx(i + 4)
    dy(i + 5) = dx(i + 5)
    dy(i + 6) = dx(i + 6)
  end do
  return
end
function ddot ( n, dx, incx, dy, incy )
!
!*******************************************************************************
!
!! DDOT forms the dot product of two vectors.
!
!
!     uses unrolled loops for increments equal to one.
!     jack dongarra, linpack, 3/11/78.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
  double precision ddot
  double precision dx(*),dy(*),dtemp
  integer i,incx,incy,ix,iy,m,n
!
  ddot = 0.0D+00
  dtemp = 0.0D+00
  if ( n <= 0)return
  if ( incx == 1.and.incy == 1)go to 20
!
!        code for unequal increments or equal increments
!          not equal to 1
!
  ix = 1
  iy = 1
  if ( incx < 0)ix = (-n+1)*incx + 1
  if ( incy < 0)iy = (-n+1)*incy + 1
  do i = 1,n
    dtemp = dtemp + dx(ix)*dy(iy)
    ix = ix + incx
    iy = iy + incy
  end do
  ddot = dtemp
  return
!
!        code for both increments equal to 1
!
!
!        clean-up loop
!
   20 m = mod(n,5)
  if (  m  ==  0 ) go to 40
  do i = 1,m
    dtemp = dtemp + dx(i)*dy(i)
  end do
  if (  n  <  5 ) go to 60
   40 continue
  do i = m+1, n, 5
    dtemp = dtemp + dx(i)*dy(i) + dx(i + 1)*dy(i + 1) + &
        dx(i + 2)*dy(i + 2) + dx(i + 3)*dy(i + 3) + dx(i + 4)*dy(i + 4)
  end do
   60 ddot = dtemp
  return
end
function dmach ( job )
!
!*******************************************************************************
!
!! DMACH computes machine parameters of floating point arithmetic 
!  for use in testing only.  
!
!  not required by linpack proper.
!
!     if trouble with automatic computation of these quantities,
!     they can be set by direct assignment statements.
!     assume the computer has
!
!        b = base of arithmetic
!        t = number of base  b  digits
!        l = smallest possible exponent
!        u = largest possible exponent
!
!     then
!
!        eps = b**(1-t)
!        tiny = 100.0D+00 *b**(-l+t)
!        huge = 0.01D+00 *b**(u-t)
!
!     dmach same as smach except t, l, u apply to
!     double precision.
!
!     cmach same as smach except if complex division
!     is done by
!
!        1/(x+i*y) = (x-i*y)/(x**2+y**2)
!
!     then
!
!        tiny = sqrt(tiny)
!        huge = sqrt(huge)
!
!
!     job is 1, 2 or 3 for epsilon, tiny and huge, respectively.
!
  double precision dmach
  integer job
  double precision eps,tiny,huge,s
!
  eps = epsilon ( eps )
  s = 1.0D+00
   20 tiny = s
  s = s/16.0D+00
  if (s*1.0 .ne. 0.0D+00) go to 20
  tiny = (tiny/eps)*100.0D+00
  huge = 1.0D+00/tiny
!
  if (job  ==  1) dmach = eps
  if (job  ==  2) dmach = tiny
  if (job  ==  3) dmach = huge
  return
end
function dnrm2 ( n, x, incx )
!
!*******************************************************************************
!
!! DNRM2 returns the euclidean norm of a vector.
!
!
!     dnrm2 := sqrt( x'*x )
!
!
!
!  -- this version written on 25-october-1982.
!     modified on 14-october-1993 to inline the call to dlassq.
!     sven hammarling, nag ltd.
!
!
!     .. parameters ..
  double precision dnrm2
  integer                           incx, n
  double precision                  x( * )
!
  double precision      one         , zero
  parameter           ( one = 1.0d+0, zero = 0.0d+0 )
!
  integer               ix
  double precision      absxi, norm, scale, ssq
  intrinsic             abs, sqrt
!
  if (  n < 1 .or. incx.lt.1 )then
     norm  = zero
  else if (  n == 1 )then
     norm  = abs( x( 1 ) )
  else
     scale = zero
     ssq   = one
!
!        the following loop is equivalent to this call to the lapack
!        auxiliary routine:
!        call dlassq( n, x, incx, scale, ssq )
!
     do ix = 1, 1 + ( n - 1 )*incx, incx
        if (  x( ix ).ne.zero )then
           absxi = abs( x( ix ) )
           if (  scale < absxi )then
              ssq   = one   + ssq*( scale/absxi )**2
              scale = absxi
           else
              ssq   = ssq   +     ( absxi/scale )**2
           end if
        end if
     end do
     norm  = scale * sqrt( ssq )
  end if
  dnrm2 = norm
  return
end
subroutine drot ( n, dx, incx, dy, incy, c, s )
!
!*******************************************************************************
!
!! DROT applies a plane rotation.
!
!     jack dongarra, linpack, 3/11/78.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
  double precision dx(*),dy(*),dtemp,c,s
  integer i,incx,incy,ix,iy,n
!
  if ( n <= 0)return
  if ( incx == 1.and.incy == 1)go to 20
!
!       code for unequal increments or equal increments not equal
!         to 1
!
  ix = 1
  iy = 1
  if ( incx < 0)ix = (-n+1)*incx + 1
  if ( incy < 0)iy = (-n+1)*incy + 1
  do i = 1,n
    dtemp = c*dx(ix) + s*dy(iy)
    dy(iy) = c*dy(iy) - s*dx(ix)
    dx(ix) = dtemp
    ix = ix + incx
    iy = iy + incy
  end do
  return
!
!       code for both increments equal to 1
!
   20 continue
  do i = 1,n
    dtemp = c*dx(i) + s*dy(i)
    dy(i) = c*dy(i) - s*dx(i)
    dx(i) = dtemp
  end do
  return
end
subroutine drotg ( da, db, c, s )
!
!*******************************************************************************
!
!! DROTG constructs a Givens plane rotation.
!
!     jack dongarra, linpack, 3/11/78.
!
  double precision da,db,c,s,roe,scale,r,z
!
  roe = db
  if (  dabs(da) > dabs(db) ) roe = da
  scale = dabs(da) + dabs(db)
  if (  scale == 0.0d0 ) then
    c = 1.0d0
    s = 0.0d0
    r = 0.0d0
    z = 0.0d0
  else
    r = scale*dsqrt((da/scale)**2 + (db/scale)**2)
    r = dsign(1.0d0,roe)*r
    c = da/r
    s = db/r
    z = 1.0d0
    if (  dabs(da) > dabs(db) ) z = s
    if (  dabs(db) >= dabs(da) .and. c /= 0.0d0 ) z = 1.0d0/c
  end if
  da = r
  db = z
  return
end
subroutine drotm ( n, dx, incx, dy, incy, dparam )
!
!*******************************************************************************
!
!! DROTM applies the modified givens transformation, h, to the 2 by n matrix
!
!     (dx**t) , where **t indicates transpose. the elements of dx are in
!     (dy**t)
!
!     dx(lx+i*incx), i = 0 to n-1, where lx = 1 if incx >= 0, else
!     lx = (-incx)*n, and similarly for sy using ly and incy.
!     with dparam(1)=dflag, h has one of the following forms..
!
!     dflag=-1.d0     dflag=0.d0        dflag=1.d0     dflag=-2.d0
!
!       (dh11  dh12)    (1.d0  dh12)    (dh11  1.d0)    (1.d0  0.d0)
!     h=(          )    (          )    (          )    (          )
!       (dh21  dh22),   (dh21  1.d0),   (-1.d0 dh22),   (0.d0  1.d0).
!     see drotmg for a description of data storage in dparam.
!
  double precision dflag
  double precision dh11
  double precision dh12
  double precision dh21
  double precision dh22
  double precision dparam(5)
  double precision dx(*)
  double precision dy(*)
  integer i
  integer incx
  integer incy
  integer kx
  integer ky
  integer n
  integer nsteps
  double precision two
  double precision w
  double precision z
  double precision zero
!
!
  data zero,two/0.d0,2.d0/
!
  dflag=dparam(1)
  if ( n  <=  0 .or.(dflag+two == zero)) go to 140
      if ( .not.(incx == incy.and. incx > 0)) go to 70
!
           nsteps=n*incx
           if ( dflag) 50,10,30
   10          continue
           dh12=dparam(4)
           dh21=dparam(3)
                do 20 i=1,nsteps,incx
                w=dx(i)
                z=dy(i)
                dx(i)=w+z*dh12
                dy(i)=w*dh21+z
   20               continue
           go to 140
   30          continue
           dh11=dparam(2)
           dh22=dparam(5)
                do 40 i=1,nsteps,incx
                w=dx(i)
                z=dy(i)
                dx(i)=w*dh11+z
                dy(i)=-w+dh22*z
   40               continue
           go to 140
   50          continue
           dh11=dparam(2)
           dh12=dparam(4)
           dh21=dparam(3)
           dh22=dparam(5)
                do 60 i=1,nsteps,incx
                w=dx(i)
                z=dy(i)
                dx(i)=w*dh11+z*dh12
                dy(i)=w*dh21+z*dh22
   60               continue
           go to 140
   70     continue
      kx=1
      ky=1
      if ( incx  <  0) kx=1+(1-n)*incx
      if ( incy  <  0) ky=1+(1-n)*incy
      if ( dflag)120,80,100
   80     continue
      dh12=dparam(4)
      dh21=dparam(3)
           do 90 i=1,n
           w=dx(kx)
           z=dy(ky)
           dx(kx)=w+z*dh12
           dy(ky)=w*dh21+z
           kx=kx+incx
           ky=ky+incy
   90          continue
      go to 140
  100     continue
      dh11=dparam(2)
      dh22=dparam(5)
      do i=1,n
           w=dx(kx)
           z=dy(ky)
           dx(kx)=w*dh11+z
           dy(ky)=-w+dh22*z
           kx=kx+incx
           ky=ky+incy
      end do
      go to 140
  120     continue
      dh11=dparam(2)
      dh12=dparam(4)
      dh21=dparam(3)
      dh22=dparam(5)
      do i=1,n
           w=dx(kx)
           z=dy(ky)
           dx(kx)=w*dh11+z*dh12
           dy(ky)=w*dh21+z*dh22
           kx=kx+incx
           ky=ky+incy
      end do
  140     continue
      return
end
!!!!!!!!!!subroutine drotmg ( dd1, dd2, dx1, dy1, dparam )
!!!!!!!!!!!
!!!!!!!!!!!*******************************************************************************
!!!!!!!!!!!
!!!!!!!!!!!! DROTMG constructs the modified givens transformation matrix h which zeros
!!!!!!!!!!!     the second component of the 2-vector  (dsqrt(dd1)*dx1,dsqrt(dd2)*
!!!!!!!!!!!     dy2)**t.
!!!!!!!!!!!     with dparam(1)=dflag, h has one of the following forms..
!!!!!!!!!!!
!!!!!!!!!!!     dflag=-1.d0     dflag=0.d0        dflag=1.d0     dflag=-2.d0
!!!!!!!!!!!
!!!!!!!!!!!       (dh11  dh12)    (1.d0  dh12)    (dh11  1.d0)    (1.d0  0.d0)
!!!!!!!!!!!     h=(          )    (          )    (          )    (          )
!!!!!!!!!!!       (dh21  dh22),   (dh21  1.d0),   (-1.d0 dh22),   (0.d0  1.d0).
!!!!!!!!!!!     locations 2-4 of dparam contain dh11, dh21, dh12, and dh22
!!!!!!!!!!!     respectively. (values of 1.d0, -1.d0, or 0.d0 implied by the
!!!!!!!!!!!     value of dparam(1) are not stored in dparam.)
!!!!!!!!!!!
!!!!!!!!!!!     the values of gamsq and rgamsq set in the data statement may be
!!!!!!!!!!!     inexact.  this is ok as they are only used for testing the size
!!!!!!!!!!!     of dd1 and dd2.  all actual scaling of data is done using gam.
!!!!!!!!!!!
!!!!!!!!!!  double precision dd1
!!!!!!!!!!  double precision dd2
!!!!!!!!!!  double precision dflag
!!!!!!!!!!  double precision dh11
!!!!!!!!!!  double precision dh12
!!!!!!!!!!  double precision dh21
!!!!!!!!!!  double precision dh22
!!!!!!!!!!  double precision dp1
!!!!!!!!!!  double precision dp2
!!!!!!!!!!  double precision dparam(5)
!!!!!!!!!!  double precision dq1
!!!!!!!!!!  double precision dq2
!!!!!!!!!!  double precision dtemp
!!!!!!!!!!  double precision du
!!!!!!!!!!  double precision dx1
!!!!!!!!!!  double precision dy1
!!!!!!!!!!  double precision gam
!!!!!!!!!!  double precision gamsq
!!!!!!!!!!  integer igo
!!!!!!!!!!  double precision one
!!!!!!!!!!  double precision rgamsq
!!!!!!!!!!  double precision two
!!!!!!!!!!  double precision zero
!!!!!!!!!!!
!!!!!!!!!!  data zero,one,two /0.d0,1.d0,2.d0/
!!!!!!!!!!  data gam,gamsq,rgamsq/4096.d0,16777216.d0,5.9604645d-8/
!!!!!!!!!!!
!!!!!!!!!!  if ( .not. dd1  <  zero) go to 10
!!!!!!!!!!!       go zero-h-d-and-dx1..
!!!!!!!!!!      go to 60
!!!!!!!!!!   10 continue
!!!!!!!!!!!     case-dd1-nonnegative
!!!!!!!!!!  dp2=dd2*dy1
!!!!!!!!!!  if ( .not. dp2  ==  zero) go to 20
!!!!!!!!!!      dflag=-two
!!!!!!!!!!      go to 260
!!!!!!!!!!!     regular-case..
!!!!!!!!!!   20 continue
!!!!!!!!!!  dp1=dd1*dx1
!!!!!!!!!!  dq2=dp2*dy1
!!!!!!!!!!  dq1=dp1*dx1
!!!!!!!!!!  if ( .not. dabs(dq1) > dabs(dq2)) go to 40
!!!!!!!!!!      dh21=-dy1/dx1
!!!!!!!!!!      dh12=dp2/dp1
!!!!!!!!!!      du=one-dh12*dh21
!!!!!!!!!!      if ( .not. du  <=  zero) go to 30
!!!!!!!!!!!         go zero-h-d-and-dx1..
!!!!!!!!!!           go to 60
!!!!!!!!!!   30     continue
!!!!!!!!!!           dflag=zero
!!!!!!!!!!           dd1=dd1/du
!!!!!!!!!!           dd2=dd2/du
!!!!!!!!!!           dx1=dx1*du
!!!!!!!!!!!         go scale-check..
!!!!!!!!!!           go to 100
!!!!!!!!!!   40 continue
!!!!!!!!!!      if ( .not. dq2  <  zero) go to 50
!!!!!!!!!!!         go zero-h-d-and-dx1..
!!!!!!!!!!           go to 60
!!!!!!!!!!   50     continue
!!!!!!!!!!           dflag=one
!!!!!!!!!!           dh11=dp1/dp2
!!!!!!!!!!           dh22=dx1/dy1
!!!!!!!!!!           du=one+dh11*dh22
!!!!!!!!!!           dtemp=dd2/du
!!!!!!!!!!           dd2=dd1/du
!!!!!!!!!!           dd1=dtemp
!!!!!!!!!!           dx1=dy1*du
!!!!!!!!!!!         go scale-check
!!!!!!!!!!           go to 100
!!!!!!!!!!!     procedure..zero-h-d-and-dx1..
!!!!!!!!!!   60 continue
!!!!!!!!!!      dflag=-one
!!!!!!!!!!      dh11=zero
!!!!!!!!!!      dh12=zero
!!!!!!!!!!      dh21=zero
!!!!!!!!!!      dh22=zero
!!!!!!!!!!      dd1=zero
!!!!!!!!!!      dd2=zero
!!!!!!!!!!      dx1=zero
!!!!!!!!!!!         return..
!!!!!!!!!!      go to 220
!!!!!!!!!!!     procedure..fix-h..
!!!!!!!!!!   70 continue
!!!!!!!!!!  if ( .not. dflag >= zero) go to 90
!!!!!!!!!!      if ( .not. dflag  ==  zero) go to 80
!!!!!!!!!!      dh11=one
!!!!!!!!!!      dh22=one
!!!!!!!!!!      dflag=-one
!!!!!!!!!!      go to 90
!!!!!!!!!!   80     continue
!!!!!!!!!!      dh21=-one
!!!!!!!!!!      dh12=one
!!!!!!!!!!      dflag=-one
!!!!!!!!!!   90 continue
!!!!!!!!!!  go to igo,(120,150,180,210)
!!!!!!!!!!!     procedure..scale-check
!!!!!!!!!!  100 continue
!!!!!!!!!!  110     continue
!!!!!!!!!!      if ( .not. dd1  <=  rgamsq) go to 130
!!!!!!!!!!           if ( dd1  ==  zero) go to 160
!!!!!!!!!!           assign 120 to igo
!!!!!!!!!!!              fix-h..
!!!!!!!!!!           go to 70
!!!!!!!!!!  120          continue
!!!!!!!!!!           dd1=dd1*gam**2
!!!!!!!!!!           dx1=dx1/gam
!!!!!!!!!!           dh11=dh11/gam
!!!!!!!!!!           dh12=dh12/gam
!!!!!!!!!!      go to 110
!!!!!!!!!!  130 continue
!!!!!!!!!!  140     continue
!!!!!!!!!!      if ( .not. dd1 >= gamsq) go to 160
!!!!!!!!!!           assign 150 to igo
!!!!!!!!!!!              fix-h..
!!!!!!!!!!           go to 70
!!!!!!!!!!  150          continue
!!!!!!!!!!           dd1=dd1/gam**2
!!!!!!!!!!           dx1=dx1*gam
!!!!!!!!!!           dh11=dh11*gam
!!!!!!!!!!           dh12=dh12*gam
!!!!!!!!!!      go to 140
!!!!!!!!!!  160 continue
!!!!!!!!!!  170     continue
!!!!!!!!!!      if ( .not. dabs(dd2)  <=  rgamsq) go to 190
!!!!!!!!!!           if ( dd2  ==  zero) go to 220
!!!!!!!!!!           assign 180 to igo
!!!!!!!!!!!              fix-h..
!!!!!!!!!!           go to 70
!!!!!!!!!!  180          continue
!!!!!!!!!!           dd2=dd2*gam**2
!!!!!!!!!!           dh21=dh21/gam
!!!!!!!!!!           dh22=dh22/gam
!!!!!!!!!!      go to 170
!!!!!!!!!!  190 continue
!!!!!!!!!!  200     continue
!!!!!!!!!!      if ( .not. dabs(dd2) >= gamsq) go to 220
!!!!!!!!!!           assign 210 to igo
!!!!!!!!!!!              fix-h..
!!!!!!!!!!           go to 70
!!!!!!!!!!  210          continue
!!!!!!!!!!           dd2=dd2/gam**2
!!!!!!!!!!           dh21=dh21*gam
!!!!!!!!!!           dh22=dh22*gam
!!!!!!!!!!      go to 200
!!!!!!!!!!  220 continue
!!!!!!!!!!      if ( dflag)250,230,240
!!!!!!!!!!  230     continue
!!!!!!!!!!           dparam(3)=dh21
!!!!!!!!!!           dparam(4)=dh12
!!!!!!!!!!           go to 260
!!!!!!!!!!  240     continue
!!!!!!!!!!           dparam(2)=dh11
!!!!!!!!!!           dparam(5)=dh22
!!!!!!!!!!           go to 260
!!!!!!!!!!  250     continue
!!!!!!!!!!           dparam(2)=dh11
!!!!!!!!!!           dparam(3)=dh21
!!!!!!!!!!           dparam(4)=dh12
!!!!!!!!!!           dparam(5)=dh22
!!!!!!!!!!  260 continue
!!!!!!!!!!      dparam(1)=dflag
!!!!!!!!!!  return
!!!!!!!!!!end
subroutine dscal ( n, da, dx, incx )
!
!*******************************************************************************
!
!! DSCAL scales a vector by a constant.
!     uses unrolled loops for increment equal to one.
!     jack dongarra, linpack, 3/11/78.
!     modified 3/93 to return if incx  <=  0.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
  double precision da,dx(*)
  integer i,incx,m,n,nincx
!
  if (  n <= 0 .or. incx <= 0 )return
  if ( incx == 1)go to 20
!
!        code for increment not equal to 1
!
  nincx = n*incx
  do i = 1,nincx,incx
    dx(i) = da*dx(i)
  end do
  return
!
!        code for increment equal to 1
!
!
!        clean-up loop
!
   20 continue
  m = mod(n,5)
  if (  m  ==  0 ) go to 40
  dx(1:m) = da*dx(1:m)
  if (  n  <  5 ) return
40 continue
  do i = m+1,n,5
    dx(i) = da*dx(i)
    dx(i + 1) = da*dx(i + 1)
    dx(i + 2) = da*dx(i + 2)
    dx(i + 3) = da*dx(i + 3)
    dx(i + 4) = da*dx(i + 4)
  end do
  return
end
subroutine dswap ( n, x, incx, y, incy )
!
!*******************************************************************************
!
!! DSWAP interchanges two vectors.
!
!
!  Modified:
!
!    08 April 1999
!
!  Reference:
!
!    Lawson, Hanson, Kincaid, Krogh,
!    Basic Linear Algebra Subprograms for Fortran Usage,
!    Algorithm 539,
!    ACM Transactions on Mathematical Software,
!    Volume 5, Number 3, September 1979, pages 308-323.
!
!  Parameters:
!
!    Input, integer N, the number of entries in the vectors.
!
!    Input/output, double precision X(*), one of the vectors to swap.
!
!    Input, integer INCX, the increment between successive entries of X.
!
!    Input/output, double precision Y(*), one of the vectors to swap.
!
!    Input, integer INCY, the increment between successive elements of Y.
!
  integer i
  integer incx
  integer incy
  integer ix
  integer iy
  integer m
  integer n
  double precision temp
  double precision x(*)
  double precision y(*)
!
  if ( n <= 0 ) then
  else if ( incx == 1 .and. incy == 1 ) then
    m = mod ( n, 3 )
    do i = 1, m
      temp = x(i)
      x(i) = y(i)
      y(i) = temp
    end do
    do i = m+1, n, 3
      temp = x(i)
      x(i) = y(i)
      y(i) = temp
      temp = x(i+1)
      x(i+1) = y(i+1)
      y(i+1) = temp
      temp = x(i+2)
      x(i+2) = y(i+2)
      y(i+2) = temp
    end do
  else
    if ( incx >= 0 ) then
      ix = 1
    else
      ix = ( - n + 1 ) * incx + 1
    end if
    if ( incy >= 0 ) then
      iy = 1
    else
      iy = ( - n + 1 ) * incy + 1
    end if
    do i = 1, n
      temp = x(ix)
      x(ix) = y(iy)
      y(iy) = temp
      ix = ix + incx
      iy = iy + incy
    end do
  end if
  return
end
function idamax ( n, dx, incx )
!
!*******************************************************************************
!
!! IDAMAX finds the index of element having max. absolute value.
!
!     jack dongarra, linpack, 3/11/78.
!     modified 3/93 to return if incx  <=  0.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
  integer idamax
  double precision dmax
  double precision dx(*)
  integer i,incx,ix,n
!
  idamax = 0
  if (  n < 1 .or. incx <= 0 ) return
  idamax = 1
  if ( n == 1)return
  if ( incx == 1)go to 20
!
!        code for increment not equal to 1
!
  ix = 1
  dmax = dabs(dx(1))
  ix = ix + incx
  do i = 2,n
     if ( dabs(dx(ix)) <= dmax) go to 5
     idamax = i
     dmax = dabs(dx(ix))
    5    ix = ix + incx
  end do
  return
!
!        code for increment equal to 1
!
   20 continue
  dmax = abs ( dx(1) )
  do i = 2,n
    if ( abs ( dx(i) ) > dmax ) then
      idamax = i
      dmax = abs ( dx(i) )
    end if
  end do
  return
end
function lsame ( CA, CB )
!
!*******************************************************************************
!
!! LSAME returns .TRUE. if CA is the same letter as CB regardless of case.
!
!
!  Arguments
!
!  CA      (input) CHARACTER*1
!  CB      (input) CHARACTER*1
!          CA and CB specify the single characters to be compared.
!
  logical lsame
  CHARACTER CA, CB
  INTEGER INTA, INTB, ZCODE
!
!  Test if the characters are equal
!
  LSAME = ( CA == CB )
  IF ( LSAME ) then
    RETURN
  end if
!
!  Now test for equivalence if both characters are alphabetic.
!
  ZCODE = ICHAR ( 'Z' )
!
!  Use 'Z' rather than 'A' so that ASCII can be detected on Prime
!  machines, on which ICHAR returns a value with bit 8 set.
!  ICHAR('A') on Prime machines returns 193 which is the same as
!  ICHAR('A') on an EBCDIC machine.
!
  INTA = ICHAR ( CA )
  INTB = ICHAR ( CB )
  IF ( ZCODE == 90 .OR. ZCODE == 122 ) THEN
!
!  ASCII is assumed - ZCODE is the ASCII code of either lower or
!  upper case 'Z'.
!
    IF ( INTA >= 97 .AND. INTA <= 122 ) INTA = INTA - 32
    IF ( INTB >= 97 .AND. INTB <= 122 ) INTB = INTB - 32
  ELSE IF ( ZCODE == 233 .OR. ZCODE == 169 ) THEN
!
!  EBCDIC is assumed - ZCODE is the EBCDIC code of either lower or
!  upper case 'Z'.
!
    IF ( INTA >= 129 .AND. INTA <= 137 .OR. &
         INTA >= 145 .AND. INTA <= 153 .OR. &
         INTA >= 162 .AND. INTA <= 169 ) INTA = INTA + 64
    IF ( INTB >= 129 .AND. INTB <= 137 .OR. &
         INTB >= 145 .AND. INTB <= 153 .OR. &
         INTB >= 162 .AND. INTB <= 169 ) INTB = INTB + 64
  ELSE IF ( ZCODE == 218 .OR. ZCODE == 250 ) THEN
!
!  ASCII is assumed, on Prime machines - ZCODE is the ASCII code
!  plus 128 of either lower or upper case 'Z'.
!
    IF ( INTA >= 225 .AND. INTA <= 250 ) INTA = INTA - 32
    IF ( INTB >= 225 .AND. INTB <= 250 ) INTB = INTB - 32
  END IF
  LSAME = ( INTA == INTB )
  return
end
SUBROUTINE XERBLA ( SRNAME, INFO )
!
!*******************************************************************************
!
!! XERBLA is an error handler for the LAPACK routines.
!
!
!  Arguments
!
!  SRNAME  (input) CHARACTER*6
!          The name of the routine which called XERBLA.
!
!  INFO    (input) INTEGER
!          The position of the invalid parameter in the parameter list
!          of the calling routine.
!
  CHARACTER ( len = 6 ) SRNAME
  INTEGER INFO
!
  WRITE ( *, FMT = 9999 )SRNAME, INFO
  STOP
 9999 FORMAT( ' ** On entry to ', A6, ' parameter number ', I2, ' had ', &
  'an illegal value' )
END

module module_blas_etc

  use precision_kinds, only: dp
  implicit none

contains

  pure subroutine daxpy(n,da,dx,incx,dy,incy)
    implicit none
    integer, intent(in) :: n, incx, incy
    real(dp), intent(in) :: da, dx(n)
    real(dp), intent(inout) :: dy(n)
    if (n<=0) then
      return
    else
      dy = dy + da*dx
    end if
  end subroutine daxpy

  pure subroutine dcopy(n,dx,incx,dy,incy)
    implicit none
    integer, intent(in) :: n, incx, incy
    real(dp), intent(in) :: dx(n)
    real(dp), intent(inout) :: dy(n)
    if (n<=0) then
      return
    else
      dy = dx
    end if
  end subroutine dcopy

  pure function ddot(n,dx,incx,dy,incy)
    implicit none
    real(dp) :: ddot
    integer, intent(in) :: n, incx, incy
    real(dp), intent(in) :: dx(n),dy(n)
    real(dp) :: dtemp
    integer i,ix,iy,m,mp1
    ddot = 0.0d0
    dtemp = 0.0d0
    if(n.le.0)return
    if(incx.eq.1.and.incy.eq.1) go to 20
    ix = 1
    iy = 1
    if(incx.lt.0)ix = (-n+1)*incx + 1
    if(incy.lt.0)iy = (-n+1)*incy + 1
    do i = 1,n
      dtemp = dtemp + dx(ix)*dy(iy)
      ix = ix + incx
      iy = iy + incy
    end do
    ddot = dtemp
    return
    20 m = mod(n,5)
    if( m .eq. 0 ) go to 40
    do i = 1,m
      dtemp = dtemp + dx(i)*dy(i)
    end do
    if( n .lt. 5 ) go to 60
    40 mp1 = m + 1
    do i = mp1,n,5
      dtemp = dtemp + dx(i)*dy(i) + dx(i + 1)*dy(i + 1) +  dx(i + 2)*dy(i + 2) + dx(i + 3)*dy(i + 3) + dx(i + 4)*dy(i + 4)
    end do
    60 ddot = dtemp
    return
  end function ddot

  subroutine  dscal(n,da,dx,incx)
    implicit none
    real(dp) da,dx(*)
    integer i,incx,m,mp1,n,nincx
    if( n.le.0 .or. incx.le.0 )return
    if(incx.eq.1)go to 20
    nincx = n*incx
    do 10 i = 1,nincx,incx
      dx(i) = da*dx(i)
    10 continue
      return
    20 m = mod(n,5)
    if( m .eq. 0 ) go to 40
    do 30 i = 1,m
      dx(i) = da*dx(i)
    30 continue
    if( n .lt. 5 ) return
    40 mp1 = m + 1
    do 50 i = mp1,n,5
      dx(i) = da*dx(i)
      dx(i + 1) = da*dx(i + 1)
      dx(i + 2) = da*dx(i + 2)
      dx(i + 3) = da*dx(i + 3)
      dx(i + 4) = da*dx(i + 4)
    50 continue
    return
  end subroutine  dscal

  subroutine dpofa(a,lda,n,info)
    implicit none
    integer lda,n,info
    real(dp) a(lda,*)
    real(dp) t!,ddot
    real(dp) s
    integer j,jm1,k
    do 30 j = 1, n
      info = j
      s = 0.0d0
      jm1 = j - 1
      if (jm1 .lt. 1) go to 20
      do k = 1, jm1
        t = a(k,j) - ddot(k-1,a(1,k),1,a(1,j),1)
        t = t/a(k,k)
                 a(k,j) = t
                 s = s + t*t
              end do
     20       continue
              s = a(j,j) - s
              if (s .le. 0.0d0) return
              a(j,j) = sqrt(s)
     30    continue
           info = 0
  end subroutine dpofa

  pure subroutine dtrsl(t,ldt,n,b,job,info)
    !
    ! linpack's dtrsl is the time-consuming routine of lbfgs
    !
    implicit none
    integer, intent(in) :: ldt,n,job
    integer, intent(out) :: info
    real(dp), intent(in) :: t(ldt,n)
    real(dp), intent(inout) :: b(n)
    real(dp) temp!,ddot
    integer case,j,jj
    do info = 1, n
      if (t(info,info) .eq. 0.0d0) return
    end do
    info = 0
    case = 1
    if (mod(job,10) .ne. 0) case = 2
    if (mod(job,100)/10 .ne. 0) case = case + 2
    go to (20,50,80,110), case
    20  continue
    b(1) = b(1)/t(1,1)
    if (n .lt. 2) go to 40
    do j = 2, n
      temp = -b(j-1)
      call daxpy(n-j+1,temp,t(j,j-1),1,b(j),1)
      b(j) = b(j)/t(j,j)
    end do
    40       continue
    go to 140
    50    continue
    b(n) = b(n)/t(n,n)
    if (n .lt. 2) go to 70
    do jj = 2, n
      j = n - jj + 1
      temp = -b(j+1)
      call daxpy(j,temp,t(1,j+1),1,b(1),1)
      b(j) = b(j)/t(j,j)
    end do
    70       continue
    go to 140
    80    continue
    b(n) = b(n)/t(n,n)
    if (n .lt. 2) go to 100
    do jj = 2, n
      j = n - jj + 1
      b(j) = b(j) - ddot(jj-1,t(j+1,j),1,b(j+1),1)
      b(j) = b(j)/t(j,j)
    end do
    100       continue
    go to 140
    110    continue
    b(1) = b(1)/t(1,1)
    if (n .lt. 2) go to 130
    do j = 2, n
      b(j) = b(j) - ddot(j-1,t(1,j),1,b(1),1)
      b(j) = b(j)/t(j,j)
    end do
    130       continue
    140    continue
    150 continue
  end subroutine dtrsl

end module module_blas_etc

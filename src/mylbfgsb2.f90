module module_blas_etc
  use iso_c_binding, only: c_double, c_float
  implicit none
contains

  pure subroutine daxpy(n,da,dx,incx,dy,incy)
    implicit none
    integer, intent(in) :: n, incx, incy
    real(c_double), intent(in) :: da, dx(n)
    real(c_double), intent(inout) :: dy(n)
    if (n<=0) then
      return
    else
      dy = dy + da*dx
    end if
  end subroutine daxpy

  pure subroutine dcopy(n,dx,incx,dy,incy)
    implicit none
    integer, intent(in) :: n, incx, incy
    real(c_double), intent(in) :: dx(n)
    real(c_double), intent(inout) :: dy(n)
    if (n<=0) then
      return
    else
      dy = dx
    end if
  end subroutine dcopy

  pure function ddot(n,dx,incx,dy,incy)
    implicit none
    real(c_double) :: ddot
    integer, intent(in) :: n, incx, incy
    real(c_double), intent(in) :: dx(n),dy(n)
    real(c_double) :: dtemp
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
        real(c_double) da,dx(*)
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
        real(c_double) a(lda,*)

        real(c_double) t!,ddot
        real(c_double) s
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
              if (s .le. 0.0d0) go to 40
              a(j,j) = sqrt(s)
     30    continue
           info = 0
     40 continue
        return
        end subroutine dpofa


        subroutine dtrsl(t,ldt,n,b,job,info)
        implicit none
        integer, intent(in) :: ldt,n,job
        integer, intent(out) :: info
        real(c_double), intent(in) :: t(ldt,n)
        real(c_double), intent(out) :: b(n)
        real(c_double) temp!,ddot
        integer case,j,jj
           do 10 info = 1, n
              if (t(info,info) .eq. 0.0d0) go to 150
     10    continue
           info = 0
           case = 1
           if (mod(job,10) .ne. 0) case = 2
           if (mod(job,100)/10 .ne. 0) case = case + 2
           go to (20,50,80,110), case
     20    continue
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
              do 60 jj = 2, n
                 j = n - jj + 1
                 temp = -b(j+1)
                 call daxpy(j,temp,t(1,j+1),1,b(1),1)
                 b(j) = b(j)/t(j,j)
     60       continue
     70       continue
           go to 140
     80    continue
              b(n) = b(n)/t(n,n)
              if (n .lt. 2) go to 100
              do 90 jj = 2, n
                 j = n - jj + 1
                 b(j) = b(j) - ddot(jj-1,t(j+1,j),1,b(j+1),1)
                 b(j) = b(j)/t(j,j)
     90       continue
    100       continue
           go to 140
    110    continue
              b(1) = b(1)/t(1,1)
              if (n .lt. 2) go to 130
              do 120 j = 2, n
                 b(j) = b(j) - ddot(j-1,t(1,j),1,b(1),1)
                 b(j) = b(j)/t(j,j)
    120       continue
    130       continue
    140    continue
    150 continue
        return
        end subroutine dtrsl

end module module_blas_etc




module module_lbfgs

  use iso_c_binding, only: c_double, c_float
  use module_blas_etc

contains

    ! subroutine setulb(n, m, x, l, u, nbd, f, g, factr, pgtol, wa, iwa, task, iprint, csave, lsave, isave, dsave) MAX 9 MAI 2016 SINCE WE ALWAYS DO UNCONSTRAINED MINIMIZATION
      subroutine setulb(n, m, x,            f, g, factr, pgtol, wa, iwa, task, iprint, csave, lsave, isave, dsave)
        implicit none
      character(len=60)     task, csave
      logical          lsave(4)
      integer          n, m, iprint, nbd(n), iwa(3*n), isave(44)
      real(c_double) :: f, factr, pgtol, x(n), l(n), u(n), g(n), wa(2*m*n + 5*n + 11*m*m + 8*m), dsave(29)


      integer   lws,lr,lz,lt,ld,lxp,lwa,lwy,lsy,lss,lwt,lwn,lsnd

      nbd=0
      l=0.
      u=0.

      if (task .eq. 'START') then
         isave(1)  = m*n
         isave(2)  = m**2
         isave(3)  = 4*m**2
         isave(4)  = 1                      ! ws      m*n
         isave(5)  = isave(4)  + isave(1)   ! wy      m*n
         isave(6)  = isave(5)  + isave(1)   ! wsy     m**2
         isave(7)  = isave(6)  + isave(2)   ! wss     m**2
         isave(8)  = isave(7)  + isave(2)   ! wt      m**2
         isave(9)  = isave(8)  + isave(2)   ! wn      4*m**2
         isave(10) = isave(9)  + isave(3)   ! wsnd    4*m**2
         isave(11) = isave(10) + isave(3)   ! wz      n
         isave(12) = isave(11) + n          ! wr      n
         isave(13) = isave(12) + n          ! wd      n
         isave(14) = isave(13) + n          ! wt      n
         isave(15) = isave(14) + n          ! wxp     n
         isave(16) = isave(15) + n          ! wa      8*m
      endif
      lws  = isave(4)
      lwy  = isave(5)
      lsy  = isave(6)
      lss  = isave(7)
      lwt  = isave(8)
      lwn  = isave(9)
      lsnd = isave(10)
      lz   = isave(11)
      lr   = isave(12)
      ld   = isave(13)
      lt   = isave(14)
      lxp  = isave(15)
      lwa  = isave(16)

      call mainlb(n,m,x,l,u,nbd,f,g,factr,pgtol, wa(lws),wa(lwy),wa(lsy),wa(lss), wa(lwt),&
       wa(lwn),wa(lsnd),wa(lz),wa(lr),wa(ld),wa(lt),wa(lxp),wa(lwa), iwa(1),iwa(n+1),iwa(2*n+1),task,iprint,&
       csave,lsave,isave(22),dsave)

    end subroutine setulb


      subroutine mainlb(n, m, x, l, u, nbd, f, g, factr, pgtol, ws, wy,&
                      sy, ss, wt, wn, snd, z, r, d, t, xp, wa,&
                      index, iwhere, indx2, task,&
                       iprint, csave, lsave, isave, dsave)
      implicit none
      character(len=60)     task, csave
      logical          lsave(4)
      integer          n, m, iprint, nbd(n), index(n), iwhere(n), indx2(n), isave(23)
      real(c_double) f, factr, pgtol, x(n), l(n), u(n), g(n), z(n), r(n), d(n), t(n), xp(n), wa(8*m), &
                      ws(n, m), wy(n, m), sy(m, m), ss(m, m),  wt(m, m), wn(2*m, 2*m), snd(2*m, 2*m), dsave(29)


      logical          prjctd,cnstnd,boxed,updatd,wrk
      character(len=3)      word
      integer          i,k,nintol,itfile,iback,nskip,head,col,iter,itail,iupdat,nseg,nfgv,info,ifun,iword,nfree,nact,ileave,nenter
      real(c_double) theta,fold,dr,rr,tol,xstep,sbgnrm,ddum,dnorm,dtd,epsmch,cpu1,cpu2,cachyt,sbtime,lnscht,time1,time2,&
                      gd,gdold,stp,stpmx,time!,ddot
      real(c_double) one,zero
      parameter        (one=1.0d0,zero=0.0d0)

      if (task .eq. 'START') then

         epsmch = epsilon(one)

         call cpu_time(time1)

         col    = 0
         head   = 1
         theta  = one
         iupdat = 0
         updatd = .false.
         iback  = 0
         itail  = 0
         iword  = 0
         nact   = 0
         ileave = 0
         nenter = 0
         fold   = zero
         dnorm  = zero
         cpu1   = zero
         gd     = zero
         stpmx  = zero
         sbgnrm = zero
         stp    = zero
         gdold  = zero
         dtd    = zero

         iter   = 0
         nfgv   = 0
         nseg   = 0
         nintol = 0
         nskip  = 0
         nfree  = n
         ifun   = 0
         tol = factr*epsmch

         cachyt = 0
         sbtime = 0
         lnscht = 0

         word = '---'

         info = 0

         itfile = 8
         if (iprint .ge. 1) then
            open (8, file = 'output/iterate.dat', status = 'unknown')
         endif


         call errclb(n,m,factr,l,u,nbd,task,info,k)
         if (task(1:5) .eq. 'ERROR') then
            call prn3lb(n,x,f,task,iprint,info,itfile, iter,nfgv,nintol,nskip,nact,sbgnrm,&
                       zero,nseg,word,iback,stp,xstep,k,&
                       cachyt,sbtime,lnscht)
            return
         endif

         call prn1lb(n,m,l,u,x,iprint,itfile,epsmch)


         call active(n,l,u,nbd,x,iwhere,iprint,prjctd,cnstnd,boxed)


      else

         prjctd = lsave(1)
         cnstnd = lsave(2)
         boxed  = lsave(3)
         updatd = lsave(4)

         nintol = isave(1)
         itfile = isave(3)
         iback  = isave(4)
         nskip  = isave(5)
         head   = isave(6)
         col    = isave(7)
         itail  = isave(8)
         iter   = isave(9)
         iupdat = isave(10)
         nseg   = isave(12)
         nfgv   = isave(13)
         info   = isave(14)
         ifun   = isave(15)
         iword  = isave(16)
         nfree  = isave(17)
         nact   = isave(18)
         ileave = isave(19)
         nenter = isave(20)

         theta  = dsave(1)
         fold   = dsave(2)
         tol    = dsave(3)
         dnorm  = dsave(4)
         epsmch = dsave(5)
         cpu1   = dsave(6)
         cachyt = dsave(7)
         sbtime = dsave(8)
         lnscht = dsave(9)
         time1  = dsave(10)
         gd     = dsave(11)
         stpmx  = dsave(12)
         sbgnrm = dsave(13)
         stp    = dsave(14)
         gdold  = dsave(15)
         dtd    = dsave(16)


         if (task(1:5) .eq. 'FG_LN') goto 666
         if (task(1:5) .eq. 'NEW_X') goto 777
         if (task(1:5) .eq. 'FG_ST') goto 111
         if (task(1:4) .eq. 'STOP') then
            if (task(7:9) .eq. 'CPU') then
               call dcopy(n,t,1,x,1)
               call dcopy(n,r,1,g,1)
               f = fold
            endif
            goto 999
         endif
      endif


      task = 'FG_START'
      goto 1000
 111  continue
      nfgv = 1


      call projgr(n,l,u,nbd,x,g,sbgnrm)

      if (iprint .ge. 1) then
         write (6,1002) iter,f,sbgnrm
         write (itfile,1003) iter,nfgv,sbgnrm,f
      endif
      if (sbgnrm .le. pgtol) then
         task = 'CONVERGENCE: NORM_OF_PROJECTED_GRADIENT_<=_PGTOL'
         goto 999
      endif


 222  continue
      if (iprint .ge. 99) write (6,1001) iter + 1
      iword = -1
      if (.not. cnstnd .and. col .gt. 0) then
         call dcopy(n,x,1,z,1)
         wrk = updatd
         nseg = 0
         goto 333
      endif


      call cpu_time(cpu1)
      call cauchy(n,x,l,u,nbd,g,indx2,iwhere,t,d,z, m,wy,ws,sy,wt,theta,col,head, wa(1),wa(2*m+1),wa(4*m+1),wa(6*m+1),nseg,&
                 iprint, sbgnrm, info, epsmch)
      if (info .ne. 0) then
         if(iprint .ge. 1) write (6, 1005)
         info   = 0
         col    = 0
         head   = 1
         theta  = one
         iupdat = 0
         updatd = .false.
         call cpu_time(cpu2)
         cachyt = cachyt + cpu2 - cpu1
         goto 222
      endif
      call cpu_time(cpu2)
      cachyt = cachyt + cpu2 - cpu1
      nintol = nintol + nseg


      call freev(n,nfree,index,nenter,ileave,indx2,iwhere,wrk,updatd,cnstnd,iprint,iter)
      nact = n - nfree

 333  continue


      if (nfree .eq. 0 .or. col .eq. 0) goto 555

      call cpu_time(cpu1)


      if (wrk) call formk(n,nfree,index,nenter,ileave,indx2,iupdat,updatd,wn,snd,m,ws,wy,sy,theta,col,head,info)
      if (info .ne. 0) then
         if(iprint .ge. 1) write (6, 1006)
         info   = 0
         col    = 0
         head   = 1
         theta  = one
         iupdat = 0
         updatd = .false.
         call cpu_time(cpu2)
         sbtime = sbtime + cpu2 - cpu1
         goto 222
      endif

      call cmprlb(n,m,x,g,ws,wy,sy,wt,z,r,wa,index,theta,col,head,nfree,cnstnd,info)
      if (info .ne. 0) goto 444


      call subsm( n, m, nfree, index, l, u, nbd, z, r, xp, ws, wy, theta, x, g, col, head, iword, wa, wn, iprint, info)
 444  continue
      if (info .ne. 0) then
         if(iprint .ge. 1) write (6, 1005)
         info   = 0
         col    = 0
         head   = 1
         theta  = one
         iupdat = 0
         updatd = .false.
         call cpu_time(cpu2)
         sbtime = sbtime + cpu2 - cpu1
         goto 222
      endif

      call cpu_time(cpu2)
      sbtime = sbtime + cpu2 - cpu1
 555  continue


      do 40 i = 1, n
         d(i) = z(i) - x(i)
  40  continue
      call cpu_time(cpu1)
 666  continue
      call lnsrlb(n,l,u,nbd,x,f,fold,gd,gdold,g,d,r,t,z,stp,dnorm, dtd,xstep,stpmx,iter,ifun,iback,nfgv,info,task,&
                 boxed,cnstnd,csave,isave(22),dsave(17))
      if (info .ne. 0 .or. iback .ge. 20) then
         call dcopy(n,t,1,x,1)
         call dcopy(n,r,1,g,1)
         f = fold
         if (col .eq. 0) then
            if (info .eq. 0) then
               info = -9
               nfgv = nfgv - 1
               ifun = ifun - 1
               iback = iback - 1
            endif
            task = 'ABNORMAL_TERMINATION_IN_LNSRCH'
            iter = iter + 1
            goto 999
         else
            if(iprint .ge. 1) write (6, 1008)
            if (info .eq. 0) nfgv = nfgv - 1
            info   = 0
            col    = 0
            head   = 1
            theta  = one
            iupdat = 0
            updatd = .false.
            task   = 'RESTART_FROM_LNSRCH'
            call cpu_time(cpu2)
            lnscht = lnscht + cpu2 - cpu1
            goto 222
         endif
      else if (task(1:5) .eq. 'FG_LN') then
         goto 1000
      else
         call cpu_time(cpu2)
         lnscht = lnscht + cpu2 - cpu1
         iter = iter + 1


         call projgr(n,l,u,nbd,x,g,sbgnrm)


         call prn2lb(n,x,f,g,iprint,itfile,iter,nfgv,nact,sbgnrm,nseg,word,iword,iback,stp,xstep)
         goto 1000
      endif
 777  continue


      if (sbgnrm .le. pgtol) then
         task = 'CONVERGENCE: NORM_OF_PROJECTED_GRADIENT_<=_PGTOL'
         goto 999
      endif

      ddum = max(abs(fold), abs(f), one)
      if ((fold - f) .le. tol*ddum) then
         task = 'CONVERGENCE: REL_REDUCTION_OF_F_<=_FACTR*EPSMCH'
         if (iback .ge. 10) info = -5
         goto 999
      endif


      do 42 i = 1, n
         r(i) = g(i) - r(i)
  42  continue
      rr = ddot(n,r,1,r,1)
      if (stp .eq. one) then
         dr = gd - gdold
         ddum = -gdold
      else
         dr = (gd - gdold)*stp
         call dscal(n,stp,d,1)
         ddum = -gdold*stp
      endif

      if (dr .le. epsmch*ddum) then
         nskip = nskip + 1
         updatd = .false.
         if (iprint .ge. 1) write (6,1004) dr, ddum
         goto 888
      endif


      updatd = .true.
      iupdat = iupdat + 1


      call matupd(n,m,ws,wy,sy,ss,d,r,itail, iupdat,col,head,theta,rr,dr,stp,dtd)


      call formt(m,wt,sy,ss,col,theta,info)

      if (info .ne. 0) then
         if(iprint .ge. 1) write (6, 1007)
         info = 0
         col = 0
         head = 1
         theta = one
         iupdat = 0
         updatd = .false.
         goto 222
      endif


 888  continue


      goto 222
 999  continue
      call cpu_time(time2)
      time = time2 - time1
      call prn3lb(n,x,f,task,iprint,info,itfile, iter,nfgv,nintol,nskip,nact,sbgnrm, time,nseg,word,iback,stp,xstep,k,&
                 cachyt,sbtime,lnscht)
 1000 continue


      lsave(1)  = prjctd
      lsave(2)  = cnstnd
      lsave(3)  = boxed
      lsave(4)  = updatd

      isave(1)  = nintol
      isave(3)  = itfile
      isave(4)  = iback
      isave(5)  = nskip
      isave(6)  = head
      isave(7)  = col
      isave(8)  = itail
      isave(9)  = iter
      isave(10) = iupdat
      isave(12) = nseg
      isave(13) = nfgv
      isave(14) = info
      isave(15) = ifun
      isave(16) = iword
      isave(17) = nfree
      isave(18) = nact
      isave(19) = ileave
      isave(20) = nenter

      dsave(1)  = theta
      dsave(2)  = fold
      dsave(3)  = tol
      dsave(4)  = dnorm
      dsave(5)  = epsmch
      dsave(6)  = cpu1
      dsave(7)  = cachyt
      dsave(8)  = sbtime
      dsave(9)  = lnscht
      dsave(10) = time1
      dsave(11) = gd
      dsave(12) = stpmx
      dsave(13) = sbgnrm
      dsave(14) = stp
      dsave(15) = gdold
      dsave(16) = dtd

 1001 format (//,'ITERATION ',i5)
 1002 format (/,'At iterate',i5,4x,'f= ',1p,d12.5,4x,'|proj g|= ',1p,d12.5)
 1003 format (2(1x,i4),5x,'-',5x,'-',3x,'-',5x,'-',5x,'-',8x,'-',3x, 1p,2(1x,d10.3))
 1004 format ('  ys=',1p,e10.3,'  -gs=',1p,e10.3,' BFGS update SKIPPED')
 1005 format (/,' Singular triangular system detected;',/,'   refresh the lbfgs memory and restart the iteration.')
 1006 format (/,' Nonpositive definiteness in Cholesky factorization in formk;',/,&
    '   refresh the lbfgs memory and restart the iteration.')
 1007 format (/,&
     ' Nonpositive definiteness in Cholesky factorization in formt;',/,&
     '   refresh the lbfgs memory and restart the iteration.')
 1008 format (/,&
     ' Bad direction in the line search;',/,&
     '   refresh the lbfgs memory and restart the iteration.')

      return

      end


      subroutine active(n, l, u, nbd, x, iwhere, iprint, prjctd, cnstnd, boxed)
        implicit none
      logical          prjctd, cnstnd, boxed
      integer          n, iprint, nbd(n), iwhere(n)
      real(c_double) x(n), l(n), u(n)


      integer          nbdd,i
      real(c_double) zero
      parameter        (zero=0.0d0)


      nbdd = 0
      prjctd = .false.
      cnstnd = .false.
      boxed = .true.


      do 10 i = 1, n
         if (nbd(i) .gt. 0) then
            if (nbd(i) .le. 2 .and. x(i) .le. l(i)) then
               if (x(i) .lt. l(i)) then
                  prjctd = .true.
                  x(i) = l(i)
               endif
               nbdd = nbdd + 1
            else if (nbd(i) .ge. 2 .and. x(i) .ge. u(i)) then
               if (x(i) .gt. u(i)) then
                  prjctd = .true.
                  x(i) = u(i)
               endif
               nbdd = nbdd + 1
            endif
         endif
  10  continue


      do 20 i = 1, n
         if (nbd(i) .ne. 2) boxed = .false.
         if (nbd(i) .eq. 0) then
            iwhere(i) = -1

         else
            cnstnd = .true.
            if (nbd(i) .eq. 2 .and. u(i) - l(i) .le. zero) then
               iwhere(i) = 3
            else
               iwhere(i) = 0
            endif
         endif
  20  continue

      if (iprint .ge. 0) then
         if (prjctd) write (6,*) 'The initial X is infeasible.  Restart with its projection.'
         if (.not. cnstnd) write (6,*) 'This problem is unconstrained.'
      endif

      if (iprint .gt. 0) write (6,1001) nbdd

 1001 format (/,'At X0 ',i9,' variables are exactly at the bounds')

      return

      end


      subroutine bmv(m, sy, wt, col, v, p, info)
        implicit none
      integer m, col, info
      real(c_double) sy(m, m), wt(m, m), v(2*col), p(2*col)


      integer          i,k,i2
      real(c_double) sum

      if (col .eq. 0) return

      p(col + 1) = v(col + 1)
      do 20 i = 2, col
         i2 = col + i
         sum = 0.0d0
         do 10 k = 1, i - 1
            sum = sum + sy(i,k)*v(k)/sy(k,k)
  10     continue
         p(i2) = v(i2) + sum
  20  continue
      call dtrsl(wt,m,col,p(col+1),11,info)
      if (info .ne. 0) return

      do 30 i = 1, col
         p(i) = v(i)/sqrt(sy(i,i))
  30  continue

      call dtrsl(wt,m,col,p(col+1),01,info)
      if (info .ne. 0) return

      do 40 i = 1, col
         p(i) = -p(i)/sqrt(sy(i,i))
  40  continue
      do 60 i = 1, col
         sum = 0.d0
         do 50 k = i + 1, col
            sum = sum + sy(k,i)*p(col+k)/sy(i,i)
  50     continue
         p(i) = p(i) + sum
  60  continue

      return

      end


      subroutine cauchy(n, x, l, u, nbd, g, iorder, iwhere, t, d, xcp, m, wy, ws, sy, wt, theta, col, head, p, c, wbp,&
                       v, nseg, iprint, sbgnrm, info, epsmch)
      implicit none
      integer          n, m, head, col, nseg, iprint, info,  nbd(n), iorder(n), iwhere(n)
      real(c_double) theta, epsmch,&
                      x(n), l(n), u(n), g(n), t(n), d(n), xcp(n),&
                     wy(n, col), ws(n, col), sy(m, m),&
                     wt(m, m), p(2*m), c(2*m), wbp(2*m), v(2*m)


      logical          xlower,xupper,bnded
      integer          i,j,col2,nfree,nbreak,pointr,    ibp,nleft,ibkmin,iter
      real(c_double) f1,f2,dt,dtm,tsum,dibp,zibp,dibp2,bkmin,tu,tl,wmc,wmp,wmw,tj,tj0,neggi,sbgnrm,f2_org!,ddot
      real(c_double) one,zero
      parameter        (one=1.0d0,zero=0.0d0)


      if (sbgnrm .le. zero) then
         if (iprint .ge. 0) write (6,*) 'Subgnorm = 0.  GCP = X.'
         call dcopy(n,x,1,xcp,1)
         return
      endif
      bnded = .true.
      nfree = n + 1
      nbreak = 0
      ibkmin = 0
      bkmin = zero
      col2 = 2*col
      f1 = zero
      if (iprint .ge. 99) write (6,3010)


      do 20 i = 1, col2
         p(i) = zero
  20  continue


      do 50 i = 1, n
         neggi = -g(i)
         if (iwhere(i) .ne. 3 .and. iwhere(i) .ne. -1) then
            if (nbd(i) .le. 2) tl = x(i) - l(i)
            if (nbd(i) .ge. 2) tu = u(i) - x(i)

            xlower = nbd(i) .le. 2 .and. tl .le. zero
            xupper = nbd(i) .ge. 2 .and. tu .le. zero

            iwhere(i) = 0
            if (xlower) then
               if (neggi .le. zero) iwhere(i) = 1
            else if (xupper) then
               if (neggi .ge. zero) iwhere(i) = 2
            else
               if (abs(neggi) .le. zero) iwhere(i) = -3
            endif
         endif
         pointr = head
         if (iwhere(i) .ne. 0 .and. iwhere(i) .ne. -1) then
            d(i) = zero
         else
            d(i) = neggi
            f1 = f1 - neggi*neggi
            do 40 j = 1, col
               p(j) = p(j) +  wy(i,pointr)* neggi
               p(col + j) = p(col + j) + ws(i,pointr)*neggi
               pointr = mod(pointr,m) + 1
  40        continue
            if (nbd(i) .le. 2 .and. nbd(i) .ne. 0  .and. neggi .lt. zero) then
               nbreak = nbreak + 1
               iorder(nbreak) = i
               t(nbreak) = tl/(-neggi)
               if (nbreak .eq. 1 .or. t(nbreak) .lt. bkmin) then
                  bkmin = t(nbreak)
                  ibkmin = nbreak
               endif
            else if (nbd(i) .ge. 2 .and. neggi .gt. zero) then
               nbreak = nbreak + 1
               iorder(nbreak) = i
               t(nbreak) = tu/neggi
               if (nbreak .eq. 1 .or. t(nbreak) .lt. bkmin) then
                  bkmin = t(nbreak)
                  ibkmin = nbreak
               endif
            else
               nfree = nfree - 1
               iorder(nfree) = i
               if (abs(neggi) .gt. zero) bnded = .false.
            endif
         endif
  50  continue


      if (theta .ne. one) then
         call dscal(col,theta,p(col+1),1)
      endif


      call dcopy(n,x,1,xcp,1)

      if (nbreak .eq. 0 .and. nfree .eq. n + 1) then
         if (iprint .gt. 100) write (6,1010) (xcp(i), i = 1, n)
         return
      endif


      do 60 j = 1, col2
         c(j) = zero
  60  continue


      f2 =  -theta*f1
      f2_org  =  f2
      if (col .gt. 0) then
         call bmv(m,sy,wt,col,p,v,info)
         if (info .ne. 0) return
         f2 = f2 - ddot(col2,v,1,p,1)
      endif
      dtm = -f1/f2
      tsum = zero
      nseg = 1
      if (iprint .ge. 99)  write (6,*) 'There are ',nbreak,'  breakpoints '


      if (nbreak .eq. 0) goto 888

      nleft = nbreak
      iter = 1


      tj = zero


 777  continue


      tj0 = tj
      if (iter .eq. 1) then
         tj = bkmin
         ibp = iorder(ibkmin)
      else
         if (iter .eq. 2) then
            if (ibkmin .ne. nbreak) then
               t(ibkmin) = t(nbreak)
               iorder(ibkmin) = iorder(nbreak)
            endif
         endif
         call hpsolb(nleft,t,iorder,iter-2)
         tj = t(nleft)
         ibp = iorder(nleft)
      endif

      dt = tj - tj0

      if (dt .ne. zero .and. iprint .ge. 100) then
         write (6,4011) nseg,f1,f2
         write (6,5010) dt
         write (6,6010) dtm
      endif


      if (dtm .lt. dt) goto 888


      tsum = tsum + dt
      nleft = nleft - 1
      iter = iter + 1
      dibp = d(ibp)
      d(ibp) = zero
      if (dibp .gt. zero) then
         zibp = u(ibp) - x(ibp)
         xcp(ibp) = u(ibp)
         iwhere(ibp) = 2
      else
         zibp = l(ibp) - x(ibp)
         xcp(ibp) = l(ibp)
         iwhere(ibp) = 1
      endif
      if (iprint .ge. 100) write (6,*) 'Variable  ',ibp,'  is fixed.'
      if (nleft .eq. 0 .and. nbreak .eq. n) then
         dtm = dt
         goto 999
      endif


      nseg = nseg + 1
      dibp2 = dibp**2


      f1 = f1 + dt*f2 + dibp2 - theta*dibp*zibp
      f2 = f2 - theta*dibp2

      if (col .gt. 0) then
         call daxpy(col2,dt,p,1,c,1)

         pointr = head
         do 70 j = 1,col
            wbp(j) = wy(ibp,pointr)
            wbp(col + j) = theta*ws(ibp,pointr)
            pointr = mod(pointr,m) + 1
  70     continue

         call bmv(m,sy,wt,col,wbp,v,info)
         if (info .ne. 0) return
         wmc = ddot(col2,c,1,v,1)
         wmp = ddot(col2,p,1,v,1)
         wmw = ddot(col2,wbp,1,v,1)

         call daxpy(col2,-dibp,wbp,1,p,1)

         f1 = f1 + dibp*wmc
         f2 = f2 + 2.0d0*dibp*wmp - dibp2*wmw
      endif

      f2 = max(epsmch*f2_org,f2)
      if (nleft .gt. 0) then
         dtm = -f1/f2
         goto 777
      else if(bnded) then
         f1 = zero
         f2 = zero
         dtm = zero
      else
         dtm = -f1/f2
      endif


 888  continue
      if (iprint .ge. 99) then
         write (6,*)
         write (6,*) 'GCP found in this segment'
         write (6,4010) nseg,f1,f2
         write (6,6010) dtm
      endif
      if (dtm .le. zero) dtm = zero
      tsum = tsum + dtm


      call daxpy(n,tsum,d,1,xcp,1)

 999  continue


      if (col .gt. 0) call daxpy(col2,dtm,p,1,c,1)
      if (iprint .gt. 100) write (6,1010) (xcp(i),i = 1,n)
      if (iprint .ge. 99) write (6,2010)

 1010 format ('Cauchy X =  ',/,(4x,1p,6(1x,d11.4)))
 2010 format (/,'---------------- exit CAUCHY----------------------',/)
 3010 format (/,'---------------- CAUCHY entered-------------------')
 4010 format ('Piece    ',i3,' --f1, f2 at start point ',1p,2(1x,d11.4))
 4011 format (/,'Piece    ',i3,' --f1, f2 at start point ',  1p,2(1x,d11.4))
 5010 format ('Distance to the next break point =  ',1p,d11.4)
 6010 format ('Distance to the stationary point =  ',1p,d11.4)

      return

      end


      subroutine cmprlb(n, m, x, g, ws, wy, sy, wt, z, r, wa, index,  theta, col, head, nfree, cnstnd, info)
        implicit none
      logical          cnstnd
      integer          n, m, col, head, nfree, info, index(n)
      real(c_double) theta, x(n), g(n), z(n), r(n), wa(4*m), ws(n, m), wy(n, m), sy(m, m), wt(m, m)


      integer          i,j,k,pointr
      real(c_double) a1,a2

      if (.not. cnstnd .and. col .gt. 0) then
         do 26 i = 1, n
            r(i) = -g(i)
  26     continue
      else
         do 30 i = 1, nfree
            k = index(i)
            r(i) = -theta*(z(k) - x(k)) - g(k)
  30     continue
         call bmv(m,sy,wt,col,wa(2*m+1),wa(1),info)
         if (info .ne. 0) then
            info = -8
            return
         endif
         pointr = head
         do 34 j = 1, col
            a1 = wa(j)
            a2 = theta*wa(col + j)
            do 32 i = 1, nfree
               k = index(i)
               r(i) = r(i) + wy(k,pointr)*a1 + ws(k,pointr)*a2
  32        continue
            pointr = mod(pointr,m) + 1
  34     continue
      endif

      return

      end


      subroutine errclb(n, m, factr, l, u, nbd, task, info, k)
        implicit none
      character(len=60)     task
      integer          n, m, info, k, nbd(n)
      real(c_double) factr, l(n), u(n)


      integer          i
      real(c_double) one,zero
      parameter        (one=1.0d0,zero=0.0d0)


      if (n .le. 0) task = 'ERROR: N .LE. 0'
      if (m .le. 0) task = 'ERROR: M .LE. 0'
      if (factr .lt. zero) task = 'ERROR: FACTR .LT. 0'


      do 10 i = 1, n
         if (nbd(i) .lt. 0 .or. nbd(i) .gt. 3) then
            task = 'ERROR: INVALID NBD'
            info = -6
            k = i
         endif
         if (nbd(i) .eq. 2) then
            if (l(i) .gt. u(i)) then
               task = 'ERROR: NO FEASIBLE SOLUTION'
               info = -7
               k = i
            endif
         endif
  10  continue

      return

      end


      subroutine formk(n, nsub, ind, nenter, ileave, indx2, iupdat, updatd, wn, wn1, m, ws, wy, sy, theta, col, head, info)
        implicit none
      integer          n, nsub, m, col, head, nenter, ileave, iupdat,  info, ind(n), indx2(n)
      real(c_double) theta, wn(2*m, 2*m), wn1(2*m, 2*m),  ws(n, m), wy(n, m), sy(m, m)
      logical          updatd


      integer          m2,ipntr,jpntr,iy,is,jy,js,is1,js1,k1,i,k, col2,pbegin,pend,dbegin,dend,upcl
      real(c_double) temp1,temp2,temp3,temp4!,ddot
      real(c_double) one,zero
      parameter        (one=1.0d0,zero=0.0d0)


      if (updatd) then
         if (iupdat .gt. m) then
            do 10 jy = 1, m - 1
               js = m + jy
               call dcopy(m-jy,wn1(jy+1,jy+1),1,wn1(jy,jy),1)
               call dcopy(m-jy,wn1(js+1,js+1),1,wn1(js,js),1)
               call dcopy(m-1,wn1(m+2,jy+1),1,wn1(m+1,jy),1)
  10        continue
         endif

         pbegin = 1
         pend = nsub
         dbegin = nsub + 1
         dend = n
         iy = col
         is = m + col
         ipntr = head + col - 1
         if (ipntr .gt. m) ipntr = ipntr - m
         jpntr = head
         do 20 jy = 1, col
            js = m + jy
            temp1 = zero
            temp2 = zero
            temp3 = zero
            do 15 k = pbegin, pend
               k1 = ind(k)
               temp1 = temp1 + wy(k1,ipntr)*wy(k1,jpntr)
  15        continue
            do 16 k = dbegin, dend
               k1 = ind(k)
               temp2 = temp2 + ws(k1,ipntr)*ws(k1,jpntr)
               temp3 = temp3 + ws(k1,ipntr)*wy(k1,jpntr)
  16        continue
            wn1(iy,jy) = temp1
            wn1(is,js) = temp2
            wn1(is,jy) = temp3
            jpntr = mod(jpntr,m) + 1
  20     continue

         jy = col
         jpntr = head + col - 1
         if (jpntr .gt. m) jpntr = jpntr - m
         ipntr = head
         do 30 i = 1, col
            is = m + i
            temp3 = zero
            do 25 k = pbegin, pend
               k1 = ind(k)
               temp3 = temp3 + ws(k1,ipntr)*wy(k1,jpntr)
  25        continue
            ipntr = mod(ipntr,m) + 1
            wn1(is,jy) = temp3
  30     continue
         upcl = col - 1
      else
         upcl = col
      endif

      ipntr = head
      do 45 iy = 1, upcl
         is = m + iy
         jpntr = head
         do 40 jy = 1, iy
            js = m + jy
            temp1 = zero
            temp2 = zero
            temp3 = zero
            temp4 = zero
            do 35 k = 1, nenter
               k1 = indx2(k)
               temp1 = temp1 + wy(k1,ipntr)*wy(k1,jpntr)
               temp2 = temp2 + ws(k1,ipntr)*ws(k1,jpntr)
  35        continue
            do 36 k = ileave, n
               k1 = indx2(k)
               temp3 = temp3 + wy(k1,ipntr)*wy(k1,jpntr)
               temp4 = temp4 + ws(k1,ipntr)*ws(k1,jpntr)
  36        continue
            wn1(iy,jy) = wn1(iy,jy) + temp1 - temp3
            wn1(is,js) = wn1(is,js) - temp2 + temp4
            jpntr = mod(jpntr,m) + 1
  40     continue
         ipntr = mod(ipntr,m) + 1
  45  continue

      ipntr = head
      do 60 is = m + 1, m + upcl
         jpntr = head
         do 55 jy = 1, upcl
            temp1 = zero
            temp3 = zero
            do 50 k = 1, nenter
               k1 = indx2(k)
               temp1 = temp1 + ws(k1,ipntr)*wy(k1,jpntr)
  50        continue
            do 51 k = ileave, n
               k1 = indx2(k)
               temp3 = temp3 + ws(k1,ipntr)*wy(k1,jpntr)
  51        continue
         if (is .le. jy + m) then
               wn1(is,jy) = wn1(is,jy) + temp1 - temp3
            else
               wn1(is,jy) = wn1(is,jy) - temp1 + temp3
            endif
            jpntr = mod(jpntr,m) + 1
  55     continue
         ipntr = mod(ipntr,m) + 1
  60  continue


      m2 = 2*m
      do 70 iy = 1, col
         is = col + iy
         is1 = m + iy
         do 65 jy = 1, iy
            js = col + jy
            js1 = m + jy
            wn(jy,iy) = wn1(iy,jy)/theta
            wn(js,is) = wn1(is1,js1)*theta
  65     continue
         do 66 jy = 1, iy - 1
            wn(jy,is) = -wn1(is1,jy)
  66     continue
         do 67 jy = iy, col
            wn(jy,is) = wn1(is1,jy)
  67     continue
         wn(iy,iy) = wn(iy,iy) + sy(iy,iy)
  70  continue


      call dpofa(wn,m2,col,info)
      if (info .ne. 0) then
         info = -1
         return
      endif
      col2 = 2*col
      do  js = col+1 ,col2
         call dtrsl(wn,m2,col,wn(1,js),11,info)
      end do



      do is = col+1, col2
        do js = is, col2
               wn(is,js) = wn(is,js) + ddot(col,wn(1,is),1,wn(1,js),1)
        end do
      end do


      call dpofa(wn(col+1,col+1),m2,col,info)
      if (info .ne. 0) then
         info = -2
         return
      endif

      return

      end


      subroutine formt(m, wt, sy, ss, col, theta, info)
        implicit none
      integer          m, col, info
      real(c_double) theta, wt(m, m), sy(m, m), ss(m, m)


      integer          i,j,k,k1
      real(c_double) ddum
      real(c_double) zero
      parameter        (zero=0.0d0)



      do  j = 1, col
         wt(1,j) = theta*ss(1,j)
      end do
      do 55 i = 2, col
         do  j = i, col
            k1 = min(i,j) - 1
            ddum  = zero
            do k = 1, k1
               ddum  = ddum + sy(i,k)*sy(j,k)/sy(k,k)
            end do
            wt(i,j) = ddum + theta*ss(i,j)
        end do
  55  continue


      call dpofa(wt,m,col,info)
      if (info .ne. 0) then
         info = -3
      endif

      return

      end


      subroutine freev(n, nfree, index, nenter, ileave, indx2,  iwhere, wrk, updatd, cnstnd, iprint, iter)
        implicit none
      integer n, nfree, nenter, ileave, iprint, iter,  index(n), indx2(n), iwhere(n)
      logical wrk, updatd, cnstnd


      integer iact,i,k

      nenter = 0
      ileave = n + 1
      if (iter .gt. 0 .and. cnstnd) then
         do 20 i = 1, nfree
            k = index(i)


            if (iwhere(k) .gt. 0) then
               ileave = ileave - 1
               indx2(ileave) = k
               if (iprint .ge. 100) write (6,*) 'Variable ',k,' leaves the set of free variables'
            endif
  20     continue
         do 22 i = 1 + nfree, n
            k = index(i)
            if (iwhere(k) .le. 0) then
               nenter = nenter + 1
               indx2(nenter) = k
               if (iprint .ge. 100) write (6,*) 'Variable ',k,' enters the set of free variables'
            endif
  22     continue
         if (iprint .ge. 99) write (6,*)   n+1-ileave,' variables leave; ',nenter,' variables enter'
      endif
      wrk = (ileave .lt. n+1) .or. (nenter .gt. 0) .or. updatd


      nfree = 0
      iact = n + 1
      do 24 i = 1, n
         if (iwhere(i) .le. 0) then
            nfree = nfree + 1
            index(nfree) = i
         else
            iact = iact - 1
            index(iact) = i
         endif
  24  continue
      if (iprint .ge. 99) write (6,*)    nfree,' variables are free at GCP ',iter + 1

      return

      end


      subroutine hpsolb(n, t, iorder, iheap)
        implicit none
      integer          iheap, n, iorder(n)
      real(c_double) t(n)


      integer          i,j,k,indxin,indxou
      real(c_double) ddum,out

      if (iheap .eq. 0) then


         do 20 k = 2, n
            ddum  = t(k)
            indxin = iorder(k)

            i = k
   10       continue
            if (i.gt.1) then
               j = i/2
               if (ddum .lt. t(j)) then
                  t(i) = t(j)
                  iorder(i) = iorder(j)
                  i = j
                  goto 10
               endif
            endif
            t(i) = ddum
            iorder(i) = indxin
   20    continue
      endif


      if (n .gt. 1) then
         i = 1
         out = t(1)
         indxou = iorder(1)
         ddum  = t(n)
         indxin  = iorder(n)

   30    continue
         j = i+i
         if (j .le. n-1) then
            if (t(j+1) .lt. t(j)) j = j+1
            if (t(j) .lt. ddum ) then
               t(i) = t(j)
               iorder(i) = iorder(j)
               i = j
               goto 30
            endif
         endif
         t(i) = ddum
         iorder(i) = indxin


         t(n) = out
         iorder(n) = indxou
      endif

      return

      end


      subroutine lnsrlb(n, l, u, nbd, x, f, fold, gd, gdold, g, d, r, t,&
                       z, stp, dnorm, dtd, xstep, stpmx, iter, ifun,&
                       iback, nfgv, info, task, boxed, cnstnd, csave,&
                       isave, dsave)
      implicit none
      character(len=60)     task, csave
      logical          boxed, cnstnd
      integer          n, iter, ifun, iback, nfgv, info,  nbd(n), isave(2)
      real(c_double) f, fold, gd, gdold, stp, dnorm, dtd, xstep, stpmx, x(n), l(n), u(n), g(n), d(n), r(n), t(n), z(n), dsave(13)

      integer          i
      double           precision a1,a2!,ddot
      real(c_double) one,zero,big
      parameter        (one=1.0d0,zero=0.0d0,big=1.0d+10)
      real(c_double) ftol,gtol,xtol
      parameter        (ftol=1.0d-3,gtol=0.9d0,xtol=0.1d0)

      if (task(1:5) .eq. 'FG_LN') goto 556

      dtd = ddot(n,d,1,d,1)
      dnorm = sqrt(dtd)


      stpmx = big
      if (cnstnd) then
         if (iter .eq. 0) then
            stpmx = one
         else
            do 43 i = 1, n
               a1 = d(i)
               if (nbd(i) .ne. 0) then
                  if (a1 .lt. zero .and. nbd(i) .le. 2) then
                     a2 = l(i) - x(i)
                     if (a2 .ge. zero) then
                        stpmx = zero
                     else if (a1*stpmx .lt. a2) then
                        stpmx = a2/a1
                     endif
                  else if (a1 .gt. zero .and. nbd(i) .ge. 2) then
                     a2 = u(i) - x(i)
                     if (a2 .le. zero) then
                        stpmx = zero
                     else if (a1*stpmx .gt. a2) then
                        stpmx = a2/a1
                     endif
                  endif
               endif
  43        continue
         endif
      endif

      if (iter .eq. 0 .and. .not. boxed) then
         stp = min(one/dnorm, stpmx)
      else
         stp = one
      endif

      call dcopy(n,x,1,t,1)
      call dcopy(n,g,1,r,1)
      fold = f
      ifun = 0
      iback = 0
      csave = 'START'
 556  continue
      gd = ddot(n,g,1,d,1)
      if (ifun .eq. 0) then
         gdold=gd
         if (gd .ge. zero) then
            write(6,*)' ascent direction in projection gd = ', gd
            info = -4
            return
         endif
      endif

      call dcsrch(f,gd,stp,ftol,gtol,xtol,zero,stpmx,csave,isave,dsave)

      xstep = stp*dnorm
      if (csave(1:4) .ne. 'CONV' .and. csave(1:4) .ne. 'WARN') then
         task = 'FG_LNSRCH'
         ifun = ifun + 1
         nfgv = nfgv + 1
         iback = ifun - 1
         if (stp .eq. one) then
            call dcopy(n,z,1,x,1)
         else
            do 41 i = 1, n
               x(i) = stp*d(i) + t(i)
  41        continue
         endif
      else
         task = 'NEW_X'
      endif

      return

      end


      subroutine matupd(n, m, ws, wy, sy, ss, d, r, itail,  iupdat, col, head, theta, rr, dr, stp, dtd)
        implicit none
      integer          n, m, itail, iupdat, col, head
      real(c_double) theta, rr, dr, stp, dtd, d(n), r(n), ws(n, m), wy(n, m), sy(m, m), ss(m, m)


      integer          j,pointr
      ! real(c_double) ddot
      real(c_double) one
      parameter        (one=1.0d0)


      if (iupdat .le. m) then
         col = iupdat
         itail = mod(head+iupdat-2,m) + 1
      else
         itail = mod(itail,m) + 1
         head = mod(head,m) + 1
      endif


      call dcopy(n,d,1,ws(1,itail),1)
      call dcopy(n,r,1,wy(1,itail),1)


      theta = rr/dr


      if (iupdat .gt. m) then
         do 50 j = 1, col - 1
            call dcopy(j,ss(2,j+1),1,ss(1,j),1)
            call dcopy(col-j,sy(j+1,j+1),1,sy(j,j),1)
  50     continue
      endif
      pointr = head
      do 51 j = 1, col - 1
         sy(col,j) = ddot(n,d,1,wy(1,pointr),1)
         ss(j,col) = ddot(n,ws(1,pointr),1,d,1)
         pointr = mod(pointr,m) + 1
  51  continue
      if (stp .eq. one) then
         ss(col,col) = dtd
      else
         ss(col,col) = stp*stp*dtd
      endif
      sy(col,col) = dr

      return

      end


      subroutine prn1lb(n, m, l, u, x, iprint, itfile, epsmch)
        implicit none
      integer n, m, iprint, itfile
      real(c_double) epsmch, x(n), l(n), u(n)


      integer i

      if (iprint .ge. 0) then
         write (6,7001) epsmch
         write (6,*) 'N = ',n,'    M = ',m
         if (iprint .ge. 1) then
            write (itfile,2001) epsmch
            write (itfile,*)'N = ',n,'    M = ',m
            write (itfile,9001)
            if (iprint .gt. 100) then
               write (6,1004) 'L =',(l(i),i = 1,n)
               write (6,1004) 'X0 =',(x(i),i = 1,n)
               write (6,1004) 'U =',(u(i),i = 1,n)
            endif
         endif
      endif

 1004 format (/,a4, 1p, 6(1x,d11.4),/,(4x,1p,6(1x,d11.4)))
 2001 format ('RUNNING THE L-BFGS-B CODE',/,/,&
      'it    = iteration number',/,&
      'nf    = number of function evaluations',/,&
      'nseg  = number of segments explored during the Cauchy search',/,&
      'nact  = number of active bounds at the generalized Cauchy point'&
      ,/,&
      'sub   = manner in which the subspace minimization terminated:'&
      ,/,'        con = converged, bnd = a bound was reached',/,&
      'itls  = number of iterations performed in the line search',/,&
      'stepl = step length used',/,&
      'tstep = norm of the displacement (total step)',/,&
      'projg = norm of the projected gradient',/,&
      'f     = function value',/,/,&
      '           * * *',/,/,&
      'Machine precision =',1p,d10.3)
 7001 format ('RUNNING THE L-BFGS-B CODE',/,/,'           * * *',/,/,'Machine precision =',1p,d10.3)
 9001 format (/,3x,'it',3x,'nf',2x,'nseg',2x,'nact',2x,'sub',2x,'itls',   2x,'stepl',4x,'tstep',5x,'projg',8x,'f')

      return

      end


      subroutine prn2lb(n, x, f, g, iprint, itfile, iter, nfgv, nact, sbgnrm, nseg, word, iword, iback, stp, xstep)
        implicit none
      character(len=3)      word
      integer          n, iprint, itfile, iter, nfgv, nact, nseg, iword, iback
      real(c_double) f, sbgnrm, stp, xstep, x(n), g(n)


      integer i,imod

      if (iword .eq. 0) then
         word = 'con'
      else if (iword .eq. 1) then
         word = 'bnd'
      else if (iword .eq. 5) then
         word = 'TNT'
      else
         word = '---'
      endif
      if (iprint .ge. 99) then
         write (6,*) 'LINE SEARCH',iback,' times; norm of step = ',xstep
         write (6,2001) iter,f,sbgnrm
         if (iprint .gt. 100) then
            write (6,1004) 'X =',(x(i), i = 1, n)
            write (6,1004) 'G =',(g(i), i = 1, n)
         endif
      else if (iprint .gt. 0) then
         imod = mod(iter,iprint)
         if (imod .eq. 0) write (6,2001) iter,f,sbgnrm
      endif
      if (iprint .ge. 1) write (itfile,3001) iter,nfgv,nseg,nact,word,iback,stp,xstep,sbgnrm,f

 1004 format (/,a4, 1p, 6(1x,d11.4),/,(4x,1p,6(1x,d11.4)))
 2001 format (/,'At iterate',i5,4x,'f= ',1p,d12.5,4x,'|proj g|= ',1p,d12.5)
 3001 format(2(1x,i4),2(1x,i5),2x,a3,1x,i4,1p,2(2x,d7.1),1p,2(1x,d10.3))

      return

      end


      subroutine prn3lb(n, x, f, task, iprint, info, itfile,&
                       iter, nfgv, nintol, nskip, nact, sbgnrm,&
                       time, nseg, word, iback, stp, xstep, k,&
                       cachyt, sbtime, lnscht)
        implicit none
      character(len=60)     task
      character(len=3)      word
      integer          n, iprint, info, itfile, iter, nfgv, nintol,  nskip, nact, nseg, iback, k
      real(c_double) f, sbgnrm, time, stp, xstep, cachyt, sbtime, lnscht, x(n)


      integer i

      if (task(1:5) .eq. 'ERROR') goto 999

      if (iprint .ge. 0) then
         write (6,3003)
         write (6,3004)
         write(6,3005) n,iter,nfgv,nintol,nskip,nact,sbgnrm,f
         if (iprint .ge. 100) then
            write (6,1004) 'X =',(x(i),i = 1,n)
         endif
         if (iprint .ge. 1) write (6,*) ' F =',f
      endif
 999  continue
      if (iprint .ge. 0) then
         write (6,3009) task
         if (info .ne. 0) then
            if (info .eq. -1) write (6,9011)
            if (info .eq. -2) write (6,9012)
            if (info .eq. -3) write (6,9013)
            if (info .eq. -4) write (6,9014)
            if (info .eq. -5) write (6,9015)
            if (info .eq. -6) write (6,*)' Input nbd(',k,') is invalid.'
            if (info .eq. -7) write (6,*)' l(',k,') > u(',k,').  No feasible solution.'
            if (info .eq. -8) write (6,9018)
            if (info .eq. -9) write (6,9019)
         endif
         if (iprint .ge. 1) write (6,3007) cachyt,sbtime,lnscht
         write (6,3008) time
         if (iprint .ge. 1) then
            if (info .eq. -4 .or. info .eq. -9) then
               write (itfile,3002) iter,nfgv,nseg,nact,word,iback,stp,xstep
            endif
            write (itfile,3009) task
            if (info .ne. 0) then
               if (info .eq. -1) write (itfile,9011)
               if (info .eq. -2) write (itfile,9012)
               if (info .eq. -3) write (itfile,9013)
               if (info .eq. -4) write (itfile,9014)
               if (info .eq. -5) write (itfile,9015)
               if (info .eq. -8) write (itfile,9018)
               if (info .eq. -9) write (itfile,9019)
            endif
            write (itfile,3008) time
         endif
      endif

 1004 format (/,a4, 1p, 6(1x,d11.4),/,(4x,1p,6(1x,d11.4)))
 3002 format(2(1x,i4),2(1x,i5),2x,a3,1x,i4,1p,2(2x,d7.1),6x,'-',10x,'-')
 3003 format (/,&
      '           * * *',/,/,&
      'Tit   = total number of iterations',/,&
      'Tnf   = total number of function evaluations',/,&
      'Tnint = total number of segments explored during',&
                ' Cauchy searches',/,&
      'Skip  = number of BFGS updates skipped',/,&
      'Nact  = number of active bounds at final generalized',&
               ' Cauchy point',/,&
      'Projg = norm of the final projected gradient',/,&
      'F     = final function value',/,/,&
      '           * * *')
 3004 format (/,3x,'N',4x,'Tit',5x,'Tnf',2x,'Tnint',2x, 'Skip',2x,'Nact',5x,'Projg',8x,'F')
 3005 format (i5,2(1x,i6),(1x,i6),(2x,i4),(1x,i5),1p,2(2x,d10.3))
 3007 format (/,' Cauchy                time',1p,e10.3,' seconds.',/&
            ' Subspace minimization time',1p,e10.3,' seconds.',/&
             ' Line search           time',1p,e10.3,' seconds.')
 3008 format (/,' Total User time',1p,e10.3,' seconds.',/)
 3009 format (/,a60)
 9011 format (/,' Matrix in 1st Cholesky factorization in formk is not Pos. Def.')
 9012 format (/,' Matrix in 2st Cholesky factorization in formk is not Pos. Def.')
 9013 format (/,' Matrix in the Cholesky factorization in formt is not Pos. Def.')
 9014 format (/,&
     ' Derivative >= 0, backtracking line search impossible.',/,&
     '   Previous x, f and g restored.',/,&
     ' Possible causes: 1 error in function or gradient evaluation;',/,&
     '                  2 rounding errors dominate computation.')
 9015 format (/,&
     ' Warning:  more than 10 function and gradient',/,&
     '   evaluations in the last line search.  Termination',/,&
     '   may possibly be caused by a bad search direction.')
 9018 format (/,' The triangular system is singular.')
 9019 format (/,&
     ' Line search cannot locate an adequate point after 20 function',/&
     ,'  and gradient evaluations.  Previous x, f and g restored.',/,&
     ' Possible causes: 1 error in function or gradient evaluation;',/,&
     '                  2 rounding error dominate computation.')

      return

      end


      subroutine projgr(n, l, u, nbd, x, g, sbgnrm)
        implicit none
      integer          n, nbd(n)
      real(c_double) sbgnrm, x(n), l(n), u(n), g(n)


      integer i
      real(c_double) gi
      real(c_double) one,zero
      parameter        (one=1.0d0,zero=0.0d0)

      sbgnrm = zero
      do 15 i = 1, n
        gi = g(i)
        if (nbd(i) .ne. 0) then
           if (gi .lt. zero) then
              if (nbd(i) .ge. 2) gi = max((x(i)-u(i)),gi)
           else
              if (nbd(i) .le. 2) gi = min((x(i)-l(i)),gi)
           endif
        endif
        sbgnrm = max(sbgnrm,abs(gi))
  15  continue

      return

      end


      subroutine subsm ( n, m, nsub, ind, l, u, nbd, x, d, xp, ws, wy, theta, xx, gg, col, head, iword, wv, wn, iprint, info )
      implicit none
      integer          n, m, nsub, col, head, iword, iprint, info, ind(nsub), nbd(n)
      real(c_double) theta, l(n), u(n), x(n), d(n), xp(n), xx(n), gg(n), ws(n, m), wy(n, m), wv(2*m), wn(2*m, 2*m)




      integer          pointr,m2,col2,ibd,jy,js,i,j,k
      real(c_double) alpha, xk, dk, temp1, temp2
      real(c_double) one,zero
      parameter        (one=1.0d0,zero=0.0d0)
      real(c_double) dd_p

      if (nsub .le. 0) return
      if (iprint .ge. 99) write (6,1001)


      pointr = head
      do 20 i = 1, col
         temp1 = zero
         temp2 = zero
         do 10 j = 1, nsub
            k = ind(j)
            temp1 = temp1 + wy(k,pointr)*d(j)
            temp2 = temp2 + ws(k,pointr)*d(j)
  10     continue
         wv(i) = temp1
         wv(col + i) = theta*temp2
         pointr = mod(pointr,m) + 1
  20  continue


      m2 = 2*m
      col2 = 2*col
      call dtrsl(wn,m2,col2,wv,11,info)
      if (info .ne. 0) return
      do 25 i = 1, col
         wv(i) = -wv(i)
  25     continue
      call dtrsl(wn,m2,col2,wv,01,info)
      if (info .ne. 0) return


      pointr = head
      do 40 jy = 1, col
         js = col + jy
         do 30 i = 1, nsub
            k = ind(i)
            d(i) = d(i) + wy(k,pointr)*wv(jy)/theta   + ws(k,pointr)*wv(js)
  30     continue
         pointr = mod(pointr,m) + 1
  40  continue

      call dscal( nsub, one/theta, d, 1 )

      iword = 0

      call dcopy ( n, x, 1, xp, 1 )
      do 50 i=1, nsub
         k  = ind(i)
         dk = d(i)
         xk = x(k)
         if ( nbd(k) .ne. 0 ) then
            if ( nbd(k).eq.1 ) then          ! lower bounds only
               x(k) = max( l(k), xk + dk )
               if ( x(k).eq.l(k) ) iword = 1
            else
               if ( nbd(k).eq.2 ) then       ! upper and lower bounds
                  xk   = max( l(k), xk + dk )
                  x(k) = min( u(k), xk )
                  if ( x(k).eq.l(k) .or. x(k).eq.u(k) ) iword = 1
               else
                  if ( nbd(k).eq.3 ) then    ! upper bounds only
                     x(k) = min( u(k), xk + dk )
                     if ( x(k).eq.u(k) ) iword = 1
                  end if
               end if
            end if
         else                                ! free variables
            x(k) = xk + dk
         end if
 50   continue
      if ( iword.eq.0 ) then
         go to 911
      end if
      dd_p = zero
      do 55 i=1, n
         dd_p  = dd_p + (x(i) - xx(i))*gg(i)
 55   continue
      if ( dd_p .gt.zero ) then
         call dcopy( n, xp, 1, x, 1 )
         write(6,*) ' Positive dir derivative in projection '
         write(6,*) ' Using the backtracking step '
      else
         go to 911
      endif
      alpha = one
      temp1 = alpha
      ibd   = 0
      do 60 i = 1, nsub
         k = ind(i)
         dk = d(i)
         if (nbd(k) .ne. 0) then
            if (dk .lt. zero .and. nbd(k) .le. 2) then
               temp2 = l(k) - x(k)
               if (temp2 .ge. zero) then
                  temp1 = zero
               else if (dk*alpha .lt. temp2) then
                  temp1 = temp2/dk
               endif
            else if (dk .gt. zero .and. nbd(k) .ge. 2) then
               temp2 = u(k) - x(k)
               if (temp2 .le. zero) then
                  temp1 = zero
               else if (dk*alpha .gt. temp2) then
                  temp1 = temp2/dk
               endif
            endif
            if (temp1 .lt. alpha) then
               alpha = temp1
               ibd = i
            endif
         endif
 60   continue

      if (alpha .lt. one) then
         dk = d(ibd)
         k = ind(ibd)
         if (dk .gt. zero) then
            x(k) = u(k)
            d(ibd) = zero
         else if (dk .lt. zero) then
            x(k) = l(k)
            d(ibd) = zero
         endif
      endif
      do 70 i = 1, nsub
         k    = ind(i)
         x(k) = x(k) + alpha*d(i)
 70   continue
 911  continue

      if (iprint .ge. 99) write (6,1004)

 1001 format (/,'----------------SUBSM entered-----------------',/)
 1004 format (/,'----------------exit SUBSM --------------------',/)

      return

      end

      subroutine dcsrch(f,g,stp,ftol,gtol,xtol,stpmin,stpmax, task,isave,dsave)
        implicit none
      character*(*) task
      integer isave(2)
      real(c_double) f,g,stp,ftol,gtol,xtol,stpmin,stpmax
      real(c_double) dsave(13)
      real(c_double) zero,p5,p66
      parameter(zero=0.0d0,p5=0.5d0,p66=0.66d0)
      real(c_double) xtrapl,xtrapu
      parameter(xtrapl=1.1d0,xtrapu=4.0d0)

      logical brackt
      integer stage
      real(c_double) finit,ftest,fm,fx,fxm,fy,fym,ginit,gtest, gm,gx,gxm,gy,gym,stx,sty,stmin,stmax,width,width1


      if (task(1:5) .eq. 'START') then


         if (stp .lt. stpmin) task = 'ERROR: STP .LT. STPMIN'
         if (stp .gt. stpmax) task = 'ERROR: STP .GT. STPMAX'
         if (g .ge. zero) task = 'ERROR: INITIAL G .GE. ZERO'
         if (ftol .lt. zero) task = 'ERROR: FTOL .LT. ZERO'
         if (gtol .lt. zero) task = 'ERROR: GTOL .LT. ZERO'
         if (xtol .lt. zero) task = 'ERROR: XTOL .LT. ZERO'
         if (stpmin .lt. zero) task = 'ERROR: STPMIN .LT. ZERO'
         if (stpmax .lt. stpmin) task = 'ERROR: STPMAX .LT. STPMIN'


         if (task(1:5) .eq. 'ERROR') return


         brackt = .false.
         stage = 1
         finit = f
         ginit = g
         gtest = ftol*ginit
         width = stpmax - stpmin
         width1 = width/p5


         stx = zero
         fx = finit
         gx = ginit
         sty = zero
         fy = finit
         gy = ginit
         stmin = zero
         stmax = stp + xtrapu*stp
         task = 'FG'

         goto 1000

      else


         if (isave(1) .eq. 1) then
            brackt = .true.
         else
            brackt = .false.
         endif
         stage = isave(2)
         ginit = dsave(1)
         gtest = dsave(2)
         gx = dsave(3)
         gy = dsave(4)
         finit = dsave(5)
         fx = dsave(6)
         fy = dsave(7)
         stx = dsave(8)
         sty = dsave(9)
         stmin = dsave(10)
         stmax = dsave(11)
         width = dsave(12)
         width1 = dsave(13)

      endif


      ftest = finit + stp*gtest
      if (stage .eq. 1 .and. f .le. ftest .and. g .ge. zero)   stage = 2


      if (brackt .and. (stp .le. stmin .or. stp .ge. stmax))  task = 'WARNING: ROUNDING ERRORS PREVENT PROGRESS'
      if (brackt .and. stmax - stmin .le. xtol*stmax)          task = 'WARNING: XTOL TEST SATISFIED'
      if (stp .eq. stpmax .and. f .le. ftest .and. g .le. gtest)   task = 'WARNING: STP = STPMAX'
      if (stp .eq. stpmin .and. (f .gt. ftest .or. g .ge. gtest))   task = 'WARNING: STP = STPMIN'


      if (f .le. ftest .and. abs(g) .le. gtol*(-ginit))   task = 'CONVERGENCE'


      if (task(1:4) .eq. 'WARN' .or. task(1:4) .eq. 'CONV') goto 1000


      if (stage .eq. 1 .and. f .le. fx .and. f .gt. ftest) then


         fm = f - stp*gtest
         fxm = fx - stx*gtest
         fym = fy - sty*gtest
         gm = g - gtest
         gxm = gx - gtest
         gym = gy - gtest


         call dcstep(stx,fxm,gxm,sty,fym,gym,stp,fm,gm, brackt,stmin,stmax)


         fx = fxm + stx*gtest
         fy = fym + sty*gtest
         gx = gxm + gtest
         gy = gym + gtest

      else


        call dcstep(stx,fx,gx,sty,fy,gy,stp,f,g, brackt,stmin,stmax)

      endif


      if (brackt) then
         if (abs(sty-stx) .ge. p66*width1) stp = stx + p5*(sty - stx)
         width1 = width
         width = abs(sty-stx)
      endif


      if (brackt) then
         stmin = min(stx,sty)
         stmax = max(stx,sty)
      else
         stmin = stp + xtrapl*(stp - stx)
         stmax = stp + xtrapu*(stp - stx)
      endif


      stp = max(stp,stpmin)
      stp = min(stp,stpmax)


      if (brackt .and. (stp .le. stmin .or. stp .ge. stmax)   .or. (brackt .and. stmax-stmin .le. xtol*stmax)) stp = stx


      task = 'FG'

 1000 continue


      if (brackt) then
         isave(1) = 1
      else
         isave(1) = 0
      endif
      isave(2) = stage
      dsave(1) =  ginit
      dsave(2) =  gtest
      dsave(3) =  gx
      dsave(4) =  gy
      dsave(5) =  finit
      dsave(6) =  fx
      dsave(7) =  fy
      dsave(8) =  stx
      dsave(9) =  sty
      dsave(10) = stmin
      dsave(11) = stmax
      dsave(12) = width
      dsave(13) = width1

      return
      end


      subroutine dcstep(stx,fx,dx,sty,fy,dy,stp,fp,dp,brackt,  stpmin,stpmax)
        implicit none
      logical brackt
      real(c_double) stx,fx,dx,sty,fy,dy,stp,fp,dp,stpmin,stpmax

      real(c_double) zero,p66,two,three
      parameter(zero=0.0d0,p66=0.66d0,two=2.0d0,three=3.0d0)

      real(c_double) gamma,p,q,r,s,sgnd,stpc,stpf,stpq,theta

      sgnd = dp*(dx/abs(dx))


      if (fp .gt. fx) then
         theta = three*(fx - fp)/(stp - stx) + dx + dp
         s = max(abs(theta),abs(dx),abs(dp))
         gamma = s*sqrt((theta/s)**2 - (dx/s)*(dp/s))
         if (stp .lt. stx) gamma = -gamma
         p = (gamma - dx) + theta
         q = ((gamma - dx) + gamma) + dp
         r = p/q
         stpc = stx + r*(stp - stx)
         stpq = stx + ((dx/((fx - fp)/(stp - stx) + dx))/two)*(stp - stx)
         if (abs(stpc-stx) .lt. abs(stpq-stx)) then
            stpf = stpc
         else
            stpf = stpc + (stpq - stpc)/two
         endif
         brackt = .true.


      else if (sgnd .lt. zero) then
         theta = three*(fx - fp)/(stp - stx) + dx + dp
         s = max(abs(theta),abs(dx),abs(dp))
         gamma = s*sqrt((theta/s)**2 - (dx/s)*(dp/s))
         if (stp .gt. stx) gamma = -gamma
         p = (gamma - dp) + theta
         q = ((gamma - dp) + gamma) + dx
         r = p/q
         stpc = stp + r*(stx - stp)
         stpq = stp + (dp/(dp - dx))*(stx - stp)
         if (abs(stpc-stp) .gt. abs(stpq-stp)) then
            stpf = stpc
         else
            stpf = stpq
         endif
         brackt = .true.


      else if (abs(dp) .lt. abs(dx)) then


         theta = three*(fx - fp)/(stp - stx) + dx + dp
         s = max(abs(theta),abs(dx),abs(dp))


         gamma = s*sqrt(max(zero,(theta/s)**2-(dx/s)*(dp/s)))
         if (stp .gt. stx) gamma = -gamma
         p = (gamma - dp) + theta
         q = (gamma + (dx - dp)) + gamma
         r = p/q
         if (r .lt. zero .and. gamma .ne. zero) then
            stpc = stp + r*(stx - stp)
         else if (stp .gt. stx) then
            stpc = stpmax
         else
            stpc = stpmin
         endif
         stpq = stp + (dp/(dp - dx))*(stx - stp)

         if (brackt) then


            if (abs(stpc-stp) .lt. abs(stpq-stp)) then
               stpf = stpc
            else
               stpf = stpq
            endif
            if (stp .gt. stx) then
               stpf = min(stp+p66*(sty-stp),stpf)
            else
               stpf = max(stp+p66*(sty-stp),stpf)
            endif
         else


            if (abs(stpc-stp) .gt. abs(stpq-stp)) then
               stpf = stpc
            else
               stpf = stpq
            endif
            stpf = min(stpmax,stpf)
            stpf = max(stpmin,stpf)
         endif


      else
         if (brackt) then
            theta = three*(fp - fy)/(sty - stp) + dy + dp
            s = max(abs(theta),abs(dy),abs(dp))
            gamma = s*sqrt((theta/s)**2 - (dy/s)*(dp/s))
            if (stp .gt. sty) gamma = -gamma
            p = (gamma - dp) + theta
            q = ((gamma - dp) + gamma) + dy
            r = p/q
            stpc = stp + r*(sty - stp)
            stpf = stpc
         else if (stp .gt. stx) then
            stpf = stpmax
         else
            stpf = stpmin
         endif
      endif


      if (fp .gt. fx) then
         sty = stp
         fy = fp
         dy = dp
      else
         if (sgnd .lt. zero) then
            sty = stx
            fy = fx
            dy = dx
         endif
         stx = stp
         fx = fp
         dx = dp
      endif


      stp = stpf

      return
      end
end module module_lbfgs

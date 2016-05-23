MODULE mathematica
    use precision_kinds, only: dp
    ! This module implements several usefull functions of Mathematica
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: chop, TriLinearInterpolation, UTest_TrilinearInterpolation, floorNode, UTest_floorNode, ceilingNode, &
    UTest_ceilingNode, distToFloorNode, UTest_distToFloorNode, deduce_optimal_histogram_properties, &
    splint, spline, akima_spline, fact

    real(dp), parameter :: fact(0:34) = [ 1.0000000000000000     ,&
                                          1.0000000000000000     ,&
                                          2.0000000000000000     ,&
                                          6.0000000000000000     ,&
                                          24.000000000000000     ,&
                                          120.00000000000000     ,&
                                          720.00000000000000     ,&
                                          5040.0000000000000     ,&
                                          40320.000000000000     ,&
                                          362880.00000000000     ,&
                                          3628800.0000000000     ,&
                                          39916800.000000000     ,&
                                          479001600.00000000     ,&
                                          6227020800.0000000     ,&
                                          87178291200.000000     ,&
                                          1307674368000.0000     ,&
                                          20922789888000.000     ,&
                                          355687428096000.00     ,&
                                          6402373705728000.0     ,&
                                          1.2164510040883200E+017,&
                                          2.4329020081766400E+018,&
                                          5.1090942171709440E+019,&
                                          1.1240007277776077E+021,&
                                          2.5852016738884978E+022,&
                                          6.2044840173323941E+023,&
                                          1.5511210043330986E+025,&
                                          4.0329146112660565E+026,&
                                          1.0888869450418352E+028,&
                                          3.0488834461171384E+029,&
                                          8.8417619937397008E+030,&
                                          2.6525285981219103E+032,&
                                          8.2228386541779224E+033,&
                                          2.6313083693369352E+035,&
                                          8.6833176188118859E+036,&
                                          2.9523279903960412E+038 ]


CONTAINS

    !===============================================================================================================================
    PURE FUNCTION chop(x,delta)
        ! see http://reference.wolfram.com/mathematica/ref/Chop.html
        ! It replaces numbers smaller in absolute magnitude than delta by 0.
        ! chop uses a default tolerance of 10._dp**(-10)
        IMPLICIT NONE
        REAL(dp) :: chop
        REAL(dp), INTENT(IN) :: x
        REAL(dp), OPTIONAL, INTENT(IN) :: delta
        REAL(dp), PARAMETER :: defaultdelta=10._dp**(-10)
        REAL(dp) :: d
        IF (PRESENT(delta)) THEN
            d=delta
        ELSE
            d=defaultdelta
        END IF
        IF (abs(x)<=d) THEN
            chop=0._dp
        ELSE
            chop=x
        END IF
    END FUNCTION chop


    !===============================================================================================================================
    PURE FUNCTION TriLinearInterpolation (cube,x)
        ! Returns the value at position x(1:3) within the cube. The value is known at each corner of the cube.
        ! It is thus an interpolation of value at the corners the cube to a point inside the cube.
        IMPLICIT NONE
        REAL(dp), INTENT(IN) :: cube(0:1,0:1,0:1), x(1:3)
        REAL(dp) :: TriLinearInterpolation
        IF( ALL(cube==cube(0,0,0)) ) THEN ! homogeneous case
            TrilinearInterpolation = cube(0,0,0)
        ELSE
            TriLinearInterpolation = cube(0,0,0) * (1._dp-x(1)) * (1._dp-x(2)) * (1._dp-x(3)) &
            +cube(1,0,0) * x(1) * (1._dp-x(2)) * (1._dp-x(3)) &
            +cube(0,1,0) * (1._dp-x(1)) * x(2) * (1._dp-x(3)) &
            +cube(0,0,1) * (1._dp-x(1)) * (1._dp-x(2)) * x(3) &
            +cube(1,0,1) * x(1) * (1._dp-x(2)) * x(3) &
            +cube(0,1,1) * (1._dp-x(1)) * x(2) * x(3) &
            +cube(1,1,0) * x(1) * x(2) * (1._dp-x(3)) &
            +cube(1,1,1) * x(1) * x(2) * x(3)
        END IF
    END FUNCTION TriLinearInterpolation


    !===============================================================================================================================
    SUBROUTINE UTest_TrilinearInterpolation
        ! Tests the pure function TriLinearInterpolation where result is known:
        ! - if the point is one of the corners
        ! - if it is on the center of the cube
        ! Then it tests that no answer is higher or lower than the maximum or minimum value of any corner.
        IMPLICIT NONE
        REAL(dp) :: x(1:3)
        REAL(dp), PARAMETER :: z=0._dp, o=1.0_dp
        LOGICAL, SAVE :: alreadydone=.FALSE.
        REAL(dp) :: cube(0:1,0:1,0:1), t0, t1
        IF (alreadydone) RETURN
        CALL CPU_TIME(t0)
        CALL RANDOM_NUMBER(cube)
        cube = cube * 1000._dp
        IF( TriLinearInterpolation(cube,[z,z,z]) /= cube(0,0,0) .OR.&
        TriLinearInterpolation(cube,[o,z,z]) /= cube(1,0,0) .OR.&
        TriLinearInterpolation(cube,[z,o,z]) /= cube(0,1,0) .OR.&
        TriLinearInterpolation(cube,[z,z,o]) /= cube(0,0,1) .OR.&
        TriLinearInterpolation(cube,[o,o,z]) /= cube(1,1,0) .OR.&
        TriLinearInterpolation(cube,[o,z,o]) /= cube(1,0,1) .OR.&
        TriLinearInterpolation(cube,[z,o,o]) /= cube(0,1,1) .OR.&
        TriLinearInterpolation(cube,[o,o,o]) /= cube(1,1,1) .OR.&
        ABS(TriLinearInterpolation(cube,[o,o,o]/2._dp) - SUM(cube)/8._dp )>EPSILON(1.0_dp) ) THEN
        STOP "Problem detected in UTest_TriLinearInterpolation"
    END IF
    CALL CPU_TIME(t1)
    DO WHILE(t1-t0<0.005_dp) ! Test for 5 ms
        CALL RANDOM_NUMBER(x)
        CALL RANDOM_NUMBER(cube)
        IF ( TriLinearInterpolation(cube,x) < MINVAL(cube) &
        .OR. TriLinearInterpolation(cube,x) > MAXVAL(cube) ) STOP "Problem detected in UTest_TriLinearInterpolation"
        CALL CPU_TIME(t1)
    END DO
    alreadydone = .TRUE.
END SUBROUTINE UTest_TrilinearInterpolation


!===============================================================================================================================
PURE FUNCTION floorNode(gridnode,gridlen,x,pbc)
    IMPLICIT NONE
    INTEGER :: floorNode(3)
    INTEGER, INTENT(IN) :: gridnode(3)
    REAL(dp), INTENT(IN) :: gridlen(3), x(3)
    LOGICAL, INTENT(IN) :: pbc ! periodic boundary counditions
    if(pbc) floorNode = FLOOR(MODULO(x,gridlen)/(gridlen/REAL(gridnode))) +1
END FUNCTION floorNode


!===============================================================================================================================
SUBROUTINE UTest_floorNode
    IMPLICIT NONE
    LOGICAL, SAVE :: alreadydone=.FALSE.
    REAL(dp), PARAMETER :: z=0._dp, o=1.0_dp
    REAL(dp) :: gridlen(3)
    IF (alreadydone) RETURN
    CALL RANDOM_NUMBER(gridlen)
    IF( ANY( floorNode([1,1,1],gridlen*100,[z,z,z],.TRUE.) /=[1,1,1]) ) STOP "problem 1 in UTest_floorNode"
    IF( ANY( floorNode(INT(gridlen*1000)*10,gridlen*100,[z,z,z],.TRUE.) /=[1,1,1]) ) STOP "problem 2 in UTest_floorNode"
    IF( ANY( floorNode([100,1,1],[50._dp,o,o],[51._dp,z,z],.TRUE.) /=[3,1,1]) ) STOP "problem 3 in UTest_floorNode"
    IF( ANY( floorNode([100,1,1],[50._dp,o,o],[49.999_dp,z,z],.TRUE.) /=[100,1,1]) ) STOP "problem 4 in UTest_floorNode"
    IF( ANY( floorNode([100,1,1],[50._dp,o,o],[50._dp,z,z],.TRUE.) /=[1,1,1]) ) STOP "problem 5 in UTest_floorNode"
    alreadydone=.TRUE.
END SUBROUTINE UTest_floorNode


!===============================================================================================================================
PURE FUNCTION ceilingNode(gridnode,gridlen,x,pbc)
    IMPLICIT NONE
    INTEGER :: ceilingNode(3)
    INTEGER, INTENT(IN) :: gridnode(3)
    REAL(dp), INTENT(IN) :: gridlen(3), x(3)
    LOGICAL, INTENT(IN) :: pbc ! periodic boundary counditions
    ceilingNode = MODULO( floorNode(gridnode,gridlen,x,pbc)  ,gridnode) +1
END FUNCTION ceilingNode


!===============================================================================================================================
SUBROUTINE UTest_ceilingNode
    IMPLICIT NONE
    LOGICAL, SAVE :: alreadydone=.FALSE.
    REAL(dp), PARAMETER :: z=0._dp
    INTEGER :: gridnode(3)
    REAL(dp) :: gridlen(3), x(3)
    IF (alreadydone) RETURN

    ! Test 1, for x at origin, i.e., obviously in first bin, i.e., in [1,1,1]. Ceiling Node should be [2,2,2].
    gridlen=[100._dp,100._dp,100._dp]
    gridnode = [100,100,100]
    x = [z,z,z]
    IF ( ANY( ceilingNode(gridnode,gridlen,x,pbc=.TRUE.) /= [2,2,2] ) ) THEN
        STOP "Test 1 of UTest_ceilingNode failed"
    END IF

    ! Test 2, for x at end of supercell, i.e., at gridlen, x is gridlen==[0,0,0], so ceilingNode should be the same as Test 1
    x = gridlen
    IF ( ANY( ceilingNode(gridnode,gridlen,x,pbc=.TRUE.) /= [2,2,2] ) ) THEN
        STOP "Test 2 of UTest_ceilingNode failed"
    END IF

    ! Test 3, for x just after the origin, it is obviously in first bin again, so ceiling should be [2,2,2] again
    x = EPSILON(1.0_dp)
    IF ( ANY( ceilingNode(gridnode,gridlen,x,pbc=.TRUE.) /= [2,2,2] ) ) THEN
        STOP "Test 3 of UTest_ceilingNode failed"
    END IF

    ! Test 4: For x just below gridlen, ceiling should be gridnode
    x = gridlen - gridlen*EPSILON(1.0_dp)
    IF ( ANY( ceilingNode(gridnode,gridlen,x,pbc=.TRUE.) /= [1,1,1] ) ) THEN
        PRINT*,ceilingNode(gridnode,gridlen,x,pbc=.TRUE.)
        STOP "Test 4 of UTest_ceilingNode failed"
    END IF

    alreadydone=.TRUE.
END SUBROUTINE UTest_ceilingNode


!===============================================================================================================================
PURE FUNCTION distToFloorNode(gridnode,gridlen,x,pbc)
    ! Given a grid (number of nodes per direction and length in Angstroms per direction),
    ! returns the distance to floor node in  grid units, i.e., in dx.
    ! 0._dp <= distToFloorNode < 1._dp
    IMPLICIT NONE
    REAL(dp) :: distToFloorNode(3), xfloor(3)
    INTEGER, INTENT(IN) :: gridnode(3)
    REAL(dp), INTENT(IN) :: gridlen(3), x(3)
    LOGICAL, INTENT(IN) :: pbc ! periodic boundary counditions
    if( pbc) xfloor = ABS(  x/(gridlen/REAL(gridnode)) - FLOOR(x/(gridlen/REAL(gridnode)))  )
    distToFloorNode = xfloor
    !~         distToFloorNode = MIN(&
    !~                                 xfloor,&
    !~                                 ABS(1._dp-xfloor)&
    !~                             )
END FUNCTION distToFloorNode


!===============================================================================================================================
SUBROUTINE UTest_distToFloorNode
    implicit none
    REAL(dp), PARAMETER :: z=0._dp, o=1.0_dp

    IF( ANY(   distToFloorNode([10,10,10],[10._dp,10._dp,10._dp],[z,z,z],pbc=.TRUE.) /= [z,z,z] )) THEN
        STOP "UTest_distToFloorNode: Test 1 Failed."
    END IF

    IF( ANY(   distToFloorNode([10,10,10],[10._dp,10._dp,10._dp],[1._dp,1._dp,1._dp],pbc=.TRUE.) /= [z,z,z] )) THEN
        STOP "UTest_distToFloorNode: Test 2 Failed."
    END IF

    IF( ANY(   distToFloorNode([10,10,10],[10._dp,10._dp,10._dp],[2._dp,2._dp,2._dp],pbc=.TRUE.) /= [z,z,z]     )) THEN
        STOP "UTest_distToFloorNode: Test 3 Failed."
    END IF

    IF( ANY(   distToFloorNode([10,10,10],[10._dp,10._dp,10._dp],[10._dp,10._dp,10._dp],pbc=.TRUE.) /= [z,z,z]  )) THEN
        STOP "UTest_distToFloorNode: Test 4 Failed."
    END IF

    IF( ANY(   distToFloorNode([10,10,10],[10._dp,10._dp,10._dp],[11._dp,11._dp,11._dp],pbc=.TRUE.) /= [z,z,z]  )) THEN
        STOP "UTest_distToFloorNode: Test 5 Failed."
    END IF

    IF( ANY(   distToFloorNode([10,10,10],[10._dp,10._dp,10._dp],[2._dp,3._dp,4._dp],pbc=.TRUE.) /= [z,z,z]     )) THEN
        STOP "UTest_distToFloorNode: Test 6 Failed."
    END IF

    IF( ANY(   distToFloorNode([10,10,10],[10._dp,10._dp,10._dp],[10._dp,10._dp,10._dp],.TRUE.) /= [z,z,z]  )) THEN
        STOP "UTest_distToFloorNode: Test 7 Failed."
    END IF

    IF( ANY(   distToFloorNode([10,10,10],[10._dp,10._dp,10._dp],[0.5_dp,0.5_dp,0.5_dp],.TRUE.) /= [o,o,o]/2.0_dp  )) THEN
        STOP "UTest_distToFloorNode: Test 8 Failed."
    END IF

    IF( ANY(   distToFloorNode([10,10,10],[100._dp,100._dp,100._dp],[50._dp,50._dp,50._dp],.TRUE.) /= [z,z,z]   )) THEN
        STOP "UTest_distToFloorNode: Test 9 Failed."
    END IF

    IF( ANY(   distToFloorNode([10,10,10],[100._dp,100._dp,100._dp],[55._dp,55._dp,55._dp],.TRUE.) /= [0.5_dp,0.5_dp,0.5_dp] ))&
    THEN
    STOP "UTest_distToFloorNode: Test 10 Failed."
END IF

IF( ANY(   distToFloorNode([10,10,10],[100._dp,100._dp,100._dp],[15._dp,15._dp,15._dp],.TRUE.) /= [0.5_dp,0.5_dp,0.5_dp] ))&
THEN
STOP "UTest_distToFloorNode: Test 11 Failed."
END IF

IF( ANY(   distToFloorNode([10,10,10],[100._dp,100._dp,100._dp],[51._dp,51._dp,51._dp],.TRUE.) - [0.1_dp,0.1_dp,0.1_dp] &
> EPSILON(1.0_dp)*[100._dp,100._dp,100._dp])) THEN
STOP "UTest_distToFloorNode: Test 12 Failed."
END IF

IF( ANY(   distToFloorNode([10,10,10],[10._dp,10._dp,10._dp],[13._dp,13._dp,13._dp],.TRUE.) - [z,z,z] &
> EPSILON(1.0_dp)*[10._dp,10._dp,10._dp] )) THEN
STOP "UTest_distToFloorNode: Test 13 Failed."
END IF

IF( ANY(   distToFloorNode([10,10,10],[100._dp,100._dp,100._dp],[13._dp,13._dp,13._dp],.TRUE.) - [0.3_dp,0.3_dp,0.3_dp] &
> EPSILON(1.0_dp)*[100._dp,100._dp,100._dp] )) THEN
STOP "UTest_distToFloorNode: Test 14 Failed."
END IF

IF( ANY(   distToFloorNode([100,100,100],[10._dp,10._dp,10._dp],[13._dp,13._dp,13._dp],.TRUE.) - [z,z,z] &
> EPSILON(1.0_dp)*[10._dp,10._dp,10._dp] )) THEN
STOP "UTest_distToFloorNode: Test 15 Failed."
END IF

IF( ANY(   distToFloorNode([100,100,100],[10._dp,10._dp,10._dp],[13.1_dp,13.1_dp,13.1_dp],.TRUE.) - [z,z,z] &
> EPSILON(1.0_dp)*[10._dp,10._dp,10._dp] )) THEN
STOP "UTest_distToFloorNode: Test 16 Failed."
END IF

IF( ANY(   distToFloorNode([100,100,100],[10._dp,10._dp,10._dp],[13.31_dp,13.31_dp,13.31_dp],.TRUE.) - [o,o,o]/10._dp &
> EPSILON(1.0_dp)*[10._dp,10._dp,10._dp] )) THEN
STOP "UTest_distToFloorNode: Test 17 Failed."
END IF

IF( ANY(   distToFloorNode([100,100,100],[10._dp,10._dp,10._dp],[19.99_dp,19.99_dp,19.99_dp],.TRUE.) - [o,o,o]*0.9_dp &
> EPSILON(1.0_dp)*[10._dp,10._dp,10._dp] )) THEN
STOP "UTest_distToFloorNode: Test 18 Failed."
END IF

END SUBROUTINE UTest_distToFloorNode

!=================================================================================================================================
pure function factorial(n) ! computes the factorial of any integer n, i.e., n!
    implicit none
    integer, intent(in) :: n
    integer :: i, factorial
    select case (n)
    case (0)
        factorial = 1
    case default
        factorial = product([(i, i=1,n)])
    end select
end function factorial

!=================================================================================================================================
pure subroutine deduce_optimal_histogram_properties( n, maxrange, nbins, binwidth)
    implicit none
    integer, intent(in)   :: n ! total number of points to be histogramed
    real(dp), intent(in)  :: maxrange ! maximum range of the histogram (e.g., r max for g(r))
    integer, intent(out)  :: nbins ! number of bins
    real(dp), intent(out) :: binwidth ! width of a bin
    nbins    = ceiling( 2*real(n)**(1._dp/3._dp) ) ! Rice Rule, see http://en.wikipedia.org/wiki/Histogram
    nbins    = nbins/2 ! @Max sans lien avec la Rice Rule mais pragmatique. Mon nombre de bins est trop grand.
    binwidth = maxrange/real(nbins,dp) ! Width of each bin of the histogram
end subroutine deduce_optimal_histogram_properties

!=================================================================================================================================
subroutine spline(x,y,n,yp1,ypn,y2)
    implicit none
    integer, intent(in) :: n
    real(dp), intent(in) :: yp1,ypn,x(n),y(n)
    real(dp), intent(out) :: y2(n)
    ! Given arrays x(1:n) and y(1:n) containing a tabulated function, i.e., yi = f(xi), with
    ! x1 < x2 < ... < xN , and given values yp1 and ypn for the first derivative of the interpolating function
    ! at points 1 and n, respectively, this routine returns an array y2(1:n) of
    ! length n which contains the second derivatives of the interpolating function at the tabulated
    ! points xi. If yp1 and/or ypn are equal to 10^30 or larger, the routine is signaled to set
    ! the corresponding boundary condition for a natural spline, with zero second derivative on
    ! that boundary.
    ! see http://www.haoli.org/nr/bookfpdf/f3-3.pdf
    integer :: i,k
    real(dp) :: p, qn, sig, un, u(n)
    if(yp1>0.99e30 .or. yp1==huge(1._dp)) then
        y2(1)=0._dp
        u(1)=0._dp
    else
        y2(1)=-0.5
        u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
    end if
    ! this next loop induces IEEE underflow and denormal
    do i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2
        y2(i)=(sig-1)/p
        u(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
    end do
    if (ypn>0.99e30 .or. ypn==huge(1._dp)) then
        qn=0.
        un=0.
    else
        qn=0.5
        un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
    end if
    y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
    do k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
    end do
end subroutine spline

!=================================================================================================================================
pure subroutine splint(xa,ya,y2a,n,x,y)
    implicit none
    INTEGER, intent(in) :: n
    REAL(dp), intent(in) :: x,xa(n),y2a(n),ya(n)
    real(dp), intent(out) :: y
    !Given the arrays xa(1:n) and ya(1:n) of length n, which tabulate a function (with the
    !  xai ’s in order), and given the array y2a(1:n), which is the output from spline above,
    !  and given a value of x, this routine returns a cubic-spline interpolated value y
    INTEGER :: k,khi,klo
    REAL(dp) :: a,b,h
    klo=1
    khi=n
    do while (khi-klo>1)
        k=(khi+klo)/2
        if(xa(k).gt.x)then
            khi=k
        else
            klo=k
        end if
    end do
    h=xa(khi)-xa(klo)
    a=(xa(khi)-x)/h
    b=(x-xa(klo))/h
    y=a*ya(klo)+b*ya(khi)+ ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.
end subroutine splint

!-----------------------------------------------------------------------
! Cubic spline interpolation
!   [Reference:] Akima, H., 1970: J. ACM, 17, 589-602.
!-----------------------------------------------------------------------
pure subroutine akima_spline(ndim,x,y,n,x5,y5)
    implicit none
    integer,intent(in) :: ndim         ! number of grid points
    real(dp),intent(in) :: x(ndim) ! coordinate
    real(dp),intent(in) :: y(ndim) ! variable
    integer,intent(in) :: n            ! number of targets
    real(dp),intent(in) :: x5(n)   ! target coordinates
    real(dp),intent(out) :: y5(n)  ! target values
    integer :: i,j,m
    real(dp) :: dydx(5), ddydx(4), t(2), dx21, dx
    real(dp) :: wk

    tgt: do j=1,n
        do i=1,ndim
            if(x5(j) == x(i)) then
                y5(j) = y(i)
                cycle tgt
            end if
            if(x5(j) < x(i)) exit
        end do
        !       i-3   i-2   i-1    i    i+1   i+2
        !     ---+-----+-----+---*-+-----+-----+---
        !dydx       1     2     3     4     5
        !ddydx         1     2     3     4
        !t                   1     2
        if(i==2) then
            do m=3,5
                dydx(m) = (y(i-3+m)-y(i-4+m)) / (x(i-3+m)-x(i-4+m))
            end do
            dydx(2) = 2.0_dp*dydx(3) - dydx(4)
            dydx(1) = 2.0_dp*dydx(2) - dydx(3)
        else if(i==3) then
            do m=2,5
                dydx(m) = (y(i-3+m)-y(i-4+m)) / (x(i-3+m)-x(i-4+m))
            end do
            dydx(1) = 2.0_dp*dydx(2) - dydx(3)
        else if(i==ndim) then
            do m=1,3
                dydx(m) = (y(i-3+m)-y(i-4+m)) / (x(i-3+m)-x(i-4+m))
            end do
            dydx(4) = 2.0_dp*dydx(3) - dydx(2)
            dydx(5) = 2.0_dp*dydx(4) - dydx(3)
        else if(i==ndim-1) then
            do m=1,4
                dydx(m) = (y(i-3+m)-y(i-4+m)) / (x(i-3+m)-x(i-4+m))
            end do
            dydx(5) = 2.0_dp*dydx(4) - dydx(3)
        else
            do m=1,5
                dydx(m) = (y(i-3+m)-y(i-4+m)) / (x(i-3+m)-x(i-4+m))
            end do
        end if
        do m=1,4
            ddydx(m) = abs(dydx(m+1) - dydx(m))
        end do
        do m=1,2
            wk = ddydx(m+2) + ddydx(m)
            if(wk == 0) then
                t(m) = 0.0_dp
            else
                t(m) = (ddydx(m+2)*dydx(m+1)+ddydx(m)*dydx(m+2))/wk
            end if
        end do
        dx21 = x(i)-x(i-1)
        dx = x5(j) - x(i-1)
        y5(j) = y(i-1) &
        & + dx*t(1) &
        & + dx*dx*(3.0_dp*dydx(3)-2.0_dp*t(1)-t(2))/dx21 &
        & + dx*dx*dx*(t(1)+t(2)-2.0_dp*dydx(3))/dx21/dx21
    end do tgt

end subroutine akima_spline

END MODULE mathematica

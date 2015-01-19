!===================================================================================================================================
MODULE mathematica
!===================================================================================================================================
! This module implements several usefull functions of Mathematica
    USE precision_kinds     ,ONLY:dp,i2b
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: chop, TriLinearInterpolation, UTest_TrilinearInterpolation, floorNode, UTest_floorNode, ceilingNode, &
        UTest_ceilingNode, distToFloorNode, UTest_distToFloorNode, factorial, deduce_optimal_histogram_properties, &
        splint, spline

    CONTAINS

    !===============================================================================================================================
    PURE FUNCTION chop(x,delta)
    !===============================================================================================================================
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


    !===============================================================================================================================
    PURE FUNCTION TriLinearInterpolation (cube,x)
    !===============================================================================================================================
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


    !===============================================================================================================================
    SUBROUTINE UTest_TrilinearInterpolation
    !===============================================================================================================================
    ! Tests the pure function TriLinearInterpolation where result is known:
    ! - if the point is one of the corners
    ! - if it is on the center of the cube
    ! Then it tests that no answer is higher or lower than the maximum or minimum value of any corner.
        IMPLICIT NONE
        REAL(dp) :: A(0:1,0:1,0:1), x(1:3)
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


    !===============================================================================================================================
    PURE FUNCTION floorNode(gridnode,gridlen,x,pbc)
    !===============================================================================================================================
        IMPLICIT NONE
        INTEGER(i2b) :: floorNode(3)
        INTEGER(i2b), INTENT(IN) :: gridnode(3)
        REAL(dp), INTENT(IN) :: gridlen(3), x(3)
        LOGICAL, INTENT(IN) :: pbc ! periodic boundary counditions
        floorNode = FLOOR(MODULO(x,gridlen)/(gridlen/REAL(gridnode))) +1
    END FUNCTION floorNode
    !===============================================================================================================================


    !===============================================================================================================================
    SUBROUTINE UTest_floorNode
    !===============================================================================================================================
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


    !===============================================================================================================================
    PURE FUNCTION ceilingNode(gridnode,gridlen,x,pbc)
    !===============================================================================================================================
        IMPLICIT NONE
        INTEGER(i2b) :: ceilingNode(3)
        INTEGER(i2b), INTENT(IN) :: gridnode(3)
        REAL(dp), INTENT(IN) :: gridlen(3), x(3)
        LOGICAL, INTENT(IN) :: pbc ! periodic boundary counditions
        ceilingNode = MODULO( floorNode(gridnode,gridlen,x,pbc)  ,gridnode) +1
    END FUNCTION ceilingNode
    !===============================================================================================================================


    !===============================================================================================================================
    SUBROUTINE UTest_ceilingNode
    !===============================================================================================================================
        IMPLICIT NONE
        LOGICAL, SAVE :: alreadydone=.FALSE.
        REAL(dp), PARAMETER :: z=0._dp, o=1.0_dp
        INTEGER(i2b) :: gridnode(3)
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


    !===============================================================================================================================
    PURE FUNCTION distToFloorNode(gridnode,gridlen,x,pbc)
    !===============================================================================================================================
    ! Given a grid (number of nodes per direction and length in Angstroms per direction),
    ! returns the distance to floor node in  grid units, i.e., in dx.
    ! 0._dp <= distToFloorNode < 1._dp
        IMPLICIT NONE
        REAL(dp) :: distToFloorNode(3), xfloor(3)
        INTEGER(i2b), INTENT(IN) :: gridnode(3)
        REAL(dp), INTENT(IN) :: gridlen(3), x(3)
        LOGICAL, INTENT(IN) :: pbc ! periodic boundary counditions
        xfloor = ABS(  x/(gridlen/REAL(gridnode)) - FLOOR(x/(gridlen/REAL(gridnode)))  )
        distToFloorNode = xfloor
!~         distToFloorNode = MIN(&
!~                                 xfloor,&
!~                                 ABS(1._dp-xfloor)&
!~                             )
    END FUNCTION distToFloorNode
    !===============================================================================================================================


    !===============================================================================================================================
    SUBROUTINE UTest_distToFloorNode
    !===============================================================================================================================
        REAL(dp)     :: xfloor(3)
        LOGICAL      :: pbc ! periodic boundary counditions
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
    !===============================================================================================================================

  !=================================================================================================================================
  pure function factorial(n) ! computes the factorial of any integer n, i.e., n!
  !=================================================================================================================================
    use precision_kinds, only: i2b
    implicit none
    integer(i2b), intent(in) :: n
    integer(i2b) :: i, factorial
    select case (n)
    case (0)
      factorial = 1
    case default
      factorial = product([(i, i=1,n)])
    end select
  end function factorial
  !=================================================================================================================================

  !=================================================================================================================================
  pure subroutine deduce_optimal_histogram_properties( n, maxrange, nbins, binwidth)
  !=================================================================================================================================
    implicit none
    integer, intent(in)   :: n ! total number of points to be histogramed
    real(dp), intent(in)  :: maxrange ! maximum range of the histogram (e.g., r max for g(r))
    integer, intent(out)  :: nbins ! number of bins
    real(dp), intent(out) :: binwidth ! width of a bin
    nbins    = ceiling( 2*real(n)**(1._dp/3._dp) ) ! Rice Rule, see http://en.wikipedia.org/wiki/Histogram
    binwidth = maxrange/real(nbins,dp) ! Width of each bin of the histogram
  end subroutine deduce_optimal_histogram_properties
  !=================================================================================================================================

  !=================================================================================================================================
  pure subroutine spline(x,y,n,yp1,ypn,y2)
  !=================================================================================================================================
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
    if(yp1>0.99e30) then
      y2(1)=0._dp
      u(1)=0._dp
    else
      y2(1)=-0.5
      u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
    end if
    do i=2,n-1
      sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
      p=sig*y2(i-1)+2
      y2(i)=(sig-1)/p
      u(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
    end do
    if (ypn>0.99e30) then
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

  !=================================================================================================================================
  pure subroutine splint(xa,ya,y2a,n,x,y)
  !=================================================================================================================================
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
    1 if (khi-klo>1) then
      k=(khi+klo)/2
      if(xa(k).gt.x)then
        khi=k
      else
        klo=k
      end if
      goto 1
    end if
    h=xa(khi)-xa(klo)
    a=(xa(khi)-x)/h
    b=(x-xa(klo))/h
    y=a*ya(klo)+b*ya(khi)+ ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.
  end subroutine splint
  !=================================================================================================================================

END MODULE
!===================================================================================================================================

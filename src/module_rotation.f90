module module_rotation
    use iso_c_binding, only: dp=>c_double
    implicit none
    private
    public :: angle, thetaofq, phiofq, rotation_matrix_between_complex_spherical_harmonics_lu
    real(dp), private :: epsdp=epsilon(1._dp)
contains

    pure function thetaofq(qx,qy,qz)
        implicit none
        real(dp), intent(in) :: qx, qy, qz
        real(dp) :: thetaofq
        if (qx**2+qy**2+qz**2 > epsdp) then
            thetaofq=acos(qz/sqrt(qx**2+qy**2+qz**2))
        else
            thetaofq=0._dp
        end if
    end function thetaofq

    pure function phiofq(qx,qy,qz)
        implicit none
        real(dp), intent(in) :: qx, qy
        real(dp), intent(in), optional :: qz
        real(dp) :: phiofq
        if (qx**2+qy**2>epsdp) then
            phiofq=angle(qx,qy)
        else
          phiofq=0._dp
        end if
    end function phiofq


    pure function angle (x,y)
        !
        ! Returns the oriented angle between (0,x) and (x,y). This angle is called phi
        !
        implicit none
        real(dp), parameter :: twopi=2._dp*acos(-1._dp)
        real(dp), intent(in) :: x,y
        real(dp) :: angle
        real(dp) :: xx
        if (x==0._dp .and. y==0._dp) then
            angle = 0._dp
        else
            xx = acos( x/sqrt(x**2 + y**2) )
            if (y>=0._dp) then
                angle = xx
            else
                angle = twopi - xx
            end if
        end if
    end function angle

    pure function rotation_matrix_from_lab_to_q_frame_lu (q) result (rmat)
        implicit none
        real(dp), intent(in) :: q(3)
        real(dp) :: rmat(3,3), rmat1(3), rmat2(3), rmat3(3)
        if (q(1)==0._dp .and. q(2)==0._dp .and. q(3)==0._dp) then
            rmat3 = (/0._dp,0._dp,1._dp/)                            ! theta definied as zero.
        else
            rmat3 = q
        end if
        if (rmat3(1)/=0._dp .or. rmat3(2)/=0._dp) then       ! if rmat3 is along with axis z, the GSH is null, and we don't carre about phi.
            rmat2 = cross_product((/0._dp,0._dp,1._dp/),rmat3) ! in the MDFT definition of Omega, the rotation axes are z-y-z.
        else
            rmat2 = cross_product(rmat3,(/1._dp,0._dp,0._dp/)) ! cross product of rmat3 and axis x gives axis y, phi definied as zero.
        end if
        rmat3 = rmat3 / norm2(rmat3)
        rmat2 = rmat2 / norm2(rmat2)       ! to avoid round up error if rmat3 is so closed to z.
        rmat1 = cross_product(rmat2,rmat3)
        rmat(:,1) = rmat1
        rmat(:,2) = rmat2
        rmat(:,3) = rmat3
    end function rotation_matrix_from_lab_to_q_frame_lu


    pure function cross_product(a,b)
        !
        ! Cross product of two vectors
        !
        implicit none
        real(dp) :: cross_product(3)
        real(dp), intent(in) :: a(3), b(3)
        cross_product(1) = a(2)*b(3)-a(3)*b(2)
        cross_product(2) = a(3)*b(1)-a(1)*b(3)
        cross_product(3) = a(1)*b(2)-a(2)*b(1)
    end function cross_product

     function rotation_matrix_between_complex_spherical_harmonics_lu(q,mmax) result(R)

        use precision_kinds, only : i2b,dp
        implicit none

        integer(i2b), intent(in) :: mmax
        real(dp), parameter :: rac2=sqrt(2._dp)
        real(dp), dimension(3), intent(in) :: q
        complex(dp), dimension(0:mmax,-mmax:mmax,-mmax:mmax) :: R
        real(dp), dimension(3) :: rmat1,rmat2,rmat3
        real(dp), dimension(3,3) :: rmat
        real(dp), dimension(0:mmax,-mmax:mmax,-mmax:mmax) :: f,g
        integer(i2b) :: l,l1,m,m1,m1min,k
        real(dp), dimension(0:mmax,-mmax:mmax,-mmax:mmax) :: a,b,c,d
        real(dp), dimension(-1:2*mmax+1) :: rac
        real(dp), parameter :: zerodp=0._dp
        complex(dp), parameter :: zeroc=complex(0._dp,0._dp)

        a=zerodp
        b=zerodp
        c=zerodp
        d=zerodp
        f=zerodp
        g=zerodp
        R=zeroc

        ! prepare coefficients
        rac(-1)=0._dp
        do k=0,2*mmax+1
            rac(k)=sqrt(real(k,dp))
        end do

        if (mmax /= 0) then
            do l=1,mmax; do m=-l,l
                do m1=-l+1,l-1
                    a(l,m,m1) = rac(l+m)*rac(l-m)/(rac(l+m1)*rac(l-m1))
                    b(l,m,m1) = rac(l+m)*rac(l+m-1)/(rac2*rac(l+m1)*rac(l-m1))
                end do
                m1 = l
                c(l,m,m1) = rac2*rac(l+m)*rac(l-m)/(rac(l+m1)*rac(l+m1-1))
                d(l,m,m1) = rac(l+m)*rac(l+m-1)/(rac(l+m1)*rac(l+m1-1))
            end do; end do
        end if


        ! mmax = 0
        f(0,0,0) = 1._dp
        g(0,0,0) = 0._dp

        if (mmax == 0) then
            R = cmplx(f,g,dp)
            return
        end if

        ! build q-frame
        if (q(1)==0._dp .and. q(2)==0._dp .and. q(3)==0._dp) then
            rmat3 = (/0._dp,0._dp,1._dp/)                            ! theta definied as zero.
        else
            rmat3 = q
        end if
        if (rmat3(1)/=0._dp .or. rmat3(2)/=0._dp) then       ! if rmat3 is along with axis z, the GSH is null, and we don't carre about phi.
            rmat2 = cross_product((/0._dp,0._dp,1._dp/),rmat3) ! in the MDFT definition of Omega, the rotation axes are z-y-z.
        else ! qx=qy=0
            if (rmat3(3)>=0) then
                rmat2 = cross_product(rmat3,(/1._dp,0._dp,0._dp/)) ! cross product of rmat3 and axis x gives axis y, phi definied as zero.
            else
                rmat2 = -cross_product(rmat3,(/1._dp,0._dp,0._dp/))
            end if
        end if
        rmat3 = rmat3 / norm2(rmat3)
        rmat2 = rmat2 / norm2(rmat2)       ! to avoid round up error if rmat3 is so closed to z.
        rmat1 = cross_product(rmat2,rmat3)
        rmat(:,1) = rmat1
        rmat(:,2) = rmat2
        rmat(:,3) = rmat3

        ! mmax = 1
        f(1,-1,-1) = (rmat(2,2)+rmat(1,1))/2._dp
        f(1,-1, 0) = rmat(1,3)/rac2
        f(1,-1,+1) = (rmat(2,2)-rmat(1,1))/2._dp
        f(1, 0,-1) = rmat(3,1)/rac2
        f(1, 0, 0) = rmat(3,3)
        f(1, 0,+1) = -rmat(3,1)/rac2
        f(1,+1,-1) = (rmat(2,2)-rmat(1,1))/2._dp
        f(1,+1, 0) = -rmat(1,3)/rac2
        f(1,+1,+1) = (rmat(2,2)+rmat(1,1))/2._dp

        g(1,-1,-1) = (rmat(2,1)-rmat(1,2))/2._dp
        g(1,-1, 0) = rmat(2,3)/rac2
        g(1,-1,+1) = (-rmat(2,1)-rmat(1,2))/2._dp
        g(1, 0,-1) = -rmat(3,2)/rac2
        g(1, 0, 0) = 0._dp
        g(1, 0,+1) = -rmat(3,2)/rac2
        g(1,+1,-1) = (rmat(2,1)+rmat(1,2))/2._dp
        g(1,+1, 0) = rmat(2,3)/rac2
        g(1,+1,+1) = (rmat(1,2)-rmat(2,1))/2._dp

        if (mmax == 1) then
            R = cmplx(f,g,dp)
            return
        end if

        ! mmax > 1
        do l=2,mmax
            l1 = l-1
            do m=-l,l
                m1min = 0
                if (m>0) m1min = 1
                do m1=m1min,l-1
                    if( m==-l) then
                        f(l,m,m1)=b(l,-m,m1)*(f(1,-1,0)*f(l1,m+1,m1)-g(1,-1,0)*g(l1,m+1,m1))
                        g(l,m,m1)=b(l,-m,m1)*(f(1,-1,0)*g(l1,m+1,m1)+g(1,-1,0)*f(l1,m+1,m1))
                    else if( m==l) then
                        f(l,m,m1)=b(l,m,m1)*(f(1, 1,0)*f(l1,m-1,m1)-g(1, 1,0)*g(l1,m-1,m1))
                        g(l,m,m1)=b(l,m,m1)*(f(1, 1,0)*g(l1,m-1,m1)+g(1, 1,0)*f(l1,m-1,m1))
                    else
                        f(l,m,m1)=a(l, m,m1)*(f(1, 0,0)*f(l1,m  ,m1))+   &
                        b(l, m,m1)*(f(1, 1,0)*f(l1,m-1,m1)-g(1, 1,0)*g(l1,m-1,m1))+   &
                        b(l,-m,m1)*(f(1,-1,0)*f(l1,m+1,m1)-g(1,-1,0)*g(l1,m+1,m1))
                        g(l,m,m1)=a(l, m,m1)*(f(1, 0,0)*g(l1,m  ,m1))+   &
                        b(l, m,m1)*(f(1, 1,0)*g(l1,m-1,m1)+g(1, 1,0)*f(l1,m-1,m1))+   &
                        b(l,-m,m1)*(f(1,-1,0)*g(l1,m+1,m1)+g(1,-1,0)*f(l1,m+1,m1))
                    end if
                    f(l,-m,-m1)=(-1)**(m+m1)*f(l,m,m1)
                    g(l,-m,-m1)=-(-1)**(m+m1)*g(l,m,m1)

                    ! f(l,m,m1) = a(l,m,m1)*(f(1,0,0)*f(l1,m,m1)-g(1,0,0)*g(l1,m,m1))+         &
                    ! b(l,m,m1)*(f(1,+1,0)*f(l1,m-1,m1)-g(1,+1,0)*g(l1,m-1,m1))+   &
                    ! b(l,-m,m1)*(f(1,-1,0)*f(l1,m+1,m1)-g(1,-1,0)*g(l1,m+1,m1))
                    ! g(l,m,m1) = a(l,m,m1)*(f(1,0,0)*g(l1,m,m1)+g(1,0,0)*f(l1,m,m1))+         &
                    ! b(l,m,m1)*(f(1,+1,0)*g(l1,m-1,m1)+g(1,+1,0)*f(l1,m-1,m1))+   &
                    ! b(l,-m,m1)*(f(1,-1,0)*g(l1,m+1,m1)+g(1,-1,0)*f(l1,m+1,m1))
                    ! f(l,-m,-m1) = (-1)**(m+m1)*f(l,m,m1)
                    ! g(l,-m,-m1) = -(-1)**(m+m1)*g(l,m,m1)
                end do
                m1 = l
                if( m==-l) then
                    f(l,m,m1)=d(l,-m,m1)*(f(1,-1,+1)*f(l1,m+1,m1-1)-g(1,-1,+1)*g(l1,m+1,m1-1))
                    g(l,m,m1)=d(l,-m,m1)*(f(1,-1,+1)*g(l1,m+1,m1-1)+g(1,-1,+1)*f(l1,m+1,m1-1))
                else if( m==l) then
                    f(l,m,m1)=d(l,m,m1)*(f(1,+1,+1)*f(l1,m-1,m1-1)-g(1,+1,+1)*g(l1,m-1,m1-1))
                    g(l,m,m1)=d(l,m,m1)*(f(1,+1,+1)*g(l1,m-1,m1-1)+g(1,+1,+1)*f(l1,m-1,m1-1))
                else
                    f(l,m,m1)=c(l,m,m1)*(f(1,0,+1)*f(l1,m,m1-1)-g(1,0,+1)*g(l1,m,m1-1))+   &
                        d(l,m,m1)*(f(1,+1,+1)*f(l1,m-1,m1-1)-g(1,+1,+1)*g(l1,m-1,m1-1))+   &
                        d(l,-m,m1)*(f(1,-1,+1)*f(l1,m+1,m1-1)-g(1,-1,+1)*g(l1,m+1,m1-1))

                    g(l,m,m1)=c(l,m,m1)*(f(1,0,+1)*g(l1,m,m1-1)+g(1,0,+1)*f(l1,m,m1-1))+   &
                        d(l,m,m1)*(f(1,+1,+1)*g(l1,m-1,m1-1)+g(1,+1,+1)*f(l1,m-1,m1-1))+   &
                        d(l,-m,m1)*(f(1,-1,+1)*g(l1,m+1,m1-1)+g(1,-1,+1)*f(l1,m+1,m1-1))
                end if
                f(l,-m,-m1)=(-1)**(m+m1)*f(l,m,m1)
                g(l,-m,-m1)=-(-1)**(m+m1)*g(l,m,m1)

                ! f(l,m,m1) = c(l,m,m1)*(f(1,0,+1)*f(l1,m,m1-1)-g(1,0,+1)*g(l1,m,m1-1))+         &
                ! d(l,m,m1)*(f(1,+1,+1)*f(l1,m-1,m1-1)-g(1,+1,+1)*g(l1,m-1,m1-1))+   &
                ! d(l,-m,m1)*(f(1,-1,+1)*f(l1,m+1,m1-1)-g(1,-1,+1)*g(l1,m+1,m1-1))
                ! g(l,m,m1) = c(l,m,m1)*(f(1,0,+1)*g(l1,m,m1-1)+g(1,0,+1)*f(l1,m,m1-1))+         &
                ! d(l,m,m1)*(f(1,+1,+1)*g(l1,m-1,m1-1)+g(1,+1,+1)*f(l1,m-1,m1-1))+   &
                ! d(l,-m,m1)*(f(1,-1,+1)*g(l1,m+1,m1-1)+g(1,-1,+1)*f(l1,m+1,m1-1))
                ! f(l,-m,-m1) = (-1)**(m+m1)*f(l,m,m1)
                ! g(l,-m,-m1) = -(-1)**(m+m1)*g(l,m,m1)
            end do
        end do

        R = cmplx(f,g,dp)

    end function rotation_matrix_between_complex_spherical_harmonics_lu




end module module_rotation

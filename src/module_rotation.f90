module module_rotation
    use precision_kinds, only: dp
    implicit none
    private
    public :: angle, thetaofq, phiofq, rotation_matrix_between_complex_spherical_harmonics_lu, init, omega_prime_fonction_de_omega
    real(dp), private :: epsdp=epsilon(1._dp)
    integer, parameter, private :: mmax_max = 6 ! we will never have mmax > 5, actually we do for ACN
    real(dp), dimension(2:mmax_max,-mmax_max:mmax_max,-mmax_max:mmax_max), private :: a, b
    real(dp), dimension(2:mmax_max,-mmax_max:mmax_max), private :: c, d

contains

    subroutine init
        implicit none
        real(dp), parameter :: sqrt2=sqrt(2._dp)
        real(dp), parameter :: sqrtof(-1:2*mmax_max +1) = [0._dp, 0._dp, 1._dp, sqrt(2._dp), sqrt(3._dp), sqrt(4._dp), sqrt(5._dp),&
                                             sqrt(6._dp), sqrt(7._dp), sqrt(8._dp), sqrt(9._dp), sqrt(10._dp), sqrt(11._dp),sqrt(12._dp),sqrt(13._dp) ]
        integer :: l, m, m1
        a = 0._dp
        b = 0._dp
        c = 0._dp
        d = 0._dp
        do l=2,mmax_max
          do m=-l,l
            do m1=-l+1,l-1
                a(l,m,m1) = sqrtof(l+m)*sqrtof(l-m)/(sqrtof(l+m1)*sqrtof(l-m1))
                b(l,m,m1) = sqrtof(l+m)*sqrtof(l+m-1)/(sqrt2*sqrtof(l+m1)*sqrtof(l-m1))
            end do
            m1 = l
            c(l,m) = sqrt2*sqrtof(l+m)*sqrtof(l-m)/(sqrtof(l+m1)*sqrtof(l+m1-1))
            d(l,m) = sqrtof(l+m)*sqrtof(l+m-1)/(sqrtof(l+m1)*sqrtof(l+m1-1))
          end do
        end do
    end subroutine init

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

    pure function phiofq(qx,qy)
        implicit none
        real(dp), intent(in) :: qx, qy
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

    ! pure function rotation_matrix_from_lab_to_q_frame_lu (q) result (rmat)
    !     implicit none
    !     real(dp), intent(in) :: q(3)
    !     real(dp) :: rmat(3,3)
    !     real(dp) :: rmat1(3), rmat2(3), rmat3(3)
    !     real(dp) :: n2rmat3
    !     if (q(1)==0._dp .and. q(2)==0._dp .and. q(3)==0._dp) then
    !         rmat3 = (/0._dp,0._dp,1._dp/)                            ! theta definied as zero.
    !         n2rmat3 = 1._dp
    !     else
    !         rmat3 = q
    !         n2rmat3 = sqrt(q(1)**2+q(2)**2+q(3)**2)
    !     end if
    !     if (rmat3(1)/=0._dp .or. rmat3(2)/=0._dp) then       ! if rmat3 is along with axis z, the GSH is null, and we don't carre about phi.
    !         rmat2 = cross_product((/0._dp,0._dp,1._dp/),rmat3) ! in the MDFT definition of Omega, the rotation axes are z-y-z.
    !     else
    !         rmat2 = cross_product(rmat3,(/1._dp,0._dp,0._dp/)) ! cross product of rmat3 and axis x gives axis y, phi definied as zero.
    !     end if
    !     rmat3 = rmat3 / n2rmat3
    !     rmat2 = rmat2 / norm2(rmat2)       ! to avoid round up error if rmat3 is so closed to z.
    !     rmat1 = cross_product(rmat2,rmat3)
    !     rmat(:,1) = rmat1
    !     rmat(:,2) = rmat2
    !     rmat(:,3) = rmat3
    ! end function rotation_matrix_from_lab_to_q_frame_lu


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

    pure function rotation_matrix_between_complex_spherical_harmonics_lu (mmax, q) result(R)
        implicit none
        integer, intent(in) :: mmax
        real(dp), intent(in) :: q(3)
        complex(dp), dimension(0:mmax,-mmax:mmax,-mmax:mmax) :: R
        real(dp), dimension(3) :: rmat1,rmat2,rmat3
        real(dp), dimension(3,3) :: rmat
        real(dp), dimension(0:mmax,-mmax:mmax,-mmax:mmax) :: f, g
        integer :: l, l1, m, m1, m1min
        real(dp), parameter :: sqrt2 = sqrt(2._dp)

        if (mmax == 0) then
            R = (1._dp,0._dp)
            return
        end if

        ! R = cmplx(f,g)
        f = 0._dp
        g = 0._dp
        R = (0._dp,0._dp)

        ! m = 0
        f(0,0,0) = 1._dp
        g(0,0,0) = 0._dp

        !
        ! Build q-frame XYZ
        !

        !
        ! Start by aligning the new Z along q
        !
        if (q(1)==0._dp .and. q(2)==0._dp .and. q(3)==0._dp) then
            rmat3 = (/0._dp,0._dp,1._dp/)                            ! theta definied as zero.
        else
            rmat3 = q
        end if
        !
        ! Then deal with X and Y
        !
        if (rmat3(1)/=0._dp .or. rmat3(2)/=0._dp) then       ! if rmat3 is along with axis z, the GSH is null, and we don't carre about phi.
            rmat2 = cross_product((/0._dp,0._dp,1._dp/),rmat3) ! in the MDFT definition of Omega, the rotation axes are z-y-z.
        else ! qx=qy=0
            if (rmat3(3)>=0) then
                rmat2 = cross_product(rmat3,(/1._dp,0._dp,0._dp/)) ! cross product of rmat3 and axis x gives axis y, phi definied as zero.
            else
                rmat2 = -cross_product(rmat3,(/1._dp,0._dp,0._dp/))
            end if
        end if
        !
        ! Normalize
        !
        rmat3 = rmat3 / norm2(rmat3)
        rmat2 = rmat2 / norm2(rmat2)       ! to avoid round up error if rmat3 is so closed to z.
        rmat1 = cross_product(rmat2,rmat3)
        rmat(:,1) = rmat1
        rmat(:,2) = rmat2
        rmat(:,3) = rmat3

        ! m = 1
        f(1,-1,-1) = (rmat(2,2)+rmat(1,1))*0.5_dp
        f(1,-1, 0) = rmat(1,3)/sqrt2
        f(1,-1,+1) = (rmat(2,2)-rmat(1,1))*0.5_dp
        f(1, 0,-1) = rmat(3,1)/sqrt2
        f(1, 0, 0) = rmat(3,3)
        f(1, 0,+1) = -rmat(3,1)/sqrt2
        f(1,+1,-1) = (rmat(2,2)-rmat(1,1))*0.5_dp
        f(1,+1, 0) = -rmat(1,3)/sqrt2
        f(1,+1,+1) = (rmat(2,2)+rmat(1,1))*0.5_dp

        g(1,-1,-1) = (rmat(2,1)-rmat(1,2))*0.5_dp
        g(1,-1, 0) = rmat(2,3)/sqrt2
        g(1,-1,+1) = (-rmat(2,1)-rmat(1,2))*0.5_dp
        g(1, 0,-1) = -rmat(3,2)/sqrt2
        g(1, 0, 0) = 0._dp
        g(1, 0,+1) = -rmat(3,2)/sqrt2
        g(1,+1,-1) = (rmat(2,1)+rmat(1,2))*0.5_dp
        g(1,+1, 0) = rmat(2,3)/sqrt2
        g(1,+1,+1) = (rmat(1,2)-rmat(2,1))*0.5_dp

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
                        f(l,m,m1) = b(l,-m,m1)*(f(1,-1,0)*f(l1,m+1,m1)-g(1,-1,0)*g(l1,m+1,m1))
                        g(l,m,m1) = b(l,-m,m1)*(f(1,-1,0)*g(l1,m+1,m1)+g(1,-1,0)*f(l1,m+1,m1))
                    else if( m==l) then
                        f(l,m,m1) = b(l,m,m1)*(f(1, 1,0)*f(l1,m-1,m1)-g(1, 1,0)*g(l1,m-1,m1))
                        g(l,m,m1) = b(l,m,m1)*(f(1, 1,0)*g(l1,m-1,m1)+g(1, 1,0)*f(l1,m-1,m1))
                    else
                        f(l,m,m1) = a(l, m,m1)*(f(1, 0,0)*f(l1,m  ,m1))+   &
                                    b(l, m,m1)*(f(1, 1,0)*f(l1,m-1,m1)-g(1, 1,0)*g(l1,m-1,m1))+   &
                                    b(l,-m,m1)*(f(1,-1,0)*f(l1,m+1,m1)-g(1,-1,0)*g(l1,m+1,m1))
                        g(l,m,m1) = a(l, m,m1)*(f(1, 0,0)*g(l1,m  ,m1))+   &
                                    b(l, m,m1)*(f(1, 1,0)*g(l1,m-1,m1)+g(1, 1,0)*f(l1,m-1,m1))+   &
                                    b(l,-m,m1)*(f(1,-1,0)*g(l1,m+1,m1)+g(1,-1,0)*f(l1,m+1,m1))
                    end if
                    f(l,-m,-m1)=(-1)**(m+m1)*f(l,m,m1)
                    g(l,-m,-m1)=-(-1)**(m+m1)*g(l,m,m1)
                end do

                m1 = l
                if( m==-l) then
                    f(l,m,m1)=d(l,-m)*(f(1,-1,+1)*f(l1,m+1,m1-1)-g(1,-1,+1)*g(l1,m+1,m1-1))
                    g(l,m,m1)=d(l,-m)*(f(1,-1,+1)*g(l1,m+1,m1-1)+g(1,-1,+1)*f(l1,m+1,m1-1))
                else if( m==l) then
                    f(l,m,m1)=d(l,m)*(f(1,+1,+1)*f(l1,m-1,m1-1)-g(1,+1,+1)*g(l1,m-1,m1-1))
                    g(l,m,m1)=d(l,m)*(f(1,+1,+1)*g(l1,m-1,m1-1)+g(1,+1,+1)*f(l1,m-1,m1-1))
                else
                    f(l,m,m1)=c(l,m) *(f(1,0,+1)*f(l1,m,m1-1)-g(1,0,+1)*g(l1,m,m1-1))+   &
                              d(l,m) *(f(1,+1,+1)*f(l1,m-1,m1-1)-g(1,+1,+1)*g(l1,m-1,m1-1))+   &
                              d(l,-m)*(f(1,-1,+1)*f(l1,m+1,m1-1)-g(1,-1,+1)*g(l1,m+1,m1-1))

                    g(l,m,m1)=c(l,m)*(f(1,0,+1)*g(l1,m,m1-1)+g(1,0,+1)*f(l1,m,m1-1))+   &
                              d(l,m)*(f(1,+1,+1)*g(l1,m-1,m1-1)+g(1,+1,+1)*f(l1,m-1,m1-1))+   &
                              d(l,-m)*(f(1,-1,+1)*g(l1,m+1,m1-1)+g(1,-1,+1)*f(l1,m+1,m1-1))
                end if
                f(l,-m,-m1)=(-1)**(m+m1)*f(l,m,m1)
                g(l,-m,-m1)=-(-1)**(m+m1)*g(l,m,m1)
            end do
        end do
        R = cmplx(f,g,dp)
    end function rotation_matrix_between_complex_spherical_harmonics_lu


    pure function omega_prime_fonction_de_omega( q, omega) result (omegaPrime)
        use precision_kinds, only: dp
        implicit none
        real(dp), intent(in) :: q(3), omega(3)
        real(dp) :: omegaPrime(3)
        real(dp) :: theta, phi, psi
        real(dp) :: cosTheta, sinTheta, cosPsiPrimeMinPsi, sinPsiPrimeMinPsi, psiPrimeMinPsi
        real(dp) :: thetaPrime, phiPrime, psiPrime, cosThetaPrime, sinThetaPrime, cosPhiPrime, sinPhiPrime
        real(dp) :: thetaQ, phiQ, qx, qy, qz, cosThetaQ, sinThetaQ, cosPhiMinPhiQ, sinPhiMinPhiQ
        real(dp), parameter :: epsdp = epsilon(1._dp)

        theta = omega(1)
        phi = omega(2)
        psi = omega(3)

        cosTheta = cos(theta)
        sinTheta = sin(theta)

        thetaQ = thetaOfQ( q(1), q(2), q(3) )
        cosThetaQ = cos(thetaQ)
        sinThetaQ = sin(thetaQ)
        phiQ = phiOfQ( q(1), q(2) )
        cosPhiMinPhiQ = cos(phi-phiQ)
        sinPhiMinPhiQ = sin(phi-phiQ)

        cosThetaPrime = cosThetaQ * cosTheta + sinThetaQ * sinTheta * cosPhiMinPhiQ
        sinThetaPrime = sqrt( 1._dp - costhetaPrime**2 )
        thetaPrime = acos(cosThetaPrime)

        cosPhiPrime = (cosThetaQ*sinTheta*cosPhiMinPhiQ-sinThetaQ*cosTheta) / sinThetaPrime
        sinPhiPrime = sinTheta * sinPhiMinPhiQ / sinThetaPrime
        phiPrime = acos(cosPhiPrime)
        if(sinPhiPrime < 0) phiPrime = -phiPrime

        cosPsiPrimeMinPsi = (-sinThetaQ*cosTheta*cosPhiMinPhiQ+cosThetaQ*sinTheta)/sinThetaPrime
        sinPsiPrimeMinPsi = -sinThetaQ*sinPhiMinPhiQ/sinThetaPrime
        psiPrimeMinPsi = acos(cosPsiPrimeMinPsi)
        if(sinPsiPrimeMinPsi < 0) psiPrimeMinPsi = -psiPrimeMinPsi
        psiPrime = psiPrimeMinPsi + psi

        omegaPrime(1) = thetaPrime
        omegaPrime(2) = phiPrime
        omegaPrime(3) = psiPrime

    end function omega_prime_fonction_de_omega

end module module_rotation

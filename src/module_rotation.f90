module module_rotation
    use iso_c_binding, only: dp=>c_double
    implicit none
    private
    public :: rotation_matrix_from_lab_to_q_frame, harm_sph, rotation_matrix_between_complex_spherical_harmonics, angle,&
    thetaofq, phiofq, rotation_matrix_between_complex_spherical_harmonics_lu
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

    function rotation_matrix_from_lab_to_q_frame_lu (q) result (rmat)
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

    function rotation_matrix_from_lab_to_q_frame(q) result(R)
        !
        ! Given the coordinates q(3) in a frame (let's call it lab frame),
        ! this function produces the rotation matrix that
        ! transforms any coordinates in lab frame to some frame
        ! whose (0,0,1) vector is aligned with q.
        ! In other words, transform a vector in lab frame to a vector in intermolecular (or q) frame.
        !
        implicit none
        real(dp), intent(in) :: q(3)
        real(dp) :: R1(3), R2(3), R3(3), R(3,3)
        if (norm2(q)<=epsilon(1._dp)) then
            ! R(1,:)=[0,1,0]
            ! R(2,:)=[1,0,0]
            ! R(3,:)=[0,0,1] va pas sur psi
            R1=[1,0,0]
            R2=[0,1,0]
            R3=[0,0,1]
            R(1,:)=[1,0,0]
            R(2,:)=[0,1,0]
            R(3,:)=[0,0,1]
        else
            R3 = q/NORM2(q) ! The rotation matrix R will induce z to be aligned to vector q
            IF( ABS(R3(1)) < 0.99 ) THEN ! TODO COMPRENDRE POURQUOI JE PEUX PAS FIARE AVEC 010 au lieu de 100 PAR DEFAUT
                R1 = cross_product( REAL([1,0,0],dp) , R3 ) ! we dont care about R2 and R1, as long as they are orthogonal to R3.
            ELSE
                R1 = cross_product( REAL([0,1,0],dp) , R3 )
            END IF
            R1 = R1/NORM2(R1) ! normalization
            R2 = cross_product( R3, R1 )
            R(1,:) = R1
            R(2,:) = R2
            R(3,:) = R3
        end if
        ! block
        !     real(dp) :: theta,phi,psi
        !     real(dp) :: cos_theta,sin_theta,cos_phi,sin_phi,cos_psi,sin_psi
        !     theta=thetaofq(q(1),q(2),q(3))
        !     phi=phiofq(q(1),q(2))
        !     psi=0._dp
        !     cos_theta=cos(theta)
        !     sin_theta=sin(theta)
        !     cos_phi=cos(phi)
        !     sin_phi=sin(phi)
        !     cos_psi=cos(psi)
        !     sin_psi=sin(psi)
        !     R(1,1) = cos_theta*cos_phi*cos_psi-sin_phi*sin_psi
        !     R(1,2) = cos_theta*sin_phi*cos_psi+cos_phi*sin_psi
        !     R(1,3) =-sin_theta*cos_psi
        !     R(2,1) =-cos_theta*cos_phi*sin_psi-sin_phi*cos_psi
        !     R(2,2) =-cos_theta*sin_phi*sin_psi+cos_phi*cos_psi
        !     R(2,3) =sin_theta*sin_psi
        !     R(3,1) = sin_theta*cos_phi
        !     R(3,2) = sin_theta*sin_phi
        !     R(3,3) = cos_theta
        ! ! print*,"rotamat",norm2(R2-R(1,:))
        ! end block
    end function rotation_matrix_from_lab_to_q_frame

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

    function harm_sph( m, mu, mup, beta ) ! Luc's luc72p143
        !
        ! R^m_{mu,mup}(\beta)
        ! beta is the angle in radian
        !
        use mathematica, only: fact
        implicit none
        REAL(dp) :: harm_sph
        INTEGER, intent(in) :: m, mu, mup
        REAL(dp), intent(in) :: beta
        REAL(dp) :: beta0, x, cc, ss, pm, pm1, pm2
        INTEGER :: mu0, l, it
        if( abs(mu)>m .or. abs(mup)>m ) then
            harm_sph = 0
            return
        end if
        if( m==0 ) then
            harm_sph = 1
            return
        end if
        if( mu==0 .or. mup==0 ) then
            mu0 = mu
            beta0 = beta
            if(mu==0) then            ! je mets le 0 en second
                mu0 = mup
                beta0 = -beta
            end if
            x = 1                 ! si mu negatif, ca vaut (-1)**mu * valeur pour -mu
            if(mu0<0) then
                x = (-1)**mu0
                mu0 = -mu0
            end if
            cc = cos(beta0)           ! plutôt a partir d'une formule de recurrence stable
            pm1 = 0.               ! des polynomes de legendre associes pl,m
            pm = 1.                ! luc73p96
            do l = mu0+1,m
                pm2 = pm1
                pm1 = pm
                pm = (cc*REAL(2*l-1,dp)*pm1-REAL(l+mu0-1,dp)*pm2)/REAL(l-mu0,dp)
            end do
            harm_sph = x*(-1.)**mu0*sqrt(fact(m-mu0)/fact(m+mu0))&
            *fact(2*mu0)/(2.**mu0*fact(mu0))*sin(beta0)**mu0*pm
        else                        ! donc mu et mup non nuls, utiliser betement la formule de wigner
            harm_sph = 0
            cc = cos(0.5*beta)
            ss = sin(0.5*beta)
            do it = max(0,mu-mup),min(m+mu,m-mup)
                harm_sph = harm_sph+(-1.)**it/(fact(m+mu-it)*fact(m-mup-it)*fact(it)*&
                fact(it-mu+mup))*cc**(2*m+mu-mup-2*it)*ss**(2*it-mu+mup)
            end do
            harm_sph = sqrt(fact(m+mu)*fact(m-mu)*fact(m+mup)*fact(m-mup))*harm_sph
        end if
    end function harm_sph


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

        ! prepare coefficients
        rac(-1)=0._dp
        do k=0,2*mmax+1
            rac(k)=sqrt(real(k,dp))
        end do
        a = 0._dp; b = 0._dp; c = 0._dp; d = 0._dp
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
        else
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
                    f(l,m,m1)=c(l,m,m1)*(f(1,0,+1)*f(l1,m,m1-1)-g(1,0,+1)*g(l1,m,m1-1))+         &
                    d(l,m,m1)*(f(1,+1,+1)*f(l1,m-1,m1-1)-g(1,+1,+1)*g(l1,m-1,m1-1))+   &
                    d(l,-m,m1)*(f(1,-1,+1)*f(l1,m+1,m1-1)-g(1,-1,+1)*g(l1,m+1,m1-1))
                    g(l,m,m1)=c(l,m,m1)*(f(1,0,+1)*g(l1,m,m1-1)+g(1,0,+1)*f(l1,m,m1-1))+         &
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


    function rotation_matrix_between_complex_spherical_harmonics (q, lmax) result (R)
        !
        ! Given a column vector of real values q of dimension 3,
        ! this function returns the rotation matrix between complex spherical harmonics.
        ! The algorithm is from Choi et al., J. Chem. Phys. 111, 8825 (1999)
        ! http://dx.doi.org/10.1063/1.480229
        ! see eq 6.1 and following and 7.3 and following
        ! We use Choi's notation
        ! R^l_{mup,khi}(q)
        !
        ! Luc Belloni wrote the original fortran code (luc85p17)
        ! Maximilien Levesque rewrote it from scratch thanks to Luc's explanations
        ! Conditions for m==-l and m==l by ML
        ! All bugs due to ML
        !
        use precision_kinds, only: dp
        implicit none
        !
        integer, intent(in) :: lmax
        real(dp), intent(in) :: q(3)
        complex(dp) :: R(0:lmax,-lmax:lmax,-lmax:lmax)
        real(dp), dimension(0:lmax,-lmax:lmax,-lmax:lmax) :: fi, gi, a, b, c, d
        real(dp) :: rot(3,3)
        real(dp), parameter :: tsqrt2=sqrt(2._dp)
        real(dp) :: tsqrt(-1:2*lmax+1), adenom, cdenom, anumerator, bnumerator
        integer :: l,l1,m,m1,m1min
        !
        ! Init R, fi and gi to zero
        !
        fi = 0
        gi = 0
        R = complex(0._dp,0._dp)
        !
        ! l = 0
        !
        fi(0,0,0) = 1
        gi(0,0,0) = 0
        IF( lmax == 0 ) THEN
            R = complex(1._dp,0._dp)
            RETURN
        ELSE
            !
            ! We will need the rotation matrix from lab to q frame
            !
            Rot = rotation_matrix_from_lab_to_q_frame(q)
            if (norm2(rot)<=epsilon(1._dp))  error stop "oijouhskueghfuyfgjcfgr"
            where (rot<=epsilon(1._dp)) rot=0._dp
        END IF
        !
        ! l = 1
        !
        fi(1,-1,-1) = (Rot(2,2)+Rot(1,1))/2._dp
        fi(1,-1, 0) = Rot(1,3)/tsqrt2
        fi(1,-1, 1) = (Rot(2,2)-Rot(1,1))/2._dp
        fi(1, 0,-1) = Rot(3,1)/tsqrt2
        fi(1, 0, 0) = Rot(3,3)
        fi(1, 0, 1) = -Rot(3,1)/tsqrt2
        fi(1, 1,-1) = (Rot(2,2)-Rot(1,1))/2._dp
        fi(1, 1, 0) = -Rot(1,3)/tsqrt2
        fi(1, 1, 1) = (Rot(2,2)+Rot(1,1))/2._dp
        gi(1,-1,-1) = (Rot(2,1)-Rot(1,2))/2._dp
        gi(1,-1, 0) = Rot(2,3)/tsqrt2
        gi(1,-1, 1) = (-Rot(2,1)-Rot(1,2))/2._dp
        gi(1, 0,-1) = -Rot(3,2)/tsqrt2
        gi(1, 0, 0) = 0._dp
        gi(1, 0, 1) = -Rot(3,2)/tsqrt2
        gi(1, 1,-1) = (Rot(2,1)+Rot(1,2))/2._dp
        gi(1, 1, 0) = Rot(2,3)/tsqrt2
        gi(1, 1, 1) = (Rot(1,2)-Rot(2,1))/2._dp
        IF( lmax == 1 ) THEN
            R = cmplx(fi,gi,dp)
            RETURN
        ELSE
            !
            ! Tabulate sqrt
            !
            tsqrt(-1)=0._dp ! to avoid problems later
            tsqrt(0)=0._dp
            tsqrt(1)=1._dp
            tsqrt(2)=tsqrt2
            do l=3,2*lmax+1
                tsqrt(l) = SQRT(REAL(l,dp))
            end do
            !
            ! Coefficients needed for the recursion
            ! see eq 6.1 to 6.5 of Choi et al.
            !
            do l=1,lmax
                a(l,-l,:)=0._dp
                b(l,-l,:)=0._dp
                c(l,-l,:)=0._dp
                d(l,-l,:)=0._dp
                do m=-l+1,l
                    anumerator=tsqrt(l+m)*tsqrt(l-m)
                    bnumerator=tsqrt(l+m)*tsqrt(l+m-1)
                    do m1=-l+1,l-1
                        adenom = 1._dp/(tsqrt(l+m1)*tsqrt(l-m1))
                        a(l,m,m1)=anumerator*adenom
                        b(l,m,m1)=bnumerator*adenom/tsqrt2
                    end do
                    m1=l
                    cdenom=1/(tsqrt(l+m1)*tsqrt(l+m1-1))
                    c(l,m,m1)=tsqrt2*anumerator*cdenom
                    d(l,m,m1)=       bnumerator*cdenom
                end do
            end do
        END IF
        !
        ! l > 1    Use recursion
        !
        do l=2,lmax
            l1=l-1
            do m=-l,l
                m1min=0
                IF(m>0) m1min=1
                do m1=m1min,l-1
                    if( m==-l) then
                        fi(l,m,m1)=b(l,-m,m1)*(fi(1,-1,0)*fi(l1,m+1,m1)-gi(1,-1,0)*gi(l1,m+1,m1))
                        gi(l,m,m1)=b(l,-m,m1)*(fi(1,-1,0)*gi(l1,m+1,m1)+gi(1,-1,0)*fi(l1,m+1,m1))
                    else if( m==l) then
                        fi(l,m,m1)=b(l,m,m1)*(fi(1, 1,0)*fi(l1,m-1,m1)-gi(1, 1,0)*gi(l1,m-1,m1))
                        gi(l,m,m1)=b(l,m,m1)*(fi(1, 1,0)*gi(l1,m-1,m1)+gi(1, 1,0)*fi(l1,m-1,m1))
                    else
                        fi(l,m,m1)=a(l, m,m1)*(fi(1, 0,0)*fi(l1,m  ,m1))+   &
                        b(l, m,m1)*(fi(1, 1,0)*fi(l1,m-1,m1)-gi(1, 1,0)*gi(l1,m-1,m1))+   &
                        b(l,-m,m1)*(fi(1,-1,0)*fi(l1,m+1,m1)-gi(1,-1,0)*gi(l1,m+1,m1))
                        gi(l,m,m1)=a(l, m,m1)*(fi(1, 0,0)*gi(l1,m  ,m1))+   &
                        b(l, m,m1)*(fi(1, 1,0)*gi(l1,m-1,m1)+gi(1, 1,0)*fi(l1,m-1,m1))+   &
                        b(l,-m,m1)*(fi(1,-1,0)*gi(l1,m+1,m1)+gi(1,-1,0)*fi(l1,m+1,m1))
                    end if
                    fi(l,-m,-m1)=(-1)**(m+m1)*fi(l,m,m1)
                    gi(l,-m,-m1)=-(-1)**(m+m1)*gi(l,m,m1)
                end do
                m1=l
                if( m==-l) then
                    fi(l,m,m1)=d(l,-m,m1)*(fi(1,-1,+1)*fi(l1,m+1,m1-1)-gi(1,-1,+1)*gi(l1,m+1,m1-1))
                    gi(l,m,m1)=d(l,-m,m1)*(fi(1,-1,+1)*gi(l1,m+1,m1-1)+gi(1,-1,+1)*fi(l1,m+1,m1-1))
                else if( m==l) then
                    fi(l,m,m1)=d(l,m,m1)*(fi(1,+1,+1)*fi(l1,m-1,m1-1)-gi(1,+1,+1)*gi(l1,m-1,m1-1))
                    gi(l,m,m1)=d(l,m,m1)*(fi(1,+1,+1)*gi(l1,m-1,m1-1)+gi(1,+1,+1)*fi(l1,m-1,m1-1))
                else
                    fi(l,m,m1)=c(l,m,m1)*(fi(1,0,+1)*fi(l1,m,m1-1)-gi(1,0,+1)*gi(l1,m,m1-1))+         &
                    d(l,m,m1)*(fi(1,+1,+1)*fi(l1,m-1,m1-1)-gi(1,+1,+1)*gi(l1,m-1,m1-1))+   &
                    d(l,-m,m1)*(fi(1,-1,+1)*fi(l1,m+1,m1-1)-gi(1,-1,+1)*gi(l1,m+1,m1-1))
                    gi(l,m,m1)=c(l,m,m1)*(fi(1,0,+1)*gi(l1,m,m1-1)+gi(1,0,+1)*fi(l1,m,m1-1))+         &
                    d(l,m,m1)*(fi(1,+1,+1)*gi(l1,m-1,m1-1)+gi(1,+1,+1)*fi(l1,m-1,m1-1))+   &
                    d(l,-m,m1)*(fi(1,-1,+1)*gi(l1,m+1,m1-1)+gi(1,-1,+1)*fi(l1,m+1,m1-1))
                end if
                fi(l,-m,-m1)=(-1)**(m+m1)*fi(l,m,m1)
                gi(l,-m,-m1)=-(-1)**(m+m1)*gi(l,m,m1)
            end do
        end do
        !
        ! fi and gi are the REAL and imaginary part of our rotation matrix
        !
        R = cmplx(fi,gi,dp)
    end function rotation_matrix_between_complex_spherical_harmonics


end module module_rotation

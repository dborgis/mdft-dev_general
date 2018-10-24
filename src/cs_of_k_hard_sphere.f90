! This subroutine computes the direct correlation function of a hard sphere fluid
!it was added back in Oct 2018
subroutine cs_of_k_hard_sphere

   use precision_kinds ,only: i2b, dp
   use module_solvent, only: solvent
   use module_input           ,only: input_line, n_linesInFile, verbose,getinput
   use constants       ,only: fourpi, pi, zerodp
   use module_dcf             ,only:  c_s, c_s_hs
   use mathematica     ,only: spline, splint
   use hardspheres     ,only: hs, read_hard_sphere_radius
   use module_grid, only : grid ,dq

   implicit none

   real(dp) :: phi ! excess free energy density
   real(dp) :: n0, n1, n2, n3 ! weighted densities of homogeneous fluid
   real(dp) :: w0, w1, w2, w3 ! weight functions of homogeneous fluid
   real(dp) :: kR,R,x_loc,y_loc,k,coskR,sinkR
   real(dp), dimension (0:3,0:3) :: d2phi ! second partial derivative of phi wrt ni nj
   integer(i2b) :: i,j,nb_k, ios
   logical :: PY,CS ! the one which is true is the right equation of state
   character(2)  :: hs_functional

   if( allocated(c_s_hs%x) ) then
     return
     stop "cs_of_k_hard_sphere.f90 strange c_s_hs is already existing..."
   end if

   ! type HS should be allocated. It contains radius etc.
   if( .not. allocated(hs) ) then
     call read_hard_sphere_radius
   end if
  
   !I juste compute cs for 1 specie of solvent for now
   if( solvent(1)%nspec/=1 ) then
     stop "cs_of_k_hard_sphere works only for 1 solvent species"
   end if

   !This is quite ugly, but it is a inheritance and Ido not want to spend time changing it
   hs_functional=getinput%char('hs_functional')
   if (hs_functional=="CS") then
     CS = .true.
     PY = .false.
   else
     CS = .false.
     PY = .true.
   end if

   if( all([CS,PY].eqv..false.) ) then
     stop "in cs_of_k_hard_sphere I could not detect if you want PY or CS functions"
   end if

   ! Here we could generate as many points as wanted. In order to be coherent with the excess functional
   ! we will use as much nk than qmaxnecessary in module_energy_cproj_mrso
   nb_k = int(norm2([maxval(grid%kx(1:grid%nx)), maxval(grid%ky(1:grid%ny)), maxval(grid%kz(1:grid%nz/2+1))])/dq)+1
   allocate( c_s_hs%x(nb_k),  source=zerodp) ! k
   allocate( c_s_hs%y(nb_k),  source=zerodp) ! c(k) hard sphere
   allocate( c_s_hs%y2(nb_k), source=zerodp) ! dc(k)/dk
   R=hs(1)%R
   do i = 1, nb_k
     k=(i-1)*dq
     !print*, k
     c_s_hs%x(i) = k
     !R = hs(1)%R
     kR = k*R
     sinkR = sin(kR)
     coskR = cos(kR)
     if (abs(k)<=epsilon(1._dp)) then ! k=0
       w0 = 1.
       w1 = R
       w2 = FourPi*R**2
       w3 = FourPi/3.*R**3

       ! weighted densities of homogeneous fluid
       n0 = solvent(1)%n0 * w0
       n1 = solvent(1)%n0 * w1
       n2 = solvent(1)%n0 * w2
       n3 = solvent(1)%n0 * w3

     ELSE ! k > 0
       w0 = coskR + kR*sinkR/2.
       w1 = (sinkR + kR*coskR)/(2.*k)
       w2 = (FourPi*R*sinkR)/k
       w3 = FourPi*(sinkR-kR*coskR)/k**3
       ! n0 etc that are used below are thoses for the k=0
     END IF

     if( abs(n3-1)<=epsilon(1.0_dp) ) then
       stop "n3-1=0 in cs_of_k_hard_sphere.f90:69"
     end if
     ! expression of phi_exc depends obviously of the choice of the excess functional
     if ( PY ) then
       ! second partial derivative of phi wrt n_\alpha n_\beta
       d2phi(0,0) = 0
       d2phi(0,1) = 0
       d2phi(0,2) = 0
       d2phi(0,3) = 1./(1.-n3)

       d2phi(1,0) = 0
       d2phi(1,1) = 0
       d2phi(1,2) = 1./(1.-n3)
       d2phi(1,3) = n2/(1.-n3)**2

       d2phi(2,0) = 0
       d2phi(2,1) = d2phi(1,2)
       d2phi(2,2) = n2/(4*pi*(1-n3)**2)
       d2phi(2,3) = n1/(1.-n3)**2 +n2**2/(4*pi*(1-n3)**3)

       d2phi(3,0) = d2phi(0,3)
       d2phi(3,1) = d2phi(1,3)
       d2phi(3,2) = d2phi(2,3)
       d2phi(3,3) = (2.*n1*n2)/(1.-n3)**3 + n0/(1.-n3)**2 + n2**3/(4.*Pi*(1.-n3)**4)
     ELSE IF ( CS ) then ! CS

       ! second partial derivative of phi wrt n_\alpha n_\beta
       d2phi(0,0) = 0.d0
       d2phi(0,1) = 0.d0
       d2phi(0,2) = 0.d0
       d2phi(0,3) = 1.d0/(1-n3)

       d2phi(1,0) = d2phi(0,1)
       d2phi(1,1) = 0.d0
       d2phi(1,2) = 1.d0/(1.d0 - n3)
       d2phi(1,3) = n2/(1.d0 - n3)**2

       d2phi(2,0) = d2phi(0,2)
       d2phi(2,1) = d2phi(1,2)
       d2phi(2,2) = n2/(6.d0*(1.d0 - n3)**2*n3*Pi) + (n2*Log(1.d0 - n3))/(6.d0*n3**2*Pi) ! simplified by Mathematica
       d2phi(2,3) = (n3*(n2**2*(2.d0 - 5.d0*n3 + n3**2) + 12.d0*n1*(-1.d0 + n3)*n3**2*Pi) &
                    - 2.d0*n2**2*(-1.d0 + n3)**3*Log(1-n3))/(12.d0*(-1.d0 + n3)**3*n3**3*Pi)  ! simplified by Mathematica

       d2phi(3,0) = d2phi(0,3)
       d2phi(3,1) = d2phi(1,3)
       d2phi(3,2) = d2phi(2,3)
       d2phi(3,3) = (n3*(n2**3*(6.d0 - 21.d0*n3 + 26.d0*n3**2 - 5.d0*n3**3) - 72.d0*n1*n2*(-1.d0 + n3)*n3**3*Pi &
                  + 36.d0*n0*(-1.d0 + n3)**2*n3**3*Pi) + 6.d0*n2**3*(-1.d0 + n3)**4*Log(1.d0 - n3))/(36.d0*(-1.d0 + n3)**4*n3**4*Pi) ! simplified by Mathematica
     END IF ! only PY or CS for now

     ! direct correlation function
     c_s_hs%y(i) =    -( d2phi(0,0)*w0*w0 &
                       + d2phi(0,1)*w0*w1 &
                       + d2phi(0,2)*w0*w2 &
                       + d2phi(0,3)*w0*w3 &
                       + d2phi(1,0)*w1*w0 &
                       + d2phi(1,1)*w1*w1 &
                       + d2phi(1,2)*w1*w2 &
                       + d2phi(1,3)*w1*w3 &
                       + d2phi(2,0)*w2*w0 &
                       + d2phi(2,1)*w2*w1 &
                       + d2phi(2,2)*w2*w2 &
                       + d2phi(2,3)*w2*w3 &
                       + d2phi(3,0)*w3*w0 &
                       + d2phi(3,1)*w3*w1 &
                       + d2phi(3,2)*w3*w2 &
                       + d2phi(3,3)*w3*w3 )
   end do

   open(99,file="output/cs_hs_tabulated.dat")
     do i=1,size(c_s_hs%x)
       write(99,*) c_s_hs%x(i), c_s_hs%y(i)
     end do
   close(99)

   ! Prepare the future use of splines by calculating the second order derivative (y2) at each point
   !call spline( x=c_s_hs%x, y=c_s_hs%y, n=size(c_s_hs%x), yp1=huge(1._dp), ypn=huge(1._dp), y2=c_s_hs%y2)

end subroutine cs_of_k_hard_sphere

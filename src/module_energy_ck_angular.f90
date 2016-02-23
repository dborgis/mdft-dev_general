module module_energy_ck_angular

  use precision_kinds, only: dp
  use iso_c_binding

  implicit none
  private
  !
  ! FFTW3 header - modern (fortran 2003) version. Expects iso_c_binding
  !
  include 'fftw3.f03'
  type(c_ptr) :: plan_forward
  type(c_ptr) :: plan_backward
  type ck_type
    logical :: is_initiated = .false.
    real(dp) :: dk
    integer :: nk, npsi, nphi, ncostheta, molrotsymorder
    real(dp) :: dphi, dpsi, dcostheta
  end type ck_type
  type(ck_type) :: ck
  type angleind_type
      integer :: costheta, psi
      real(dp) :: phi
  end type angleind_type
  type ang_type
    real(dp), allocatable :: phi(:,:,:,:)
    integer, allocatable :: psi(:,:,:,:), costheta(:,:,:,:)
  end type ang_type
  type(ang_type) :: ang
  type(angleind_type), allocatable, dimension(:,:,:,:) :: angleInd ! integer table of correspondence omega(k,Omega)
  complex(dp), allocatable :: ck_angular(:,:,:,:,:,:)
  complex(dp), allocatable :: delta_rho_k(:,:,:,:) ! delta_rho(omega,k)
  real(dp), allocatable :: in_forward(:,:,:), out_backward(:,:,:)
  complex(dp), allocatable :: out_forward(:,:,:), in_backward(:,:,:)
  real(dp), allocatable :: y(:,:,:,:) ! gamma(omega,r)
  public :: energy_ck_angular

contains




subroutine energy_ck_angular (ff, df)
  use omp_lib
  use precision_kinds, only: dp
  use module_grid, only: grid
  use module_thermo, only: thermo
  use module_solvent, only: solvent
  implicit none
  real(dp), intent(out) :: ff
  real(dp), intent(inout) :: df(:,:,:,:,:) ! x y z o s
  integer :: nx, ny, nz, no, mrso, io, io1, io2, ix, iy, iz, ik, iphi12,ipsi1,ipsi2,icostheta1,icostheta2
  real(dp), parameter :: zerodp=0._dp, twopi=2._dp*acos(-1._dp)

  real(dp) :: rho0, dv, drho, kT, phi12, xi, vexc, phi1
  complex(dp) :: my_ck_angular
  real :: t(10)

  nx=grid%nx
  ny=grid%ny
  nz=grid%nz
  no=grid%no
  mrso=grid%molrotsymorder
  ff=zerodp
  rho0=solvent(1)%rho0
  dv=grid%dv
  kT=thermo%kbT

  call cpu_time(t(1))
  if (ck%is_initiated .eqv. .false.) then
    call read_ck_angular
    allocate( in_forward(nx,ny,nz), source=0._dp)
    allocate( out_forward(nx/2+1,ny,nz), source=complex(0._dp,0._dp))
    allocate( in_backward(nx/2+1,ny,nz), source=complex(0._dp,0._dp))
    allocate( out_backward(nx,ny,nz), source=0._dp)
    call dfftw_plan_dft_r2c_3d(  plan_forward,  nx,ny,nz, in_forward , out_forward , FFTW_MEASURE )
    call dfftw_plan_dft_c2r_3d(  plan_backward, nx,ny,nz, in_backward, out_backward, FFTW_MEASURE )
    allocate (delta_rho_k(no,nx/2+1,ny,nz), source=complex(0._dp,0._dp))
    allocate ( y(no,nx,ny,nz), source=0._dp)
    call gen_angleind
  end if
  call cpu_time(t(2))

  !
  ! Fourier transform delta rho => delta_rho_k
  !
  do io=1,no
    in_forward = rho0*(solvent(1)%xi(io,:,:,:)**2 -1._dp) ! delta_rho
    call dfftw_execute_dft_r2c( plan_forward, in_forward, out_forward )
    delta_rho_k(io,:,:,:) = out_forward*grid%w(io) ! NOTE THAT WE INCLUDE THE ANGULAR WEIGHT HERE
  end do
  call cpu_time(t(3))

  !
  ! Calcul de gamma == y
  !
  call cpu_time(t(4))


  do io1=1,no

    in_backward = complex(0._dp,0._dp)
!$OMP PARALLEL DO DEFAULT(FIRSTPRIVATE) SHARED(ck_angular,ang,grid,delta_rho_k,ck), REDUCTION(+:in_backward)
    do iz=1,nz
      do iy=1,ny
        do ix=1,nx/2+1
          ik = int( sqrt(grid%kx(ix)**2+grid%ky(iy)**2+grid%kz(iz)**2) /ck%dk +0.5) +1

          ipsi1 = ang%psi(io1,ix,iy,iz)
          icostheta1 = ang%costheta(io1,ix,iy,iz)
          phi1 = ang%phi(io1,ix,iy,iz)

          do io2=1,no
            phi12 = modulo( phi1 - ang%phi(io2,ix,iy,iz) ,twopi )! phi1 - phi2
            iphi12 = mod( int(phi12*ck%nphi/twopi), ck%nphi ) +1
            ipsi2 = ang%psi(io2,ix,iy,iz)
            icostheta2 = ang%costheta(io2,ix,iy,iz)
            my_ck_angular = ck_angular(ipsi2,iphi12,icostheta2,ipsi1,icostheta1,ik)
            in_backward(ix,iy,iz) = in_backward(ix,iy,iz) + my_ck_angular*delta_rho_k(io2,ix,iy,iz)
          end do

        end do
      end do
    end do
!$OMP END PARALLEL DO

    call dfftw_execute_dft_c2r ( plan_backward , in_backward, out_backward )
    y(io1,:,:,:) = out_backward/real(nx*ny*nz,dp)
  end do

  call cpu_time(t(5))


  ! excess free energy
  ff = 0._dp
  do iz=1,nz
    do iy=1,ny
      do ix=1,nx
        do io=1,no
          xi = solvent(1)%xi(io,ix,iy,iz)
          vexc = -kT*grid%w(io)*y(io,ix,iy,iz)
          ff = ff + vexc*rho0*(xi**2-1._dp) /2._dp *dv
          df(io,ix,iy,iz,1) = df(io,ix,iy,iz,1) +vexc*2._dp*rho0*xi
        end do
      end do
    end do
  end do
  call cpu_time(t(6))

  print*,
  print*,
  print*, "init =>", t(2)-t(1)
  print*, "fft deltarho", t(3)-t(2)
  print*, "coeur =", t(5)-t(4)
  print*, "calcul ff df", t(6)-t(5)
  print*,
  print*,
end subroutine energy_ck_angular





subroutine gen_angleind
  use module_grid, only: grid
  use module_rotation, only: angle
  implicit none
  integer :: ix,iy,iz,io
  real(dp), parameter :: twopi=2._dp*acos(-1._dp)
  real(dp) :: costheta_k, sintheta_k, cosphi_k, sinphi_k, phi_k,w_kz,w_kx,w_ky,u_kz,v_kz,cos_value,phi_value,psi_value,ksq

  allocate ( angleInd(grid%no,grid%nx/2+1,grid%ny,grid%nz) )
  angleInd%costheta = 0
  angleInd%phi      = 0._dp
  angleInd%psi      = 0

  do concurrent (ix=1:grid%nx/2+1, iy=1:grid%ny, iz=1:grid%nz)

    ksq = grid%kx(ix)**2 +grid%ky(iy)**2 +grid%kz(iz)**2

    if ( ksq > epsilon(1._dp) ) then
      costheta_k = grid%kz(iz)/sqrt(ksq)
      sintheta_k = sqrt(1._dp - costheta_k**2)
    else
      costheta_k = 1.0_dp
      sintheta_k = 0._dp
    end if

    if ( grid%kx(ix)**2 + grid%ky(iy)**2 > epsilon(1._dp) ) then
      phi_k = angle( grid%kx(ix) , grid%ky(iy) )
      cosphi_k = cos(phi_k)
      sinphi_k = sin(phi_k)
    else
      cosphi_k = 1._dp
      sinphi_k = 0._dp
    end if

    do io=1,grid%no
      w_kz =   grid%rotxz(io)*sintheta_k*cosphi_k &
             + grid%rotyz(io)*sintheta_k*sinphi_k &
             + grid%rotzz(io)*costheta_k
      w_kx =   grid%rotxz(io)*costheta_k*cosphi_k &
             + grid%rotyz(io)*costheta_k*sinphi_k &
             - grid%rotzz(io)*sintheta_k
      w_ky = - grid%rotxz(io)*sinphi_k &
             + grid%rotyz(io)*cosphi_k
      u_kz =   grid%rotxx(io)*sintheta_k*cosphi_k &
             + grid%rotyx(io)*sintheta_k*sinphi_k &
             + grid%rotzx(io)*costheta_k
      !   u_kx =   grid%rotxx(io)*costheta_k*cosphi_k &
      !          + grid%rotyx(io)*costheta_k*sinphi_k &
      !          - grid%rotzx(io)*sintheta_k
      !   u_ky = - grid%rotxx(io)*sinphi_k &
      !          + grid%rotyx(io)*cosphi_k
      v_kz =   grid%rotxy(io)*sintheta_k*cosphi_k &
             + grid%rotyy(io)*sintheta_k*sinphi_k &
             + grid%rotzy(io)*costheta_k
      !   v_kx =   grid%rotxy(io)*costheta_k*cosphi_k &
      !          + grid%rotyy(io)*costheta_k*sinphi_k &
      !          - grid%rotzy(io)*sintheta_k
      !   v_ky = - grid%rotxy(io)*sinphi_k &
      !          + grid%rotyy(io)*cosphi_k

      ! Calculate angles omega
      cos_value = w_kz
      phi_value = angle( w_kx , w_ky )
      psi_value = MODULO(angle(-u_kz,v_kz), twopi/grid%molRotSymOrder)

      ! Real omega values for interpolation
      ! IF (karim) THEN
      !   angleVal(l,m,n,o,p)%costheta = cos_value
      !   angleVal(l,m,n,o,p)%phi      = phi_value
      !   angleVal(l,m,n,o,p)%psi      = psi_value
      !
      !   ! Index omega(k,Omega)
      ! ELSE
      angleInd(io,ix,iy,iz)%costheta = MIN(INT((1._dp + cos_value)*ck%ncostheta/2._dp) + 1, ck%ncostheta)
      angleInd(io,ix,iy,iz)%phi      = phi_value
      angleInd(io,ix,iy,iz)%psi      = MOD(INT(psi_value*ck%npsi*grid%molRotSymOrder/twopi), ck%npsi) + 1
      ! END IF
    end do ! io
  end do ! concurrent ix,iy,iz


  ! Check final nomega
  if ( any(angleind%costheta <= 0._dp)) error stop "angleind%costheta is somewhere negative"
  if ( any(angleind%psi <= 0._dp)) error stop "angleind%psi is somewhere negative"
  if ( any(angleind%costheta > ck%ncostheta)) error stop "angleind%costheta is somewhere > ck%ncostheta"
  if ( any(angleInd%psi > ck%npsi)) error stop "angleind%psi is somewhere > ck%npsi"


  ! move in new arrays
  allocate (ang%costheta(grid%no,grid%nx/2+1,grid%ny,grid%nz), source=angleind%costheta)
  allocate (ang%phi(grid%no,grid%nx/2+1,grid%ny,grid%nz), source=angleind%phi)
  allocate (ang%psi(grid%no,grid%nx/2+1,grid%ny,grid%nz), source=angleind%psi)
  deallocate (angleind)

end subroutine gen_angleind








subroutine read_ck_angular
  use module_grid, only: grid
  implicit none
  logical :: ok
  integer, parameter :: myunit=33
  real(dp), parameter :: twopi=2._dp*acos(-1._dp)
  character(len=len("./input/ck_angular.in")), parameter :: filename="./input/ck_angular.in"

  if (ck%is_initiated .eqv. .true.) then
    print*, "error in read_ck_angular in module_energy_ck_angular"
    print*, "ck seems already initiated"
    error stop
  end if

  inquire(file=filename, exist=ok)
  if (.not.ok) error stop "fatal error. File input/ck_angular.in not found"

  open(myunit, file=filename, FORM='unformatted', action="read")
  read(myunit) ck%nk, ck%dk, ck%ncostheta, ck%nphi, ck%npsi, ck%molrotsymorder  ! Note that psi is from 0 to pi for water, while no other symetry is taken into account
  print*,
  print*,"==== ck_angular ==="
  print*,"ck%nk             =",ck%nk
  print*,"ck%dk             =",ck%dk
  print*,"ck%ncostheta      =",ck%ncostheta
  print*,"ck%nphi           =",ck%nphi
  print*,"ck%npsi           =",ck%npsi
  print*,"ck%molrotsymorder =",ck%molrotsymorder
  print*,"==== ck_angular ==="
  print*,

  if (ck%molrotsymorder /= grid%molrotsymorder) then
    print*, "llkjiuklhfsdkhukhfswilsfewjl le fichier ck_angular.in n'a pas molrotsymorder=2"
    print*, "ck%molrotsymorder de ck_angular.in =", ck%molrotsymorder
    print*, "while grid%molrotsymorder =", grid%molrotsymorder
    error stop "987IUHKJN" ! unique tag to simplify bash greps to buggy source
  end if

  allocate( ck_angular(ck%npsi,ck%npsi,ck%nphi,ck%ncostheta,ck%ncostheta,ck%nk)  ,source=complex(0._dp,0._dp) )

  rewind(myunit)
  ! what follows is here just because ck_angular.in is a single precision binary. One can't read directly the complex(dp)
  block
    complex :: ck_angular_sp(ck%npsi,ck%npsi,ck%nphi,ck%ncostheta,ck%ncostheta,ck%nk)
    read(myunit) ck%nk, ck%dk, ck%ncostheta, ck%nphi, ck%npsi, ck%molrotsymorder, ck_angular_sp
    ck_angular = ck_angular_sp
  end block


  block
    integer :: ipsi1,ipsi2,icostheta1,icostheta2,iphi12,ik
    complex(dp) :: ck_tmp(ck%npsi,ck%nphi,ck%ncostheta,ck%npsi,ck%ncostheta,ck%nk)
    do concurrent(ipsi1=1:ck%npsi, ipsi2=1:ck%npsi, iphi12=1:ck%nphi, icostheta1=1:ck%ncostheta, icostheta2=1:ck%ncostheta,&
       ik=1:ck%nk)
      ck_tmp(ipsi2,iphi12,icostheta2,ipsi1,icostheta1,ik) = ck_angular(ipsi1,ipsi2,iphi12,icostheta1,icostheta2,ik)
    end do
    deallocate (ck_angular)
    allocate (ck_angular(ck%npsi,ck%nphi,ck%ncostheta,ck%npsi,ck%ncostheta,ck%nk), source=ck_tmp)
  end block


  ck%dphi=twopi/ck%nphi
  ck%dpsi=twopi/(ck%npsi*ck%molrotsymorder)
  ck%dcostheta=2._dp/ck%ncostheta

  ck%is_initiated = .true.

end subroutine



end module module_energy_ck_angular

SUBROUTINE energy_threebody_faster
USE precision_kinds, only:i2b, dp
use input,only : input_line , input_log
use constants, only: twopi
use quadrature, only : angGrid, molRotGrid
use system, only: nfft1 , nfft2 , nfft3 , deltaV , rho_0 , sig_mol , sig_solv , Lx , Ly , Lz ,&
&    id_mol, x_mol , y_mol , z_mol , kbT , nb_species, nb_solute_sites, deltax, deltay, deltaz&
& , lambda1_mol , lambda2_mol, deltaV,n_0
USE minimizer, ONLY:cg_vect,dF,FF
use fft, only : fftw3

IMPLICIT NONE
real(dp), parameter :: rmin1 = 1.5_dp, rsw1 = 2.0_dp, rmin2 = 2.25_dp, rsw2 = 2.5_dp, rmax2 = 5.0_dp, d_w = 1.9_dp
integer(i2b)::icg
integer(i2b) :: i,j,k,o,p,n, i1, j1, k1
real(dp) :: time0,time1
real(dp), allocatable, dimension(:,:,:)::rho
complex(dp), allocatable, dimension(:,:,:)::rho_k
real(dp) ::rk2,xk2,yk2,zk2,r,x,y,z
real(dp) ::  DHxx_ijk,DHyy_ijk,DHzz_ijk,DHxy_ijk,DHxz_ijk,DHyz_ijk,DH0_ijk,DHx_ijk,DHy_ijk,DHz_ijk
real(dp), dimension(nfft1,nfft2,nfft3) :: Gxx,Gyy,Gzz,Gxy,Gxz,Gyz,Gx,Gy,Gz,G0 
real(dp), dimension(nfft1,nfft2,nfft3,nb_solute_sites) :: DHxx,DHyy,DHzz,DHxy,DHxz,DHyz,DH0,DHx,DHy,DHz 
real(dp), dimension(nb_solute_sites) :: Hxx,Hyy,Hzz,Hxy,Hxz,Hyz,Hx,Hy,Hz,H0
real(dp), dimension(nfft1,nfft2,nfft3) :: Fxx,Fyy,Fzz,Fxy,Fxz,Fyz,Fx,Fy,Fz,F0
complex(dp), dimension(nfft1/2+1,nfft2,nfft3) :: Fxx_k,Fyy_k,Fzz_k,Fxy_k,Fxz_k,Fyz_k,Fx_k,Fy_k,Fz_k,F0_k
real(dp), dimension(nfft1,nfft2,nfft3) :: Axx,Ayy,Azz,Axy,Axz,Ayz,Ax,Ay,Az,A0
complex(dp), dimension(nfft1/2+1,nfft2,nfft3) :: Axx_k,Ayy_k,Azz_k,Axy_k,Axz_k,Ayz_k,Ax_k,Ay_k,Az_k,A0_k
real(dp) :: fk1,rmax1,fk2,rho_temp,psi
integer(i2b)::nmax1x,nmax1y,nmax1z,ix,iy,iz,nmax2x,nmax2y,nmax2z
real(dp)::deltaVk,Hxxpreviousstep
real(dp)::fw,f_ww,F3B1,F3B2,costheta0
real(dp)::rb
complex(dp), dimension(nfft1/2+1,nfft2,nfft3) ::Gxx_k,Gyy_k,Gzz_k,Gxy_k,Gxz_k,Gyz_k,Gx_k,Gy_k,Gz_k,G0_k , function_rho_0k
real(dp), dimension(nfft1,nfft2,nfft3) :: FGxx,FGyy,FGzz,FGxy,FGxz,FGyz,FGx,FGy,FGz,FG0, function_rho_0
real(dp) :: lambda_w , F3B_ww, rmax_w!lambda parameter for water water interaction
!real(dp), dimension(nfft1,nfft2,nfft3) :: FAxx,FAyy,FAzz,FAxy,FAyz,FAxz,FAx,FAy,FAz,FA0
!integer(kind=i2B) ::nmax_wx, nmax_wy, nmax_wz ! nmax for water water interactions along x y z
F3B1=0.0_dp
F3B2=0.0_dp
costheta0=-1.0_dp/3.0_dp
deltaVk=(twopi)**3/(Lx*Ly*Lz)
!lambda_w=20.0_dp
! check if user wants to use this part of the functional
do i = 1 , size ( input_line )
  j = len ( 'threebody' )
  if ( input_line (i) ( 1 : j ) == 'threebody' .and. input_line ( i ) ( j + 4 : j + 4 ) == 'F' ) return
END DO
if (.not. input_log('F3B_new')) return
!> start timer
call cpu_time(time0)
! get density
allocate(rho(nfft1,nfft2,nfft3))
rho=0.0_dp
icg=0
do i=1,nfft1
  do j=1,nfft2
    do k=1,nfft3
      do o=1,angGrid%n_angles
        do p=1, molRotGrid%n_angles
          icg=icg+1
          rho(i,j,k) = rho(i,j,k) + rho_0*angGrid%weight(o)*molRotGrid%weight(p)*cg_vect(icg)**2
        END DO
      END DO
    END DO
  END DO
END DO
!rho_temp=rho
!rho=0.0_dp
      Axx=0.0_dp
      Ayy=0.0_dp
      Azz=0.0_dp
      Axy=0.0_dp
      Axz=0.0_dp
      Ayz=0.0_dp
      Ax=0.0_dp
      Ay=0.0_dp
      Az=0.0_dp
      A0=0.0_dp
      Axx_k=0.0_dp
      Ayy_k=0.0_dp
      Azz_k=0.0_dp
      Axy_k=0.0_dp
      Axz_k=0.0_dp
      Ayz_k=0.0_dp
      Ax_k=0.0_dp
      Ay_k=0.0_dp
      Az_k=0.0_dp
      A0_K=0.0_dp
do i=1, nfft1
  i1=i-1
  if ( i > nfft1/2 ) i1=i-1-nfft1
  do j=1, nfft2
     j1=j-1
     if (j > nfft2/2) j1=j-1-nfft2
    do k=1, nfft3
       k1=k-1
       if (k>nfft3/2) k1=k-1-nfft3
      x=i1*deltax
      y=j1*deltay
      z=k1*deltaz
      r=sqrt(x**2+y**2+z**2)
      fw=f_ww(r,rmin2,rsw2,rmax2)
        if (r/=0.0_dp) then
      A0(i,j,k)=fw
        Axx(i,j,k)=fw*x**2/(r**2)
        Ayy(i,j,k)=fw*y**2/(r**2)
        Azz(i,j,k)=fw*z**2/(r**2)
        Axy(i,j,k)=fw*x*y/(r**2)
        Axz(i,j,k)=fw*x*z/(r**2)
        Ayz(i,j,k)=fw*z*y/(r**2)
        Ax(i,j,k)=fw*x/r
        Ay(i,j,k)=fw*y/r
        Az(i,j,k)=fw*z/r
       END IF
    END DO     !k,z
  END DO       !j,y
END DO         !i,x
!stop
allocate(rho_k(nfft1/2+1,nfft2,nfft3))
fftw3%in_forward=rho
call dfftw_execute(fftw3%plan_forward)
rho_k=fftw3%out_forward*deltaV
      fftw3%in_forward=Axx
      call dfftw_execute(fftw3%plan_forward)
      Axx_k=fftw3%out_forward*deltaV
      fftw3%in_forward=Ayy
      call dfftw_execute(fftw3%plan_forward)
      Ayy_k=fftw3%out_forward*deltaV
      fftw3%in_forward=Azz
      call dfftw_execute(fftw3%plan_forward)
      Azz_k=fftw3%out_forward*deltaV
      fftw3%in_forward=Axy
      call dfftw_execute(fftw3%plan_forward)
      Axy_k=fftw3%out_forward*deltaV
      fftw3%in_forward=Axz
      call dfftw_execute(fftw3%plan_forward)
      Axz_k=fftw3%out_forward*deltaV
      fftw3%in_forward=Ayz
      call dfftw_execute(fftw3%plan_forward)
      Ayz_k=fftw3%out_forward*deltaV
      fftw3%in_forward=Ax
      call dfftw_execute(fftw3%plan_forward)
      Ax_k=fftw3%out_forward*deltaV
      fftw3%in_forward=Ay
      call dfftw_execute(fftw3%plan_forward)
      Ay_k=fftw3%out_forward*deltaV
      fftw3%in_forward=Az
      call dfftw_execute(fftw3%plan_forward)
      Az_k=fftw3%out_forward*deltaV
      fftw3%in_forward=A0
      call dfftw_execute(fftw3%plan_forward)
      A0_k=fftw3%out_forward*deltaV
      Fxx_k=Axx_k*rho_k
      Fyy_k=Ayy_k*rho_k
      Fzz_k=Azz_k*rho_k
      Fxy_k=Axy_k*rho_k
      Fxz_k=Axz_k*rho_k
      Fyz_k=Ayz_k*rho_k
      Fx_k=Ax_k*rho_k
      Fy_k=Ay_k*rho_k
      Fz_k=Az_k*rho_k
      F0_k=A0_k*rho_k
      fftw3%in_backward=Fxx_k
      call dfftw_execute(fftw3%plan_backward)
      Fxx=fftw3%out_backward*deltaVk/(twopi)**3
      fftw3%in_backward=Fyy_k
      call dfftw_execute(fftw3%plan_backward)
      Fyy=fftw3%out_backward*deltaVk/(twopi)**3
      fftw3%in_backward=Fzz_k
      call dfftw_execute(fftw3%plan_backward)
      Fzz=fftw3%out_backward*deltaVk/(twopi)**3
      fftw3%in_backward=Fxy_k
      call dfftw_execute(fftw3%plan_backward)
      Fxy=fftw3%out_backward*deltaVk/(twopi)**3
      fftw3%in_backward=Fxz_k
      call dfftw_execute(fftw3%plan_backward)
      Fxz=fftw3%out_backward*deltaVk/(twopi)**3
      fftw3%in_backward=Fyz_k
      call dfftw_execute(fftw3%plan_backward)
      Fyz=fftw3%out_backward*deltaVk/(twopi)**3
      fftw3%in_backward=Fx_k
      call dfftw_execute(fftw3%plan_backward)
      Fx=fftw3%out_backward*deltaVk/(twopi)**3
      fftw3%in_backward=Fy_k
      call dfftw_execute(fftw3%plan_backward)
      Fy=fftw3%out_backward*deltaVk/(twopi)**3
      fftw3%in_backward=Fz_k
      call dfftw_execute(fftw3%plan_backward)
      Fz=fftw3%out_backward*deltaVk/(twopi)**3
      fftw3%in_backward=F0_k
      call dfftw_execute(fftw3%plan_backward)
      F0=fftw3%out_backward*deltaVk/(twopi)**3
!      function_rho_0=n_0
!      in_forward=function_rho_0
!      call dfftw_execute(plan_forward)
!      function_rho_0k=out_forward
!      in_backward=(rho_k-function_rho_0k)*Axx_k
!      call dfftw_execute(plan_backward)
!      FAxx=out_backward*deltaVk/(twopi)**3
!      in_backward=(rho_k-function_rho_0k)*Ayy_k
!      call dfftw_execute(plan_backward)
!      FAyy=out_backward*deltaVk/(twopi)**3
!      in_backward=(rho_k-function_rho_0k)*Azz_k
!      call dfftw_execute(plan_backward)
!      FAzz=out_backward*deltaVk/(twopi)**3
!      in_backward=(rho_k-function_rho_0k)*Axy_k
!      call dfftw_execute(plan_backward)
!      FAxy=out_backward*deltaVk/(twopi)**3
!      in_backward=(rho_k-function_rho_0k)*Axz_k
!      call dfftw_execute(plan_backward)
!      FAxz=out_backward*deltaVk/(twopi)**3
!      in_backward=(rho_k-function_rho_0k)*Ayz_k
!      call dfftw_execute(plan_backward)
!      FAyz=out_backward*deltaVk/(twopi)**3
!      in_backward=(rho_k-function_rho_0k)*Ax_k
!      call dfftw_execute(plan_backward)
!      FAx=out_backward*deltaVk/(twopi)**3
!      in_backward=(rho_k-function_rho_0k)*Ay_k
!      call dfftw_execute(plan_backward)
!      FAy=out_backward*deltaVk/(twopi)**3
!      in_backward=(rho_k-function_rho_0k)*Az_k
!      call dfftw_execute(plan_backward)
!      FAz=out_backward*deltaVk/(twopi)**3
!      in_backward=(rho_k-function_rho_0k)*A0_k
!      call dfftw_execute(plan_backward)
!      FA0=out_backward*deltaVk/(twopi)**3
!  
Hxx=0.0_dp
Hyy=0.0_dp
Hzz=0.0_dp
Hxy=0.0_dp
Hxz=0.0_dp
Hyz=0.0_dp
Hx=0.0_dp
Hy=0.0_dp
Hz=0.0_dp
H0=0.0_dp
Gxx=0.0_dp
Gyy=0.0_dp
Gzz=0.0_dp
Gxy=0.0_dp
Gxz=0.0_dp
Gyz=0.0_dp
Gx=0.0_dp
Gy=0.0_dp
Gz=0.0_dp
G0=0.0_dp
do n=1, nb_solute_sites
if (lambda1_mol(n)/=0.0_dp) then
rmax1=0.5_dp*(sig_mol(id_mol(n))+sig_solv(1)) + d_w 
nmax1x = int(rmax1/deltax)
nmax1y = int(rmax1/deltay)
nmax1z = int(rmax1/deltaz)
ix = int(x_mol(n)/deltax) + 1
iy = int(y_mol(n)/deltay) + 1
iz = int(z_mol(n)/deltaz) + 1
  do i=ix-nmax1x, ix+nmax1x
    do j=iy-nmax1y, iy+nmax1y
      do k=iz-nmax1z,iz+nmax1z
        rk2=sqrt((x_mol(n)-(i-1)*deltax)**2+(y_mol(n)-(j-1)*deltay)**2+(z_mol(n)-(k-1)*deltaz)**2)
        xk2=-(x_mol(n)-(i-1)*deltax)
        yk2=-(y_mol(n)-(j-1)*deltay)
        zk2=-(z_mol(n)-(k-1)*deltaz)
        fk1=f_ww(rk2,rmin1,rsw1,rmax1)
       
        if (rk2 /= 0.0_dp )then
          rho_temp=rho(i,j,k)
          DHxx_ijk=deltaV*fk1*xk2**2/(rk2**2)
          DHyy_ijk=deltaV*fk1*yk2**2/(rk2**2)
          DHzz_ijk=deltaV*fk1*zk2**2/(rk2**2)
          DHxy_ijk=deltaV*fk1*xk2*yk2/(rk2**2)
          DHxz_ijk=deltaV*fk1*xk2*zk2/(rk2**2)
          DHyz_ijk=deltaV*fk1*yk2*zk2/(rk2**2)
          DHx_ijk=deltaV*fk1*xk2/rk2
          DHy_ijk=deltaV*fk1*yk2/rk2
          DHz_ijk=deltaV*fk1*zk2/rk2
          DH0_ijk=fk1*deltaV 
          Hxx(n)=Hxx(n)+DHxx_ijk*rho_temp
          Hyy(n)=Hyy(n)+DHyy_ijk*rho_temp
          Hzz(n)=Hzz(n)+DHzz_ijk*rho_temp
          Hxy(n)=Hxy(n)+DHxy_ijk*rho_temp
          Hxz(n)=Hxz(n)+DHxz_ijk*rho_temp
          Hyz(n)=Hyz(n)+DHyz_ijk*rho_temp
          Hx(n)=Hx(n)+DHx_ijk*rho_temp
          Hy(n)=Hy(n)+DHy_ijk*rho_temp
          Hz(n)=Hz(n)+DHz_ijk*rho_temp
          H0(n)=H0(n)+DH0_ijk*rho_temp
          
          DHxx(i,j,k,n)=DHxx_ijk
          DHyy(i,j,k,n)=DHyy_ijk
          DHzz(i,j,k,n)=DHzz_ijk
          DHxy(i,j,k,n)=DHxy_ijk
          DHxz(i,j,k,n)=DHxz_ijk
          DHx(i,j,k,n)= DHx_ijk 
          DHy(i,j,k,n)= DHy_ijk 
          DHz(i,j,k,n)= DHz_ijk
          DHyz(i,j,k,n)=DHyz_ijk
          DH0(i,j,k,n)= DH0_ijk         
        END IF
      END DO
    END DO
  END DO
END IF
END DO
do n=1, nb_solute_sites
if (lambda1_mol(n)/=0.0_dp) then
nmax2x = int(rmax1/deltax)
nmax2y = int(rmax1/deltay)
nmax2z = int(rmax1/deltaz)
ix = int(x_mol(n)/deltax) + 1
iy = int(y_mol(n)/deltay) + 1
iz = int(z_mol(n)/deltaz) + 1
  do i=ix-nmax2x, ix+nmax2x+1
    do j=iy-nmax2y, iy+nmax2y+1
      do k=iz-nmax2z,iz+nmax2z+1
      rk2=sqrt((x_mol(n)-(i-1)*deltax)**2+(y_mol(n)-(j-1)*deltay)**2+(z_mol(n)-(k-1)*deltaz)**2)
      xk2=-(x_mol(n)-(i-1)*deltax)
      yk2=-(y_mol(n)-(j-1)*deltay)
      zk2=-(z_mol(n)-(k-1)*deltaz)
      fk2=f_ww(rk2,rmin1,rsw1,rmax1)
          G0(i,j,k)=G0(i,j,k)+lambda2_mol(n)*fk2
!        if (rk2 /= 0.0_dp)then
          Gxx(i,j,k)=Gxx(i,j,k)+lambda2_mol(n)*fk2*xk2**2/(rk2**2)
          Gyy(i,j,k)=Gyy(i,j,k)+lambda2_mol(n)*fk2*yk2**2/(rk2**2)
          Gzz(i,j,k)=Gzz(i,j,k)+lambda2_mol(n)*fk2*zk2**2/(rk2**2)
          Gxy(i,j,k)=Gxy(i,j,k)+lambda2_mol(n)*fk2*xk2*yk2/(rk2**2)
          Gxz(i,j,k)=Gxz(i,j,k)+lambda2_mol(n)*fk2*xk2*zk2/(rk2**2)
          Gyz(i,j,k)=Gyz(i,j,k)+lambda2_mol(n)*fk2*yk2*zk2/(rk2**2)
          Gx(i,j,k)=Gx(i,j,k)+lambda2_mol(n)*fk2*(xk2)/rk2
          Gy(i,j,k)=Gy(i,j,k)+lambda2_mol(n)*fk2*(yk2)/rk2
          Gz(i,j,k)=Gz(i,j,k)+lambda2_mol(n)*fk2*(zk2)/rk2
!if (Gx(i,j,k)/=0.0_dp) then
!print*, i,j,j,Gx(i,j,k)
!END IF
!        END IF
      END DO
    END DO
  END DO
END IF
END DO
Gxx_k=0.0_dp
Gyy_k=0.0_dp
Gzz_k=0.0_dp
Gxy_k=0.0_dp
Gyz_k=0.0_dp
Gxz_k=0.0_dp
Gx_k=0.0_dp
Gy_k=0.0_dp
Gz_k=0.0_dp
G0_k=0.0_dp
FGxx=0.0_dp
FGyy=0.0_dp
FGzz=0.0_dp
FGxy=0.0_dp
FGyz=0.0_dp
FGxz=0.0_dp
FGx=0.0_dp
FGy=0.0_dp
FGz=0.0_dp
FG0=0.0_dp
fftw3%in_forward=Gxx*rho
call dfftw_execute(fftw3%plan_forward)
Gxx_k=fftw3%out_forward*deltaV
fftw3%in_forward=Gyy*rho
call dfftw_execute(fftw3%plan_forward)
Gyy_k=fftw3%out_forward*deltaV
fftw3%in_forward=Gzz*rho
call dfftw_execute(fftw3%plan_forward)
Gzz_k=fftw3%out_forward*deltaV
fftw3%in_forward=Gxy*rho
call dfftw_execute(fftw3%plan_forward)
Gxy_k=fftw3%out_forward*deltaV
fftw3%in_forward=Gxz*rho
call dfftw_execute(fftw3%plan_forward)
Gxz_k=fftw3%out_forward*deltaV
fftw3%in_forward=Gyz*rho
call dfftw_execute(fftw3%plan_forward)
Gyz_k=fftw3%out_forward*deltaV
fftw3%in_forward=Gx*rho
call dfftw_execute(fftw3%plan_forward)
Gx_k=fftw3%out_forward*deltaV
fftw3%in_forward=Gy*rho
call dfftw_execute(fftw3%plan_forward)
Gy_k=fftw3%out_forward*deltaV
fftw3%in_forward=Gz*rho
call dfftw_execute(fftw3%plan_forward)
Gz_k=fftw3%out_forward*deltaV
fftw3%in_forward=G0*rho
call dfftw_execute(fftw3%plan_forward)
G0_k=fftw3%out_forward*deltaV
Gxx_k=Axx_k*Gxx_k
Gyy_k=Ayy_k*Gyy_k
Gzz_k=Azz_k*Gzz_k
Gxy_k=Axy_k*Gxy_k
Gxz_k=Axz_k*Gxz_k
Gyz_k=Ayz_k*Gyz_k
Gx_k=Ax_k*Gx_k
Gy_k=Ay_k*Gy_k
Gz_k=Az_k*Gz_k
G0_k=A0_k*G0_k
fftw3%in_backward=Gxx_k
call dfftw_execute(fftw3%plan_backward)
FGxx=fftw3%out_backward*deltaVk/(twopi)**3
fftw3%in_backward=Gyy_k
call dfftw_execute(fftw3%plan_backward)
FGyy=fftw3%out_backward*deltaVk/(twopi)**3
fftw3%in_backward=Gzz_k
call dfftw_execute(fftw3%plan_backward)
FGzz=fftw3%out_backward*deltaVk/(twopi)**3
fftw3%in_backward=Gxy_k
call dfftw_execute(fftw3%plan_backward)
FGxy=fftw3%out_backward*deltaVk/(twopi)**3
fftw3%in_backward=Gxz_k
call dfftw_execute(fftw3%plan_backward)
FGxz=fftw3%out_backward*deltaVk/(twopi)**3
fftw3%in_backward=Gyz_k
call dfftw_execute(fftw3%plan_backward)
FGyz=fftw3%out_backward*deltaVk/(twopi)**3
fftw3%in_backward=Gx_k
call dfftw_execute(fftw3%plan_backward)
FGx=fftw3%out_backward*deltaVk/(twopi)**3
fftw3%in_backward=Gy_k
call dfftw_execute(fftw3%plan_backward)
FGy=fftw3%out_backward*deltaVk/(twopi)**3
fftw3%in_backward=Gz_k
call dfftw_execute(fftw3%plan_backward)
FGz=fftw3%out_backward*deltaVk/(twopi)**3
fftw3%in_backward=G0_k
call dfftw_execute(fftw3%plan_backward)
FG0=fftw3%out_backward*deltaVk/(twopi)**3
!print*,maxval(aimag(G0_k)),maxval(aimag(Gxx_k)),maxval(aimag(Gyy_k)),maxval(aimag(Gzz_k)),maxval(aimag(Gxy_k))
!print*,maxval(aimag(Gyz_k)),maxval(aimag(Gxz_k)),maxval(aimag(Gx_k)),maxval(aimag(Gy_k)),maxval(aimag(Gz_k))
!print*,maxval(real(G0_k)),maxval(real(Gxx_k)),maxval(real(Gyy_k)),maxval(real(Gzz_k)),maxval(real(Gxy_k))
!print*,maxval(real(Gyz_k)),maxval(real(Gxz_k)),maxval(real(Gx_k)),maxval(real(Gy_k)),maxval(real(Gz_k))
!print*,maxval(G0),maxval(Gxx),maxval(Gyy),maxval(Gzz),maxval(Gxy),maxval(Gyz),maxval(Gxz),maxval(Gx),maxval(Gy),maxval(Gz)
!print*,maxval(FG0),maxval(FGxx),maxval(FGyy),maxval(FGzz),maxval(FGxy),maxval(FGyz),maxval(FGxz),maxval(FGx),maxval(FGy),maxval(FGz)
!stop
do n=1, nb_solute_sites
  F3B1=F3B1+lambda1_mol(n)*kBT*0.5_dp*((Hxx(n))**2+(Hyy(n))**2+(Hzz(n))**2+&
            2.0_dp*(Hxy(n))**2+2.0_dp*(Hxz(n))**2+2.0_dp*(Hyz(n))**2&
           -2.0_dp*costheta0*((Hx(n))**2+(Hy(n))**2+(Hz(n))**2)+costheta0**2*(H0(n))**2)
END DO
do i=ix-nmax2x, ix+nmax2x
  do j=iy-nmax2y, iy+nmax2y
    do k=iz-nmax2z,iz+nmax2z
       F3B2=F3B2+(kBT*0.5_dp*(Fxx(i,j,k)*Gxx(i,j,k)+Fyy(i,j,k)*Gyy(i,j,k)+Fzz(i,j,k)*Gzz(i,j,k)&
                +2.0_dp*Fxy(i,j,k)*Gxy(i,j,k)+ 2.0_dp*Fxz(i,j,k)*Gxz(i,j,k)+ 2.0_dp*Fyz(i,j,k)*Gyz(i,j,k))&
                - 2.0_dp*costheta0*kBT*0.5_dp*(Fx(i,j,k)*Gx(i,j,k)+Fy(i,j,k)*Gy(i,j,k)+Fz(i,j,k)*Gz(i,j,k))&
                +costheta0**2*F0(i,j,k)*G0(i,j,k)*kBT*0.5_dp)*deltaV*rho(i,j,k)
    
    END DO
  END DO
END DO
!F3B_ww=0.0_dp
!do i=1,nfft1
!  do j=1, nfft2
!    do k=1, nfft3
!      F3B_ww=F3B_ww+0.5_dp*kbT*deltaV*rho(i,j,k)*lambda_w*(FAxx(i,j,k)**2+FAyy(i,j,k)**2+FAzz(i,j,k)**2+&
!             2.0_dp*(FAxy(i,j,k)**2+FAxz(i,j,k)**2+FAyz(i,j,k)**2)-2.0_dp*costheta0*(FAx(i,j,k)**2+FAy(i,j,k)**2+FAz(i,j,k)**2)&
!             +costheta0**2*FA0(i,j,k)**2)
!    END DO
!  END DO
!END DO
print*,'F3B_ww = ', F3B_ww
icg=0
open(12,file='output/dF_2S_new.dat')
do i=1,nfft1
  do j=1, nfft2
    do k=1, nfft3
      do o=1,angGrid%n_angles
        do p=1,molRotGrid%n_angles
          icg=icg+1    
          psi=cg_vect(icg)
          dF(icg)=dF(icg)+kBT*psi*deltaV*angGrid%weight(o)*molRotGrid%weight(p)*rho_0*(&
                  (Fxx(i,j,k)*Gxx(i,j,k)+Fyy(i,j,k)*Gyy(i,j,k)+Fzz(i,j,k)*Gzz(i,j,k)&
                  +2.0_dp*Fxy(i,j,k)*Gxy(i,j,k)+ 2.0_dp*Fxz(i,j,k)*Gxz(i,j,k)+ 2.0_dp*Fyz(i,j,k)*Gyz(i,j,k))&
                  -2.0_dp*costheta0*(Fx(i,j,k)*Gx(i,j,k)+Fy(i,j,k)*Gy(i,j,k)+Fz(i,j,k)*Gz(i,j,k))&
                  +costheta0**2*F0(i,j,k)*G0(i,j,k))+&
                  kBT*psi*angGrid%weight(o)*molRotGrid%weight(p)*&
                  rho_0*deltaV*(FGxx(i,j,k)+FGyy(i,j,k)+FGzz(i,j,k)+2.0_dp*FGxy(i,j,k)&
                 +2.0_dp*FGxz(i,j,k)+2.0_dp*FGyz(i,j,k)+2.0_dp*costheta0*(FGx(i,j,k)+FGy(i,j,k)+FGz(i,j,k))+costheta0**2*FG0(i,j,k))
!          dF(icg)=dF(icg)+kbT*deltaV*angGrid%weight(o)*molRotGrid%weight(p)*rho_0*lambda_w*psi*((Fxx(i,j,k)**2+Fyy(i,j,k)**2+Fzz(i,j,k)**2+&
!                  2.0_dp*(Fxy(i,j,k)**2+Fxz(i,j,k)**2+Fyz(i,j,k)**2)-2.0_dp*costheta0*(Fx(i,j,k)**2+Fy(i,j,k)**2+Fz(i,j,k)**2)&
!                  +costheta0**2*F0(i,j,k)**2)&
!                  +2.0_dp*(FAxx(i,j,k)+FAyy(i,j,k)+FAzz(i,j,k)+2.0_dp*FAxy(i,j,k)+2.0_dp*FAyz(i,j,k)+2.0_dp*FAxz(i,j,k)&
!                  -2.0_dp*costheta0*(FAx(i,j,k)+FAy(i,j,k)+FAz(i,j,k))+costheta0**2*FA0(i,j,k) ))
!                  
!          dF(icg)=dF(icg)+kbT*deltaV*angGrid%weight(o)*molRotGrid%weight(p)*rho_0*lambda_w*psi*((FAxx(i,j,k)**2+FAyy(i,j,k)**2+FAzz(i,j,k)**2+&
!                 2.0_dp*(FAxy(i,j,k)**2+FAxz(i,j,k)**2+FAyz(i,j,k)**2)-2.0_dp*costheta0*(FAx(i,j,k)**2+FAy(i,j,k)**2+FAz(i,j,k)**2)&
!                  +costheta0**2*FA0(i,j,k)**2)&
!                  +2.0_dp*(FAxx(i,j,k)*Axx(i,j,k)+FAyy(i,j,k)*Ayy(i,j,k)+FAzz(i,j,k)*Azz(i,j,k)+2.0_dp*FAxy(i,j,k)*Axy(i,j,k)+&
!                  2.0_dp*FAyz(i,j,k)*Ayz(i,j,k)+2.0_dp*FAxz(i,j,k)*Axz(i,j,k)&
!              -2.0_dp*costheta0*(FAx(i,j,k)*Ax(i,j,k)+FAy(i,j,k)*Ay(i,j,k)+FAz(i,j,k)*Az(i,j,k))+costheta0**2*FA0(i,j,k)*A0(i,j,k)))
                  
          do n=1,nb_solute_sites
            dF(icg)=dF(icg)+lambda1_mol(n)*kBT*psi*angGrid%weight(o)*&
            molRotGrid%weight(p)*((Hxx(n)*DHxx(i,j,k,n)+Hyy(n)*DHyy(i,j,k,n)+&
Hzz(n)*DHzz(i,j,k,n)&
            +2.0_dp*Hxy(n)*DHxy(i,j,k,n)+2.0_dp*Hxz(n)*DHxz(i,j,k,n)+2.0_dp*Hyz(n)*DHyz(i,j,k,n))&
      -2.0_dp*costheta0*(Hx(n)*DHx(i,j,k,n)+Hy(n)*DHy(i,j,k,n)+Hz(n)*DHz(i,j,k,n))+costheta0**2*H0(n)*DH0(i,j,k,n))*rho_0*2.0_dp
          END DO
         
        END DO
      END DO
    END DO
  END DO
END DO
close(12)
call cpu_time(time1)
print*, 'F3B1=', F3B1
print*, 'F3B2=', F3B2
print*, 'in' , time1-time0, 'sec'
FF=FF+F3B2+F3B1!+F3B_ww
END SUBROUTINE

function f_ww( r , rmin, rsw, rmax )
    USE precision_kinds, only: dp,i2b
    IMPLICIT NONE
    real(dp) :: f_ww, r, rmin, rsw, rmax
    real(dp) , parameter :: gam = 2.0_dp/3.0_dp
    real(dp) :: deltar, exp_term, Switch
    deltar = rsw-rmin 
    if(r > rmin .and. r<rmax) then
        exp_term = exp(gam*rmax/(r-rmax))
        if(r<rsw) then
            Switch = (r-rmin)**2*(-2.0_dp*(r-rmin)/deltar+3.0_dp)/deltar**2
        ELSE
            Switch = 1.0_dp
        END IF
        f_ww = Switch*exp_term
    ELSE
        f_ww = 0.0_dp
    END IF
    RETURN
end function f_ww

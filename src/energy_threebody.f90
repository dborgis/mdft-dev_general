!-----------3-BODY TERM-------------------------------
! This routine works for molecules with several specific sites
! giving rise to H-bonding (O, OH, N, NH, ...)
! Corresponding parameters are specified in tables lambda1_mol, lambda2_mol which
! have dimension nb_solute_sites
! See the paper by Borgis et al. JCP 134 194102 (2011), a_1 = rmax1, a_2 = rmax2
subroutine energy_threebody
use precision_kinds, only: dp,i2b
use system, only: nfft1 , nfft2 , nfft3 , deltaV , rho_0_multispec , sig_mol , sig_solv , Lx , Ly , Lz ,&
&    id_mol, x_mol , y_mol , z_mol , kbT , nb_species, nb_solute_sites, deltax, deltay, deltaz&
& , lambda1_mol , lambda2_mol, nb_psi
use cg, only: cg_vect,FF,dF
use quadrature, only: weight , weight_psi, angGrid
use constants, only: fourpi
use input , only : input_line, input_log
implicit none
! 3-BODY PARAMETERS
real(dp), parameter :: rmin1 = 1.5_dp, rsw1 = 2.0_dp, rmin2 = 2.25_dp, rsw2 = 2.5_dp, rmax2 = 5.0_dp, d_w = 1.9_dp
real(dp), parameter :: cos_theta0 = -1.0_dp/3.0_dp   ! cos of the 3 body angle one would like (for instance cos_theta0=cos(109.5°) in water
real(dp) :: rmax1 ! distance constraints from the solute site to the 2nd shell
real(dp), allocatable, dimension(:,:,:) :: rho
real(dp), allocatable, dimension(:,:,:) :: V3 ! dF3/drho
real(dp), allocatable, dimension(:,:,:) :: V3_fs ! V3_cor
real(dp) :: F3 ! total energy due to 3body interactions
real(dp) :: F3_fs ! first shell part of F3
real(dp) :: fw1, fw2, fw12, f_w
real(dp) :: energy_3b
real(dp) :: fact_n !>@var integration constant for the double sum over all positions = fact**2
real(dp) :: x1, y1, z1, x2, y2, z2 ! coordinates for first and second integrations
real(dp) :: x12, y12, z12, r12 ! projections and distance between first and second shell
real(dp) :: r1, r2 ! distance between solvent site and unique solute
real(dp) :: cos_theta ! cos of the effective angle between 3bodies
integer(i2b) :: i1, j1, k1, i2, j2, k2, ix, iy, iz, n!>@var dummy
!real(dp) :: a40, a3, b3, rho_s   ! 3-body parameters
integer(i2b) :: nmax1x, nmax1y, nmax1z, nmax2x, nmax2y, nmax2z
real(dp) :: psi ! psi**2 = rho(r)
integer(i2b) :: i , j , k , icg , o , p! dummy for loops
real(dp) :: time0,time1 ! timer steps
!real(dp),  dimension(nb_solute_sites) :: lambda1_mol, lambda2_mol
! check if user wants to use this part of the functional
if (.not. input_log('threebody')  ) return
if (.not. input_log('F3B_old')) return
!> start timer
call cpu_time(time0)
!******************A CHANGER******************************************
!> Should be read in input/dft.in and solute.in. Here for test purpose only.
!  lambda1_mol and lambda1_mol should have dimension nb_solute_sites
! FOR CH3OH only 
!lambda1_mol = 0.0_dp
!lambda2_mol = 0.0_dp
!lambda1_mol(2)=30.0_dp
!lambda2_mol(2)=30.0_dp
!lambda1_mol(4)=45.0_dp
!lambda2_mol(4)=45.0_dp
!**********************************************************************
fact_n=(deltaV*rho_0_multispec(1) )**2 !CARE TODO HERE ONLY VALID FOR NB_SPECIES =1
if ( nb_species /= 1 ) then
  print *, 'error in energy_thhreebody.f90. You try to compute it but it is only valid for 1 component fluid'
  stop
end if
! get density
allocate(rho(nfft1,nfft2,nfft3))
rho=0.0_dp
icg=0
do i=1,nfft1
  do j=1,nfft2
    do k=1,nfft3
      do o=1,angGrid%n_angles
        do p=1, nb_psi
          icg=icg+1
          rho(i,j,k) = rho(i,j,k) + weight(o)*weight_psi(p)*cg_vect(icg)**2
        end do
      end do
    end do
  end do
end do
!> First shell
allocate(V3_fs(nfft1,nfft2,nfft3))
V3_fs = 0.0_dp ! array
F3_fs = 0.0_dp ! scalar
!> Second shell
allocate(V3(nfft1,nfft2,nfft3))
V3 = 0.0_dp ! array
F3 = 0.0_dp
!>define grid parameters for water-water
nmax2x = int(rmax2/deltax)
nmax2y = int(rmax2/deltay)
nmax2z = int(rmax2/deltaz)
!Loop over sites
DO n = 1, nb_solute_sites
if(lambda1_mol(n) > 1.0_dp) THEN
rmax1=0.5_dp*(sig_mol(id_mol(n))+sig_solv(1)) + d_w
nmax1x = int(rmax1/deltax)
nmax1y = int(rmax1/deltay)
nmax1z = int(rmax1/deltaz)
ix = int(x_mol(n)/deltax) + 1
iy = int(y_mol(n)/deltay) + 1
iz = int(z_mol(n)/deltaz) + 1
do k1 = iz - nmax1z, iz + nmax1z +1
    z1 = (k1-1)*deltaz - z_mol(n)  !TODO PBC. Here solute is supposed at center of BOX
      do j1 = iy - nmax1y, iy + nmax1y +1
         y1 = (j1-1)*deltay - y_mol(n)
         do i1 = ix-nmax1x,ix+nmax1x + 1
            x1 = (i1-1)*deltax - x_mol(n)
      r1 = sqrt(x1**2+y1**2+z1**2)
      if (r1>rmin1 .and. r1<rmax1) then
        
        fw1 = f_w(r1, rmin1, rsw1, rmax1)
!  First shell contribution
      do k2 = iz - nmax1z, iz + nmax1z +1
         z2 = (k2-1)*deltaz - z_mol(n)  !TODO PBC
          do j2 = iy - nmax1y, iy + nmax1y +1
            y2 = (j2-1)*deltay - y_mol(n)
             do i2 = ix-nmax1x,ix+nmax1x + 1
               x2 = (i2-1)*deltax - x_mol(n)
               r2 = sqrt(x2**2+y2**2+z2**2)
 
              if(r2>rmin1 .and. r2<rmax1) then
                cos_theta = (x1*x2+y1*y2+z1*z2)/(r1*r2)
                fw2 = f_w(r2, rmin1, rsw1, rmax1)
                energy_3b = 0.5_dp*fact_n*kBT*lambda1_mol(n)*(cos_theta-cos_theta0)**2*fw1*fw2
                F3_fs =F3_fs + energy_3b*rho(i2,j2,k2)*rho(i1,j1,k1)
                V3_fs(i1,j1,k1) = V3_fs(i1,j1,k1) + energy_3b*rho(i2,j2,k2)
                V3_fs(i2,j2,k2) = V3_fs(i2,j2,k2) + energy_3b*rho(i1,j1,k1)
              endif
            end do
          end do
        end do
   
! Second shell contribution
        ! second loop over all positions
        do i2 = i1-nmax2x,i1+nmax2x
          x12 = (i2-i1)*deltax
          do j2 = j1-nmax2y,j1+nmax2y
            y12 = (j2-j1)*deltay
            do k2 = k1-nmax2z,k1+nmax2z
              z12 = (k2-k1)*deltaz
              r12 = sqrt(x12**2+y12**2+z12**2)
              
              if(r12>rmin2 .and. r12<rmax2) then
              
                cos_theta = -(x1*x12+y1*y12+z1*z12)/r1/r12
                fw12 = f_w(r12, rmin2, rsw2, rmax2)
            
                energy_3b = 0.5_dp*fact_n*kBT*lambda2_mol(n)*(cos_theta-cos_theta0)**2*fw1*fw12
                                  
                F3 =F3 + energy_3b*rho(i2,j2,k2)*rho(i1,j1,k1)
                V3(i1,j1,k1) = V3(i1,j1,k1) + energy_3b*rho(i2,j2,k2)
                V3(i2,j2,k2) = V3(i2,j2,k2) + energy_3b*rho(i1,j1,k1)
                
              end if
            end do
          end do
        end do
          
      end if  
    end do
  end do
end do
END IF
END DO ! end loop over sites
!> The detail is useless for computing the whole energy and its gradient
!write(*,*) 'first shell=',F3_fs,'  second shell =', F3
print*, 'first_shell = ', F3_fs , 'second_shell', F3   !Rajouter par Guillaume le 19/11/12 à enlever!!!
F3 = F3 + F3_fs
V3 = V3+V3_fs!V3_fs + V3
! add gradient of the three body term to total gradient
open (11, file='output/dF_2S_old.dat')
icg = 0
do i = 1 , nfft1
  do j = 1 , nfft2
    do k = 1 , nfft3
      do o = 1, angGrid%n_angles
        do p=1, nb_psi
          icg = icg + 1
          psi = CG_vect ( icg )
          dF ( icg ) = dF ( icg ) + 2.0_dp * psi * weight ( o ) * weight_psi(p)* (V3 ( i , j , k ))
!           write(11,*) icg,  2.0_dp * psi * weight ( o ) * weight_psi(p)* (V3 ( i , j , k )-V3_fs(i,j,k))
        end do
      end do
    end do
  end do
end do
close(11)
! Add the threebody contributions to total energy FF
FF = FF + F3
! Warn user
call cpu_time ( time1 )
write(*,*) 'F3body      = ',F3,'computed in (sec)',time1-time0
end subroutine energy_threebody
function f_w( r , rmin, rsw, rmax )
use precision_kinds, only: dp,i2b
implicit none
real(dp) :: f_w, r, rmin, rsw, rmax
real(dp) , parameter :: gam = 2.0_dp/3.0_dp
real(dp) :: deltar, exp_term, Switch
deltar = rsw-rmin 
if(r > rmin .and. r<rmax) then
                
exp_term = exp(gam*rmax/(r-rmax))
   if(r<rsw) then
     Switch = (r-rmin)**2*(-2.0_dp*(r-rmin)/deltar+3.0_dp)/deltar**2
     else
     Switch = 1.0_dp
   end if
f_w = Switch*exp_term
else
f_w = 0.0_dp
end if
return
end function f_w

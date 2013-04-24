!-----------3-BODY TERM-------------------------------
! Apply only for one-site solute at this stage, or with respect to oxygen of water
! Needs to be generalized 

! in the paper by Borgis et al. JCP 134 194102 (2011), a_1 = rmax1, a_2 = rmax2
subroutine compute_3body_term(rho,F3,V3)
use precision_kinds, only: dp,i2b
use system, only: nfft1,nfft2,nfft3,deltaV,rho_0,sig_mol,sig_solv,Lx,Ly,Lz,x_mol,y_mol,z_mol,kbT
implicit none

real(kind=dp), dimension(nfft1,nfft2,nfft3) :: rho
real(kind=dp), dimension(nfft1,nfft2,nfft3) :: V3 ! grad
real(kind=dp), dimension(nfft1,nfft2,nfft3) :: V3_fs ! V3_cor
real(kind=dp) :: F3 ! total energy due to 3body interactions
real(kind=dp) :: F3_fs ! F3_cor
real(kind=dp) :: S1, S2, exp1, exp2
real(kind=dp) :: energy_3b
real(kind=dp) :: fact_n !>@var integration constant for the double sum over all positions = fact**2
real(kind=dp) :: x1, y1, z1, x2, y2, z2 ! coordinates for first and second integrations
real(kind=dp) :: x12, y12, z12, r12 ! projections and distance between first and second shell
real(kind=dp) :: r1, r2 ! distance between solvent site and unique solute
real(kind=dp) :: cos_theta ! cos of the effective angle between 3bodies
real(kind=dp) :: cos_theta0 ! cos of the 3 body angle one would like (for instance cos_theta0=cos(109.5°) in water
integer(kind=i2b) :: i1, j1, k1, i2, j2, k2 !>@var dummy
real(kind=dp) :: gam ! empirical, based on 
real(kind=dp) :: gam1 ! dummy 1st shell
real(kind=dp) :: gam2 ! dummy 2nd shell
real(kind=dp) :: lambda, lambda_fs
real(kind=dp) :: rmin1, rmax1 ! distance constraints from the solute site to the first shell
real(kind=dp) :: rmin2, rmax2 ! distance constraints from the solute site to the 2nd shell
!real(kind=dp) :: a40, a3, b3, rho_s   ! 3-body parameters
real(kind=dp) :: rsw1, rsw2, deltar1, deltar2  ! 3-body parameters
integer(kind=i2b) :: nmin1, nmax1, nmin2, nmax2
integer(kind=i2b) :: nf1,nf2,nf3
real(kind=dp) :: delta_r



!> Should be read in input/dft.in. Here for test purpose only.
lambda_fs=45.0_dp
lambda=30.0_dp
gam=2.0_dp/3.0_dp
rmin1=1.5_dp
rsw1=2.0_dp !2.0
rmax1=0.5_dp*(sig_mol(1)+sig_solv(1))+1.9_dp
rmin2=2.25_dp
rsw2=2.5_dp
rmax2=4.75_dp

!ORIGINAL VALUES
!lambda_fs=100.0_dp
!lambda=100.0_dp
!gam=2.0_dp/3.0_dp
!rmin1=1.5_dp
!rsw1=2.0_dp
!rmax1=0.5_dp*(sig_mol(1)+sig_solv(1))+1.9_dp
!rmin2=2.25_dp
!rsw2=2.5_dp
!rmax2=4.75_dp




nf1=nfft1/2
nf2=nfft2/2
nf3=nfft3/2
fact_n=(deltaV*rho_0)**2
gam1=gam*rmax1
gam2=gam*rmax2
delta_r = Lx/nfft1 ! valid only if cubic box OR if number of point per angstrom is constant
deltar1 = rsw1-rmin1
deltar2 = rsw2-rmin2
nmax1 = int(rmax1/delta_r)
nmin1 = int(rmin1/delta_r)
nmax2 = int(rmax2/delta_r)
nmin2 = int(rmin2/delta_r)
cos_theta0 = -1.0_dp/3.0_dp ! =cos(109.5°)





!> First shell
V3_fs = 0.0_dp ! array
F3_fs = 0.0_dp ! scalar
do i1 = nf1+1-nmax1,nf1+1+nmax1
  x1 = (i1-1)*delta_r - x_mol(1)
  do j1 = nf2+1-nmax1,nf2+1+nmax1
    y1 = (j1-1)*delta_r - y_mol(1)
    do k1 = nf3+1-nmax1, nf3+1+nmax1
      z1 = (k1-1)*delta_r - z_mol(1)

      r1 = sqrt(x1**2+y1**2+z1**2)

      if (r1<rmax1 .and. r1>rmin1) then
        exp1 = exp(gam1/(r1-rmax1))
        if(r1<rsw1) then
          S1 = (r1-rmin1)**2*(-2.0_dp*(r1-rmin1)/deltar1+3.0_dp)/deltar1**2
        else
          S1 = 1.0_dp
        end if

        do i2 = nf1+1-nmax1,nf1+1+nmax1
          x2 = (i2-1)*delta_r-x_mol(1)
          do j2 = nf2+1-nmax1,nf2+1+nmax1
            y2 = (j2-1)*delta_r-y_mol(1)
            do k2 = nf3+1-nmax1,nf3+1+nmax1
              z2 = (k2-1)*delta_r-z_mol(1)

              r2 = sqrt(x2**2+y2**2+z2**2)
 
              if(r2<rmax1 .and. r2>rmin1) then
                cos_theta = (x1*x2+y1*y2+z1*z2)/(r1*r2)
                exp2 = exp(gam1/(r2-rmax1))
                if(r2<rsw1) then
                  S2 = (r2-rmin1)**2*(-2.0_dp*(r2-rmin1)/deltar1+3.0_dp)/deltar1**2
                else
                  S2 = 1.0_dp
                end if
                energy_3b = 0.5_dp*fact_n*kBT*lambda_fs*(cos_theta-cos_theta0)**2*S1*S2*exp2*exp1
                F3_fs =F3_fs + energy_3b*rho(i2,j2,k2)*rho(i1,j1,k1)
                V3_fs(i1,j1,k1) = V3_fs(i1,j1,k1) + energy_3b*rho(i2,j2,k2)
                V3_fs(i2,j2,k2) = V3_fs(i2,j2,k2) + energy_3b*rho(i1,j1,k1)
              endif

            end do
          end do
        end do
             
      end if  

    end do
  end do
end do




! 2nd shell
F3 = 0.0_dp
V3 = 0.0_dp ! array
do i1 = nf1+1-nmax1,nf1+1+nmax1
  x1 = (i1-1)*delta_r - x_mol(1)
  do j1 = nf2+1-nmax1,nf2+1+nmax1
    y1 = (j1-1)*delta_r - y_mol(1)
    do k1 = nf3+1-nmax1, nf3+1+nmax1
      z1 = (k1-1)*delta_r - z_mol(1)

      r1 = sqrt(x1**2+y1**2+z1**2)

      if (r1<rmax1 .and. r1>rmin1) then ! rmin1 and rmax1 are readen in input/dft.in
        exp1 = exp(gam1/(r1-rmax1))
        if(r1<rsw1) then
          S1 = (r1-rmin1)**2*(-2.0_dp*(r1-rmin1)/deltar1+3.0_dp)/deltar1**2
        else
          S1 = 1.0_dp
        end if

        ! second loop over all positions
        do i2 = i1-nmax2,i1+nmax2
          x12 = (i2-i1)*delta_r
          do j2 = j1-nmax2,j1+nmax2
            y12 = (j2-j1)*delta_r
            do k2 = k1-nmax2,k1+nmax2
              z12 = (k2-k1)*delta_r

              r12 = sqrt(x12**2+y12**2+z12**2)

              if(r12<=rmax2 .and. r12>rmin2) then
                cos_theta = -(x1*x12+y1*y12+z1*z12)/r1/r12
                exp2 = exp(gam2/(r12-rmax2))
                if(r12<rsw2) then
                  S2 = (r12-rmin2)**2*(-2.0_dp*(r12-rmin2)/deltar2+3.0_dp)/deltar2**2
                else
                  S2 = 1.0_dp
                end if
                energy_3b = 0.5_dp*fact_n*kBT*lambda*(cos_theta-cos_theta0)**2*exp2*exp1*S1*S2
                F3 =F3 + energy_3b*rho(i2,j2,k2)*rho(i1,j1,k1)
                V3(i1,j1,k1) = V3(i1,j1,k1) + energy_3b*rho(i2,j2,k2)
                V3(i2,j2,k2) = V3(i2,j2,k2) + energy_3b*rho(i1,j1,k1)
              endif

            end do
          end do
        end do
             
      end if
    end do
  end do
end do



! Correction to reduce height of first pic
!allocate (V3_cor(nfft1,nfft2,nfft3))
!V3_cor = 0.0_dp ! array
!write(*,*) 'a3, rho_s=',a3,rho_s
!F3_cor =0.0_dp
!do i=1,nfft1
!  do j=1,nfft2
!    do  k=1,nfft3
!      ! if(rho(i,j,k).ge.fourpi*rho_s) then
!      F3_cor = F3_cor + a3*kBT*(rho(i,j,k)-fourpi*rho_s)**3
!      V3_cor(i,j,k) = 3.0*a3*kBT*(rho(i,j,k)-fourpi*rho_s)**2
!      ! end if
!    end do
!  end do
!end do
!write(*,*) 'Fcor=', F3_cor


!> The detail is useless for computing the whole energy and its gradient
F3 = F3 + F3_fs
V3 = V3 + V3_fs

end subroutine compute_3body_term

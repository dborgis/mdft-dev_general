SUBROUTINE energy_threebody_faster (F3B1,F3B2)
    
    USE precision_kinds ,ONLY: i2b, dp
    USE input           ,ONLY: verbose, input_dp
    USE constants       ,ONLY: twopi,zeroC
    USE quadrature      ,ONLY: angGrid, molRotGrid
    USE system          ,ONLY: thermocond, spaceGrid, soluteSite, solventSite, solvent
    USE minimizer       ,ONLY: cg_vect,dF,FF
    USE fft             ,ONLY: fftw3,kproj
    
    IMPLICIT NONE
    REAL(dp), INTENT(OUT)                       :: F3B1, F3B2
    REAL(dp), PARAMETER                         :: rmin1=1.5_dp, rsw1=2.0_dp, rmin2=2.25_dp, rsw2=2.5_dp, rmax2=5.0_dp, d_w=1.9_dp
    INTEGER(i2b)                                :: icg,i,j,k,o,p,n,i1,j1,k1,nfft1,nfft2,nfft3,nb_solute_sites
    REAL(dp)                                    :: rk2,xk2,yk2,zk2,r,x,y,z,deltaVk,rb,fw,deltaV,time0,time1,deltax,deltay,deltaz
    REAL(dp)                                    :: fk1,rmax1,fk2,rho_temp,psi
    REAL(dp), ALLOCATABLE, DIMENSION(:,:,:)     :: Gxx,Gyy,Gzz,Gxy,Gxz,Gyz,Gx,Gy,Gz,G0    
    REAL(dp), ALLOCATABLE, DIMENSION(:,:,:)     :: Fxx,Fyy,Fzz,Fxy,Fxz,Fyz,Fx,Fy,Fz,F0 
    REAL(dp), ALLOCATABLE, DIMENSION(:,:,:)     :: Axx,Ayy,Azz,Axy,Axz,Ayz,Ax,Ay,Az,A0
    REAL(dp), ALLOCATABLE, DIMENSION(:,:,:,:)   :: DHxx,DHyy,DHzz,DHxy,DHxz,DHyz,DH0,DHx,DHy,DHz 
    REAL(dp), ALLOCATABLE, DIMENSION(:)         :: Hxx,Hyy,Hzz,Hxy,Hxz,Hyz,Hx,Hy,Hz,H0
    REAL(dp), ALLOCATABLE, DIMENSION(:,:,:)     :: rho
    COMPLEX(dp), ALLOCATABLE, DIMENSION(:,:,:)  :: rho_k
    COMPLEX(dp), ALLOCATABLE, DIMENSION(:,:,:)  :: Axx_k,Ayy_k,Azz_k,Axy_k,Axz_k,Ayz_k,Ax_k,Ay_k,Az_k,A0_k
    INTEGER(i2b)                                :: nmax1x,nmax1y,nmax1z,ix,iy,iz,nmax2x,nmax2y,nmax2z
    REAL(dp), PARAMETER                         :: costheta0 = -1.0_dp/3.0_dp
    COMPLEX(dp) ,ALLOCATABLE, DIMENSION(:,:,:)  :: Gxx_k,Gyy_k,Gzz_k,Gxy_k,Gxz_k,Gyz_k,Gx_k,Gy_k,Gz_k,G0_k , function_rho_0k
    REAL(dp)    ,ALLOCATABLE, DIMENSION(:,:,:)  :: FGxx,FGyy,FGzz,FGxy,FGxz,FGyz,FGx,FGy,FGz,FG0
    REAL(dp)                                    :: lambda_w , F3B_ww, rmax_w !lambda parameter for water water interaction
    REAL(dp)    ,ALLOCATABLE, DIMENSION(:,:,:)  :: FAxx,FAyy,FAzz,FAxy,FAyz,FAxz,FAx,FAy,FAz,FA0
    
    !integer(kind=i2B) ::nmax_wx, nmax_wy, nmax_wz ! nmax for water water interactions along x y z
    nb_solute_sites = SIZE(soluteSite)
    deltaVk=(twopi)**3/PRODUCT(spaceGrid%length)
    lambda_w=input_dp ('lambda_solvent')!5.0_dp
    ! check if user wants to use this part of the functional

    CALL CPU_TIME(time0)
    ! TODO : add test that no lambda2(n) zero if lambda1(n)==0

    deltaV =spaceGrid%dV
    deltax =spaceGrid%dl(1)
    deltay =spaceGrid%dl(2)
    deltaz =spaceGrid%dl(3)
    nfft1 =spaceGrid%n_nodes(1)
    nfft2 =spaceGrid%n_nodes(2)
    nfft3 =spaceGrid%n_nodes(3)

    ALLOCATE(rho(nfft1,nfft2,nfft3), SOURCE=0._dp)
    icg=0
    DO i=1,nfft1
      DO j=1,nfft2
        DO k=1,nfft3
          DO o=1,angGrid%n_angles
            DO p=1, molRotGrid%n_angles
                icg=icg+1
                rho(i,j,k) = rho(i,j,k) + solvent(1)%rho0 * angGrid%weight(o) * molRotGrid%weight(p) * cg_vect(icg)**2
            END DO
          END DO
        END DO
      END DO
    END DO
    !rho_temp=rho
    !rho=0.0_dp


    fftw3%in_forward=rho
    CALL dfftw_execute(fftw3%plan_forward)
    ALLOCATE (rho_k(nfft1/2+1,nfft2,nfft3)  ,SOURCE=fftw3%out_forward*deltaV)

    fftw3%in_forward=solvent(1)%n0
    CALL dfftw_execute(fftw3%plan_forward)
    function_rho_0k=fftw3%out_forward


    ALLOCATE(Axx(nfft1,nfft2,nfft3), SOURCE=0.0_dp)
    ALLOCATE(Ayy(nfft1,nfft2,nfft3), SOURCE=0.0_dp)
    ALLOCATE(Azz(nfft1,nfft2,nfft3), SOURCE=0.0_dp)
    ALLOCATE(Axy(nfft1,nfft2,nfft3), SOURCE=0.0_dp)
    ALLOCATE(Axz(nfft1,nfft2,nfft3), SOURCE=0.0_dp)
    ALLOCATE(Ayz(nfft1,nfft2,nfft3), SOURCE=0.0_dp)
    ALLOCATE( Ax(nfft1,nfft2,nfft3), SOURCE=0.0_dp)
    ALLOCATE( Ay(nfft1,nfft2,nfft3), SOURCE=0.0_dp)
    ALLOCATE( Az(nfft1,nfft2,nfft3), SOURCE=0.0_dp)
    ALLOCATE( A0(nfft1,nfft2,nfft3), SOURCE=0.0_dp)


!!!!!!!

    BLOCK
        REAL(dp) :: dl(3),x,y,z,r,r2
        INTEGER(i2B) :: i,j,k
        dl = spaceGrid%length * spaceGrid%dl / twopi
        DO CONCURRENT (i=1:nfft1, j=1:nfft2, k=1:nfft3)
            x = kproj(1,i) * dl(1)
            y = kproj(2,j) * dl(2)
            z = kproj(3,k) * dl(3)
            r = SQRT(x**2+y**2+z**2)
            IF (r/=0.0_dp) THEN
                fw = f_ww(r,rmin2,rsw2,rmax2)
                r2 = r**2
                A0 (i,j,k) = fw
                Axx(i,j,k) = fw*x*x/r2
                Ayy(i,j,k) = fw*y*y/r2
                Azz(i,j,k) = fw*z*z/r2
                Axy(i,j,k) = fw*x*y/r2
                Axz(i,j,k) = fw*x*z/r2
                Ayz(i,j,k) = fw*z*y/r2
                Ax (i,j,k) = fw*x/r
                Ay (i,j,k) = fw*y/r
                Az (i,j,k) = fw*z/r
           !ELSE zero in all tables, i.e., they keep their initial value
            END IF
        END DO
    END BLOCK

!!!!!!!

    fftw3%in_forward=Axx
    CALL dfftw_execute(fftw3%plan_forward)
    Axx_k=fftw3%out_forward*deltaV
    fftw3%in_backward=Axx_k*rho_k
    CALL dfftw_execute(fftw3%plan_backward)
    ALLOCATE(Fxx(nfft1,nfft2,nfft3), SOURCE=fftw3%out_backward*deltaVk/(twopi)**3)
    fftw3%in_backward=(rho_k-function_rho_0k)*Axx_k
    CALL dfftw_execute(fftw3%plan_backward)
    ALLOCATE(FAxx(nfft1,nfft2,nfft3), SOURCE=fftw3%out_backward*deltaVk/(twopi)**3)

    fftw3%in_forward=Ayy
    CALL dfftw_execute(fftw3%plan_forward)
    Ayy_k=fftw3%out_forward*deltaV
    fftw3%in_backward=Ayy_k*rho_k
    CALL dfftw_execute(fftw3%plan_backward)
    ALLOCATE(Fyy(nfft1,nfft2,nfft3), SOURCE=fftw3%out_backward*deltaVk/(twopi)**3)
    fftw3%in_backward=(rho_k-function_rho_0k)*Ayy_k
    CALL dfftw_execute(fftw3%plan_backward)
    ALLOCATE(FAyy(nfft1,nfft2,nfft3), SOURCE=fftw3%out_backward*deltaVk/(twopi)**3)
    
    fftw3%in_forward=Azz
    CALL dfftw_execute(fftw3%plan_forward)
    Azz_k=fftw3%out_forward*deltaV
    fftw3%in_backward=Azz_k*rho_k
    CALL dfftw_execute(fftw3%plan_backward)
    ALLOCATE(Fzz(nfft1,nfft2,nfft3), SOURCE=fftw3%out_backward*deltaVk/(twopi)**3)
    fftw3%in_backward=(rho_k-function_rho_0k)*Azz_k
    CALL dfftw_execute(fftw3%plan_backward)
    ALLOCATE(FAzz(nfft1,nfft2,nfft3), SOURCE=fftw3%out_backward*deltaVk/(twopi)**3)
    
    fftw3%in_forward=Axy
    CALL dfftw_execute(fftw3%plan_forward)
    Axy_k=fftw3%out_forward*deltaV
    fftw3%in_backward=Axy_k*rho_k
    CALL dfftw_execute(fftw3%plan_backward)
    ALLOCATE(Fxy(nfft1,nfft2,nfft3), SOURCE=fftw3%out_backward*deltaVk/(twopi)**3)
    fftw3%in_backward=(rho_k-function_rho_0k)*Axy_k
    CALL dfftw_execute(fftw3%plan_backward)
    ALLOCATE(FAxy(nfft1,nfft2,nfft3), SOURCE=fftw3%out_backward*deltaVk/(twopi)**3)

    fftw3%in_forward=Axz
    CALL dfftw_execute(fftw3%plan_forward)
    Axz_k=fftw3%out_forward*deltaV
    fftw3%in_backward=Axz_k*rho_k
    CALL dfftw_execute(fftw3%plan_backward)
    ALLOCATE(Fxz(nfft1,nfft2,nfft3), SOURCE=fftw3%out_backward*deltaVk/(twopi)**3)
    fftw3%in_backward=(rho_k-function_rho_0k)*Axz_k
    CALL dfftw_execute(fftw3%plan_backward)
    ALLOCATE(FAxz(nfft1,nfft2,nfft3), SOURCE=fftw3%out_backward*deltaVk/(twopi)**3)
    
    fftw3%in_forward=Ayz
    CALL dfftw_execute(fftw3%plan_forward)
    Ayz_k=fftw3%out_forward*deltaV
    fftw3%in_backward=Ayz_k*rho_k
    CALL dfftw_execute(fftw3%plan_backward)
    ALLOCATE(Fyz(nfft1,nfft2,nfft3), SOURCE=fftw3%out_backward*deltaVk/(twopi)**3)
    fftw3%in_backward=(rho_k-function_rho_0k)*Ayz_k
    CALL dfftw_execute(fftw3%plan_backward)
    ALLOCATE(FAyz(nfft1,nfft2,nfft3), SOURCE=fftw3%out_backward*deltaVk/(twopi)**3)

    fftw3%in_forward=Ax
    CALL dfftw_execute(fftw3%plan_forward)
    Ax_k=fftw3%out_forward*deltaV
    fftw3%in_backward=Ax_k*rho_k
    CALL dfftw_execute(fftw3%plan_backward)
    ALLOCATE( Fx(nfft1,nfft2,nfft3), SOURCE=fftw3%out_backward*deltaVk/(twopi)**3)
    fftw3%in_backward=(rho_k-function_rho_0k)*Ax_k
    CALL dfftw_execute(fftw3%plan_backward)
    ALLOCATE(FAx(nfft1,nfft2,nfft3), SOURCE=fftw3%out_backward*deltaVk/(twopi)**3)
    
    fftw3%in_forward=Ay
    CALL dfftw_execute(fftw3%plan_forward)
    Ay_k=fftw3%out_forward*deltaV
    fftw3%in_backward=Ay_k*rho_k
    CALL dfftw_execute(fftw3%plan_backward)
    ALLOCATE( Fy(nfft1,nfft2,nfft3), SOURCE=fftw3%out_backward*deltaVk/(twopi)**3)
    fftw3%in_backward=(rho_k-function_rho_0k)*Ay_k
    CALL dfftw_execute(fftw3%plan_backward)
    ALLOCATE(FAy(nfft1,nfft2,nfft3), SOURCE=fftw3%out_backward*deltaVk/(twopi)**3)
    
    fftw3%in_forward=Az
    CALL dfftw_execute(fftw3%plan_forward)
    Az_k=fftw3%out_forward*deltaV
    fftw3%in_backward=Az_k*rho_k
    CALL dfftw_execute(fftw3%plan_backward)
    ALLOCATE( Fz(nfft1,nfft2,nfft3), SOURCE=fftw3%out_backward*deltaVk/(twopi)**3)
    fftw3%in_backward=(rho_k-function_rho_0k)*Az_k
    CALL dfftw_execute(fftw3%plan_backward)
    ALLOCATE(FAz(nfft1,nfft2,nfft3), SOURCE=fftw3%out_backward*deltaVk/(twopi)**3)
    
    fftw3%in_forward=A0
    CALL dfftw_execute(fftw3%plan_forward)
    A0_k=fftw3%out_forward*deltaV
    fftw3%in_backward=A0_k*rho_k
    CALL dfftw_execute(fftw3%plan_backward)
    ALLOCATE( F0(nfft1,nfft2,nfft3), SOURCE=fftw3%out_backward*deltaVk/(twopi)**3)
    fftw3%in_backward=(rho_k-function_rho_0k)*A0_k
    CALL dfftw_execute(fftw3%plan_backward)
    ALLOCATE(FA0(nfft1,nfft2,nfft3), SOURCE=fftw3%out_backward*deltaVk/(twopi)**3)


    
!!!!!!!!!

    ALLOCATE ( DHxx (nfft1,nfft2,nfft3,nb_solute_sites) ,SOURCE=0._dp)
    ALLOCATE ( DHyy (nfft1,nfft2,nfft3,nb_solute_sites) ,SOURCE=0._dp)
    ALLOCATE ( DHzz (nfft1,nfft2,nfft3,nb_solute_sites) ,SOURCE=0._dp)
    ALLOCATE ( DHxy (nfft1,nfft2,nfft3,nb_solute_sites) ,SOURCE=0._dp)
    ALLOCATE ( DHxz (nfft1,nfft2,nfft3,nb_solute_sites) ,SOURCE=0._dp)
    ALLOCATE ( DHyz (nfft1,nfft2,nfft3,nb_solute_sites) ,SOURCE=0._dp)
    ALLOCATE ( DHx  (nfft1,nfft2,nfft3,nb_solute_sites) ,SOURCE=0._dp)
    ALLOCATE ( DHy  (nfft1,nfft2,nfft3,nb_solute_sites) ,SOURCE=0._dp)
    ALLOCATE ( DHz  (nfft1,nfft2,nfft3,nb_solute_sites) ,SOURCE=0._dp)
    ALLOCATE ( DH0  (nfft1,nfft2,nfft3,nb_solute_sites) ,SOURCE=0._dp)

    DO CONCURRENT (n=1:nb_solute_sites, soluteSite(n)%lambda1/=0._dp)
    
        rmax1=0.5_dp*(soluteSite(n)%sig + solventSite(1)%sig) + d_w 
        nmax1x = int(rmax1/deltax)
        nmax1y = int(rmax1/deltay)
        nmax1z = int(rmax1/deltaz)
        ix = int(soluteSite(n)%r(1)/deltax) + 1
        iy = int(soluteSite(n)%r(2)/deltay) + 1
        iz = int(soluteSite(n)%r(3)/deltaz) + 1

        DO k=iz-nmax1z,iz+nmax1z
            zk2 = -(soluteSite(n)%r(3)-(k-1)*deltaz)
            DO j=iy-nmax1y, iy+nmax1y
                yk2 = -(soluteSite(n)%r(2)-(j-1)*deltay)
                DO i=ix-nmax1x, ix+nmax1x
                    xk2 = -(soluteSite(n)%r(1)-(i-1)*deltax)
                    rk2 = SQRT( xk2**2 + yk2**2 + zk2**2 )
                    fk1= deltaV*f_ww(rk2,rmin1,rsw1,rmax1)
                    IF (rk2 /= 0.0_dp .AND. fk1/=0._dp) THEN
                        rho_temp=rho(i,j,k)
                        DHxx(i,j,k,n)= fk1*xk2**2/(rk2**2)
                        DHyy(i,j,k,n)= fk1*yk2**2/(rk2**2)
                        DHzz(i,j,k,n)= fk1*zk2**2/(rk2**2)
                        DHxy(i,j,k,n)= fk1*xk2*yk2/(rk2**2)
                        DHxz(i,j,k,n)= fk1*xk2*zk2/(rk2**2)
                        DHyz(i,j,k,n)= fk1*yk2*zk2/(rk2**2)
                        DHx(i,j,k,n) = fk1*xk2/rk2
                        DHy(i,j,k,n) = fk1*yk2/rk2
                        DHz(i,j,k,n) = fk1*zk2/rk2
                        DH0(i,j,k,n) = fk1
                    END IF
                END DO
            END DO
        END DO
    END DO

    ALLOCATE ( Hxx (nb_solute_sites) ,SOURCE=0._dp)
    ALLOCATE ( Hyy (nb_solute_sites) ,SOURCE=0._dp)
    ALLOCATE ( Hzz (nb_solute_sites) ,SOURCE=0._dp)
    ALLOCATE ( Hxy (nb_solute_sites) ,SOURCE=0._dp)
    ALLOCATE ( Hxz (nb_solute_sites) ,SOURCE=0._dp)
    ALLOCATE ( Hyz (nb_solute_sites) ,SOURCE=0._dp)
    ALLOCATE ( Hx  (nb_solute_sites) ,SOURCE=0._dp)
    ALLOCATE ( Hy  (nb_solute_sites) ,SOURCE=0._dp)
    ALLOCATE ( Hz  (nb_solute_sites) ,SOURCE=0._dp)
    ALLOCATE ( H0  (nb_solute_sites) ,SOURCE=0._dp) ! USE ARRAY CONSTRUCTOR INSTEAD
    DO CONCURRENT (n=1:nb_solute_sites)
        Hxx(n)= SUM(  DHxx(:,:,:,n)*rho )
        Hyy(n)= SUM(  DHyy(:,:,:,n)*rho )
        Hzz(n)= SUM(  DHzz(:,:,:,n)*rho )
        Hxy(n)= SUM(  DHxy(:,:,:,n)*rho )
        Hxz(n)= SUM(  DHxz(:,:,:,n)*rho )
        Hyz(n)= SUM(  DHyz(:,:,:,n)*rho )
        Hx(n) = SUM(  DHx(:,:,:,n) *rho )
        Hy(n) = SUM(  DHy(:,:,:,n) *rho )
        Hz(n) = SUM(  DHz(:,:,:,n) *rho )
        H0(n) = SUM(  DH0(:,:,:,n) *rho )
    END DO


    ALLOCATE(Gxx(nfft1,nfft2,nfft3), SOURCE=0.0_dp)
    ALLOCATE(Gyy(nfft1,nfft2,nfft3), SOURCE=0.0_dp)
    ALLOCATE(Gzz(nfft1,nfft2,nfft3), SOURCE=0.0_dp)
    ALLOCATE(Gxy(nfft1,nfft2,nfft3), SOURCE=0.0_dp)
    ALLOCATE(Gxz(nfft1,nfft2,nfft3), SOURCE=0.0_dp)
    ALLOCATE(Gyz(nfft1,nfft2,nfft3), SOURCE=0.0_dp)
    ALLOCATE( Gx(nfft1,nfft2,nfft3), SOURCE=0.0_dp)
    ALLOCATE( Gy(nfft1,nfft2,nfft3), SOURCE=0.0_dp)
    ALLOCATE( Gz(nfft1,nfft2,nfft3), SOURCE=0.0_dp)
    ALLOCATE( G0(nfft1,nfft2,nfft3), SOURCE=0.0_dp)

    
    DO CONCURRENT ( n=1:nb_solute_sites , soluteSite(n)%lambda2/=0._dp )
        nmax2x = int(rmax1/deltax)
        nmax2y = int(rmax1/deltay)
        nmax2z = int(rmax1/deltaz)
        ix = int(soluteSite(n)%r(1)/deltax) + 1
        iy = int(soluteSite(n)%r(2)/deltay) + 1
        iz = int(soluteSite(n)%r(3)/deltaz) + 1
        DO k=iz-nmax2z,iz+nmax2z+1
            zk2 = -(soluteSite(n)%r(3)-(k-1)*deltaz)
            DO j=iy-nmax2y, iy+nmax2y+1
                yk2 = -(soluteSite(n)%r(2)-(j-1)*deltay)
                DO i=ix-nmax2x, ix+nmax2x+1
                    xk2 = -(soluteSite(n)%r(1)-(i-1)*deltax)
                    rk2 = SQRT( xk2**2 + yk2**2 + zk2**2 )
                    fk2 = soluteSite(n)%lambda2 * f_ww(rk2,rmin1,rsw1,rmax1)
                    G0(i,j,k) = G0(i,j,k) + fk2                    
                    IF (rk2/=0.0_dp .AND. fk2/=0._dp) THEN
                        Gxx(i,j,k)= Gxx(i,j,k)+ fk2*xk2**2 /(rk2**2)
                        Gyy(i,j,k)= Gyy(i,j,k)+ fk2*yk2**2 /(rk2**2)
                        Gzz(i,j,k)= Gzz(i,j,k)+ fk2*zk2**2 /(rk2**2)
                        Gxy(i,j,k)= Gxy(i,j,k)+ fk2*xk2*yk2/(rk2**2)
                        Gxz(i,j,k)= Gxz(i,j,k)+ fk2*xk2*zk2/(rk2**2)
                        Gyz(i,j,k)= Gyz(i,j,k)+ fk2*yk2*zk2/(rk2**2)
                        Gx (i,j,k)= Gx (i,j,k)+ fk2*xk2/rk2
                        Gy (i,j,k)= Gy (i,j,k)+ fk2*yk2/rk2
                        Gz (i,j,k)= Gz (i,j,k)+ fk2*zk2/rk2
                    END IF
                END DO
            END DO
        END DO
    END DO


    fftw3%in_forward=Gxx*rho
    CALL dfftw_execute(fftw3%plan_forward)
    ALLOCATE(Gxx_k(nfft1/2+1,nfft2,nfft3), SOURCE=fftw3%out_forward*deltaV)
    CALL dothestuff( FGxx , Axx_k , Gxx_k )
        
    fftw3%in_forward=Gyy*rho
    CALL dfftw_execute(fftw3%plan_forward)
    ALLOCATE(Gyy_k(nfft1/2+1,nfft2,nfft3), SOURCE=fftw3%out_forward*deltaV)
    CALL dothestuff( FGyy , Ayy_k , Gyy_k )
    
    fftw3%in_forward=Gzz*rho
    CALL dfftw_execute(fftw3%plan_forward)
    ALLOCATE(Gzz_k(nfft1/2+1,nfft2,nfft3), SOURCE=fftw3%out_forward*deltaV)
    CALL dothestuff( FGzz , Azz_k , Gzz_k )
    
    fftw3%in_forward=Gxy*rho
    CALL dfftw_execute(fftw3%plan_forward)
    ALLOCATE(Gxy_k(nfft1/2+1,nfft2,nfft3), SOURCE=fftw3%out_forward*deltaV)
    CALL dothestuff( FGxy , Axy_k , Gxy_k )
    
    fftw3%in_forward=Gxz*rho
    CALL dfftw_execute(fftw3%plan_forward)
    ALLOCATE(Gxz_k(nfft1/2+1,nfft2,nfft3), SOURCE=fftw3%out_forward*deltaV)
    CALL dothestuff( FGxz , Axz_k , Gxz_k )
    
    fftw3%in_forward=Gyz*rho
    CALL dfftw_execute(fftw3%plan_forward)
    ALLOCATE(Gyz_k(nfft1/2+1,nfft2,nfft3), SOURCE=fftw3%out_forward*deltaV)
    CALL dothestuff( FGyz , Ayz_k , Gyz_k )
    
    fftw3%in_forward=Gx*rho
    CALL dfftw_execute(fftw3%plan_forward)
    ALLOCATE( Gx_k(nfft1/2+1,nfft2,nfft3), SOURCE=fftw3%out_forward*deltaV)
    CALL dothestuff( FGx  ,  Ax_k , Gx_k  )
    
    fftw3%in_forward=Gy*rho
    CALL dfftw_execute(fftw3%plan_forward)
    ALLOCATE( Gy_k(nfft1/2+1,nfft2,nfft3), SOURCE=fftw3%out_forward*deltaV)
    CALL dothestuff( FGy  ,  Ay_k , Gy_k  )
    
    fftw3%in_forward=Gz*rho
    CALL dfftw_execute(fftw3%plan_forward)
    ALLOCATE( Gz_k(nfft1/2+1,nfft2,nfft3), SOURCE=fftw3%out_forward*deltaV)
    CALL dothestuff( FGz  ,  Az_k , Gz_k  )
    
    fftw3%in_forward=G0*rho
    CALL dfftw_execute(fftw3%plan_forward)
    ALLOCATE( G0_k(nfft1/2+1,nfft2,nfft3), SOURCE=fftw3%out_forward*deltaV)
    CALL dothestuff( FG0  ,  A0_k , G0_k  )


    F3B1 = thermocond%kbT/2._dp* SUM (soluteSite%lambda1*( Hxx**2+Hyy**2+Hzz**2 +2.0_dp*Hxy**2+2.0_dp*Hxz**2+2.0_dp*Hyz**2&
                                          -2.0_dp*costheta0*(Hx**2+Hy**2+Hz**2)+costheta0**2*H0**2   ))
    

    F3B2 = thermocond%kbT/2._dp*deltaV*SUM(rho*( (Fxx*Gxx+Fyy*Gyy+Fzz*Gzz+2._dp*Fxy*Gxy+2._dp*Fxz*Gxz+2._dp*Fyz*Gyz)   &
                                       -2._dp*costheta0*(Fx*Gx+Fy*Gy+Fz*Gz)      &
                                       +costheta0**2*F0*G0   ))


    CALL compute_water_water_3body_term (F3B_ww)

BLOCK
    REAL(dp) :: opweight
    icg=0
    DO i=1,nfft1
        DO j=1, nfft2
            DO k=1, nfft3
                DO o=1,angGrid%n_angles
                    DO p=1,molRotGrid%n_angles
                        
                        opweight = angGrid%weight(o)*molRotGrid%weight(p)
                        
                        icg=icg+1    
                        psi=cg_vect(icg)
                        
                        dF(icg)=dF(icg)+thermocond%kbT*psi*deltaV*opweight*solvent(1)%rho0*(&
                      (Fxx(i,j,k)*Gxx(i,j,k)+Fyy(i,j,k)*Gyy(i,j,k)+Fzz(i,j,k)*Gzz(i,j,k)&
                      +2.0_dp*Fxy(i,j,k)*Gxy(i,j,k)+ 2.0_dp*Fxz(i,j,k)*Gxz(i,j,k)+ 2.0_dp*Fyz(i,j,k)*Gyz(i,j,k))&
                      -2.0_dp*costheta0*(Fx(i,j,k)*Gx(i,j,k)+Fy(i,j,k)*Gy(i,j,k)+Fz(i,j,k)*Gz(i,j,k))&
                      +costheta0**2*F0(i,j,k)*G0(i,j,k))+&
                      thermocond%kbT*psi*opweight*&
                      solvent(1)%rho0*deltaV*(FGxx(i,j,k)+FGyy(i,j,k)+FGzz(i,j,k)+2.0_dp*FGxy(i,j,k)&
                     +2.0_dp*FGxz(i,j,k)+2.0_dp*FGyz(i,j,k)+2.0_dp*costheta0*(FGx(i,j,k)+FGy(i,j,k)+FGz(i,j,k))+&
                     costheta0**2*FG0(i,j,k))
             
              dF(icg)=dF(icg)+thermocond%kbT*deltaV*opweight*solvent(1)%rho0*lambda_w*psi*((FAxx(i,j,k)**2+FAyy(i,j,k)**2+&
              FAzz(i,j,k)**2+2.0_dp*(FAxy(i,j,k)**2+FAxz(i,j,k)**2+FAyz(i,j,k)**2)-2.0_dp*costheta0*(FAx(i,j,k)**2+FAy(i,j,k)**2&
              +FAz(i,j,k)**2)+costheta0**2*FA0(i,j,k)**2)&
                      +2.0_dp*(FAxx(i,j,k)*Axx(i,j,k)+FAyy(i,j,k)*Ayy(i,j,k)+FAzz(i,j,k)*Azz(i,j,k)+2.0_dp*FAxy(i,j,k)*Axy(i,j,k)+&
                      2.0_dp*FAyz(i,j,k)*Ayz(i,j,k)+2.0_dp*FAxz(i,j,k)*Axz(i,j,k)&
              -2.0_dp*costheta0*(FAx(i,j,k)*Ax(i,j,k)+FAy(i,j,k)*Ay(i,j,k)+FAz(i,j,k)*Az(i,j,k))+costheta0**2*FA0(i,j,k)*A0(i,j,k)))
                      
                    dF(icg)=dF(icg) +SUM(soluteSite%lambda1*thermocond%kbT*psi*opweight*solvent(1)%rho0*2.0_dp*(&
                        (Hxx(:)*DHxx(i,j,k,:)+Hyy(:)*DHyy(i,j,k,:)+Hzz(:)*DHzz(i,j,k,:)&
                        +2.0_dp*Hxy(:)*DHxy(i,j,k,:)+2.0_dp*Hxz(:)*DHxz(i,j,k,:)+2.0_dp*Hyz(:)*DHyz(i,j,k,:))&
                        -2.0_dp*costheta0*(Hx(:)*DHx(i,j,k,:)+Hy(:)*DHy(i,j,k,:)+Hz(:)*DHz(i,j,k,:))&
                        +costheta0**2*H0(:)*DH0(i,j,k,:)))
!~                         DO n=1,nb_solute_sites
!~                             dF(icg)=dF(icg)+soluteSite(n)%lambda1*thermocond%kbT*psi*opweight*((Hxx(n)*DHxx(i,j,k,n)+Hyy(n)*DHyy(i,j,k,n)+&
!~                                 Hzz(n)*DHzz(i,j,k,n)&
!~                                 +2.0_dp*Hxy(n)*DHxy(i,j,k,n)+2.0_dp*Hxz(n)*DHxz(i,j,k,n)+2.0_dp*Hyz(n)*DHyz(i,j,k,n))&
!~           -2.0_dp*costheta0*(Hx(n)*DHx(i,j,k,n)+Hy(n)*DHy(i,j,k,n)+Hz(n)*DHz(i,j,k,n))+costheta0**2*H0(n)*DH0(i,j,k,n))*solvent(1)%rho0*2.0_dp
!~                         END DO
             
                    END DO
                END DO
            END DO
        END DO
    END DO
END BLOCK

    FF = FF + F3B2 + F3B1 + F3B_ww

    CALL CPU_TIME (time1)
    
    IF (verbose) THEN
        WRITE(*,'(''    F3B1               = '',f11.3,'' in '',I5,'' sec'')') F3B1 , NINT(time1-time0)
        WRITE(*,'(''    F3B2               = '',f11.3,'' in '',I5,'' sec'')') F3B2 , NINT(time1-time0)
    END IF

    DEALLOCATE(FAxx)
    DEALLOCATE(FAyy)
    DEALLOCATE(FAzz)
    DEALLOCATE(FAxy)
    DEALLOCATE(FAxz)
    DEALLOCATE(FAyz)
    DEALLOCATE( FAx)
    DEALLOCATE( FAy)
    DEALLOCATE( FAz)
    DEALLOCATE( FA0)
    DEALLOCATE(Axx)
    DEALLOCATE(Ayy)
    DEALLOCATE(Azz)
    DEALLOCATE(Axy)
    DEALLOCATE(Axz)
    DEALLOCATE(Ayz)
    DEALLOCATE( Ax)
    DEALLOCATE( Ay)
    DEALLOCATE( Az)
    DEALLOCATE( A0)
    DEALLOCATE(Gxx)
    DEALLOCATE(Gyy)
    DEALLOCATE(Gzz)
    DEALLOCATE(Gxy)
    DEALLOCATE(Gxz)
    DEALLOCATE(Gyz)
    DEALLOCATE( Gx)
    DEALLOCATE( Gy)
    DEALLOCATE( Gz)
    DEALLOCATE( G0)


    CONTAINS
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    SUBROUTINE dothestuff (A,B,C)
        USE system  ,ONLY:  spaceGrid
        INTEGER(i2b)    ::  i,j,k
        REAL(dp), ALLOCATABLE, INTENT(INOUT)    :: A(:,:,:)
        COMPLEX(dp), ALLOCATABLE, INTENT(INOUT) :: B(:,:,:), C(:,:,:)
        i=spaceGrid%n_nodes(1)
        j=spaceGrid%n_nodes(2)
        k=spaceGrid%n_nodes(3)
        ALLOCATE( A(i,j,k) ,SOURCE=0.0_dp )
        fftw3%in_backward=B*C
        DEALLOCATE(B)
        DEALLOCATE(C)
        CALL dfftw_execute ( fftw3%plan_backward )
        A=fftw3%out_backward*deltaVk/(twopi)**3
    END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    PURE FUNCTION f_ww (r,rmin,rsw,rmax)
        USE precision_kinds, only: dp,i2b
        IMPLICIT NONE
        REAL(dp)             :: f_ww
        REAL(dp), INTENT(IN) :: r, rmin, rsw, rmax
        REAL(dp), PARAMETER  :: gam = 2.0_dp/3.0_dp
        REAL(dp)             :: deltar, exp_term, Switch
        deltar = rsw-rmin
        IF (r > rmin .and. r<rmax) THEN
            exp_term = exp(gam*rmax/(r-rmax))
            IF (r<rsw) THEN
                Switch = (r-rmin)**2*(-2.0_dp*(r-rmin)/deltar+3.0_dp)/deltar**2
            ELSE
                Switch = 1.0_dp
            END IF
            f_ww = Switch*exp_term
        ELSE
            f_ww = 0.0_dp
        END IF
    END FUNCTION f_ww

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! water water 3 body term
    PURE SUBROUTINE compute_water_water_3body_term (F3B_ww)
        REAL(dp), INTENT(OUT) :: F3B_ww
        IF (lambda_w/=0._dp) THEN
            F3B_ww= 0.5_dp*thermocond%kbT*deltaV*lambda_w*SUM(rho*(FAxx**2+FAyy**2+FAzz**2+2._dp*(FAxy**2+FAxz**2+FAyz**2)&
                -2._dp*costheta0*(FAx**2+FAy**2+FAz**2) +costheta0**2*FA0**2))
        ELSE
            F3B_ww=0._dp
        END IF
    END SUBROUTINE compute_water_water_3body_term

END SUBROUTINE energy_threebody_faster




    
    
    
    


!     Last change:  LB   27 May 2015    5:43 pm
!     Last change:  LB   13 Jan 2012    3:13 pm
subroutine accumu_mnmunukhi
implicit REAL(8) (a-h,o-z)
REAL * 8 iacc,iacc_o
common/nombre/n_atomes
common/positi/xx(10000),yy(10000),zz(10000)        ! positions cellule centrale
common/numeri/p,num
COMMON/nb_proj/mnmax,ialpmax
PARAMETER(ialpmaxx=1226)
PARAMETER(n_omega=8)
COMMON/proj_l/mm(ialpmaxx),nn(ialpmaxx),ll(ialpmaxx),mumu(ialpmaxx),nunu(ialpmaxx)
COMMON/proj_khi/mm1(ialpmaxx),nn1(ialpmaxx),khi1(ialpmaxx),mumu1(ialpmaxx),nunu1(ialpmaxx)
common/cumul/iacc(ialpmaxx,2000),gr(ialpmaxx,2000),grint(2000)
common/cumul_omega/iacc_o(n_omega,2000),gr_o(n_omega,2000),dcosbet,dphi,coeff_a,coeff_b           ! pour g(r,omega) luc85p176,190
common/diamet/sigma_a,sigma2,rcutoff_a,rcut2,roh_a,theta_d
common/conc/conc_h2o,xl_a
COMMON/effica/aaa,ichoixaaa,aaa_bis(80)                         ! differentes aaa dans aaa_bis luc82p82
DIMENSION k_bis(10)
COMMON/densite_a_2_corps/i_rho2
DIMENSION rac(-1:2*mnmax+1)
REAL(8), DIMENSION(0:mnmax,-mnmax:mnmax,-mnmax:mnmax):: a,b,c,d,fi,gi,fj,gj
REAL(8),DIMENSION(3,3):: ri,rj,rij,rmati,rmatj
!
interface
    FUNCTION vector(r1,r2)
        implicit REAL(8) (a-h,o-z)
        REAL(8), DIMENSION(1:3):: r1,r2,vector
    end function
END interface
!
!
!
!        projections gmnmunukhi(r) accumulees dans iacc  pour H2O
!        luc85p23 et Choi et al. J. Chem. Phys. 111, 8825 (1999)
!
!        version light pour illustration pour Lu, Maximilien, Danielluc89p178
!
cosdphi2=COS(dphi/2.d0)
sindphi2=sin(dphi/2.d0)
costetra=1.d0/SQRT(3.d0)
!
n_h2o=n_atomes/3
roh=roh_a/xl_a
pi=4.d0*ATAN(1.d0)
theta=pi/180.d0*theta_d
theta2=theta/2.d0
cthet2=COS(theta2); sthet2=SIN(theta2)
rhh=2.d0*roh*sthet2
rohh=2.d0*roh*cthet2
rac2=sqrt(2.d0)
rac(-1)=0.                    ! pour eviter des pbs
do k=0,2*mnmax+1
rac(k)=SQRT(dble(k))                           ! racine(k) en tableau
end do
!
!                                 coeff pour recurrence
!
do l=1,mnmax
 do m=-l,l
  do m1=-l+1,l-1
  a(l,m,m1)=rac(l+m)*rac(l-m)/(rac(l+m1)*rac(l-m1))              ! coeff a,b,c,d lmm'
  b(l,m,m1)=rac(l+m)*rac(l+m-1)/(rac(2)*rac(l+m1)*rac(l-m1))
  end do
 m1=l
 c(l,m,m1)=rac(2)*rac(l+m)*rac(l-m)/(rac(l+m1)*rac(l+m1-1))
 d(l,m,m1)=rac(l+m)*rac(l+m-1)/(rac(l+m1)*rac(l+m1-1))
 end do
end do


!
!               en r, je construis les gmnmunukhi(r)
!
    do i_h2o=1,n_h2o             ! molecule i
!
        io=3*(i_h2o-1)+1
        ih1=io+1
        ih2=io+2
!          construire les coord des vecteurs lies a i_h2o dans repere fixe

        ri(3,1)=(xx(ih1)+xx(ih2)-2.d0*xx(io))/rohh        ! axe principal z suivant O - H1+H2
        ri(3,2)=(yy(ih1)+yy(ih2)-2.d0*yy(io))/rohh
        ri(3,3)=(zz(ih1)+zz(ih2)-2.d0*zz(io))/rohh
        ri(1,1)=(xx(ih2)-xx(ih1))/rhh                      ! axe x suivant H1H2
        ri(1,2)=(yy(ih2)-yy(ih1))/rhh
        ri(1,3)=(zz(ih2)-zz(ih1))/rhh
        ri(2,1:3)=vector(ri(3,1:3),ri(1,1:3))                  ! 3eme axe, perpendiculaire aux premiers
!
      do j_h2o=i_h2o+1,n_h2o
!
        jo=3*(j_h2o-1)+1
        jh1=jo+1
        jh2=jo+2
        xij=xx(jo)-xx(io)                               ! vecteur Oi - Oj
        yij=yy(jo)-yy(io)
        zij=zz(jo)-zz(io)
        if(xij.gt.0.5) xij=xij-1.
        if(xij.lt.-0.5) xij=xij+1.
        if(yij.gt.0.5) yij=yij-1.
        if(yij.lt.-0.5) yij=yij+1.
        if(zij.gt.0.5) zij=zij-1.
        if(zij.lt.-0.5) zij=zij+1.
         r2=xij**2+yij**2+zij**2
         IF(r2.gt.rmax2_p) cycle
         r=sqrt(r2)
         kr=int(r/p)+1
!
!          construire les coord des vecteurs du repere intermoleculaire ij dans repere fixe
!
!           pour Lu...
!
!        voila comment je construis la matrice de rotation qui fait passer du repere fixe a un repere lie a r12
!        pour un repere lie a q, c'est pareil!
!
         rij(3,1)=xij/r; rij(3,2)=yij/r; rij(3,3)=zij/r                ! vecteur unitaire rij
         IF(ABS(rij(3,1))<0.99) then                                ! puis x*rij (sauf si rij suivant x)
         rij(1,1:3)=vector((/1.d0,0.d0,0.d0/),rij(3,1:3))
         else                                              ! dans ce cas y*rij
         rij(1,1:3)=vector((/0.d0,1.d0,0.d0/),rij(3,1:3))
         endif
         rij1_norm=SQRT(SUM(rij(1,:)**2))
         rij(1,:)=rij(1,:)/rij1_norm
         rij(2,1:3)=vector(rij(3,1:3),rij(1,1:3))
!          construire les coord des vecteurs lies a j_h2o dans repere fixe
        rj(3,1)=(xx(jh1)+xx(jh2)-2.d0*xx(jo))/rohh        ! axe principal z suivant O - H1+H2
        rj(3,2)=(yy(jh1)+yy(jh2)-2.d0*yy(jo))/rohh
        rj(3,3)=(zz(jh1)+zz(jh2)-2.d0*zz(jo))/rohh
        rj(1,1)=(xx(jh2)-xx(jh1))/rhh                      ! axe x suivant H1H2
        rj(1,2)=(yy(jh2)-yy(jh1))/rhh
        rj(1,3)=(zz(jh2)-zz(jh1))/rhh
        rj(2,1:3)=vector(rj(3,1:3),rj(1,1:3))               ! 3eme axe, perpendiculaire aux premiers
!
!         construire alors les matrices de Messiah Rlmm' pour i
!
fi(0,0,0)=1.d0; gi(0,0,0)=0.          ! facile!
!         coord des vecteurs de i dans le repere intermol
        do k=1,3
         do l=1,3
         rmati(k,l)=DOT_PRODUCT(rij(k,:),ri(l,:))
         end do
        end do
!
!         pour Lu...
!        j'ai donc construit la matrice rotation "rmati" qui fait passer du repere local a celui lie a la molecule i
!        et je m'apprete a en deduire les Rm a commencer par R1
!        Si l'on veut faire la meme chose pour la rotation qui fait passer du repere fixe au repere lie a r12 (ou q pour Lu...),
!        il suffit de faire la meme chose qui suit mais en partant de la matrice rotation "rij"!
!
!         matrice Messiah Rlmm' pour l=1  = Flmm' + i Glmm'
        fi(1,-1,-1)=(rmati(2,2)+rmati(1,1))/2.d0
        fi(1,-1,0)=rmati(1,3)/rac2
        fi(1,-1,+1)=(rmati(2,2)-rmati(1,1))/2.d0
        fi(1,0,-1)=rmati(3,1)/rac2
        fi(1,0,0)=rmati(3,3)
        fi(1,0,+1)=-rmati(3,1)/rac2
        fi(1,+1,-1)=(rmati(2,2)-rmati(1,1))/2.d0
        fi(1,+1,0)=-rmati(1,3)/rac2
        fi(1,+1,+1)=(rmati(2,2)+rmati(1,1))/2.d0
        gi(1,-1,-1)=(rmati(2,1)-rmati(1,2))/2.d0
        gi(1,-1,0)=rmati(2,3)/rac2
        gi(1,-1,+1)=(-rmati(2,1)-rmati(1,2))/2.d0
        gi(1,0,-1)=-rmati(3,2)/rac2
        gi(1,0,0)=0.
        gi(1,0,+1)=-rmati(3,2)/rac2
        gi(1,+1,-1)=(rmati(2,1)+rmati(1,2))/2.d0
        gi(1,+1,0)=rmati(2,3)/rac2
        gi(1,+1,+1)=(rmati(1,2)-rmati(2,1))/2.d0
!        matrice Messiah Rlmm' pour l>1 par recurrence Choi... luc85p17
do l=2,mnmax
l1=l-1
 do m=-l,l
 m1min=0; IF(m>0) m1min=1
  do m1=m1min,l-1
  fi(l,m,m1)=a(l,m,m1)*(fi(1,0,0)*fi(l1,m,m1)-gi(1,0,0)*gi(l1,m,m1))+         &
            b(l,m,m1)*(fi(1,+1,0)*fi(l1,m-1,m1)-gi(1,+1,0)*gi(l1,m-1,m1))+   &
            b(l,-m,m1)*(fi(1,-1,0)*fi(l1,m+1,m1)-gi(1,-1,0)*gi(l1,m+1,m1))
  gi(l,m,m1)=a(l,m,m1)*(fi(1,0,0)*gi(l1,m,m1)+gi(1,0,0)*fi(l1,m,m1))+         &
            b(l,m,m1)*(fi(1,+1,0)*gi(l1,m-1,m1)+gi(1,+1,0)*fi(l1,m-1,m1))+   &
            b(l,-m,m1)*(fi(1,-1,0)*gi(l1,m+1,m1)+gi(1,-1,0)*fi(l1,m+1,m1))
  fi(l,-m,-m1)=(-1)**(m+m1)*fi(l,m,m1)
  gi(l,-m,-m1)=-(-1)**(m+m1)*gi(l,m,m1)
  end do
 m1=l
 fi(l,m,m1)=c(l,m,m1)*(fi(1,0,+1)*fi(l1,m,m1-1)-gi(1,0,+1)*gi(l1,m,m1-1))+         &
           d(l,m,m1)*(fi(1,+1,+1)*fi(l1,m-1,m1-1)-gi(1,+1,+1)*gi(l1,m-1,m1-1))+   &
           d(l,-m,m1)*(fi(1,-1,+1)*fi(l1,m+1,m1-1)-gi(1,-1,+1)*gi(l1,m+1,m1-1))
 gi(l,m,m1)=c(l,m,m1)*(fi(1,0,+1)*gi(l1,m,m1-1)+gi(1,0,+1)*fi(l1,m,m1-1))+         &
           d(l,m,m1)*(fi(1,+1,+1)*gi(l1,m-1,m1-1)+gi(1,+1,+1)*fi(l1,m-1,m1-1))+   &
           d(l,-m,m1)*(fi(1,-1,+1)*gi(l1,m+1,m1-1)+gi(1,-1,+1)*fi(l1,m+1,m1-1))
 fi(l,-m,-m1)=(-1)**(m+m1)*fi(l,m,m1)
 gi(l,-m,-m1)=-(-1)**(m+m1)*gi(l,m,m1)
 end do
end do
!
!           idem pour j
!
fj(0,0,0)=1.d0; gj(0,0,0)=0.          ! facile!
!         coord des vecteurs de j dans le repere intermol
        do k=1,3
         do l=1,3
         rmatj(k,l)=DOT_PRODUCT(rij(k,:),rj(l,:))
         end do
        end do
!PRINT*, 'j interm',((rmatj(k,l),l=1,3),k=1,3)
!         matrice Messiah Rlmm' pour l=1  = Flmm' + i Glmm'
        fj(1,-1,-1)=(rmatj(2,2)+rmatj(1,1))/2.d0
        fj(1,-1,0)=rmatj(1,3)/rac2
        fj(1,-1,+1)=(rmatj(2,2)-rmatj(1,1))/2.d0
        fj(1,0,-1)=rmatj(3,1)/rac2
        fj(1,0,0)=rmatj(3,3)
        fj(1,0,+1)=-rmatj(3,1)/rac2
        fj(1,+1,-1)=(rmatj(2,2)-rmatj(1,1))/2.d0
        fj(1,+1,0)=-rmatj(1,3)/rac2
        fj(1,+1,+1)=(rmatj(2,2)+rmatj(1,1))/2.d0
        gj(1,-1,-1)=(rmatj(2,1)-rmatj(1,2))/2.d0
        gj(1,-1,0)=rmatj(2,3)/rac2
        gj(1,-1,+1)=(-rmatj(2,1)-rmatj(1,2))/2.d0
        gj(1,0,-1)=-rmatj(3,2)/rac2
        gj(1,0,0)=0.
        gj(1,0,+1)=-rmatj(3,2)/rac2
        gj(1,+1,-1)=(rmatj(2,1)+rmatj(1,2))/2.d0
        gj(1,+1,0)=rmatj(2,3)/rac2
        gj(1,+1,+1)=(rmatj(1,2)-rmatj(2,1))/2.d0
!        matrice Messiah Rlmm' pour l>1 par recurrence Choi... luc85p17
do l=2,mnmax
l1=l-1
 do m=-l,l
 m1min=0; IF(m>0) m1min=1
  do m1=m1min,l-1
  fj(l,m,m1)=a(l,m,m1)*(fj(1,0,0)*fj(l1,m,m1)-gj(1,0,0)*gj(l1,m,m1))+         &
            b(l,m,m1)*(fj(1,+1,0)*fj(l1,m-1,m1)-gj(1,+1,0)*gj(l1,m-1,m1))+   &
            b(l,-m,m1)*(fj(1,-1,0)*fj(l1,m+1,m1)-gj(1,-1,0)*gj(l1,m+1,m1))
  gj(l,m,m1)=a(l,m,m1)*(fj(1,0,0)*gj(l1,m,m1)+gj(1,0,0)*fj(l1,m,m1))+         &
            b(l,m,m1)*(fj(1,+1,0)*gj(l1,m-1,m1)+gj(1,+1,0)*fj(l1,m-1,m1))+   &
            b(l,-m,m1)*(fj(1,-1,0)*gj(l1,m+1,m1)+gj(1,-1,0)*fj(l1,m+1,m1))
  fj(l,-m,-m1)=(-1)**(m+m1)*fj(l,m,m1)
  gj(l,-m,-m1)=-(-1)**(m+m1)*gj(l,m,m1)
  end do
 m1=l
 fj(l,m,m1)=c(l,m,m1)*(fj(1,0,+1)*fj(l1,m,m1-1)-gj(1,0,+1)*gj(l1,m,m1-1))+         &
           d(l,m,m1)*(fj(1,+1,+1)*fj(l1,m-1,m1-1)-gj(1,+1,+1)*gj(l1,m-1,m1-1))+   &
           d(l,-m,m1)*(fj(1,-1,+1)*fj(l1,m+1,m1-1)-gj(1,-1,+1)*gj(l1,m+1,m1-1))
 gj(l,m,m1)=c(l,m,m1)*(fj(1,0,+1)*gj(l1,m,m1-1)+gj(1,0,+1)*fj(l1,m,m1-1))+         &
           d(l,m,m1)*(fj(1,+1,+1)*gj(l1,m-1,m1-1)+gj(1,+1,+1)*fj(l1,m-1,m1-1))+   &
           d(l,-m,m1)*(fj(1,-1,+1)*gj(l1,m+1,m1-1)+gj(1,-1,+1)*fj(l1,m+1,m1-1))
 fj(l,-m,-m1)=(-1)**(m+m1)*fj(l,m,m1)
 gj(l,-m,-m1)=-(-1)**(m+m1)*gj(l,m,m1)
 end do
end do
!
!        construire alors Reelle(Rmkhimu(i)*Rn-khinu(j)) pour tous les ialp de la liste
!
!phitest1=0.
do ialp=1,ialpmax
m=mm1(ialp); n=nn1(ialp); mu=mumu1(ialp); nu=nunu1(ialp); khi=khi1(ialp)
phi=fi(m,khi,mu)*fj(n,-khi,nu)-gi(m,khi,mu)*gj(n,-khi,nu)
IF(m/=n.or.mu/=nu) phi=(phi+(-1)**(m+n)*(fi(n,khi,nu)*fj(m,-khi,mu)-gi(n,khi,nu)*gj(m,-khi,mu)))/2.d0
phi=rac(2*m+1)*rac(2*n+1)*phi
!
IF(i_rho2==2) phi=phi/(vol*grint(kr))              ! j'accumule eventuellement PHIij/(deltaV*V)  luc88p138
!
iacc(ialp,kr)=iacc(ialp,kr)+phi
IF(ialp==ialp_aaa.and.kr==k_aaa) aaa=aaa+phi
IF(ialp<=7) then
do kk=1,k_bis_max
IF(kr==k_bis(kk)) aaa_bis(10*ialp+kk)=aaa_bis(10*ialp+kk)+phi
end do
endif
END do
!
      END do                            ! fin j
!
    END do                            ! fin i
!
!PRINT*, 'phimnlmunu(ij) fute = ',phitest1
!                                              on normalise aaa pour que ce soit g
      IF(ichoixaaa>100) then
      den=n_h2o**2/2.*4.*pi/3.*((k_aaa*p)**3-((k_aaa-1)*p)**3)
      aaa=aaa/den
      endif
      do kk=1,k_bis_max
      den=n_h2o**2/2.*4.*pi/3.*((k_bis(kk)*p)**3-((k_bis(kk)-1)*p)**3)
      aaa_bis(10*(/1,2,3,4,5,6,7/)+kk)=aaa_bis(10*(/1,2,3,4,5,6,7/)+kk)/den
      end do
!
!
      aaa_bis(80)=aaa           ! je mets la variable aaa d'avant a la fin du tableau aaa_bis
!
      return
      end
!
!
function vector(r1,r2)
IMPLICIT REAL(8) (a-h,o-z)
DIMENSION r1(3),r2(3),vector(3)
!
!        fait le produit vectoriel r = r1 * r2
!        luc85 page 23
!
x3=r1(1); y3=r1(2); z3=r1(3)              ! r peut eventuellement remplacer en memoire r1
vector(1)=y3*r2(3)-z3*r2(2)
vector(2)=z3*r2(1)-x3*r2(3)
vector(3)=x3*r2(2)-y3*r2(1)
end
!
function harm_sph(m,mu,mup,beta)
implicit real(8) (a-h,o-z)
COMMON/facto/fac(0:101)
!
!          pour harmonique spherique Rm,mu,mup(Omega=omega,beta,phi)
!                                    =exp(-i*omega)*r(m)mu,mup*exp(-i*phi)
!          calcul l'element mu,mup de la matrice r(m) en fonction de l'angle beta
!          formule de Wigner dans Messiah eq.72 p922 betement
!          luc72p143
!          si mu ou mup nul, formule de recurrence stable
!
IF(ABS(mu)>m.or.ABS(mup)>m) THEN; harm_sph=0.; RETURN; endif
!
IF(mu==0.or.mup==0) THEN                      ! mu=0 ou mup=0
!
mu0=mu; beta0=beta
IF(mu==0) THEN; mu0=mup; beta0=-beta; ENDIF       ! je mets le 0 en second
!          si mu negatif, ca vaut (-1)**mu * valeur pour -mu
x=1.d0;
IF(mu0<0) THEN; x=(-1)**mu0; mu0=-mu0; endif
cc=COS(beta0)                            ! plutôt a partir d'une formule de recurrence stable
pm1=0.                                  ! des polynomes de Legendre associes Pl,m
pm=1.d0                                 ! luc73p96
do l=mu0+1,m
pm2=pm1
pm1=pm
pm=(cc*dble(2*l-1)*pm1-dble(l+mu0-1)*pm2)/DBLE(l-mu0)
end do
harm_sph=x*(-1)**mu0*SQRT(fac(m-mu0)/fac(m+mu0))*fac(2*mu0)/(2.d0**mu0*fac(mu0))*SIN(beta0)**mu0*pm
!
        else                  !   donc mu et mup non nuls, utiliser betement la formule de Wigner
!
harm_sph=0.
cc=COS(0.5d0*beta)
ss=SIN(0.5d0*beta)
do it=MAX(0,mu-mup),MIN(m+mu,m-mup)
harm_sph=harm_sph+(-1)**it/(fac(m+mu-it)*fac(m-mup-it)*fac(it)*fac(it-mu+mup))*  &
                           cc**(2*m+mu-mup-2*it)*ss**(2*it-mu+mup)
end do
harm_sph=SQRT(fac(m+mu)*fac(m-mu)*fac(m+mup)*fac(m-mup))*harm_sph
        endif
end


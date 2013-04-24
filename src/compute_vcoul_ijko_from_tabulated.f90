subroutine compute_vcoul_ijko_from_tabulated(nrgrid,drgrid,nb_id_solv,nb_id_mol,Rotxx,Rotxy,Rotxz,Rotyx,Rotyy,&
                                        Rotyz,Rotzx,Rotzy,Rotzz,tabulated_coulsq )


use precision_kinds, only: dp,i2b

use system, only: nfft1,nfft2,nfft3,deltax,deltay,deltaz,nb_omega,nb_psi,id_solv,id_mol,x_solv,y_solv,z_solv,x_mol,y_mol,z_mol,&
                        beta,nb_solute_sites,nb_solvent_sites,chg_mol,chg_solv,Lx,Ly,Lz,Rc , nb_species

use constants , only : fourpi , qfact

use external_potential , only : Vext_q


implicit none


integer(kind=i2b),intent(in) :: nrgrid ! total number of radial grid nodes
real(kind=dp), intent(in) :: drgrid ! distance between two radial grid nodes
integer(kind=i2b), intent(in) :: nb_id_solv,nb_id_mol ! number of different kinds of solvent sites or solute sites
real(kind=dp), dimension(nb_omega,nb_psi), intent(in) :: Rotxx,Rotxy,Rotxz,Rotyx,Rotyy,Rotyz,Rotzx,Rotzy,Rotzz
real(kind=dp), dimension(nb_id_solv,nb_id_mol,0:nrgrid), intent(in) :: tabulated_coulsq

integer(kind=i2b) :: i,j,k,o,p,m,n ! dummy
real(kind=dp),dimension(nb_solvent_sites,nb_psi,nb_omega) :: xmod,ymod,zmod
real(kind=dp) :: x_grid,y_grid,z_grid ! coordinates of grid nodes
real(kind=dp) :: x_m,y_m,z_m ! solvent sites coordinates
real(kind=dp) :: x_nm,y_nm,z_nm ! coordinate of vecteur solute-solvent
real(kind=dp) :: r_nm2 ! norm**2 of vector x_nm;y_nm;z_nm
real(kind=dp) :: V_psi ! dummy
integer(kind=i2b) :: ids,idm ! id of the solvent site or the solute site that is being used
real(kind=dp) :: time0,time1 ! timer start and end
real(kind=dp), dimension(:,:,:), allocatable :: temparray
character(50):: filename
integer(kind=i2b) :: px,py,pz
real(kind=dp) :: qfactcc ! qfact*chg_solv()*chg_mol()
real(kind=dp) :: rc2 ! charge pseudo radius **2   == Rc**2
real(kind=dp) :: tempVcoul ! temporary Vcoul(i,j,k,o)
integer(kind=i2b) :: ou_on_en_est ! where we are in the loop to follow on screen. goes from nfft3 to 0 when finished

integer ( kind = i2b ) :: species ! dummy for loops over species


print*,'!!!!!!!WARNING!!!!!!!WARNING!!!!!!!WARNING!!!!!!!WARNING!!!!!!!WARNING!!!!!!!WARNING!!!!!!!WARNING!!!!!!!WARNING!!!!!'
print*,'This routine use Rc which is defined in dft.in know what you are doing'
write(*,*)'Rc = ',Rc 
print*,'!!!!!!!WARNING!!!!!!!WARNING!!!!!!!WARNING!!!!!!!WARNING!!!!!!!WARNING!!!!!!!WARNING!!!!!!!WARNING!!!!!!!WARNING!!!!!'



!> initiate
call cpu_time(time0)
Vext_q = 0.0_dp
px=0
py=0
pz=0
Rc2=Rc**2
ou_on_en_est=nfft3


! test if all solutes have zero charge then don't waste your time : go to end of subroutine

if (minval(chg_mol)==0.0_dp .and. maxval(chg_mol)==0.0_dp) then   ! solute is not charged

  write(*,*)'solute has no charge'

  go to 999 ! end of subroutine

else if (minval(chg_solv)==0.0_dp .and. maxval(chg_solv)==0.0_dp) then   ! solvent is not charged

  write(*,*)'solvent has no charge'

  go to 999 ! end of subroutine

end if



!> tabulate Rotij(omega,psi)*k_solv(a) for speeding up
do o=1,nb_omega
  do p=1,nb_psi
    do m=1,nb_solvent_sites
      xmod(m,p,o)= Rotxx(o,p)*x_solv(m) + Rotxy(o,p)*y_solv(m) + Rotxz(o,p)*z_solv(m)
      ymod(m,p,o)= Rotyx(o,p)*x_solv(m) + Rotyy(o,p)*y_solv(m) + Rotyz(o,p)*z_solv(m)
      zmod(m,p,o)= Rotzx(o,p)*x_solv(m) + Rotzy(o,p)*y_solv(m) + Rotzz(o,p)*z_solv(m)
    end do
  end do
end do


!$OMP PARALLEL DO &
!$OMP DEFAULT(private) & ! variables are all private by default, but the ones in next SHARED(..)
!$OMP SHARED(nfft3,nfft2,nfft1,ou_on_en_est,deltaz,deltay,deltax,nb_omega,nb_psi,nb_solvent_sites, &
!$OMP id_solv,chg_solv,xmod,ymod,zmod,nb_solute_sites,id_mol,chg_mol,x_mol,y_mol,z_mol,Lx,Ly,Lz, &
!$OMP Rc2,beta,Vcoul) &
!$OMP SCHEDULE(DYNAMIC) ! dynamic allocation of work over nodes (no node is sleeping while others work)
do species = 1 , nb_species

do k=1,nfft3
  ou_on_en_est=ou_on_en_est-1
  write(*,*)ou_on_en_est
  z_grid=real(k-1,kind=dp)*deltaz

!  if(z_grid>12.0_dp .and. z_grid<30.0_dp) then
!    Vcoul(:,:,k,:)=100.0_dp
!    cycle
!  end if

  do j=1,nfft2
    y_grid=real(j-1,kind=dp)*deltay

    do i=1,nfft1
      x_grid=real(i-1,kind=dp)*deltax

      do o=1,nb_omega
        tempVcoul=0.0_dp

ploop:  do p=1,nb_psi
          V_psi=0.0_dp

          do m=1,nb_solvent_sites
            ids=id_solv(m)
            if (chg_solv(ids)==0.0_dp) cycle
            x_m = x_grid + xmod(m,p,o)
            y_m = y_grid + ymod(m,p,o)
            z_m = z_grid + zmod(m,p,o)
            do n=1,nb_solute_sites
              idm=id_mol(n)
              if (chg_mol(idm)==0.0_dp) cycle
              qfactcc=qfact*chg_solv(ids)*chg_mol(idm)
           !   write(*,*) qfactcc, qfact
        !       consider every period in every direction => 3*3*3=27 solutes
              do px=-1,1
                x_nm= x_m-(x_mol(n)+px*Lx)
               do py=-1,1
                  y_nm= y_m-(y_mol(n)+py*Ly)
                  do pz=-1,1
                    z_nm= z_m-(z_mol(n)+pz*Lz)
                    r_nm2 = x_nm**2+y_nm**2+z_nm**2
!                    i_r = int(r_nm2/drgrid+0.5_dp) ! FOR TABULATION USE ONLY
                    if(r_nm2<Rc2) then
                      V_psi = huge(1.0_dp)
                      cycle ploop ! sous-entendu tempVcoul=tempVcoul+exp(-beta*huge)=tempVcoul+0
                    else
                      V_psi = V_psi+qfactcc/sqrt(r_nm2)!tabulated_coulsq(ids,idm,i_r)
!                      V_psi = V_psi+tabulated_coulsq(ids,idm,i_r)
                    end if
               end do
              end do
            end do
            end do ! solute
          end do ! solvent

        !  tempVcoul=tempVcoul+exp(-beta*V_psi)
       

  !      if (tempVcoul<1.e-10_dp) then
          Vext_q(i,j,k,o,p,species)=V_psi
  !      else
  !        Vext_q(i,j,k,o,1)=-log(tempVcoul/nb_psi)/beta
  !      end if

        ! Limit maximum value
        if (Vext_q(i,j,k,o,p,species)>100.0_dp) Vext_q(i,j,k,o,p,species)=100.0_dp

!!!!!        if (Vcoul(i,j,k,o)<=-10.0_dp) then
!!!!!          write(*,*)'problem in compute_vcoul_ijko_from_tabulated.f90'
!!!!!          write(*,*)'reduced coordinate of problematic grid point : ',x_grid/Lx,y_grid/Ly,z_grid/Lz
!!!!!          write(*,*)'Vcoul(i,j,k,o) = ',Vcoul(i,j,k,o)
!!!!!        else if (Vcoul(i,j,k,o)>100.0_dp) then
!!!!!          Vcoul(i,j,k,o)=100.0_dp
!!!!!        end if
         end do ploop ! psi
      end do ! omega
    end do ! i
  end do ! j
end do ! k

end do ! species
!$OMP END PARALLEL DO



999 continue


write(*,*) 'minval(Vcoul)= ' , minval( Vext_q (:,:,:,:,:,:) )

write(*,*) 'maxval(Vcoul)= ' , maxval( Vext_q (:,:,:,:,:,:) )

call cpu_time(time1)

write(*,*)'time for compute_vcoul_ijko_from_tabulated ' , time1 - time0

!> Get the external potential over orientations and print it

allocate ( temparray ( nfft1 , nfft2 , nfft3 ) )

call mean_over_orientations ( Vext_q (:,:,:,:,:,1) , temparray )

temparray = temparray / fourpi

filename = 'output/Vcoul.cube'

call write_to_cube_file ( temparray , filename )

filename = 'output/Vcoul_along-z.dat'

call compute_z_density ( temparray , filename )

deallocate(temparray)




end subroutine compute_vcoul_ijko_from_tabulated

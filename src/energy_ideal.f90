! This subroutine computes the ideal part of the free energy functional.
! It sum over species of n(log(n)-1)
subroutine energy_ideal
use precision_kinds , only : i2b , dp
use cg , only : cg_vect , FF , dF
use system, only : nfft1, nfft2, nfft3, nb_omega, deltaV, rho_0, kBT, nb_species, rho_0_multispec, mole_fraction, nb_psi,&
n_0_multispec
use quadrature , only : weight , weight_psi,sym_order
use input, only : input_log, input_char
use constants, only : fourpi
implicit none
real(dp):: Fideal, Fid_lin ! ideal free energy, linearized ideal free energy 
real (dp) :: Fid_lin_temp, dFid_lin_temp !dummy for temporary store linearized ideal free energy and the corresponding gradient
integer(i2b):: icg , i , j , k , o , p! dummy for loops
integer(i2b):: species ! dummy between 1 and nb_species
real(dp):: psi ! dummy for cg_vext(i)
real(dp):: rho, rhon ! local density
real(dp):: logrho ! dummy for log(rho)
real(dp):: time0 , time1 ! timesteps
real(dp), dimension(nfft1,nfft2,nfft3) :: rho_n  !one-particle number density
! init timer
call cpu_time(time0)
! init Fideal to zero and its gradient
Fideal = 0.0_dp
print*, 'BETA' ,1/kbt, rho_0_multispec
! the following is a dummy for speeding up loops
!get rho_n, currently implemented only for 1 species
if (input_log('Linearize_entropy') .and. trim(adjustl(input_char('if_Linearize_entropy'))) == '1' ) then
  if (nb_species/=1 ) then
    print*, 'the linearized ideal F is only implemented for 1 specie for the moment'
    stop
  end if
end if
icg=0
Fid_lin=0.0_dp
if (input_log('Linearize_entropy').and. trim(adjustl(input_char('if_Linearize_entropy'))) == '1' ) then
  do i=1,nfft1
    do j=1,nfft2
      do k=1,nfft3
        rhon=0.0_dp
        do o=1,nb_omega
          do p=1,nb_psi
            icg=icg+1
            rhon=rhon+cg_vect(icg)**2*weight(o)*weight_psi(p)/(fourpi**2/(sym_order*2.0_dp))
          end do
        end do
        rho_n(i,j,k)=rhon
        if (rhon <=1.0_dp) then
            Fid_lin_temp=0.0_dp
        else  
            Fid_lin_temp=-(rhon*Log(rhon)-rhon+1.0_dp-0.5_dp*(rhon-1.0_dp)**2)
        end if
        Fid_lin=Fid_lin+Fid_lin_temp
      end do
    end do
  end do
end if
!put linearized ideal free energy in kJ/mol
Fid_lin=Fid_lin*kBT*n_0_multispec(1)*deltaV
print*, Fid_lin
print*, minval(rho_n), maxval(rho_n)
icg = 0
if (input_log('Linearize_entropy')  ) then
  if (trim(adjustl(input_char('if_Linearize_entropy'))) == '1') then
    do species = 1 , nb_species
      do i = 1 , nfft1
        do j = 1 , nfft2
          do k = 1 , nfft3        
            do o = 1 , nb_omega    
              do p=1 , nb_psi
                icg = icg + 1
                psi = CG_vect ( icg )
                rhon=rho_n(i,j,k)
                if (rhon<=1.0_dp) then
                   dFid_lin_temp=0.0_dp
                else
                   dFid_lin_temp=-KbT*(Log(rhon)-rhon+1.0_dp)
                end if
                if ( psi <= 0.0_dp ) then ! <= because sometimes comes -0.0_dp
                  Fideal = Fideal + weight ( o ) *weight_psi(p)* rho_0_multispec ( species ) * mole_fraction ( species ) ! lim xlogx= 0 when x->0
                  dF (icg) = dF ( icg ) + 0.0_dp ! lim x logx = 0.0
                else
                  rho = psi ** 2
                  logrho = log ( rho )
                  Fideal = Fideal + weight ( o ) * weight_psi (p) * rho_0_multispec ( species ) * mole_fraction ( species ) &
                                    * ( rho * logrho - rho + 1.0_dp )
                  dF (icg) = dF ( icg ) + 2.0_dp * psi * weight ( o ) * weight_psi(p) * DeltaV * rho_0_multispec ( species )&
                                    *( kBT * logrho + dFid_lin_temp )
                end if
              end do ! nb_psi
            end do ! nb_omega
          end do ! nfft3
        end do ! nfft2
      end do ! nfft1
    end do ! nb_species
  else if (trim(adjustl(input_char('if_Linearize_entropy'))) == '2') then
    Fid_lin=0.0_dp
    do species = 1 , nb_species
      do i = 1 , nfft1
        do j = 1 , nfft2
          do k = 1 , nfft3        
            do o = 1 , nb_omega    
              do p=1 , nb_psi
                icg = icg + 1
                psi = CG_vect ( icg )
                rho=psi**2
                logrho = log ( rho )
                if (rho<=1.0_dp) then
                   dFid_lin_temp=0.0_dp
                   Fid_lin_temp=0.0_dp
                else
                   dFid_lin_temp=-KbT*(Logrho-rho+1.0_dp)
                   Fid_lin_temp=-KbT*(rho*Logrho-rho+1.0_dp-0.5_dp*(rho-1.0_dp)**2)
                end if
                if ( psi <= 0.0_dp ) then ! <= because sometimes comes -0.0_dp
                  Fideal = Fideal + weight ( o ) *weight_psi(p)* rho_0_multispec ( species ) * mole_fraction ( species ) ! lim xlogx= 0 when x->0
                  dF (icg) = dF ( icg ) + 0.0_dp ! lim x logx = 0.0
                else
                  Fideal = Fideal + weight ( o ) * weight_psi (p) * rho_0_multispec ( species ) * mole_fraction ( species ) &
                                    * ( rho * logrho - rho + 1.0_dp )
                  Fid_lin=Fid_lin+weight(o)*weight_psi(p)*DeltaV * rho_0_multispec ( species )*Fid_lin_temp
                  dF (icg) = dF ( icg ) + 2.0_dp * psi * weight ( o ) * weight_psi(p) * DeltaV * rho_0_multispec ( species )&
                                    *( kBT * logrho + dFid_lin_temp )
                end if
              end do ! nb_psi
            end do ! nb_omega
          end do ! nfft3
        end do ! nfft2
      end do ! nfft1
    end do ! nb_species
  else
    print*, 'error in energy_ideal, you specified that you wanna Linearize the entropy but you did not choose HOW:1 for n 2 for rho'
    stop
  end if
else
  do species = 1 , nb_species
    do i = 1 , nfft1
      do j = 1 , nfft2
        do k = 1 , nfft3
          do o = 1 , nb_omega
            do p=1 , nb_psi
              icg = icg + 1
              psi = CG_vect ( icg )
              if ( psi <= 0.0_dp ) then ! <= because sometimes comes -0.0_dp
                Fideal = Fideal + weight ( o ) *weight_psi(p)* rho_0_multispec ( species ) * mole_fraction ( species ) ! lim xlogx= 0 when x->0
                dF (icg) = dF ( icg ) + 0.0_dp ! lim x logx = 0.0
            else
                rho = psi ** 2
                logrho = log ( rho )
                Fideal = Fideal + weight ( o ) * weight_psi (p) * rho_0_multispec ( species ) * mole_fraction ( species ) &
                                  * ( rho * logrho - rho + 1.0_dp )
                dF (icg) = dF ( icg ) + 2.0_dp * psi * weight ( o ) * weight_psi(p) * DeltaV * rho_0_multispec ( species )&
                                  * kBT * logrho
            end if
            end do ! nb_psi
          end do ! nb_omega
        end do ! nfft3
      end do ! nfft2
    end do ! nfft1
  end do ! nb_species
end if
! integration factor
Fideal = Fideal * kBT * DeltaV
! conclude
print*, 'Fid_lin =' , Fid_lin
FF = FF + Fideal +Fid_lin
! warn user
call cpu_time ( time1 )
write (*,*) 'ideal       = ' , Fideal , 'computed in (sec)' , time1 - time0
end subroutine energy_ideal

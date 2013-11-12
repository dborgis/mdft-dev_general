! This subroutine computes the ideal part of the free energy functional.
! It sum over species of n(log(n)-1)

subroutine energy_ideal

    USE precision_kinds, ONLY: i2b, dp
    USE cg, ONLY: cg_vect, FF, dF
    USE system, ONLY: nfft1, nfft2, nfft3, lx, ly, lz, rho_0, kBT, nb_species, rho_0_multispec, mole_fraction, &
                        n_0_multispec, spaceGrid
    USE quadrature, ONLY: sym_order, angGrid, molRotGrid
    USE input, ONLY: input_log, input_char
    USE constants, ONLY: fourpi

    IMPLICIT NONE
    
    REAL(dp):: Fideal, Fid_lin ! ideal free energy, linearized ideal free energy 
    REAL(dp) :: Fid_lin_temp, dFid_lin_temp !dummy for temporary store linearized ideal free energy and the corresponding gradient
    INTEGER(i2b):: icg , i , j , k , o , p, s! dummy for loops
    INTEGER(i2b):: species ! dummy between 1 and nb_species
    REAL(dp):: psi, tmp ! dummy for cg_vext(i)
    REAL(dp):: rho, rhon ! local density
    REAL(dp), POINTER :: deltaV => spaceGrid%dv
    REAL(dp):: logrho ! dummy for log(rho)
    REAL(dp):: time0 , time1 ! timesteps
    REAL(dp), dimension(nfft1,nfft2,nfft3) :: rho_n  !one-particle number density

    CALL CPU_TIME (time0) ! init timer

    Fideal = 0.0_dp! init Fideal to zero and its gradient

    icg=0
    Fid_lin=0.0_dp
    IF (input_log('Linearize_entropy').and. trim(adjustl(input_char('if_Linearize_entropy'))) == '1' ) then
        IF (nb_species/=1 ) STOP 'the linearized ideal F is only implemented for 1 specie for the moment'
        do i=1,nfft1
            do j=1,nfft2
                do k=1,nfft3
                    rhon=0.0_dp
                    do o=1,angGrid%n_angles
                        do p=1,molRotGrid%n_angles
                            icg=icg+1
                            rhon=rhon+cg_vect(icg)**2*angGrid%weight(o)*molRotGrid%weight(p)/(fourpi**2/(sym_order*2.0_dp))
                        end do
                    end do
                    rho_n(i,j,k)=rhon
                    if (rhon <=1.0_dp) then
                        Fid_lin_temp=0.0_dp
                    else  
                        Fid_lin_temp=-(rhon*Log(rhon)-rhon+1.0_dp-0.5_dp*(rhon-1.0_dp)**2)
                    end if
                    Fid_lin = Fid_lin+Fid_lin_temp
                end do
            end do
        end do
    END IF
    Fid_lin = Fid_lin*kBT*n_0_multispec(1)*spaceGrid%dv ! convert to kJ/mol



    icg = 0
    IF (input_log('Linearize_entropy')  ) THEN
    
    
        IF (trim(adjustl(input_char('if_Linearize_entropy'))) == '1') then
            do species = 1 , nb_species
                do i = 1 , nfft1
                    do j = 1 , nfft2
                        do k = 1 , nfft3        
                            do o = 1 , angGrid%n_angles    
                                do p=1 , molRotGrid%n_angles
                                    icg = icg + 1
                                    psi = CG_vect ( icg )
                                    rhon = rho_n(i,j,k)
                                    if (rhon<=1.0_dp) then
                                        dFid_lin_temp=0.0_dp
                                    else
                                        dFid_lin_temp = -KbT*(Log(rhon)-rhon+1.0_dp)
                                    end if
                                    rho = psi ** 2
                                    Fideal = Fideal + Fideal_local(o,p,s,rho)
                                    dF(icg) =  dF(icg) + dFideal_local(o,p,s,psi,dFid_lin_temp)
                                end do
                            end do
                        end do
                    end do
                end do
            end do

        ELSE IF (trim(adjustl(input_char('if_Linearize_entropy'))) == '2') then
            Fid_lin = 0.0_dp
            do s = 1 , nb_species
                do i = 1 , nfft1
                    do j = 1 , nfft2
                        do k = 1 , nfft3        
                            do o = 1 , angGrid%n_angles    
                                do p=1 , molRotGrid%n_angles
                                    icg = icg + 1
                                    psi = CG_vect ( icg )
                                    rho=psi**2
                                    logrho = log ( rho )
                                    IF ( rho > 0._dp ) THEN
                                        dFid_lin_temp=-KbT*(Logrho-rho+1.0_dp)
                                        Fid_lin_temp=-KbT*(rho*Logrho-rho+1.0_dp-0.5_dp*(rho-1.0_dp)**2)
                                    END IF
                                    Fideal = Fideal + Fideal_local (o,p,s,rho)
                                    IF ( psi /= 0.0_dp ) THEN
                                        Fid_lin = Fid_lin + prefactor(o,p,s) *DeltaV *Fid_lin_temp
                                        dF(icg) = dF(icg) + dFideal_local (o,p,s,psi,dFid_lin_temp)
                                    END IF
                                end do
                            end do
                        end do
                    end do
                end do
            end do
            
        ELSE
            STOP 'error in energy_ideal, you want to linearize the entropy but you did not choose HOW:1 for n 2 for rho'
        END IF
    
    ELSE
        do s = 1 , nb_species
            do i = 1 , nfft1
                do j = 1 , nfft2
                    do k = 1 , nfft3
                        do o = 1 , angGrid%n_angles
                            do p = 1 , molRotGrid%n_angles
                                icg = icg + 1
                                psi = CG_vect (icg)
                                rho = psi**2
                                Fideal = Fideal + Fideal_local (o,p,s,rho)
                                dF (icg) = dF (icg) + dFideal_local (o,p,s,psi,0._dp)
                            end do
                        end do
                    end do
                end do
            end do
        end do
    END IF
    Fideal = Fideal * kBT * spaceGrid%dv ! integration factor

    FF = FF + Fideal + Fid_lin
    
    CALL CPU_TIME (time1)
    WRITE(*,*) 'Fideal      = ' , Fideal , 'computed in (sec)' , time1 - time0


    CONTAINS


    
    PURE FUNCTION dFideal_local (o,p,s,psi,toadd)
        INTEGER(i2b), INTENT(IN) :: o,p,s
        REAL(dp), INTENT(IN) :: psi, toadd
        REAL(dp) :: dFideal_local
        IF (psi /= 0._dp) THEN
            dFideal_local = 2.0_dp * psi * prefactor(o,p,s) * spaceGrid%dv * ( kBT*LOG(psi**2) + toadd )
        ELSE
            dFideal_local = 0._dp
        END IF
    END FUNCTION dFideal_local

    
    
    PURE FUNCTION Fideal_local (o,p,s,rho)
        INTEGER(i2b), INTENT(IN) :: o,p,s
        REAL(dp), INTENT(IN) :: rho
        REAL(dp) :: Fideal_local
        IF (rho /= 0._dp) THEN
            Fideal_local = prefactor(o,p,s) * mole_fraction(s) * (rho*LOG(rho)-rho+1.0_dp)
        ELSE IF ( rho == 0._dp ) THEN
            Fideal_local = prefactor(o,p,s) * mole_fraction(s)
        END IF
    END FUNCTION Fideal_local


    PURE FUNCTION prefactor (o,p,s)
        INTEGER(i2b), INTENT(IN) :: o,p,s
        REAL(dp) :: prefactor
        prefactor = angGrid%weight(o) * molRotGrid%weight(p) * rho_0_multispec(s)
    END FUNCTION


END SUBROUTINE energy_ideal

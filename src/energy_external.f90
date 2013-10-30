! this subroutine compute the external part of the free energy functional
subroutine energy_external

    use precision_kinds , only : dp , i2b
    use system , only : lx,ly,lz, nfft1, nfft2, nfft3, nb_species, rhoBulk=>rho_0_multispec
    use quadrature , only : angGrid, molRotGrid
    use cg , only : CG_vect , FF , dF
    use external_potential , only : Vext_total
    
    implicit none
    real(dp):: Fext , Fext_q! external part of the total energy
    integer(i2b):: icg , i , j , k , o , p! dummy for loops
    real(dp):: psi !dummy
    real(dp):: rho ! dummy = psi**2
    real(dp):: wdfve, deltaV
    real(dp):: time0,time1! timers
    integer(i2b):: spec  ! dummy between 1 and nb_species
    
    call cpu_time ( time0 )

    Fext = 0.0_dp ; Fext_q = 0.0_dp

    ! Compute Fext and dFext
    ! F_{ext}[\rho(\vec{r},\vec{\Omega})]=\int d \vec{r} d \vec{\Omega} V_{ext}(\vec{r},\vec{\Omega})\rho(\vec{r},\vec{\Omega})
    deltaV = lx*ly*lz/real(nfft1*nfft2*nfft3,dp) ! volume integration
    icg = 0
    do spec = 1 , nb_species
        do i = 1 , nfft1
            do j = 1 , nfft2
                do k = 1 , nfft3
                    do o = 1 , angGrid%n_angles
                        do p=1 , molRotGrid%n_angles            
                            icg = icg + 1
                            psi = cg_vect ( icg )
                            rho = psi ** 2
                            wdfve = angGrid%weight(o) * molRotGrid%weight(p) * DeltaV * rhoBulk(spec) * Vext_total(i,j,k,o,p,spec)
                            Fext = Fext + rho * wdfve
                            dF(icg) = dF(icg) + 2.0_dp*psi*wdfve
                        end do
                    end do
                end do
            end do
        end do
    end do

    FF = FF + Fext

    call cpu_time(time1)
    write (*,*) 'external    = ' , Fext , 'computed in (sec)' , time1 - time0

end subroutine energy_external

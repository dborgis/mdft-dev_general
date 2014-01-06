!This routine compute : -The solvent molecular charge density, which can be used into Vcoul_from_solvent_charge_density.f90 to 
!evaluate the electrostatic potential.
!                       -The solvent molecular polarization (from Ranieriet al : J. Chem. Phys. 98 (11) 1993) which ca be used
!into energy_polarization_myway.f90 to compute the (multipolar) polarization Free energy.

SUBROUTINE chargeDensityAndMolecularPolarizationOfASolventMoleculeAtOrigin (Rotxx,Rotxy,Rotxz,Rotyx,Rotyy,Rotyz,Rotzx,Rotzy,Rotzz)

    USE precision_kinds, ONLY: i2b, dp
    USE constants, ONLY: i_complex, twopi
    USE system, ONLY: chg_solv, x_solv, y_solv, z_solv, nb_solvent_sites, id_solv,&
                      sigma_k,molec_polarx_k, molec_polary_k, molec_polarz_k, nb_species, spaceGrid
    USE quadrature, ONLY: angGrid, molRotGrid
    USE fft, ONLY: kx, ky, kz, k2
    
    IMPLICIT NONE
    
    INTEGER(i2b) :: i, j, k, o, p, n, nf1, s, nfft1, nfft2, nfft3
    REAL(dp), DIMENSION(angGrid%n_angles,molRotGrid%n_angles), INTENT(IN) :: Rotxx,Rotxy,Rotxz,Rotyx,Rotyy,Rotyz,Rotzx,Rotzy,Rotzz
    REAL(dp) :: xmod, ymod, zmod, deltaVk, Rc, Lx, Ly, Lz, chg_, kr
    COMPLEX(dp), PARAMETER :: zeroC=(0.0_dp,0.0_dp)
    COMPLEX(dp) :: fac

    Lx = spaceGrid%length(1)
    Ly = spaceGrid%length(2)
    Lz = spaceGrid%length(3)
    deltaVk = twopi**3/(Lx*Ly*Lz)
    nfft1 = spaceGrid%n_nodes(1)
    nfft2 = spaceGrid%n_nodes(2)
    nfft3 = spaceGrid%n_nodes(3)
    nf1 = nfft1/2

    !~ real (dp), dimension(nfft1,nfft2,nfft3)::molecpolarx,molecpolary,molecpolarz
    !            ====================================================
    !            !    	Initialization				!
    !            !							!
    !            ====================================================
    Rc=0.5_dp
    IF (Rc/=0.0_dp) THEN
        PRINT*,'WARNING: you convolute molecular Charge Density and POLARIZATION with a Gaussian be sure that is what you want!!!'
    END IF
    
    ALLOCATE (sigma_k (nf1+1, nfft2, nfft3, angGrid%n_angles, molRotGrid%n_angles,nb_species), SOURCE=zeroC)
    ALLOCATE (molec_polarx_k (nf1+1 , nfft2 , nfft3 , angGrid%n_angles , molRotGrid%n_angles,nb_species), SOURCE=zeroC)
    ALLOCATE (molec_polary_k (nf1+1 , nfft2 , nfft3 , angGrid%n_angles , molRotGrid%n_angles,nb_species), SOURCE=zeroC)
    ALLOCATE (molec_polarz_k (nf1+1 , nfft2 , nfft3 , angGrid%n_angles , molRotGrid%n_angles,nb_species), SOURCE=zeroC)

    !            ====================================================
    !            !    Compute sigma and molecular polarization   	!
    !            !							                        !
    !            ====================================================
    DO s= 1, nb_species
        DO i= 1, spaceGrid%n_nodes(1)/2 +1
            DO j= 1, spaceGrid%n_nodes(2)
                DO k= 1, spaceGrid%n_nodes(3)
                    DO o= 1, angGrid%n_angles
                        DO p= 1, molRotGrid%n_angles
                            DO n= 1, nb_solvent_sites
                                chg_ = chg_solv(id_solv(n))
                                IF ( chg_ == 0._dp ) CYCLE
                                xmod= Rotxx(o,p)*x_solv(n) + Rotxy(o,p)*y_solv(n) + Rotxz(o,p)*z_solv(n)
                                ymod= Rotyx(o,p)*x_solv(n) + Rotyy(o,p)*y_solv(n) + Rotyz(o,p)*z_solv(n)   
                                zmod= Rotzx(o,p)*x_solv(n) + Rotzy(o,p)*y_solv(n) + Rotzz(o,p)*z_solv(n)  
                                sigma_k (i,j,k,o,p,s) = sigma_k(i,j,k,o,p,s) + chg_*&
                                    EXP(-i_complex*(kx(i)*xmod+ky(j)*ymod+ kz(k)*zmod ))*EXP(-(Rc**2*k2(i,j,k))/2)
                                kr = xmod*kx(i)+ymod*ky(j)+zmod*kz(k)
                                IF ( kr==0.0_dp ) THEN
                                    molec_polarx_k (i,j,k,o,p,s) = molec_polarx_k(i,j,k,o,p,s) + chg_*xmod
                                    molec_polary_k (i,j,k,o,p,s) = molec_polary_k(i,j,k,o,p,s) + chg_*ymod
                                    molec_polarz_k (i,j,k,o,p,s) = molec_polarz_k(i,j,k,o,p,s) + chg_*zmod
                                ELSE
                                    fac = 1._dp/kr*(EXP(i_complex*kr)-1)*(-i_complex) *EXP(-(Rc**2*k2(i,j,k))/2)
                                    molec_polarx_k (i,j,k,o,p,s) = molec_polarx_k(i,j,k,o,p,s) + chg_*xmod*fac
                                    molec_polary_k (i,j,k,o,p,s) = molec_polary_k(i,j,k,o,p,s) + chg_*ymod*fac
                                    molec_polarz_k (i,j,k,o,p,s) = molec_polarz_k(i,j,k,o,p,s) + chg_*zmod*fac
                                END IF
                            END DO
                        END DO
                    END DO
                END DO
            END DO
        END DO
    END DO

END SUBROUTINE chargeDensityAndMolecularPolarizationOfASolventMoleculeAtOrigin

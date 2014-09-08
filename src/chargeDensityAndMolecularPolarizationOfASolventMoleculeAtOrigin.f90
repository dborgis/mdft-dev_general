!This routine compute : -The solvent molecular charge density, which can be used into Vcoul_from_solvent_charge_density.f90 to 
!evaluate the electrostatic potential.
!                       -The solvent molecular polarization (from Ranieriet al : J. Chem. Phys. 98 (11) 1993) that can be used
!into energy_polarization_..._.f90 to compute the (multipolar) polarization Free energy.

SUBROUTINE chargeDensityAndMolecularPolarizationOfASolventMoleculeAtOrigin (Rotxx,Rotxy,Rotxz,Rotyx,Rotyy,Rotyz,Rotzx,Rotzy,Rotzz)

    USE precision_kinds     ,ONLY: i2b, dp
    USE constants           ,ONLY: iC=>i_complex, zeroC
    USE system              ,ONLY: solvent, sigma_k, molec_polarx_k, molec_polary_k, molec_polarz_k, nb_species, spaceGrid
    USE quadrature          ,ONLY: angGrid, molRotGrid
    USE fft                 ,ONLY: kx, ky, kz, k2

    IMPLICIT NONE
    
    REAL(dp), DIMENSION(angGrid%n_angles,molRotGrid%n_angles), INTENT(IN) :: Rotxx,Rotxy,Rotxz,Rotyx,Rotyy,Rotyz,Rotzx,Rotzy,Rotzz
    INTEGER(i2b) :: i, j, k, o, p, n, nf1, s, nfft1, nfft2, nfft3
    REAL(dp)     :: xmod, ymod, zmod, Rc, Lx, Ly, Lz, kr,angle_number
    COMPLEX(dp)  :: fac

    Lx = spaceGrid%length(1)
    Ly = spaceGrid%length(2)
    Lz = spaceGrid%length(3)
    nfft1 = spaceGrid%n_nodes(1)
    nfft2 = spaceGrid%n_nodes(2)
    nfft3 = spaceGrid%n_nodes(3)
    nf1 = nfft1/2

    !            ====================================================
    !            !    	Initialization				!
    !            !							!
    !            ====================================================
    Rc=0.5_dp
    IF (Rc/=0.0_dp) THEN
        PRINT*,'WARNING: you convolute the solvent molecular charge density and the solvent molecular polarization &
                with a Gaussian of radius 0.5. Angstroms in chargeDensity...'
    END IF
    
    ALLOCATE (sigma_k        (nf1+1, nfft2, nfft3, angGrid%n_angles, molRotGrid%n_angles, nb_species), SOURCE=zeroC)
    ALLOCATE (molec_polarx_k (nf1+1, nfft2, nfft3, angGrid%n_angles, molRotGrid%n_angles, nb_species), SOURCE=zeroC)
    ALLOCATE (molec_polary_k (nf1+1, nfft2, nfft3, angGrid%n_angles, molRotGrid%n_angles, nb_species), SOURCE=zeroC)
    ALLOCATE (molec_polarz_k (nf1+1, nfft2, nfft3, angGrid%n_angles, molRotGrid%n_angles, nb_species), SOURCE=zeroC)

    !            ====================================================
    !            !    Compute sigma and molecular polarization   	!
    !            !							                        !
    !            ====================================================

    DO CONCURRENT ( i=1:nfft1/2+1, j=1:nfft2, k=1:nfft3, o=1:angGrid%n_angles, p=1:molRotGrid%n_angles, &
        n=1:SIZE(solvent(1)%site), s=1:nb_species, (solvent(1)%site(n)%q/=0._dp) )

        xmod= Rotxx(o,p)*solvent(1)%site(n)%r(1) + Rotxy(o,p)*solvent(1)%site(n)%r(2) + Rotxz(o,p)*solvent(1)%site(n)%r(3)
        ymod= Rotyx(o,p)*solvent(1)%site(n)%r(1) + Rotyy(o,p)*solvent(1)%site(n)%r(2) + Rotyz(o,p)*solvent(1)%site(n)%r(3)   
        zmod= Rotzx(o,p)*solvent(1)%site(n)%r(1) + Rotzy(o,p)*solvent(1)%site(n)%r(2) + Rotzz(o,p)*solvent(1)%site(n)%r(3)  
        kr = xmod*kx(i)+ymod*ky(j)+zmod*kz(k)

        sigma_k (i,j,k,o,p,s) = sigma_k(i,j,k,o,p,s) + solvent(1)%site(n)%q *EXP(-iC*kr) *EXP(-Rc**2*k2(i,j,k)/2)

        IF ( kr==0.0_dp ) THEN
            molec_polarx_k (i,j,k,o,p,s) = molec_polarx_k(i,j,k,o,p,s) + solvent(1)%site(n)%q*xmod
            molec_polary_k (i,j,k,o,p,s) = molec_polary_k(i,j,k,o,p,s) + solvent(1)%site(n)%q*ymod
            molec_polarz_k (i,j,k,o,p,s) = molec_polarz_k(i,j,k,o,p,s) + solvent(1)%site(n)%q*zmod
        ELSE
            fac = 1._dp/kr*(EXP(iC*kr)-1)*(-iC) *EXP(-(Rc**2*k2(i,j,k))/2)
            molec_polarx_k (i,j,k,o,p,s) = molec_polarx_k(i,j,k,o,p,s) + solvent(1)%site(n)%q*xmod*fac
            molec_polary_k (i,j,k,o,p,s) = molec_polary_k(i,j,k,o,p,s) + solvent(1)%site(n)%q*ymod*fac
            molec_polarz_k (i,j,k,o,p,s) = molec_polarz_k(i,j,k,o,p,s) + solvent(1)%site(n)%q*zmod*fac
        END IF

    END DO

    CALL normalizeMolecularDensity

    CONTAINS
        
        SUBROUTINE normalizeMolecularDensity
            USE constants ,ONLY: epsN
            REAL(dp)          :: angle_number
            COMPLEX(dp)       :: molxk, molyk,molzk
            INTEGER(i2b)      :: i,j,k,s
            i =angGrid%n_angles * molRotGrid%n_angles
            IF (i<1) STOP "I found no angle in normalizeMolecularDensity"
            angle_number =REAL(i,dp)
            DO CONCURRENT (i=1:nfft1/2+1, j=1:nfft2, k=1:nfft3, s=1:nb_species)
                molxk =SUM( molec_polarx_k (i,j,k,:,:,s) )
                molyk =SUM( molec_polary_k (i,j,k,:,:,s) )        
                molzk =SUM( molec_polarz_k (i,j,k,:,:,s) )
                molec_polarx_k (i,j,k,:,:,s) =molec_polarx_k (i,j,k,:,:,s)-molxk/angle_number
                molec_polary_k (i,j,k,:,:,s) =molec_polary_k (i,j,k,:,:,s)-molyk/angle_number
                molec_polarz_k (i,j,k,:,:,s) =molec_polarz_k (i,j,k,:,:,s)-molzk/angle_number
            END DO
        END SUBROUTINE normalizeMolecularDensity

END SUBROUTINE chargeDensityAndMolecularPolarizationOfASolventMoleculeAtOrigin

!This routine compute : -The solvent molecular charge density, which can be used into Vcoul_from_solvent_charge_density.f90 to
!evaluate the electrostatic potential.
!                       -The solvent molecular polarization (from Ranieriet al : J. Chem. Phys. 98 (11) 1993) that can be used
!into energy_polarization_..._.f90 to compute the (multipolar) polarization Free energy.

SUBROUTINE chargeDensityAndMolecularPolarizationOfASolventMoleculeAtOrigin (Rotxx,Rotxy,Rotxz,Rotyx,Rotyy,Rotyz,Rotzx,Rotzy,Rotzz)

    USE precision_kinds     ,ONLY: i2b, dp
    USE constants           ,ONLY: iC=>i_complex, zeroC
    USE system              ,ONLY: solvent, spaceGrid
    USE quadrature          ,ONLY: angGrid, molRotGrid
    USE fft                 ,ONLY: kx, ky, kz, k2

    IMPLICIT NONE

    REAL(dp), DIMENSION(angGrid%n_angles,molRotGrid%n_angles), INTENT(IN) :: Rotxx,Rotxy,Rotxz,Rotyx,Rotyy,Rotyz,Rotzx,Rotzy,Rotzz
    integer, pointer :: nfft1 => spacegrid%n_nodes(1), nfft2 => spacegrid%n_nodes(2), nfft3 => spacegrid%n_nodes(3)
    INTEGER(i2b) :: i, j, k, o, p, n, s
    REAL(dp)     :: r(3), Rc, kr, kvec(3), prefactor
    COMPLEX(dp)  :: fac

    !            ====================================================
    !            !    	Initialization				!
    !            !							!
    !            ====================================================
    Rc=0.0_dp
    PRINT*,"WARNING YOU CHANGED SIGMA_K RADIUS FROM 0.5 TO 0"
    IF (Rc/=0.0_dp) THEN
        PRINT*,'WARNING: you convolute the solvent molecular charge density and the solvent molecular polarization &
                &with a Gaussian of radius 0.5. Angstroms in chargeDensity...'
    END IF

    do concurrent( s=1:size(solvent) )
        allocate( solvent(s)%sigma_k(nfft1/2+1, nfft2, nfft3, angGrid%n_angles, molRotGrid%n_angles), SOURCE=zeroC)
        allocate( solvent(s)%molec_polar_k(3,nfft1/2+1, nfft2, nfft3, angGrid%n_angles, molRotGrid%n_angles), SOURCE=zeroC)
    end do

    !            ====================================================
    !            !    Compute sigma and molecular polarization   	!
    !            !							                        !
    !            ====================================================

    do concurrent ( i=1:nfft1/2+1, j=1:nfft2, k=1:nfft3, s=1:size(solvent) )

        kvec = [ kx(i), ky(j), kz(k) ]
        prefactor =  exp(-Rc**2*k2(i,j,k)/2)

        do concurrent ( o=1:angGrid%n_angles, p=1:molRotGrid%n_angles, &
                        n=1:SIZE(solvent(s)%site), (abs(solvent(s)%site(n)%q)>=epsilon(1.0_dp)) )

            r(1) = dot_product(   [Rotxx(o,p),Rotxy(o,p),Rotxz(o,p)]  ,  solvent(s)%site(n)%r  )
            r(2) = dot_product(   [Rotyx(o,p),Rotyy(o,p),Rotyz(o,p)]  ,  solvent(s)%site(n)%r  )
            r(3) = dot_product(   [Rotzx(o,p),Rotzy(o,p),Rotzz(o,p)]  ,  solvent(s)%site(n)%r  )

            kr = dot_product( kvec, r )

            solvent(s)%sigma_k(i,j,k,o,p) = solvent(s)%sigma_k(i,j,k,o,p) + solvent(s)%site(n)%q *EXP(-iC*kr) *prefactor

            if ( abs(kr)<=epsilon(1.0_dp) ) THEN
                solvent(s)%molec_polar_k(:,i,j,k,o,p) = solvent(s)%molec_polar_k(:,i,j,k,o,p) + solvent(s)%site(n)%q *r
            else
                fac = (EXP(iC*kr)-1)/kr*(-iC) *prefactor
                solvent(s)%molec_polar_k(:,i,j,k,o,p) = solvent(s)%molec_polar_k(:,i,j,k,o,p) + fac*solvent(s)%site(n)%q *r
            end if

        end do
    end do

    CALL normalizeMolecularDensity

    CONTAINS

        SUBROUTINE normalizeMolecularDensity
            REAL(dp)          :: angle_number
            INTEGER(i2b)      :: i,j,k,s,d
            angle_number = angGrid%n_angles * molRotGrid%n_angles
            IF (angle_number<1) STOP "I found no angle in normalizeMolecularDensity"
            DO CONCURRENT (i=1:nfft1/2+1, j=1:nfft2, k=1:nfft3, s=1:size(solvent), d=1:3)
                solvent(s)%molec_polar_k(d,i,j,k,:,:) = solvent(s)%molec_polar_k(d,i,j,k,:,:)  &
                                                       -SUM( solvent(s)%molec_polar_k(d,i,j,k,:,:) ) /real(angle_number,dp)
            END DO
        END SUBROUTINE normalizeMolecularDensity

END SUBROUTINE chargeDensityAndMolecularPolarizationOfASolventMoleculeAtOrigin

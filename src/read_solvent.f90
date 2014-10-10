!===================================================================================================================================
SUBROUTINE read_solvent
!===================================================================================================================================
! Read solvent atomic positions, charge, and lennard jones values in solvent.in
! charge in electron units, sigma in Angstroms, epsilon in KJ/mol.

    USE precision_kinds, ONLY: i2b , dp, sp
    USE system, ONLY: nb_solvent_sites, solvent
    use mathematica, only: chop
    
    IMPLICIT NONE

    INTEGER(i2b) :: n, ios, i, j, k, l, s

    OPEN(5, FILE= 'input/solvent.in', STATUS= 'old', IOSTAT= ios )! open input/solvent.in and check if it is readable
        IF ( ios/=0 ) STOP 'ERROR: solvent.in can not be opened.'
        READ (5,*) ! first line is made of comments
        READ (5,*) nb_solvent_sites
        allocate (solvent(1)%site(nb_solvent_sites), stat=ios)
            if (ios /= 0) stop "ERROR: wrong allocate of solvent%site in read_solvent.f90"
        READ(5,*) ! comment line
        DO n = 1 , size(solvent(1)%site)
            READ(5,*) ios, solvent(1)%site(n)%q, solvent(1)%site(n)%sig, solvent(1)%site(n)%eps, solvent(1)%site(n)%r
        END DO
    CLOSE(5)

    !... compute monopole, dipole, quadrupole, octupole and hexadecapole of each solvent species
    !... 1 Debye (D)  = 3.33564095 x10-30 C·m (= -0.20819435 e-·Å) 
    do concurrent (s=1:size(solvent))
        !... monopole = net charge
        solvent(s)%monopole = chop(sum( solvent(s)%site%q ))

        !... dipole
        do concurrent (i=1:3)
            solvent(s)%dipole(i) = chop(sum( solvent(s)%site%q * solvent(s)%site%r(i) ))
        end do

        !... quadrupole
        do concurrent (i=1:3, j=1:3)
            solvent(s)%quadrupole(i,j) = chop(sum( solvent(s)%site%q * solvent(s)%site%r(i) * solvent(s)%site%r(j) ))
        end do

        !... octupole
        do concurrent (i=1:3, j=1:3, k=1:3)
            solvent(s)%octupole(i,j,k) = chop(&
                sum( solvent(s)%site%q * solvent(s)%site%r(i) * solvent(s)%site%r(j) * solvent(s)%site%r(k) ))
        end do

        !... hexadecapole
        do concurrent (i=1:3, j=1:3, k=1:3, l=1:3)
            solvent(s)%hexadecapole(i,j,k,l) = chop( sum( solvent(s)%site%q * &
                solvent(s)%site%r(i) * solvent(s)%site%r(j) * solvent(s)%site%r(k) * solvent(s)%site%r(l) ))
        end do
    end do

END SUBROUTINE read_solvent
!===================================================================================================================================

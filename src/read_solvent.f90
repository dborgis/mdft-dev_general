!===================================================================================================================================
SUBROUTINE read_solvent
!===================================================================================================================================
! Read solvent atomic positions, charge, and lennard jones values in solvent.in
! charge in electron units, sigma in Angstroms, epsilon in KJ/mol.

    USE precision_kinds, ONLY: i2b , dp
    USE system, ONLY: nb_solvent_sites, solventSite, solvent
    IMPLICIT NONE

    INTEGER(i2b) :: n, ios
    
    OPEN(5, FILE= 'input/solvent.in', STATUS= 'old', IOSTAT= ios )! open input/solvent.in and check if it is readable
        IF ( ios/=0 ) STOP 'ERROR: solvent.in can not be opened.'
        READ (5,*) ! pass first line that are comments
        READ (5,*) nb_solvent_sites
        ALLOCATE( solventSite (nb_solvent_sites), stat=ios)
        if (ios /= 0) stop "ERROR: wrong allocate of solventSite in read_solvent.f90"
        allocate (solvent(1)%site(nb_solvent_sites), stat=ios)
        if (ios /= 0) stop "ERROR: wrong allocate of solvent%site in read_solvent.f90"
        ! read rest of the file
        ! check the total charge of the solvent to WARN USER if it is not neutral
        ! init total charge of the solvent to zero
        READ(5,*) ! comment line
        DO n = 1 , nb_solvent_sites
            READ(5,*) ios, solventSite(n)%q, solventSite(n)%sig, solventSite(n)%eps, solventSite(n)%r
            solvent(1)%site(n)%q   = solventSite(n)%q
            solvent(1)%site(n)%sig = solventSite(n)%sig
            solvent(1)%site(n)%eps = solventSite(n)%eps
            solvent(1)%site(n)%r   = solventSite(n)%r
        END DO
    CLOSE(5)

    IF( SUM(solventSite%q) /= 0._dp ) PRINT*,"WARNING: your solvent has a net charge of",SUM(solventSite%q)

END SUBROUTINE read_solvent
!===================================================================================================================================

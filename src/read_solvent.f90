!===================================================================================================================================
SUBROUTINE read_solvent
!===================================================================================================================================
! Read solvent atomic positions, charge, and lennard jones values in solvent.in
! charge in electron units, sigma in Angstroms, epsilon in KJ/mol.

    USE precision_kinds, ONLY: i2b , dp
    USE system, ONLY: nb_solvent_sites, solventSite
    IMPLICIT NONE

    INTEGER(i2b) :: i, n, ios
    
    OPEN(5, FILE= 'input/solvent.in', STATUS= 'old', IOSTAT= ios )! open input/solvent.in and check if it is readable
        IF ( ios/=0 ) STOP 'solvent.in cannot be opened ! => STOP !'
        READ (5,*) ! pass first line of comments
        READ (5,*) nb_solvent_sites
        ALLOCATE( solventSite (nb_solvent_sites) )
        ! read rest of the file
        ! check the total charge of the solvent to WARN USER if it is not neutral
        ! init total charge of the solvent to zero
        READ(5,*) ! comment line
        DO n = 1 , nb_solvent_sites
            READ(5,*) i, solventSite(n)%q, solventSite(n)%sig, solventSite(n)%eps, solventSite(n)%r
        END DO
        ! close input/solvent.in
    CLOSE(5)

    IF( SUM(solventSite%q) /= 0._dp ) PRINT*,"YOUR SOLVENT WEAR A TOTAL CHARGE OF ",SUM(solventSite%q),". YOU ARE NOW AWARE..."

END SUBROUTINE read_solvent
!===================================================================================================================================

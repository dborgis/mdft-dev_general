!> read solvent atomic positions, charge, and lennard jones values in solvent.in
!! charge in electron units, sigma in Angstroms, epsilon in KJ/mol.
SUBROUTINE read_solvent

    USE precision_kinds, ONLY: i2b , dp
    USE system, ONLY: nb_solvent_sites, chg_solv , sig_solv , eps_solv , id_solv, solventSite
    IMPLICIT NONE

    INTEGER(i2b) :: i, n, ios
    INTEGER(i2b) :: nb_id_solv ! number of types of solvent sites
    
    OPEN(5, FILE= 'input/solvent.in', STATUS= 'old', IOSTAT= ios )! open input/solvent.in and check if it is readable
        IF ( ios/=0 ) STOP 'solvent.in cannot be opened ! => STOP !'
        READ (5,*) ! pass first line of comments
        READ (5,*) nb_solvent_sites, nb_id_solv ! Second line is the total number of atom sites of the solvent AND the total number of different types of atoms
        ALLOCATE( solventSite (nb_solvent_sites) )
        ALLOCATE ( id_solv ( nb_solvent_sites ) ) ! from solvent_site to id for instance id_solv(1)=1, id_solv(2)=2 and id_solv(3)=2 for OH2
        ! read rest of the file
        ! check the total charge of the solvent to WARN USER if it is not neutral
        ! init total charge of the solvent to zero
        READ(5,*) ! comment line
        DO n = 1 , nb_solvent_sites
            READ(5,*) id_solv(n), solventSite(n)%q, solventSite(n)%sig, solventSite(n)%eps, solventSite(n)%r
        END DO
        ! close input/solvent.in
    CLOSE(5)

    ALLOCATE ( chg_solv ( nb_id_solv ) )
    ALLOCATE ( sig_solv ( nb_id_solv ) )
    ALLOCATE ( eps_solv ( nb_id_solv ) )
    DO n = 1, SIZE(solventSite)
        i = id_solv(n)
        chg_solv(i) = solventSite(n)%q
        sig_solv(i) = solventSite(n)%sig
        eps_solv(i) = solventSite(n)%eps
    END DO

    IF( SUM(solventSite%q) /= 0._dp ) PRINT*,"YOUR SOLVENT WEAR A TOTAL CHARGE OF ",SUM(solventSite%q),". YOU ARE NOW AWARE..."

END SUBROUTINE read_solvent

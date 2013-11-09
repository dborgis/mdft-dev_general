!> read solvent atomic positions, charge, and lennard jones values in solvent.in
!! charge in electron units, sigma in Angstroms, epsilon in KJ/mol.
SUBROUTINE read_solvent

    USE precision_kinds, ONLY: i2b , dp
    USE system, ONLY: nb_solvent_sites , x_solv , y_solv , z_solv , chg_solv , sig_solv , eps_solv , id_solv, solventSite
    IMPLICIT NONE

    integer(i2b):: i, n
    integer(i2b):: stat ! status du fichier ouvert
    integer(i2b):: nb_id_solv ! number of types of solvent sites
    real(dp):: total_charge ! charge of the whole system (should be zero in almost every case)
    
    OPEN(5, file= 'input/solvent.in', status= 'old', iostat= stat )! open input/solvent.in and check if it is readable
        IF ( stat/=0 ) STOP 'solvent.in cannot be opened ! => STOP !'
        read (5,*) ! pass first line of comments
        read (5,*) nb_solvent_sites , nb_id_solv ! Second line is the total number of atom sites of the solvent AND the total number of different types of atoms
        
        ! allocate accordingly
        allocate ( x_solv ( nb_solvent_sites ) )
        allocate ( y_solv ( nb_solvent_sites ) )
        allocate ( z_solv ( nb_solvent_sites ) )
        allocate ( id_solv ( nb_solvent_sites ) ) ! from solvent_site to id for instance id_solv(1)=1, id_solv(2)=2 and id_solv(3)=2 for OH2
        allocate ( chg_solv ( nb_id_solv ) )
        allocate ( sig_solv ( nb_id_solv ) )
        allocate ( eps_solv ( nb_id_solv ) )
        
        ! read rest of the file
        ! check the total charge of the solvent to WARN USER if it is not neutral
        ! init total charge of the solvent to zero
        total_charge = 0.0_dp
        read(5,*) ! comment line
        do n = 1 , nb_solvent_sites
            read (5,*) id_solv ( n ) , chg_solv ( id_solv ( n ) ) , sig_solv ( id_solv ( n ) ) , &
                    eps_solv ( id_solv ( n ) ) , x_solv ( n ) , y_solv ( n ) , z_solv ( n )
            total_charge = total_charge + chg_solv ( id_solv ( n ) )
        end do
        ! close input/solvent.in
    CLOSE(5)

    ALLOCATE( solventSite (nb_solvent_sites) )
    DO n = 1, nb_solvent_sites
        i = id_solv(n)
        solventSite(n)%q = chg_solv(i)
        solventSite(n)%sig = sig_solv(i)
        solventSite(n)%eps = eps_solv(i)
        solventSite(n)%r(1) = x_solv(n)
        solventSite(n)%r(2) = y_solv(n)
        solventSite(n)%r(3) = z_solv(n)
    END DO

    IF( SUM(solventSite%q) /= 0._dp ) PRINT*,"YOUR SOLVENT WEAR A TOTAL CHARGE OF ",total_charge,". BE SURE YOU KNOW WHAT YOU'RE DOING."

END SUBROUTINE read_solvent

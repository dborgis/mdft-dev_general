!> read solute atomic positions, charge, and lennard jones values in solute.in
!! charge in electron units, sigma in Angstroms, epsilon in KJ/mol.
SUBROUTINE read_solute

USE precision_kinds,only : i2b,dp
use system,only : nb_solute_sites , x_mol , y_mol , z_mol , chg_mol , sig_mol , eps_mol , atomic_nbr , id_mol , Lx , Ly , Lz&
&, lambda1_mol , lambda2_mol, soluteSite
use input,only : input_line
use periodic_table,only : init_periodic_table , ptable

IMPLICIT NONE

integer(i2b):: n, i
integer(i2b):: stat ! status du fichier ouvert
integer(i2b):: nb_id_mol ! number of different kinds of site (ie two LJ sites with same epsilon and sigma, but only diff positions)
! init module periodic_table so that all informations are available

    CALL init_periodic_table
    ! print *, ptable ( 1 ) % name
    ! open and test if input/solute.in is ok

    OPEN (5, FILE='input/solute.in', STATUS='old', IOSTAT=stat)
        IF (stat /= 0) THEN
            PRINT*,'solute.in cannot be opened ! => STOP !'
            STOP
        END IF

        READ (5,*) ! comment line
        READ (5,*) nb_solute_sites, nb_id_mol ! total number of atom sites of the solute AND the total number of different types of atoms
        ALLOCATE(soluteSite(nb_solute_sites))
        ALLOCATE(x_mol(nb_solute_sites))
        ALLOCATE(y_mol(nb_solute_sites))
        ALLOCATE(z_mol(nb_solute_sites))
        ALLOCATE(id_mol(nb_solute_sites)) ! from solute_site to id for instance id_mol(1)=1, id_mol(2)=2 and id_mol(3)=2 for OH2
        ALLOCATE(atomic_nbr(nb_solute_sites))
        ALLOCATE(chg_mol(nb_id_mol))
        ALLOCATE(sig_mol(nb_id_mol))
        ALLOCATE(eps_mol(nb_id_mol))
        ALLOCATE(lambda1_mol(nb_solute_sites))
        ALLOCATE(lambda2_mol(nb_solute_sites))

        READ (5,*) ! comment line
            DO n = 1, SIZE(soluteSite)
                READ(5,*) id_mol(n), soluteSite(n)%q, soluteSite(n)%sig, soluteSite(n)%eps, &
                        soluteSite(n)%lambda1, soluteSite(n)%lambda2, soluteSite(n)%r, soluteSite(n)%Z
!~ READ(5,*) id_mol(n), chg_mol(id_mol(n)) , sig_mol(id_mol(n)) , eps_mol(id_mol(n)) ,lambda1_mol(n) , lambda2_mol(n), x_mol(n) , y_mol(n) , z_mol(n) ,atomic_nbr(n)
            END DO
        CLOSE (5)

    CALL translate_to_center_of_supercell_if_needed ! if user wants all the sites to be translated to the center of the box, ie by Lx/2, Ly/2, Lz/2
    CALL print_supercell_xsf ! Print periodic XSF file to be read by VMD or equivalent
    CALL assure_coo_inside_cell ! check if cartesian coordinates read in input/solute.in are in the supercell



    ! As a first step toward removing all x_mol etc, I make them as pointers to our new derived type
    x_mol = soluteSite%r(1)
    y_mol = soluteSite%r(2)
    z_mol = soluteSite%r(3)
    atomic_nbr = soluteSite%Z
    lambda1_mol = soluteSite%lambda1
    lambda2_mol = soluteSite%lambda2
    DO n = 1, nb_solute_sites
        i = id_mol(n)
        chg_mol(i) = soluteSite(n)%q
        sig_mol(i) = soluteSite(n)%sig
        eps_mol(i) = soluteSite(n)%eps
    END DO


    CONTAINS

    ! if user asks for it (tag 'translate_solute_to_center'), add Lx/2, Ly/2, Lz/2 to all solute coordinates
    SUBROUTINE translate_to_center_of_supercell_if_needed
        USE input, ONLY: input_log
        USE system, ONLY: spaceGrid
        IF (input_log( 'translate_solute_to_center' )) THEN
            soluteSite%r(1) = soluteSite%r(1) + spaceGrid%length(1)/2.0_dp
            soluteSite%r(2) = soluteSite%r(2) + spaceGrid%length(2)/2.0_dp
            soluteSite%r(3) = soluteSite%r(3) + spaceGrid%length(3)/2.0_dp
        END IF
    END SUBROUTINE translate_to_center_of_supercell_if_needed



    SUBROUTINE assure_coo_inside_cell
        USE SYSTEM, ONLY: spaceGrid
        INTEGER(i2b) :: i 
        ! check if some positions are out of the supercell
        !j is a test tag. We loop over this test until every atom is in the box.
        ! This allows for instance, if a site is two boxes too far to still be ok.
        DO CONCURRENT ( i=1:SIZE(soluteSite) )
            soluteSite(i)%r(1) = MODULO ( soluteSite(i)%r(1) , spaceGrid%length(1) )
            soluteSite(i)%r(2) = MODULO ( soluteSite(i)%r(2) , spaceGrid%length(2) )
            soluteSite(i)%r(3) = MODULO ( soluteSite(i)%r(3) , spaceGrid%length(3) )
        END DO
    END SUBROUTINE assure_coo_inside_cell


END SUBROUTINE read_solute

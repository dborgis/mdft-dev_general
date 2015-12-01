! This subroutine writes an intent(in) array called 'array' to 'filename' in the file format .cube that handles 3-dimensional data.
! The standard .cube file can be read by vmd for instance. It is meant to vizualize periodic supercells.
! see for instance http://www.ks.uiuc.edu/Research/vmd/plugins/molfile/cubeplugin.html
! or http://paulbourke.net/dataformats/cube/
! It is a very easy format to deal with because it's very natural to write and read.
! Inputs needed : array to write, filename, number of points in each direction, length of each supercell vector, atomic number of the
! atoms and their coordinates. The total number of atoms is determined as the number of different coordinates ie as size(solute%site)
! output : a file called 'filename' containing 3dimensional data in .cube format
! written by Maximilien Levesque, while in post doc at Ecole Normale Superieure, Paris in Daniel Borgis's theoretical chemistry group
! 20110919  Maximilien Levesque, clean version for Virginie M.
SUBROUTINE write_to_cube_file ( array , filename )

    use module_solute, only: solute
    use module_grid, only: grid
    use precision_kinds, ONLY: i2b, dp

    IMPLICIT NONE

    INTEGER(i2b) :: i , j , k , n ! dummy variables
    CHARACTER(50), INTENT(IN) :: filename ! filename of .cube file. For example : "density.cube"
    REAL(dp), INTENT(IN), DIMENSION(grid%n_nodes(1), grid%n_nodes(2),grid%n_nodes(3)) :: array
    REAL(dp), PARAMETER :: angtobohr = 1.889725989_dp ! 1 Bohr = 1.889725989 Ang. Necessary because of VMD understanding of lengths

    ! define formats for writing to file

    ! first open the file you want to print array in. It's a formatted file.
    open ( 10 , file = filename , form = 'formatted' )
        write ( 10 , * ) ' CPMD CUBE FILE.' ! default text
        write ( 10 , * ) ' OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z' ! default text, for remembering
        ! write the total number of sites and a default text nobody knows it meaning
        write ( 10 , * ) SIZE(solute%site) ,' 0.0 0.0 0.0 '
        ! write primary vectors
        write ( 10 , * ) grid%n_nodes(1),grid%dl(1)*angtobohr,' 0.0 0.0'
        write ( 10 , * ) grid%n_nodes(2),' 0.0 ',grid%dl(2)*angtobohr, ' 0.0'
        write ( 10 , * ) grid%n_nodes(3),' 0.0 0.0 ',grid%dl(3)*angtobohr
        ! write the atoms and their coordinates in Bohr
        DO n = 1, SIZE(solute%site)! , dim = 1 )
            WRITE(10,*) solute%site(n)%Z, ' 0.0 ' , &
                    solute%site(n)%r(1)*angtobohr, solute%site(n)%r(2)*angtobohr, solute%site(n)%r(3)*angtobohr
        END DO
        ! write .cube file. One value per line. As said earlier, run over x then y then z. The one varying the most rapidly is z.
        do i = 1 , grid%n_nodes(1)
            do j = 1 , grid%n_nodes(2)
                do k = 1 , grid%n_nodes(3)
                    WRITE(10,*) real(array(i,j,k))
                end do
            end do
        end do
    close(10)
end subroutine

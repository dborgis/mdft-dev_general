! This subroutine writes an intent(in) array called 'array' to 'filename' in the file format .cube that handles 3-dimensional data.
! The standard .cube file can be read by vmd for instance. It is meant to vizualize periodic supercells.
! see for instance http://www.ks.uiuc.edu/Research/vmd/plugins/molfile/cubeplugin.html
! or http://paulbourke.net/dataformats/cube/
! It is a very easy format to deal with because it's very natural to write and read.
! Inputs needed : array to write, filename, number of points in each direction, length of each supercell vector, atomic number of the
! atoms and their coordinates. The total number of atoms is determined as the number of different coordinates ie as size(x_mol)
! output : a file called 'filename' containing 3dimensional data in .cube format
! written by Maximilien Levesque, while in post doc at Ecole Normale Superieure, Paris in Daniel Borgis's theoretical chemistry group
! 20110919  Maximilien Levesque, clean version for Virginie M.
SUBROUTINE write_to_cube_file ( array , filename )
    
    USE system, ONLY: spaceGrid, soluteSite
    USE precision_kinds, ONLY: i2b, dp
    
    IMPLICIT NONE
    
    INTEGER(i2b) :: i , j , k , n ! dummy variables
    CHARACTER(50), INTENT(IN) :: filename ! filename of .cube file. For example : "density.cube"
    REAL(dp), INTENT(IN), DIMENSION(spaceGrid%n_nodes(1), spaceGrid%n_nodes(2),spaceGrid%n_nodes(3)) :: array
    REAL(dp), PARAMETER :: angtobohr = 1.889725989_dp ! 1 Bohr = 1.889725989 Ang. Necessary because of VMD understanding of lengths

    ! define formats for writing to file
    104 FORMAT ( xI3 , xA , 3(xxF10.5) )
    102 FORMAT ( 1(xF10.5) )

    ! first open the file you want to print array in. It's a formatted file.
    open ( 10 , file = filename , form = 'formatted' )
        write ( 10 , * ) ' CPMD CUBE FILE.' ! default text
        write ( 10 , * ) ' OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z' ! default text, for remembering
        ! write the total number of sites and a default text nobody knows it meaning
        write ( 10 , * ) SIZE(soluteSite) ,' 0.0 0.0 0.0 '
        ! write primary vectors
        write ( 10 , * ) spaceGrid%n_nodes(1),spaceGrid%dl(1)*angtobohr,' 0.0 0.0'
        write ( 10 , * ) spaceGrid%n_nodes(2),' 0.0 ',spaceGrid%dl(2)*angtobohr, ' 0.0'
        write ( 10 , * ) spaceGrid%n_nodes(3),' 0.0 0.0 ',spaceGrid%dl(3)*angtobohr
        ! write the atoms and their coordinates in Bohr
        DO n = 1, SIZE(soluteSite)! , dim = 1 )
            WRITE(10,104) soluteSite(n)%Z, ' 0.0 ' , &
                    soluteSite(n)%r(1)*angtobohr, soluteSite(n)%r(2)*angtobohr, soluteSite(n)%r(3)*angtobohr
        END DO
        ! write .cube file. One value per line. As said earlier, run over x then y then z. The one varying the most rapidly is z.
        do i = 1 , spaceGrid%n_nodes(1)
            do j = 1 , spaceGrid%n_nodes(2)
                do k = 1 , spaceGrid%n_nodes(3)
                    WRITE(10,102) array(i,j,k)
                end do
            end do
        end do
    close(10)
end subroutine

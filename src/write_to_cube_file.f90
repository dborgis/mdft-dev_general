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

subroutine write_to_cube_file ( array , filename )

    use module_solute, only: solute
    use module_grid, only: grid
    use precision_kinds, only: dp

    implicit none

    integer :: i, j, k, n ! dummy variables
    character(50), intent(in) :: filename ! filename of .cube file. For example : "density.cube"
    real(dp), intent(in), dimension(grid%nx, grid%ny, grid%nz) :: array
    real(dp), parameter :: angtobohr = 1.889725989_dp ! 1 Bohr = 1.889725989 Ang. Necessary because of VMD understanding of lengths

    ! define formats for writing to file

    ! first open the file you want to print array in. It's a formatted file.
    open( 10 , file = filename , form = 'formatted' )
        write(10,*) ' CPMD CUBE FILE.' ! default text
        write(10,*) ' OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z' ! default text, for remembering
        ! write the total number of sites and a default text nobody knows it meaning
        write(10,*) size(solute%site) ,' 0.0 0.0 0.0 '
        ! write primary vectors
        write(10,*) grid%nx, grid%dx*angtobohr,' 0.0 0.0'
        write(10,*) grid%ny,' 0.0 ', grid%dy*angtobohr, ' 0.0'
        write(10,*) grid%nz,' 0.0 0.0 ', grid%dz*angtobohr
        ! write the atoms and their coordinates in Bohr
        do n = 1, size(solute%site)! , dim = 1 )
            write(10,*) solute%site(n)%Z, ' 0.0 ' , &
              solute%site(n)%r(1)*angtobohr, solute%site(n)%r(2)*angtobohr, solute%site(n)%r(3)*angtobohr
        end do
        ! write .cube file. One value per line. As said earlier, run over x then y then z. The one varying the most rapidly is z.
        do i=1, grid%nx
            do j=1, grid%ny
                do k=1, grid%nz
                    write(10,*) real(array(i,j,k))
                end do
            end do
        end do
    close(10)

end subroutine

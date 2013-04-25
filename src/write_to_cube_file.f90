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

subroutine write_to_cube_file ( array , filename )

use system , only: nfft1 , nfft2 , nfft3 , Lx , Ly , Lz , atomic_nbr , x_mol , y_mol , z_mol

! nfft1, nfft2 and nfft3 are the number of points (resolution) in directions x y and z

! Lx, Ly and Lz are the length of the supercell primary vectors. They are in Angstroms and thus need to be multiplied by angtobohr to be in bohr.

! atomic_nbr is the atomic number (periodic classification) of the atoms

! x_mol, y_mol and z_mol are the x, y and z coordinates of the atoms

use precision_kinds , only : dp ! defines the precision associated to simple and double

implicit none

integer :: i , j , k , n ! dummy variables

character ( 50 ) , intent(in) :: filename ! filename of .cube file. For example : "density.cube"

real(dp), intent(in) , dimension ( nfft1 , nfft2 , nfft3 ) :: array ! array printed in .cube file

real(dp), parameter :: angtobohr = 1.889725989_dp ! 1 Bohr = 1.889725989 Ang. Necessary because of VMD understanding of lengths

! define formats for writing to file



104 format ( xI3 , xA , 3(xxF10.5) )

102 format ( 1(xF10.5) )

! first open the file you want to print array in. It's a formatted file.

open ( 10 , file = filename , form = 'formatted' )

write ( 10 , * ) ' CPMD CUBE FILE.' ! default text

write ( 10 , * ) ' OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z' ! default text, for remembering

! write the total number of sites and a default text nobody knows it meaning

write ( 10 , * ) size ( x_mol ) ,' 0.0 0.0 0.0 ' ! 0 0 0 ou Lx/2 Ly/2 Lz/2 ?  Size(x_mol) is the total number of atoms

! write primary vectors

write ( 10 , * ) nfft1 , Lx / real ( nfft1 , kind = dp ) * angtobohr , ' 0.0 0.0'

write ( 10 , * ) nfft2 , ' 0.0 ' , Ly / real ( nfft2 , kind = dp ) * angtobohr , ' 0.0'

write ( 10 , * ) nfft3 , ' 0.0 0.0 ' , Lz / real ( nfft3 , kind = dp ) * angtobohr

! write the atoms and their coordinates in Bohr

do n = 1 , size ( x_mol , dim = 1 )

  write ( 10 , 104 ) atomic_nbr ( n ) , ' 0.0 ' , x_mol ( n ) * angtobohr , y_mol ( n ) * angtobohr , z_mol ( n ) * angtobohr

end do

! write .cube file. One value per line. As said earlier, run over x then y then z. The one varying the most rapidly is z.

do i = 1 , nfft1

  do j = 1 , nfft2
  
    do k = 1 , nfft3

      write ( 10 , 102 ) array ( i , j , k )

    end do

  end do

end do

! close the .cube file called filename

close(10)

! warn user

write ( * , * ) filename , ' written'

end subroutine

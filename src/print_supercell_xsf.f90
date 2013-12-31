!> Print an XSF file of the supercell for visualisation in VMD for instance.
!! Type vmd --xsf output/supercell.xsf to visualise it.
SUBROUTINE print_supercell_xsf
    USE precision_kinds, ONLY: i2b
    USE system, ONLY: spaceGrid, soluteSite
    IMPLICIT NONE
    integer(i2b) :: i
    
    open(5,file='output/supercell.xsf')
        100 format (xA)
        101 format (3(xF10.5))
        102 format (xI5,xI1)
        103 format (xI3,3(xxF10.5))
        write(5,100)'# this is the specification file of the supercell'
        write(5,100)'# lines beginning with # are commented. There cannot be comment lines within the sections'
        write(5,100)'# XSF format specifications can be found on the XCrySDen website http://www.xcrysden.org/doc/XSF.html'
        write(5,100)'# I strongly recommends to read this documentation.'
        write(5,*)
        write(5,100)'# for periodic structures one has to begin with word CRYSTAL'
        write(5,100)'CRYSTAL'
        write(5,100)
        write(5,100)'# Then one needs to specify the lattice vectors'
        write(5,100)'# specification of PRIMVEC (in ANGSTROMS) like:'
        write(5,100)'#         ax, ay, az    (first lattice vector)'
        write(5,100)'#         bx, by, bz    (second lattice vector)'
        write(5,100)'#         cx, cy, cz    (third lattice vector)'
        write(5,100)'# pay attention to vectors as they are written in horizontal way which is quite unusual'
        write(5,100)'# for now only orthorhombic structures are allowed (free norms of lattice vectors, all angles are 90 degrees)'
        write(5,100)'PRIMVEC'
        write(5,101) spaceGrid%length(1), 0., 0.
        write(5,101) 0., spaceGrid%length(2), 0.
        write(5,101) 0., 0., spaceGrid%length(3)
        write(5,*)
        write(5,100)'# Then one needs to specify the atoms belonging to the unit cell. '
        write(5,100)'# First number stands for number of atoms in the primitive cell (2 in this case).'
        write(5,100)'# The second number is always 1 for PRIMCOORD coordinates.'
        write(5,100)'# in angstroms and cartesian coordinates'
        write(5,100)'PRIMCOORD'
        write(5,102) SIZE(soluteSite), 1
        do i = 1, SIZE(soluteSite)
            write(5,103) soluteSite(i)%Z, soluteSite(i)%r
        END DO
    close(5)
END SUBROUTINE print_supercell_xsf

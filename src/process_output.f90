! This subroutinen calculates every output asked by user.
subroutine process_output

    ! defines precision of reals and intergers
    use precision_kinds , only: dp , i2b
    use system, only: nfft1 , nfft2 , nfft3 , nb_solute_sites , nb_species
    use input, only: input_line
    use solute_geometry, only: soluteIsPlanar => isPlanar, soluteIsLinear => isLinear
    
    implicit none
    
    real(dp), dimension( nfft1 , nfft2 , nfft3 , nb_species ) :: neq ! equilibrium density(position)
    real(dp), dimension( nfft1 , nfft2 , nfft3 , nb_species ) :: Px , Py , Pz ! equilibrium polarization(position)
    character(50):: filename
    real(dp), allocatable , dimension ( : , : , : , : ) :: temparray
    integer(i2b):: i , j ! dummy

    ! print output/density.bin which has the content of cg_vect 
    call print_cg_vect ! output/density.bin
    ! Get the final density (position) from the last minimizer step.
    call get_final_density ( neq )
    !> Write density in .cube file
    filename = 'output/density.cube'
    call write_to_cube_file ( neq ( : , : , : , 1 ) , filename ) ! TODO for now only write for the first species
    ! get the final polarization, if needed
    ! look for tag polarization in input
    do i = 1 , size ( input_line )
        j = len ( 'polarization' )
        if ( input_line (i) (1:j) == 'polarization' .and. input_line (i) (j+4:j+4) == 'T' ) then
            call get_final_polarization ( Px , Py , Pz )
            ! Write polarization in .cube file
            filename = 'output/polarization.cube'
            allocate ( temparray ( nfft1 , nfft2 , nfft3 , nb_species ) )
            temparray = sqrt ( Px ** 2 + Py ** 2 + Pz ** 2 )
            call write_to_cube_file ( temparray ( : , : , : , 1 ) , filename ) ! TODO for now only write for the first species
            deallocate ( temparray )
            exit
        end if
    end do

    !> Get radial distribution functions
    ! TODO for now it's only if nb_solute_sites is not too big because else mmalloc CRASH
    filename = 'output/g.rdf'
    call compute_rdf ( neq , filename )
    ! If calculation is for hard sphere fluid in presence of a hard wall compute profile perp wall
    ! TODO: DONT HAVE TIME TO WRITE THE TEST TODAY
    filename = 'output/z_density.dat'
    call compute_z_density ( neq ( : , : , : , 1 ) , filename ) ! TODO for now only write for the first species
    !> Check if solute is linear

    if( soluteIsLinear() ) then
        ! nothing for now
    end if
    
    if( soluteIsPlanar() ) then
        filename = 'output/planardensity.out'
        call compute_planar_density ( neq ( : , : , : , 1 ) , filename ) ! TODO for now only write for the first species
    end if

    contains

        subroutine print_cg_vect
!~             use cg , only : cg_vect
!~             if ( .not. allocated ( cg_vect ) ) then
!~                 print *, 'cg_vect is not allocated in subroutine print_cg_vect in process_output.f90. STOP.'
!~                 stop
!~             end if
!~             open( unit = 10 , file = 'output/density.bin' , form = 'unformatted' )
!~             write ( 10 ) cg_vect
!~             print *, ' output/density.bin                                 written'
        end subroutine print_cg_vect

end subroutine process_output

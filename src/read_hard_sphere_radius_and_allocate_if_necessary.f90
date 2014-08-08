!===================================================================================================================================
! Read hard sphere radius
!===================================================================================================================================
! Read hard sphere radius of every constituant of the fluid. If necessary.
! It begins by checking if it has already been allocated. If not it allocates it.
! Then, it reads every line of input_line which contains the inputs for the tag 'hard_sphere_radius"
! It then reads line after line the hard sphere radius of each constituant
!===================================================================================================================================
  SUBROUTINE read_hard_sphere_radius_and_allocate_if_necessary
    
    USE precision_kinds  ,ONLY: dp, i2b
    USE input            ,ONLY: input_line
    USE system           ,ONLY: n_0_multispec, thermocond, nb_species, spaceGrid, solventSite
    USE hardspheres      ,ONLY: hs
    
    IMPLICIT NONE
    
    INTEGER(i2b) :: i,j,s
    REAL(dp)     :: d_wca ! optimal diameter for hard spheres in the case of lennard jones perturbation as defined by Verlet and Weis, Phys Rev A 1972

    ! if allocation of radius has not been not earlier (no reason for now but perhaps later someone will want to implement it), allocate it
    IF (.NOT. ALLOCATED(hs) ) THEN ! init all radius to zero
        ALLOCATE( hs(nb_species) )
        hs%R = 0._dp
    END IF
    
    ! read it in dft.in
    ! Get hard sphere radius
    ! two possibilities : 1/ pure hard spheres : read radius in input/dft.in    2/ perturbated HS by a lennard jones
    ! next few lines is a test to know if we have
    ! Algo : one looks for line containing 'hard_sphere_radius'. Next nb_species lines contain the radius of each hard sphere species.
    ! If the radius is negative, the user means that the radius has to be calculated using the week chandler anderson model 
    ! by reading the lennard jones sigma and epsilon values accordingly to integer species between 1 and nb_species
    DO i= 1, SIZE(input_line)
        j = LEN( 'hard_sphere_radius' )
        IF ( input_line (i) (1:j) == 'hard_sphere_radius' ) THEN
            DO s = 1 , nb_species
                READ( input_line(i+s) ,*) hs(s)%R
                IF ( hs(s)%R <= 0.0_dp ) THEN ! check if one radius is negative, ie if one has to compute wca diameter
                    CALL compute_wca_diameter ( n_0_multispec(s) , thermocond%T, solventSite(s)%sig , solventSite(s)%eps , d_wca )
                    hs(s)%R = d_wca / 2.0_dp
                END IF
            END DO
            EXIT
        END IF
    END DO
    
    ! check if hard sphere radius is coherent with the grid
    BLOCK
        REAL(dp) :: dmin
        dmin = MINVAL ( spaceGrid%length/REAL(spaceGrid%n_nodes,dp) )
        DO s = 1, nb_species
            IF ( hs(s)%R <= dmin ) THEN
                PRINT*,'Radius of hard sphere ' , s , ' is ' , hs(s)%R ,&
                    ' and is smaller than the distance between two nfft grid nodes ' , dmin
                STOP
            END IF
        END DO
    END BLOCK
    
END SUBROUTINE read_hard_sphere_radius_and_allocate_if_necessary

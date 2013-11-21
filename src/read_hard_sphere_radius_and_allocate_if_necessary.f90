!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Read hard sphere radius
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Read hard sphere radius of every constituant of the fluid. If necessary.
! It begins by checking if it has already been allocated. If not it allocates it.
! Then, it reads every line of input_line which contains the inputs for the tag 'hard_sphere_radius"
! It then reads line after line the hard sphere radius of each constituant
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE read_hard_sphere_radius_and_allocate_if_necessary
    USE precision_kinds , only : dp , i2b
    use input , only : input_line
    use system , only : radius , n_0_multispec , temp , sig_solv , eps_solv , Lx , Ly , Lz , nfft1 , nfft2 , nfft3 , nb_species
    IMPLICIT NONE
    integer(i2b):: i , j , species ! dummy
    real(dp):: d_wca !>@var optimal diameter for hard spheres in the case of lennard jones perturbation as defined by Verlet and Weis, Phys Rev A 1972
    real(dp):: distance_between_grid_nodes ! minimum distance between two nodes of the nfft grid = min ( Lx/nfft1 , Ly/nfft2 , Lz/nfft3 )
    ! if allocation of radius has not been not earlier (no reason for now but perhaps later someone will want to implement it), allocate it
    if ( .not. allocated ( radius ) ) then
      allocate ( radius ( nb_species ) )
      radius = 0.0_dp
    END IF
    ! read it in dft.in
    ! Get hard sphere radius
    ! two possibilities : 1/ pure hard spheres : read radius in input/dft.in    2/ perturbated HS by a lennard jones
    ! next few lines is a test to know if we have
    ! Algo : one looks for line containing 'hard_sphere_radius'. Next nb_species lines contain the radius of each hard sphere species.
    ! If the radius is negative, the user means that the radius has to be calculated using the week chandler anderson model 
    ! by reading the lennard jones sigma and epsilon values accordingly to integer species between 1 and nb_species
    do i = 1 , size ( input_line )
      j = len ( 'hard_sphere_radius' )
      if ( input_line (i) (1:j) == 'hard_sphere_radius' ) then
        do species = 1 , nb_species
          read ( input_line ( i + species ) , * ) radius ( species )
          ! check if one radius is negative, ie if one has to compute wca diameter
          if ( radius ( species ) <= 0.0_dp ) then
            call compute_wca_diameter ( n_0_multispec ( species ) , temp, sig_solv ( species ) , eps_solv ( species ) , d_wca )
            radius ( species ) = d_wca / 2.0_dp
          END IF
          write (*,*) 'Hard sphere radius (' , species , ') = ' , radius ( species )
        END DO
        exit
      END IF
    END DO
    ! chheck if hard sphere radius is coherent with the grid as it should never be small than the distance between two grid nodes
    distance_between_grid_nodes = min ( Lx / real ( nfft1 , dp ) , Ly / real ( nfft2 , dp ) , Lz / real ( nfft3 / dp ) )
    do species = 1 , nb_species
      if ( radius ( species ) <= distance_between_grid_nodes ) then
        write (*,*) 'Radius of hard sphere ' , species , ' is ' , radius ( species ) ,&
 ' and is smaller than the distance between two nfft grid nodes ' , distance_between_grid_nodes
        stop
      END IF
    END DO
  END SUBROUTINE read_hard_sphere_radius_and_allocate_if_necessary

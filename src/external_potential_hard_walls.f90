! this subroutine computes the external potential created by several (0 to infty) hard walls
! The hard walls "plans" have coordinates read in input/dft.in
! For a brief description of the coordinates of a plan, see wikipedia: 
! Step 1/ we read how many walls are to be found in the supercell.
! Step 2/ we read their coordinates
! Step 3/ we read their thinkness (2*radius) in dft.in
! Step 4/ we compute the external potential created by the walls, which depends upon the fluid radius
subroutine external_potential_hard_walls
! precision_kinds defines the precision of the real, double, integer variables.
use precision_kinds , only : i2b , dp
! i2b = integer simple precision_kind
! dp = double precision real
! contains input/dft.in put in input_line
use input , only : input_line
! input_line (:) contains all lines in dft.in
use system , only : nfft1 , nfft2 , nfft3 , Lz , radius , nb_species , deltax , deltay , deltaz
! nfft1 = number of grid nodes in X direction
! Lz = Length (Angstroms) of the supercell in X direction
use external_potential , only : Vext_total
use quadrature, only: angGrid, molRotGrid
! Vext_total = external potential as used in the total free energy calculation
implicit none
integer(i2b):: i , j , k , wall ! dummy
real(dp):: dplan ! local distance between point M (grid point) and plan defining wall
integer(i2b):: species ! dummy between 1 and nb_species
integer(i2b):: number_of_hard_walls ! number of hard walls in the supercell
real (dp) , allocatable , dimension ( : ) :: thickness ! thickness of each wall (thickness = 2radius. Don't mix them up)
real (dp) , allocatable , dimension ( : , : ) :: normal_vec ! normal vector of each plan defining wall
real (dp) , allocatable , dimension ( : ) :: norm2_normal_vec ! norm of normal_vec
real (dp) , allocatable , dimension ( : , : ) :: OA ! A is a point of coordinates OA(1,2,3) which is in the plan normal to normal_vec
real (dp) , allocatable , dimension ( : ) :: dot_product_normal_vec_OA ! dummy
real (dp) , dimension ( 3 ) :: OM ! coordinates of grid points
real (dp) :: OMx , OMy , OMz ! projections of OM
real (dp) , parameter :: infty = huge ( 1.0_dp )
! read how many hard walls are wanted in the supercell
do i = 1 , size ( input_line )
  j = len ( 'hard_wall_number' )
  if ( input_line (i) (1:j) == 'hard_wall_number' ) then
    read ( input_line (i) (j+4:j+5) , * ) number_of_hard_walls
    if ( number_of_hard_walls == 0 ) return ! no hard wall contribution to external potential. go back to init_external_potential
!    call test_number_of_hard_walls ( number_of_hard_walls ) ! test if the total number of hard walls in the supercell is physical
    allocate ( thickness ( number_of_hard_walls ) ) ! thickness of each wall (thickness = 2radius. Don't mix them up)
    allocate ( normal_vec ( 3 , number_of_hard_walls ) ) ! normal vector of each plan defining wall
    allocate ( norm2_normal_vec ( number_of_hard_walls ) ) ! norm of normal_vec. dummy
    allocate ( OA ( 3 , number_of_hard_walls ) ) ! A is a point of coordinates OA(1,2,3) which is in the plan normal to normal_vec
    allocate ( dot_product_normal_vec_OA ( number_of_hard_walls ) ) ! dummy
    do wall = 1 , number_of_hard_walls
      read ( input_line ( i + wall ) , * ) thickness ( wall ) , normal_vec ( 1 , wall ) , normal_vec ( 2 , wall ) , & 
            normal_vec ( 3 , wall ) , OA ( 1 , wall ) , OA ( 2 , wall ) , OA ( 3 , wall ) ! for each wall, read thickness, normal vector coordinates and point in plan
      norm2_normal_vec ( wall ) = norm2 ( normal_vec ( : , wall ) ) ! pretabulate norm of normal vec
      dot_product_normal_vec_OA ( wall ) = dot_product ( - normal_vec ( : , wall ) , OA ( : , wall ) )
    end do
    exit ! out of loop over i
  end if
end do
! be sure Vext_total is allocated
if ( .not. allocated ( Vext_total ) ) then
  allocate ( Vext_total ( nfft1 , nfft2 , nfft3 , angGrid%n_angles , molRotGrid%n_angles , nb_species ) )
  Vext_total = 0.0_dp
end if
! be sure radius(:) (solvent radius) is already computed
if ( .not. allocated ( radius ) ) then
  write (*,*) 'radius is not allocated. critial stop in external_potential_hard_walls.f90'
  allocate ( radius ( nb_species ) ) 
  radius ( : ) = 0.0_dp
end if
! potential
do k = 1 , nfft3
  OMz = real ( k - 1 , dp ) * deltaz
do j = 1 , nfft2
  OMy = real ( j - 1 , dp ) * deltay
do i = 1 , nfft1
  OMx = real ( i - 1 , dp ) * deltax
  OM = (/ OMx , OMy , OMz /)
  do wall = 1 , number_of_hard_walls ! for each wall
    dplan = abs ( dot_product ( normal_vec ( : , wall ) , OM ) + dot_product_normal_vec_OA ( wall ) ) / norm2_normal_vec ( wall ) ! compute distance between grid point and plan
    do species = 1 , nb_species
      if ( dplan <= 0.5_dp * thickness ( wall ) + radius ( species ) ) Vext_total ( i , j , k , : , : , species ) = infty
    end do ! species
  end do ! walls
end do
end do
end do
contains
!===================================================================================================================================
  ! this subroutine tests if the number of hard walls is physical. For now, only 1 or 2 is implemented.
!===================================================================================================================================
  subroutine test_number_of_hard_walls ( number_of_hard_walls )
    implicit none
    integer , intent(in) :: number_of_hard_walls
    if ( number_of_hard_walls > 2 ) then
      print *, 'the number of hard walls wanted by the user should not be > 2.'
      print *, 'that is not implemented for now'
      stop
    else if ( number_of_hard_walls < 0 ) then
      print *, 'number of hard walls is negative. That is non-sense.'
      stop
    end if
  end subroutine test_number_of_hard_walls
!===================================================================================================================================
end subroutine external_potential_hard_walls

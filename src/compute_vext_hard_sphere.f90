!> Compute external potential in a basic way
!> Solvation of a hard sphere in a hard sphere fluid
!This has been reintroduced in the dev code in October 2018
module module_Vext_hard_sphere

implicit none
private
public :: compute_vext_hard_sphere


contains

  subroutine compute_vext_hard_sphere 
  
   use module_grid, only: grid
   use precision_kinds,only : i2b , dp
   use module_input,only : getinput
   use module_solute, only : solute
   use module_solvent, only : solvent
   use hardspheres ,ONLY: hs, compute_hard_spheres_parameters
  
   IMPLICIT NONE
   integer(i2b):: i , j , k !> @var dummy
   integer(i2b):: s ! dummy between 1 and solvent(1)%nspec
   integer(i2b):: solute_site ! dummy between 1 and
   real(dp):: x_grid , y_grid , z_grid !> @var coordinates of grid mesh nodes
   real(dp):: x_nm2 , y_nm2 , z_nm2 , r_nm2  ! distance between solute and grid point
   real(dp):: hard_sphere_solute_radius ! radius of the hard sphere solute. tag in dft.in
   real(dp):: sum_of_solute_and_solvent_radius ! sum of solute and solvent radius
   INTEGER(i2b) :: nx, ny, nz
   REAL(dp) :: lx, ly, lz, deltax, deltay, deltaz
  
       lx= grid%lx
       ly= grid%ly
       lz= grid%lz
       deltax= grid%dx
       deltay= grid%dy
       deltaz= grid%dz
       nx= grid%nx
       ny= grid%ny
       nz= grid%nz
  
  
   
   ! look for tag 'hard_sphere_solute_radius' which gives the radius of the sphere for the external potential
   hard_sphere_solute_radius=getinput%dp('radius_of_hard_sphere_solute',defaultvalue=0.0_dp)
   if  (hard_sphere_solute_radius==0.0) return


   !check if radius is allocated. if not it may be because you wand an hard sphere solute but no hardsphere solvent, in that case
   !radius=0, if you want an hard sphere solvent there is a bug
   if ( .not. allocated ( hs ) ) then
       call compute_hard_spheres_parameters
   END IF
  
  
   ! compute the sum of solute and solvent radius
   do solute_site = 1 , solute%nsite
     do s = 1 , solvent(1)%nspec
       sum_of_solute_and_solvent_radius = ( hard_sphere_solute_radius + HS(s)%R ) ** 2
       do k = 1 , nz
         z_grid = real( k - 1 , dp ) * deltaz
         z_nm2 = ( z_grid - solute%site(solute_site)%r(3) ) ** 2
         do j = 1 , ny
           y_grid = real ( j - 1 , dp ) * deltay
           y_nm2 = ( y_grid - solute%site(solute_site)%r(2) ) ** 2
           do i = 1 , nx
             x_grid = real ( i - 1 , dp ) * deltax
             x_nm2 = ( x_grid - solute%site(solute_site)%r(1) ) ** 2
             r_nm2 = x_nm2 + y_nm2 + z_nm2
             ! test if distance is smaller than sum of radius of solvent and solute
             if ( r_nm2 <= sum_of_solute_and_solvent_radius ) then
               ! we're inside the hard potential. vext is very high
               solvent(s)%vext(:,i,j,k) = huge ( 1.0_dp ) ! no dependancy over Omega here
             END IF
           END DO
         END DO
       END DO
     END DO
   END DO
  
   
   end subroutine compute_vext_hard_sphere 
 
 END module module_Vext_hard_sphere

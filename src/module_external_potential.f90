! Here is a module dedicated to the calculation and the use of the external potential
module external_potential
USE precision_kinds,only : dp, i2b
IMPLICIT NONE
!real(dp), allocatable, dimension(:,:,:,:) :: Vext ! external potential of the whole solute on a position and angular grid
!real(dp), allocatable, dimension(:,:,:,:) :: Vlj ! lennard jones potential of the whole solute on a position and angular grid
!real(dp), allocatable, dimension(:,:,:,:) :: Vcoul ! electrostatic potential of the whole solute on a position and angular grid
! we want the outer loop (slowest varying) to be over species, so nb_species is the last rank.
! the arrays are thus defined as  array ( nfft1 , nfft2 , nfft3 , angGrid%n_angles , nb_species )
! the most efficient loops are thus over
! species
!   omega
!     nfft3
!       nfft2
!         nfft1
real(dp), allocatable , dimension ( : , : , : , : , : , : ) :: Vext_total ! external potential as the sum of all external potentials (LJ + charge + ... )
real(dp), allocatable , dimension ( : , : , : , : , : , : ) :: Vext_lj
real(dp), allocatable , dimension ( : , : , : , : , : , : ) :: Vext_q ! ( nfft1 , nfft2 , nfft3 , angGrid%n_angles , nb_species )
real(dp), allocatable , dimension ( : , : , : ) :: V_c ! electrostatic potential calculated from poisson equation
real(dp), allocatable , dimension ( : , : , : , : , : ) :: Vext_hard ! hard potential
real(dp), allocatable , dimension ( : , : , : , : ) :: Vext_hard_core ! hard core potential based on van der walls radius
!Calculation of charge density
real( dp ) , allocatable , dimension ( : , : , : ) ::  q_charge
integer (i2b) , allocatable , dimension (: , : , : ) :: x_charge, y_charge, z_charge
integer (i2b) :: nb_of_interpolation
real (dp) :: Fcoul



contains 

    ! Deallocates everything properly
    SUBROUTINE deallocate_everything_external_potential
        IMPLICIT NONE
        if ( allocated ( Vext_total ) ) deallocate ( Vext_total )
        if ( allocated ( Vext_lj ) ) deallocate ( Vext_lj )
        if ( allocated ( Vext_q ) ) deallocate ( Vext_q )
        if ( allocated ( x_charge ) ) deallocate ( x_charge )
        if ( allocated ( y_charge ) ) deallocate ( y_charge )
        if ( allocated ( z_charge ) ) deallocate ( z_charge )
        if ( allocated ( q_charge ) ) deallocate ( q_charge )
    END SUBROUTINE deallocate_everything_external_potential

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! see 
    PURE SUBROUTINE vextdef0
        IMPLICIT NONE
    END SUBROUTINE vextdef0
    

end module external_potential

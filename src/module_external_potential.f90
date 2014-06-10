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
real(dp), allocatable , dimension (:,:,:,:) :: Vext_hard_core ! hard core potential
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

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! see manual/external_potentials.pdf
    SUBROUTINE vextdef0
        USE precision_kinds ,ONLY: dp
        USE system          ,ONLY: grd=>spaceGrid, nb_species
        USE quadrature      ,ONLY: angrd=>angGrid,&
                                   mrgrd=>molRotGrid
        IMPLICIT NONE
        REAL(dp), PARAMETER :: radius=1.0_dp 
        
        CALL allocate_vext_hard_core_if_necessary
        CALL disk
        CALL wall
        
        CONTAINS
            
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            
            SUBROUTINE allocate_vext_hard_core_if_necessary
                USE constants, ONLY: zero
                IF (.NOT. ALLOCATED(vext_hard_core)) THEN
                    ALLOCATE(   vext_hard_core( grd%n_nodes(1),grd%n_nodes(2),grd%n_nodes(3),nb_species)  ,SOURCE=zero )
                END IF
            END SUBROUTINE allocate_vext_hard_core_if_necessary
            
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            
            SUBROUTINE wall
                REAL(dp) :: d, distNode2Center, coo(2), l(2)
                INTEGER(i2b) :: i,j,n(2)
                l(:) = grd%length(1:2)
                n(:) = grd%n_nodes(1:2)
                coo = l/REAL(n,dp) + [radius,l(2)/2.]
                vext_hard_core(1,:,:,:) = HUGE(1.0_dp) ! line of wall at farther left
                ! find all nodes that are inside a disk of radius RADIUS. The disk is at the center of the plane (Lx,Ly)
                DO CONCURRENT (i=1:n(1),j=1:n(2))
                    d=ABS(REAL(j-1,dp)*l(2)/n(2)-coo(2))
                    IF (d>radius+l(1)/n(1) .AND. REAL(i,dp)*l(1)/n(1)<radius) vext_hard_core(i,j,:,:)=HUGE(1._dp)
                END DO
            END SUBROUTINE wall
            
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            
            SUBROUTINE disk
                REAL(dp) :: distNode2Center, coo(2), l(2)
                INTEGER(i2b) :: i,j,n(2)
                l(:) = grd%length(1:2)
                n(:) = grd%n_nodes(1:2)
                coo = l/2._dp! - [0,0]*l(1)/n(1) (if you want to translate it toward wall) ! The cylinder is at the center of the supercell plane {x,y}
                ! find all nodes that are inside a disk of radius RADIUS. The disk is at the center of the plane (Lx,Ly)
                DO i=1,n(1); DO j=1,n(2)!DO CONCURRENT (i=1:n(1),j=1:n(2))
                    distNode2Center = SQRT(SUM(    (    REAL([i,j]-1,dp)*l/n - coo    )**2    ))
                    IF (distNode2Center<=radius) vext_hard_core(i,j,:,:)=HUGE(1._dp)
                END DO; END DO
            END SUBROUTINE disk

    END SUBROUTINE vextdef0
    

end module external_potential

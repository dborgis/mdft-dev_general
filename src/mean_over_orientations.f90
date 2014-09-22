!> This SUBROUTINE do the legendre integration over all orientations of any array (density, vext, ...)
SUBROUTINE mean_over_orientations(arrayin,arrayout)
USE precision_kinds, only: dp,i2b
use system, only: spaceGrid
use quadrature, only: angGrid, molRotGrid
IMPLICIT NONE
real(dp), dimension(spaceGrid%n_nodes(1),spaceGrid%n_nodes(2),spaceGrid%n_nodes(3),angGrid%n_angles,molRotGrid%n_angles),&
        intent(in) :: arrayin ! input array
real(dp), dimension(spaceGrid%n_nodes(1),spaceGrid%n_nodes(2),spaceGrid%n_nodes(3)), intent(out) :: arrayout ! output array
integer(i2b) :: n , p
arrayout=0.0_dp
do n=1,angGrid%n_angles
    do p=1,molRotGrid%n_angles
        arrayout(:,:,:)=arrayout(:,:,:)+arrayin(:,:,:,n,p)*angGrid%weight(n)*molRotGrid%weight(p)
    END DO
END DO
END SUBROUTINE mean_over_orientations

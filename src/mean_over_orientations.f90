!> This SUBROUTINE do the legendre integration over all orientations of any array (density, vext, ...)
SUBROUTINE mean_over_orientations ( arrayin , arrayout )
    use precision_kinds, only: dp
    use module_grid, only: grid
    implicit none
    real(dp), intent(in) :: arrayin(:,:,:,:)  ! x, y, z, orientation
    real(dp), intent(out) :: arrayout(:,:,:)   ! <x, y, z>_orientations
    ! real(dp), dimension(grid%n_nodes(1),grid%n_nodes(2),grid%n_nodes(3),angGrid%n_angles,molRotGrid%n_angles),&
    !         intent(in) :: arrayin ! input array
    ! real(dp), dimension(grid%n_nodes(1),grid%n_nodes(2),grid%n_nodes(3)), intent(out) :: arrayout ! output array
    integer :: io

    arrayout = 0.0_dp
    do io = 1, GRID%no
        arrayout = arrayout + arrayin(:,:,:,io)*GRID%w(io)
    end do

    ! arrayout = sum( arrayin * grid%w(:))
    ! do n=1,angGrid%n_angles
    !     do p=1,molRotGrid%n_angles
    !         arrayout(:,:,:)=arrayout(:,:,:)+arrayin(:,:,:,n,p)*angGrid%weight(n)*molRotGrid%weight(p)
    !     END DO
    ! END DO

END SUBROUTINE mean_over_orientations

! SUBROUTINE get_final_polarization ( Px , Py , Pz )
!
!     use precision_kinds, only: dp, i2b
!     use module_solvent, only: solvent
!     use module_grid, only: grid
!
!     IMPLICIT NONE
!     INTEGER(i2b) :: i, j, k, io, s
!     REAL(dp) :: x, local_Px, local_Py, local_Pz
!     REAL(dp), dimension(:,:,:,:), intent(out) :: Px, Py, Pz ! equilibrium polarization(r)
!     integer :: nx, ny, nz, no, ns
!     real(dp), parameter :: zerodp = 0._dp
!
!
!     Px = zerodp
!     Py = zerodp
!     Pz = zerodp
!
!     DO s =1,solvent(1)%nspec
!         DO i =1,grid%nx
!             DO j =1,grid%ny
!                 DO k =1,grid%nz
!                     local_Px = 0.0_dp
!                     local_Py = 0.0_dp
!                     local_Pz = 0.0_dp
!                     DO io =1,grid%no
!                         x = solvent(s)%density(i,j,k,io)
!                         local_Px = local_Px + grid%omx(io) * grid%w(io) * x
!                         local_Py = local_Py + grid%omy(io) * grid%w(io) * x
!                         local_Pz = local_Pz + grid%omz(io) * grid%w(io) * x
!                     END DO
!                     Px(i,j,k,s) = local_Px
!                     Py(i,j,k,s) = local_Py
!                     Pz(i,j,k,s) = local_Pz
!                 END DO
!             END DO
!         END DO
!     END DO
!
! END SUBROUTINE get_final_polarization

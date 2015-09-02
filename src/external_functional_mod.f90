!module external_functional_mod
!    use iso_c_binding, only: dp => c_double
!    implicit none
!    private
!    public :: get
!contains
!    subroutine get(vext,gr,den,fext)
!        use external_potential_mod, only: external_potential
!        use grid_mod, only: grid
!        use density_mod, only: density
!        use functional_mod, only: functional
!        type(external_potential), intent(in) :: vext
!        type(grid), intent(in) :: gr
!        type(density), intent(in) :: den
!        type(functional), intent(out) :: fext
!        integer :: io, ix, iy, iz
!        fext%energy = 0._dp
!        fext%grad = 0._dp
!        do iz = 1, grid%nz
!            do iy = 1, grid%ny
!                do ix = 1, grid%nx
!                    do io = 1, grid%no
!                        fext%energy = fext%energy + den%rho(io,ix,iy,iz) * grid%tw(io) * grid%dv * vext%tot(io,ix,iy,iz)
!                        fext%grad(io,ix,iy,iz) = grid%tw(io) * vext%tot(io,ix,iy,iz)
!                    end do
!                end do
!            end do
!        end do
!    end subroutine
!end module

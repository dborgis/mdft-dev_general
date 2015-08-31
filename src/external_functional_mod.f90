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
!        do iz = 1, gr%nz
!            do iy = 1, gr%ny
!                do ix = 1, gr%nx
!                    do io = 1, gr%no
!                        fext%energy = fext%energy + den%rho(io,ix,iy,iz) * gr%tw(io) * gr%dv * vext%tot(io,ix,iy,iz)
!                        fext%grad(io,ix,iy,iz) = gr%tw(io) * vext%tot(io,ix,iy,iz)
!                    end do
!                end do
!            end do
!        end do
!    end subroutine
!end module

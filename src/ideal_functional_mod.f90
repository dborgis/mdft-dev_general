module ideal_functional_mod
    use iso_c_binding, only: dp => c_double
    use functional_mod, only: functional
    use grid_mod, only: grid
    use density_mod, only: density
    implicit none
    private
    public :: get
contains
    subroutine get(gr,den,fid)
        use input, only: input_dp
        implicit none
        type(grid), intent(in) :: gr
        type(density), intent(in) :: den
        type(functional), intent(out) :: fid
        integer :: ix, iy, iz, io
        real(dp) :: rho, wo, rho_b, kT
        real(dp), parameter :: eps = epsilon(1._dp)
        real(dp), parameter :: boltz = 1.3806488e-23_dp
        real(dp), parameter :: navo = 6.02214129e23_dp
        kT = input_dp("temperature", defaultvalue=300._dp) * Boltz * Navo /1000._dp
        if (.not. fid%is_built) call fid%build(gr)
        fid%energy = 0._dp
        fid%grad = 0._dp
        rho_b = den%bulk
        do concurrent (io=1:gr%no, ix=1:gr%nx, iy=1:gr%ny, iz=1:gr%nz)
            rho = den%rho(io,ix,iy,iz)
            if (rho < eps) then
                fid%energy = gr%tw(io)*rho_b
                fid%grad(io,ix,iy,iz) = 0._dp
            else
                fid%energy = fid%energy + rho*log(rho/rho_b)-rho+rho_b
                fid%grad(io,ix,iy,iz) = gr%tw(io)*log(rho/rho_b)
                print*, fid%grad(io,ix,iy,iz)
            end if
        end do
        fid%energy = fid%energy * kT * gr%dv
        fid%grad   = fid%grad   * kT * gr%dv
    end subroutine
end module

module module_energy_ideal_and_external
    implicit none
    private
    public :: energy_ideal_and_external
contains

    ! This SUBROUTINE computes the ideal part of the free energy functional.
    subroutine energy_ideal_and_external (fid, fext, df)

        use precision_kinds, only: dp
        use module_thermo, only: thermo
        use module_solvent, only: solvent, print_solvent_not_allocated
        use module_grid, only: grid
        use module_input, only: getinput

        implicit none

        real(dp), intent(out) :: fid, fext
        real(dp), intent(out) :: df(:,:,:,:,:)
        integer :: is, io, ns, ix, iy, iz
        real(dp) :: x, x0, vext, dv, kT, mu, w, dfid
        real(dp), parameter :: zerodp = 0._dp
        real(dp), parameter :: epsdp = epsilon(0._dp)

        if (.not.allocated(solvent)) call print_solvent_not_allocated ("Look at subroutine energy_ideal_and_external")

        ns = solvent(1)%nspec
        kT = thermo%kbT
        dv = grid%dv
        mu = getinput%dp( 'imposed_chempot', defaultvalue=0._dp)
        if ( mu/=0._dp) stop "mu /=0 in module_energy_ideal_and_external. That's implemented but for now shutitoff"
        if ( ns/=1 .AND. mu/=0._dp) STOP "Imposing a chemical potential is valid only for single-species solvent"
        !
        ! fid = integrate over whole space of   x.log(x/x0)-x+x0 = Int[x.(log(x/x0)-1)] + Int[x0]
        !
        ! fext = Int[x*(vext-mu)]
        !
        fid = zerodp
        fext = zerodp
        df = zerodp
        do is=1,ns
            x0 = solvent(is)%rho0 ! bulk density
            do iz=1,grid%nz
                do iy=1,grid%ny
                    do ix=1,grid%nx
                        do io=1,grid%no
                            x = solvent(is)%density(io,ix,iy,iz)
                            if (x>epsdp) then
                                fid = fid + kT*dv*grid%w(io)*(x*log(x/x0)-x+x0)
                                dfid = kT*dv*grid%w(io)*log(x/x0)
                            else if (x<=epsdp) then
                                fid = fid + kT*dv*grid%w(io)*x0 ! par continuite. Plot x*log(x/x0)-x+x0 dans Mathematica pour un x non nul mais arbitrairement petit pour t'en convaincre.
                                dfid = 0._dp
                            end if
                            fext = fext + dv*grid%w(io)*(solvent(is)%vext(io,ix,iy,iz) - mu)*x
                            df(io,ix,iy,iz,is) = df(io,ix,iy,iz,is) + dfid + dv*grid%w(io)*(solvent(is)%vext(io,ix,iy,iz) - mu)
                        end do
                    end do
                end do
            end do
        end do

    end subroutine energy_ideal_and_external
end module module_energy_ideal_and_external

module module_energy_ideal_and_external

    implicit none
    private
    public :: energy_ideal_and_external

contains

    !
    ! This subroutine computes the ideal and external parts of the free energy functional.
    ! We compute Fid and Fext together because they're one-body potential related.
    !

    subroutine energy_ideal_and_external (fid, fext, df)

        use precision_kinds, only: dp
        use module_thermo, only: thermo
        use module_solvent, only: solvent, print_solvent_not_allocated
        use module_grid, only: grid
        use module_input, only: getinput

        implicit none

        real(dp), intent(out) :: fid, fext
        real(dp), intent(inout) :: df(:,:,:,:,:)
        integer :: is, io, ns, ix, iy, iz
        real(dp) :: xi, rho0, dv, kT, mu, rho
        real(dp), parameter :: zerodp = 0._dp
        real(dp), parameter :: epsdp = epsilon(0._dp)

        ns = solvent(1)%nspec ! number of solvent species
        kT = thermo%kbT
        dv = grid%dv
        mu = getinput%dp( 'imposed_chempot', defaultvalue=0._dp )
        if ( mu/=0._dp) stop "mu /=0 in module_energy_ideal_and_external. That's implemented but for now shutitoff"
        if ( ns/=1 .AND. mu/=0._dp) STOP "Imposing a chemical potential is valid only for single-species solvent"

        !
        ! fid = integrate over whole space of   x.log(x/x0)-x+x0 = Int[x.(log(x/x0)-1)] + Int[x0]
        !
        ! fext = Int[x*(vext-mu)]
        !
        fid = zerodp
        fext = zerodp
        do is=1,ns
          rho0 = solvent(is)%rho0 ! bulk density
          do iz=1,grid%nz
            do iy=1,grid%ny
              do ix=1,grid%nx
                do io=1,grid%no

                  xi = solvent(is)%xi(io,ix,iy,iz)
                  rho = xi**2 *rho0

                  if (rho > epsdp) then
                    fid = fid + kT*dv*grid%w(io)*( rho*log(rho/rho0)-rho+rho0 )
                    fext = fext + dv*grid%w(io)*rho*solvent(is)%vext(io,ix,iy,iz)
                    df(io,ix,iy,iz,is) = df(io,ix,iy,iz,is) &
                      + grid%w(io)*2._dp*rho0*xi*(kT*log(rho/rho0)+solvent(is)%vext(io,ix,iy,iz))
                  else
                    ! limit of Fid when rho->0 is not 0, it is:
                    fid = fid + kT*dv*grid%w(io)*  rho0
                  end if

                end do
              end do
            end do
          end do
        end do

    end subroutine energy_ideal_and_external
end module module_energy_ideal_and_external

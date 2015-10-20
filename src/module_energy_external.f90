module module_energy_external
    implicit none
    private
    public : energy_external
contains

    ! Compute the external contribution to free energy functional
    pure subroutine energy_external (fext, dfext)
        use precision_kinds, only: dp
        use module_solvent, only: solvent
        use module_input, only: getinput
        use module_grid, only: grid
        implicit none
        real(dp), intent(out) :: fext, dfext(:,:,:,:,:)
        integer :: ix, iy, iz, io, is
        real(dp) :: mu ! chemical potential

        ! Impose a chemical potential
        mu = getinput%dp( 'imposed_chempot', defaultvalue=0._dp) ! 0 means that you don't impose any chemical potential (delta_chempot is built this way)
        if ( solvent(1)%nspec/=1 .AND. imposedChemPot/=0._dp) STOP "Imposing a chemical potential is valid only for single-species solvent"

        ! F_{ext}[\rho(\vec{r},\vec{\Omega})]=\int d \vec{r} d \vec{\Omega} V_{ext}(\vec{r},\vec{\Omega})\rho(\vec{r},\vec{\Omega})

        ! I know that dfext is constant during the whole program but I don't want to keep it in memory
        fext = 0._dp
        do is=1,solvent(1)%nspec
            do io=1,grid%no
                do iz=1,grid%nz
                    do iy=1,grid%ny
                        do ix=1,grid%nx
                            fext                  = grid%w(io) * (solvent(is)%vext(ix,iy,iz,io) -mu) * solvent(s)%density(ix,iy,iz,io)
                            dfext(ix,iy,iz,io,is) = grid%w(io) * (solvent(is)%vext(ix,iy,iz,io) -mu)
                        end do
                    end do
                end do
            end do
        end do
    end subroutine energy_external
end module module_energy_external

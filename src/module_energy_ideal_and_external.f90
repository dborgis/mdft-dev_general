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
        real(dp), intent(inout) :: df(:,:,:,:,:)
        integer :: is, io, ns, ix, iy, iz
        real(dp) :: x, x0, vext, volume, dv, kT, mu, w, dfid
        real(dp), parameter :: zerodp = 0._dp
        real(dp), parameter :: epsdp = epsilon(0._dp)

        if (.not.allocated(solvent)) call print_solvent_not_allocated ("Look at subroutine energy_ideal_and_external")

        ns = solvent(1)%nspec
        kT = thermo%kbT
        dv = grid%dv
        volume = grid%v ! volume
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
        do is=1,ns
            x0 = solvent(is)%rho0 ! bulk density
            do io=1,grid%no
                w = grid%w(io) ! weight of the orientation
                do iz=1,grid%nz
                    do iy=1,grid%ny
                        do ix=1,grid%nx
                            x = solvent(is)%density(ix,iy,iz,io)
                            vext = solvent(is)%vext(ix,iy,iz,io) - mu
                            ! if (x<1e-10 .and. vext<36._dp) then
                            !     print*, x, vext
                            ! end if
                            if (x>epsdp) then
                                fid = fid + kT*dv*w*(x*log(x/x0)-x+x0)
                                dfid = kT*dv*w*log(x/x0)
                            else if (x<=epsdp) then
                                fid = fid + kT*dv*w*x0 ! par continuite. Plot x*log(x/x0)-x+x0 dans Mathematica pour un x non nul mais arbitrairement petit pour t'en convaincre.
                                dfid = 0._dp
                            end if
                            fext = fext + dv*w*vext*x
                            df(ix,iy,iz,io,is) = df(ix,iy,iz,io,is) + dfid + dv*w*vext



                            ! if (vext>=solvent(is)%vext_threeshold ) then
                            !     fid = fid + kT*(x*log(x/x0)-x+x0)*dv*w
                            !     fext = fext + x*vext*dv*w
                            !     df (ix,iy,iz,io,is) = df (ix,iy,iz,io,is) + w*dv*(kT*log(x/x0)+vext)
                            ! else ! we take the limit, which is well-behaved
                            !     fid = fid + kT*(x*log(x/x0)-x+x0)*dv*w
                            !     fext = fext + x*vext*dv*w
                            !     df (ix,iy,iz,io,is) = df (ix,iy,iz,io,is) + w*(kT*log(x/x0)+vext)
                            ! end if

                            ! fid = fid + kT*(x*log(x/x0)-x+x0)*dv*w
                            ! fext = fext + x*vext*dv*w
                            ! df (ix,iy,iz,io,is) = df (ix,iy,iz,io,is) + w*(kT*log(x/x0)+vext)

                        end do
                    end do
                end do
            end do
        end do

    end subroutine energy_ideal_and_external
end module module_energy_ideal_and_external

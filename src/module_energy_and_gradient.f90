module module_energy_and_gradient
    use precision_kinds, only: dp
    implicit none
    private
    type f_type
        real(dp) :: id = 0._dp,&
        ext = 0._dp,&
        exc_cs = 0._dp,&
        exc_cdeltacd = 0._dp,&
        exc_cproj = 0._dp,&
        exc_fmt = 0._dp,&
        exc_wca = 0._dp,&
        exc_3b = 0._dp,&
        exc_dipolar = 0._dp,&
        exc_multipolar_without_coupling_to_density = 0._dp,&
        exc_multipolar_with_coupling_to_density = 0._dp,&
        exc_hydro = 0._dp,&
        exc_nn_cs_plus_nbar = 0._dp
    end type
    type (f_type), protected :: ff
    public :: energy_and_gradient
contains

    subroutine energy_and_gradient (f, df)

        ! In this SUBROUTINE one calls the different parts of the total energy
        ! first is computed the radial_part, then ... blabla.
        ! this SUBROUTINE is the one called by the minimization stuff
        ! for computing the total energy and associated gradient
        ! FF is the TOTAL ENERGY of the system, it is thus the functional of the density that is minimized by solver
        ! dF_new is the gradient of FF with respect to all coordinates. Remember it is of the kind dF_new ( number of variables over density (ie angles etc))

        use precision_kinds, only: dp
        use module_solvent, only: solvent
        use module_grid, only: grid
        use module_energy_ideal_and_external, only: energy_ideal_and_external
        use module_energy_cs, only: energy_cs
        use module_energy_cdeltacd, only: energy_cdeltacd
        use module_energy_cproj, only: energy_cproj

        implicit none

        real(dp), intent(out) :: f
        real(dp), intent(out) :: df (grid%nx, grid%ny, grid%nz, grid%no, solvent(1)%nspec)
        real(dp) :: df_loc (grid%nx, grid%ny, grid%nz, grid%no, solvent(1)%nspec)
        real(dp), parameter :: zerodp=0._dp
        real :: t(10)
        integer :: s, ns

        if (.not. allocated(solvent)) then
            print*, "in energy_and_gradient, solvent()% is not allocated"
            error stop
        end if
        ns = solvent(1)%nspec

        f  = zerodp
        df = zerodp
        df_loc = zerodp

        print*,

        do s=1,solvent(1)%nspec
            if (solvent(s)%do%id_and_ext) then
                call cpu_time(t(1))
                call energy_ideal_and_external (ff%id, ff%ext, df_loc)
                call cpu_time(t(2))
                print*, "ff%ext          =", ff%ext
                print*, "ff%id           =", ff%id, " in",t(2)-t(1),"sec"
                f=f+ff%id+ff%ext
                df = df + df_loc
            end if
            if (solvent(s)%do%exc_cs) then
                call cpu_time(t(3))
                call energy_cs (ff%exc_cs, df_loc)
                call cpu_time(t(4))
                print*, "ff%exc_cs       =", ff%exc_cs,       "and norm2@df_exc_cs      =",norm2(df_loc),"in",t(4)-t(3),"sec"
                f = f + ff%exc_cs
                df = df + df_loc
            end if
            if (solvent(s)%do%exc_cdeltacd) then
                call cpu_time(t(5))
                call energy_cdeltacd (ff%exc_cdeltacd, df_loc)
                call cpu_time(t(6))
                print*, "ff%exc_cdeltacd =", ff%exc_cdeltacd, "and norm2@df_exc_cdeltacd=",norm2(df_loc),"in",t(6)-t(5),"sec"
                f = f + ff%exc_cdeltacd
                df = df + df_loc
            end if
            if (solvent(s)%do%exc_cproj) then
                call cpu_time(t(7))
                call energy_cproj (ff%exc_cproj, df_loc)
                call cpu_time(t(8))
                print*, "ff%exc_cproj    =", ff%exc_cproj,    "and norm2@df_exc_cproj   =",norm2(df_loc),"in",t(8)-t(7),"sec"
                ! f = f + ff%exc_cproj
                ! df = df + df_loc
            end if
            ! if (solvent(s)%do%exc_fmt) call energy_fmt (ff%exc_fmt, df)
            ! if (solvent(s)%do%wca) call lennard_jones_perturbation_to_hard_spheres (ff%exc_wca, df)
            ! if (solvent(s)%do%exc_multipolar_without_coupling_to_density) &
            !         call energy_polarization_multi (ff%exc_multipolar_without_coupling_to_density, df)
            ! if (solvent(s)%do%exc_multipolar_with_coupling_to_density) &
            !         call energy_polarization_multi_with_nccoupling (ff%exc_multipolar_with_coupling_to_density, df)
            ! if (solvent(s)%do%exc_hydro) call energy_hydro (ff%exc_hydro, df)
            ! if (solvent(s)%do%exc_nn_cs_plus_nbar) call energy_nn_cs_plus_nbar (ff%exc_nn_cs_plus_nbar, df)
            ! if (solvent(s)%do%exc_3b) call energy_threebody_faster (ff%exc_3d, df)
        end do

        print*, "____________"
        print*, "TOTAL (FF), |df| =", f, norm2(df)


    end subroutine energy_and_gradient
end module module_energy_and_gradient

module module_energy_and_gradient
  use iso_c_binding, only: c_float
  use precision_kinds, only: dp
  implicit none
  private
  type f_type
    real(dp) :: id = 0._dp,&
                ext = 0._dp,&
                exc_cs = 0._dp,&
                exc_cdeltacd = 0._dp,&
                exc_cproj = 0._dp,&
                exc_ck_angular = 0._dp,&
                exc_fmt = 0._dp,&
                exc_wca = 0._dp,&
                exc_3b = 0._dp,&
                exc_dipolar = 0._dp,&
                exc_multipolar_without_coupling_to_density = 0._dp,&
                exc_multipolar_with_coupling_to_density = 0._dp,&
                exc_hydro = 0._dp,&
                exc_nn_cs_plus_nbar = 0._dp, &
                tot=0._dp,&
                pscheme_correction=-999._dp,&
                pbc_correction=-999._dp
    integer :: ieval=0
  end type
  type (f_type), public :: ff
  public :: energy_and_gradient




contains




subroutine energy_and_gradient (f, df)

    ! In this subroutine one calls the different parts of the total energy
    ! This subroutine is the one called by the minimization stuff
    ! for computing the total energy and associated gradient
    ! FF is the TOTAL ENERGY of the system, it is thus the functional of the density that is minimized by solver
    ! dF_new is the gradient of FF with respect to all coordinates. Remember it is of the kind dF_new ( number of variables over density (ie angles etc))

    use iso_c_binding, only: c_float
    use precision_kinds, only: dp
    use module_solvent, only: solvent
    use module_grid, only: grid
    use module_energy_ideal_and_external, only: energy_ideal_and_external
    ! use module_energy_cs, only: energy_cs
    ! use module_energy_cdeltacd, only: energy_cdeltacd
    use module_energy_cproj_mrso, only: energy_cproj_mrso
    ! use module_energy_cproj_no_symetry, only: energy_cproj_no_symetry
    ! use module_energy_ck_angular, only: energy_ck_angular
    ! use module_energy_luc, only: energy_luc
    ! use module_energy_luc_fast, only: energy_luc_fast
    use module_input, only: getinput
    ! use module_energy_cproj_slow, only: energy_cproj_slow
    use module_lbfgs_nocedal_mdft, only:lbfgsb

    implicit none

    real(dp), intent(out) :: f
    real(dp), intent(out), optional :: df (grid%no, grid%nx, grid%ny, grid%nz, solvent(1)%nspec)
    real(dp), parameter :: zerodp=0._dp
    real(c_float) :: t(10)
    real(dp) :: fold
    integer :: ns, s

    ff%ieval = ff%ieval +1

    if (.not. allocated(solvent)) then
      print*, "in energy_and_gradient, solvent()% is not allocated"
      error stop
    end if
    ns = solvent(1)%nspec

    fold=f
    f  = zerodp
    if(present(df)) df = zerodp
    s=1


    !
    ! Ideal and external free energy functionals
    !
    if (solvent(s)%do%id_and_ext) then
      if(present(df)) then
        call cpu_time(t(1))
        call energy_ideal_and_external (ff%id, ff%ext, df)
        call cpu_time(t(2))
        ! print*, "ff%ext            =", real(ff%ext)
        ! print*, "ff%id             =", real(ff%id), " in",t(2)-t(1),"sec"
      else
        call energy_ideal_and_external (ff%id, ff%ext)
      end if
      f = f +ff%id +ff%ext
    end if

    !
    ! purely radial component, with c_s
    !
    ! if (solvent(s)%do%exc_cs) then
    !     call cpu_time(t(3))
    !     call energy_cs (ff%exc_cs, df)
    !     call cpu_time(t(4))
    !     print*, "ff%exc_cs         =", real(ff%exc_cs), " in",t(4)-t(3),"sec"
    !     f = f + ff%exc_cs
    ! end if

    !
    ! with c_d and c_delta
    !
    ! if (solvent(s)%do%exc_cdeltacd) then
    !     call cpu_time(t(5))
    !     call energy_cdeltacd (ff%exc_cdeltacd, df)
    !     call cpu_time(t(6))
    !     print*, "ff%exc_cdeltacd   =", real(ff%exc_cdeltacd), " in",t(6)-t(5),"sec"
    !     f = f + ff%exc_cdeltacd
    ! end if


    !
    ! with Luc's routine
    !
    if (solvent(s)%do%exc_cproj) then
        if(present(df)) then
          call cpu_time(t(5))
          call energy_cproj_mrso( ff%exc_cproj, df, print_timers=.false.)
          call cpu_time(t(6))
        !   print*, "ff%exc_cproj_mrso =", real(ff%exc_cproj), " in",t(6)-t(5),"sec"
        else
          call energy_cproj_mrso( ff%exc_cproj, print_timers=.false.)
        end if
        f = f + ff%exc_cproj
    end if



    ! if (solvent(s)%do%exc_cproj) then
    !     call cpu_time(t(7))
    !     ! call energy_cproj_mrso (ff%exc_cproj, df)
    !     ! call energy_cproj_no_symetry (ff%exc_cproj, df)
    !     ! call energy_luc ( ff%exc_cproj , df )
    !     call energy_luc_fast ( ff%exc_cproj, df)
    !     call cpu_time(t(8))
    !     print*, "ff%exc_cproj      =", ff%exc_cproj,   "in",t(8)-t(7),"sec"
    !     ! stop "energy_and_gradient after call to energy_luc"
    !     ! print*, "ff%exc_cproj - (ff%exc_cs+ff%exc_cdeltacd) =", ff%exc_cproj-(ff%exc_cs + ff%exc_cdeltacd)
    !     f = f + ff%exc_cproj
    !     ! stop
    ! end if

    ! if (solvent(s)%do%exc_ck_angular) then
    !     call cpu_time(t(9))
    !     call energy_ck_angular (ff%exc_ck_angular, df)
    !     call cpu_time(t(10))
    !     print*, "ff%exc_ck_angular =", ff%exc_ck_angular,"in",t(10)-t(9),"sec"
    !     f = f + ff%exc_ck_angular
    ! end if

    !
    ! adhoc corrections to the solvation free energy (Hunenberger, pressure etc.)
    !
    if (.not. getinput%log('direct_sum', defaultvalue=.false.)) then
        call typeB_corrections
        f = f + ff%pbc_correction
        call typeC_corrections
        f = f + ff%pscheme_correction
    else
        ff%pbc_correction=0._dp
        ff%pscheme_correction=0._dp
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

    ff%tot = f



    block
        logical, save :: printheader = .true.
        if(printheader) then
            write(*,'(A5,11A14)') "#eval","Ftot","Fext","Fid","Fexc","Cpbc","Cpsch","relF","pgtol","Ttot","Text+id","Texc"
            printheader = .false.
        end if
    end block
    block
        real(dp) :: reldf, Texc, Ttot, Textid, pgtol
        Texc = t(6)-t(5)
        Textid = t(2)-t(1)
        Ttot = Texc+Textid
        pgtol = real(maxval(df))
        reldf = (fold-f)/maxval([abs(fold),abs(f),1._dp])
        write(*,"(I5,11F14.4)") ff%ieval, ff%tot, ff%ext, ff%id, ff%exc_cproj, ff%pbc_correction, ff%pscheme_correction, reldf, pgtol, Ttot, Textid, Texc
    end block

! if(present(df)) then
!     print*, "Δf/f =", real((fold-f)/maxval([abs(fold),abs(f),1._dp])) ,"target=",lbfgsb%factr*epsilon(1._dp)
!     print*, "pgtol=", real(maxval(df)),                                "target=",lbfgsb%pgtol
! end if

end subroutine energy_and_gradient





  subroutine typeB_corrections
    !
    !   Type B correction of Hunenberger : finite size and periodicity (\propto q_{I}^{2}/L)
    !   see J. Chem. Phys. 124, 224501 (2006), eq. 32 with R_i=0 (ionic radius = 0)
    !   see also Hunenberger and McCammon, JCP 110, 1856 (1999), doi: 10.1063/1.477873
    !
    use precision_kinds, only: dp
    use module_solute, only: solute
    use module_grid, only: grid
    implicit none
    real(dp) :: solute_net_charge, L
    real(dp), parameter :: dielectric_constant_spce=71._dp
    solute_net_charge = sum (solute%site%q)
    if (.not. all(grid%length==grid%length(1))) then
      print*, "The grid is not cubic."
      print*, "The periodic boundary conditions correction is intended for cubic cells."
      print*, "We use the average length sum(len)/3."
      L = sum(grid%length)/3._dp
    else
      L = grid%lx
    end if
    if (L<=epsilon(1._dp)) then
      error stop "sherY6S%hx6YYUJ"
    end if
    ! The dielectric constant of SPC/E water is 71. See Kusalik and Svishchev, "The Spatial Structure in Liquid Water", Science 265, 1219 (1994) doi:10.1126/science.265.5176.1219
    ! The original SPC/E paper by Berendsen does not provide this information.
    ! see https://www.wolframalpha.com/input/?i=-2.837297*(electron+charge)%5E2%2F(4*pi*vacuum+permittivity*2*angstroms)+to+kJ%2Fmol
    ff%pbc_correction = -1971.01_dp*solute_net_charge**2/L*(1-1._dp/dielectric_constant_spce)
  end subroutine typeB_corrections




  subroutine typeC_corrections
    !
    ! Hunenberger's corrections
    !
    !... because we use PBC/LS with P-scheme electrostatics in MDFT.
    ! See Kastenholz and Hunenberger, JCP 124, 124106 (2006), page 224501-8, equations 35, 36 and 37 with the radius of the ion, R_i = 0
    ! "To be applied if the solvent molecule is rigid and involves a single van der Waals interaction site M,
    ! and that any scheme relying on molecular-cutoff truncation refers to this specific site for applying the truncation."
    ! R_i should not be 0 in reality (see Hunenberger's paper) but (i) the contribution is small for small ions.
    ! For bigger ions (ie charged molecules), what should one do ?
    !
      use precision_kinds, only: dp
      use module_solvent, only: solvent
      use module_solute, only: solute
      implicit none
      real(dp) :: solute_net_charge ! net charge of the solute
      real(dp) :: gamma ! trace of the quadrupole moment. Should be 0.848 e.nm² for SPCE and 0.820 for SPC water.
      gamma = solvent(1)%quadrupole(1,1)+solvent(1)%quadrupole(2,2)+solvent(1)%quadrupole(3,3) ! quadrupole moment trace
      solute_net_charge = sum(solute%site%q)
      ff%pscheme_correction = -gamma*solvent(1)%n0*2909.857_dp*solute_net_charge ! in kJ/mol
  end subroutine typeC_corrections


end module module_energy_and_gradient

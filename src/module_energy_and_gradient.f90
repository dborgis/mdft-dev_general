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

    use iso_c_binding, only: c_sizeof
    use precision_kinds, only: dp
    use module_solvent, only: solvent
    use module_grid, only: grid
    use module_energy_ideal_and_external, only: energy_ideal_and_external
    use module_energy_cs, only: energy_cs
    use module_energy_cdeltacd, only: energy_cdeltacd
    use module_energy_cproj_mrso, only: energy_cproj_mrso
    use module_energy_cproj_no_symetry, only: energy_cproj_no_symetry
    use module_energy_ck_angular, only: energy_ck_angular
    use module_energy_luc, only: energy_luc
    use module_energy_luc_fast, only: energy_luc_fast
    use module_input, only: getinput
    ! use module_energy_cproj_slow, only: energy_cproj_slow

    implicit none

    real(dp), intent(out) :: f
    real(dp), intent(out) :: df (grid%no, grid%nx, grid%ny, grid%nz, solvent(1)%nspec)
    ! real(dp) :: df_trash (grid%no, grid%nx, grid%ny, grid%nz, solvent(1)%nspec)

    real(dp), parameter :: zerodp=0._dp
    real :: t(10)
    real(dp) :: fold
    integer :: s, ns

    if (.not. allocated(solvent)) then
      print*, "in energy_and_gradient, solvent()% is not allocated"
      error stop
    end if
    ns = solvent(1)%nspec

    if (allocated (solvent(1)%vextq)) then
      do s=1,ns
        deallocate (solvent(s)%vextq)
      end do
    end if

    fold=f
    f  = zerodp
    df = zerodp
    ! df_trash = zerodp

    print*,


    do s=1,solvent(1)%nspec

      !
      ! Ideal and external free energy functionals
      !
      if (solvent(s)%do%id_and_ext) then
        call cpu_time(t(1))
        call energy_ideal_and_external (ff%id, ff%ext, df)
        call cpu_time(t(2))
        print*, "ff%ext            =", ff%ext
        print*, "ff%id             =", ff%id, " in",t(2)-t(1),"sec"
        f = f +ff%id +ff%ext
      end if

      !
      ! purely radial component, with c_s
      !
      if (solvent(s)%do%exc_cs) then
        call cpu_time(t(3))
        call energy_cs (ff%exc_cs, df)
        call cpu_time(t(4))
        print*, "ff%exc_cs         =", ff%exc_cs, " in",t(4)-t(3),"sec"
        f = f + ff%exc_cs
      end if

      !
      ! with c_d and c_delta
      !
      if (solvent(s)%do%exc_cdeltacd) then
        call cpu_time(t(5))
        call energy_cdeltacd (ff%exc_cdeltacd, df)
        call cpu_time(t(6))
        print*, "ff%exc_cdeltacd   =", ff%exc_cdeltacd, "in",t(6)-t(5),"sec"
        f = f + ff%exc_cdeltacd
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
        if( ff%pscheme_correction==-999._dp) call corrections
        print*, "ff%pscheme corr   =", ff%pscheme_correction
        f = f + ff%pscheme_correction
        print*, "ff%pbc correction =", ff%pbc_correction
        f = f + ff%pbc_correction
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

    ff%tot = f

    print*, "-------------------------------------------------------------------------------"
    print*, "TOTAL (FF) =", real(f), "|   Δf/f =", real((fold-ff%tot)/&
             maxval([abs(fold),abs(f),1._dp])), "|  l2@df=",real(norm2(df))
    print*, "-------------------------------------------------------------------------------"
    print*,

  end subroutine energy_and_gradient









  subroutine corrections
    implicit none
    !
    ! Hunenberger's corrections
    !
    !... We use P-scheme instead of M-scheme for the electrostatics in MDFT.
    ! See Kastenholz and Hunenberger, JCP 124, 124106 (2006), page 224501-8, equations 35, 36 and 37 with the radius of the ion, R_i = 0
    ! "To be applied if the solvent molecule is rigid and involves a single van der Waals interaction site M,
    ! and that any scheme relying on molecular-cutoff truncation refers to this specific site for applying the truncation."
    ! R_i should not be 0 in reality (see Hunenberger's paper) but (i) the contribution is small for small ions.
    ! For bigger ions (ie charged molecules), what should one do ?
    !
    block
      use module_solvent, only: solvent
      use module_solute, only: solute
      double precision :: solute_net_charge ! net charge of the solute
      double precision :: gamma ! trace of the quadrupole moment
      gamma = solvent(1)%quadrupole(1,1)+solvent(1)%quadrupole(2,2)+solvent(1)%quadrupole(3,3) ! quadrupole moment trace
      solute_net_charge = sum(solute%site%q)
      ff%pscheme_correction = -gamma*solvent(1)%n0*2.909857E3*solute_net_charge ! in kJ/mol
      open(79,file="output/Pscheme_correction")
      write(79,*) ff%pscheme_correction
      close(79)
    end block

    !
    !   Type B correction of Hunenberger (due to lattice sums with periodic boundary conditions)
    !   see J. Chem. Phys. 124, 224501 (2006), eq. 32 with R_i=0 (ionic radius = 0)
    !
    block
      use module_solute, only: solute
      use module_grid, only: grid
      double precision :: solute_net_charge, L
      solute_net_charge = sum (solute%site%q)
      if (.not. all(grid%length==grid%length(1))) then
        print*, "The grid is not cubic."
        print*, "The periodic boundary conditions correction is intended for cubic cells."
        print*, "We use the average length sum(len)/3."
        L = sum(grid%length)/3._dp
      else
        L = grid%length(1)
      end if
      if (L<=epsilon(1._dp)) then
        error stop "sherY6S%hx6YYUJ"
      end if
      ff%pbc_correction = -1949.0466_dp*solute_net_charge**2/L
      open(79,file="output/PBC_correction")
      write(79,*) ff%pbc_correction
      close(79)
    end block

  end subroutine corrections


end module module_energy_and_gradient

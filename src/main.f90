! Molecular (classical) density functional theory code.
! main program : skeleton of the structure.

program main

    use iso_c_binding, only: dp => c_double
!    use grid_mod, only: grid
!    use density_mod, only: density
!    use ideal_functional_mod, only: fid_get => get
!    use functional_mod, only: functional
!    use external_potential_mod, only: external_potential

    implicit none
    integer(8) :: count0, count1, count_rate
    real :: mdft_wholetime


!    type(grid) :: gr
!    type(density) :: d
!    type(functional) :: ff, fid
!    type(external_potential) :: vext
!
!    call gr%build()
!    call d%build(gr)
!    call ff%build(gr)
!    call fid%build(gr)
!    call random_number(d%rho )
!    call fid_get(gr,d,fid)
!    ff%energy = ff%energy + fid%energy
!    ff%grad   = ff%grad   + fid%grad
!    ! call fext_get(gr,d,fext)
!
!    ! print*, ideal_functional_get(gr,d)
!stop "ok for the new main"
!    ! print*,fid%energy,"kj/mol"
!
!
!

    call system_clock (count0, count_rate)

    CALL init_simu
    CALL find_equilibrium_density
    CALL process_output

    write(*,'(A)')"=="

    call system_clock (count1)
    mdft_wholetime = (count1-count0)/real(count_rate)

    if( mdft_wholetime < 5*60 ) then ! less than 5 minutes
        write(*,'(A,F12.2,A)')"MDFT finished with status OK. CPU time",mdft_wholetime," sec."
    else if( mdft_wholetime < 5*60*60 ) then ! less than 5 h
        write(*,'(A,F12.2,A)')"MDFT finished with status OK. CPU time",mdft_wholetime/60.," min."
    else
        write(*,'(A,F12.2,A)')"MDFT finished with status OK. CPU time",mdft_wholetime/60./60.," hours."
    end if

end program

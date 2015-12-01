! Molecular (classical) density functional theory code.
! main program : skeleton of the structure.

program main

    use precision_kinds, only: dp
    use module_postprocessing, only: init_postprocessing
    use module_init_simu, only: init_simu

    implicit none
    integer(8) :: count0, count1, count_rate
    real :: mdft_wholetime

    call system_clock (count0, count_rate)




    call init_simu
    call energy_minimization
    call init_postprocessing

    print*,
    print*,
    print*,"---"

    call system_clock (count1)
    mdft_wholetime = real(count1-count0)/real(count_rate)
    if( mdft_wholetime < 5*60 ) then ! less than 5 minutes
        write(*,'(A,F12.2,A)') "MDFT finished in",mdft_wholetime," sec."
    else if( mdft_wholetime < 5*60*60 ) then ! less than 5 h
        write(*,'(A,F12.2,A)') "MDFT finished in",mdft_wholetime/60.," min."
    else
        write(*,'(A,F12.2,A)') "MDFT finished in",mdft_wholetime/60./60.," hours."
    end if

end program

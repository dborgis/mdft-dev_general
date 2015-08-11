! Molecular (classical) density functional theory code.
! main program : skeleton of the structure.

program main

    implicit none
    integer(8) :: count0, count1, count_rate
    real :: mdft_wholetime

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

! Molecular (classical) density functional theory code.
! main program : skeleton of the structure.

program main

    use precision_kinds, only: dp
    use module_postprocessing, only: init_postprocessing
    use module_init_simu, only: init_simu

    implicit none
    integer(8) :: count0, count1, count_rate
    real :: mdft_execution_time

    call system_clock (count0, count_rate)
    if (count_rate == 0) error stop "Bug in main.f90, count_rate==0"





    call init_simu
    call energy_minimization
    call init_postprocessing





    print*
    print*
    print*,"---"

    call system_clock (count1)

    block
        real :: five_min = 5*60
        real :: five_hours = 5*60*60
        real :: mdft_execution_time
        mdft_execution_time = real(count1-count0)/real(count_rate)
        if(      mdft_execution_time < five_min ) then
            write(*,'(A,F12.2,A)') "MDFT finished in", mdft_execution_time," sec."
        else if( mdft_execution_time < five_hours ) then
            write(*,'(A,F12.2,A)') "MDFT finished in",mdft_execution_time/60.," min."
        else
            write(*,'(A,F12.2,A)') "MDFT finished in",mdft_execution_time/60./60.," hours."
        end if
    end block

end program

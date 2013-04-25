! Molecular (classical) density functional theory code.
! main program : skeleton of the structure.
program mdft

    use precision_kinds , only : dp ! definition of the precision kinds
    implicit none
    real(dp) :: time0, time1 ! time steps

    call cpu_time ( time0 ) ! init timer
    call print_header ! package name, day & time
    call init_simu ! initialization of simulation. read, allocate etc.
    call find_equilibrium_density ! DFT part! Minimize functional of the density(position, orientation)
    call process_output ! process results
    call close_simu
    call cpu_time ( time1 ) ! close timer and tell user runtime
    print*,'Total execution time =' , time1 - time0

end program mdft

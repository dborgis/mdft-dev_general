! Molecular (classical) density functional theory code.

! main program : skeleton of the structure.

program mdft

  use precision_kinds , only : dp ! definition of the precision kinds

  implicit none
  real ( kind = dp ) :: time0 , time1 ! time steps


  ! print header
  call print_header

  ! init timer
  call cpu_time ( time0 )

  ! init simu
  call init_simu

  ! Minimize functional of the density(position, orientation)
  call find_equilibrium_density

  ! process results
  call process_output

  ! close simulation properly (deallocate etc) (may be a bit useless as gfortran closes everything as properly as possible by default)
  call close_simu

  ! close timer and inform user
  call cpu_time ( time1 )
  print*,'Total execution time =' , time1 - time0

end program mdft

! Molecular (classical) density functional theory code.
! main program : skeleton of the structure.
PROGRAM mdft

    USE precision_kinds ,ONLY: dp

    IMPLICIT NONE
    REAL(dp) :: time0, time1, cput

    CALL CPU_TIME (time0)

    CALL init_simu
    CALL find_equilibrium_density
    CALL process_output

    write(*,'(A)')"=="

    CALL CPU_TIME (time1)
    cput = time1-time0
    if( cput < 5*60 ) then ! less than 5 minutes
      write(*,'(A,F12.2,A)')"MDFT finished with status OK. CPU time",cput," sec."
    else if( cput < 5*60*60 ) then ! less than 5 h
      write(*,'(A,F12.2,A)')"MDFT finished with status OK. CPU time",cput/60.," min."
    else
      write(*,'(A,F12.2,A)')"MDFT finished with status OK. CPU time",cput/60./60.," hours."
    end if

END PROGRAM mdft

! Molecular (classical) density functional theory code.
! main program : skeleton of the structure.
PROGRAM mdft

    USE precision_kinds ,ONLY: dp
    
    IMPLICIT NONE
    REAL(dp) :: time0, time1

    CALL CPU_TIME (time0)

    CALL init_simu
    CALL find_equilibrium_density
    CALL process_output

    CALL CPU_TIME (time1)
    PRINT*,'Execution time =',time1-time0

END PROGRAM mdft

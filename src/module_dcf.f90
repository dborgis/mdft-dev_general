MODULE dcf

    USE precision_kinds, ONLY: dp, i2b
    USE input, ONLY: input_log, n_linesInFile

    IMPLICIT NONE


    CONTAINS
    
    
    SUBROUTINE init
        IF (input_log('read_ck_or_chi')) CALL read_ck! If calculation based on direct correlation functions read them in kspace in cs.in, cd.in, cdelta.in. depends on tag read_ck = T    or read_ck = F  in input/dft.in
        IF (input_log('bridge_hard_sphere')) CALL cs_of_k_hard_sphere! in case of bridge calculation, one also need the direct correlation function c2 of the hard sphere.
        IF (input_log('hard_sphere_fluid')) CALL compute_hard_spheres_parameters!> If calculation based on Fundamental Measure Theory read FMT parameters and compute weight functions etc
    END SUBROUTINE init
    
    
    
!~     SUBROUTINE readDensityDensityCorrelationFunction (cs)
!~     
!~     END SUBROUTINE readDensityDensityCorrelationFunction
!~ 
!~ 
!~     SUBROUTINE readPolarizationPolarizationCorrelationFunction (c_delta,c_d)
!~     
!~     END SUBROUTINE readPolarizationPolarizationCorrelationFunction
!~     
!~     
!~     SUBROUTINE readDielectricSusceptibilities (chi_l,chi_t)
!~     
!~     END SUBROUTINE readDielectricSusceptibilities




END MODULE dcf

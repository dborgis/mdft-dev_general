! Subroutine to initiate simulation parameters
SUBROUTINE init_simu

    USE quadrature  ,ONLY: prepare_quadrature => init
    USE fft         ,ONLY: prepare_fft => init
    USE dcf         ,ONLY: init_dcf => init
    USE input       ,ONLY: input_log, input_char
    
    IMPLICIT NONE

    CALL print_header ! package name, day & time
    CALL put_input_in_character_array ! more powerfull way of reading dft.in : put each line in an array
    CALL allocate_from_input ! TODO this should be removed and done only when needed ! read dft.in, solute.in and solvent.in
    CALL prepare_fft
    CALL print_input_to_output_folder! Print input parameters found in input files to output folder
    CALL read_solvent ! Read solvent atomic positions, charge and Lennard-Jones param
    CALL read_solute ! Read solute atomic positions, charge and Lennard-Jones param
    CALL prepare_quadrature ! prepare numerical integration (for angles)
    CALL init_dcf
    IF (input_log('hard_sphere_fluid')) CALL compute_hard_spheres_parameters ! If calculation based on Fundamental Measure Theory read FMT parameters and compute weight functions etc
    CALL prepare_minimizer ! Prepare minimization ! allocate tables and computes precision and so on
    CALL init_external_potential
    CALL init_solvent_polarization
    CALL init_density

END SUBROUTINE init_simu

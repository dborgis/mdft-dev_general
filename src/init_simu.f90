! Subroutine to initiate simulation parameters
SUBROUTINE init_simu

    USE input, ONLY: input_log
    USE quadrature, ONLY: prepare_quadrature => init
    USE fft, ONLY: prepare_fft => init
    IMPLICIT NONE

    CALL print_header ! package name, day & time
    CALL put_input_in_character_array! more powerfull way of reading dft.in : put each line in an array
    CALL allocate_from_input ! TODO this should be removed and done only when needed ! read dft.in, solute.in and solvent.in
    CALL prepare_fft
    CALL read_solvent! Read solvent atomic positions, charge and Lennard-Jones param
    CALL read_solute! Read solute atomic positions, charge and Lennard-Jones param
    CALL soluteChargeDensityFromSoluteChargeCoordinates
    CALL print_input_to_output_folder! Print input parameters found in input files to output folder
    CALL prepare_quadrature ! prepare numerical integration (for angles)
    IF (input_log('read_ck_or_chi')) CALL read_ck! If calculation based on direct correlation functions read them in kspace in cs.in, cd.in, cdelta.in. depends on tag read_ck = T    or read_ck = F  in input/dft.in
    IF (input_log('bridge_hard_sphere')) CALL cs_of_k_hard_sphere! in case of bridge calculation, one also need the direct correlation function c2 of the hard sphere.
    IF (input_log('hard_sphere_fluid')) CALL compute_hard_spheres_parameters!> If calculation based on Fundamental Measure Theory read FMT parameters and compute weight functions etc
    CALL prepare_minimizer! Prepare minimization ! allocate tables and computes precision and so on
    CALL init_external_potential! Compute external potential ( electrostatic, lennard-jones, yukawa, hard wall, ... )
    CALL init_density! Compute initial density(position, orientation) -> cg_vect(:) ! Compute rho at first iteration put directly in CG_vect ! need cg_vect being defined (in prepare minimizer)! vcoul is now useless, deallocate it. only vext is used in the free energy functional

END SUBROUTINE init_simu

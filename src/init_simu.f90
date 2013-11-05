! Subroutine to initiate simulation parameters
subroutine init_simu

    use precision_kinds , only: i2b
    use input, only: input_log, input_char! contains all input file dft.in in a character array
    use quadrature, only : prepare_quadrature => init
    USE fft, ONLY: prepare_fft => init
    implicit none

    call put_input_in_character_array! more powerfull way of reading dft.in : put each line in an array
    call allocate_from_input ! TODO this should be removed and done only when needed ! read dft.in, solute.in and solvent.in
    call prepare_fft
    call read_solvent! Read solvent atomic positions, charge and Lennard-Jones param
    call read_solute! Read solute atomic positions, charge and Lennard-Jones param
    call charge_density_from_point_charge_positions! Get charge density rho_c from point charge positions
    call print_input_to_output_folder! Print input parameters found in input files to output folder
    call prepare_quadrature ! prepare numerical integration (for angles)
    if (input_log('read_ck_or_chi')) call read_ck! If calculation based on direct correlation functions read them in kspace in cs.in, cd.in, cdelta.in. depends on tag read_ck = T    or read_ck = F  in input/dft.in
    if (input_log('bridge_hard_sphere')) call cs_of_k_hard_sphere! in case of bridge calculation, one also need the direct correlation function c2 of the hard sphere.
    if (input_log('hard_sphere_fluid')) call compute_hard_spheres_parameters!> If calculation based on Fundamental Measure Theory read FMT parameters and compute weight functions etc
    call prepare_minimizer! Prepare minimization ! allocate tables and computes precision and so on
    call init_external_potential! Compute external potential ( electrostatic, lennard-jones, yukawa, hard wall, ... )
    call init_density! Compute initial density(position, orientation) -> cg_vect(:) ! Compute rho at first iteration put directly in CG_vect ! need cg_vect being defined (in prepare minimizer)! vcoul is now useless, deallocate it. only vext is used in the free energy functional

end subroutine init_simu

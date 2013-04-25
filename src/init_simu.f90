! Subroutine to initiate simulation parameters
subroutine init_simu
use precision_kinds , only: i2b
use input , only : input_line , input_log,input_char! contains all input file dft.in in a character array
implicit none
integer(i2b):: i , j ! dummy
! more powerfull way of reading dft.in : put each line in an array
call put_input_in_character_array
! read dft.in, solute.in and solvent.in
call allocate_from_input ! TODO this should be removed and done only when needed
! initialize FFTW3 plans
call init_fftw3_plans
! tabulate the value of norm(\vec{k}) and k2 = norm_k**2 in order to speed up program
call tabulate_normk_and_k2
! Read solvent atomic positions, charge and Lennard-Jones param
call read_solvent
! Read solute atomic positions, charge and Lennard-Jones param
call read_solute
! Get charge density rho_c from point charge positions
call charge_density_from_point_charge_positions
! Print input parameters found in input files to output folder
call print_input_to_output_folder
call get_psi_integration_roots_and_weights
!Read which quadrature is wanted and read the corresponding weights the angles
if (trim(adjustl(input_char('quadrature')))=='L') then
   call get_lebedev_integration_roots_and_weights  ! Lebedev integration for angular integrations 
   print*, 'quadrature = Lebedev'
end if
if (trim(adjustl(input_char('quadrature')))=='GL') then
   call get_gauss_legendre_integration_roots_and_weights  !Gauss Legendre integration for angular integrations 
   print*, 'quadrature = Gauss-Legendre'
end if
! If calculation based on direct correlation functions read them in kspace in cs.in, cd.in, cdelta.in
! test if cs cdelta cd have to be read.
! depends on tag read_ck = T    or read_ck = F  in input/dft.in
if (input_log('read_ck_or_chi')) call read_ck ! go read cs cdelta cd in input/
! in case of bridge calculation, one also need the direct correlation function c2 of the hard sphere.
if (input_log('bridge_hard_sphere')) call cs_of_k_hard_sphere
!> If calculation based on Fundamental Measure Theory read FMT parameters and compute weight functions etc
! look for tag 'hard_sphere_fluid'
if (input_log('hard_sphere_fluid')) call compute_hard_spheres_parameters
! Prepare minimization
call prepare_minimizer ! allocate tables and computes precision and so on
! Compute external potential ( electrostatic, lennard-jones, yukawa, hard wall, ... )
call init_external_potential
! Compute initial density(position, orientation) -> cg_vect(:)
call init_density ! Compute rho at first iteration put directly in CG_vect ! need cg_vect being defined (in prepare minimizer)
! vcoul is now useless, deallocate it. only vext is used in the free energy functional
end subroutine init_simu

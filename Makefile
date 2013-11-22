# modules directory
MODDIR = mod

# sources directory
SRCDIR = src

# The compiler
FC = gfortran
# flags for debugging or for maximum performance, comment as necessary
# flags forall (e.g. look for system .mod files, required in gfortran)
FCFLAGS = -J $(MODDIR)
# libraries needed for linking, unused in the examples
LDFLAGS = -lfftw3


#FLAGS = -g -fopenmp -lfftw3 -lfftw3f -lfftw3l -Jmodules
#-g Turns on code generating options for debugging
#-fopenmp enables openmp parallel compiling

#FLAGS = -lfftw3 -J $(MODDIR)

DEBUGFLAGS = -g -Wall -fbounds-check -fcheck=all -fcheck-array-temporaries -Warray-temporaries -Wconversion -pg -Wunused-parameter
#-g turns on debugging
#-p turns on profiling

OPTIM = -O3 -march=native
#-Ofast -funroll-loops -march=native
# FOR BACKUP : -march=native -O3 -ffast-math -funroll-loops   VERY AGRESSIVE
# -fopenmp for OPENMP support

EXE = mdft


OBJS = $(SRCDIR)/module_precision_kinds.f90 \
       $(SRCDIR)/module_constants.f90 \
       $(SRCDIR)/module_input.f90 \
       $(SRCDIR)/module_system.f90 \
       $(SRCDIR)/module_quadrature.f90 \
       $(SRCDIR)/module_fft.f90 \
       $(SRCDIR)/module_external_potential.f90 \
       $(SRCDIR)/module_minimizer.f90 \
       $(SRCDIR)/module_periodic_table.f90 \
       $(SRCDIR)/module_mod_lj.f90 \
       $(SRCDIR)/module_solute_geometry.f90 \
       $(SRCDIR)/module_bfgs.f90 \
       $(SRCDIR)/main.f90 \
       $(SRCDIR)/dblas1.f90 \
       $(SRCDIR)/trilinear_interpolation_v.f90 \
       $(SRCDIR)/get_charge_factor.f90 \
       $(SRCDIR)/get_charge_density.f90\
       $(SRCDIR)/get_charge_density_k.f90 \
       $(SRCDIR)/find_equilibrium_density.f90 \
       $(SRCDIR)/close_simu.f90 \
       $(SRCDIR)/print_supercell_xsf.f90 \
       $(SRCDIR)/process_output.f90 \
       $(SRCDIR)/get_final_density.f90 \
       $(SRCDIR)/get_final_polarization.f90 \
       $(SRCDIR)/compute_z_density.f90 \
       $(SRCDIR)/steepest_descent_main.f90 \
       $(SRCDIR)/mean_over_orientations.f90 \
       $(SRCDIR)/init_simu.f90 \
       $(SRCDIR)/soluteChargeDensityFromSoluteChargeCoordinates.f90 \
       $(SRCDIR)/convert_coordinate_into_icg.f90\
       $(SRCDIR)/poissonSolver.f90 \
       $(SRCDIR)/vext_q_from_v_c.f90 \
       $(SRCDIR)/read_hard_sphere_radius_and_allocate_if_necessary.f90 \
       $(SRCDIR)/cs_of_k_hard_sphere.f90 \
       $(SRCDIR)/init_density.f90 \
       $(SRCDIR)/init_fftw3_plans.f90 \
       $(SRCDIR)/init_external_potential.f90 \
       $(SRCDIR)/put_input_in_character_array.f90 \
       $(SRCDIR)/print_input_in_output_folder.f90 \
       $(SRCDIR)/get_gauss_legendre_integration_roots_and_weights.f90 \
       $(SRCDIR)/get_lebedev_integration_roots_and_weights.f90 \
       $(SRCDIR)/read_ck.f90 \
       $(SRCDIR)/read_solvent.f90 \
       $(SRCDIR)/read_solute.f90 \
       $(SRCDIR)/allocate_from_input.f90 \
       $(SRCDIR)/compute_hard_spheres_parameters.f90 \
       $(SRCDIR)/compute_wca_diameter.f90 \
       $(SRCDIR)/prepare_minimizer.f90 \
       $(SRCDIR)/compute_rdf.f90 \
       $(SRCDIR)/compute_planar_density.f90 \
       $(SRCDIR)/write_to_cube_file.f90 \
       $(SRCDIR)/vext_total_sum.f90 \
       $(SRCDIR)/compute_vcoul_as_sum_of_pointcharges.f90 \
       $(SRCDIR)/compute_vext_hard_sphere.f90 \
       $(SRCDIR)/external_potential_hard_walls.f90 \
       $(SRCDIR)/compute_vext_hard_cylinder.f90 \
       $(SRCDIR)/compute_purely_repulsive_potential.f90 \
       $(SRCDIR)/compute_vext_perso.f90 \
       $(SRCDIR)/energy_and_gradient.f90 \
       $(SRCDIR)/energy_hard_sphere_fmt.f90 \
       $(SRCDIR)/energy_cs_hard_sphere.f90 \
       $(SRCDIR)/energy_external.f90 \
       $(SRCDIR)/energy_ideal.f90 \
       $(SRCDIR)/energy_polarization.f90 \
       $(SRCDIR)/energy_polarization_myway.f90 \
       $(SRCDIR)/energy_threebody.f90 \
       $(SRCDIR)/energy_threebody_faster.f90 \
       $(SRCDIR)/lennard_jones_perturbation_to_hard_spheres.f90 \
       $(SRCDIR)/cs_from_dcf.f90 \
       $(SRCDIR)/cs_plus_hydro.f90 \
       $(SRCDIR)/energy_hydro.f90 \
       $(SRCDIR)/print_header.f90




# symbol '@' in front of a line makes it silent. Otherwise it is printed in terminal when called

 all: $(OBJS)
	 $(FC) $(FCFLAGS) $(LDFLAGS) -o $(EXE) $(OBJS)

 optim: $(OBJS)
	 $(FC) $(FCFLAGS) $(LDFLAGS) $(OPTIM) -o $(EXE) $(OBJS)

 debug: $(OBJS)
	 $(FC) $(FCFLAGS) $(LDFLAGS) $(DEBUGFLAGS) -o $(EXE) $(OBJS)

 clean:
	rm -vf gmon.out $(EXE) $(MODDIR)/* 

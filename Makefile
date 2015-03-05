# ——————————————— Program property ———————————————

EXE = mdft
MKFLAGS = $(NULLSTRING)

# ——————————————— Variables of locations ———————————————

MODDIR = mod
SRCDIR = src
OBJDIR = obj

# _______________ Libraries and other folders __________

FFTW_INCLUDES  = -I/usr/local/include -I/usr/include
FFTW_LIBRARIES = -L/usr/local/lib -L/usr/lib -lfftw3 -lfftw3_threads
NLOPT_INCLUDES = -I/usr/local/include
NLOPT_LIBRARIES = -lnlopt -lm

# ——————————————— Fortran compiler ———————————————

FC = gfortran

# ——————————————— Compiling options ———————————————

FCFLAGS = -J$(MODDIR) -I$(MODDIR) $(FFTW_INCLUDES) $(NLOPT_INCLUDES) -Wfatal-errors -fdiagnostics-color=auto
LDFLAGS = $(FFTW_LIBRARIES) $(NLOPT_LIBRARIES)

# For POINCARE:
# FCFLAGS = -J $(MODDIR) -I $(MODDIR) -I $(FFTW_INC_DIR) -Wfatal-errors # -fdiagnostics-color=auto
# LDFLAGS = -L $(FFTW_LIB_DIR) -lfftw3_threads -lfftw3

DEBUG = -Og -g -fimplicit-none -fbacktrace -pedantic -fwhole-file -Wline-truncation -Wcharacter-truncation -Wsurprising -Waliasing -fbounds-check -pg -frecursive -fcheck=all -Wall 
# -g turns on debugging
# -p turns on profiling
# -Wextra turns on extra warning. It is extremely verbose.
# -fcheck-array-temporaries -Warray-temporaries -Wconversion -Wimplicit-interface

OPTIM = -O3 -march=native -ffast-math -funroll-loops -ftree-vectorizer-verbose=6
# FOR BACKUP : -march=native -O3 -ffast-math -funroll-loops   VERY AGRESSIVE
# -fopenmp for OPENMP support

# ——————————————— Files to compile ———————————————

FOBJ = $(OBJDIR)/allocate_from_input.o\
	$(OBJDIR)/chargeDensityAndMolecularPolarizationOfASolventMoleculeAtOrigin.o\
	$(OBJDIR)/compute_planar_density.o\
	$(OBJDIR)/compute_purely_repulsive_potential.o\
	$(OBJDIR)/output_rdf.o\
	$(OBJDIR)/compute_vcoul_as_sum_of_pointcharges.o\
	$(OBJDIR)/compute_vext_hard_cylinder.o\
	$(OBJDIR)/compute_vext_hard_sphere.o\
	$(OBJDIR)/compute_vext_perso.o\
	$(OBJDIR)/compute_wca_diameter.o\
	$(OBJDIR)/compute_z_density.o\
	$(OBJDIR)/convert_coordinate_into_icg.o\
	$(OBJDIR)/cs_of_k_hard_sphere.o\
	$(OBJDIR)/energy_and_gradient.o\
	$(OBJDIR)/adhoc_corrections_to_gsolv.o\
	$(OBJDIR)/energy_ck_angular.o\
	$(OBJDIR)/energy_external.o\
	$(OBJDIR)/energy_fmt.o\
	$(OBJDIR)/energy_hydro.o\
	$(OBJDIR)/energy_ideal.o\
	$(OBJDIR)/energy_nn_cs_plus_nbar.o\
	$(OBJDIR)/energy_polarization_dipol.o\
	$(OBJDIR)/energy_polarization_multi.o\
	$(OBJDIR)/energy_polarization_multi_with_nccoupling.o\
	$(OBJDIR)/energy_threebody_faster.o\
	$(OBJDIR)/external_potential_hard_walls.o\
	$(OBJDIR)/find_equilibrium_density.o\
	$(OBJDIR)/get_final_density.o\
	$(OBJDIR)/get_final_polarization.o\
	$(OBJDIR)/get_gauss_legendre_integration_roots_and_weights.o\
	$(OBJDIR)/get_lebedev_integration_roots_and_weights.o\
	$(OBJDIR)/init_density.o\
	$(OBJDIR)/init_external_potential.o\
	$(OBJDIR)/init_fftw3_plans.o\
	$(OBJDIR)/init_simu.o\
	$(OBJDIR)/init_solvent_polarization.o\
	$(OBJDIR)/lennard_jones_perturbation_to_hard_spheres.o\
	$(OBJDIR)/main.o\
	$(OBJDIR)/mean_over_orientations.o\
	$(OBJDIR)/module_constants.o\
	$(OBJDIR)/module_dcf.o\
	$(OBJDIR)/module_external_potential.o\
	$(OBJDIR)/module_fft.o\
	$(OBJDIR)/module_fastPoissonSolver.o\
	$(OBJDIR)/module_hardspheres.o\
	$(OBJDIR)/module_input.o\
	$(OBJDIR)/module_mathematica.o\
	$(OBJDIR)/module_minimizer.o\
	$(OBJDIR)/module_mod_lj.o\
	$(OBJDIR)/module_periodic_table.o\
	$(OBJDIR)/module_precision_kinds.o\
	$(OBJDIR)/module_quadrature.o\
	$(OBJDIR)/module_solute_geometry.o\
	$(OBJDIR)/module_system.o\
	$(OBJDIR)/print_header.o\
	$(OBJDIR)/print_input_in_output_folder.o\
	$(OBJDIR)/print_supercell_xsf.o\
	$(OBJDIR)/process_output.o\
	$(OBJDIR)/put_input_in_character_array.o\
	$(OBJDIR)/read_solute.o\
	$(OBJDIR)/read_solvent.o\
	$(OBJDIR)/soluteChargeDensityFromSoluteChargeCoordinates.o\
	$(OBJDIR)/vext_total_sum.o\
	$(OBJDIR)/write_to_cube_file.o\
	$(OBJDIR)/finalize_simu.o\
	$(OBJDIR)/output_g-r-theta.o\
	$(OBJDIR)/output_gsitesite.o\
	$(OBJDIR)/module_time.o\
	$(OBJDIR)/energy_cs.o

# ——————————————— Global rules ———————————————

all: $(EXE)
	@ $(FC) $(FCFLAGS) $(MKFLAGS) -o $(EXE) $(FOBJ) $(LDFLAGS)

optim: MKFLAGS = $(OPTIM)

optim: $(EXE)
	@ $(FC) $(FCFLAGS) $(MKFLAGS) -o $(EXE) $(FOBJ) $(LDFLAGS)

debug: MKFLAGS = $(DEBUG)

debug: $(EXE)
	$(FC) $(FCFLAGS) $(MKFLAGS) -o $(EXE) $(FOBJ) $(LDFLAGS)

clean:
	-@ rm -vf gmon.out $(EXE) $(MODDIR)/* $(OBJDIR)/* >/dev/null 2>/dev/null

# ——————————————— Pattern rules ———————————————

$(OBJDIR)/%.o : $(SRCDIR)/%.f90
	@ $(FC) $(FCFLAGS) $(MKFLAGS) -c $< -o $@

# For GNU make, *.f90 cannot be compiled automatically.

# ——————————————— Dependence rules ———————————————

$(EXE): $(FOBJ)

$(OBJDIR)/allocate_from_input.o:\
	$(SRCDIR)/allocate_from_input.f90\
	$(OBJDIR)/module_quadrature.o\
	$(OBJDIR)/module_constants.o\
	$(OBJDIR)/module_system.o\
	$(OBJDIR)/module_input.o\
	$(OBJDIR)/module_precision_kinds.o

$(OBJDIR)/chargeDensityAndMolecularPolarizationOfASolventMoleculeAtOrigin.o:\
	$(SRCDIR)/chargeDensityAndMolecularPolarizationOfASolventMoleculeAtOrigin.f90\
	$(OBJDIR)/module_fft.o\
	$(OBJDIR)/module_quadrature.o\
	$(OBJDIR)/module_system.o\
	$(OBJDIR)/module_constants.o\
	$(OBJDIR)/module_mathematica.o\
	$(OBJDIR)/module_precision_kinds.o

$(OBJDIR)/compute_planar_density.o:\
	$(SRCDIR)/compute_planar_density.f90\
	$(OBJDIR)/module_system.o\
	$(OBJDIR)/module_precision_kinds.o

$(OBJDIR)/compute_purely_repulsive_potential.o:\
	$(SRCDIR)/compute_purely_repulsive_potential.f90\
	$(OBJDIR)/module_constants.o\
	$(OBJDIR)/module_quadrature.o\
	$(OBJDIR)/module_external_potential.o\
	$(OBJDIR)/module_system.o\
	$(OBJDIR)/module_input.o\
	$(OBJDIR)/module_precision_kinds.o

$(OBJDIR)/output_rdf.o:\
	$(SRCDIR)/output_rdf.f90\
	$(OBJDIR)/module_system.o\
	$(OBJDIR)/module_precision_kinds.o\
	$(OBJDIR)/module_mathematica.o

$(OBJDIR)/compute_vcoul_as_sum_of_pointcharges.o:\
	$(SRCDIR)/compute_vcoul_as_sum_of_pointcharges.f90\
	$(OBJDIR)/module_input.o\
	$(OBJDIR)/module_quadrature.o\
	$(OBJDIR)/module_external_potential.o\
	$(OBJDIR)/module_constants.o\
	$(OBJDIR)/module_system.o\
	$(OBJDIR)/module_precision_kinds.o

$(OBJDIR)/compute_vext_hard_cylinder.o:\
	$(SRCDIR)/compute_vext_hard_cylinder.f90\
	$(OBJDIR)/module_hardspheres.o\
	$(OBJDIR)/module_quadrature.o\
	$(OBJDIR)/module_external_potential.o\
	$(OBJDIR)/module_system.o\
	$(OBJDIR)/module_input.o\
	$(OBJDIR)/module_precision_kinds.o

$(OBJDIR)/compute_vext_hard_sphere.o:\
	$(SRCDIR)/compute_vext_hard_sphere.f90\
	$(OBJDIR)/module_hardspheres.o\
	$(OBJDIR)/module_quadrature.o\
	$(OBJDIR)/module_external_potential.o\
	$(OBJDIR)/module_system.o\
	$(OBJDIR)/module_input.o\
	$(OBJDIR)/module_precision_kinds.o

$(OBJDIR)/compute_vext_perso.o:\
	$(SRCDIR)/compute_vext_perso.f90\
	$(OBJDIR)/module_hardspheres.o\
	$(OBJDIR)/module_external_potential.o\
	$(OBJDIR)/module_system.o\
	$(OBJDIR)/module_input.o\
	$(OBJDIR)/module_precision_kinds.o

$(OBJDIR)/compute_wca_diameter.o:\
	$(SRCDIR)/compute_wca_diameter.f90\
	$(OBJDIR)/module_constants.o\
	$(OBJDIR)/module_precision_kinds.o

$(OBJDIR)/compute_z_density.o:\
	$(SRCDIR)/compute_z_density.f90\
	$(OBJDIR)/module_mathematica.o\
	$(OBJDIR)/module_system.o\
	$(OBJDIR)/module_precision_kinds.o

$(OBJDIR)/convert_coordinate_into_icg.o:\
	$(SRCDIR)/convert_coordinate_into_icg.f90\
	$(OBJDIR)/module_quadrature.o\
	$(OBJDIR)/module_system.o\
	$(OBJDIR)/module_precision_kinds.o

$(OBJDIR)/cs_of_k_hard_sphere.o:\
	$(SRCDIR)/cs_of_k_hard_sphere.f90\
	$(OBJDIR)/module_hardspheres.o\
	$(OBJDIR)/module_dcf.o\
	$(OBJDIR)/module_constants.o\
	$(OBJDIR)/module_input.o\
	$(OBJDIR)/module_system.o\
	$(OBJDIR)/module_precision_kinds.o

$(OBJDIR)/adhoc_corrections_to_gsolv.o:\
	$(SRCDIR)/adhoc_corrections_to_gsolv.f90\
	$(OBJDIR)/module_precision_kinds.o\
	$(OBJDIR)/module_system.o\
	$(OBJDIR)/module_minimizer.o\
	$(OBJDIR)/get_final_density.o

$(OBJDIR)/energy_and_gradient.o:\
	$(SRCDIR)/energy_and_gradient.f90\
	$(OBJDIR)/module_system.o\
	$(OBJDIR)/module_minimizer.o\
	$(OBJDIR)/module_input.o\
	$(OBJDIR)/module_precision_kinds.o

$(OBJDIR)/energy_ck_angular.o:\
	$(SRCDIR)/energy_ck_angular.f90\
	$(OBJDIR)/module_dcf.o\
	$(OBJDIR)/module_input.o\
	$(OBJDIR)/module_fft.o\
	$(OBJDIR)/module_constants.o\
	$(OBJDIR)/module_minimizer.o\
	$(OBJDIR)/module_quadrature.o\
	$(OBJDIR)/module_system.o\
	$(OBJDIR)/module_precision_kinds.o

$(OBJDIR)/energy_external.o:\
	$(SRCDIR)/energy_external.f90\
	$(OBJDIR)/module_external_potential.o\
	$(OBJDIR)/module_minimizer.o\
	$(OBJDIR)/module_quadrature.o\
	$(OBJDIR)/module_system.o\
	$(OBJDIR)/module_precision_kinds.o\
	$(OBJDIR)/module_input.o

$(OBJDIR)/energy_fmt.o:\
	$(SRCDIR)/energy_fmt.f90\
	$(OBJDIR)/module_hardspheres.o\
	$(OBJDIR)/module_input.o\
	$(OBJDIR)/module_fft.o\
	$(OBJDIR)/module_constants.o\
	$(OBJDIR)/module_minimizer.o\
	$(OBJDIR)/module_quadrature.o\
	$(OBJDIR)/module_system.o\
	$(OBJDIR)/module_precision_kinds.o

$(OBJDIR)/energy_hydro.o:\
	$(SRCDIR)/energy_hydro.f90\
	$(OBJDIR)/module_input.o\
	$(OBJDIR)/module_fft.o\
	$(OBJDIR)/module_quadrature.o\
	$(OBJDIR)/module_minimizer.o\
	$(OBJDIR)/module_constants.o\
	$(OBJDIR)/module_dcf.o\
	$(OBJDIR)/module_system.o\
	$(OBJDIR)/module_precision_kinds.o\
	$(OBJDIR)/module_mathematica.o

$(OBJDIR)/energy_ideal.o:\
	$(SRCDIR)/energy_ideal.f90\
	$(OBJDIR)/module_quadrature.o\
	$(OBJDIR)/module_system.o\
	$(OBJDIR)/module_minimizer.o\
	$(OBJDIR)/module_precision_kinds.o

$(OBJDIR)/energy_nn_cs_plus_nbar.o:\
	$(SRCDIR)/energy_nn_cs_plus_nbar.f90\
	$(OBJDIR)/module_input.o\
	$(OBJDIR)/module_fft.o\
	$(OBJDIR)/module_quadrature.o\
	$(OBJDIR)/module_minimizer.o\
	$(OBJDIR)/module_constants.o\
	$(OBJDIR)/module_dcf.o\
	$(OBJDIR)/module_system.o\
	$(OBJDIR)/module_mathematica.o\
	$(OBJDIR)/module_precision_kinds.o

$(OBJDIR)/energy_polarization_dipol.o:\
	$(SRCDIR)/energy_polarization_dipol.f90\
	$(OBJDIR)/module_dcf.o\
	$(OBJDIR)/module_input.o\
	$(OBJDIR)/module_fft.o\
	$(OBJDIR)/module_constants.o\
	$(OBJDIR)/module_minimizer.o\
	$(OBJDIR)/module_quadrature.o\
	$(OBJDIR)/module_system.o\
	$(OBJDIR)/module_mathematica.o\
	$(OBJDIR)/module_precision_kinds.o

$(OBJDIR)/energy_polarization_multi.o:\
	$(SRCDIR)/energy_polarization_multi.f90\
	$(OBJDIR)/module_input.o\
	$(OBJDIR)/module_fft.o\
	$(OBJDIR)/module_constants.o\
	$(OBJDIR)/module_minimizer.o\
	$(OBJDIR)/module_quadrature.o\
	$(OBJDIR)/module_dcf.o\
	$(OBJDIR)/module_system.o\
	$(OBJDIR)/module_precision_kinds.o

$(OBJDIR)/energy_polarization_multi_with_nccoupling.o:\
	$(SRCDIR)/energy_polarization_multi_with_nccoupling.f90\
	$(OBJDIR)/module_input.o\
	$(OBJDIR)/module_fft.o\
	$(OBJDIR)/module_constants.o\
	$(OBJDIR)/module_minimizer.o\
	$(OBJDIR)/module_quadrature.o\
	$(OBJDIR)/module_dcf.o\
	$(OBJDIR)/module_system.o\
	$(OBJDIR)/module_precision_kinds.o

$(OBJDIR)/energy_threebody_faster.o:\
	$(SRCDIR)/energy_threebody_faster.f90\
	$(OBJDIR)/module_fft.o\
	$(OBJDIR)/module_minimizer.o\
	$(OBJDIR)/module_system.o\
	$(OBJDIR)/module_quadrature.o\
	$(OBJDIR)/module_constants.o\
	$(OBJDIR)/module_input.o\
	$(OBJDIR)/module_precision_kinds.o

$(OBJDIR)/external_potential_hard_walls.o:\
	$(SRCDIR)/external_potential_hard_walls.f90\
	$(OBJDIR)/module_hardspheres.o\
	$(OBJDIR)/module_quadrature.o\
	$(OBJDIR)/module_external_potential.o\
	$(OBJDIR)/module_system.o\
	$(OBJDIR)/module_input.o\
	$(OBJDIR)/module_precision_kinds.o

$(OBJDIR)/find_equilibrium_density.o:\
	$(SRCDIR)/find_equilibrium_density.f90\
	$(OBJDIR)/module_minimizer.o\
	$(OBJDIR)/module_precision_kinds.o\
	$(OBJDIR)/adhoc_corrections_to_gsolv.o

$(OBJDIR)/get_final_density.o:\
	$(SRCDIR)/get_final_density.f90\
	$(OBJDIR)/module_input.o\
	$(OBJDIR)/module_fft.o\
	$(OBJDIR)/module_quadrature.o\
	$(OBJDIR)/module_minimizer.o\
	$(OBJDIR)/module_constants.o\
	$(OBJDIR)/module_system.o\
	$(OBJDIR)/module_precision_kinds.o

$(OBJDIR)/get_final_polarization.o:\
	$(SRCDIR)/get_final_polarization.f90\
	$(OBJDIR)/module_input.o\
	$(OBJDIR)/module_quadrature.o\
	$(OBJDIR)/module_minimizer.o\
	$(OBJDIR)/module_constants.o\
	$(OBJDIR)/module_system.o\
	$(OBJDIR)/module_precision_kinds.o

$(OBJDIR)/get_gauss_legendre_integration_roots_and_weights.o:\
	$(SRCDIR)/get_gauss_legendre_integration_roots_and_weights.f90\
	$(OBJDIR)/module_precision_kinds.o

$(OBJDIR)/get_lebedev_integration_roots_and_weights.o:\
	$(SRCDIR)/get_lebedev_integration_roots_and_weights.f90\
	$(OBJDIR)/module_precision_kinds.o

$(OBJDIR)/init_density.o:\
	$(SRCDIR)/init_density.f90\
	$(OBJDIR)/module_mathematica.o\
	$(OBJDIR)/module_input.o\
	$(OBJDIR)/module_external_potential.o\
	$(OBJDIR)/module_minimizer.o\
	$(OBJDIR)/module_quadrature.o\
	$(OBJDIR)/module_system.o\
	$(OBJDIR)/module_precision_kinds.o

$(OBJDIR)/init_external_potential.o:\
	$(SRCDIR)/init_external_potential.f90\
	$(OBJDIR)/module_constants.o\
	$(OBJDIR)/module_quadrature.o\
	$(OBJDIR)/module_mod_lj.o\
	$(OBJDIR)/module_external_potential.o\
	$(OBJDIR)/module_system.o\
	$(OBJDIR)/module_input.o\
	$(OBJDIR)/module_precision_kinds.o\
	$(OBJDIR)/module_fastPoissonSolver.o

$(OBJDIR)/init_fftw3_plans.o:\
	$(SRCDIR)/init_fftw3_plans.f90\
	$(OBJDIR)/module_system.o\
	$(OBJDIR)/module_fft.o

$(OBJDIR)/init_simu.o:\
	$(SRCDIR)/init_simu.f90\
	$(OBJDIR)/module_input.o\
	$(OBJDIR)/module_dcf.o\
	$(OBJDIR)/module_fft.o\
	$(OBJDIR)/module_quadrature.o\
	$(OBJDIR)/module_minimizer.o\
	$(OBJDIR)/module_hardspheres.o

$(OBJDIR)/init_solvent_polarization.o:\
	$(SRCDIR)/init_solvent_polarization.f90\
	$(OBJDIR)/module_input.o\
	$(OBJDIR)/module_quadrature.o

$(OBJDIR)/lennard_jones_perturbation_to_hard_spheres.o:\
	$(SRCDIR)/lennard_jones_perturbation_to_hard_spheres.f90\
	$(OBJDIR)/module_fft.o\
	$(OBJDIR)/module_constants.o\
	$(OBJDIR)/module_precision_kinds.o\
	$(OBJDIR)/module_minimizer.o\
	$(OBJDIR)/module_quadrature.o\
	$(OBJDIR)/module_system.o

$(OBJDIR)/main.o:\
	$(SRCDIR)/main.f90\
	$(OBJDIR)/module_precision_kinds.o

$(OBJDIR)/mean_over_orientations.o:\
	$(SRCDIR)/mean_over_orientations.f90\
	$(OBJDIR)/module_quadrature.o\
	$(OBJDIR)/module_system.o\
	$(OBJDIR)/module_precision_kinds.o

$(OBJDIR)/module_constants.o:\
	$(SRCDIR)/module_constants.f90\
	$(OBJDIR)/module_precision_kinds.o

$(OBJDIR)/module_dcf.o:\
	$(SRCDIR)/module_dcf.f90\
	$(OBJDIR)/module_system.o\
	$(OBJDIR)/module_constants.o\
	$(OBJDIR)/module_fft.o\
	$(OBJDIR)/module_quadrature.o\
	$(OBJDIR)/module_input.o\
	$(OBJDIR)/module_precision_kinds.o\
	$(OBJDIR)/module_mathematica.o

$(OBJDIR)/module_external_potential.o:\
	$(SRCDIR)/module_external_potential.f90\
	$(OBJDIR)/module_constants.o\
	$(OBJDIR)/module_system.o\
	$(OBJDIR)/module_precision_kinds.o

$(OBJDIR)/module_fft.o:\
	$(SRCDIR)/module_fft.f90\
	$(OBJDIR)/module_system.o\
	$(OBJDIR)/module_precision_kinds.o\
	$(OBJDIR)/module_input.o

$(OBJDIR)/module_hardspheres.o:\
	$(SRCDIR)/module_hardspheres.f90\
	$(OBJDIR)/module_fft.o\
	$(OBJDIR)/module_system.o\
	$(OBJDIR)/module_constants.o\
	$(OBJDIR)/module_precision_kinds.o\
	$(OBJDIR)/module_dcf.o

$(OBJDIR)/module_input.o:\
	$(SRCDIR)/module_input.f90\
	$(OBJDIR)/module_precision_kinds.o

$(OBJDIR)/module_minimizer.o:\
	$(SRCDIR)/module_minimizer.f90\
	$(OBJDIR)/module_precision_kinds.o\
	$(OBJDIR)/module_system.o\
	$(OBJDIR)/module_quadrature.o

$(OBJDIR)/module_mod_lj.o:\
	$(SRCDIR)/module_mod_lj.f90\
	$(OBJDIR)/module_quadrature.o\
	$(OBJDIR)/module_input.o\
	$(OBJDIR)/module_constants.o\
	$(OBJDIR)/module_system.o\
	$(OBJDIR)/module_precision_kinds.o\
	$(OBJDIR)/module_external_potential.o

$(OBJDIR)/module_periodic_table.o:\
	$(SRCDIR)/module_periodic_table.f90\
	$(OBJDIR)/module_precision_kinds.o

$(OBJDIR)/module_precision_kinds.o:\
	$(SRCDIR)/module_precision_kinds.f90

$(OBJDIR)/module_quadrature.o:\
	$(SRCDIR)/module_quadrature.f90\
	$(OBJDIR)/module_input.o\
	$(OBJDIR)/module_constants.o\
	$(OBJDIR)/module_precision_kinds.o

$(OBJDIR)/module_solute_geometry.o:\
	$(SRCDIR)/module_solute_geometry.f90\
	$(OBJDIR)/module_system.o\
	$(OBJDIR)/module_precision_kinds.o

$(OBJDIR)/module_system.o:\
	$(SRCDIR)/module_system.f90\
	$(OBJDIR)/module_precision_kinds.o

$(OBJDIR)/print_header.o:\
	$(SRCDIR)/print_header.f90

$(OBJDIR)/print_input_in_output_folder.o:\
	$(SRCDIR)/print_input_in_output_folder.f90\
	$(OBJDIR)/module_precision_kinds.o\
	$(OBJDIR)/module_input.o

$(OBJDIR)/print_supercell_xsf.o:\
	$(SRCDIR)/print_supercell_xsf.f90\
	$(OBJDIR)/module_system.o\
	$(OBJDIR)/module_precision_kinds.o

$(OBJDIR)/process_output.o:\
	$(SRCDIR)/process_output.f90\
	$(OBJDIR)/module_minimizer.o\
	$(OBJDIR)/module_constants.o\
	$(OBJDIR)/module_solute_geometry.o\
	$(OBJDIR)/module_input.o\
	$(OBJDIR)/module_system.o\
	$(OBJDIR)/module_precision_kinds.o\
	$(OBJDIR)/module_hardspheres.o

$(OBJDIR)/put_input_in_character_array.o:\
	$(SRCDIR)/put_input_in_character_array.f90\
	$(OBJDIR)/module_input.o\
	$(OBJDIR)/module_precision_kinds.o

$(OBJDIR)/read_solute.o:\
	$(SRCDIR)/read_solute.f90\
	$(OBJDIR)/module_periodic_table.o\
	$(OBJDIR)/module_input.o\
	$(OBJDIR)/module_system.o\
	$(OBJDIR)/module_precision_kinds.o

$(OBJDIR)/read_solvent.o:\
	$(SRCDIR)/read_solvent.f90\
	$(OBJDIR)/module_system.o\
	$(OBJDIR)/module_mathematica.o\
	$(OBJDIR)/module_precision_kinds.o

$(OBJDIR)/soluteChargeDensityFromSoluteChargeCoordinates.o:\
	$(SRCDIR)/soluteChargeDensityFromSoluteChargeCoordinates.f90\
	$(OBJDIR)/module_input.o\
	$(OBJDIR)/module_system.o\
	$(OBJDIR)/module_precision_kinds.o

$(OBJDIR)/vext_total_sum.o:\
	$(SRCDIR)/vext_total_sum.f90\
	$(OBJDIR)/module_input.o\
	$(OBJDIR)/module_external_potential.o\
	$(OBJDIR)/module_constants.o\
	$(OBJDIR)/module_quadrature.o\
	$(OBJDIR)/module_system.o\
	$(OBJDIR)/module_precision_kinds.o

$(OBJDIR)/write_to_cube_file.o:\
	$(SRCDIR)/write_to_cube_file.f90\
	$(OBJDIR)/module_precision_kinds.o\
	$(OBJDIR)/module_system.o

$(OBJDIR)/module_fastPoissonSolver.o:\
	$(SRCDIR)/module_fastPoissonSolver.f90\
	$(OBJDIR)/module_precision_kinds.o\
	$(OBJDIR)/module_quadrature.o\
	$(OBJDIR)/module_system.o

$(OBJDIR)/finalize_simu.o:\
	$(SRCDIR)/finalize_simu.f90\
	$(OBJDIR)/module_fft.o

$(OBJDIR)/module_mathematica.o:\
	$(SRCDIR)/module_mathematica.f90\
	$(OBJDIR)/module_precision_kinds.o

$(OBJDIR)/output_g-r-theta.o:\
	$(SRCDIR)/output_g-r-theta.f90\
	$(OBJDIR)/module_precision_kinds.o\
	$(OBJDIR)/module_system.o\
	$(OBJDIR)/module_quadrature.o\
	$(OBJDIR)/module_constants.o\
	$(OBJDIR)/module_minimizer.o\
	$(OBJDIR)/module_mathematica.o

$(OBJDIR)/output_gsitesite.o:\
	$(SRCDIR)/output_gsitesite.f90\
	$(OBJDIR)/module_precision_kinds.o\
	$(OBJDIR)/module_system.o\
	$(OBJDIR)/module_quadrature.o\
	$(OBJDIR)/module_constants.o\
	$(OBJDIR)/module_minimizer.o\
	$(OBJDIR)/module_mathematica.o

$(OBJDIR)/module_time.o:\
	$(SRCDIR)/module_time.f90

$(OBJDIR)/energy_cs.o:\
	$(SRCDIR)/energy_cs.f90\
	$(OBJDIR)/module_precision_kinds.o\
	$(OBJDIR)/module_system.o\
	$(OBJDIR)/module_quadrature.o\
	$(OBJDIR)/module_minimizer.o\
	$(OBJDIR)/module_fft.o\
	$(OBJDIR)/module_dcf.o\
	$(OBJDIR)/module_mathematica.o

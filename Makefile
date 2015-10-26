# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.0

# Default target executed when no arguments are given to make.
default_target: all
.PHONY : default_target

# Allow only one "make -f Makefile2" at a time, but pass parallelism.
.NOTPARALLEL:
.PHONY : .NOTPARALLEL

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list

# Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/levesque/Recherche/src/MDFT-dev

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/levesque/Recherche/src/MDFT-dev

#=============================================================================
# Targets provided globally by CMake.

# Special rule for the target edit_cache
edit_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake cache editor..."
	/usr/bin/cmake-gui -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : edit_cache

# Special rule for the target edit_cache
edit_cache/fast: edit_cache
.PHONY : edit_cache/fast

# Special rule for the target install
install: preinstall
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Install the project..."
	/usr/bin/cmake -P cmake_install.cmake
.PHONY : install

# Special rule for the target install
install/fast: preinstall/fast
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Install the project..."
	/usr/bin/cmake -P cmake_install.cmake
.PHONY : install/fast

# Special rule for the target install/local
install/local: preinstall
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Installing only the local directory..."
	/usr/bin/cmake -DCMAKE_INSTALL_LOCAL_ONLY=1 -P cmake_install.cmake
.PHONY : install/local

# Special rule for the target install/local
install/local/fast: install/local
.PHONY : install/local/fast

# Special rule for the target install/strip
install/strip: preinstall
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Installing the project stripped..."
	/usr/bin/cmake -DCMAKE_INSTALL_DO_STRIP=1 -P cmake_install.cmake
.PHONY : install/strip

# Special rule for the target install/strip
install/strip/fast: install/strip
.PHONY : install/strip/fast

# Special rule for the target list_install_components
list_install_components:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Only default component available"
.PHONY : list_install_components

# Special rule for the target list_install_components
list_install_components/fast: list_install_components
.PHONY : list_install_components/fast

# Special rule for the target rebuild_cache
rebuild_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake to regenerate build system..."
	/usr/bin/cmake -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : rebuild_cache

# Special rule for the target rebuild_cache
rebuild_cache/fast: rebuild_cache
.PHONY : rebuild_cache/fast

# The main all target
all: cmake_check_build_system
	$(CMAKE_COMMAND) -E cmake_progress_start /home/levesque/Recherche/src/MDFT-dev/CMakeFiles /home/levesque/Recherche/src/MDFT-dev/CMakeFiles/progress.marks
	$(MAKE) -f CMakeFiles/Makefile2 all
	$(CMAKE_COMMAND) -E cmake_progress_start /home/levesque/Recherche/src/MDFT-dev/CMakeFiles 0
.PHONY : all

# The main clean target
clean:
	$(MAKE) -f CMakeFiles/Makefile2 clean
.PHONY : clean

# The main clean target
clean/fast: clean
.PHONY : clean/fast

# Prepare targets for installation.
preinstall: all
	$(MAKE) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall

# Prepare targets for installation.
preinstall/fast:
	$(MAKE) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall/fast

# clear depends
depend:
	$(CMAKE_COMMAND) -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 1
.PHONY : depend

#=============================================================================
# Target rules for targets named mdft-dev

# Build rule for target.
mdft-dev: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 mdft-dev
.PHONY : mdft-dev

# fast build rule for target.
mdft-dev/fast:
	$(MAKE) -f CMakeFiles/mdft-dev.dir/build.make CMakeFiles/mdft-dev.dir/build
.PHONY : mdft-dev/fast

src/compute_purely_repulsive_potential.o: src/compute_purely_repulsive_potential.f90.o
.PHONY : src/compute_purely_repulsive_potential.o

# target to build an object file
src/compute_purely_repulsive_potential.f90.o:
	$(MAKE) -f CMakeFiles/mdft-dev.dir/build.make CMakeFiles/mdft-dev.dir/src/compute_purely_repulsive_potential.f90.o
.PHONY : src/compute_purely_repulsive_potential.f90.o

src/compute_vcoul_as_sum_of_pointcharges.o: src/compute_vcoul_as_sum_of_pointcharges.f90.o
.PHONY : src/compute_vcoul_as_sum_of_pointcharges.o

# target to build an object file
src/compute_vcoul_as_sum_of_pointcharges.f90.o:
	$(MAKE) -f CMakeFiles/mdft-dev.dir/build.make CMakeFiles/mdft-dev.dir/src/compute_vcoul_as_sum_of_pointcharges.f90.o
.PHONY : src/compute_vcoul_as_sum_of_pointcharges.f90.o

src/compute_vext_hard_cylinder.o: src/compute_vext_hard_cylinder.f90.o
.PHONY : src/compute_vext_hard_cylinder.o

# target to build an object file
src/compute_vext_hard_cylinder.f90.o:
	$(MAKE) -f CMakeFiles/mdft-dev.dir/build.make CMakeFiles/mdft-dev.dir/src/compute_vext_hard_cylinder.f90.o
.PHONY : src/compute_vext_hard_cylinder.f90.o

src/compute_vext_hard_sphere.o: src/compute_vext_hard_sphere.f90.o
.PHONY : src/compute_vext_hard_sphere.o

# target to build an object file
src/compute_vext_hard_sphere.f90.o:
	$(MAKE) -f CMakeFiles/mdft-dev.dir/build.make CMakeFiles/mdft-dev.dir/src/compute_vext_hard_sphere.f90.o
.PHONY : src/compute_vext_hard_sphere.f90.o

src/compute_vext_perso.o: src/compute_vext_perso.f90.o
.PHONY : src/compute_vext_perso.o

# target to build an object file
src/compute_vext_perso.f90.o:
	$(MAKE) -f CMakeFiles/mdft-dev.dir/build.make CMakeFiles/mdft-dev.dir/src/compute_vext_perso.f90.o
.PHONY : src/compute_vext_perso.f90.o

src/compute_wca_diameter.o: src/compute_wca_diameter.f90.o
.PHONY : src/compute_wca_diameter.o

# target to build an object file
src/compute_wca_diameter.f90.o:
	$(MAKE) -f CMakeFiles/mdft-dev.dir/build.make CMakeFiles/mdft-dev.dir/src/compute_wca_diameter.f90.o
.PHONY : src/compute_wca_diameter.f90.o

src/compute_z_density.o: src/compute_z_density.f90.o
.PHONY : src/compute_z_density.o

# target to build an object file
src/compute_z_density.f90.o:
	$(MAKE) -f CMakeFiles/mdft-dev.dir/build.make CMakeFiles/mdft-dev.dir/src/compute_z_density.f90.o
.PHONY : src/compute_z_density.f90.o

src/convert_coordinate_into_icg.o: src/convert_coordinate_into_icg.f90.o
.PHONY : src/convert_coordinate_into_icg.o

# target to build an object file
src/convert_coordinate_into_icg.f90.o:
	$(MAKE) -f CMakeFiles/mdft-dev.dir/build.make CMakeFiles/mdft-dev.dir/src/convert_coordinate_into_icg.f90.o
.PHONY : src/convert_coordinate_into_icg.f90.o

src/cs_of_k_hard_sphere.o: src/cs_of_k_hard_sphere.f90.o
.PHONY : src/cs_of_k_hard_sphere.o

# target to build an object file
src/cs_of_k_hard_sphere.f90.o:
	$(MAKE) -f CMakeFiles/mdft-dev.dir/build.make CMakeFiles/mdft-dev.dir/src/cs_of_k_hard_sphere.f90.o
.PHONY : src/cs_of_k_hard_sphere.f90.o

src/energy_and_gradient.o: src/energy_and_gradient.f90.o
.PHONY : src/energy_and_gradient.o

# target to build an object file
src/energy_and_gradient.f90.o:
	$(MAKE) -f CMakeFiles/mdft-dev.dir/build.make CMakeFiles/mdft-dev.dir/src/energy_and_gradient.f90.o
.PHONY : src/energy_and_gradient.f90.o

src/energy_ck_angular.o: src/energy_ck_angular.f90.o
.PHONY : src/energy_ck_angular.o

# target to build an object file
src/energy_ck_angular.f90.o:
	$(MAKE) -f CMakeFiles/mdft-dev.dir/build.make CMakeFiles/mdft-dev.dir/src/energy_ck_angular.f90.o
.PHONY : src/energy_ck_angular.f90.o

src/energy_hydro.o: src/energy_hydro.f90.o
.PHONY : src/energy_hydro.o

# target to build an object file
src/energy_hydro.f90.o:
	$(MAKE) -f CMakeFiles/mdft-dev.dir/build.make CMakeFiles/mdft-dev.dir/src/energy_hydro.f90.o
.PHONY : src/energy_hydro.f90.o

src/energy_ideal_and_external.o: src/energy_ideal_and_external.f90.o
.PHONY : src/energy_ideal_and_external.o

# target to build an object file
src/energy_ideal_and_external.f90.o:
	$(MAKE) -f CMakeFiles/mdft-dev.dir/build.make CMakeFiles/mdft-dev.dir/src/energy_ideal_and_external.f90.o
.PHONY : src/energy_ideal_and_external.f90.o

src/energy_minimization.o: src/energy_minimization.f90.o
.PHONY : src/energy_minimization.o

# target to build an object file
src/energy_minimization.f90.o:
	$(MAKE) -f CMakeFiles/mdft-dev.dir/build.make CMakeFiles/mdft-dev.dir/src/energy_minimization.f90.o
.PHONY : src/energy_minimization.f90.o

src/energy_nn_cs_plus_nbar.o: src/energy_nn_cs_plus_nbar.f90.o
.PHONY : src/energy_nn_cs_plus_nbar.o

# target to build an object file
src/energy_nn_cs_plus_nbar.f90.o:
	$(MAKE) -f CMakeFiles/mdft-dev.dir/build.make CMakeFiles/mdft-dev.dir/src/energy_nn_cs_plus_nbar.f90.o
.PHONY : src/energy_nn_cs_plus_nbar.f90.o

src/energy_polarization_dipol.o: src/energy_polarization_dipol.f90.o
.PHONY : src/energy_polarization_dipol.o

# target to build an object file
src/energy_polarization_dipol.f90.o:
	$(MAKE) -f CMakeFiles/mdft-dev.dir/build.make CMakeFiles/mdft-dev.dir/src/energy_polarization_dipol.f90.o
.PHONY : src/energy_polarization_dipol.f90.o

src/energy_polarization_multi.o: src/energy_polarization_multi.f90.o
.PHONY : src/energy_polarization_multi.o

# target to build an object file
src/energy_polarization_multi.f90.o:
	$(MAKE) -f CMakeFiles/mdft-dev.dir/build.make CMakeFiles/mdft-dev.dir/src/energy_polarization_multi.f90.o
.PHONY : src/energy_polarization_multi.f90.o

src/energy_polarization_multi_with_nccoupling.o: src/energy_polarization_multi_with_nccoupling.f90.o
.PHONY : src/energy_polarization_multi_with_nccoupling.o

# target to build an object file
src/energy_polarization_multi_with_nccoupling.f90.o:
	$(MAKE) -f CMakeFiles/mdft-dev.dir/build.make CMakeFiles/mdft-dev.dir/src/energy_polarization_multi_with_nccoupling.f90.o
.PHONY : src/energy_polarization_multi_with_nccoupling.f90.o

src/energy_threebody_faster.o: src/energy_threebody_faster.f90.o
.PHONY : src/energy_threebody_faster.o

# target to build an object file
src/energy_threebody_faster.f90.o:
	$(MAKE) -f CMakeFiles/mdft-dev.dir/build.make CMakeFiles/mdft-dev.dir/src/energy_threebody_faster.f90.o
.PHONY : src/energy_threebody_faster.f90.o

src/external_potential_hard_walls.o: src/external_potential_hard_walls.f90.o
.PHONY : src/external_potential_hard_walls.o

# target to build an object file
src/external_potential_hard_walls.f90.o:
	$(MAKE) -f CMakeFiles/mdft-dev.dir/build.make CMakeFiles/mdft-dev.dir/src/external_potential_hard_walls.f90.o
.PHONY : src/external_potential_hard_walls.f90.o

src/get_final_polarization.o: src/get_final_polarization.f90.o
.PHONY : src/get_final_polarization.o

# target to build an object file
src/get_final_polarization.f90.o:
	$(MAKE) -f CMakeFiles/mdft-dev.dir/build.make CMakeFiles/mdft-dev.dir/src/get_final_polarization.f90.o
.PHONY : src/get_final_polarization.f90.o

src/init_external_potential.o: src/init_external_potential.f90.o
.PHONY : src/init_external_potential.o

# target to build an object file
src/init_external_potential.f90.o:
	$(MAKE) -f CMakeFiles/mdft-dev.dir/build.make CMakeFiles/mdft-dev.dir/src/init_external_potential.f90.o
.PHONY : src/init_external_potential.f90.o

src/main.o: src/main.f90.o
.PHONY : src/main.o

# target to build an object file
src/main.f90.o:
	$(MAKE) -f CMakeFiles/mdft-dev.dir/build.make CMakeFiles/mdft-dev.dir/src/main.f90.o
.PHONY : src/main.f90.o

src/module_constants.o: src/module_constants.f90.o
.PHONY : src/module_constants.o

# target to build an object file
src/module_constants.f90.o:
	$(MAKE) -f CMakeFiles/mdft-dev.dir/build.make CMakeFiles/mdft-dev.dir/src/module_constants.f90.o
.PHONY : src/module_constants.f90.o

src/module_dcf.o: src/module_dcf.f90.o
.PHONY : src/module_dcf.o

# target to build an object file
src/module_dcf.f90.o:
	$(MAKE) -f CMakeFiles/mdft-dev.dir/build.make CMakeFiles/mdft-dev.dir/src/module_dcf.f90.o
.PHONY : src/module_dcf.f90.o

src/module_density.o: src/module_density.f90.o
.PHONY : src/module_density.o

# target to build an object file
src/module_density.f90.o:
	$(MAKE) -f CMakeFiles/mdft-dev.dir/build.make CMakeFiles/mdft-dev.dir/src/module_density.f90.o
.PHONY : src/module_density.f90.o

src/module_external_potential.o: src/module_external_potential.f90.o
.PHONY : src/module_external_potential.o

# target to build an object file
src/module_external_potential.f90.o:
	$(MAKE) -f CMakeFiles/mdft-dev.dir/build.make CMakeFiles/mdft-dev.dir/src/module_external_potential.f90.o
.PHONY : src/module_external_potential.f90.o

src/module_fastPoissonSolver.o: src/module_fastPoissonSolver.f90.o
.PHONY : src/module_fastPoissonSolver.o

# target to build an object file
src/module_fastPoissonSolver.f90.o:
	$(MAKE) -f CMakeFiles/mdft-dev.dir/build.make CMakeFiles/mdft-dev.dir/src/module_fastPoissonSolver.f90.o
.PHONY : src/module_fastPoissonSolver.f90.o

src/module_fft.o: src/module_fft.f90.o
.PHONY : src/module_fft.o

# target to build an object file
src/module_fft.f90.o:
	$(MAKE) -f CMakeFiles/mdft-dev.dir/build.make CMakeFiles/mdft-dev.dir/src/module_fft.f90.o
.PHONY : src/module_fft.f90.o

src/module_geometry.o: src/module_geometry.f90.o
.PHONY : src/module_geometry.o

# target to build an object file
src/module_geometry.f90.o:
	$(MAKE) -f CMakeFiles/mdft-dev.dir/build.make CMakeFiles/mdft-dev.dir/src/module_geometry.f90.o
.PHONY : src/module_geometry.f90.o

src/module_grid.o: src/module_grid.f90.o
.PHONY : src/module_grid.o

# target to build an object file
src/module_grid.f90.o:
	$(MAKE) -f CMakeFiles/mdft-dev.dir/build.make CMakeFiles/mdft-dev.dir/src/module_grid.f90.o
.PHONY : src/module_grid.f90.o

src/module_hardspheres.o: src/module_hardspheres.f90.o
.PHONY : src/module_hardspheres.o

# target to build an object file
src/module_hardspheres.f90.o:
	$(MAKE) -f CMakeFiles/mdft-dev.dir/build.make CMakeFiles/mdft-dev.dir/src/module_hardspheres.f90.o
.PHONY : src/module_hardspheres.f90.o

src/module_init_simu.o: src/module_init_simu.f90.o
.PHONY : src/module_init_simu.o

# target to build an object file
src/module_init_simu.f90.o:
	$(MAKE) -f CMakeFiles/mdft-dev.dir/build.make CMakeFiles/mdft-dev.dir/src/module_init_simu.f90.o
.PHONY : src/module_init_simu.f90.o

src/module_input.o: src/module_input.f90.o
.PHONY : src/module_input.o

# target to build an object file
src/module_input.f90.o:
	$(MAKE) -f CMakeFiles/mdft-dev.dir/build.make CMakeFiles/mdft-dev.dir/src/module_input.f90.o
.PHONY : src/module_input.f90.o

src/module_lennardjones.o: src/module_lennardjones.f90.o
.PHONY : src/module_lennardjones.o

# target to build an object file
src/module_lennardjones.f90.o:
	$(MAKE) -f CMakeFiles/mdft-dev.dir/build.make CMakeFiles/mdft-dev.dir/src/module_lennardjones.f90.o
.PHONY : src/module_lennardjones.f90.o

src/module_mathematica.o: src/module_mathematica.f90.o
.PHONY : src/module_mathematica.o

# target to build an object file
src/module_mathematica.f90.o:
	$(MAKE) -f CMakeFiles/mdft-dev.dir/build.make CMakeFiles/mdft-dev.dir/src/module_mathematica.f90.o
.PHONY : src/module_mathematica.f90.o

src/module_minimizer.o: src/module_minimizer.f90.o
.PHONY : src/module_minimizer.o

# target to build an object file
src/module_minimizer.f90.o:
	$(MAKE) -f CMakeFiles/mdft-dev.dir/build.make CMakeFiles/mdft-dev.dir/src/module_minimizer.f90.o
.PHONY : src/module_minimizer.f90.o

src/module_periodic_table.o: src/module_periodic_table.f90.o
.PHONY : src/module_periodic_table.o

# target to build an object file
src/module_periodic_table.f90.o:
	$(MAKE) -f CMakeFiles/mdft-dev.dir/build.make CMakeFiles/mdft-dev.dir/src/module_periodic_table.f90.o
.PHONY : src/module_periodic_table.f90.o

src/module_postprocessing.o: src/module_postprocessing.f90.o
.PHONY : src/module_postprocessing.o

# target to build an object file
src/module_postprocessing.f90.o:
	$(MAKE) -f CMakeFiles/mdft-dev.dir/build.make CMakeFiles/mdft-dev.dir/src/module_postprocessing.f90.o
.PHONY : src/module_postprocessing.f90.o

src/module_precision_kinds.o: src/module_precision_kinds.f90.o
.PHONY : src/module_precision_kinds.o

# target to build an object file
src/module_precision_kinds.f90.o:
	$(MAKE) -f CMakeFiles/mdft-dev.dir/build.make CMakeFiles/mdft-dev.dir/src/module_precision_kinds.f90.o
.PHONY : src/module_precision_kinds.f90.o

src/module_quadrature.o: src/module_quadrature.f90.o
.PHONY : src/module_quadrature.o

# target to build an object file
src/module_quadrature.f90.o:
	$(MAKE) -f CMakeFiles/mdft-dev.dir/build.make CMakeFiles/mdft-dev.dir/src/module_quadrature.f90.o
.PHONY : src/module_quadrature.f90.o

src/module_solute.o: src/module_solute.f90.o
.PHONY : src/module_solute.o

# target to build an object file
src/module_solute.f90.o:
	$(MAKE) -f CMakeFiles/mdft-dev.dir/build.make CMakeFiles/mdft-dev.dir/src/module_solute.f90.o
.PHONY : src/module_solute.f90.o

src/module_solvent.o: src/module_solvent.f90.o
.PHONY : src/module_solvent.o

# target to build an object file
src/module_solvent.f90.o:
	$(MAKE) -f CMakeFiles/mdft-dev.dir/build.make CMakeFiles/mdft-dev.dir/src/module_solvent.f90.o
.PHONY : src/module_solvent.f90.o

src/module_system.o: src/module_system.f90.o
.PHONY : src/module_system.o

# target to build an object file
src/module_system.f90.o:
	$(MAKE) -f CMakeFiles/mdft-dev.dir/build.make CMakeFiles/mdft-dev.dir/src/module_system.f90.o
.PHONY : src/module_system.f90.o

src/module_time.o: src/module_time.f90.o
.PHONY : src/module_time.o

# target to build an object file
src/module_time.f90.o:
	$(MAKE) -f CMakeFiles/mdft-dev.dir/build.make CMakeFiles/mdft-dev.dir/src/module_time.f90.o
.PHONY : src/module_time.f90.o

src/mylbfgsb.o: src/mylbfgsb.f.o
.PHONY : src/mylbfgsb.o

# target to build an object file
src/mylbfgsb.f.o:
	$(MAKE) -f CMakeFiles/mdft-dev.dir/build.make CMakeFiles/mdft-dev.dir/src/mylbfgsb.f.o
.PHONY : src/mylbfgsb.f.o

src/output_g-r-theta.o: src/output_g-r-theta.f90.o
.PHONY : src/output_g-r-theta.o

# target to build an object file
src/output_g-r-theta.f90.o:
	$(MAKE) -f CMakeFiles/mdft-dev.dir/build.make CMakeFiles/mdft-dev.dir/src/output_g-r-theta.f90.o
.PHONY : src/output_g-r-theta.f90.o

src/output_gsitesite.o: src/output_gsitesite.f90.o
.PHONY : src/output_gsitesite.o

# target to build an object file
src/output_gsitesite.f90.o:
	$(MAKE) -f CMakeFiles/mdft-dev.dir/build.make CMakeFiles/mdft-dev.dir/src/output_gsitesite.f90.o
.PHONY : src/output_gsitesite.f90.o

src/output_rdf.o: src/output_rdf.f90.o
.PHONY : src/output_rdf.o

# target to build an object file
src/output_rdf.f90.o:
	$(MAKE) -f CMakeFiles/mdft-dev.dir/build.make CMakeFiles/mdft-dev.dir/src/output_rdf.f90.o
.PHONY : src/output_rdf.f90.o

src/vext_total_sum.o: src/vext_total_sum.f90.o
.PHONY : src/vext_total_sum.o

# target to build an object file
src/vext_total_sum.f90.o:
	$(MAKE) -f CMakeFiles/mdft-dev.dir/build.make CMakeFiles/mdft-dev.dir/src/vext_total_sum.f90.o
.PHONY : src/vext_total_sum.f90.o

src/write_to_cube_file.o: src/write_to_cube_file.f90.o
.PHONY : src/write_to_cube_file.o

# target to build an object file
src/write_to_cube_file.f90.o:
	$(MAKE) -f CMakeFiles/mdft-dev.dir/build.make CMakeFiles/mdft-dev.dir/src/write_to_cube_file.f90.o
.PHONY : src/write_to_cube_file.f90.o

# Help Target
help:
	@echo "The following are some of the valid targets for this Makefile:"
	@echo "... all (the default if no target is provided)"
	@echo "... clean"
	@echo "... depend"
	@echo "... edit_cache"
	@echo "... install"
	@echo "... install/local"
	@echo "... install/strip"
	@echo "... list_install_components"
	@echo "... mdft-dev"
	@echo "... rebuild_cache"
	@echo "... src/compute_purely_repulsive_potential.o"
	@echo "... src/compute_vcoul_as_sum_of_pointcharges.o"
	@echo "... src/compute_vext_hard_cylinder.o"
	@echo "... src/compute_vext_hard_sphere.o"
	@echo "... src/compute_vext_perso.o"
	@echo "... src/compute_wca_diameter.o"
	@echo "... src/compute_z_density.o"
	@echo "... src/convert_coordinate_into_icg.o"
	@echo "... src/cs_of_k_hard_sphere.o"
	@echo "... src/energy_and_gradient.o"
	@echo "... src/energy_ck_angular.o"
	@echo "... src/energy_hydro.o"
	@echo "... src/energy_ideal_and_external.o"
	@echo "... src/energy_minimization.o"
	@echo "... src/energy_nn_cs_plus_nbar.o"
	@echo "... src/energy_polarization_dipol.o"
	@echo "... src/energy_polarization_multi.o"
	@echo "... src/energy_polarization_multi_with_nccoupling.o"
	@echo "... src/energy_threebody_faster.o"
	@echo "... src/external_potential_hard_walls.o"
	@echo "... src/get_final_polarization.o"
	@echo "... src/init_external_potential.o"
	@echo "... src/main.o"
	@echo "... src/module_constants.o"
	@echo "... src/module_dcf.o"
	@echo "... src/module_density.o"
	@echo "... src/module_external_potential.o"
	@echo "... src/module_fastPoissonSolver.o"
	@echo "... src/module_fft.o"
	@echo "... src/module_geometry.o"
	@echo "... src/module_grid.o"
	@echo "... src/module_hardspheres.o"
	@echo "... src/module_init_simu.o"
	@echo "... src/module_input.o"
	@echo "... src/module_lennardjones.o"
	@echo "... src/module_mathematica.o"
	@echo "... src/module_minimizer.o"
	@echo "... src/module_periodic_table.o"
	@echo "... src/module_postprocessing.o"
	@echo "... src/module_precision_kinds.o"
	@echo "... src/module_quadrature.o"
	@echo "... src/module_solute.o"
	@echo "... src/module_solvent.o"
	@echo "... src/module_system.o"
	@echo "... src/module_time.o"
	@echo "... src/mylbfgsb.o"
	@echo "... src/output_g-r-theta.o"
	@echo "... src/output_gsitesite.o"
	@echo "... src/output_rdf.o"
	@echo "... src/vext_total_sum.o"
	@echo "... src/write_to_cube_file.o"
.PHONY : help



#=============================================================================
# Special targets to cleanup operation of make.

# Special rule to run CMake to check the build system integrity.
# No rule that depends on this can have commands that come from listfiles
# because they might be regenerated.
cmake_check_build_system:
	$(CMAKE_COMMAND) -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 0
.PHONY : cmake_check_build_system


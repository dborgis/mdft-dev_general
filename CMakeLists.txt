#
# Some help to write CMakeLists.cmake for fortran project:
#
# - http://www.cmake.org/Wiki/CMakeForFortranExample
# - http://www.ngssc.se/courses/specialized-courses/advanced-programming/make_cmake.pdf
# - https://github.com/SethMMorton/cmake_fortran_template
#

# require cmake 2.8.0
CMAKE_MINIMUM_REQUIRED(VERSION 2.8.0 FATAL_ERROR)

##################
project("mdft-dev")
##################

# Fortran is disabled by default, so we need to wake-up
# note: string given to enable_language is case-sensitive
enable_language(Fortran)

## guard against in-source builds
#if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
#  message(FATAL_ERROR "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there. You may need to remove CMakeCache.txt. ")
#endif()

# Build type :
# - RELEASE (with optimization flags)
# - DEBUG   (with debug flags)
# make sure that the default is RELEASE
if (NOT CMAKE_BUILD_TYPE)
    set (CMAKE_BUILD_TYPE  RELEASE  CACHE STRING
      "Choose the type of build, options are: None Debug Release."
      FORCE)
endif (NOT CMAKE_BUILD_TYPE)

# allow cmake to use custom modules located in cmake subdir
# (e.g. FindFFTW)
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

# external libraries dependencies
find_package(FFTW REQUIRED)

#
# check compiler (gfortran, ifort, ...)
# set FFLAGS depend on the compiler
# set where fortran module files will be located
#
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
message(STATUS "Using compiler ${Fortran_COMPILER_NAME}")

# Have the .mod files placed in the lib folder
SET(CMAKE_Fortran_MODULE_DIRECTORY "${PROJECT_BINARY_DIR}/mod")

if (Fortran_COMPILER_NAME MATCHES "gfortran.*")
  # GNU fortran compiler
  set (CMAKE_Fortran_FLAGS_RELEASE "-march=native -ffast-math -funroll-all-loops -fimplicit-none -ffree-line-length-none -O3")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-g -fimplicit-none -fbacktrace -pedantic -fwhole-file -Wline-truncation -Wcharacter-truncation -Wsurprising -Waliasing -fbounds-check -pg -frecursive -fcheck=all -Wall -ffpe-trap=zero,underflow,overflow")
elseif (Fortran_COMPILER_NAME MATCHES "ifort.*")
  # INTEL fortran compiler
  set (CMAKE_Fortran_FLAGS_RELEASE "-O3")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
elseif (Fortran_COMPILER_NAME MATCHES "pgf90")
  # PGI fortran compiler
  set (CMAKE_Fortran_FLAGS_RELEASE "-O3")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
endif()

#############################
# Now build mdft executable #
#############################

# set source files list
set(mdft_SRCS
  src/module_rotation.f90
  src/module_energy_cproj_mrso.f90
  src/module_energy_cdeltacd.f90
  src/module_thermo.f90
  src/module_debug.f90
  src/compute_vext_hard_cylinder.f90
  src/compute_vext_hard_sphere.f90
  src/compute_wca_diameter.f90
  src/compute_z_density.f90
  src/convert_coordinate_into_icg.f90
  src/cs_of_k_hard_sphere.f90
  src/module_energy_and_gradient.f90
  src/energy_ck_angular.f90
  src/module_energy_cs.f90
#  src/energy_fmt.f90
  src/energy_hydro.f90
  src/module_energy_ideal_and_external.f90
  src/energy_minimization.f90
  src/energy_nn_cs_plus_nbar.f90
  src/energy_polarization_multi.f90
  src/energy_polarization_multi_with_nccoupling.f90
  src/energy_threebody_faster.f90
  src/get_final_polarization.f90
  src/main.f90
  src/module_constants.f90
  src/module_dcf.f90
  src/module_density.f90
  src/module_vext.f90
  src/module_fastpoissonsolver.f90
  src/module_fft.f90
  src/module_geometry.f90
  src/module_grid.f90
  src/module_hardspheres.f90
  src/module_init_simu.f90
  src/module_input.f90
  src/module_mathematica.f90
  src/module_minimizer.f90
  src/module_periodic_table.f90
  src/module_postprocessing.f90
  src/module_precision_kinds.f90
  src/module_solute.f90
  src/module_solvent.f90
  src/module_system.f90
  src/module_time.f90
  src/mylbfgsb.f
  src/output_g-r-theta.f90
  src/output_gsitesite.f90
  src/output_rdf.f90
  src/write_to_cube_file.f90
)

if (UNIX)
  message(STATUS "Building on a UNIX-like platform")
  add_executable(mdft-dev ${mdft_SRCS})
  set_target_properties(mdft-dev PROPERTIES COMPILE_FLAGS "-I${FFTW_INCLUDES}")
  target_link_libraries(mdft-dev ${FFTW_LIBRARIES} m)
elseif(WIN32)
  message(STATUS "WIN32 platform not yet supported...")
endif()


#####################################
# Tell how to install this executable
#####################################

IF(WIN32)
    SET(CMAKE_INSTALL_PREFIX "C:\\Program Files")
ELSE()
    SET(CMAKE_INSTALL_PREFIX /usr/local)
ENDIF(WIN32)
INSTALL(TARGETS ${FOOEXE} RUNTIME DESTINATION bin)

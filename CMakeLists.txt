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

#option to build FFTW in same time as MDFT
option(BUILD_FFTW "build fftw when compiling MDFT" OFF)

## guard against in-source builds
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  message(FATAL_ERROR "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there. You may need to remove CMakeCache.txt. Follow instructions in README.md.")
endif()

# Build type :
# - RELEASE (with optimization flags, without openmp)
# - DEBUG   (with debug flags, without openmp)
# make sure that the default is RELEASE
if (NOT CMAKE_BUILD_TYPE)
    set (CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are: None Debug Release."
      FORCE)
endif (NOT CMAKE_BUILD_TYPE)

# allow cmake to use custom modules located in cmake subdir
# (e.g. FindFFTW)
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
#find_package(MPI REQUIRED)

#
# check compiler (gfortran, ifort, ...)
# set FFLAGS depend on the compiler
# set where fortran module files will be located
#
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
#SET (CMAKE_Fortran_COMPILER  gfortran)
message(STATUS "Using compiler ${Fortran_COMPILER_NAME}")

# Have the .mod files placed in the lib folder
SET(CMAKE_Fortran_MODULE_DIRECTORY "${PROJECT_BINARY_DIR}/mod")

if (Fortran_COMPILER_NAME MATCHES "gfortran.*")
  set (CMAKE_Fortran_FLAGS_RELEASE "-mcmodel=large -march=native -ffree-line-length-none -O3 -ffast-math -pedantic -std=f2008 -fopenmp")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-mcmodel=large -march=native -ffree-line-length-none -O3             -pedantic -std=f2008 -ffpe-trap=zero,underflow,overflow -Wall -fcheck=all -g")
  #-ftree-vectorizer-verbose=10 -fopt-info
  #set (CMAKE_Fortran_FLAGS_DEBUG   "-g -fbacktrace -pedantic -fwhole-file -Wline-truncation -Wcharacter-truncation -Wsurprising -Waliasing -fbounds-check -pg -frecursive -fcheck=all -Wall -ffpe-trap=zero,underflow,overflow")
elseif (Fortran_COMPILER_NAME MATCHES "ifort.*")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O3 -xhost")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
elseif (Fortran_COMPILER_NAME MATCHES "pgf90")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O3")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
endif()


#add_subdirectory(vendor/bpp bpp EXCLUDE_FROM_ALL)
#
#bpp_preprocess(SRC_Bpp src/module_energy_cproj_mrso.F90.bpp src/module_rotation.F90.bpp)



#############################
# Now build mdft executable #
#############################

# set source files list
set(mdft_SRCS
  src/module_blas_etc.f90
  src/module_lecture.f90
#  ${SRC_Bpp}
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
  src/module_energy_ck_angular.f90
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
  src/module_periodic_table.f90
  src/module_postprocessing.f90
  src/module_precision_kinds.f90
  src/module_solute.f90
  src/module_solvent.f90
  src/module_system.f90
  src/module_time.f90
  src/module_lbfgs_nocedal_mdft.f90
  src/output_g-r-theta.f90
  src/output_gsitesite.f90
  src/output_rdf.f90
  src/write_to_cube_file.f90
  src/module_wigner_d.f90
  src/module_orientation_projection_transform.f90
  src/module_lennardjones.f90
#  src/module_mpi.F90
)


#
# Create a symlink to the input folder
#
EXECUTE_PROCESS(COMMAND "${CMAKE_COMMAND}" -E create_symlink "${CMAKE_SOURCE_DIR}/input" "${CMAKE_BINARY_DIR}/input")

#
# Build FFTW
#
include(ExternalProject)
if(BUILD_FFTW)
	ExternalProject_Add(
	  fftw
	  URL ${PROJECT_SOURCE_DIR}/src/fftw-3.3.4.tar.gz
	  PREFIX ${CMAKE_CURRENT_BINARY_DIR}/fftw3-cmake
	  CONFIGURE_COMMAND <SOURCE_DIR>/configure --prefix=${CMAKE_CURRENT_BINARY_DIR}/libs/fftw3  --enable-openmp --enable-threads --enable-sse2 --enable-avx
	  BUILD_COMMAND make -j4
	  INSTALL_COMMAND make install -j4
	)
	ExternalProject_Add(
	  fftwf
	  URL ${PROJECT_SOURCE_DIR}/src/fftw-3.3.4.tar.gz
	  PREFIX ${CMAKE_CURRENT_BINARY_DIR}/fftw3-cmake
	  CONFIGURE_COMMAND <SOURCE_DIR>/configure --prefix=${CMAKE_CURRENT_BINARY_DIR}/libs/fftw3 --enable-float --enable-openmp --enable-threads --enable-sse --enable-sse2 --enable-avx
	  BUILD_COMMAND make -j4
	  INSTALL_COMMAND make install -j4
	)
	SET(FFTW_INCLUDES ${CMAKE_CURRENT_BINARY_DIR}/libs/fftw3/include)
	SET(FFTW_LIBRARIES ${CMAKE_CURRENT_BINARY_DIR}/libs/fftw3/lib/libfftw3.a;${CMAKE_CURRENT_BINARY_DIR}/libs/fftw3/lib/libfftw3f.a)
else()
	find_package(FFTW REQUIRED)
endif()



if (UNIX)
  message(STATUS "Building on a UNIX-like platform")


  include_directories(${FFTW_INCLUDES})
  add_executable(mdft-dev ${mdft_SRCS})

  message(STATUS "Compiling with fftw ${FFTW_LIBRARIES}")
  target_link_libraries(mdft-dev ${FFTW_LIBRARIES})

  if(BUILD_FFTW) 
    add_dependencies(mdft-dev fftw fftwf)
  endif()

#  if("${MPI_Fortran_HAS_MODULE}")
#        set_property(TARGET  mdft-dev
#                APPEND PROPERTY
#                COMPILE_DEFINITIONS "MPI_FORTRAN_HAS_MODULE")
#  endif()

#  include_directories(${MPI_Fortran_INCLUDE_PATH} ${MPI_C_INCLUDE_PATH})
#  set_property(TARGET  mdft-dev APPEND_STRING PROPERTY COMPILE_FLAGS
#        " ${MPI_Fortran_COMPILE_FLAGS}  ${MPI_C_COMPILE_FLAGS} ")
#  set_property(TARGET  mdft-dev APPEND_STRING PROPERTY LINK_FLAGS
#        " ${MPI_Fortran_LINK_FLAGS} ${MPI_C_LINK_FLAGS} ")
#  target_link_libraries(mdft-dev ${MPI_Fortran_LIBRARIES} ${MPI_C_LIBRARIES})

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

################################################################################
#
# \file      tests/regression/CMakeLists.txt
# \copyright 2012-2015 J. Bakosi,
#            2016-2018 Los Alamos National Security, LLC.,
#            2019-2021 Triad National Security, LLC.
#            2022-2023 J. Bakosi
#            All rights reserved. See the LICENSE file for details.
# \brief     Cmake code common to all regression testing
#
################################################################################

# Set paths required for running tests
set(TEST_RUNNER ${CMAKE_SOURCE_DIR}/../cmake/test_runner.cmake)

# Find gmsh executable
find_package(Gmsh QUIET)
if(NOT GMSH_FOUND)
  message(STATUS
         "Gmsh not found, will not configure meshconv tests with gmsh output")
endif()

# Find valgrind executable
set(ENABLE_VALGRIND false CACHE BOOL "Run all regression tests using valgrind")
if (ENABLE_VALGRIND)
  find_program(VALGRIND valgrind)
  if (VALGRIND)
    set(USE_VALGRIND true)
    message(STATUS "Will run all regression tests using valgrind")
  else()
    message(WARNING "Valgrind was enabled but NOT found")
  endif()
endif()

# Include function used to add regression tests
include(add_regression_test)

# Include regression test directories

# UnitTest regression tests
message(STATUS "Adding regression tests for ${UNITTEST_EXECUTABLE}")
add_subdirectory(unittest/exec)

# MeshConv regression tests
message(STATUS "Adding regression tests for ${MESHCONV_EXECUTABLE}")
if(GMSH_FOUND)
  add_subdirectory(meshconv/gmsh_output)
endif()
add_subdirectory(meshconv/netgen_output)
add_subdirectory(meshconv/exo_output)
add_subdirectory(meshconv/exec)

# Inciter regression tests
message(STATUS "Adding regression tests for ${INCITER_EXECUTABLE}")
add_subdirectory(inciter/exec)
add_subdirectory(inciter/asynclogic)
add_subdirectory(inciter/restart)
add_subdirectory(inciter/compflow/Euler/VorticalFlow)
add_subdirectory(inciter/compflow/Euler/NonlinearEnergyGrowth)
add_subdirectory(inciter/compflow/Euler/TaylorGreen)
add_subdirectory(inciter/compflow/Euler/RayleighTaylor)
add_subdirectory(inciter/compflow/Euler/Sod)
add_subdirectory(inciter/compflow/Euler/Sedov)
add_subdirectory(inciter/compflow/Euler/Pipe)
add_subdirectory(inciter/transport/slot_cyl)
add_subdirectory(inciter/mesh_refinement/t0ref)
#add_subdirectory(inciter/mesh_refinement/dtref)

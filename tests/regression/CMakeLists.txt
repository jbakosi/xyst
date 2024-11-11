################################################################################
#
# \file      tests/regression/CMakeLists.txt
# \copyright 2012-2015 J. Bakosi,
#            2016-2018 Los Alamos National Security, LLC.,
#            2019-2021 Triad National Security, LLC.,
#            2022-2024 J. Bakosi
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
    set(sup "${PROJECT_SOURCE_DIR}/../tests/regression/valgrind.sup")
    #set(VALGRIND_CMD "${VALGRIND} --tool=memcheck --suppressions=${sup}")
    set(VALGRIND_CMD "${VALGRIND} --tool=memcheck --suppressions=${sup} --gen-suppressions=all --log-file=${sup}.log")
    message(STATUS "Will run all regression tests prefixed by '${VALGRIND_CMD}'")
  else()
    message(WARNING "Valgrind was enabled but NOT found")
  endif()
endif()

add_subdirectory(numdiff)

# configure helper executable to test signals from executables
add_executable(signal_wrapper signal_wrapper.cpp)
set(SIGNAL_WRAPPER "${CMAKE_CURRENT_BINARY_DIR}/signal_wrapper")
add_dependencies(${MESHCONV_EXECUTABLE} signal_wrapper)

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
add_subdirectory(inciter/asynclogic)
add_subdirectory(inciter/exec)
add_subdirectory(inciter/href/t0ref)
#add_subdirectory(inciter/href/dtref)
add_subdirectory(inciter/restart)
add_subdirectory(inciter/RieCG/VorticalFlow)
add_subdirectory(inciter/RieCG/NonlinearEnergyGrowth)
add_subdirectory(inciter/RieCG/TaylorGreen)
add_subdirectory(inciter/RieCG/RayleighTaylor)
add_subdirectory(inciter/RieCG/Sod)
add_subdirectory(inciter/RieCG/Sedov)
add_subdirectory(inciter/RieCG/Pipe)
add_subdirectory(inciter/RieCG/SlotCyl)
add_subdirectory(inciter/RieCG/Canyon)
add_subdirectory(inciter/KozCG/VorticalFlow)
add_subdirectory(inciter/KozCG/TaylorGreen)
add_subdirectory(inciter/KozCG/RayleighTaylor)
add_subdirectory(inciter/KozCG/NonlinearEnergyGrowth)
add_subdirectory(inciter/KozCG/Sod)
add_subdirectory(inciter/KozCG/SlotCyl)
add_subdirectory(inciter/ZalCG/Sod)
add_subdirectory(inciter/ZalCG/Sedov)
add_subdirectory(inciter/ZalCG/SlotCyl)
add_subdirectory(inciter/ZalCG/Bump)
add_subdirectory(inciter/LaxCG/Bump)
add_subdirectory(inciter/ChoCG/Poisson)
add_subdirectory(inciter/ChoCG/Sphere)
add_subdirectory(inciter/ChoCG/Poiseuille)
add_subdirectory(inciter/ChoCG/Lid)
add_subdirectory(inciter/LohCG/Poiseuille)
add_subdirectory(inciter/LohCG/Lid)

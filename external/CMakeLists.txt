# Build Xyst external libraries

cmake_minimum_required(VERSION 3.1.0 FATAL_ERROR)

project(xyst-external C CXX)
include(ExternalProject)

# Set cmake modules directory
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../cmake")

# Disallow in-source builds
include(DisallowInSourceBuilds)
# Set CMAKE_BUILD_TYPE
include(BuildType)
# Find MPI compilers
include(MPICompilers)

MESSAGE(STATUS "Underlying C++ compiler: " ${UNDERLYING_CXX_COMPILER})

set(BUILD_SHARED_LIBS ON
    CACHE BOOL "Build shared libraries. Possible values: ON | OFF")

# Set compiler id
string(TOLOWER "${CMAKE_CXX_COMPILER_ID}" COMP)

option(SMP "Build Charm++ in SMP mode" OFF)
if(SMP)
  set(smp "-smp")
endif()

# Set install directory based on compiler ID if was not given on command line
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX "${PROJECT_SOURCE_DIR}/install/${COMP}${smp}${rndq}"
      CACHE STRING "Third-party libraries install prefix." FORCE)
  message(STATUS "CMAKE_INSTALL_PREFIX not specified, setting to ${CMAKE_INSTALL_PREFIX}")
else()
  message(STATUS "CMAKE_INSTALL_PREFIX: " ${CMAKE_INSTALL_PREFIX})
endif()

ExternalProject_Add(
  charm
  PREFIX charm
  GIT_REPOSITORY https://github.com/UIUC-PPL/charm
  GIT_TAG v7.0.0-rc2
  PATCH_COMMAND patch -p1 -t -i ${PROJECT_SOURCE_DIR}/../cmake/charm-fortran.patch
  CMAKE_ARGS -D TARGET=charm
             -D NETWORK=mpi
             -D SMP=${SMP}
             -D ENABLE_FORTRAN=off
             -D CMAKE_C_COMPILER:PATH=${MPI_C_COMPILER}
             -D CMAKE_CXX_COMPILER:PATH=${MPI_CXX_COMPILER}
             -D CMAKE_INSTALL_PREFIX:PATH=${CMAKE_INSTALL_PREFIX}
  LOG_DOWNLOAD 1
  LOG_CONFIGURE 1
  LOG_BUILD 1
  LOG_INSTALL 1
)

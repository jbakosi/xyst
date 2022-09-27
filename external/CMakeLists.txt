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

option(SMP "Build Charm++ in SMP mode")
if(SMP)
  set(smp "-smp")
endif()

# Set install directory based on compiler ID if was not given on command line
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX "${PROJECT_SOURCE_DIR}/install/${COMP}${smp}"
      CACHE STRING "Third-party libraries install prefix." FORCE)
  message(STATUS "CMAKE_INSTALL_PREFIX not specified, setting to ${CMAKE_INSTALL_PREFIX}")
else()
  message(STATUS "CMAKE_INSTALL_PREFIX: " ${CMAKE_INSTALL_PREFIX})
endif()

set(HDF5_PREFER_PARALLEL true)
if(NOT BUILD_SHARED_LIBS)
  set(HDF5_USE_STATIC_LIBRARIES true)
endif()
find_package(HDF5 COMPONENTS C HL REQUIRED)
find_package(NetCDF REQUIRED)

ExternalProject_Add(
  trilinos
  PREFIX trilinos
  GIT_REPOSITORY https://github.com/trilinos/Trilinos.git
  GIT_TAG trilinos-release-12-12-1
  PATCH_COMMAND patch -p1 -t -i ${PROJECT_SOURCE_DIR}/../cmake/exodus_gcc10.patch
  CMAKE_ARGS -D Trilinos_ENABLE_ALL_PACKAGES:BOOL=OFF
             -D CMAKE_CXX_STANDARD:STRING=17
             -D TPL_ENABLE_MPI:BOOL=ON
             -D BUILD_SHARED_LIBS:BOOL=ON
             -D Trilinos_ENABLE_ALL_OPTIONAL_PACKAGES=OFF
             -D Trilinos_ENABLE_Zoltan2:BOOL=ON
             -D Trilinos_ENABLE_SEACASExodus:BOOL=ON
             -D Trilinos_ENABLE_SEACASExodiff:BOOL=ON
             -D Trilinos_ENABLE_Fortran=OFF
             -D Trilinos_ENABLE_DEBUG=OFF
             -D TPL_ENABLE_Netcdf:BOOL=ON
             -D Netcdf_LIBRARY_DIRS:PATH=${NETCDF_LIBRARY_DIRS}
             -D TPL_Netcdf_INCLUDE_DIRS:PATH=${NETCDF_INCLUDES}
             -D CMAKE_INSTALL_RPATH:PATH=${CMAKE_INSTALL_PREFIX}/lib
             -D CMAKE_C_COMPILER:PATH=${MPI_C_COMPILER}
             -D CMAKE_CXX_COMPILER:PATH=${MPI_CXX_COMPILER}
             -D CMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
             -D CMAKE_INSTALL_PREFIX:PATH=${CMAKE_INSTALL_PREFIX}
  LOG_DOWNLOAD 1
  LOG_PATCH 1
  LOG_CONFIGURE 1
  LOG_BUILD 1
  LOG_INSTALL 1
)

ExternalProject_Add(
  charm
  PREFIX charm
  GIT_REPOSITORY https://github.com/UIUC-PPL/charm
  GIT_TAG v7.0.0-rc2
  PATCH_COMMAND patch -p1 -t -i ${PROJECT_SOURCE_DIR}/../cmake/charm-fortran.patch
  CMAKE_ARGS -D TARGET=LIBS
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

ExternalProject_Add(
  tut
  PREFIX tut
  GIT_REPOSITORY https://github.com/mrzechonek/tut-framework
  CONFIGURE_COMMAND ""
  BUILD_COMMAND ""
  INSTALL_COMMAND ${CMAKE_COMMAND} -E copy_directory
                  ${PROJECT_BINARY_DIR}/tut/src/tut/include/tut
                  ${CMAKE_INSTALL_PREFIX}/include/tut
  LOG_DOWNLOAD 1
  LOG_INSTALL 1
)

ExternalProject_Add(
  numdiff
  PREFIX numdiff
  GIT_REPOSITORY https://codeberg.org/xyst/ndiff
  CMAKE_ARGS -D CMAKE_C_COMPILER:PATH=${CMAKE_C_COMPILER}
             -D CMAKE_CXX_COMPILER:PATH=${CMAKE_CXX_COMPILER}
             -D CMAKE_INSTALL_PREFIX:PATH=${CMAKE_INSTALL_PREFIX}
  LOG_DOWNLOAD 1
  LOG_CONFIGURE 1
  LOG_BUILD 1
  LOG_INSTALL 1
)

ExternalProject_Add(
  highwayhash
  PREFIX highwayhash
  GIT_REPOSITORY https://github.com/google/highwayhash
  CONFIGURE_COMMAND ""
  BUILD_COMMAND ""
  INSTALL_COMMAND ${CMAKE_COMMAND} -E copy_directory
                  ${PROJECT_BINARY_DIR}/highwayhash/src/highwayhash/highwayhash
                  ${CMAKE_INSTALL_PREFIX}/include/highwayhash
  LOG_DOWNLOAD 1
  LOG_INSTALL 1
)

ExternalProject_Add(
  brigand 
  PREFIX brigand
  GIT_REPOSITORY https://github.com/edouarda/brigand
  CONFIGURE_COMMAND ""
  BUILD_COMMAND ""
  INSTALL_COMMAND ${CMAKE_COMMAND} -E copy_directory
                  ${PROJECT_BINARY_DIR}/brigand/src/brigand/include/brigand
                  ${CMAKE_INSTALL_PREFIX}/include/brigand
  LOG_DOWNLOAD 1
  LOG_INSTALL 1
)

ExternalProject_Add(
  testu01 
  PREFIX testu01
  GIT_REPOSITORY https://codeberg.org/xyst/TestU01-threadsafe
  CMAKE_ARGS -D CMAKE_C_COMPILER:PATH=${CMAKE_C_COMPILER}
             -D CMAKE_CXX_COMPILER:PATH=${CMAKE_CXX_COMPILER}
             -D CMAKE_INSTALL_PREFIX:PATH=${CMAKE_INSTALL_PREFIX}/testu01
  LOG_DOWNLOAD 1
  LOG_CONFIGURE 1
  LOG_BUILD 1
  LOG_INSTALL 1
)

ExternalProject_Add(
  mcss
  PREFIX mcss
  GIT_REPOSITORY https://github.com/mosra/m.css
  GIT_TAG 91ff035
  CONFIGURE_COMMAND ""
  BUILD_COMMAND ""
  INSTALL_COMMAND ${CMAKE_COMMAND} -E copy_directory
                  ${PROJECT_BINARY_DIR}/mcss/src/mcss
                  ${CMAKE_INSTALL_PREFIX}/m.css
  LOG_DOWNLOAD 1
  LOG_INSTALL 1
)

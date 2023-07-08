################################################################################
#
# \file      cmake/FindHDF5AMPI.cmake
# \copyright 2012-2015 J. Bakosi,
#            2016-2018 Los Alamos National Security, LLC.,
#            2019-2021 Triad National Security, LLC.
#            2022-2023 J. Bakosi
#            All rights reserved. See the LICENSE file for details.
# \brief     Find HDF5
#
################################################################################

# Optionally, set HDF5_INSTALL_DIR before calling find_package to a custom path.

if(HDF5_INCLUDE_DIRS AND HDF5_C_LIBRARIES AND HDF5_HL_LIBRARIES)
  set(HDF5AMPI_FIND_QUIETLY TRUE)
endif()

find_library(HDF5_C_LIBRARY NAMES libhdf5.a
             PATHS ${HDF5_INSTALL_DIR}
                   ${CMAKE_BINARY_DIR}/hdf5/install
             PATH_SUFFIXES lib)

set(HDF5_C_LIBRARIES ${HDF5_C_LIBRARY})

find_library(HDF5_HL_LIBRARY NAMES libhdf5_hl.a
             PATHS ${HDF5_INSTALL_DIR}
                   ${CMAKE_BINARY_DIR}/hdf5/install
             PATH_SUFFIXES lib)

set(HDF5_HL_LIBRARIES ${HDF5_HL_LIBRARY})

find_path(HDF5_INCLUDE_DIR NAMES hdf5.h
                           PATHS ${HDF5_INSTALL_DIR}
                                 ${CMAKE_BINARY_DIR}/hdf5/install
                           PATH_SUFFIXES include)

set(HDF5_INCLUDE_DIRS ${HDF5_INCLUDE_DIR})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(HDF5AMPI DEFAULT_MSG HDF5_INCLUDE_DIRS HDF5_C_LIBRARIES HDF5_HL_LIBRARIES)

mark_as_advanced(HDF5_INCLUDE_DIRS HDF5_C_LIBRARIES HDF5_HL_LIBRARIES)

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

if(HDF5_INCLUDE_DIR AND HDF5_INCLUDE_DIRS AND
   HDF5_LIBRARIES AND HDF5_C_LIBRARY AND HDF5_HL_LIBRARY)
  set(HDF5AMPI_FIND_QUIETLY TRUE)
endif()

find_library(HDF5_C_LIBRARY NAMES libhdf5.a
             PATHS ${HDF5_INSTALL_DIR}
                   ${CMAKE_BINARY_DIR}/hdf5/install
             PATH_SUFFIXES lib)

find_library(HDF5_HL_LIBRARY NAMES libhdf5_hl.a
             PATHS ${HDF5_INSTALL_DIR}
                   ${CMAKE_BINARY_DIR}/hdf5/install
             PATH_SUFFIXES lib)

set(HDF5_LIBRARIES ${HDF5_HL_LIBRARY} ${HDF5_C_LIBRARY})

find_path(HDF5_INCLUDE_DIR NAMES hdf5.h
                           PATHS ${HDF5_INSTALL_DIR}
                                 ${CMAKE_BINARY_DIR}/hdf5/install
                           PATH_SUFFIXES include)

set(HDF5_INCLUDE_DIRS ${HDF5_INCLUDE_DIR})

include (FindPackageHandleStandardArgs)
find_package_handle_standard_args(HDF5AMPI DEFAULT_MSG HDF5_INCLUDE_DIR HDF5_INCLUDE_DIRS HDF5_LIBRARIES HDF5_C_LIBRARY HDF5_HL_LIBRARY)

mark_as_advanced(HDF5_INCLUDE_DIR HDF5_INCLUDE_DIRS HDF5_LIBRARIES HDF5_C_LIBRARY HDF5_HL_LIBRARY)

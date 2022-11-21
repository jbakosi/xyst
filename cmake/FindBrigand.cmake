################################################################################
#
# \file      FindBrigand.cmake
# \copyright 2012-2015 J. Bakosi,
#            2016-2018 Los Alamos National Security, LLC.,
#            2019-2021 Triad National Security, LLC.
#            2022-2023 J. Bakosi
#            All rights reserved. See the LICENSE file for details.
# \brief     Find the Brigand template metaprogramming library
#
################################################################################

# Brigand: https://github.com/edouarda/brigand
#
#  BRIGAND_FOUND - System has Brigand
#  BRIGAND_INCLUDE_DIRS - The Brigand include directory

# If already in cache, be silent
if(BRIGAND_INCLUDE_DIRS)
  set (BRIGAND_FIND_QUIETLY TRUE)
endif()

find_path(BRIGAND_INCLUDE_DIR NAMES brigand.hpp PATH_SUFFIXES brigand)

set(BRIGAND_INCLUDE_DIRS ${BRIGAND_INCLUDE_DIR})

# Handle the QUIETLY and REQUIRED arguments and set Brigand_FOUND to TRUE if
# all listed variables are TRUE.
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(Brigand REQUIRED_VARS BRIGAND_INCLUDE_DIRS)

MARK_AS_ADVANCED(BRIGAND_INCLUDE_DIRS)

################################################################################
#
# \file      FindGmsh.cmake
# \copyright 2012-2015 J. Bakosi,
#            2016-2018 Los Alamos National Security, LLC.,
#            2019-2021 Triad National Security, LLC.,
#            2022-2025 J. Bakosi
#            All rights reserved. See the LICENSE file for details.
# \brief     Find Gmsh
#
################################################################################

# Gmsh: http://gmsh.info
#
# GMSH_FOUND - System has gmsh
# GMSH_EXECUTABLE - The gmsh executable
#
# Usage:
#
# set(GMSH_ROOT "/path/to/custom/gmsh") # prefer over system
# find_package(Gmsh)

if(GMSH_EXECUTABLE)
  # Already in cache, be silent
  set (GMSH_FIND_QUIETLY TRUE)
endif()

find_program(GMSH_EXECUTABLE NAMES gmsh)

# Handle the QUIETLY and REQUIRED arguments and set GMSH_FOUND to TRUE if
# all listed variables are TRUE.
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(Gmsh DEFAULT_MSG GMSH_EXECUTABLE)

MARK_AS_ADVANCED(GMSH_EXECUTABLE)

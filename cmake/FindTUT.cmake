################################################################################
#
# \file      FindTUT.cmake
# \copyright 2012-2015 J. Bakosi,
#            2016-2018 Los Alamos National Security, LLC.,
#            2019-2021 Triad National Security, LLC.
#            All rights reserved. See the LICENSE file for details.
# \brief     Find the Template Unit Test library headers
#
################################################################################

# TUT: http://mrzechonek.github.io/tut-framework/
#
#  TUT_FOUND        - True if TUT is found
#  TUT_INCLUDE_DIRS - TUT include files directories

# Look for the header file
FIND_PATH(TUT_INCLUDE_DIR NAMES tut.hpp PREFIX tut)

set(TUT_INCLUDE_DIRS ${TUT_INCLUDE_DIR})

# Handle the QUIETLY and REQUIRED arguments and set TUT_FOUND to TRUE if 
# all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(TUT DEFAULT_MSG TUT_INCLUDE_DIRS)

MARK_AS_ADVANCED(TUT_INCLUDE_DIRS)

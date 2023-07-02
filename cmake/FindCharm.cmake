################################################################################
#
# \file      FindCharm.cmake
# \copyright 2012-2015 J. Bakosi,
#            2016-2018 Los Alamos National Security, LLC.,
#            2019-2021 Triad National Security, LLC.
#            2022-2023 J. Bakosi
#            All rights reserved. See the LICENSE file for details.
# \brief     Find Charm++
#
################################################################################

# Charm++: http://charmplusplus.org
#
#  CHARM_FOUND        - True if the charmc compiler wrapper was found
#  CHARM_INCLUDE_DIRS - Charm++ include files paths
#  CHARM_COMPILER     - Charmc compiler wrapper
#  CHARM_RUN          - Charmrun executable runner
#
#  Set CHARM_INSTALL_DIR before calling find_package to a path to add an
#  additional search path.

function(_GET_CHARMINC _OUT_INC _charmc)
  file(STRINGS ${_charmc} _contents REGEX "^CHARMINC=")
  if(_contents)
    string(REGEX REPLACE "^CHARMINC=\"(.*)\"" "\\1" ${_OUT_INC} "${_contents}")
    set(${_OUT_INC} ${${_OUT_INC}} PARENT_SCOPE)
  else()
    message(FATAL_ERROR "file ${_charmc} does not exist")
  endif()
endfunction()

# If already in cache, be silent
if (CHARM_INCLUDE_DIRS AND CHARM_COMPILER AND CHARM_RUN)
  set (CHARM_FIND_QUIETLY TRUE)
endif()

FIND_PROGRAM(CHARM_COMPILER
  NAMES charmc
  PATHS ${CHARM_INSTALL_DIR}
        ${CMAKE_BINARY_DIR}/charm/install
  PATH_SUFFIXES bin
)

FIND_PROGRAM(CHARM_RUN
  NAMES charmrun
  PATHS ${CHARM_INSTALL_DIR}
        ${CMAKE_BINARY_DIR}/charm/install
  PATH_SUFFIXES bin
)

if(CHARM_COMPILER)
  _GET_CHARMINC(HINTS_CHARMINC ${CHARM_COMPILER})
endif()

FIND_PATH(CHARM_INCLUDE_DIR NAMES charm.h
                            HINTS ${HINTS_CHARMINC}
                                  ${CHARM_INSTALL_DIR}/include
                                  ${CMAKE_BINARY_DIR}/charm/install/include
                            PATH_SUFFIXES charm)

if(CHARM_INCLUDE_DIR)
  set(CHARM_INCLUDE_DIRS ${CHARM_INCLUDE_DIR})
else()
  set(CHARM_INCLUDE_DIRS "")
endif()

# Handle the QUIETLY and REQUIRED arguments and set CHARM_FOUND to TRUE if all
# listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(Charm DEFAULT_MSG CHARM_COMPILER
                                  CHARM_INCLUDE_DIRS CHARM_RUN)

if(CHARM_COMPILER)
  set(CMAKE_REQUIRED_QUIET 1)
  include(CheckIncludeFiles)
  CHECK_INCLUDE_FILES("${CHARM_INCLUDE_DIR}/conv-autoconfig.h"
                      HAVE_CHARM_CONV_MACH_OPT)
  if (HAVE_CHARM_CONV_MACH_OPT)
    include(CheckSymbolExists)
    CHECK_SYMBOL_EXISTS(CMK_SMP "${CHARM_INCLUDE_DIR}/conv-autoconfig.h" SMP)
    if (SMP)
      message(STATUS "Charm++ built in SMP mode")
    else()
      message(STATUS "Charm++ built in non-SMP mode")
    endif()
  endif()
endif()

MARK_AS_ADVANCED(CHARM_COMPILER CHARM_INCLUDE_DIRS CHARM_RUN SMP)

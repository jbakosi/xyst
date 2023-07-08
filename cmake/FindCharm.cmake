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

# Optionally, set CHARM_INSTALL_DIR before calling find_package to a custom path.

function(_GET_CHARMINC _OUT_INC _charmc)
  file(STRINGS ${_charmc} _contents REGEX "^CHARMINC=")
  if(_contents)
    string(REGEX REPLACE "^CHARMINC=\"(.*)\"" "\\1" ${_OUT_INC} "${_contents}")
    set(${_OUT_INC} ${${_OUT_INC}} PARENT_SCOPE)
  else()
    message(FATAL_ERROR "file ${_charmc} does not exist")
  endif()
endfunction()

# find out if Charm++ was built in SMP mode
function(GET_CHARM_DEF RES SEARCH_DEF)
  file(STRINGS ${CHARM_INCLUDE_DIR}/conv-autoconfig.h _contents
       REGEX ".*#define ${SEARCH_DEF}[ \t]+")
  if(_contents)
    string(REGEX REPLACE ".*#define ${SEARCH_DEF}[ \t]+(on|ON|off|OFF|0|1|true|TRUE|false|FALSE).*" "\\1" ${RES} "${_contents}")
    string(TOLOWER ${${RES}} ${RES})
    if (${RES} STREQUAL "1" OR
        ${RES} STREQUAL "on" OR
        ${RES} STREQUAL "true")
      set(${RES} true PARENT_SCOPE)
    else()
      set(${RES} false PARENT_SCOPE)
    endif()
 endif()
endfunction()

find_program(CHARM_COMPILER
  NAMES charmc
  PATHS ${CHARM_INSTALL_DIR}
        ${CMAKE_BINARY_DIR}/charm/install
  PATH_SUFFIXES bin
)

find_program(CHARM_RUN
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
                                  ${CHARM_INSTALL_DIR}
                                  ${CMAKE_BINARY_DIR}/charm/install
                            PATH_SUFFIXES include)

if(CHARM_INCLUDE_DIR)
  set(CHARM_INCLUDE_DIRS ${CHARM_INCLUDE_DIR})
else()
  set(CHARM_INCLUDE_DIRS "")
endif()

if(CHARM_COMPILER)
  GET_CHARM_DEF(SMP CMK_SMP)
  if (SMP)
    message(STATUS "Charm++ built in SMP mode")
  else()
    message(STATUS "Charm++ built in non-SMP mode")
  endif()

  GET_CHARM_DEF(RNDQ CMK_RANDOMIZED_MSGQ)
  if (RNDQ)
    message(STATUS "Charm++ built with randomized message queues")
  endif()
endif()

find_program(AMPI_C_COMPILER
  NAMES ampicc
  PATHS ${CHARM_INSTALL_DIR}
        ${CMAKE_BINARY_DIR}/charm/install
  PATH_SUFFIXES bin
)

find_program(AMPI_CXX_COMPILER
  NAMES ampicxx
  PATHS ${CHARM_INSTALL_DIR}
        ${CMAKE_BINARY_DIR}/charm/install
  PATH_SUFFIXES bin
)

find_program(AMPI_RUN
  NAMES ampirun
  PATHS ${CHARM_INSTALL_DIR}
        ${CMAKE_BINARY_DIR}/charm/install
  PATH_SUFFIXES bin
)

if(AMPI_C_COMPILER AND AMPI_CXX_COMPILER AND AMPI_RUN)
  set(AMPI_FOUND true)
  message(STATUS "Charm++ built with AMPI")
  FIND_PATH(MPI_INCLUDE_DIR
            NAMES mpi.h
            PATHS ${HINTS_CHARMINC}
                  ${CHARM_INSTALL_DIR}/include/ampi
                  ${CMAKE_BINARY_DIR}/charm/install/include/ampi)
  set(MPI_C_INCLUDE_DIRS ${MPI_INCLUDE_DIR} ${CHARM_INCLUDE_DIR})
  set(MPI_CXX_INCLUDE_DIRS ${MPI_INCLUDE_DIR} ${CHARM_INCLUDE_DIR})
endif()

# Handle the QUIETLY and REQUIRED arguments and set CHARM_FOUND to TRUE if all
# listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(Charm DEFAULT_MSG CHARM_COMPILER CHARM_INCLUDE_DIRS CHARM_RUN)

MARK_AS_ADVANCED(CHARM_COMPILER CHARM_INCLUDE_DIRS CHARM_RUN SMP AMPI_FOUND AMPI_C_COMPILER AMPI_CXX_COMPILER AMPI_RUN MPI_C_INCLUDE_DIRS MPI_CXX_INCLUDE_DIRS)

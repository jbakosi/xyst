################################################################################
#
# \file      CodeCoverage.cmake
# \copyright 2012-2015 J. Bakosi,
#            2016-2018 Los Alamos National Security, LLC.,
#            2019-2021 Triad National Security, LLC.
#            2022-2023 J. Bakosi
#            All rights reserved. See the LICENSE file for details.
# \brief     Setup target for code coverage analysis
#
################################################################################

# ##############################################################################
# Function to add code coverage target for all individual coverage targets
# included. This is very similar to setup_target_for_coverage(), defined above,
# but compiles a coverage report for all the individual test suites. The
# TESTRUNNER_ARGS are hard-coded.
#
# setup_target_for_all_coverage( <path> <targetname> [DEPENDS dep1 dep2 ... ] )
#
# Mandatory arguments:
# --------------------
#
# path - Path to prepend to where the report is generated:
# <path>${targetname}/index.html.
#
# targetname - The name of the code coverage target. The HTML report on code
# coverage is generated at the path <path>/${targetname}/index.html.
#
# unittestrunner - Command line of the unittest harness runner.
#
# unittestrunner_ncpus_arg - Command line argument specifying the number of cpus
# argument of the unittest harness runner.
#
# Optional arguments:
# -------------------
#
# DEPENDS dep1 dep2 ... - Optional dependencies added to test coverage target.
# Default: "". Here all dependencies should be given that should be covered by
# the test suite the coverage is being setup for, as well as those that are
# required for successfully building the tests and the test runner.
#
# ##############################################################################
FUNCTION(SETUP_TARGET_FOR_ALL_COVERAGE path targetname unittestrunner
                                       unittestrunner_ncpus_arg)

  set(multiValueArgs DEPENDS)
  cmake_parse_arguments(ARG "${options}" "${oneValueArgs}" "${multiValueArgs}"
                        ${ARGN})

  IF(NOT FASTCOV)
    MESSAGE(FATAL_ERROR "fastcov not found! Aborting...")
  ENDIF()

  IF(NOT GENHTML)
    MESSAGE(FATAL_ERROR "genhtml not found! Aborting...")
  ENDIF()

  IF(NOT SED)
    MESSAGE(FATAL_ERROR "sed not found! Aborting...")
  ENDIF()

  # Set shorcut for output: path/target
  set(OUTPUT ${path}/${targetname})
  file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/${path})

  # Setup code coverage target
  ADD_CUSTOM_TARGET(${targetname}
    # Zero coverage counters
    COMMAND ${FASTCOV} --zerocounters
    # Run all test suites
    COMMAND ${unittestrunner} ${unittestrunner_ncpus_arg} ${PROCESSOR_COUNT} Main/${UNITTEST_EXECUTABLE} -v
    COMMAND ${CMAKE_CTEST_COMMAND} -j${PROCESSOR_COUNT} -LE extreme
    # Process gcov output for genhtml
    COMMAND ${FASTCOV} --branch-coverage --exceptional-branch-coverage --lcov -o ${OUTPUT}.info --exclude tests/ c++/ include/ boost/ brigand/ charm/ decl.h def.h openmpi pegtl exodus/ tut/ highwayhash/ zoltan/ moduleinit
    # Copy over report customization files for genhtml
    COMMAND ${CMAKE_COMMAND} -E copy
            ${CMAKE_SOURCE_DIR}/../doc/xyst.gcov.css
            ${CMAKE_BINARY_DIR}
    COMMAND ${CMAKE_COMMAND} -E copy
            ${CMAKE_SOURCE_DIR}/../doc/xyst.lcov.prolog
            ${CMAKE_BINARY_DIR}
    # Generate HTML report
    COMMAND ${GENHTML} --legend --rc genhtml_branch_coverage=1 --demangle-cpp --css-file xyst.gcov.css --ignore-errors source --html-prolog xyst.lcov.prolog -o ${OUTPUT} ${OUTPUT}.info
    # Customize page headers in generated html
    COMMAND find ${OUTPUT} -type f -exec ${SED} -i "s/LCOV - code coverage report/Xyst test code coverage report/g" {} +
    COMMAND find ${OUTPUT} -type f -exec ${SED} -i "s/<td class=\"headerItem\">Test:<\\/td>/<td class=\"headerItem\">Commit:<\\/td>/g" {} +
    COMMAND find ${OUTPUT} -type f -exec ${SED} -i "s/test_coverage.info/<a href=\"https:\\/\\/codeberg.org\\/xyst\\/xyst\\/commit\\/${GIT_HEAD_SHA1}\">${GIT_HEAD_SHA1}<\\/a>/g" {} +
    COMMAND find ${OUTPUT} -type f -exec ${SED} -i "s/\\/home\\/jbakosi\\/code\\/xyst\\/build\\/gnu\\/Main/build\\/Main/g" {} +
    # Cleanup intermediate data
    COMMAND ${CMAKE_COMMAND} -E remove ${OUTPUT}.info
    # Set work directory for target
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    # Echo what is being done
    COMMENT "Test code coverage report"
    VERBATIM USES_TERMINAL
  )

  # Make test coverage target dependent on optional dependencies passed in using
  # keyword DEPENDS
  add_dependencies(${targetname} ${ARG_DEPENDS})

  # Output code coverage target enabled
  message(STATUS "Enabling code coverage target '${targetname}' generating full coverage, dependencies {${ARG_DEPENDS}}, report at ${OUTPUT}/index.html")

ENDFUNCTION()

################################################################################
#
# \file      CppCheck.cmake
# \copyright 2012-2015 J. Bakosi,
#            2016-2018 Los Alamos National Security, LLC.,
#            2019-2021 Triad National Security, LLC.,
#            2022-2024 J. Bakosi
#            All rights reserved. See the LICENSE file for details.
# \brief     Setup target for code coverage analysis
#
################################################################################

find_program( CPPCHECK cppcheck )
find_program( CPPCHECK_HTMLREPORT cppcheck-htmlreport )

if(CPPCHECK AND CPPCHECK_HTMLREPORT)

  find_program( FILEFIND find )
  find_program( SED sed )

  if(FILEFIND AND SED)
    file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/doc/html/${CMAKE_BUILD_TYPE}/cppcheck)
    ADD_CUSTOM_TARGET(cppcheck-xml
      # Run cppcheck static analysis
      COMMAND ${CPPCHECK} --inline-suppr --enable=all --force
              --xml --xml-version=2 #-j${PROCESSOR_COUNT}
              -I${PROJECT_SOURCE_DIR}/Base
              -I${PROJECT_SOURCE_DIR}/Control
              -I${PROJECT_SOURCE_DIR}/NoWarning
              -I${PROJECT_BINARY_DIR}/Main
              -i${PROJECT_SOURCE_DIR}/brigand
              -i${PROJECT_SOURCE_DIR}/highwayhash
              -i${PROJECT_SOURCE_DIR}/exodus
              -i${PROJECT_SOURCE_DIR}/tut
              -i${PROJECT_SOURCE_DIR}/zoltan
              ${CMAKE_CURRENT_SOURCE_DIR}
              ${CMAKE_CURRENT_SOURCE_DIR}/../tests/unit
              2> doc/html/${CMAKE_BUILD_TYPE}/cppcheck/cppcheck-report.xml
      # Generate html output
      COMMAND ${CPPCHECK_HTMLREPORT}
              --file=doc/html/${CMAKE_BUILD_TYPE}/cppcheck/cppcheck-report.xml
              --report-dir=doc/html/${CMAKE_BUILD_TYPE}/cppcheck --source-dir=.
      # Customize page headers in generated html
      COMMAND ${FILEFIND} doc/html/${CMAKE_BUILD_TYPE}/cppcheck -type f -exec ${SED} -i "s/project name/Xyst/g" {} +
      # Set work directory for target
      WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
      # Echo what is being done
      COMMENT "Cppcheck-xml static analysis report"
      VERBATIM USES_TERMINAL
    )
    # Output code coverage target enabled
    message(STATUS "Enabling cppcheck static analysis target 'cppcheck-xml', report at ${CMAKE_BINARY_DIR}/doc/html/${CMAKE_BUILD_TYPE}/cppcheck/index.html")
  endif()

endif()

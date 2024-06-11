project(UnitTest CXX)

# Include function for adding Charm++ modules
include(charm)

add_library(UnitTest
            TUTSuite.cpp
            TUTTest.cpp
            QuietCerr.cpp)

target_include_directories(UnitTest PUBLIC
                           ${PROJECT_SOURCE_DIR}/..
                           ${PROJECT_SOURCE_DIR}/../Base
                           ${PROJECT_SOURCE_DIR}/../Control
                           ${PROJECT_SOURCE_DIR}/../Main
                           ${PROJECT_SOURCE_DIR}/../UnitTest
                           ${PROJECT_SOURCE_DIR}/../Inciter
                           ${PROJECT_SOURCE_DIR}/../Mesh
                           ${PROJECT_SOURCE_DIR}/../LinearSolver
                           ${CHARM_INCLUDE_DIRS}
                           ${PROJECT_BINARY_DIR}/../Base
                           ${PROJECT_BINARY_DIR}/../Main
                           ${PROJECT_BINARY_DIR}/../IO
                           ${PROJECT_BINARY_DIR}/../Inciter
                           ${PROJECT_BINARY_DIR}/../LinearSolver
                           ${PROJECT_BINARY_DIR}/../UnitTest)

addCharmModule( "charmtimer" "UnitTest" )
addCharmModule( "migrated" "UnitTest" )
addCharmModule( "cgreceiver" "UnitTest" )
addCharmModule( "quietcerr" "UnitTest" )
addCharmModule( "tutsuite" "UnitTest" )
addCharmModule( "tuttest" "UnitTest" )

# Add extra dependency of UnitTest on unittestCharmModule. This is required as a
# dependency of UnitTest, TUTSuite, refers to the main Charm++ proxy defined in
# the Charm++ module unittest (in Main/UnitTest.cpp).
add_dependencies("UnitTest" "unittestCharmModule")

# Add extra dependencies of UnitTest on those CharmModules required for testing
# Inciter.
add_dependencies("UnitTest" "discretizationCharmModule")

set_target_properties(UnitTest PROPERTIES LIBRARY_OUTPUT_NAME xyst_unittest)

INSTALL(TARGETS UnitTest
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Runtime
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Runtime
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development
)

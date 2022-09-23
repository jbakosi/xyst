project(UnitTest CXX)

# Include function for adding Charm++ modules
include(charm)

add_library(UnitTest
            Assessment.cpp
            TUTSuite.cpp
            TUTTest.cpp
            QuietCerr.cpp
            MPIRunner.cpp)

target_include_directories(UnitTest PUBLIC
                           ${PROJECT_SOURCE_DIR}/..
                           ${PROJECT_SOURCE_DIR}/../Base
                           ${PROJECT_SOURCE_DIR}/../Control
                           ${PROJECT_SOURCE_DIR}/../Main
                           ${PROJECT_SOURCE_DIR}/../UnitTest
                           ${PROJECT_SOURCE_DIR}/../Inciter
                           ${PROJECT_SOURCE_DIR}/../Mesh
                           ${PROJECT_SOURCE_DIR}/../LinearSolver
                           ${PROJECT_SOURCE_DIR}/../PDE
                           ${TUT_INCLUDE_DIRS}
                           ${HIGHWAYHASH_INCLUDE_DIRS}
                           ${BRIGAND_INCLUDE_DIRS}
                           ${CHARM_INCLUDE_DIRS}
                           ${PEGTL_INCLUDE_DIRS}
                           ${MPI_CXX_INCLUDE_DIRS}
                           ${MPI_CXX_INCLUDE_PATH}
                           ${PROJECT_BINARY_DIR}/../Base
                           ${PROJECT_BINARY_DIR}/../Main
                           ${PROJECT_BINARY_DIR}/../IO
                           ${PROJECT_BINARY_DIR}/../LinSys
                           ${PROJECT_BINARY_DIR}/../LoadBalance
                           ${PROJECT_BINARY_DIR}/../LinearSolver
                           ${PROJECT_BINARY_DIR}/../Inciter
                           ${PROJECT_BINARY_DIR}/../UnitTest)

addCharmModule( "charmchild" "UnitTest" )
addCharmModule( "charmtimer" "UnitTest" )
addCharmModule( "testarray" "UnitTest" )
addCharmModule( "cgreceiver" "UnitTest" )
addCharmModule( "migrated_base" "UnitTest" )
addCharmModule( "quietcerr" "UnitTest" )

addCharmModule( "migrated_inciter" "UnitTest" )
addCharmModule( "tutsuiteinciter" "UnitTest" )
addCharmModule( "tuttestinciter" "UnitTest" )
addCharmModule( "mpirunnerinciter" "UnitTest" )
add_dependencies( "tutsuiteinciterCharmModule" "mpirunnerinciterCharmModule" )

# Add extra dependency of UnitTest on unittestCharmModule. This is required as a
# dependency of UnitTest, TUTSuite, refers to the main Charm++ proxy defined in
# the Charm++ module unittest (in Main/UnitTest.cpp).
add_dependencies("UnitTest" "unittestinciterCharmModule")

# Add extra dependencies of UnitTest on those CharmModules that are required for
# testing Inciter.
add_dependencies("UnitTest" "diagcgCharmModule"
                            "alecgCharmModule"
                            "distfctCharmModule"
                            "aleCharmModule"
                            "conjugategradientsCharmModule"
                            "discretizationCharmModule"
                            "transporterCharmModule")

# Add extra dependency of TUTSuite on charm modules testing advanced array
# element placement using maps. This is required since TUTSuite spawns chare
# arrays testing advanced array element placement maps.
add_dependencies( "tutsuiteinciterCharmModule"
                  "linearmapCharmModule" "unsmeshmapCharmModule" )

set_target_properties(UnitTest PROPERTIES LIBRARY_OUTPUT_NAME xyst_unittest)

INSTALL(TARGETS UnitTest
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Runtime
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Runtime
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development
)

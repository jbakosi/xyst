project(Inciter CXX)

# Include function for adding Charm++ modules
include(charm)

add_library(Inciter
            Transporter.cpp
            Partitioner.cpp
            Discretization.cpp
            Refiner.cpp
            Sorter.cpp
            RieCG.cpp
            ZalCG.cpp
            KozCG.cpp
            DiagReducer.cpp
            GraphReducer.cpp
            PartsReducer.cpp
            NodeDiagnostics.cpp
            IntegralReducer.cpp)

target_include_directories(Inciter PUBLIC
                           ${PROJECT_SOURCE_DIR}/..
                           ${PROJECT_SOURCE_DIR}/../exodus/include
                           ${PROJECT_SOURCE_DIR}/../Base
                           ${PROJECT_SOURCE_DIR}/../Mesh
                           ${PROJECT_SOURCE_DIR}/../IO
                           ${PROJECT_SOURCE_DIR}/../Control
                           ${PROJECT_SOURCE_DIR}/../Main
                           ${PROJECT_SOURCE_DIR}/../Partition
                           ${PROJECT_SOURCE_DIR}/../Statistics
                           ${PROJECT_SOURCE_DIR}/../Inciter
                           ${PROJECT_SOURCE_DIR}/../Physics
                           ${CHARM_INCLUDE_DIRS}
                           ${MPI_C_INCLUDE_DIRS}
                           ${NETCDF_INCLUDE_DIRS}
                           ${PROJECT_BINARY_DIR}/../Inciter
                           ${PROJECT_BINARY_DIR}/../Base
                           ${PROJECT_BINARY_DIR}/../IO
                           ${PROJECT_BINARY_DIR}/../Mesh
                           ${PROJECT_BINARY_DIR}/../Main)

add_library(MeshRefinement
            PUPAMR.cpp
            AMR/mesh_adapter.cpp
            AMR/util.cpp
            AMR/edge.cpp
            AMR/Error.cpp)

target_include_directories(MeshRefinement PUBLIC
                           ${PROJECT_SOURCE_DIR}/..
                           ${PROJECT_SOURCE_DIR}/../Base
                           ${PROJECT_SOURCE_DIR}/../Control
                           ${PROJECT_SOURCE_DIR}/../Mesh
                           ${PROJECT_SOURCE_DIR}/../Inciter
                           ${CHARM_INCLUDE_DIRS}
                           ${TPL_INCLUDE_DIR}
                           ${PROJECT_BINARY_DIR}/../Main)

addCharmModule( "transporter" "Inciter" )
addCharmModule( "partitioner" "Inciter" )
addCharmModule( "discretization" "Inciter" )
addCharmModule( "refiner" "Inciter" )
addCharmModule( "sorter" "Inciter" )
addCharmModule( "riecg" "Inciter" )
addCharmModule( "zalcg" "Inciter" )
addCharmModule( "kozcg" "Inciter" )

# Add extra dependency of Inciter on inciterCharmModule. This is required as
# Inciter refers to the main Charm++ proxy defined in the Charm++ module
# inciter (in Main/Inciter.cpp).
add_dependencies( "Inciter" "inciterCharmModule" )

# Add extra dependency of Transporter charm module on Partitioner charm module.
# This is required so that partitioner.decl.h and partitioner.def.h are
# generated before Transporter including those.
add_dependencies( "transporterCharmModule" "partitionerCharmModule" )

# Add extra dependency of Discretization charm module on MeshWriter charm
# module. This is so meshwriter.decl.h and meshwriter.def.h are generated
# before Discretization which needs meshwriter's type information.
add_dependencies( "discretizationCharmModule" "meshwriterCharmModule" )

set_target_properties(Inciter PROPERTIES LIBRARY_OUTPUT_NAME xyst_inciter)
set_target_properties(MeshRefinement PROPERTIES LIBRARY_OUTPUT_NAME
                      xyst_mesh_refinement)

INSTALL(TARGETS Inciter
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Runtime
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Runtime
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development
)

INSTALL(TARGETS MeshRefinement
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Runtime
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Runtime
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development
)

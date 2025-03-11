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
            LaxCG.cpp
            ZalCG.cpp
            KozCG.cpp
            ChoCG.cpp
            LohCG.cpp
            DiagReducer.cpp
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
                           ${PROJECT_SOURCE_DIR}/../LinearSolver
                           ${PROJECT_SOURCE_DIR}/../Transfer
                           ${CHARM_INCLUDE_DIRS}
                           ${MPI_C_INCLUDE_DIRS}
                           ${NETCDF_INCLUDE_DIRS}
                           ${PROJECT_BINARY_DIR}/../Inciter
                           ${PROJECT_BINARY_DIR}/../Base
                           ${PROJECT_BINARY_DIR}/../IO
                           ${PROJECT_BINARY_DIR}/../Mesh
                           ${PROJECT_BINARY_DIR}/../LinearSolver
                           ${PROJECT_BINARY_DIR}/../Transfer
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
addCharmModule( "laxcg" "Inciter" )
addCharmModule( "zalcg" "Inciter" )
addCharmModule( "kozcg" "Inciter" )
addCharmModule( "chocg" "Inciter" )
addCharmModule( "lohcg" "Inciter" )

# Add extra dependency of Inciter on inciterCharmModule
add_dependencies( "Inciter" "inciterCharmModule" )

# Add extra dependency of Transporter charm module on Partitioner charm module
add_dependencies( "transporterCharmModule" "partitionerCharmModule" )

# Add extra dependency of Discretization charm module on MeshWriter charm module
add_dependencies( "discretizationCharmModule" "meshwriterCharmModule" )

# Add extra dependency of Partitioner charm module on ConjugateGradients charm
# module
add_dependencies( "partitionerCharmModule" "conjugategradientsCharmModule" )

# Add extra dependency of Refiner charm module on ConjugateGradients charm
# module
add_dependencies( "refinerCharmModule" "conjugategradientsCharmModule" )

# Add extra dependency of Sorter charm module on ConjugateGradients charm
# module
add_dependencies( "sorterCharmModule" "conjugategradientsCharmModule" )

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

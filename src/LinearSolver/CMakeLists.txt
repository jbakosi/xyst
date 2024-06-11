project(LinearSolver CXX)

add_library(LinearSolver
            CSR.cpp
            ConjugateGradients.cpp)

target_include_directories(LinearSolver PUBLIC
                           ${PROJECT_SOURCE_DIR}
                           ${PROJECT_SOURCE_DIR}/..
                           ${PROJECT_SOURCE_DIR}/../Base
                           ${PROJECT_SOURCE_DIR}/../Mesh
                           ${CHARM_INCLUDE_DIRS}
                           ${PROJECT_BINARY_DIR}/../Main)

set_target_properties(LinearSolver PROPERTIES LIBRARY_OUTPUT_NAME
                      xyst_linearsolver)

addCharmModule( "conjugategradients" "LinearSolver" )

INSTALL(TARGETS LinearSolver
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Runtime
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Runtime
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development
)

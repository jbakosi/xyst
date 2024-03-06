project(Partition CXX)

include(charm)

add_library(Partition
            ZoltanGeom.cpp
            ZoltanGraph.cpp)

target_include_directories(Partition PUBLIC
                      ${PROJECT_SOURCE_DIR}/..
                      ${PROJECT_SOURCE_DIR}/../Base
                      ${PROJECT_SOURCE_DIR}/../Control
                      ${PROJECT_SOURCE_DIR}/../Mesh
                      ${PROJECT_SOURCE_DIR}/../zoltan/src/include
                      ${PROJECT_SOURCE_DIR}/../zoltan/src/all
                      ${PROJECT_BINARY_DIR}/../Main
                      ${CHARM_INCLUDE_DIRS}
                      ${MPI_C_INCLUDE_DIRS})

set_target_properties(Partition PROPERTIES
                      LIBRARY_OUTPUT_NAME xyst_partition)

INSTALL(TARGETS Partition
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Runtime
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Runtime
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development)

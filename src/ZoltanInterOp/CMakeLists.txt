project(ZoltanInterOp CXX)

include(charm)

add_library(ZoltanInterOp
            ZoltanInterOp.cpp
            ZoltanGeom.cpp
            #ZoltanGraph.cpp
            )

target_include_directories(ZoltanInterOp PUBLIC
                           ${PROJECT_SOURCE_DIR}/..
                           ${PROJECT_SOURCE_DIR}/../Base
                           ${PROJECT_SOURCE_DIR}/../Control
                           ${PROJECT_SOURCE_DIR}/../Mesh
                           ${PROJECT_SOURCE_DIR}/../zoltan/src/include
                           ${PROJECT_SOURCE_DIR}/../zoltan/src/all
                           ${PROJECT_BINARY_DIR}/../Main
                           ${PEGTL_INCLUDE_DIRS}
                           ${CHARM_INCLUDE_DIRS}
                           ${MPI_CXX_INCLUDE_DIRS})

set_target_properties(ZoltanInterOp PROPERTIES
                      LIBRARY_OUTPUT_NAME xyst_zoltaninterop)

add_dependencies(ZoltanInterOp charm)

INSTALL(TARGETS ZoltanInterOp
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Runtime
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Runtime
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development
)

project(Transfer CXX)

# Include function for adding Charm++ modules
include(charm)

add_library(Transfer
            Transfer.cpp
            NodeSearch.cpp)

target_include_directories(Transfer PUBLIC
                           ${PROJECT_SOURCE_DIR}
                           ${PROJECT_SOURCE_DIR}/..
                           ${PROJECT_SOURCE_DIR}/../Base
                           ${PROJECT_SOURCE_DIR}/../Mesh
                           ${PROJECT_BINARY_DIR}
                           ${PROJECT_BINARY_DIR}/../Main
                           ${CHARM_INCLUDE_DIRS})

set_target_properties(Transfer PROPERTIES LIBRARY_OUTPUT_NAME xyst_transfer)

INSTALL(TARGETS Transfer
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Runtime
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Runtime
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development)

addCharmModule( "nodesearch" "Transfer" )
addCharmModule( "transfer" "Transfer" )

# Add extra dependency of Inciter on transferCharmModule
add_dependencies( "Inciter" "transferCharmModule" )

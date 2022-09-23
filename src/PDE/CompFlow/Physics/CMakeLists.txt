add_library(CGCompFlowPhysics
            CGNavierStokes.cpp)

target_include_directories(CGCompFlowPhysics PUBLIC
                           ${PROJECT_SOURCE_DIR}/..
                           ${PROJECT_SOURCE_DIR}/../Base
                           ${PROJECT_SOURCE_DIR}/../Control
                           ${PROJECT_SOURCE_DIR}/../Inciter
                           ${PROJECT_BINARY_DIR}/../Main
                           ${PROJECT_SOURCE_DIR}/../PDE
                           ${PROJECT_SOURCE_DIR}/../PDE/CompFlow
                           ${BRIGAND_INCLUDE_DIRS}
                           ${PEGTL_INCLUDE_DIRS}
                           ${CHARM_INCLUDE_DIRS})

set_target_properties(CGCompFlowPhysics
                      PROPERTIES LIBRARY_OUTPUT_NAME xyst_cgcompflowphysics)

INSTALL(TARGETS CGCompFlowPhysics
  	RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Runtime
  	LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Runtime
 	ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development)

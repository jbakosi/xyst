project(PDE CXX)

add_library(PDE
            PDEStack.cpp
            ConfigureTransport.cpp
            ConfigureCompFlow.cpp
            ConfigureOutVar.cpp
            CGPDE.cpp
)

target_include_directories(PDE PUBLIC
                           ${PROJECT_SOURCE_DIR}
                           ${PROJECT_SOURCE_DIR}/..
                           ${PROJECT_SOURCE_DIR}/../Base
                           ${PROJECT_SOURCE_DIR}/../Control
                           ${PROJECT_SOURCE_DIR}/../Inciter
                           ${PROJECT_SOURCE_DIR}/../Mesh
                           ${PROJECT_SOURCE_DIR}/CompFlow
                           ${BRIGAND_INCLUDE_DIRS}
                           ${PEGTL_INCLUDE_DIRS}
                           ${HIGHWAYHASH_INCLUDE_DIRS}
                           ${CHARM_INCLUDE_DIRS}
                           ${MKL_INCLUDE_DIRS}
                           ${LAPACKE_INCLUDE_DIRS}
                           ${PROJECT_BINARY_DIR}/../Main)

set_target_properties(PDE PROPERTIES LIBRARY_OUTPUT_NAME xyst_pde)

INSTALL(TARGETS PDE
  	RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Runtime
  	LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Runtime
 	ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development)

# Configure targets for physics and problem types
add_subdirectory(Transport/Physics)
add_subdirectory(Transport/Problem)
add_subdirectory(CompFlow/Physics)
add_subdirectory(CompFlow/Problem)

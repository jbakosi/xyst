project(Mesh CXX)

add_library(Mesh
            DerivedData.cpp
            Gradients.cpp
            Reorder.cpp
            CommMap.cpp)

target_include_directories(Mesh PUBLIC
                           ${PROJECT_SOURCE_DIR}/..
                           ${PROJECT_SOURCE_DIR}/../Base
                           ${PROJECT_SOURCE_DIR}/../Control
                           ${PROJECT_SOURCE_DIR}/../Mesh
                           ${PROJECT_BINARY_DIR}/../Main
                           ${TPL_INCLUDE_DIR}
                           ${PEGTL_INCLUDE_DIRS}
                           ${CHARM_INCLUDE_DIRS})

set_target_properties(Mesh PROPERTIES LIBRARY_OUTPUT_NAME xyst_mesh)

add_dependencies(Mesh charm)

INSTALL(TARGETS Mesh
  	RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Runtime
  	LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Runtime
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development)

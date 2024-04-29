project(Physics CXX)

add_library(Physics
            BC.cpp
            Box.cpp
            Problems.cpp
            Riemann.cpp
            Kozak.cpp
            Zalesak.cpp
            Lax.cpp
)

target_include_directories(Physics PUBLIC
                           ${PROJECT_SOURCE_DIR}/..
                           ${PROJECT_SOURCE_DIR}/../Base
                           ${PROJECT_SOURCE_DIR}/../Mesh
                           ${PROJECT_SOURCE_DIR}/../Control
                           ${PROJECT_BINARY_DIR}/../Main
                           ${CHARM_INCLUDE_DIRS})

set_target_properties(Physics PROPERTIES LIBRARY_OUTPUT_NAME xyst_physics)

INSTALL(TARGETS Physics
  	RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Runtime
  	LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Runtime
 	ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development)

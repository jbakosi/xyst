project(Statistics CXX)

add_library(Statistics
            PDFReducer.cpp)

target_include_directories(Statistics PUBLIC
                           ${PROJECT_SOURCE_DIR}/..
                           ${PROJECT_SOURCE_DIR}/../Base
                           ${PROJECT_SOURCE_DIR}/../Control
                           ${PROJECT_BINARY_DIR}/../Main
                           ${CHARM_INCLUDE_DIRS})

set_target_properties(Statistics PROPERTIES
                      LIBRARY_OUTPUT_NAME xyst_statistics)

INSTALL(TARGETS Statistics
  	RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Runtime
  	LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Runtime
 	ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development
)

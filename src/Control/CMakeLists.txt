project(Control C CXX)

#### Inciter control ###########################################################
add_library(InciterControl
            Lua.c
            InciterInputDeck.cpp)

target_include_directories(InciterControl PUBLIC
                           ${PROJECT_SOURCE_DIR}
                           ${PROJECT_SOURCE_DIR}/..
                           ${PROJECT_SOURCE_DIR}/../Base
                           ${PROJECT_BINARY_DIR}/../Main
                           ${CHARM_INCLUDE_DIRS})

set_target_properties(InciterControl PROPERTIES
                      LIBRARY_OUTPUT_NAME xyst_incitercontrol)

INSTALL(TARGETS InciterControl
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Runtime
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Runtime
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development)

project(Control CXX)

#### Control build configuration ###############################################

include(charm)

#### Inciter control ###########################################################
add_library(InciterControl
            FileParser.cpp
            StringParser.cpp
            Inciter/InputDeck/Parser.cpp
            Inciter/CmdLine/Parser.cpp)

target_include_directories(InciterControl PUBLIC
                           ${PROJECT_SOURCE_DIR}
                           ${PROJECT_SOURCE_DIR}/..
                           ${PROJECT_SOURCE_DIR}/../Base
                           ${PROJECT_SOURCE_DIR}/../Mesh
                           ${PROJECT_SOURCE_DIR}/../PDE
                           ${PROJECT_SOURCE_DIR}/../Inciter
                           ${Boost_INCLUDE_DIRS}
                           ${PEGTL_INCLUDE_DIRS}
                           ${CHARM_INCLUDE_DIRS}
                           ${PROJECT_BINARY_DIR}/../Main)

set_target_properties(InciterControl PROPERTIES
                      LIBRARY_OUTPUT_NAME xyst_incitercontrol)

INSTALL(TARGETS InciterControl
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Runtime
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Runtime
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development)

#### UnitTest control ##########################################################
project(UnitTestControl CXX)

add_library(UnitTestControl
            StringParser.cpp
            FileParser.cpp
            UnitTest/CmdLine/Parser.cpp)

target_include_directories(UnitTestControl PUBLIC
                           ${PROJECT_SOURCE_DIR}/..
                           ${PROJECT_SOURCE_DIR}/../Base
                           ${PROJECT_SOURCE_DIR}/../Control
                           ${PROJECT_BINARY_DIR}/../Main
                           ${Boost_INCLUDE_DIRS}
                           ${PEGTL_INCLUDE_DIRS}
                           ${CHARM_INCLUDE_DIRS})

set_target_properties(UnitTestControl PROPERTIES
                      LIBRARY_OUTPUT_NAME xyst_unittestcontrol)

INSTALL(TARGETS UnitTestControl
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Runtime
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Runtime
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development)

#### MeshConv control ##########################################################
project(MeshConvControl CXX)

add_library(MeshConvControl
            StringParser.cpp
            MeshConv/CmdLine/Parser.cpp)

target_include_directories(MeshConvControl PUBLIC
                           ${PROJECT_SOURCE_DIR}/..
                           ${PROJECT_SOURCE_DIR}/../Base
                           ${PROJECT_SOURCE_DIR}/../IO
                           ${PROJECT_SOURCE_DIR}/../Control
                           ${PROJECT_BINARY_DIR}/../Main
                           ${Boost_INCLUDE_DIRS}
                           ${PEGTL_INCLUDE_DIRS}
                           ${CHARM_INCLUDE_DIRS})

set_target_properties(MeshConvControl PROPERTIES
                      LIBRARY_OUTPUT_NAME xyst_meshconvcontrol)

INSTALL(TARGETS MeshConvControl
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Runtime
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Runtime
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development)

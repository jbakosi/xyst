project(xyst CXX)

### Executables build configuration ############################################

# Basic configuration
add_definitions()

# Include function for adding Charm++ modules
include(charm)

if(AMPI)
  set(CHARM_LINKER ${AMPI_CXX_COMPILER})
else()
  set(CHARM_LINKER ${CHARM_COMPILER})
endif()

STRING(REGEX REPLACE "<CMAKE_CXX_COMPILER>"
  "${CHARM_LINKER} -module CommonLBs -module collidecharm ${EXTRA_LINK_ARGS}"
  CMAKE_CXX_LINK_EXECUTABLE "${CMAKE_CXX_LINK_EXECUTABLE}")

include(ConfigExecutable)

# Configure configuration library
add_library(Config ${PROJECT_BINARY_DIR}/XystConfig.cpp)

set_target_properties(Config PROPERTIES LIBRARY_OUTPUT_NAME xyst_config)

INSTALL(TARGETS Config
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Runtime
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Runtime
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development)

# Configure executable initialization library
add_library(Init Init.cpp)

set_target_properties(Init PROPERTIES LIBRARY_OUTPUT_NAME xyst_init)

target_include_directories(Init PUBLIC
                           ${PROJECT_SOURCE_DIR}/..
                           ${PROJECT_SOURCE_DIR}/../Base
                           ${PROJECT_SOURCE_DIR}/../Control
                           ${PROJECT_BINARY_DIR}/../Main
                           ${CHARM_INCLUDE_DIRS}
                           ${MPI_CXX_INCLUDE_DIRS})

INSTALL(TARGETS Init
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Runtime
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Runtime
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development)

# Configure executable targets

## UnitTest executable ########################################################
if (ENABLE_TESTS)
  add_executable(${UNITTEST_EXECUTABLE}
                 UnitTest.cpp
                 ../../tests/unit/Base/TestContainerUtil.cpp
                 ../../tests/unit/Base/TestData.cpp
                 ../../tests/unit/Base/TestException.cpp
                 ../../tests/unit/Base/TestFlip_map.cpp
                 ../../tests/unit/Base/TestHas.cpp
                 ../../tests/unit/Base/TestPUPUtil.cpp
                 ../../tests/unit/Base/TestReader.cpp
                 ../../tests/unit/Base/TestTaggedTuple.cpp
                 ../../tests/unit/Base/TestPrintUtil.cpp
                 ../../tests/unit/Base/TestPrintTaggedTuple.cpp
                 ../../tests/unit/Base/TestPrintTaggedTupleDeep.cpp
                 ../../tests/unit/Base/TestTimer.cpp
                 ../../tests/unit/Base/TestVector.cpp
                 ../../tests/unit/Base/TestWriter.cpp
                 ../../tests/unit/Base/TestLoadDistributor.cpp
                 ../../tests/unit/Inciter/AMR/TestError.cpp
                 ../../tests/unit/IO/TestExodusIIMeshReader.cpp
                 ../../tests/unit/IO/TestMesh.cpp
                 ../../tests/unit/Mesh/TestAround.cpp
                 ../../tests/unit/Mesh/TestDerivedData.cpp
                 ../../tests/unit/Mesh/TestGradients.cpp
                 ../../tests/unit/Mesh/TestReorder.cpp
                 ../../tests/unit/LinearSolver/TestCSR.cpp
                 ../../tests/unit/LinearSolver/TestConjugateGradients.cpp
  )

  target_include_directories(${UNITTEST_EXECUTABLE} PUBLIC
                             ${PROJECT_SOURCE_DIR}/..
                             ${PROJECT_SOURCE_DIR}/../LinearSolver
                             ${PROJECT_SOURCE_DIR}/../UnitTest
                             ${PROJECT_SOURCE_DIR}/../IO
                             ${PROJECT_BINARY_DIR}/../LinearSolver
                             ${PROJECT_BINARY_DIR}/../UnitTest
                             ${PROJECT_BINARY_DIR}/../IO)

  config_executable(${UNITTEST_EXECUTABLE})

  target_link_libraries(${UNITTEST_EXECUTABLE}
                        Base
                        Config
                        Init
                        MeshRefinement
                        UnitTest
                        Mesh
                        NativeMeshIO
                        MeshDetect
                        LinearSolver
                        ExodusIIMeshIO
                        exodus
                        ${MPI_C_LIBRARIES}
                        ${NETCDF_LIBRARIES}
                        ${HDF5_LIBRARIES})

  # Add custom dependencies for UnitTest's main Charm++ module
  addCharmModule( "unittest" "${UNITTEST_EXECUTABLE}" )
endif()

## Inciter executable ##########################################################

add_executable(${INCITER_EXECUTABLE}
               LBSwitch.cpp
               Inciter.cpp
               XystConfig.cpp)

config_executable(${INCITER_EXECUTABLE})

target_include_directories(${INCITER_EXECUTABLE} PUBLIC
                           ${PROJECT_BINARY_DIR}/Main)

target_link_libraries(${INCITER_EXECUTABLE}
                      InciterControl
                      Inciter
                      MeshRefinement
                      Partition
                      Base
                      Config
                      Init
                      IO
                      MeshWriter
                      ExodusIIMeshIO
                      MeshDetect
                      Mesh
                      Physics
                      Statistics
                      LinearSolver
                      Transfer
                      zoltan1
                      exodus
                      ${MPI_C_LIBRARIES}
                      ${NETCDF_LIBRARIES}
                      ${HDF5_LIBRARIES})

# Add custom dependencies for Inciter's main Charm++ module
addCharmModule( "inciter" "${INCITER_EXECUTABLE}" "-I${PROJECT_BINARY_DIR}")
addCharmModule( "lbswitch" "inciterCharmModule" )

add_dependencies( "inciterCharmModule" "meshwriterCharmModule" )

### MeshConv executable ########################################################

add_executable(${MESHCONV_EXECUTABLE}
               MeshConvDriver.cpp
               MeshConv.cpp)

config_executable(${MESHCONV_EXECUTABLE})

target_include_directories(${MESHCONV_EXECUTABLE} PUBLIC
                           ${PROJECT_SOURCE_DIR}/../IO
                           ${PROJECT_BINARY_DIR}/../Base)

target_link_libraries(${MESHCONV_EXECUTABLE}
                      NativeMeshIO
                      ExodusIIMeshIO
                      MeshDetect
                      Mesh
                      Base
                      Config
                      Init
                      exodus
                      ${MPI_C_LIBRARIES}
                      ${NETCDF_LIBRARIES}
                      ${HDF5_LIBRARIES})

# Add custom dependencies for MeshConv's main Charm++ module
addCharmModule( "meshconv" "${MESHCONV_EXECUTABLE}" )

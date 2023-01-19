project(xyst CXX)

### Executables build configuration ############################################

# Basic configuration
add_definitions()

# Include function for adding Charm++ modules
include(charm)

# Link executables with the charmc wrapper
STRING(REGEX REPLACE "<CMAKE_CXX_COMPILER>"
       "${CHARM_COMPILER} -module CommonLBs ${COLLIDECHARM} ${EXTRA_LINK_ARGS} -c++ <CMAKE_CXX_COMPILER>"
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
                           ${PEGTL_INCLUDE_DIRS}
                           ${CHARM_INCLUDE_DIRS}
                           ${BRIGAND_INCLUDE_DIRS}
                           ${PROJECT_BINARY_DIR}/../Main)

add_dependencies( "Init" "charestatecollectorCharmModule" )

INSTALL(TARGETS Init
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Runtime
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Runtime
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development)

# Configure executable targets

## UnitTest executable ########################################################

add_executable(${UNITTEST_EXECUTABLE}
               UnitTestDriver.cpp
               UnitTest.cpp
               ../../tests/unit/Base/TestContainerUtil.cpp
               ../../tests/unit/Base/TestData.cpp
               ../../tests/unit/Base/TestException.cpp
               ../../tests/unit/Base/TestExceptionMPI.cpp
               ../../tests/unit/Base/TestFactory.cpp
               ../../tests/unit/Base/TestFlip_map.cpp
               ../../tests/unit/Base/TestHas.cpp
               ../../tests/unit/Base/TestPrint.cpp
               ../../tests/unit/Base/TestProcessControl.cpp
               ../../tests/unit/Base/TestPUPUtil.cpp
               ../../tests/unit/Base/TestReader.cpp
               ../../tests/unit/Base/TestPrintUtil.cpp
               ../../tests/unit/Base/TestTaggedTuple.cpp
               ../../tests/unit/Base/TestTaggedTuplePrint.cpp
               ../../tests/unit/Base/TestTaggedTupleDeepPrint.cpp
               ../../tests/unit/Base/TestTimer.cpp
               ../../tests/unit/Base/TestVector.cpp
               ../../tests/unit/Base/TestWriter.cpp
               ../../tests/unit/Control/TestFileParser.cpp
               ../../tests/unit/Control/TestStringParser.cpp
               ../../tests/unit/Control/TestSystemComponents.cpp
               ../../tests/unit/Control/TestToggle.cpp
               ../../tests/unit/Inciter/AMR/TestError.cpp
               ../../tests/unit/IO/TestExodusIIMeshReader.cpp
               ../../tests/unit/IO/TestMesh.cpp
               ../../tests/unit/IO/TestMeshReader.cpp
               ../../tests/unit/LoadBalance/TestLoadDistributor.cpp
               ../../tests/unit/Mesh/TestAround.cpp
               ../../tests/unit/Mesh/TestDerivedData.cpp
               ../../tests/unit/Mesh/TestDerivedData_MPISingle.cpp
               ../../tests/unit/Mesh/TestGradients.cpp
               ../../tests/unit/Mesh/TestReorder.cpp
)

target_include_directories(${UNITTEST_EXECUTABLE} PUBLIC
                           ${PROJECT_SOURCE_DIR}/..
                           ${PROJECT_SOURCE_DIR}/../UnitTest
                           ${PROJECT_SOURCE_DIR}/../LoadBalance
                           ${PROJECT_SOURCE_DIR}/../IO
                           ${TUT_INCLUDE_DIRS}
                           ${LAPACKE_INCLUDE_DIRS}
                           ${RANDOM123_INCLUDE_DIRS}
                           ${PROJECT_BINARY_DIR}/../UnitTest
                           ${PROJECT_BINARY_DIR}/../IO)

config_executable(${UNITTEST_EXECUTABLE})

target_link_libraries(${UNITTEST_EXECUTABLE}
                      Base
                      Config
                      Init
                      MeshRefinement
                      UnitTest
                      UnitTestControl
                      Mesh
                      MeshDetect
                      NativeMeshIO
                      ExodusIIMeshIO
                      HyperMeshIO
                      ${PUGIXML_LIBRARIES}
                      ${SEACASExodus_LIBRARIES}
                      ${MKL_CORE_LIBRARY}
                      ${MKL_INTERFACE_LIBRARY}
                      ${MKL_SEQUENTIAL_LAYER_LIBRARY}
                      ${HDF5_C_LIBRARIES}
                      Backward::Backward
)

# Add custom dependencies for UnitTest's main Charm++ module
addCharmModule( "unittestinciter" "${UNITTEST_EXECUTABLE}" )
add_dependencies( "unittestinciterCharmModule"
                  "charestatecollectorCharmModule"
                  "mpirunnerinciterCharmModule" )

## Inciter executable ##########################################################

add_executable(${INCITER_EXECUTABLE}
               InciterDriver.cpp
               InciterPrint.cpp
               LBSwitch.cpp
               Inciter.cpp)

config_executable(${INCITER_EXECUTABLE})

target_link_libraries(${INCITER_EXECUTABLE}
                      InciterControl
                      Inciter
                      MeshRefinement
                      ZoltanInterOp
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
                      ${EXAM2M_LIBRARIES}
                      ${SEACASExodus_LIBRARIES}
                      ${Zoltan2_LIBRARIES}
                      ${LAPACKE_LIBRARIES}
                      Backward::Backward
)

# Add custom dependencies for Inciter's main Charm++ module
addCharmModule( "inciter" "${INCITER_EXECUTABLE}" "-I${PROJECT_BINARY_DIR}")
addCharmModule( "lbswitch" "inciterCharmModule" )

add_dependencies( "inciterCharmModule" "charestatecollectorCharmModule" )
add_dependencies( "inciterCharmModule" "meshwriterCharmModule" )

### MeshConv executable ########################################################

add_executable(${MESHCONV_EXECUTABLE}
               MeshConvDriver.cpp
               MeshConv.cpp)

config_executable(${MESHCONV_EXECUTABLE})

target_include_directories(${MESHCONV_EXECUTABLE} PUBLIC
                           ${PROJECT_BINARY_DIR}/../Base)

target_link_libraries(${MESHCONV_EXECUTABLE}
                      NativeMeshIO
                      ExodusIIMeshIO
                      HyperMeshIO
                      MeshDetect
                      Mesh
                      MeshConvControl
                      Base
                      Config
                      Init
                      ${PUGIXML_LIBRARIES}
                      ${SEACASExodus_LIBRARIES}
                      ${HDF5_C_LIBRARIES}
                      Backward::Backward
)

# Add custom dependencies for MeshConv's main Charm++ module
addCharmModule( "meshconv" "${MESHCONV_EXECUTABLE}" )

add_dependencies( "meshconvCharmModule" "charestatecollectorCharmModule" )

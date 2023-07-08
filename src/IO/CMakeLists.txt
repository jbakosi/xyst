project(IO CXX)

add_library(IO
            PDFWriter.cpp
            DiagWriter.cpp)

target_include_directories(IO PUBLIC
                           ${PROJECT_SOURCE_DIR}/..
                           ${PROJECT_SOURCE_DIR}/../Base
                           ${PROJECT_SOURCE_DIR}/../Control
                           ${PROJECT_SOURCE_DIR}/../Statistics
                           ${PEGTL_INCLUDE_DIRS}
                           ${CHARM_INCLUDE_DIRS}
                           ${MPI_C_INCLUDE_DIRS}
                           ${NETCDF_INCLUDE_DIRS}
                           ${PROJECT_BINARY_DIR}/../Main)

set_target_properties(IO PROPERTIES LIBRARY_OUTPUT_NAME xyst_io)

INSTALL(TARGETS IO
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Runtime
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Runtime
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development
)

# Native IO libraries do not require third-party libraries
add_library(NativeMeshIO
            MeshFactory.cpp
            GmshMeshReader.cpp
            NetgenMeshReader.cpp
            ASCMeshReader.cpp
            UGRIDMeshReader.cpp
            RDGFLOMeshReader.cpp
            MeditMeshReader.cpp
            GmshMeshWriter.cpp
            NetgenMeshWriter.cpp
)

target_include_directories(NativeMeshIO PUBLIC
                           ${PROJECT_SOURCE_DIR}/..
                           ${PROJECT_SOURCE_DIR}/../Base
                           ${PROJECT_SOURCE_DIR}/../Mesh
                           ${PROJECT_SOURCE_DIR}/../Control
                           ${CHARM_INCLUDE_DIRS}
                           ${MPI_C_INCLUDE_DIRS}
                           ${PEGTL_INCLUDE_DIRS}
                           ${NETCDF_INCLUDE_DIRS}
                           ${PROJECT_BINARY_DIR}/../Main)

set_target_properties(NativeMeshIO PROPERTIES LIBRARY_OUTPUT_NAME
                      xyst_nativemeshio)

INSTALL(TARGETS NativeMeshIO
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Runtime
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Runtime
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development
)

add_library(ExodusIIMeshIO
            ExodusIIMeshReader.cpp
            ExodusIIMeshWriter.cpp)

target_include_directories(ExodusIIMeshIO PUBLIC
                           ${PROJECT_SOURCE_DIR}/..
                           ${PROJECT_SOURCE_DIR}/../Base
                           ${PROJECT_SOURCE_DIR}/../Mesh
                           ${MPI_C_INCLUDE_DIRS}
                           ${NETCDF_INCLUDE_DIRS})

set_target_properties(ExodusIIMeshIO PROPERTIES LIBRARY_OUTPUT_NAME
                      xyst_exodusiimeshio)

INSTALL(TARGETS ExodusIIMeshIO
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Runtime
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Runtime
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development
)

add_library(MeshDetect
            MeshDetect.cpp)

target_include_directories(MeshDetect PUBLIC
                           ${PROJECT_SOURCE_DIR}/..
                           ${PROJECT_SOURCE_DIR}/../Base)

set_target_properties(MeshDetect PROPERTIES LIBRARY_OUTPUT_NAME
                      xyst_meshdetect)

INSTALL(TARGETS MeshDetect
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Runtime
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Runtime
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development
)

add_library(MeshWriter
            MeshWriter.cpp)

target_include_directories(MeshWriter PUBLIC
                           ${PROJECT_SOURCE_DIR}/..
                           ${PROJECT_SOURCE_DIR}/../Base
                           ${PROJECT_SOURCE_DIR}/../Mesh
                           ${PROJECT_SOURCE_DIR}/../Control
                           ${PROJECT_BINARY_DIR}/../Main
                           ${PROJECT_BINARY_DIR}/../IO
                           ${CHARM_INCLUDE_DIRS}
                           ${PEGTL_INCLUDE_DIRS})

set_target_properties(MeshWriter PROPERTIES LIBRARY_OUTPUT_NAME
                      xyst_meshwriter)

INSTALL(TARGETS MeshWriter
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Runtime
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Runtime
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development
)

addCharmModule( "meshwriter" "MeshWriter" )

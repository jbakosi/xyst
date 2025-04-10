project(Base CXX)

# Include function for adding Charm++ modules
include(charm)

add_library(Base
            Exception.cpp
            ProcessException.cpp
            LoadDistributor.cpp
            Timer.cpp
            Reader.cpp
            Writer.cpp
)

target_include_directories(Base PUBLIC
                           ${PROJECT_SOURCE_DIR}
                           ${PROJECT_SOURCE_DIR}/..
                           ${PROJECT_SOURCE_DIR}/../Control
                           ${PROJECT_BINARY_DIR}/../Main
                           ${CHARM_INCLUDE_DIRS}
                           ${MPI_C_INCLUDE_DIRS})

set_target_properties(Base PROPERTIES LIBRARY_OUTPUT_NAME xyst_base)

INSTALL(TARGETS Base
       RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Runtime
       LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Runtime
       ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development
)

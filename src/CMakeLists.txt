# Build xyst

cmake_minimum_required(VERSION 3.1.0 FATAL_ERROR)

project(xyst C CXX)

# Set cmake modules directory
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../cmake")

# Disallow in-source builds
include(DisallowInSourceBuilds)
# Set CMAKE_BUILD_TYPE
include(BuildType)
# Find MPI compilers
include(MPICompilers)

set(BUILD_SHARED_LIBS ON CACHE BOOL "Build shared libraries. Possible values: ON | OFF")

# Set compilers, will become compiler defines
set(COMPILER ${UNDERLYING_CXX_COMPILER})
set(MPI_COMPILER ${MPI_CXX_COMPILER})
MESSAGE(STATUS "Underlying C++ compiler: " ${UNDERLYING_CXX_COMPILER})
# Set compiler id
string(TOLOWER "${CMAKE_CXX_COMPILER_ID}" COMP)

# set version
set(VERSION_MAJOR "1")
set(VERSION_MINOR "0")
set(PROJECT_VERSION ${VERSION_MAJOR}.${VERSION_MINOR})

# Set the requirement for the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
message(STATUS "Required language standard: C++${CMAKE_CXX_STANDARD}")

# Set build type (will be exported to the source)
set(BUILD_TYPE ${CMAKE_BUILD_TYPE})
string(TOLOWER "${BUILD_TYPE}" BUILD_TYPE)

# Echo install prefix
message(STATUS "CMAKE_INSTALL_PREFIX: " ${CMAKE_INSTALL_PREFIX})

macro (add_compiler_flag flag)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${flag}")
endmacro (add_compiler_flag)

if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")

  add_compiler_flag("-fcolor-diagnostics")
  add_compiler_flag("-Weverything")
  add_compiler_flag("-Wno-c++98-compat")
  add_compiler_flag("-Wno-c++98-compat-pedantic")
  add_compiler_flag("-Wno-padded")
  add_compiler_flag("-Wno-global-constructors")
  add_compiler_flag("-Wno-weak-vtables")
  add_compiler_flag("-Wno-reserved-id-macro")
  add_compiler_flag("-Wno-unknown-pragmas")
  add_compiler_flag("-Wno-exit-time-destructors")
  add_compiler_flag("-Wno-source-uses-openmp")

  # Linker flags for clang
  set(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS)      # clear link flags
  set(CMAKE_EXE_LINKER_FLAGS
    "${CMAKE_EXE_LINKER_FLAGS} -Wno-missing-prototypes -Wno-unused-parameter -Wno-unused-command-line-argument")

elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")

  add_compiler_flag("-fdiagnostics-color")
  add_compiler_flag("-Wall")
  add_compiler_flag("-Wextra")
  add_compiler_flag("-Wno-unknown-pragmas")

endif()

# Include cmake code to enable configuration for data layouts
include(ConfigureDataLayout)

# Grab current date to stick into the executables and doc (footer)
execute_process(COMMAND "date" OUTPUT_VARIABLE BUILD_DATE)
string(REGEX REPLACE "[\r\n]" "" BUILD_DATE "${BUILD_DATE}")

execute_process(COMMAND "date" "+%A, %b %d, %Y"
                OUTPUT_VARIABLE BUILD_DATE_SIMPLE)
string(REGEX REPLACE "[\r\n]" "" BUILD_DATE_SIMPLE "${BUILD_DATE_SIMPLE}")

option(SMP "Use Charm++ built in SMP mode")
if(SMP)
  set(smp "-smp")
endif()

# Set external libraries install directory
if (NOT TPL_DIR)
  set(REL_TPL_DIR ${PROJECT_SOURCE_DIR}/../external/install/${COMP}${smp})
  get_filename_component(ABS_TPL_DIR ${REL_TPL_DIR} ABSOLUTE)
  set(TPL_DIR ${ABS_TPL_DIR} CACHE
      STRING "Third-party (external) libraries install directory.")
endif()
message(STATUS "Third-party (external) libraries install directory: " ${TPL_DIR})

# Include support for multiarch path names
include(GNUInstallDirs)
# Add our TPL_DIR to modules directory
set(CMAKE_PREFIX_PATH ${TPL_DIR} ${CMAKE_PREFIX_PATH})

find_package(Charm REQUIRED)
find_package(Boost REQUIRED)
find_package(TUT REQUIRED)
find_package(PEGTL 2.0.0 REQUIRED)
find_package(NumDiff REQUIRED)
find_package(LAPACKE REQUIRED)

find_package(Backward)
if(Backward_FOUND)
  set(HAS_BACKWARD true)  # will become compiler define
endif()

find_package(HighwayHash REQUIRED)
find_package(Brigand REQUIRED)
find_package(Pstreams REQUIRED)
find_package(TUT REQUIRED)
find_package(Pugixml REQUIRED)

set(HDF5_PREFER_PARALLEL true)
if(NOT BUILD_SHARED_LIBS)
  set(HDF5_USE_STATIC_LIBRARIES true)
endif()
find_package(HDF5 COMPONENTS C HL REQUIRED)
find_package(NetCDF REQUIRED)

find_package(Zoltan2 REQUIRED)
find_package(SEACASExodus REQUIRED)
find_package(Exodiff REQUIRED)

# Find Python: required for code coverage (fastcov) and doc (m.css)
find_package(PythonInterp 3.6)

# Attempt to find tools required for code coverage analysis
include(DetectCodeCoverage)

# Set compiler flags to perform code coverage analysis
if(CODE_COVERAGE)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-arcs -ftest-coverage")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} --coverage")
endif()

# Query build hostname
site_name(HOSTNAME)

# Set runner configuration. This is used when generating code coverage and
# running regression tests using ctest. Note that ninja puts chamrun to <build>
# while make puts it into <build>/Main, but after any executable target we copy
# <build>/Main/charmrun to <build> for non-ninja builds. See also
# cmake/ConfigExecutable.cmake.

# Defaults for runner, used to run regression tests
set(RUNNER "${CHARM_RUN}" CACHE STRING
    "Parallel/serial job runner, used to run regression tests")
set(RUNNER_NCPUS_ARG "+p" CACHE STRING
    "The argument used to specify the number of PEs for the runner")
set(RUNNER_ARGS "" CACHE STRING
    "Arguments to pass to the parallel/serial job runner, e.g., to mpirun.")
set(POSTFIX_RUNNER_ARGS "" CACHE STRING
    "Arguments to pass to the parallel/serial job runner at the end of the command line.")

# Echo runner configuration
message(STATUS "Parallel/serial job runner: ${RUNNER}")
message(STATUS "Parallel job runner ncpus arg: ${RUNNER_NCPUS_ARG}")
if(RUNNER_ARGS)
  message(STATUS "Parallel/serial job runner args: ${RUNNER_ARGS}")
endif()
if(POSTFIX_RUNNER_ARGS)
  message(STATUS "Parallel/serial job runner postfix args: ${POSTFIX_RUNNER_ARGS}")
endif()

# Echo compiler flags
if(CMAKE_C_FLAGS)
  message(STATUS "C compiler flags: '${CMAKE_C_FLAGS}'")
endif()
if(CMAKE_CXX_FLAGS)
  message(STATUS "C++ compiler flags: '${CMAKE_CXX_FLAGS}'")
endif()
if(CMAKE_EXE_LINKER_FLAGS)
  message(STATUS "Executable linker flags: '${CMAKE_EXE_LINKER_FLAGS}'")
endif()

# Query regression tests directory (will be exported to the source)
get_filename_component(REGRESSION_DIR
  "${PROJECT_SOURCE_DIR}/../tests/regression" ABSOLUTE)

set(INCITER_EXECUTABLE inciter)
set(MESHCONV_EXECUTABLE meshconv)
set(UNITTEST_EXECUTABLE unittest)

# Create set of all executables tested. These are the executables that do
# "interesting things", i.e., excluding unittest
set(NONTEST_EXECUTABLES ${MESHCONV_EXECUTABLE} ${INCITER_EXECUTABLE})

set(EXECUTABLES ${MESHCONV_EXECUTABLE} ${UNITTEST_EXECUTABLE} ${INCITER_EXECUTABLE})

# Configure cmake variables to pass to the build
configure_file( "${PROJECT_SOURCE_DIR}/Main/XystConfig.hpp.in"
                "${PROJECT_BINARY_DIR}/Main/XystConfig.hpp" )
configure_file( "${PROJECT_SOURCE_DIR}/Main/XystBuildConfig.hpp.in"
                "${PROJECT_BINARY_DIR}/Main/XystBuildConfig.hpp" )
configure_file( "${PROJECT_SOURCE_DIR}/Main/XystConfig.cpp.in"
                "${PROJECT_BINARY_DIR}/Main/XystConfig.cpp" ESCAPE_QUOTES)

add_subdirectory(Base)
add_subdirectory(Control)
add_subdirectory(LoadBalance)
add_subdirectory(Statistics)
add_subdirectory(Mesh)
add_subdirectory(IO)
add_subdirectory(Physics)
add_subdirectory(Inciter)
add_subdirectory(UnitTest)
add_subdirectory(Main)

# Configure regression tests
enable_testing()
add_subdirectory(${REGRESSION_DIR}
                 ${CMAKE_BINARY_DIR}/tests
                 EXCLUDE_FROM_ALL)

# Configure building documentation and web page
find_package(Doxygen 1.8.15)
find_package(MCSS)
if (DOXYGEN_FOUND AND MCSS_FOUND AND PYTHONINTERP_FOUND AND PYGMENTS_FOUND AND
    JINJA2_FOUND AND LATEX_FOUND)

  find_program( PERL perl )
  find_program( AWK awk )
  find_program( XARGS xargs )
  find_program( GIT git )
  find_program( HEAD head )
  find_program( EXTRACT_CMD_KEYWORDS
                ${CMAKE_CURRENT_SOURCE_DIR}/../tools/extract_cmd_keywords.sh )
  find_program( EXTRACT_CTR_KEYWORDS
                ${CMAKE_CURRENT_SOURCE_DIR}/../tools/extract_ctr_keywords.sh )

  if (PERL AND AWK AND XARGS AND GIT AND HEAD)

    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/../doc/xyst.doxy
                   ${CMAKE_CURRENT_BINARY_DIR}/doc/xyst.doxy @ONLY)

    file(GLOB doc_pages CONFIGURE_DEPENDS
         RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" "../doc/pages/*.dox")
    foreach(page ${doc_pages})
      configure_file(${page} ${CMAKE_CURRENT_BINARY_DIR}/doc/${page} @ONLY)
    endforeach()

    file(COPY ../doc/images DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/doc/html)

    # Generate targets <executable>_cmd_keywords for extracting and generating
    # reference documentation pages for command line keywords
    set(executables_cmd_keywords)
    foreach(exe IN LISTS EXECUTABLES)
      add_custom_target(${exe}_cmd_keywords
                        DEPENDS ${exe} ${doc_pages}
                        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/doc
                        COMMAND ${EXTRACT_CMD_KEYWORDS} ${exe}
                        COMMENT "Extract ${exe}'s command-line keywords")
      list(APPEND executables_cmd_keywords "${exe}_cmd_keywords")
    endforeach()

    # Generate targets <executable>_ctr_keywords for extracting and generating
    # reference documentation pages for control file keywords
    set(executables_ctr_keywords)
    foreach(exe IN LISTS EXECUTABLES)
      add_custom_target(${exe}_ctr_keywords
                        DEPENDS ${exe} ${doc_pages}
                        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/doc
                        COMMAND ${EXTRACT_CTR_KEYWORDS} ${exe}
                        COMMENT "Extract ${exe}'s control file keywords")
      list(APPEND executables_ctr_keywords "${exe}_ctr_keywords")
    endforeach()

    add_custom_command(OUTPUT xml
                       WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/doc
                       DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/doc/xyst.doxy
                               ${executables_cmd_keywords}
                               ${executables_ctr_keywords}
                       COMMAND ${DOXYGEN_EXECUTABLE} xyst.doxy
                       COMMAND ${MCSS_DOX2HTML5} --no-doxygen xyst.doxy
                       VERBATIM USES_TERMINAL)

    add_custom_target(doc DEPENDS xml
                      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/doc
                      COMMENT "Documentation")

    message(STATUS "Add target 'doc' to generate documentation to ./doc/html/index.html")

  endif()

endif()

# Setup code coverage targets
if(CODE_COVERAGE)

  # Determine number of CPUs available
  if (DEFINED ENV{SLURM_NPROCS})
    set(PROCESSOR_COUNT $ENV{SLURM_NPROCS})
  else()
    include(ProcessorCount)
    ProcessorCount(PROCESSOR_COUNT)
    if(PROCESSOR_COUNT EQUAL 0)     # if ncpus could not be determined, assign 1
      set(PROCESSOR_COUNT "1")
    endif()
  endif()
  message(STATUS "Found ${PROCESSOR_COUNT} PEs")

  include(CppCheck)
  include(CodeCoverage)
  setup_target_for_all_coverage(doc/html/${CMAKE_BUILD_TYPE} test_coverage
                                ${RUNNER} ${RUNNER_NCPUS_ARG}
                                DEPENDS ${NONTEST_EXECUTABLES}
                                        ${UNITTEST_EXECUTABLE})
endif()

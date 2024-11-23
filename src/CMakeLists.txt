# Build xyst

cmake_minimum_required(VERSION 3.1.0 FATAL_ERROR)

project(xyst C CXX)

# Set cmake modules directory
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../cmake")

# Disallow in-source builds
include(DisallowInSourceBuilds)
# Set CMAKE_BUILD_TYPE
include(BuildType)

option(AMPI "Use Charm++'s AMPI instead of MPI" OFF)

if (NOT AMPI)
  find_package(MPI REQUIRED)
  set(COMPILER ${MPI_CXX_COMPILER})
endif()

set(BUILD_SHARED_LIBS ON CACHE
    BOOL "Build shared libraries. Possible values: ON | OFF")

string(TOLOWER "${CMAKE_CXX_COMPILER_ID}" COMP)

# Set the requirement for the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
message(STATUS "Language standard: C++${CMAKE_CXX_STANDARD}")

# Set build type (will be exported to the source)
set(BUILD_TYPE ${CMAKE_BUILD_TYPE})
string(TOLOWER "${BUILD_TYPE}" BUILD_TYPE)

# Echo install prefix
message(STATUS "Install prefix: " ${CMAKE_INSTALL_PREFIX})

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
  add_compiler_flag("-Wno-unknown-warning-option")
  add_compiler_flag("-Wno-reserved-identifier")
  add_compiler_flag("-Wno-reserved-id-macro")
  add_compiler_flag("-Wno-exit-time-destructors")
  add_compiler_flag("-Wno-source-uses-openmp")

  # Override default optimization level for release builds
  if(CMAKE_BUILD_TYPE MATCHES "REL")
    set(CMAKE_CXX_FLAGS_RELEASE "-O1 -DNDEBUG" CACHE STRING
        "CMAKE_CXX_FLAGS_RELEASE: Flags used by the compiler during RELEASE builds."
        FORCE)
    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O1 -g -DNDEBUG" CACHE STRING
        "CMAKE_CXX_FLAGS_RELWITHDEBINFO: Flags used by the compiler during RELWITHDEBINFO info."
        FORCE)
  endif()

  # Linker flags for clang
  set(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS)      # clear link flags
  set(CMAKE_EXE_LINKER_FLAGS
    "${CMAKE_EXE_LINKER_FLAGS} -Wno-missing-prototypes -Wno-unused-parameter -Wno-unused-command-line-argument")

elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")

  add_compiler_flag("-fdiagnostics-color")
  add_compiler_flag("-Wall")
  add_compiler_flag("-Wextra")

elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")

  add_compiler_flag("-w3")       # enable diagnostics: remarks, warnings, errors
  add_compiler_flag("-Wremarks") # display remarks and comments

  # Globally suppress the following Intel compiler remarks:
  # 7: unrecognized token
  # 68: integer conversion resulted in a change of sign
  # 82: storage class is not first
  # 111: statement is unreachable
  # 186: pointless comparison of unsigned integer with zero
  # 177: function was declared but never referenced
  # 193: zero used for undefined preprocessing identifier
  # 280: selector expression is constant
  # 304: access control not specified ("public" by default)
  # 383: value copied to temporary, reference to temporary used
  # 327: NULL reference is not allowed
  # 411: defines no constructor for initializer_list<>
  # 444: destructor for base class is not virtual
  # 593: parameter was set but never used
  # 688: "simd" is an unrecognized __declspec attribute
  # 869: parameter was never referenced
  # 981: operands are evaluated in unspecified order
  # 1011: missing return statement at end of non-void function (constexpr)
  # 1292: unknown attribute "fallthrough"
  # 1418: external function definition with no prior declaration
  # 1419: external declaration in primary source file
  # 1572: floating-point equality and inequality comparisons are unreliable
  # 1599: declaration hides variable
  # 2259: non-pointer conversion from "int" to "char" may lose significant bits
  # 2407: the initialization of class will be done before that of member
  # 2415: variable of static storage duration was declared but never referenced
  # 2547: non-system include dir will be ignored when both a system and non-system include directory were specified
  # 2557: comparison between signed and unsigned operands
  # 3180: unrecognized OpenMP #pragma
  # 3280: declaration hides member
  # 3346: dynamic exception specifications are deprecated
  # 3656: variable may be used before its value is set
  # 10120: overriding '-O3' with '-O2'
  # 10397: optimization reports are generated in *.optrpt files in the output location
  # 15518: Structure assignment was serialized
  add_compiler_flag("-diag-disable 7,68,82,111,177,186,193,280,304,327,383,411,444,593,688,869,981,1011,1292,1418,1419,1572,1599,2259,2407,2415,2547,2557,3180,3280,3346,3656,10120,10397,15518,2282")

  # default: no openmp simd (regression tests would fail otherwise)
  add_compiler_flag("-qno-openmp-simd")
  # for vectorization: disable above and enable some below
  #add_compiler_flag("-march=native")
  #add_compiler_flag("-qopt-zmm-usage=high")
  #add_compiler_flag("-vecabi=cmdtarget")
  #add_compiler_flag("-vec-threshold0")
  #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -qopt-report=5 -qopt-report-phase=loop,vec")

  # No limit on inlining. This can increase the build time.
  if(CMAKE_BUILD_TYPE MATCHES "REL") # Release | RelWithDebInfo | MinSizeRel
    #add_compiler_flag("-no-inline-max-size")
    #add_compiler_flag("-no-inline-max-total-size")
    #add_compiler_flag("-no-inline-max-per-routine")
    add_compiler_flag("-diag-disable 11074,11076")
  endif()

endif()

# Include cmake code to enable configuration for data layouts
include(ConfigureDataLayout)

# Grab current date to stick into the executables and doc (footer)
execute_process(COMMAND "date" OUTPUT_VARIABLE BUILD_DATE)
string(REGEX REPLACE "[\r\n]" "" BUILD_DATE "${BUILD_DATE}")

execute_process(COMMAND "date" "+%A, %b %d, %Y"
                OUTPUT_VARIABLE BUILD_DATE_SIMPLE)
string(REGEX REPLACE "[\r\n]" "" BUILD_DATE_SIMPLE "${BUILD_DATE_SIMPLE}")

# Include support for multiarch path names
include(GNUInstallDirs)

include(ProcessorCount)
ProcessorCount(PROCESSOR_COUNT)
if(PROCESSOR_COUNT EQUAL 0)
  set(PROCESSOR_COUNT "1")
endif()
message(STATUS "Found ${PROCESSOR_COUNT} PEs")

# Build external project: PROJECT: project name, COMMAND: cmake build code
function(build_external_project PROJECT COMMAND)
  set(dir ${CMAKE_BINARY_DIR}/${PROJECT})
  file(MAKE_DIRECTORY ${dir} ${dir}/build)
  set(CMAKE_LIST_CONTENT "${COMMAND}")
  file(WRITE ${dir}/CMakeLists.txt "${CMAKE_LIST_CONTENT}")
  execute_process(COMMAND ${CMAKE_COMMAND} ..
                  WORKING_DIRECTORY ${dir}/build)
  execute_process(COMMAND ${CMAKE_COMMAND} --build .
                  WORKING_DIRECTORY ${dir}/build)
endfunction()

option(SMP "Build Charm++ in SMP mode" OFF)
set(CHARM_OPTS "" CACHE STRING "Forward options to Charm++ cmake")

if(AMPI)
  set(CHARM_TARGET "AMPI")
  set(CHARM_NETWORK "netlrts")
  set(XYST_AMPI AMPI)   # will become compiler define
else()
  set(CHARM_TARGET "Charm++")
  set(CHARM_NETWORK "mpi")
endif()

find_package(Charm)
if(NOT CHARM_FOUND)
  set(CHARM_LOCAL "${PROJECT_SOURCE_DIR}/../charm")
  if(EXISTS ${CHARM_LOCAL})
    set(CHARM_REPOSITORY "${CHARM_LOCAL}")
  else()
    set(CHARM_REPOSITORY "https://github.com/charmplusplus/charm")
  endif()
  message(STATUS "Will build charm++ from: ${CHARM_REPOSITORY}")
  build_external_project(charm "
    cmake_minimum_required(VERSION 3.0)
    project(charm CXX)
    include(ExternalProject)
    ExternalProject_add(charm
      PREFIX charm
      GIT_REPOSITORY ${CHARM_REPOSITORY}
      GIT_TAG v7.0.0-rc2
      PATCH_COMMAND patch -t -i ${PROJECT_SOURCE_DIR}/../cmake/charm-fortran.patch
      CMAKE_ARGS -D TARGET=${CHARM_TARGET}
                 -D NETWORK=${CHARM_NETWORK}
                 -D SMP=${SMP}
                 -D ENABLE_FORTRAN=off
                 #-D SYNCFT=on
                 -D CMAKE_BUILD_TYPE=Release ${CHARM_OPTS}
                 -D CMAKE_INSTALL_PREFIX:PATH=${CMAKE_BINARY_DIR}/charm/install
      CMAKE_GENERATOR \"Unix Makefiles\"
      BUILD_COMMAND make -j${PROCESSOR_COUNT}
      LOG_DOWNLOAD 1
      LOG_CONFIGURE 1
      LOG_BUILD 1
      LOG_INSTALL 1)
  ")
  find_package(Charm REQUIRED)
endif()

if(AMPI)
  set(COMPILER ${AMPI_CXX_COMPILER})
endif()

message(STATUS "C++ compiler: " ${COMPILER})

find_package(NetCDF REQUIRED)
if (FIND_HDF5)
  find_package(HDF5 COMPONENTS HL REQUIRED)
endif()

add_subdirectory(exodus)
add_subdirectory(zoltan)

# get git sha1
include(GitCommit)
RunGitCommand(show -s "--format=%H" HEAD)
if(exit_code EQUAL 0)
  set(GIT_HEAD_SHA1 ${output})
endif()
if(NOT ${GIT_HEAD_SHA1} STREQUAL "")
  message(STATUS "Git commit: ${GIT_HEAD_SHA1}")
endif()

option(BUILD_DOC "Build docmentation" OFF)
if (BUILD_DOC)
  # Find Python: required for code coverage (fastcov) and doc (m.css)
  find_package(PythonInterp 3.6 QUIET)
  # Attempt to find tools required for code coverage analysis
  include(DetectCodeCoverage)
  # Set compiler flags to perform code coverage analysis
  if(CODE_COVERAGE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-arcs -ftest-coverage")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} --coverage")
  endif()
endif()

# Query build hostname
site_name(HOSTNAME)

option(ENABLE_TESTS "Enable tests" ON)
if (ENABLE_TESTS)
  # Set runner configuration. This is used when generating code coverage and
  # running regression tests using ctest. Note that ninja puts chamrun to
  # <build> while make puts it into <build>/Main, but after any executable
  # target we copy <build>/Main/charmrun to <build> for non-ninja builds. See
  # also cmake/ConfigExecutable.cmake.

  # Defaults for runner, used to run regression tests
  set(RUNNER "${CHARM_RUN}" CACHE STRING
      "Parallel/serial test runner, used to run regression tests")
  set(RUNNER_NCPUS_ARG "+p" CACHE STRING
      "The argument used to specify the number of PEs for the test runner")
  set(RUNNER_ARGS "" CACHE STRING
      "Arguments to pass to the parallel/serial test runner, e.g., to mpirun.")
  set(POSTFIX_RUNNER_ARGS "" CACHE STRING
      "Arguments to pass to the parallel/serial test runner at the end of the command line.")

  # Echo runner configuration
  message(STATUS "Test runner: ${RUNNER}")
  message(STATUS "Test runner ncpus arg: ${RUNNER_NCPUS_ARG}")
  if(RUNNER_ARGS)
    message(STATUS "Test runner args: ${RUNNER_ARGS}")
  endif()
  if(POSTFIX_RUNNER_ARGS)
    message(STATUS "Test runner postfix args: ${POSTFIX_RUNNER_ARGS}")
  endif()
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

if (ENABLE_TESTS)
  # Query regression tests directory (will be exported to the source)
  get_filename_component(REGRESSION_DIR
    "${PROJECT_SOURCE_DIR}/../tests/regression" ABSOLUTE)
endif()

set(INCITER_EXECUTABLE inciter)
set(MESHCONV_EXECUTABLE meshconv)
if (ENABLE_TESTS)
  set(UNITTEST_EXECUTABLE unittest)
endif()

# Decide if the Exception class will write to std::cerr. This may not always be
# thread-safe, e.g., in Charm++'s SMP mode, so it can be disabled here. The
# default is set based on whether Charm++ is built in SMP or non-SMP mode.
if (NOT EXCEPTIONS_WRITE_TO_CERR)     # set default if not set
  if (SMP)
    set(cerr false)
  else()
    set(cerr true)
  endif()
  # will become compiler define in Main/QuinoaBuildConfig.hpp
  set(EXCEPTIONS_WRITE_TO_CERR ${cerr} CACHE BOOL "Exceptions write std::cerr.")
endif()
if(EXCEPTIONS_WRITE_TO_CERR)
  message(STATUS "Exceptions write to cerr")
else()
  message(STATUS "Exceptions do not write to cerr")
endif()

# Configure cmake variables to pass to the build
configure_file( "${PROJECT_SOURCE_DIR}/Main/XystConfig.hpp.in"
                "${PROJECT_BINARY_DIR}/Main/XystConfig.hpp" )
configure_file( "${PROJECT_SOURCE_DIR}/Main/XystBuildConfig.hpp.in"
                "${PROJECT_BINARY_DIR}/Main/XystBuildConfig.hpp" )
configure_file( "${PROJECT_SOURCE_DIR}/Main/XystConfig.cpp.in"
                "${PROJECT_BINARY_DIR}/Main/XystConfig.cpp" ESCAPE_QUOTES)

add_subdirectory(Base)
add_subdirectory(Control)
add_subdirectory(Partition)
add_subdirectory(Statistics)
add_subdirectory(Mesh)
add_subdirectory(IO)
add_subdirectory(Physics)
add_subdirectory(Inciter)
add_subdirectory(LinearSolver)
if (ENABLE_TESTS)
  add_subdirectory(UnitTest)
endif()
add_subdirectory(Main)

# Configure regression tests
if (ENABLE_TESTS)
  enable_testing()
  add_subdirectory(${REGRESSION_DIR}
                   ${CMAKE_BINARY_DIR}/tests)
endif()

# Configure building documentation and web page
if (BUILD_DOC)
  find_package(MCSS QUIET)
  if (DOXYGEN_FOUND AND MCSS_FOUND AND PYTHONINTERP_FOUND AND PYGMENTS_FOUND AND
      JINJA2_FOUND AND LATEX_FOUND)

    if (Git_FOUND)

      configure_file(${CMAKE_CURRENT_SOURCE_DIR}/../doc/xyst.doxy
                     ${CMAKE_CURRENT_BINARY_DIR}/doc/xyst.doxy @ONLY)

      file(GLOB doc_pages CONFIGURE_DEPENDS
           RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" "../doc/pages/*.dox")
      foreach(page ${doc_pages})
        configure_file(${page} ${CMAKE_CURRENT_BINARY_DIR}/doc/${page} @ONLY)
      endforeach()

      file(COPY ../doc/images DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/doc/html)

      add_custom_command(OUTPUT xml
                         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/doc
                         DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/doc/xyst.doxy
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
                                  DEPENDS ${MESHCONV_EXECUTABLE}
                                          ${INCITER_EXECUTABLE}
                                          ${UNITTEST_EXECUTABLE}
                                          numdiff exodiff)
  endif()
endif()

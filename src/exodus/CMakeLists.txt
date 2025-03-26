cmake_minimum_required(VERSION 3.18)

project(exodus C)
FILE(GLOB exodus_src src/ex_*.c)
add_library(exodus ${exodus_src})
set_property(TARGET exodus PROPERTY C_STANDARD 99)
target_include_directories(exodus PUBLIC include
                                         ${MPI_CXX_INCLUDE_DIRS}
                                         ${NETCDF_INCLUDE_DIRS})

set(suplib_cpp_src
    ${CMAKE_CURRENT_SOURCE_DIR}/suplib_cpp/FileInfo.C
    ${CMAKE_CURRENT_SOURCE_DIR}/suplib_cpp/GetLongOpt.C
    ${CMAKE_CURRENT_SOURCE_DIR}/suplib_cpp/SL_tokenize.C
    ${CMAKE_CURRENT_SOURCE_DIR}/suplib_cpp/copyright.C
    ${CMAKE_CURRENT_SOURCE_DIR}/suplib_cpp/format_time.C
    ${CMAKE_CURRENT_SOURCE_DIR}/suplib_cpp/hwm.C
    ${CMAKE_CURRENT_SOURCE_DIR}/suplib_cpp/iqsort.C
    ${CMAKE_CURRENT_SOURCE_DIR}/suplib_cpp/smart_assert.C
    ${CMAKE_CURRENT_SOURCE_DIR}/suplib_cpp/sys_info.C
    ${CMAKE_CURRENT_SOURCE_DIR}/suplib_cpp/time_stamp.C
    ${CMAKE_CURRENT_SOURCE_DIR}/suplib_cpp/timer.C)
add_library(suplib_cpp ${suplib_cpp_src})
target_include_directories(suplib_cpp PUBLIC suplib_cpp)
if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
  target_compile_options(suplib_cpp PRIVATE -Wno-documentation-unknown-command
                                            -Wno-signed-enum-bitfield
                                            -Wno-switch-enum
                                            -Wno-missing-prototypes
                                            -Wno-sign-conversion
                                            -Wno-undef
                                            -Wno-covered-switch-default
                                            -Wno-unused-template
                                            -Wno-missing-noreturn
                                            -Wno-extra-semi-stmt
                                            -Wno-old-style-cast)
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  target_compile_options(suplib_cpp PRIVATE -Wno-switch-enum
                                            -Wno-sign-conversion
                                            -Wno-undef
                                            -Wno-missing-noreturn
                                            -Wno-dangling-reference
                                            -Wno-old-style-cast)
endif()

set(suplib_c_src
    ${CMAKE_CURRENT_SOURCE_DIR}/suplib_c/add_to_log.c
    ${CMAKE_CURRENT_SOURCE_DIR}/suplib_c/adler.c
    ${CMAKE_CURRENT_SOURCE_DIR}/suplib_c/copy_string.c
    ${CMAKE_CURRENT_SOURCE_DIR}/suplib_c/open_file_limit.c
    ${CMAKE_CURRENT_SOURCE_DIR}/suplib_c/term_width.c)
add_library(suplib_c ${suplib_c_src})
target_include_directories(suplib_c PUBLIC suplib_c)

INSTALL(TARGETS exodus
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Runtime
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Runtime
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development)

project(exodiff CXX)
SET(exodiff_src ${CMAKE_CURRENT_SOURCE_DIR}/exodiff/ED_SystemInterface.C
	        ${CMAKE_CURRENT_SOURCE_DIR}/exodiff/exo_entity.C
	        ${CMAKE_CURRENT_SOURCE_DIR}/exodiff/exodiff.C
	        ${CMAKE_CURRENT_SOURCE_DIR}/exodiff/Tolerance.C
	        ${CMAKE_CURRENT_SOURCE_DIR}/exodiff/face_block.C
	        ${CMAKE_CURRENT_SOURCE_DIR}/exodiff/check.C
	        ${CMAKE_CURRENT_SOURCE_DIR}/exodiff/map.C
	        ${CMAKE_CURRENT_SOURCE_DIR}/exodiff/create_file.C
	        ${CMAKE_CURRENT_SOURCE_DIR}/exodiff/node_set.C
	        ${CMAKE_CURRENT_SOURCE_DIR}/exodiff/edge_block.C
	        ${CMAKE_CURRENT_SOURCE_DIR}/exodiff/side_set.C
	        ${CMAKE_CURRENT_SOURCE_DIR}/exodiff/exoII_read.C
	        ${CMAKE_CURRENT_SOURCE_DIR}/exodiff/stringx.C
	        ${CMAKE_CURRENT_SOURCE_DIR}/exodiff/exo_block.C
	        ${CMAKE_CURRENT_SOURCE_DIR}/exodiff/util.C)
add_executable(exodiff ${exodiff_src})
target_include_directories(exodiff PUBLIC exodiff
                                          suplib_c
                                          suplib_cpp
                                          include
                                          ${NETCDF_INCLUDE_DIRS})

if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
  target_compile_options(exodiff PRIVATE -Wno-switch-enum
                                         -Wno-undef
                                         -Wno-float-equal
                                         -Wno-undefined-func-template
                                         -Wno-missing-noreturn
                                         -Wno-unreachable-code-return
                                         -Wno-documentation
                                         -Wno-invalid-noreturn
                                         -Wno-dollar-in-identifier-extension
                                         -Wno-documentation-unknown-command
                                         -Wno-old-style-cast
                                         -Wno-weak-template-vtables
                                         -Wno-zero-as-null-pointer-constant
                                         -Wno-disabled-macro-expansion
                                         -Wno-missing-prototypes
                                         -Wno-signed-enum-bitfield
                                         -Wno-missing-variable-declarations
                                         -Wno-double-promotion
                                         -Wno-extra-semi-stmt
                                         -Wno-sign-conversion
                                         -Wno-shorten-64-to-32
                                         -Wno-implicit-float-conversion
                                         -Wno-signed-enum-bitfield)
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  target_compile_options(exodiff PRIVATE -Wno-switch-enum
                                         -Wno-undef
                                         -Wno-float-equal
                                         -Wno-missing-noreturn
                                         -Wno-old-style-cast
                                         -Wno-zero-as-null-pointer-constant
                                         -Wno-cast-function-type
                                         -Wno-double-promotion
                                         -Wno-sign-conversion)
endif()
target_link_libraries(exodiff PUBLIC exodus
                                     suplib_c
                                     suplib_cpp
                                     ${NETCDF_LIBRARIES}
                                     ${HDF5_LIBRARIES}
                                     ${MPI_CXX_LIBRARIES})
set(EXODIFF_EXECUTABLE ${CMAKE_CURRENT_BINARY_DIR}/exodiff PARENT_SCOPE)

project(numdiff C)

if(NOT "${CMAKE_C_FLAGS}" MATCHES "std=c99")
   if(NOT CMAKE_C_COMPILER_ID STREQUAL "Cray")
     if(NOT (WIN32 AND CMAKE_C_COMPILER_ID STREQUAL "Intel") AND
        NOT CMAKE_C_COMPILER_ID STREQUAL "PGI")
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99")
     else()
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /Qstd=c99")
     endif()
   endif()
endif()

file(GLOB numdiff_src src/*.c)

include_directories(src)
add_executable(numdiff ${numdiff_src})
set(NUMDIFF_EXECUTABLE ${CMAKE_CURRENT_BINARY_DIR}/numdiff PARENT_SCOPE)
if(UNIX)
   target_link_libraries(numdiff m)
endif()

if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
  target_compile_options(numdiff PRIVATE -Wno-fortify-source)
endif()

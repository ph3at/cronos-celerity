cmake_minimum_required(VERSION 3.22)
project(cronos-amr-sycl)

###################
# hipSYCL Options #
###################
find_package(hipSYCL CONFIG REQUIRED)

if(NOT HIPSYCL_DEBUG_LEVEL)
  if(CMAKE_BUILD_TYPE MATCHES "Debug")
    set(HIPSYCL_DEBUG_LEVEL "3" CACHE STRING
      "Choose the debug level, options are: 0 (no debug), 1 (print errors), 2 (also print warnings), 3 (also print general information)" FORCE)
  else()
    set(HIPSYCL_DEBUG_LEVEL "2" CACHE STRING
      "Choose the debug level, options are: 0 (no debug), 1 (print errors), 2 (also print warnings), 3 (also print general information)" FORCE)
  endif()
endif()

add_compile_definitions(HIPSYCL_DEBUG_LEVEL=${HIPSYCL_DEBUG_LEVEL})

####################
# Compiler Options #
####################
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

function(set_compiler_parameters TARGET)
  set_target_properties(${TARGET} PROPERTIES
    CXX_STANDARD 20
    CXX_EXTENSIONS OFF
  )
  target_compile_options(${TARGET} PRIVATE
    $<$<CXX_COMPILER_ID:MSVC>:/D_CRT_SECURE_NO_WARNINGS /MP /W4 /permissive- /bigobj /Zi /utf-8>
    $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-Wall -Wextra -pedantic -g -march=native>
    $<$<CXX_COMPILER_ID:GNU>:-fdiagnostics-color>
    $<$<CXX_COMPILER_ID:Clang,AppleClang>:-fcolor-diagnostics>
  )
  target_link_options(${TARGET} PRIVATE
    $<$<CXX_COMPILER_ID:MSVC>:/DEBUG /ignore:4099>
  )
endfunction()

##########
# OpenMP #
##########
option(USE_OpenMP "Use OpenMP" ON)
if(USE_OpenMP AND NOT MSVC)
  find_package(OpenMP)
  if(OPENMP_FOUND)
    add_compile_definitions(PARALLEL)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  endif()
endif()

########
# TOML #
########
set(CMAKE_FOLDER "vendor/tomlplusplus")
  add_subdirectory(vendor/tomlplusplus)
unset(CMAKE_FOLDER)

#########
# Build #
#########
file(GLOB_RECURSE SRC src/*.cpp)
list(REMOVE_ITEM SRC "${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp")
file(GLOB_RECURSE HDR src/*.h)

add_executable(cronos-amr-sycl src/main.cpp ${SRC} ${HDR})
add_sycl_to_target(TARGET cronos-amr-sycl SOURCES ${SRC})
set_compiler_parameters(cronos-amr-sycl)
target_link_libraries(cronos-amr-sycl PRIVATE tomlplusplus::tomlplusplus)

###########
# Testing #
###########
set(CMAKE_FOLDER "vendor/Catch2")
    add_subdirectory(vendor/Catch2)
    include("vendor/Catch2/extras/Catch.cmake")
unset(CMAKE_FOLDER)
add_executable(tests test/test.cpp ${SRC} ${HDR})
add_sycl_to_target(TARGET tests SOURCES ${SRC})
set_compiler_parameters(tests)
target_link_libraries(tests PRIVATE Catch2::Catch2WithMain tomlplusplus::tomlplusplus)

# Automatically register tests
include(CTest)
include(Catch)
catch_discover_tests(tests)

# Copy test data and configuration to build folder
file(COPY test/test-data DESTINATION .)
file(COPY configuration DESTINATION .)

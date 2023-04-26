cmake_minimum_required(VERSION 3.22)
project(cronos-amr-celerity)

####################
# Celerity Options #
####################
find_package(Celerity 0.3.2 REQUIRED)

####################
# Compiler Options #
####################
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

function(set_compiler_parameters TARGET)
  set_target_properties(${TARGET} PROPERTIES
    CXX_STANDARD 17
    CXX_EXTENSIONS OFF
  )
  target_compile_options(${TARGET} PRIVATE
    $<$<CXX_COMPILER_ID:MSVC>:/D_CRT_SECURE_NO_WARNINGS /MP /W4 /permissive- /bigobj /Zi /utf-8>
    $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-Wall -Wextra -pedantic -g -march=native $<$<CONFIG:Debug>:-O0>>
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
list(REMOVE_ITEM SRC "${CMAKE_CURRENT_SOURCE_DIR}/src/main-sycl.cpp")
list(REMOVE_ITEM SRC "${CMAKE_CURRENT_SOURCE_DIR}/src/main-celerity.cpp")
file(GLOB_RECURSE HDR src/*.h)

set(CRONOS_AMR_SRC ${SRC})
list(APPEND CRONOS_AMR_SRC "src/main.cpp")
add_executable(cronos-amr ${CRONOS_AMR_SRC} ${HDR})
set_compiler_parameters(cronos-amr)
target_link_libraries(cronos-amr PRIVATE tomlplusplus::tomlplusplus)
add_celerity_to_target(TARGET cronos-amr SOURCES ${CRONOS_AMR_SRC})

set(CRONOS_SYCL_SRC ${SRC})
list(APPEND CRONOS_SYCL_SRC "src/main-sycl.cpp")
add_executable(cronos-sycl ${CRONOS_SYCL_SRC} ${HDR})
set_compiler_parameters(cronos-sycl)
target_link_libraries(cronos-sycl PRIVATE tomlplusplus::tomlplusplus)
add_celerity_to_target(TARGET cronos-sycl SOURCES ${CRONOS_SYCL_SRC})

set(CRONOS_CELERITY_SRC ${SRC})
list(APPEND CRONOS_CELERITY_SRC "src/main-celerity.cpp")
add_executable(cronos-celerity ${CRONOS_CELERITY_SRC} ${HDR})
set_compiler_parameters(cronos-celerity)
target_link_libraries(cronos-celerity PRIVATE tomlplusplus::tomlplusplus)
add_celerity_to_target(TARGET cronos-celerity SOURCES ${CRONOS_CELERITY_SRC})

###########
# Testing #
###########
include(FetchContent)
FetchContent_Declare(Catch2
  URL ${PROJECT_SOURCE_DIR}/vendor/Catch2
)
FetchContent_MakeAvailable(Catch2)
list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/extras)

set(TESTS_SRC ${SRC})
list(APPEND TESTS_SRC "test/test.cpp")
add_executable(tests ${TESTS_SRC} ${HDR})
set_compiler_parameters(tests)
target_link_libraries(tests PRIVATE Catch2::Catch2 tomlplusplus::tomlplusplus)
add_celerity_to_target(TARGET tests SOURCES ${TESTS_SRC})

# Automatically register tests
include(CTest)
include(Catch)
catch_discover_tests(tests)

# Copy test data and configuration to build folder
file(COPY test/test-data DESTINATION .)
file(COPY configuration DESTINATION .)

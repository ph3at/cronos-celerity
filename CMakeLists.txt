cmake_minimum_required(VERSION 3.22)
project( cronos-amr )

####################
# Compiler Options #
####################
set( CMAKE_EXPORT_COMPILE_COMMANDS ON )
set( CMAKE_CPP_COMPILER clang )
add_compile_options( -Wall -Wextra -Wpedantic -march=native )
set( CMAKE_CXX_STANDARD 20 )
set( CMAKE_CXX_STANDARD_REQUIRED ON )
IF(CMAKE_BUILD_TYPE STREQUAL "Debug")
    add_compile_options( -O0 -pg )
ELSE()
    add_compile_options( -03 )
ENDIF()

##########
# OpenMP #
##########
OPTION (USE_OpenMP "Use OpenMP" ON)
IF(USE_OpenMP)
  FIND_PACKAGE(OpenMP)
  IF(OPENMP_FOUND)
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  ENDIF()
ENDIF()

#########
# Build #
#########
file(GLOB_RECURSE SRC src/*.cpp)
list(REMOVE_ITEM SRC "${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp")
file(GLOB_RECURSE HDR src/*.h)

add_executable(cronos-amr src/main.cpp ${SRC} ${HDR})

###########
# Testing #
###########
find_package(Catch2 2 REQUIRED)
add_executable(tests test/test.cpp ${SRC} ${HDR})
target_link_libraries(tests PRIVATE Catch2::Catch2)

# Automatically register tests
include(CTest)
include(Catch)
catch_discover_tests(tests)

# Copy test data to build folder
file(COPY test/test-data DESTINATION .)

cmake_minimum_required(VERSION 3.22)
project( cronos-amr )

####################
# Compiler Options #
####################
set( CMAKE_EXPORT_COMPILE_COMMANDS ON )
set( CMAKE_CPP_COMPILER clang )
add_compile_options( -Wall -Wextra -Wpedantic -O2 -march=native -std=c++20 )


#########
# Build #
#########
file(GLOB_RECURSE SRC src/*.cpp)
file(GLOB_RECURSE HDR src/*.h)

add_executable(cronos-amr ${SRC} ${HDR})

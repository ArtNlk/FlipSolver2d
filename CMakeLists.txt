cmake_minimum_required(VERSION 3.14)

project(FlipSolver2d)

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
#set(CMAKE_CXX_FLAGS "-O3 -Wall -Wextra")

set(USE_NUMPY_LOGGING ON)
if(USE_NUMPY_LOGGING)
    add_compile_definitions(NUMPY_LOGGING)
endif()

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    add_compile_definitions(DEBUG)
endif()

add_subdirectory(FlipSolver2dLib)
include_directories(FlipSolver2dLib)

add_subdirectory(Liquid2dRender)

cmake_minimum_required(VERSION 3.14)

project(FlipSolver2d)

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(MEMORY_SANITIZE ON)
set(ALL_WARNINGS OFF)
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
if(ALL_WARNINGS)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")
endif()
if(MEMORY_SANITIZE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address -fno-omit-frame-pointer")
    set(LDFLAGS "-fsanitize=address")
endif()

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

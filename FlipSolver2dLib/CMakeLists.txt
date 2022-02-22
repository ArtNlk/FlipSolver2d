cmake_minimum_required(VERSION 3.14)

project(FlipSolver2d LANGUAGES CXX)

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(HEADERS
    flipsolver2d.h
    linearindexable2d.h
    index2d.h
    grid2d.h
    fluidcell.h
    customassert.h
    fluidgrid.h
    simsettings.h
    dynamicsparsematrix.h
    sparsematrix.h
    pcgsolver.h
    vmath.h)

set(SOURCES
    flipsolver2d.cpp
    linearindexable2d.cpp
    index2d.cpp
    grid2d.cpp
    fluidcell.cpp
    fluidgrid.cpp
    simsettings.cpp
    dynamicsparsematrix.cpp
    sparsematrix.cpp
    pcgsolver.cpp
    vmath.cpp)

add_library(FlipSolver2d STATIC
  ${HEADERS}
  ${SOURCES}
)

target_compile_definitions(FlipSolver2d PRIVATE FLIPSOLVER2D_LIBRARY)

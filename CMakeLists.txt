cmake_minimum_required(VERSION 3.14)

project(FlipSolver2d LANGUAGES CXX)

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_library(FlipSolver2d STATIC
  flipsolver2d.cpp
  flipsolver2d.h
)

target_compile_definitions(FlipSolver2d PRIVATE FLIPSOLVER2D_LIBRARY)

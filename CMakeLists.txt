cmake_minimum_required(VERSION 3.14)

project(FlipSolver2d)

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_subdirectory(FlipSolver2dLib)
include_directories(FlipSolver2dLib)

add_subdirectory(Liquid2dRender)

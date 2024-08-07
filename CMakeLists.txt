cmake_minimum_required(VERSION 3.14)

project(FlipSolver2d)

if (CMAKE_VERSION VERSION_GREATER_EQUAL "3.24.0")
        cmake_policy(SET CMP0135 NEW)
endif()
set(FETCHCONTENT_QUIET FALSE)

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(MEMORY_SANITIZE OFF)
set(THREAD_SANITIZE OFF)
set(ALL_WARNINGS OFF)
set(POOL_SINGLE_THREAD OFF)
set(VECTORIZER_DEBUG OFF)
set(FAST_MATH ON)

if(POOL_SINGLE_THREAD)
    add_compile_definitions(ENQUEUE_IS_RUN)
endif()

set(FLUID_ACCEL "sse") #sse, avx, avx2, avx512

if(FLUID_ACCEL STREQUAL "sse")
    add_compile_definitions(FLUID_SSE)
    if(UNIX)
        add_compile_options(-msse4.2)
    elseif(MSVC)
        #no sse flag for msvc
    endif()
elseif(FLUID_ACCEL STREQUAL "avx")
    add_compile_definitions(FLUID_AVX)
    if(UNIX)
        add_compile_options(-mavx)
    elseif(MSVC)
        add_compile_options(/arch:AVX)
    endif()
elseif(FLUID_ACCEL STREQUAL "avx2")
    add_compile_definitions(FLUID_AVX2)
    if(UNIX)
        add_compile_options(-mavx2)
    elseif(MSVC)
        add_compile_options(/arch:AVX2)
    endif()
elseif(FLUID_ACCEL STREQUAL "avx512")
    #add_compile_definitions(FLUID_512)
    message(WARNING "No avx512 support yet")
endif()

if(UNIX)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pedantic-errors")
    if (CMAKE_BUILD_TYPE STREQUAL "Release" OR CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
        add_compile_options(-O3)
        if(FAST_MATH)
            add_compile_options(-ffast-math)
        endif()
    elseif(CMAKE_BUILD_TYPE STREQUAL "Debug")
        #add_compile_options(-Og)
    endif()

    if(ALL_WARNINGS)
        add_compile_options(-Wall -Wextra -Werror=return-type)
    endif()

    if(MEMORY_SANITIZE)
        add_compile_options(-fsanitize=address -fno-omit-frame-pointer)
        ADD_LINK_OPTIONS(-fsanitize=address -rdynamic)
    endif()

    if(THREAD_SANITIZE)
        add_compile_options(-fsanitize=thread -fno-omit-frame-pointer)
        ADD_LINK_OPTIONS(-fsanitize=thread -rdynamic)
    endif()
    add_compile_options(-fopenmp)

elseif(MSVC)
    if (CMAKE_BUILD_TYPE STREQUAL "Release" OR CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
        add_compile_options(/O2)
        add_compile_options(/std:c++20)
    endif()

    add_compile_options(/openmp)
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
add_subdirectory(bench)

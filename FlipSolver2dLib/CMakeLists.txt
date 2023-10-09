add_subdirectory("./threading")

set(HEADERS
    linearindexable2d.h
    index2d.h
    grid2d.h
    customassert.h
    materialgrid.h
    sdfgrid.h
    linearsolver.h
    vmath.h
    staticmatrix.h
    dynamicmatrix.h
    logger.h
    geometry2d.h
    mathfuncs.h
    emitter.h
    obstacle.h
    flipsolver2d.h
    flipsmokesolver.h
    sink.h
    staggeredvelocitygrid.h
    flipfiresolver.h
    nbflipsolver.h
    markerparticlesystem.h)

set(SOURCES
    linearindexable2d.cpp
    index2d.cpp
    grid2d.cpp
    materialgrid.cpp
    sdfgrid.cpp
    linearsolver.cpp
    vmath.cpp
    staticmatrix.cpp
    dynamicmatrix.cpp
    logger.cpp
    geometry2d.cpp
    mathfuncs.cpp
    emitter.cpp
    obstacle.cpp
    flipsolver2d.cpp
    flipsmokesolver.cpp
    sink.cpp
    staggeredvelocitygrid.cpp
    flipfiresolver.cpp
    nbflipsolver.cpp
    markerparticlesystem.cpp)

set(SSE42_HEADERS
    vops_sse42.h
    linearsolver_sse42.h
    grid_sse42.h)

set(SSE42_SOURCES
    vops_sse42.cpp
    linearsolver_sse42.cpp
    grid_sse42.cpp)

set(AVX2_HEADERS
    vops_avx2.h)

set(AVX2_SOURCES
    vops_avx2.cpp)

if(FLUID_ACCEL STREQUAL "sse")
    set(SOURCES ${SOURCES} ${SSE42_SOURCES})
    set(SOURCES ${SOURCES} ${SSE42_HEADERS})
elseif(FLUID_ACCEL STREQUAL "avx")
    message("No avx sources")
elseif(FLUID_ACCEL STREQUAL "avx2")
    set(SOURCES ${SOURCES} ${AVX2_SOURCES})
    set(SOURCES ${SOURCES} ${AVX2_HEADERS})
elseif(FLUID_ACCEL STREQUAL "avx512")
    message("No avx512 sources")
endif()

if(UNIX)

elseif(MSVC)

endif()

add_compile_options("$<$<CONFIG:DEBUG>:-DDEBUG>")

add_library(FlipSolver2d STATIC
  ${HEADERS}
  ${SOURCES}
  ${SSE42_HEADERS}
  ${SSE42_SOURCES}
  ${AVX_HEADERS}
  ${AVX_SOURCES}
)
target_include_directories(FlipSolver2d PUBLIC "./threading/")
target_link_libraries(FlipSolver2d SimThreadTools)

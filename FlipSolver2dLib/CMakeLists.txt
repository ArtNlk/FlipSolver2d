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
    uppertriangularmatrix.h
    dynamicuppertriangularsparsematrix.h
    logger.h
    geometry2d.h
    mathfuncs.h
    squarematrix.h
    emitter.h
    obstacle.h
    flipsolver2d.h
    flipsmokesolver.h
    sink.h
    staggeredvelocitygrid.h
    flipfiresolver.h
    nbflipsolver.h
    hwinfo.h
    markerparticlesystem.h)

set(SOURCES
    linearindexable2d.cpp
    index2d.cpp
    grid2d.cpp
    materialgrid.cpp
    sdfgrid.cpp
    linearsolver.cpp
    vmath.cpp
    uppertriangularmatrix.cpp
    dynamicuppertriangularsparsematrix.cpp
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
    hwinfo.cpp
    markerparticlesystem.cpp)

set(SSE42_HEADERS
    vops_sse42.h
    linearsolver_sse42.h
    grid_sse42.h)

set(SSE42_SOURCES
    vops_sse42.cpp
    linearsolver_sse42.cpp
    grid_sse42.cpp)

if(UNIX)
    set_source_files_properties(${SSE42_SOURCES} PROPERTIES COMPILE_FLAGS "-mmmx -msse -msse2 -msse3 -mssse3 -msse4 -msse4.1 -msse4.2")
elseif(MSVC)

endif()

add_compile_options("$<$<CONFIG:DEBUG>:-DDEBUG>")

add_library(FlipSolver2d STATIC
  ${HEADERS}
  ${SOURCES}
  ${SSE42_HEADERS}
  ${SSE42_SOURCES}
)
target_include_directories(FlipSolver2d PUBLIC "./threading/")
target_link_libraries(FlipSolver2d SimThreadTools)

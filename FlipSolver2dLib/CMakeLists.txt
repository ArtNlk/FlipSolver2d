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
    nbflipsolver.h)

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
    nbflipsolver.cpp)

add_compile_options("$<$<CONFIG:DEBUG>:-DDEBUG>")

add_library(FlipSolver2d STATIC
  ${HEADERS}
  ${SOURCES}
)
target_include_directories(FlipSolver2d PUBLIC "./threading/")
target_link_libraries(FlipSolver2d SimThreadTools)

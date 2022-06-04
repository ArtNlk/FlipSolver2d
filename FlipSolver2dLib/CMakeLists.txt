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
    vmath.h
    uppertriangularmatrix.h
    dynamicuppertriangularsparsematrix.h
    logger.h
    geometry2d.h
    mathfuncs.h)

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
    vmath.cpp
    uppertriangularmatrix.cpp
    dynamicuppertriangularsparsematrix.cpp
    logger.cpp
    geometry2d.cpp
    mathfuncs.cpp)

add_compile_options("$<$<CONFIG:DEBUG>:-DDEBUG>")

add_library(FlipSolver2d STATIC
  ${HEADERS}
  ${SOURCES}
)

set(HEADERS
    staticmatrix.h
    dynamicmatrix.h
    matrixweights.h)

set(SOURCES
    staticmatrix.cpp
    dynamicmatrix.cpp
    matrixweights.cpp)

add_compile_options("$<$<CONFIG:DEBUG>:-DDEBUG>")

add_library(SimMatrixTools STATIC
  ${HEADERS}
  ${SOURCES}
  matrixbase.h
)

target_include_directories(SimMatrixTools PRIVATE "../")
target_include_directories(SimMatrixTools PRIVATE "../threading")
target_link_libraries(SimMatrixTools PUBLIC SimThreadTools)

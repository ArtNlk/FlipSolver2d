set(HEADERS
   threadpool.h)

set(SOURCES
    threadpool.cpp)

add_compile_options("$<$<CONFIG:DEBUG>:-DDEBUG>")

add_library(SimThreadTools STATIC
  ${HEADERS}
  ${SOURCES}
)

find_package(Threads REQUIRED)
target_link_libraries(SimThreadTools INTERFACE Threads::Threads)

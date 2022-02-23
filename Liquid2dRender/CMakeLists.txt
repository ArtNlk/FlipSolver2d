include_directories(${GLFW3_INCLUDE_DIR})
include_directories(include)
include_directories(include/glad)
include_directories(include/KHR)

set(HEADERS
    include/glad/glad.h
    liquidrenderapp.h)

set(SOURCES
    src/glad.c
    main.cpp
    liquidrenderapp.cpp)

find_package(OpenGL REQUIRED)
find_package(glfw3 REQUIRED)

add_executable(LiquidRenderApp ${HEADERS} ${SOURCES})

target_link_libraries(LiquidRenderApp FlipSolver2d)
target_link_libraries(LiquidRenderApp ${OPENGL_LIBRARIES} ${GLFW3_LIBRARY})

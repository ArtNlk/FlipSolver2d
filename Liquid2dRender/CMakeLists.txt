include_directories(${GLFW3_INCLUDE_DIR})
include_directories(include)
include_directories(include/glad)
include_directories(include/KHR)
include_directories(include/external)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(HEADERS
    include/glad/glad.h
    liquidrenderapp.h
    fluidrenderer.h
    texture2d.h
    color.h
    vertex.h)

set(SOURCES
    src/glad.c
    main.cpp
    liquidrenderapp.cpp
    fluidrenderer.cpp
    texture2d.cpp
    color.cpp
    vertex.cpp)

find_package(OpenGL REQUIRED)
find_package(glfw3 REQUIRED)
find_package(Freetype REQUIRED)

add_executable(LiquidRenderApp ${HEADERS} ${SOURCES})

target_link_libraries(LiquidRenderApp FlipSolver2d)
target_link_libraries(LiquidRenderApp ${OPENGL_LIBRARIES} glfw)
target_link_libraries(LiquidRenderApp ${FREETYPE_LIBRARIES})
if(UNIX)
target_link_libraries(LiquidRenderApp ${CMAKE_DL_LIBS})
endif()
target_include_directories(LiquidRenderApp PRIVATE ${FREETYPE_INCLUDE_DIRS})
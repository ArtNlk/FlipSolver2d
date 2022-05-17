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
    vertex.h
    textrenderer.h
    textmenurenderer.h)

set(SOURCES
    src/glad.c
    main.cpp
    liquidrenderapp.cpp
    fluidrenderer.cpp
    texture2d.cpp
    color.cpp
    vertex.cpp
    textrenderer.cpp
    textmenurenderer.cpp)

find_package(OpenGL REQUIRED)
find_package(glfw3 REQUIRED)
find_package(Freetype REQUIRED)
find_package(glm REQUIRED)

if(NOT DEFINED GLFW3_LIBRARY)
    set(GLFW3_LIBRARY glfw)
endif()

add_executable(LiquidRenderApp ${HEADERS} ${SOURCES})

target_link_libraries(LiquidRenderApp FlipSolver2d)
target_link_libraries(LiquidRenderApp ${OPENGL_LIBRARIES} ${GLFW3_LIBRARY})
target_link_libraries(LiquidRenderApp ${FREETYPE_LIBRARIES})
if(UNIX)
target_link_libraries(LiquidRenderApp ${CMAKE_DL_LIBS})
endif()
target_include_directories(LiquidRenderApp PRIVATE ${GLM_INCLUDE_DIRS})
target_link_directories(LiquidRenderApp PRIVATE ${GLM_LIBRARY_DIRS})
add_definitions(${GLM_DEFINITIONS})
target_include_directories(LiquidRenderApp PRIVATE ${FREETYPE_INCLUDE_DIRS})

include_directories(${GLFW3_INCLUDE_DIR})
include_directories(include)
include_directories(include/glad)
include_directories(include/KHR)
include_directories(include/external)

include (FindPkgConfig)
if (NOT PKG_CONFIG_FOUND)
    message (FATAL_ERROR "pkg-config not found")
endif ()

find_package(OpenGL REQUIRED)
find_package(glfw3 REQUIRED)
find_package(Freetype REQUIRED)
find_package(glm REQUIRED)
find_package(ZLIB REQUIRED)
include_directories(${ZLIB_INCLUDE_DIR})
#find_package(libpng REQUIRED)

if(NOT DEFINED GLFW3_LIBRARY)
    set(GLFW3_LIBRARY glfw)
endif()

set(HEADERS
    include/glad/glad.h
    liquidrenderapp.h
    fluidrenderer.h
    texture2d.h
    color.h
    textrenderer.h
    textmenurenderer.h
    globalcallbackhandler.h)

set(SOURCES
    src/glad.c
    main.cpp
    liquidrenderapp.cpp
    fluidrenderer.cpp
    texture2d.cpp
    color.cpp
    textrenderer.cpp
    textmenurenderer.cpp
    globalcallbackhandler.cpp)

add_executable(LiquidRenderApp ${HEADERS} ${SOURCES})

## check for libpng
pkg_check_modules(LIBPNG libpng16 REQUIRED)
if (NOT LIBPNG_FOUND)
    message(FATAL_ERROR "You don't seem to have libpng libraries installed")
else ()
    include_directories (${LIBPNG_INCLUDE_DIRS})
    target_link_directories (LiquidRenderApp PRIVATE ${LIBPNG_LIBRARY_DIRS})
    target_link_libraries (LiquidRenderApp ${LIBPNG_LIBRARIES} ${ZLIB_LIBRARY})
endif ()
message(LIBS ${GLM_LIBRARY_DIRS})
message(INCLUDES ${LIBPNG_INCLUDE_DIRS})

target_link_libraries(LiquidRenderApp FlipSolver2d)
target_link_libraries(LiquidRenderApp ${OPENGL_LIBRARIES} ${GLFW3_LIBRARY} -static-libgcc -static-libstdc++)
target_link_libraries(LiquidRenderApp ${FREETYPE_LIBRARIES})
if(UNIX)
target_link_libraries(LiquidRenderApp ${CMAKE_DL_LIBS})
endif()
target_include_directories(LiquidRenderApp PRIVATE ${GLM_INCLUDE_DIRS})
target_link_directories(LiquidRenderApp PRIVATE ${GLM_LIBRARY_DIRS})
add_definitions(${GLM_DEFINITIONS})
target_include_directories(LiquidRenderApp PRIVATE ${FREETYPE_INCLUDE_DIRS})

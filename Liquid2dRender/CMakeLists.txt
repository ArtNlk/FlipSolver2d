#include_directories(${GLFW3_INCLUDE_DIR})
include_directories(include)
include_directories(include/glad)
include_directories(include/KHR)
include_directories(include/external)

include(FetchContent)

find_package(OpenGL REQUIRED)

#GLFW
FetchContent_Declare(
        glfw
        GIT_REPOSITORY https://github.com/glfw/glfw
        GIT_TAG 3.4
        GIT_SHALLOW TRUE
        GIT_PROGRESS TRUE
)
FetchContent_GetProperties(glfw)
if(NOT glfw_POPULATED)
    FetchContent_Populate(glfw)

    set(GLFW_BUILD_EXAMPLES OFF CACHE INTERNAL "Build the GLFW example programs")
    set(GLFW_BUILD_TESTS OFF CACHE INTERNAL "Build the GLFW test programs")
    set(GLFW_BUILD_DOCS OFF CACHE INTERNAL "Build the GLFW documentation")
    set(GLFW_INSTALL OFF CACHE INTERNAL "Generate installation target")
    if(LINUX_WAYLAND)
        set(GLFW_BUILD_WAYLAND ON)
        set(GLFW_BUILD_X11 OFF)
    else()
        set(GLFW_BUILD_WAYLAND OFF)
        set(GLFW_BUILD_X11 ON)
    endif()
    add_subdirectory(${glfw_SOURCE_DIR} ${glfw_BINARY_DIR})
endif()

add_library(imgui_glfw STATIC
    ${CMAKE_SOURCE_DIR}/imgui/imgui.cpp
    ${CMAKE_SOURCE_DIR}/imgui/imgui_draw.cpp
    ${CMAKE_SOURCE_DIR}/imgui/imgui_demo.cpp
    ${CMAKE_SOURCE_DIR}/imgui/imgui_tables.cpp
    ${CMAKE_SOURCE_DIR}/imgui/imgui_widgets.cpp

    ${CMAKE_SOURCE_DIR}/imgui/backends/imgui_impl_glfw.cpp
    ${CMAKE_SOURCE_DIR}/imgui/backends/imgui_impl_opengl3.cpp
)

target_link_libraries(imgui_glfw PUBLIC glfw ${OPENGL_LIBRARIES})

target_include_directories(imgui_glfw
PUBLIC
    ${CMAKE_SOURCE_DIR}/imgui
    ${CMAKE_SOURCE_DIR}/imgui/backends
)

if(VECTORIZER_DEBUG AND UNIX)
    add_compile_options(-fopt-info-vec-missed)
elseif(VECTORIZER_DEBUG AND MSVC)
    add_compile_options(/Qvec-report:2)
endif()

#GLM
FetchContent_Declare(
        glm
        GIT_REPOSITORY https://github.com/g-truc/glm.git
        GIT_TAG bf71a834948186f4097caa076cd2663c69a10e1e
        GIT_SHALLOW TRUE
        GIT_PROGRESS TRUE
)
FetchContent_GetProperties(glm)
if(NOT glm_POPULATED)
    FetchContent_Populate(glm)
    add_subdirectory(${glm_SOURCE_DIR} ${glm_BINARY_DIR})
endif()

set(HEADERS
    include/glad/glad.h
    liquidrenderapp.h
    fluidrenderer.h
    texture2d.h
    color.h)

set(SOURCES
    src/glad.c
    main.cpp
    liquidrenderapp.cpp
    fluidrenderer.cpp
    texture2d.cpp
    color.cpp)

file(GLOB SCENES
     "scenes/*.json"
)

set(CONFIG_FILE "./config.json")

add_executable(LiquidRenderApp ${HEADERS} ${SOURCES} ${SCENES} ${CONFIG_FILE})

get_filename_component(filename ${CONFIG_FILE} NAME)
file(READ "${CONFIG_FILE}" ignored)
unset(ignored) # Drop the memory
configure_file(${CONFIG_FILE} ${CMAKE_BINARY_DIR}/Liquid2dRender/${filename} COPYONLY)

foreach(file ${SCENES})
    get_filename_component(filename ${file} NAME)
    file(READ "${file}" ignored)
    unset(ignored) # Drop the memory
    configure_file(${file} ${CMAKE_BINARY_DIR}/Liquid2dRender/scenes/${filename} COPYONLY)
endforeach()

target_link_libraries(LiquidRenderApp FlipSolver2d)
target_link_libraries(LiquidRenderApp ${OPENGL_LIBRARIES} -static-libgcc -static-libstdc++)
target_link_libraries(LiquidRenderApp glfw)
target_link_libraries(LiquidRenderApp glm)
target_link_libraries(LiquidRenderApp Utils)
target_link_libraries(LiquidRenderApp imgui_glfw)

if(UNIX)
target_link_libraries(LiquidRenderApp ${CMAKE_DL_LIBS})
endif()

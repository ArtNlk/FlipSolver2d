#include_directories(${GLFW3_INCLUDE_DIR})
include_directories(include)
include_directories(include/glad)
include_directories(include/KHR)
include_directories(include/external)

include(FetchContent)

#GLFW
FetchContent_Declare(
        glfw
        GIT_REPOSITORY https://github.com/glfw/glfw
        GIT_TAG 7482de6071d21db77a7236155da44c172a7f6c9e
)
FetchContent_GetProperties(glfw)
if(NOT glfw_POPULATED)
    FetchContent_Populate(glfw)

    set(GLFW_BUILD_EXAMPLES OFF CACHE INTERNAL "Build the GLFW example programs")
    set(GLFW_BUILD_TESTS OFF CACHE INTERNAL "Build the GLFW test programs")
    set(GLFW_BUILD_DOCS OFF CACHE INTERNAL "Build the GLFW documentation")
    set(GLFW_INSTALL OFF CACHE INTERNAL "Generate installation target")

    add_subdirectory(${glfw_SOURCE_DIR} ${glfw_BINARY_DIR})
endif()

#FREETYPE
FetchContent_Declare(
  freetype
  URL https://download.savannah.gnu.org/releases/freetype/freetype-2.12.1.tar.gz
)
FetchContent_GetProperties(freetype)
if(NOT freetype_POPULATED)
    FetchContent_Populate(freetype)

    set(FT_DISABLE_ZLIB ON CACHE INTERNAL "Use internal ZLIB")
    set(FT_DISABLE_BZIP2 ON CACHE INTERNAL "Disable BZIP2 compressed fonts")
    set(FT_DISABLE_PNG ON CACHE INTERNAL "Disable PNG fonts")
    set(FT_DISABLE_HARFBUZZ ON CACHE INTERNAL "Disable Harf buzz autohint")
    set(FT_DISABLE_BROTLI ON CACHE INTERNAL "Disable WOFF2 fonts")

    add_subdirectory(${freetype_SOURCE_DIR} ${freetype_BINARY_DIR})
endif()

#GLM
FetchContent_Declare(
        glm
        GIT_REPOSITORY https://github.com/g-truc/glm.git
        GIT_TAG bf71a834948186f4097caa076cd2663c69a10e1e
)
FetchContent_GetProperties(glm)
if(NOT glm_POPULATED)
    FetchContent_Populate(glm)
    add_subdirectory(${glm_SOURCE_DIR} ${glm_BINARY_DIR})
endif()

#NLOHMANN_JSON
FetchContent_Declare(
        nlohmann_json
        GIT_REPOSITORY https://github.com/nlohmann/json.git
        GIT_TAG bc889afb4c5bf1c0d8ee29ef35eaaf4c8bef8a5d
)
FetchContent_GetProperties(nlohmann_json)
if(NOT nlohmann_json_POPULATED)
    FetchContent_Populate(nlohmann_json)
    add_subdirectory(${nlohmann_json_SOURCE_DIR} ${nlohmann_json_BINARY_DIR})
endif()

#ZLIB
FetchContent_Declare(
        zlib
        GIT_REPOSITORY https://github.com/madler/zlib.git
        GIT_TAG 04f42ceca40f73e2978b50e93806c2a18c1281fc
)
FetchContent_GetProperties(zlib)
if(NOT zlib_POPULATED)
    FetchContent_Populate(zlib)
    add_subdirectory(${zlib_SOURCE_DIR} ${zlib_BINARY_DIR})
endif()
set(ZLIB_LIBRARY zlib)

#LIBPNG
FetchContent_Declare(
        libpng
        GIT_REPOSITORY https://github.com/glennrp/libpng.git
        GIT_TAG 07b8803110da160b158ebfef872627da6c85cbdf
)
FetchContent_GetProperties(libpng)
if(NOT libpng_POPULATED)
    FetchContent_Populate(libpng)

    set(PNG_BUILD_ZLIB ON CACHE INTERNAL "Custom zlib location")
    set(PNG_STATIC ON CACHE INTERNAL "Build static lib")
    set(PNG_SHARED OFF CACHE INTERNAL "Build shared lib")
    set(PNG_EXECUTABLES OFF CACHE INTERNAL "Build exes")
    set(PNG_TESTS OFF CACHE INTERNAL "Build tests")

    add_subdirectory(${libpng_SOURCE_DIR} ${libpng_BINARY_DIR})
endif()

find_package(OpenGL REQUIRED)

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
target_link_libraries(LiquidRenderApp freetype)
target_link_libraries(LiquidRenderApp libpng16.a)
target_link_libraries(LiquidRenderApp zlib)
target_link_libraries(LiquidRenderApp nlohmann_json::nlohmann_json)

if(UNIX)
target_link_libraries(LiquidRenderApp ${CMAKE_DL_LIBS})
endif()

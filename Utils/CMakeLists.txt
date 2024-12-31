#NLOHMANN_JSON
FetchContent_Declare(
        nlohmann_json
        GIT_REPOSITORY https://github.com/nlohmann/json.git
        GIT_TAG bc889afb4c5bf1c0d8ee29ef35eaaf4c8bef8a5d
        GIT_SHALLOW TRUE
        GIT_PROGRESS TRUE
)
FetchContent_GetProperties(nlohmann_json)
if(NOT nlohmann_json_POPULATED)
    FetchContent_Populate(nlohmann_json)
    add_subdirectory(${nlohmann_json_SOURCE_DIR} ${nlohmann_json_BINARY_DIR})
endif()

set(HEADERS
    jsonscenereader.h)

set(SOURCES
    jsonscenereader.cpp)

add_compile_options("$<$<CONFIG:DEBUG>:-DDEBUG>")

add_library(Utils STATIC
  ${HEADERS}
  ${SOURCES}
)

target_link_libraries(Utils FlipSolver2d)
target_link_libraries(Utils nlohmann_json::nlohmann_json)

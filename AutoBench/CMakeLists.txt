project(FlipSolver2d)

#OpenXLSX
FetchContent_Declare(
        openxlsx
        GIT_REPOSITORY https://github.com/troldal/OpenXLSX.git
        GIT_TAG 8df551c1c62d77f6cd7f26cd67b725eec0e303ed
        GIT_SHALLOW TRUE
        GIT_PROGRESS TRUE
        CMAKE_ARGS
          -DOPENXLSX_CREATE_DOCS=OFF
          -DOPENXLSX_BUILD_SAMPLES=OFF
          -DOPENXLSX_BUILD_TESTS=OFF
)
FetchContent_GetProperties(openxlsx)
if(NOT openxlsx_POPULATED)
    FetchContent_Populate(openxlsx)
    add_subdirectory(${openxlsx_SOURCE_DIR} ${openxlsx_BINARY_DIR})
endif()

set(SOURCES
    main.cpp
    benchmarkrunnerapplication.cpp
    argumentparser.cpp
    benchruntable.cpp)

set(HEADERS
    benchmarkrunnerapplication.h
    argumentparser.h
    benchruntable.h)

add_executable(Autobench ${HEADERS} ${SOURCES})

target_link_libraries(Autobench FlipSolver2d)
target_link_libraries(Autobench Utils)
target_link_libraries(Autobench OpenXLSX::OpenXLSX)

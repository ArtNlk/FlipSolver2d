Include(FetchContent)

FetchContent_Declare(
  Catch2
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG        v3.8.0 # or a later release
)

FetchContent_MakeAvailable(Catch2)

set(SOURCES
    AlwaysPass.cpp
    LightViscosityMatrixTest.cpp)

set(HEADERS)

add_executable(UnitTests ${SOURCES} ${HEADERS})

target_link_libraries(UnitTests Catch2::Catch2WithMain)
target_link_libraries(UnitTests FlipSolver2d)
target_link_libraries(UnitTests Eigen3::Eigen)
target_include_directories(UnitTests PUBLIC ${Eigen3_SOURCE_DIR})

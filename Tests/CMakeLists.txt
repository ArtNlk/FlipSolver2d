Include(FetchContent)

FetchContent_Declare(
  Catch2
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG        v3.8.0 # or a later release
)

FetchContent_MakeAvailable(Catch2)

set(SOURCES
    AlwaysPass.cpp)

set(HEADERS)

add_executable(UnitTests ${SOURCES} ${HEADERS})

target_link_libraries(UnitTests PRIVATE Catch2::Catch2WithMain)

project(FlipSolver2d)

find_package(OpenMP REQUIRED)

add_executable(Bench main.cpp)
target_link_libraries(Bench PRIVATE Eigen3::Eigen)
target_link_libraries(Bench PRIVATE OpenMP::OpenMP_CXX)

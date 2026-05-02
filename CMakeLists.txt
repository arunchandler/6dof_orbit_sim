cmake_minimum_required(VERSION 3.10)
project(6dof_orbit_sim CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(Eigen3 REQUIRED NO_MODULE)

add_library(6dof_orbit_sim
    src/state_conversions.cpp
    src/elements.cpp
    src/translational_dynamics.cpp
    src/attitude_dynamics.cpp
)
target_include_directories(6dof_orbit_sim PUBLIC include)
target_link_libraries(6dof_orbit_sim PUBLIC Eigen3::Eigen)

add_subdirectory(tests)
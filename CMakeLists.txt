cmake_minimum_required (VERSION 3.9)

project (MPI)

find_package(MPI REQUIRED)

set(CMAKE_C_STANDARD 99)
link_libraries(m)

add_executable(monte_carlo monte_carlo.c)
target_link_libraries(monte_carlo PUBLIC MPI::MPI_C)

cmake_minimum_required(VERSION 3.10)
project(Bayesian_Phylo)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_subdirectory(${PROJECT_SOURCE_DIR}/lib/ncl)

include_directories(${PROJECT_SOURCE_DIR}/src)

file(GLOB_RECURSE SOURCES "${PROJECT_SOURCE_DIR}/src/*.cpp")

add_executable(main ${SOURCES})

target_include_directories(main PRIVATE ${PROJECT_SOURCE_DIR}/src)

target_link_libraries(main PRIVATE ncl_shared)
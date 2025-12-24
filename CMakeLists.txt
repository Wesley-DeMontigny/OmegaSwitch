cmake_minimum_required(VERSION 3.20)
project(Bayesian_Phylo)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

option(TEST_RUN "Enabled Test Environment Run Parameters" ON)
option(SAMPLE_PRIOR "Enabled Prior Sampling Mode" OFF)
option(DEBUG "Enabled Debug Mode" OFF)
option(TIME_PROFILE "Enabled Time Profiling" OFF)

add_subdirectory(${PROJECT_SOURCE_DIR}/lib/ncl)

include_directories(/workspaces/taskflow/) # Add the location of taskflow

include_directories(${PROJECT_SOURCE_DIR}/src)

file(GLOB_RECURSE SOURCES "${PROJECT_SOURCE_DIR}/src/*.cpp")

add_executable(main ${SOURCES})

target_link_libraries(main PRIVATE ncl_static)
target_link_libraries(main PRIVATE pthread) # This was required on RedHat in my experience
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -O3")

if(DEBUG)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=leak -fno-omit-frame-pointer -g")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address -fno-omit-frame-pointer -g")
endif()

if(TIME_PROFILE)
add_compile_definitions(TIME_PROFILE=1)
endif()

if(SAMPLE_PRIOR)
add_compile_definitions(SAMPLE_PRIOR=1)
endif()

if(TEST_RUN)
add_compile_definitions(TEST_RUN=1)
endif()

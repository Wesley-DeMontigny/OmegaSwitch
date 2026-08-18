cmake_minimum_required(VERSION 3.20)
project(OmegaSwitch)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

option(TEST_RUN "Enabled Test Environment Run Parameters" OFF)
option(SAMPLE_PRIOR "Enabled Prior Sampling Mode" OFF)
option(DEBUG "Enabled Debug Mode" OFF)
option(TIME_PROFILE "Enabled Time Profiling" OFF)

add_subdirectory(${PROJECT_SOURCE_DIR}/lib/ncl)


find_package(Taskflow QUIET)

if(NOT Taskflow_FOUND)
    include(FetchContent)
    FetchContent_Declare(
        Taskflow
        GIT_REPOSITORY https://github.com/taskflow/taskflow.git
        GIT_TAG        v3.7.0
    )
    FetchContent_MakeAvailable(Taskflow)
endif()

include_directories(${PROJECT_SOURCE_DIR}/src)

file(GLOB_RECURSE SOURCES "${PROJECT_SOURCE_DIR}/src/*.cpp")

add_executable(omega_switch ${SOURCES})

target_link_libraries(omega_switch PRIVATE ncl_static)
target_link_libraries(omega_switch PRIVATE pthread) # This was required on RedHat in my experience
target_link_libraries(omega_switch PRIVATE Taskflow)

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

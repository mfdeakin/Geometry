
cmake_minimum_required (VERSION 2.8)
project(geometry)

find_package(CUDA QUIET REQUIRED)

list(APPEND CUDA_NVCC_FLAGS "-std=c++11")

file(GLOB TIMER_SRC "src/timer/*.cpp")
file(GLOB TEST_SRC "src/tests/*.cpp")
file(GLOB TIMING_INTERSECTIONS_SRC "src/timing_intersections/*.cpp")
file(GLOB CUDATEST_SRC "src/cudatest/*.cu")

file(GLOB TIMER_INC "include/timer" "*.h" "*.hpp")
file(GLOB GEOMETRY_INC "include/geometry" "*.h" "*.hpp")

include_directories(./include/geometry ./include/timer)

set(CXX_COMPILE_FLAGS "-std=c++11")

set(CMAKE_C_FLAGS "${CXX_COMPILE_FLAGS}")

set(CMAKE_CXX_FLAGS "${CXX_COMPILE_FLAGS}")

set(CUDA_64_BIT_DEVICE_CODE ON)

cuda_add_executable(cudatest ${CUDATEST_SRC} ${TIMER_SRC})
add_executable(timing_intersections ${TIMING_INTERSECTIONS_SRC} ${TIMER_SRC})
add_executable(geomtest ${TEST_SRC})

target_link_libraries(geomtest mpfr gtest pthread)
target_link_libraries(timing_intersections mpfr)
target_link_libraries(cudatest mpfr)


cmake_minimum_required (VERSION 2.8)
project(geometry)

find_package(CUDA QUIET REQUIRED)

file(GLOB TIMER_SRC "src/timer/*.cpp")
file(GLOB TEST_SRC "src/tests/*.cpp")
file(GLOB TIMING_INTERSECTIONS_SRC "src/timing_intersections/*.cpp")
file(GLOB CUDATEST_SRC "src/cudatest/*.cu")

file(GLOB TIMER_INC "include/timer" "*.h" "*.hpp")
file(GLOB GEOMETRY_INC "include/geometry" "*.h" "*.hpp")

include_directories(./include/geometry ./include/timer)

set(CUDA_64_BIT_DEVICE_CODE ON)
list(APPEND CUDA_NVCC_FLAGS "-std=c++11")

add_library(timer STATIC ${TIMER_SRC})
cuda_add_executable(cudatest ${CUDATEST_SRC})
add_executable(timing_intersections ${TIMING_INTERSECTIONS_SRC})
add_executable(geomtest ${TEST_SRC})

set_target_properties(timer timing_intersections geomtest PROPERTIES COMPILE_FLAGS "-std=c++11")

target_link_libraries(geomtest mpfr gtest pthread)
target_link_libraries(timing_intersections mpfr timer)
target_link_libraries(cudatest mpfr timer)


cmake_minimum_required (VERSION 2.6)

project (Geometry)

set(TIMER_SRC "src/timer/timer.cpp")
file(GLOB TEST_SRC "src/tests/*.cpp")
file(GLOB TIMING_INTERSECTIONS_SRC "src/timing_intersections/*.cpp")
file(GLOB CUDATEST_SRC "src/cudatest/*.cu")

set(CXX_COMPILE_FLAGS "-std=c++11 -I./include/geometry -I./include/timer -g")

set(CMAKE_CXX_FLAGS "${CXX_COMPILE_FLAGS}")

set(LD_LIBRARIES "-lmpfr -lgtest")
set(CMAKE_EXE_LINKER_FLAGS "${LD_LIBRARIES}")
add_executable(geomtest ${TEST_SRC})

add_executable(sorttime ${TIMING_INTERSECTIONS_SRC} ${TIMER_SRC})

set(CMAKE_CXX_COMPILER "nvcc")
add_executable(cudatest ${CUDATEST_SRC} ${TIMER_SRC})

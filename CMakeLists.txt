
cmake_minimum_required (VERSION 2.8)
project(geometry)

find_package(CUDA QUIET)

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
#cuda_add_executable(cudatest ${CUDATEST_SRC})
add_executable(timing_intersections ${TIMING_INTERSECTIONS_SRC})
add_executable(geomtest ${TEST_SRC})

set_target_properties(timer timing_intersections geomtest PROPERTIES COMPILE_FLAGS "-std=c++11")

target_link_libraries(geomtest mpfr gtest pthread)
target_link_libraries(timing_intersections mpfr timer)
if(CUDA_FOUND)
        execute_process(COMMAND ${CMAKE_C_COMPILER} -dumpversion OUTPUT_VARIABLE GCC_VERSION)
        string(REGEX MATCHALL "[0-9]+" GCC_VERSION_COMPONENT ${GCC_VERSION})
        list(GET GCC_VERSION_COMPONENTS 0 GCC_MAJOR)
        list(GET GCC_VERSION_COMPONENTS 1 GCC_MINOR)
        if(GCC_VERSION VERSION_LESS 5.3 OR GCC_VERSION VERSION_GREATER 5.3)
                target_link_libraries(cudatest mpfr timer)
        endif()
endif()

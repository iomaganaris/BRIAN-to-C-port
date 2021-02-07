project(simulation)
cmake_minimum_required(VERSION 3.7)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories(src)

set(COMMON_SOURCES src/constants.cpp src/simulation.cpp)

set(ORIG_SOURCES src/neuron.cpp src/synapses.cpp)

set(OPTIMIZED_SOURCES src/neuron.cpp src/synapses_opt.cpp)

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")

    set(OMP_SOURCES src/omp/neuron.cpp src/omp/synapses.cpp)
    add_executable(${PROJECT_NAME}_omp ${COMMON_SOURCES} ${OMP_SOURCES})
endif()

add_executable(${PROJECT_NAME} ${COMMON_SOURCES} ${ORIG_SOURCES})

add_executable(${PROJECT_NAME}_opt ${COMMON_SOURCES} ${OPTIMIZED_SOURCES})


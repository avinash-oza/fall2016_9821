cmake_minimum_required(VERSION 2.8)
project(cppcodetest)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(CMAKE_BUILD_TYPE, "Debug")

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

include_directories("/usr/include/eigen3")

set(SOURCE_FILES
    Main.cpp
    MatrixDecomposition.hpp
    Solvers.hpp
    Main.cpp IterativeMethods.hpp StoppingCriterionSolvers.hpp)

add_executable(cppcodetest ${SOURCE_FILES})

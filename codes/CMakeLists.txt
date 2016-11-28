cmake_minimum_required(VERSION 2.8)
project(cppcodetest)

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O3 -DEIGEN_NO_DEBUG -fopenmp")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fPIC")
set(CMAKE_BUILD_TYPE, "Debug")

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

include_directories("/usr/include/eigen3")
include_directories("/usr/include")

set(SOURCE_FILES
    Main.cpp
    MatrixDecomposition.hpp
    Solvers.hpp
    IterativeMethods.hpp StoppingCriterionSolvers.hpp Examples.hpp OptionPricers.hpp BlackScholes.hpp
        MonteCarlo.hpp RandomNumberGenerator.hpp BlackScholesOption.hpp PDESolver.hpp PDEMesh.hpp EuropeanPDESolver.hpp)

add_executable(cppcodetest ${SOURCE_FILES})

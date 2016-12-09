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
        MonteCarlo.hpp RandomNumberGenerator.hpp PDESolver.hpp PDEMesh.hpp EuropeanPDESolver.hpp AmericanPDESolver.hpp BarrierPDE.hpp BinomialTrees.hpp BarrierOption.hpp)

add_executable(cppcodetest ${SOURCE_FILES})


################################
# GTest
################################
ADD_SUBDIRECTORY (googletest)
enable_testing()
include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})



################################
# Unit Tests
################################
# Add test cpp file
add_executable( runUnitTests tests/treetests.cpp)
# Link test executable against gtest & gtest_main
target_link_libraries(runUnitTests gtest gtest_main)
add_test( runUnitTests runUnitTests )

cmake_minimum_required(VERSION 2.8)
project(cppcodetest)

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O3 -DEIGEN_NO_DEBUG -fopenmp")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x -fPIC")
set(CMAKE_BUILD_TYPE, "Debug")

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)


include_directories("/usr/include/eigen3")

set(SOURCE_FILES
    FinalMain.cpp
    MatrixDecomposition.hpp
        BarrierPDE.hpp
    Solvers.hpp
    IterativeMethods.hpp StoppingCriterionSolvers.hpp Examples.hpp OptionPricers.hpp BlackScholes.hpp
        MonteCarlo.hpp RandomNumberGenerator.hpp PDESolver.hpp PDEMesh.hpp EuropeanPDESolver.hpp AmericanPDESolver.hpp  BinomialTrees.hpp BarrierOption.hpp
        uFunctions.hpp Exams.hpp FinalMain.cpp)

add_executable(cppcodetest ${SOURCE_FILES})


################################
# GTest
################################
ADD_SUBDIRECTORY (externals/googletest)
enable_testing()
include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})



################################
# Unit Tests
################################
# Add test cpp file
add_executable( runUnitTests tests/finitedifferencetests.cpp)
# Link test executable against gtest & gtest_main
target_link_libraries(runUnitTests gtest gtest_main)
add_test( runUnitTests runUnitTests )

cmake_minimum_required(VERSION 3.6)
project(eigen_test)

include_directories("/usr/include/eigen3")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES main.cpp)
add_executable(eigen_test ${SOURCE_FILES})
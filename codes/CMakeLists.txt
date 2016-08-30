cmake_minimum_required(VERSION 3.6)
project(codes)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES
    CMakeLists.txt
    main.cpp)

include_directories("/usr/include/eigen3")

add_executable(codes ${SOURCE_FILES})
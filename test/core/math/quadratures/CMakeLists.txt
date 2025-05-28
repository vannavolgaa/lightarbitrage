cmake_minimum_required(VERSION 3.29.3)

project(test_quadratures)

set(CMAKE_CXX_STANDARD 14)

# Add subdirectories for math and utils
include_directories(../../../../src/core/math/quadratures)
include_directories(../../../../src/core/math/eigen)

add_executable(
    test 
    test.cpp )



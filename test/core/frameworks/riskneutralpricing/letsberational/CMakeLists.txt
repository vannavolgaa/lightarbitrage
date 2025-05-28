cmake_minimum_required(VERSION 3.29.3)

project(test_lets_be_rational)

set(CMAKE_CXX_STANDARD 17)

add_executable(
    test 
    test.cpp 
    ../../../../../src/core/frameworks/riskneutralpricing/letsberational/lbr.cpp
    ../../../../../src/core/frameworks/riskneutralpricing/normblack/normblack.cpp
    ../../../../../src/core/frameworks/riskneutralpricing/blackscholes/bs.cpp
    ../../../../../src/core/math/probability/probability.cpp
    ../../../../../src/core/math/optimization/optimization.cpp)
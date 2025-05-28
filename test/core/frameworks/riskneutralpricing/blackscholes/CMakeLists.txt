cmake_minimum_required(VERSION 3.29.3)

project(test_black_scholes)

set(CMAKE_CXX_STANDARD 17)


add_executable(
    test 
    test.cpp 
    ../../../../../src/core/frameworks/riskneutralpricing/blackscholes/bs.cpp
    ../../../../../src/core/frameworks/riskneutralpricing/payoff/montecarlo.cpp
    ../../../../../src/core/frameworks/riskneutralpricing/normblack/normblack.cpp
    ../../../../../src/core/math/probability/probability.cpp
    ../../../../../src/core/math/probability/simulation.cpp
    ../../../../../src/core/math/optimization/optimization.cpp)
cmake_minimum_required(VERSION 3.29.3)

project(test_interpolation_2D)

set(CMAKE_CXX_STANDARD 11)

# Add subdirectories for math and utils
include_directories(../../../../src/core/math)

# Create the executable
#add_executable(test test_probability.cpp)

# Link the math and utils libraries to the executable
#target_link_libraries(test PRIVATE probability tools)

add_executable(
    test 
    test.cpp)
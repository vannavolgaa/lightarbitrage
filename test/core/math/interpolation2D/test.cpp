#include <iostream>
#include <cassert>
#include "../../../../src/core/math/interpolation2D.h"

void test_interpolation()
{
    // Sample data for interpolation
    std::map<double, double> test_data = {
        {1.0, 2.0},
        {2.0, 3.0},
        {3.0, 5.0},
        {4.0, 7.0},
        {5.0, 11.0}
    };
    
    // Testing Linear Interpolation
    LinearInterpolation2D linear(test_data);
    
    std::cout << "Testing Linear Interpolation..." << std::endl;
    assert(linear.evaluate(1.0) == 2.0);
    assert(linear.evaluate(3.0) == 5.0);
    assert(linear.evaluate(5.0) == 11.0);
    std::cout << "Linear Interpolation Tests Passed!" << std::endl;
    
    // Testing Cubic Spline Interpolation
    CubicSpline2D cubic(test_data);
    
    std::cout << "Testing Cubic Spline Interpolation..." << std::endl;
    assert(cubic.evaluate(1.0) == 2.0);
    assert(cubic.evaluate(3.0) == 5.0);
    assert(cubic.evaluate(5.0) == 11.0);
    std::cout << "Cubic Spline Interpolation Tests Passed!" << std::endl;
    
    // Testing out-of-range handling
    try {
        linear.evaluate(0.5); // Should throw EvaluationOutOfRange
        assert(false); // Should not reach this line
    } catch (const MathError &e) {
        std::cout << e.what() << std::endl;
    }
    
    try {
        cubic.evaluate(6.0); // Should throw EvaluationOutOfRange
        assert(false);
    } catch (const MathError &e) {
        std::cout << e.what() << std::endl;
    }
}

int main()
{
    test_interpolation();
    return 0;
}

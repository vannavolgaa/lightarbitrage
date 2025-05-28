#include <iostream>
#include <cassert>
#include "../../../../src/core/math/regression.h"

void test_compute_linear_least_square() {
    // Example data
    std::vector<double> Y = {1.0, 2.0, 3.0, 4.0, 5.0};
    std::vector<std::vector<double>> X = {
        {1.0}, {2.0}, {3.0}, {4.0}, {5.0}
    };

    // Perform regression
    RegressionResult result = compute_linear_least_square(Y, X, true);

    // Expected results (for a perfect fit, slope = 1, intercept = 0)
    double expected_intercept = 0.0;
    double expected_coefficient = 1.0;

    // Assertions
    assert(std::abs(result.intercept - expected_intercept) < 1e-6);
    assert(std::abs(result.coefficients[0] - expected_coefficient) < 1e-6);
    assert(result.r_squared > 0.9999); // Perfect fit
    assert(result.rmse < 1e-6);       // Near-zero RMSE
    assert(result.mae < 1e-6);        // Near-zero MAE

    std::cout << "Test passed: compute_linear_least_square" << std::endl;
}

void test_errors() {
    try {
        std::vector<double> Y = {1.0, 2.0, 3.0};
        std::vector<std::vector<double>> X = {{1.0}, {2.0}};
        compute_linear_least_square(Y, X, true);
    } catch (const MathError& e) {
        assert(e.get_code() == MathErrorCode::REGRESSION_MISMATCH_X_Y_SIZE);
        std::cout << "Test passed: Regression mismatch error" << std::endl;
    }
}

int main() {
    test_compute_linear_least_square();
    test_errors();
    return 0;
}
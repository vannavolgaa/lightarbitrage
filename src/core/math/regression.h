#pragma once 
#include <iostream>
#include <vector>
#include "../../../src/core/math/eigen/Eigen/Dense"
#include "eigen_tools.h"
#include "../../../src/errors/core/math.h"

struct RegressionResult 
{
    double intercept; 
    std::vector<double> coefficients;
    std::vector<double> residuals;
    double r_squared;
    double rmse;
    double mae;

    double get_predicted_value(const std::vector<double>& x) const 
    {
        double predicted = intercept;
        for (size_t i = 0; i < coefficients.size(); ++i) 
        {
            predicted += coefficients[i] * x[i];
        }
        return predicted;
    };
}; 

inline void check_matrix_dimensions(const std::vector<std::vector<double>>& X, const std::vector<double>& Y) 
{
    size_t cols = X[0].size();
    for (const auto& col : X) 
    {
        if (col.size() == 0) throw MathError(MathErrorCode::REGRESSION_INVALID_X_MATRIX);
        if (col.size() != cols) throw MathError(MathErrorCode::REGRESSION_INVALID_X_MATRIX);
    }
    if (X.size() != Y.size()) throw MathError(MathErrorCode::REGRESSION_MISMATCH_X_Y_SIZE);
}

inline RegressionResult compute_linear_least_square(std::vector<double>& Y, std::vector<std::vector<double>>& X, bool use_intercept)
{
    check_matrix_dimensions(X, Y);
    if (use_intercept) for (size_t i = 0; i < X.size(); ++i) X[i].insert(X[i].begin(), 1.0); // Add intercept term
    
    MatrixXd A = get_egein_matrix_from_vectors(X);
    VectorXd b = get_eigen_vector_from_vector(Y).transpose();

    MatrixXd AtA = A.transpose() * A;
    VectorXd Atb = A.transpose() * b;
    VectorXd x = AtA.ldlt().solve(Atb);
    
    double intercept = use_intercept ? x(0) : 0.0;
    VectorXd predictions = A * x;
    VectorXd residuals = b - predictions;

    double ss_res = residuals.squaredNorm();
    double ss_tot = (b.array() - b.mean()).matrix().squaredNorm();
    double r_squared = 1.0 - (ss_res / ss_tot);
    double rmse = sqrt(ss_res / b.size());
    double mae = residuals.array().abs().mean();
    
    return {intercept, get_vector_from_eigen_vector(x.segment(1, x.size() - 1)), get_vector_from_eigen_vector(residuals), r_squared, rmse, mae};
}
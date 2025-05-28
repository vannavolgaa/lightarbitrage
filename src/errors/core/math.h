#pragma once 
#include <iostream>
#include "../../../src/errors/base.h"

enum class MathErrorCode
{
    NEWTON_RAPHSON_DERIVATIVE_CLOSE_TO_ZERO = 0,
    INTERPOLATION_OUT_OF_RANGE = 1, 
    INTERPOLATION_WRONG_VECTOR_SIZE = 2,
    NORMAL_DISTRIBUTION_INVALID_SIGMA = 3,
    MISATCH_VECTOR_SIZE_LOSS_FUNCTION = 4, 
    LOSS_ZERO_VALUE_IN_TRUE_VALUES = 5, 
    LOSS_NEGATIVE_VALUE = 6, 
    LOSS_UNKNOWN_TYPE = 7, 
    REGRESSION_INVALID_X_MATRIX = 8,
    REGRESSION_MISMATCH_X_Y_SIZE = 9, 
    UNIFORM_INVALID_BOUNDS = 10, 
    CORRELATION_MATRIX_NOT_SQUARE = 11,
    CORRELATION_MATRIX_NOT_SYMMETRIC = 12,
    CORRELATION_MATRIX_NOT_POSITIVE_DEFINITE = 13,
    CORRELATION_INVALID_COEFFICIENT = 14,
    CORRELATION_INVALID_COEFFICIENT_SIZE = 15

};

class MathError: public CoreArbitrageError<MathErrorCode>
{
    public:
        explicit MathError(const MathErrorCode code): CoreArbitrageError<MathErrorCode>(code){};
        ~MathError() override = default;
    protected: 
        std::string get_module_name() const override {return "math";};
        std::string get_error_message() const override 
        {
            switch (get_code())
            {
                case MathErrorCode::NEWTON_RAPHSON_DERIVATIVE_CLOSE_TO_ZERO:
                    return "The derivative value in the Newton-Raphson optimizer is getting close to zero which makes the opptimizer unstable.";
                case MathErrorCode::INTERPOLATION_OUT_OF_RANGE: 
                    return "It is not possible to interpolate a value out of the interpolation range.";
                case MathErrorCode::INTERPOLATION_WRONG_VECTOR_SIZE:
                    return "The vector to interpolate must have a size strictly superior to 2.";
                case MathErrorCode::NORMAL_DISTRIBUTION_INVALID_SIGMA: 
                    return "Sigma parameter cannot be non-positive."; 
                case MathErrorCode::MISATCH_VECTOR_SIZE_LOSS_FUNCTION:
                    return "The two vectors (true values and estimates) must have the same size.";
                case MathErrorCode::LOSS_ZERO_VALUE_IN_TRUE_VALUES:
                    return "True value is zero in computation (division by zero).";
                case MathErrorCode::LOSS_NEGATIVE_VALUE:
                    return "Negative value in computation.";
                case MathErrorCode::LOSS_UNKNOWN_TYPE:
                    return "Unknown loss function type.";
                case MathErrorCode::REGRESSION_INVALID_X_MATRIX:
                    return "The X matrix is invalid. The number of cols must be consistent for each row. It also cannot be empty.";
                case MathErrorCode::REGRESSION_MISMATCH_X_Y_SIZE:
                    return "The number of rows in X must match the number of elements in Y.";
                case MathErrorCode::UNIFORM_INVALID_BOUNDS:
                    return "The bounds for the uniform distribution must be valid (lower bound must be less than upper bound).";
                case MathErrorCode::CORRELATION_MATRIX_NOT_SQUARE:
                    return "The correlation matrix must be square (same number of rows and columns).";
                case MathErrorCode::CORRELATION_MATRIX_NOT_SYMMETRIC:
                    return "The correlation matrix must be symmetric (equal to its transpose).";
                case MathErrorCode::CORRELATION_MATRIX_NOT_POSITIVE_DEFINITE:
                    return "The correlation matrix must be positive definite (all eigenvalues must be positive).";
                case MathErrorCode::CORRELATION_INVALID_COEFFICIENT:
                    return "The correlation coefficient must be in the range [-1, 1].";
                case MathErrorCode::CORRELATION_INVALID_COEFFICIENT_SIZE:
                    return "The number of correlation coefficients must match the size of the correlation matrix.";
            }
            return "Unknown math error.";
        };
};
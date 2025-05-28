#pragma once
#include <iostream>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <numeric>
#include <algorithm>
#include "../../../src/errors/core/math.h" 

// Define the enum class
enum class LossType {
    MSE,
    RMSE,
    MAE,
    MAPE,
    MSLE, 
    MSPE, 
    RMPSE
};

// Helper check function
inline void check_vector_sizes(const std::vector<double>& estimate, const std::vector<double>& true_values) {
    if (estimate.size() != true_values.size()) throw MathError(MathErrorCode::MISATCH_VECTOR_SIZE_LOSS_FUNCTION);
}

// Individual loss functions
inline double mean_squared_error(const std::vector<double>& estimate, const std::vector<double>& true_values) {
    check_vector_sizes(estimate, true_values);
    double sum = 0.0;
    for (size_t i = 0; i < estimate.size(); ++i) {
        double diff = estimate[i] - true_values[i];
        sum += diff * diff;
    }
    return sum / estimate.size();
}

inline double mean_squared_percentage_error(const std::vector<double>& estimate, const std::vector<double>& true_values) {
    check_vector_sizes(estimate, true_values);
    double sum = 0.0;
    for (size_t i = 0; i < estimate.size(); ++i) {
        if (abs(true_values[i]) <= 1e-12) throw MathError(MathErrorCode::LOSS_ZERO_VALUE_IN_TRUE_VALUES);
        double diff = (estimate[i] - true_values[i]) / true_values[i];
        sum += diff * diff;
    }
    return sum / estimate.size();
}

inline double root_mean_squared_percentage_error(const std::vector<double>& estimate, const std::vector<double>& true_values) {
    return sqrt(mean_squared_percentage_error(estimate, true_values));
}


inline double root_mean_squared_error(const std::vector<double>& estimate, const std::vector<double>& true_values) {
    return std::sqrt(mean_squared_error(estimate, true_values));
}

inline double mean_absolute_error(const std::vector<double>& estimate, const std::vector<double>& true_values) {
    check_vector_sizes(estimate, true_values);
    double sum = 0.0;
    for (size_t i = 0; i < estimate.size(); ++i) {
        sum += std::abs(estimate[i] - true_values[i]);
    }
    return sum / estimate.size();
}

inline double mean_absolute_percentage_error(const std::vector<double>& estimate, const std::vector<double>& true_values) {
    check_vector_sizes(estimate, true_values);
    double sum = 0.0;
    for (size_t i = 0; i < estimate.size(); ++i) {
        if (abs(true_values[i]) <= 1e-12) throw MathError(MathErrorCode::LOSS_ZERO_VALUE_IN_TRUE_VALUES);
        sum += std::abs((true_values[i] - estimate[i]) / true_values[i]);
    }
    return 100.0 * sum / estimate.size();
}

inline double mean_squared_log_error(const std::vector<double>& estimate, const std::vector<double>& true_values) {
    check_vector_sizes(estimate, true_values);
    double sum = 0.0;
    for (size_t i = 0; i < estimate.size(); ++i) {
        if (estimate[i] < -1.0 || true_values[i] < -1.0) throw MathError(MathErrorCode::LOSS_NEGATIVE_VALUE);
        double log_est = std::log1p(std::max(estimate[i], 0.0));  // log(1 + max(0, x))
        double log_true = std::log1p(std::max(true_values[i], 0.0));
        double diff = log_est - log_true;
        sum += diff * diff;
    }
    return sum / estimate.size();
}

// General function
inline double compute_loss(const std::vector<double>& estimate, const std::vector<double>& true_values, LossType loss_type) {
    switch (loss_type) {
        case LossType::MSE:
            return mean_squared_error(estimate, true_values);
        case LossType::RMSE:
            return root_mean_squared_error(estimate, true_values);
        case LossType::MAE:
            return mean_absolute_error(estimate, true_values);
        case LossType::MAPE:
            return mean_absolute_percentage_error(estimate, true_values);
        case LossType::MSLE:
            return mean_squared_log_error(estimate, true_values);
        case LossType::MSPE:
            return mean_squared_percentage_error(estimate, true_values);
        case LossType::RMPSE:
            return root_mean_squared_percentage_error(estimate, true_values);
        default:
            throw MathError(MathErrorCode::LOSS_UNKNOWN_TYPE);
    }
}


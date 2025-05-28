#include "stats.h"

void CorrelationMatrix::check_correlation_matrix(const Eigen::MatrixXd& matrix)
{
    bool is_symmetric = (matrix - matrix.transpose()).norm() < 1e-10;   
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(matrix);
    bool is_psd = (solver.eigenvalues().array() > 0).all();

    if (matrix.rows() != matrix.cols()) throw MathError(MathErrorCode::CORRELATION_MATRIX_NOT_SQUARE);
    if (!is_symmetric) throw MathError(MathErrorCode::CORRELATION_MATRIX_NOT_SYMMETRIC);
    if (!is_psd) throw MathError(MathErrorCode::CORRELATION_MATRIX_NOT_POSITIVE_DEFINITE);
}

Eigen::MatrixXd CorrelationMatrix::construct_2x2_matrix(const double coeff)
{
    if (abs(coeff) > 1.0) throw MathError(MathErrorCode::CORRELATION_INVALID_COEFFICIENT);
    Eigen::MatrixXd matrix(2, 2);
    matrix << 1.0, coeff, coeff, 1.0;
    return matrix;
}

Eigen::MatrixXd CorrelationMatrix::construct_nxn_matrix(const std::vector<double> coeffs)
{
    size_t c = coeffs.size();
    double p_double = (1.0 + std::sqrt(1.0 + 8.0 * c)) / 2.0;
    size_t p = static_cast<size_t>(std::round(p_double));

    if (p * (p - 1) / 2 != c) throw MathError(MathErrorCode::CORRELATION_INVALID_COEFFICIENT_SIZE);

    Eigen::MatrixXd corr = Eigen::MatrixXd::Identity(p, p);

    size_t idx = 0;
    for (size_t i = 0; i < p; ++i) {
        for (size_t j = i + 1; j < p; ++j) {
            if (abs(coeffs[idx]) > 1.0) throw MathError(MathErrorCode::CORRELATION_INVALID_COEFFICIENT);
            corr(i, j) = coeffs[idx];
            corr(j, i) = coeffs[idx]; // symmetric
            ++idx;
        }
    }

    check_correlation_matrix(corr);

    return corr;
}
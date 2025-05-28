#pragma once 
#include <iostream>
#include "../../../../src/core/math/eigen_tools.h"
#include "../../../../src/errors/core/math.h"

class CorrelationMatrix
{
    public: 
        CorrelationMatrix() : corr_matrix_(Eigen::MatrixXd::Identity(1, 1)) {};
        CorrelationMatrix(const double coeff) : corr_matrix_(construct_2x2_matrix(coeff)) {};
        CorrelationMatrix(const std::vector<double>& coeffs) : corr_matrix_(construct_nxn_matrix(coeffs)) {};
        CorrelationMatrix(const Eigen::MatrixXd& matrix) : corr_matrix_(matrix) {check_correlation_matrix(matrix);};
        ~CorrelationMatrix() = default;

        Eigen::MatrixXd get_matrix() const { return corr_matrix_; }

        Eigen::MatrixXd get_cholesky_decomposition() const 
        {
            Eigen::LLT<Eigen::MatrixXd> llt(corr_matrix_);
            return llt.matrixL();
        }

        int get_dimension() const {return corr_matrix_.rows();}

    private: 
        Eigen::MatrixXd corr_matrix_;
        void check_correlation_matrix(const Eigen::MatrixXd& matrix);
        Eigen::MatrixXd construct_2x2_matrix(const double coeff); 
        Eigen::MatrixXd construct_nxn_matrix(const std::vector<double> coeffs); 
};
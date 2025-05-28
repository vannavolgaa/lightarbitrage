#pragma once 
#include <iostream>
#include <vector>
#include "../../../src/core/math/eigen/Eigen/Dense"

using namespace std;
using namespace Eigen;

inline MatrixXd get_egein_matrix_from_vectors(const std::vector<std::vector<double>>& data) {
    int rows = data.size();
    int cols = data[0].size();
    MatrixXd mat(rows, cols);

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            mat(i, j) = data[i][j];
        }
    }
    return mat;
}

inline VectorXd get_eigen_vector_from_vector(const std::vector<double>& data) {
    return Map<const VectorXd>(data.data(), data.size());
}

inline std::vector<double> get_vector_from_eigen_vector(const VectorXd& data) {
    return std::vector<double>(data.data(), data.data() + data.size());
}


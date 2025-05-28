#pragma once
#include <iostream>
#include <vector>
#include <functional>
#include "../../../src/core/math/eigen/Eigen/Eigenvalues"
#include "eigen_tools.h"

struct GaussianQuadrature
{
    const int points;
    const std::vector<double> roots;
    const std::vector<double> weights;
    GaussianQuadrature(int points_, const std::vector<double>& roots_, const std::vector<double>& weights_)
        : points(points_), roots(roots_), weights(weights_) {}
    ~GaussianQuadrature() = default;
    double integrate(std::function<double(double)> f) const
    {
        double integral = 0.0;
        for (int i = 0; i < points; i++){integral += weights[i] * f(roots[i]);}
        return integral;
    }
};

inline std::shared_ptr<GaussianQuadrature> get_gauss_laguerre_quadrature(int points)
{
    points = std::max(2, points);
    // Construct the symmetric tridiagonal matrix (Golub-Welsch algorithm)
    Eigen::MatrixXd J = Eigen::MatrixXd::Zero(points, points);
    
    // Diagonal elements: a_k = 2k + 1
    for (int k = 0; k < points; k++)
    {
        J(k, k) = 2.0 * k + 1.0;
    }
    
    // Off-diagonal elements: b_k = k
    for (int k = 0; k < points - 1; k++)
    {
        double b = k + 1.0;
        J(k, k + 1) = b;
        J(k + 1, k) = b;
    }

    // Compute eigenvalues (roots) and eigenvectors
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(J);
    Eigen::VectorXd roots = solver.eigenvalues();
    
    // Compute weights using the first component of normalized eigenvectors
    std::vector<double> weights(points);
    Eigen::MatrixXd V = solver.eigenvectors();
    for (int i = 0; i < points; i++)
    {
        weights[i] = std::pow(V(0, i), 2);
    }

    return std::make_shared<GaussianQuadrature>(points, get_vector_from_eigen_vector(roots), weights);
}


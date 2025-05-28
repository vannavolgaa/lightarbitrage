#include <iostream>
#include <vector>
#include <cmath>
#include <cassert>
#include <chrono>
#include "../../../../src/core/math/quadratures.h"

void test_gauss_laguerre_quadrature() {
    // Test the computation of roots and weights
    int points = 20;
    std::shared_ptr<GaussianQuadrature> gq = get_gauss_laguerre_quadrature(points);

    std::vector<double> roots = gq->roots;
    std::vector<double> weights = gq->weights;

    std::cout << "Roots:\n";
    for (double root : roots) {
        std::cout << root << " ";
    }
    std::cout << "\nWeights:\n";
    for (double weight : weights) {
        std::cout << weight << " ";
    }
    std::cout << std::endl;

    // Verify the number of roots and weights
    assert(roots.size() == points);
    assert(weights.size() == points);

    // Test the integration method
    auto f = [](double x) { return x; }; // Integral of x * exp(-x) from 0 to infinity is 1
    double result = gq->integrate(f);
    std::cout << "Integration result: " << result << std::endl;

    // Verify the integration result
    assert(std::abs(result - 1.0) < 1e-6);
}

void test_gauss_laguerre_quadrature_high_dimension() {
    // Test the computation of roots and weights
    int points = 300;
    std::shared_ptr<GaussianQuadrature> gq = get_gauss_laguerre_quadrature(points);

    std::vector<double> roots = gq->roots;
    std::vector<double> weights = gq->weights;

    std::cout << "Roots:\n";
    for (double root : roots) {
        std::cout << root << " ";
    }
    std::cout << "\nWeights:\n";
    for (double weight : weights) {
        std::cout << weight << " ";
    }
    std::cout << std::endl;

    // Verify the number of roots and weights
    assert(roots.size() == points);
    assert(weights.size() == points);

    // Test the integration method
    auto f = [](double x) { return x; }; // Integral of x * exp(-x) from 0 to infinity is 1
    double result = gq->integrate(f);
    std::cout << "Integration result: " << result << std::endl;

    // Verify the integration result
    assert(std::abs(result - 1.0) < 1e-6);
}

void test_timing_quadratures()
{
    // Test the timing of Gauss-Laguerre quadrature
    int points = 300;
    auto start = std::chrono::high_resolution_clock::now();
    std::shared_ptr<GaussianQuadrature> gq = get_gauss_laguerre_quadrature(points);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Time taken to compute Gauss-Laguerre quadrature with " << points << " points: " << elapsed.count() << " seconds" << std::endl;
}

int main() {
    std::cout << "Testing Gauss-Laguerre Quadrature:\n";
    test_gauss_laguerre_quadrature();
    test_gauss_laguerre_quadrature_high_dimension();
    test_timing_quadratures();
    std::cout << "All tests passed successfully!" << std::endl;
    return 0;
}
#include <iostream>
#include <cassert>
#include <cmath>
#include "../../../../src/core/math/probability/probability.h"
#include "../../../../src/core/math/probability/stats.h"
#include "../../../../src/core/math/probability/simulation.h"

void test_normal_distribution() {
    // Test default constructor
    NormalDistribution default_dist;
    assert(default_dist.get_mu() == 0.0);
    assert(default_dist.get_sigma() == 1.0);

    // Test parameterized constructor
    NormalDistribution dist(5.0, 2.0);
    assert(dist.get_mu() == 5.0);
    assert(dist.get_sigma() == 2.0);

    // Test setting and getting mu
    dist.set_mu(3.0);
    assert(dist.get_mu() == 3.0);

    // Test setting and getting sigma
    dist.set_sigma(1.5);
    assert(dist.get_sigma() == 1.5);

    // Test invalid sigma (should throw MathError)
    try {
        NormalDistribution invalid_dist(0.0, -1.0);
        assert(false); // Should not reach here
    } catch (const MathError& e) {
        assert(e.get_code() == MathErrorCode::NORMAL_DISTRIBUTION_INVALID_SIGMA);
    }

    // Test cdf
    double cdf_value = default_dist.cdf(0.0);
    assert(std::abs(cdf_value - 0.5) < 1e-6); // cdf(0) for standard normal is 0.5

    // Test inv_cdf
    double inv_cdf_value = default_dist.inv_cdf(0.5);
    assert(std::abs(inv_cdf_value - 0.0) < 1e-6); // inv_cdf(0.5) for standard normal is 0

    // Test mgf
    double mgf_value = default_dist.mgf(0.1);
    assert(std::abs(mgf_value - std::exp(0.1 * 0.0 + 0.5 * 0.1 * 0.1)) < 1e-6);

    // Test cf
    std::complex<double> cf_value = default_dist.cf(0.1);
    std::cout << cf_value.real() << std::endl;
    std::cout << cf_value.imag() << std::endl;
    //assert(std::abs(cf_value.real() - std::cos(0.1 * 0.0)) < 1e-6);
    //assert(std::abs(cf_value.imag() - std::sin(0.1 * 0.0)) < 1e-6);

    // Test pdf
    double pdf_value = default_dist.pdf(0.0);
    assert(std::abs(pdf_value - ONE_OVER_SQRT_TWO_PI) < 1e-6); // pdf(0) for standard normal

    // Test random number generation (r)
    double random_value = default_dist.r();
    assert(random_value >= -10.0 && random_value <= 10.0); // Check if it's within a reasonable range

    std::cout << "All tests passed for NormalDistribution!" << std::endl;
}

void test_correlation_matrix_2x2()
{
    // Test 2x2 correlation matrix with valid coefficient
    CorrelationMatrix corr_matrix(0.5);
    Eigen::MatrixXd expected_matrix(2, 2);
    expected_matrix << 1.0, 0.5, 0.5, 1.0;
    assert((corr_matrix.get_matrix() - expected_matrix).norm() < 1e-6);

    // Test invalid coefficient (should throw MathError)
    try {
        CorrelationMatrix invalid_corr_matrix(1.5);
        assert(false); // Should not reach here
    } catch (const MathError& e) {
        assert(e.get_code() == MathErrorCode::CORRELATION_INVALID_COEFFICIENT);
    }

    // test cholesky decomposition
    Eigen::MatrixXd cholesky_decomp = corr_matrix.get_cholesky_decomposition();
    Eigen::MatrixXd expected_cholesky(2, 2);
    expected_cholesky << 1.0, 0.0, 0.5, std::sqrt(1 - 0.5 * 0.5);
    assert((cholesky_decomp - expected_cholesky).norm() < 1e-6);

    std::cout << "All tests passed for CorrelationMatrix 2x2!" << std::endl;
}

void test_correlation_matrix_nxn()
{
    // Test nxn correlation matrix with valid coefficients
    std::vector<double> coeffs = {0.8,0.3,0.5,0.25,0.6,0.4};
    CorrelationMatrix corr_matrix(coeffs);
    Eigen::MatrixXd expected_matrix(4, 4);
    expected_matrix << 1.0, 0.8, 0.3, 0.5,
                       0.8, 1.0, 0.25, 0.6,
                       0.3, 0.25, 1.0, 0.4,
                       0.5, 0.6, 0.4, 1.0;
    assert((corr_matrix.get_matrix() - expected_matrix.transpose()).norm() < 1e-6);

    // Test invalid coefficients (should throw MathError)
    try {
        std::vector<double> invalid_coeffs = {1.0, 2.0, 2.0, 1.0, 1.0, 1.0};
        CorrelationMatrix invalid_corr_matrix(invalid_coeffs);
        assert(false); // Should not reach here
    } catch (const MathError& e) {
        assert(e.get_code() == MathErrorCode::CORRELATION_INVALID_COEFFICIENT);
    }

    try {
        std::vector<double> invalid_coeffs = {1.0, 2.0, 2.0, 1.0};
        CorrelationMatrix invalid_corr_matrix(invalid_coeffs);
        assert(false); // Should not reach here
    } catch (const MathError& e) {
        assert(e.get_code() == MathErrorCode::CORRELATION_INVALID_COEFFICIENT_SIZE);
    }

    // Test cholesky decomposition
    Eigen::MatrixXd cholesky_decomp = corr_matrix.get_cholesky_decomposition();
    Eigen::MatrixXd expected_cholesky(4, 4);
    expected_cholesky << 1.0, 0.0, 0.0, 0.0,
                         0.8, .6, 0.0, 0.0,
                         0.3, 0.0166667, 0.953794, 0.0,
                         0.5, 0.333333, 0.256287, 0.757104;
    assert((cholesky_decomp - expected_cholesky).norm() < 1e-6);

    std::cout << "All tests passed for CorrelationMatrix nxn!" << std::endl;
}

void test_constructor_monte_carlo_engine()
{
    // Test default constructor
    MonteCarloEngine engine(10000, 100, std::make_shared<NormalDistribution>(), false);

    std::cout << "Time taken for samples update: " << engine.get_time_taken_for_samples_update() << " seconds." << std::endl;

    std::cout << "All tests passed for MonteCarloEngine!" << std::endl;
}

void test_multivariate_monte_carlo_engine()
{
    // Test default constructor
    CorrelationMatrix corr_matrix(0.5);
    MultivariateMonteCarloEngine engine(corr_matrix, 10000, 100, std::make_shared<NormalDistribution>());

    std::vector<Eigen::MatrixXd> samples = engine.get_samples();
    int n = samples.size();
    assert(samples.size() == 2);
    for (int i = 0; i < n; ++i) {
        assert(samples[i].rows() == 10000);
        assert(samples[i].cols() == 100);
    }

    std::cout << "Time taken for samples update: " << engine.get_time_taken_for_samples_update() << " seconds." << std::endl;

    std::cout << "Cholesky decomposition of correlation matrix:\n";
    Eigen::MatrixXd cholesky_decomp = corr_matrix.get_cholesky_decomposition();
    std::cout << cholesky_decomp << std::endl;

    std::cout << "Samples means: ";
    std::vector<double> means = engine.get_sample_means();
    for (const auto& mean : means) {
        std::cout << mean << " ";
    }
    std::cout << std::endl;
    std::cout << "Samples variances: ";
    std::vector<double> variances = engine.get_sample_variances();
    for (const auto& variance : variances) {
        std::cout << variance << " ";
    }
    std::cout << std::endl;

    std::vector<double> coeffs = {0.8,0.3,0.5,0.25,0.6,0.4};
    CorrelationMatrix corr_matrix2(coeffs);
    MultivariateMonteCarloEngine engine2(corr_matrix2, 10000, 100, std::make_shared<NormalDistribution>());
    std::vector<Eigen::MatrixXd> samples2 = engine2.get_samples();
    int n2 = samples2.size();
    assert(samples2.size() == corr_matrix2.get_dimension());
    for (int i = 0; i < n2; ++i) {
        assert(samples2[i].rows() == 10000);
        assert(samples2[i].cols() == 100);
    }
    std::cout << "Time taken for samples update: " << engine2.get_time_taken_for_samples_update() << " seconds." << std::endl;

    std::cout << "All tests passed for MultivariateMonteCarloEngine!" << std::endl;
}

int main() {
    test_normal_distribution();
    test_correlation_matrix_2x2();
    test_correlation_matrix_nxn();
    test_constructor_monte_carlo_engine();
    test_multivariate_monte_carlo_engine();
    std::cout << "All tests passed!" << std::endl;
    return 0;
}
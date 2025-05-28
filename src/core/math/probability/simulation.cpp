#include "simulation.h"

void MonteCarloEngine::update_samples()
{
    auto start = std::chrono::high_resolution_clock::now(); 
    UniformDistribution uniform = UniformDistribution(0.0, 1.0);

    if (antithetic_)
    {
        int N2 = N_/2;
        Eigen::MatrixXd u = Eigen::MatrixXd::NullaryExpr(N2, M_,[&]() { return uniform.r(); });
        Eigen::MatrixXd samples1 = u.unaryExpr([this](double x) { return distribution_->inv_cdf(x); });
        Eigen::MatrixXd samples2 = u.unaryExpr([this](double x) { return distribution_->inv_cdf(1-x); });
        samples_.resize(N_, M_);
        samples_.topRows(N2) = samples1;
        samples_.bottomRows(N2) = samples2;

    }
    else samples_ = Eigen::MatrixXd::NullaryExpr(N_, M_, [&]() { return distribution_->inv_cdf(uniform.r()); });
    auto end = std::chrono::high_resolution_clock::now();
    time_taken_for_update_ = std::chrono::duration<double>(end - start).count();
    updated_ = true;
}

void MultivariateMonteCarloEngine::update_samples()
{
    update_samples_non_antithetic(); 
}

void MultivariateMonteCarloEngine::update_samples_non_antithetic()
{
    auto start = std::chrono::high_resolution_clock::now();
    Eigen::MatrixXd L = correlation_matrix_.get_cholesky_decomposition();
    MonteCarloEngine engine(N_*M_, L.rows(),distribution_, false);
    Eigen::MatrixXd U = engine.get_samples();
    Eigen::MatrixXd Z = U*L.transpose();
    for (int i = 0; i < Z.cols(); ++i)
    {
        Eigen::MatrixXd S = Z.col(i);
        S.resize(N_,M_);
        samples_.push_back(S);
    }
    auto end = std::chrono::high_resolution_clock::now();
    time_taken_for_update_ = std::chrono::duration<double>(end - start).count();
    updated_ = true;
}

std::vector<double> MultivariateMonteCarloEngine::get_sample_means()
{
    if (!updated_) update_samples();
    std::vector<double> means;
    for (const auto& sample : samples_)
    {
        means.push_back(sample.mean());
    }
    return means;
}

std::vector<double> MultivariateMonteCarloEngine::get_sample_variances()
{
    if (!updated_) update_samples();
    std::vector<double> variances;
    for (const auto& sample : samples_)
    {
        variances.push_back(sample.array().square().mean() - sample.mean() * sample.mean());
    }
    return variances;
}
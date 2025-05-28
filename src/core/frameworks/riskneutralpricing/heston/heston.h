#pragma once 
#include <iostream>
#include "../../../../../src/core/frameworks/riskneutralpricing/normblack/normblack.h"
#include "../../../../../src/core/math/quadratures.h"
#include "../../../../../src/core/frameworks/riskneutralpricing/payoff/montecarlo.h"
#include "../../../../../src/core/math/probability/simulation.h"
#include "../../../../../src/errors/core/frameworks.h"

class Heston 
{
    public: 
        Heston(const double kappa, const double theta, const double eta, const double rho, const double v0)
            : kappa_(kappa), theta_(theta), eta_(eta), rho_(rho), v0_(v0), quadrature_(get_gauss_laguerre_quadrature(20))
        {
            if (kappa_ < 0 or eta_<0 or v0_<0 or theta_<0 or abs(rho_)>=1) throw FrameworksError(FrameworksErrorCode::HESTON_INVALID_PARAM);
        };
        Heston(const double kappa, const double theta, const double eta, const double rho, const double v0, int quad_points)
            : kappa_(kappa), theta_(theta), eta_(eta), rho_(rho), v0_(v0), quadrature_(get_gauss_laguerre_quadrature(quad_points))
        {
            if (kappa_ < 0 or eta_<0 or v0_<0 or theta_<0 or abs(rho_)>=1) throw FrameworksError(FrameworksErrorCode::HESTON_INVALID_PARAM);
        };
        ~Heston() = default;
        
        double get_kappa() const { return kappa_; }
        double get_theta() const { return theta_; }
        double get_eta() const { return eta_; }
        double get_rho() const { return rho_; }
        double get_v0() const { return v0_; }
        std::shared_ptr<GaussianQuadrature> get_quadrature() const { return quadrature_; }
        std::complex<double> get_characteristic_function(std::complex<double> u, double T);
        bool is_feller_condition_satisfied() const { return 2.0 * kappa_ * theta_ > eta_ * eta_;};
        double get_lewis_normalized_price(const double x, const double T, const bool is_call);
        double get_lewis_price(const double F, const double K, const double T, const double r, const bool is_call);

        void set_kappa(double kappa) { kappa_ = kappa; }
        void set_theta(double theta) { theta_ = theta; }
        void set_eta(double eta) { eta_ = eta;  }
        void set_rho(double rho) { rho_ = rho; }
        void set_v0(double v0) { v0_ = v0; }
        void set_quadrature(int points) { quadrature_ = get_gauss_laguerre_quadrature(points); }

    private: 
        double kappa_;  // Speed of mean reversion
        double theta_;  // Long-term mean variance
        double eta_;    // Volatility of volatility
        double rho_;    // Correlation between asset and volatility processes
        double v0_;     // Initial variance
        std::shared_ptr<GaussianQuadrature> quadrature_;
        
}; 

enum class HestonDiscretizationMethod
{
    EULER = 0,
    MILSTEIN = 1
};

class HestonSimulation
{
    public: 
        HestonSimulation(
            const double S, const double mu, const double T, const double kappa, const double theta, const double eta, const double rho, const double v0,
            const int number_iterations, const int sample_size)
            : engine_(CorrelationMatrix(rho),number_iterations, sample_size, std::make_shared<NormalDistribution>(0.0, 1.0)), 
            model_(kappa, theta, eta, rho, v0, 2), S_(S), mu_(mu), T_(T), full_truncation_(true), discretization_method_(HestonDiscretizationMethod::MILSTEIN),
            is_updated_(false), time_taken_for_simulation_(0.0), simulated_prices(engine_.get_samples()[1]), simulated_variances(engine_.get_samples()[0])
        {
            update();
        };
        HestonSimulation(
            const double S, const double mu, const double T, const double kappa, const double theta, const double eta, const double rho, const double v0,
            const HestonDiscretizationMethod discretization_method, const int number_iterations, const int sample_size)
            : engine_(CorrelationMatrix(rho),number_iterations, sample_size, std::make_shared<NormalDistribution>(0.0, 1.0)), 
            model_(kappa, theta, eta, rho, v0, 2), S_(S), mu_(mu), T_(T), full_truncation_(true), discretization_method_(discretization_method),
            is_updated_(false), time_taken_for_simulation_(0.0), simulated_prices(engine_.get_samples()[1]), simulated_variances(engine_.get_samples()[0])
        {
            update();
        };
        ~HestonSimulation() = default;

        double get_S() const { return S_; }
        double get_mu() const { return mu_; }
        double get_T() const { return T_; }
        double get_kappa() const { return model_.get_kappa(); };
        double get_theta() const { return model_.get_theta(); };
        double get_eta() const { return model_.get_eta(); };
        double get_rho() const { return model_.get_rho(); };
        double get_v0() const { return model_.get_v0(); };
        bool is_updated() const { return is_updated_; };
        HestonDiscretizationMethod get_discretization_method() const { return discretization_method_; };
        Eigen::MatrixXd get_simulated_prices() { if (!is_updated_) update(); return simulated_prices; };
        MultivariateMonteCarloEngine get_simulation_engine() const { return engine_; };
        MonteCarloPayoffEngine get_payoff_engine() { if (!is_updated_) update(); return MonteCarloPayoffEngine(simulated_prices, T_); };
        double get_time_taken_for_simulation() const { return time_taken_for_simulation_; };
        double get_variance_from_truncation(double v) const { return full_truncation_ ? std::max(v,0.0) : std::abs(v); }
        Heston get_heston() const { return model_; }

        void set_S(double S) { S_ = S; is_updated_ = false; };
        void set_mu(double mu) { mu_ = mu; is_updated_ = false; };
        void set_T(double T) { T_ = T; is_updated_ = false; };
        void set_kappa(double kappa) { model_.set_kappa(kappa); is_updated_ = false; };
        void set_theta(double theta) { model_.set_theta(theta); is_updated_ = false; };
        void set_eta(double eta) { model_.set_eta(eta); is_updated_ = false; };
        void set_rho(double rho) { model_.set_rho(rho); is_updated_ = false; };
        void set_v0(double v0) { model_.set_v0(v0); is_updated_ = false; };
        void set_sample_size(int sample_size) { engine_.set_sample_dimension(sample_size, engine_.get_number_iterations()); is_updated_ = false; }
        void set_number_iterations(int number_iterations) { engine_.set_sample_dimension(engine_.get_sample_size(), number_iterations); is_updated_ = false; }
        void set_full_truncation(bool full_truncation) { full_truncation_ = full_truncation; is_updated_ = false; }
        void set_discretization_method(HestonDiscretizationMethod method) { discretization_method_ = method; is_updated_ = false; }
        void update();
    
    private: 
        MultivariateMonteCarloEngine engine_; 
        Heston model_; 
        double S_; // Underlying asset price
        double mu_; // Risk-free interest rate
        double T_; // Time to maturity
        bool full_truncation_;
        HestonDiscretizationMethod discretization_method_;  // Discretization method 
        bool is_updated_;
        double time_taken_for_simulation_;
        Eigen::MatrixXd simulated_prices; 
        Eigen::MatrixXd simulated_variances;
        void update_euler_milstein(bool milstein); 
        Heston get_model(double kappa, double theta, double eta, double rho, double v0) const {return Heston(kappa, theta, eta, rho, v0, 2);}
}; 

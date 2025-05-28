#pragma once 
#include <iostream>
#include "../../../../../src/core/frameworks/riskneutralpricing/normblack/normblack.h"
#include "../../../../../src/core/frameworks/riskneutralpricing/payoff/montecarlo.h"
#include "../../../../../src/core/math/optimization/optimization.h"
#include "../../../../../src/core/math/probability/simulation.h"
#include "../../../../../src/errors/core/frameworks.h"

class BlackScholes 
{
    public: 
        BlackScholes(const double S, const double K, const double T, const double r, const double q, const double sigma, const bool is_call, const bool is_future): 
        S_(S), K_(K), T_(T), r_(r), q_(q), sigma_(sigma), call_put_flag_(get_call_put_flag(is_call)), future_flag_(get_future_flag(is_future)), is_updated_(true)
        {
            update();
        };
        ~BlackScholes() = default;

        double get_S() const {return S_;};
        double get_K() const {return K_;};
        double get_T() const {return T_;};
        double get_r() const {return r_;};
        double get_q() const {return q_;};
        double get_sigma() const {return sigma_;};
        bool get_is_call() const {return call_put_flag_ == 1 ? true : false;};
        bool get_is_future() const {return future_flag_ == 0 ? true : false;};
        
        void set_S(double S) {S_ = S; is_updated_ = false;};
        void set_K(double K) {K_ = K; is_updated_ = false;};
        void set_T(double T) {T_ = T; is_updated_ = false;};
        void set_r(double r) {r_ = r; is_updated_ = false;};
        void set_q(double q) {q_ = q; is_updated_ = false;};
        void set_sigma(double sigma) {sigma_ = sigma; is_updated_ = false;};
        void set_is_call(bool is_call) {call_put_flag_ = get_call_put_flag(is_call); is_updated_ = false;};
        void set_is_future(bool is_future) {future_flag_ = get_future_flag(is_future); is_updated_ = false;};

        double get_future_price(); 

        double get_price();
        double get_delta();
        double get_gamma();
        double get_vega();
        double get_theta();
        double get_dual_delta();
        double get_dual_gamma();
        double get_lewis_price(std::shared_ptr<GaussianQuadrature> quadrature);

    private: 
        double S_; // Underlying asset price
        double K_; // Strike price
        double T_; // Time to maturity
        double r_; // Risk-free interest rate
        double q_; // Dividend yield
        double sigma_; // Volatility
        int call_put_flag_; // Call or put option
        int future_flag_; // Future or option contract
        bool is_updated_; 
        double F_; 
        double d1_; 
        double d2_; 
        double Nd1_; 
        double Nd2_;
        double nd1_; 
        double nd2_;
        void update(); 
        void check_update(){if (!is_updated_) update();};
};

class BaroneAdesiWhaley 
{
    public: 
        BaroneAdesiWhaley(const double S, const double K, const double T, const double r, const double q, const double sigma, const bool is_call, const bool is_future): 
        S_(S), K_(K), T_(T), r_(r), q_(q), sigma_(sigma), 
        call_put_flag_(get_call_put_flag(is_call)), future_flag_(get_future_flag(is_future)), is_updated_(true), max_iter_(25), epsilon_(1e-6)
        {
            update();
        };
        ~BaroneAdesiWhaley() = default;

        double get_S() const {return S_;};
        double get_K() const {return K_;};
        double get_T() const {return T_;};
        double get_r() const {return r_;};
        double get_q() const {return q_;};
        double get_sigma() const {return sigma_;};
        bool get_is_call() const {return call_put_flag_ == 1 ? true : false;};
        bool get_is_future() const {return future_flag_ == 0 ? true : false;};

        void set_S(double S) {S_ = S; is_updated_ = false;};
        void set_K(double K) {K_ = K; is_updated_ = false;};
        void set_T(double T) {T_ = T; is_updated_ = false;};
        void set_r(double r) {r_ = r; is_updated_ = false;};
        void set_q(double q) {q_ = q; is_updated_ = false;};
        void set_sigma(double sigma) {sigma_ = sigma; is_updated_ = false;};
        void set_is_call(bool is_call) {call_put_flag_ = get_call_put_flag(is_call); is_updated_ = false;};
        void set_is_future(bool is_future) {future_flag_ = get_future_flag(is_future); is_updated_ = false;};

        double get_exercise_premium();
        double get_price();
    private: 
        double S_; // Underlying asset price
        double K_; // Strike price
        double T_; // Time to maturity
        double r_; // Risk-free interest rate
        double q_; // Dividend yield
        double sigma_; // Volatility
        int call_put_flag_; // Call or put option
        int future_flag_; // Future or option contract
        bool is_updated_; 
        int max_iter_;
        double epsilon_; // Perturbation parameter
        double qq_inf_; 
        double qq_; 
        double Se0_;
        double Se_; 
        void update(); 
        void check_update(){if (!is_updated_) update();};
        BlackScholes get_bs(double s) const;
        double target_function(double s); 
        double derivative_target_function(double s);
        double get_optimal_exercise_price();
        
};

class BlackScholesSimulation
{
    public: 
        BlackScholesSimulation(
            const double S, const double mu, const double sigma, const double T, 
            const int number_iterations, const int sample_size, bool antithetic)
            : engine_(number_iterations, sample_size, std::make_shared<NormalDistribution>(0.0, 1.0), antithetic),
            S_(S), mu_(mu), sigma_(sigma), T_(T), is_updated_(false), simulated_prices(engine_.get_samples()), 
            time_taken_for_simulation_(0.0)
        {
            update();
        };
        ~BlackScholesSimulation() = default;

        double get_S() const {return S_;};
        double get_mu() const {return mu_;};
        double get_sigma() const {return sigma_;};
        double get_T() const {return T_;};
        bool is_updated() const {return is_updated_;};
        Eigen::MatrixXd get_simulated_prices() {if (!is_updated_) update(); return simulated_prices;};
        double get_time_taken_for_simulation() const {return time_taken_for_simulation_;};
        MonteCarloEngine get_simulation_engine() const {return engine_;};
        MonteCarloPayoffEngine get_payoff_engine() {if (!is_updated_) update(); return MonteCarloPayoffEngine(simulated_prices, T_);};

        void set_S(double S) {S_ = S; is_updated_=false;};
        void set_mu(double mu) {mu_ = mu; is_updated_=false;};
        void set_sigma(double sigma) {sigma_ = sigma; is_updated_=false;};
        void set_T(double T) {T_ = T; is_updated_=false;};
        void set_sample_size(int sample_size) {engine_.set_sample_dimension(sample_size, engine_.get_number_iterations()); is_updated_=false;};
        void set_number_iterations(int number_iterations) {engine_.set_sample_dimension(engine_.get_sample_size(), number_iterations); is_updated_=false;};

    private: 
        MonteCarloEngine engine_; 
        double S_; // Underlying asset price
        double mu_; // Risk-free interest rate
        double sigma_; // Volatility
        double T_; // Time to maturity
        bool is_updated_;
        Eigen::MatrixXd simulated_prices; 
        double time_taken_for_simulation_;
        void update();

};
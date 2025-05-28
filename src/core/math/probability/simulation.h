#pragma once 
#include <iostream>
#include "../../../../src/core/math/eigen_tools.h"
#include "../../../../src/core/math/probability/probability.h"
#include "stats.h"

class MonteCarloEngine
{
    public: 
        MonteCarloEngine(const int number_iterations, const int sample_size, const std::shared_ptr<ProbabilityDistribution>& distribution, bool antithetic)
            : M_(int(max(1,abs(sample_size)))), N_(get_corrected_number_iterations(number_iterations)), 
            distribution_(distribution), antithetic_(antithetic), updated_(false), samples_(Eigen::MatrixXd::Zero(N_, M_)), 
            time_taken_for_update_(0.0)
            {update_samples();};
        MonteCarloEngine(const int number_iterations, const std::shared_ptr<ProbabilityDistribution>& distribution, bool antithetic)
            : M_(1), N_(get_corrected_number_iterations(number_iterations)), distribution_(distribution),antithetic_(antithetic), 
            updated_(false), samples_(Eigen::MatrixXd::Zero(N_, M_)), time_taken_for_update_(0.0)
            {update_samples();};
        virtual ~MonteCarloEngine() = default;

        Eigen::MatrixXd& get_samples(){if (!updated_) { update_samples(); }; return samples_;}
        double get_sample_means() {if (!updated_) { update_samples(); }; return samples_.mean();}
        double get_sample_variances() {if (!updated_) { update_samples(); }; return samples_.array().square().mean() - samples_.mean() * samples_.mean();}
        int get_sample_size() const { return M_; }
        int get_number_iterations() const { return N_; }
        std::shared_ptr<ProbabilityDistribution> get_distribution() const { return distribution_; }
        bool is_updated() const { return updated_; }
        bool is_antithetic() const { return antithetic_; }
        double get_time_taken_for_samples_update() const { return time_taken_for_update_; }

        void set_sample_dimension(int sample_size, int number_iterations) { M_ = sample_size; N_=number_iterations, samples_ = Eigen::MatrixXd::Zero(N_, M_); updated_ = false;}
        void set_distribution(const std::shared_ptr<ProbabilityDistribution>& distribution) { distribution_ = distribution; updated_ = false; }
        void set_antithetic(bool antithetic) { antithetic_ = antithetic; updated_ = false;}
        void update_samples(); 
    
    protected:
        int get_corrected_number_iterations(int N) const {return max(int(antithetic_ ? (N % 2 != 0) ? ++N : N  : abs(N)), 1);};

    private: 
        int M_;
        int N_;
        std::shared_ptr<ProbabilityDistribution> distribution_;
        bool antithetic_; // Not used in this implementation, but can be added later
        bool updated_; 
        Eigen::MatrixXd samples_; 
        double time_taken_for_update_; // Not used in this implementation, but can be added later
        
};

class MultivariateMonteCarloEngine 
{
    public: 
        MultivariateMonteCarloEngine(const CorrelationMatrix& correlation_matrix, const int number_iterations, const int sample_size, const std::shared_ptr<ProbabilityDistribution>& distribution)
            : M_(int(max(1,abs(sample_size)))), N_(get_corrected_number_iterations(number_iterations)), distribution_(distribution), 
            correlation_matrix_(correlation_matrix), updated_(false), samples_({}), time_taken_for_update_(0.0)
            {update_samples();};
        virtual ~MultivariateMonteCarloEngine() = default;

        std::vector<Eigen::MatrixXd> get_samples(){if (!updated_) { update_samples(); }; return samples_;}
        std::vector<double> get_sample_means();
        std::vector<double> get_sample_variances();
        CorrelationMatrix get_correlation_matrix() const { return correlation_matrix_; };
        int get_sample_size() const { return M_; }
        int get_number_iterations() const { return N_; }
        bool is_updated() const { return updated_; }
        double get_time_taken_for_samples_update() const { return time_taken_for_update_; }
        //bool is_antithetic() const { return antithetic_; }
        std::shared_ptr<ProbabilityDistribution> get_distribution() const { return distribution_; }

        void set_sample_dimension(int sample_size, int number_iterations){ M_ = sample_size; N_=number_iterations, updated_ = false;}
        void set_distribution(const std::shared_ptr<ProbabilityDistribution>& distribution) {distribution_ = distribution; updated_ = false; }
        //void set_antithetic(bool antithetic) {antithetic_ = antithetic; updated_ = false;}
        void set_correlation_matrix(const CorrelationMatrix& correlation_matrix){correlation_matrix_ = correlation_matrix; updated_ = false;}
        void update_samples();
    
    protected: 
        int get_corrected_number_iterations(int N) const {return max(abs(N), 1);};

    private: 
        int M_; 
        int N_; 
        std::shared_ptr<ProbabilityDistribution> distribution_;
        //bool antithetic_; // Not used in this implementation, but can be added later
        CorrelationMatrix correlation_matrix_;
        bool updated_; 
        std::vector<Eigen::MatrixXd> samples_; 
        double time_taken_for_update_;
        void update_samples_non_antithetic();
        //void update_samples_antithetic();
};


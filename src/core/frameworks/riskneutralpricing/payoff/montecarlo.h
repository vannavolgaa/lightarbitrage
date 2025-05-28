#pragma once 
#include <iostream>
#include "../../../../../src/core/math/eigen_tools.h"

class MonteCarloPayoffEngine
{
    public: 
        MonteCarloPayoffEngine(std::vector<std::vector<double>> prices, const double T): prices_(get_egein_matrix_from_vectors(prices)), T_(T){}; // Default constructor is not allowed
        MonteCarloPayoffEngine(Eigen::MatrixXd prices, const double T): prices_(prices), T_(T) {};
        ~MonteCarloPayoffEngine() = default;

        double european_vanilla_price(double r, double K, bool is_call) const; 

        
    protected: 
        Eigen::VectorXd get_terminal_prices() const {return prices_.col(prices_.cols() - 1);}

    private: 
        Eigen::MatrixXd prices_;
        double T_; 
        
};
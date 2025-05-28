#pragma once 
#include <iostream>
#include <cfloat>
#include <tuple>
#include "../../../../../src/core/frameworks/riskneutralpricing/normblack/normblack.h"
#include "../../../../../src/core/math/optimization/optimization.h"
#include "../../../../../src/errors/core/frameworks.h"

namespace LetsBeRational
{
    
    double get_rational_cubic_interpolate(double x, double x0, double x1,double y0, double y1,double dy0, double dy1,double r); 

    double get_r_right_side(double x0, double x1,double y0, double y1,double dy0, double dy1, double ddy1);

    double get_r_left_side(double x0, double x1,double y0, double y1,double dy0, double dy1, double ddy0);

    double get_right_asymptotic(double sigma); 

    double get_right_asymptotic_first_derivative(double x, double sigma);

    double get_right_asymptotic_second_derivative(double x, double sigma);

    double get_left_asymptotic_z(double x, double sigma); 

    double get_left_asymptotic(double x, double sigma); 

    double get_left_asymptotic_first_derivative(double x, double sigma);

    double get_left_asymptotic_second_derivative(double x, double sigma);

    std::tuple<double, double, double, double, double, double, double, double, double, double, bool> get_initial_data(double beta, double x, double is_call);

    double get_initial_guess_normalized_volatility(double beta, double x, bool is_call); 

    std::function<double(double)> get_target_function(double beta, double x, double b_l, double b_u); 

    std::function<double(double)> get_target_function_first_derivative(double beta, double x, double b_l, double b_u); 

    double get_newton_normalized_volatility(double beta, double x, bool is_call);

    double get_newton_black_scholes_implied_volatility(double price, double S, double K, double T, double r, double q, bool is_call, bool is_future); 

};

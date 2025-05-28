#pragma once 
#include <iostream>
#include "../../../../../src/core/math/probability/probability.h"
#include "../../../../../src/core/math/erf_cody.h"
#include "../../../../../src/core/math/quadratures.h"

using namespace constants; 

// Simplified version of the Black Price (undisounted and normalized). This frameworks is presented in "Let's be rational" by Mark Joshi.
// Some code written here has been directly sourced the official github (https://github.com/vollib/lets_be_rational/)
// From this framework we also provide the Lewis Formula expressed in undiscounted and normalized manner for further frameworks development. 

inline int get_future_flag(bool is_future){return is_future ? 0 : 1;}

inline int get_call_put_flag(bool is_call){return is_call ? 1 : -1;}

inline NormalDistribution stdnorm(0.0,1.0);

double get_lewis_normalized_black_price_gauss_laguerre(
    double x, bool is_call,
    std::function<std::complex<double>(std::complex<double>)> cf, 
    std::shared_ptr<GaussianQuadrature> quadrature);

double get_lewis_normalized_black_price_gauss_laguerre(
    double x, bool is_call,
    std::function<std::complex<double>(std::complex<double>)> cf, 
    int points);

namespace NormalizedBlack
{
    constexpr double H_LARGE = -10.0;
    constexpr double T_SMALL = 0.21;
    // Code source from https://github.com/vollib/lets_be_rational/
    double get_call_price_region_1(double h, double t); 
    // Code source from https://github.com/vollib/lets_be_rational/
    double get_call_price_region_2(double h, double t); 
    double get_call_price_region_4(double h, double t); 
    double get_call_price_region_3(double h, double t); 
    double get_intrisic(double x, bool is_call);
    double get_call_price(double x, double normalized_sigma);
    double get_vega(double x, double normalized_sigma);
    double get_volga(double x, double normalized_sigma);
    double get_price(double x, double normalized_sigma, bool is_call);
    std::complex<double> get_characteristic_function(double u, double x, double normalized_sigma);
    double get_lewis_price(double x, double normalized_sigma, bool is_call, std::shared_ptr<GaussianQuadrature> quadrature);

}

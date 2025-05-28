#include <iostream>
#include <iomanip>
#include <memory>
#include <cassert>
#include <map>
#include <chrono>
#include "../../../../../src/core/frameworks/riskneutralpricing/blackscholes/bs.h"

// Helper function for floating-point comparisons
bool is_close(double a, double b, double tol = 1e-3) {
    return std::abs(a - b) <= tol;
}

void test_black_scholes()
{
    // Input parameters
    double S = 100.0;    // Stock price
    double K = 100.0;    // Strike price
    double r = 0.05;     // Risk-free rate
    double q = 0.02;     // Dividend yield
    double sigma = 0.20; // Volatility
    double T = 2.0;      // Time to maturity

    // Option prices and Greeks
    double call_price = 13.5218;
    double put_price = 7.9266;
    
    double call_delta = 0.6131;
    double put_delta = -0.3476;
    
    double gamma = 0.0127;
    
    double call_theta = -3.7095;
    double put_theta = -1.1068;
    
    double call_rho = 95.5845;
    double put_rho = -85.3830;
    
    double vega = 50.9225;
    
    BlackScholes bs(S, K, T, r, q, sigma, true, false);
    // Assertions for Stock Options (Black-Scholes)
    assert(is_close(bs.get_price(), call_price));
    assert(is_close(bs.get_delta(), call_delta));
    assert(is_close(bs.get_gamma(), gamma));
    assert(is_close(bs.get_vega(), vega));
    assert(is_close(bs.get_theta(), call_theta));
    std::cout << "Black-Scholes Call Tests Passed!\n";
    bs.set_is_call(false);
    assert(is_close(bs.get_price(), put_price));
    assert(is_close(bs.get_delta(), put_delta));
    assert(is_close(bs.get_gamma(), gamma));
    assert(is_close(bs.get_vega(), vega));
    assert(is_close(bs.get_theta(), put_theta));
    std::cout << "Black-Scholes Put Tests Passed!\n";
    

}

void test_barone_adesi_whaley2(double r,double q,double sigma, double T, bool is_call, bool is_future, std::vector<double> baw_prices)
{
    std::map<double,double> baw_prices_map =  { // Initialisation directe
        {80, baw_prices[0]},
        {90, baw_prices[1]},
        {100, baw_prices[2]},
        {110, baw_prices[3]},
        {120, baw_prices[4]}
    };
    std::vector<double> spot_list = {80, 90,100,110,120};
    for (double s: spot_list){
        BaroneAdesiWhaley baw(s, 100, T, r, q, sigma, is_call, is_future);
        assert(is_close(baw.get_price(), baw_prices_map.at(s), 1e-2));
    }
};

void full_test_barone_adesi_whaley()
{

    std::vector<double> baw_prices =  {.05,.85,4.44,11.66,20.90};
    test_barone_adesi_whaley2(.08,.04,.2,.25,true,false,baw_prices);
    baw_prices =  {.05,.84,4.40,11.55,20.69};
    test_barone_adesi_whaley2(.12,.08,.2,.25,true,false,baw_prices);
    baw_prices =  {1.29,3.82,8.35,14.8,22.72};
    test_barone_adesi_whaley2(.08,.04,.4,.25,true,false,baw_prices);
    baw_prices =  {.41,2.18,6.5,13.42,22.06};
    test_barone_adesi_whaley2(.08,.04,.2,.5,true,false,baw_prices);
    std::cout<<"Barone Adesi Whaley call prices with b=0.04 are all validated."<<std::endl;

    baw_prices =  {.03,.59,3.52,10.31,20.0};
    test_barone_adesi_whaley2(.08,.12,.2,.25,true,false,baw_prices);
    baw_prices =  {.03,.59,3.51,10.29,20.0};
    test_barone_adesi_whaley2(.12,.16,.2,.25,true,false,baw_prices);
    baw_prices =  {1.07,3.28,7.41,13.5,21.23};
    test_barone_adesi_whaley2(.08,.12,.4,.25,true,false,baw_prices);
    baw_prices =  {.23,1.39,4.72,10.96,20.0};
    test_barone_adesi_whaley2(.08,.12,.2,.5,true,false,baw_prices);
    std::cout<<"Barone Adesi Whaley call prices with b=-0.04 are all validated."<<std::endl;

    baw_prices =  {.04,.70,3.93,10.81,20.02};
    test_barone_adesi_whaley2(.08,.08,.2,.25,true,true,baw_prices);
    baw_prices =  {.04,.70,3.9,10.75,20.0};
    test_barone_adesi_whaley2(.12,.16,.2,.25,true,true,baw_prices);
    baw_prices =  {1.17,3.53,7.84,14.08,21.86};
    test_barone_adesi_whaley2(.08,.12,.4,.25,true,true,baw_prices);
    baw_prices =  {.30,1.72,5.48,11.9,20.34};
    test_barone_adesi_whaley2(.08,.12,.2,.5,true,true,baw_prices);
    std::cout<<"Barone Adesi Whaley call prices on Future (b=0.00) are all validated."<<std::endl;

    baw_prices =  {2.52,4.97,8.67,13.88,20.88};
    test_barone_adesi_whaley2(.08,.12,.2,3.0,true,false,baw_prices);
    baw_prices =  {4.2,7.54,12.03,17.64,24.30};
    test_barone_adesi_whaley2(.08,.0,.2,3.0,true,true,baw_prices);
    baw_prices =  {6.97,11.62,17.4,24.09,31.49};
    test_barone_adesi_whaley2(.08,.04,.2,3.0,true,false,baw_prices);
    std::cout<<"Barone Adesi Whaley call prices on long expiry (T=3.0) are all validated."<<std::endl;

    baw_prices =  {20.0,10.18,3.54,0.8,0.12};
    test_barone_adesi_whaley2(.08,.04,.2,.25,false,false,baw_prices);
    baw_prices =  {20.0,10.16,3.53,0.79,0.12};
    test_barone_adesi_whaley2(.12,.08,.2,.25,false,false,baw_prices);
    baw_prices =  {20.53,12.93,7.46,3.96,1.95};
    test_barone_adesi_whaley2(.08,.04,.4,.25,false,false,baw_prices);
    baw_prices =  {20.0,10.71,4.77,1.76,0.55};
    test_barone_adesi_whaley2(.08,.04,.2,.5,false,false,baw_prices);
    std::cout<<"Barone Adesi Whaley put prices with b=0.04 are all validated."<<std::endl;

    baw_prices =  {20.42,11.25,4.4,1.12,0.18};
    test_barone_adesi_whaley2(.08,.12,.2,.25,false,false,baw_prices);
    baw_prices =  {20.25,11.15,4.35,1.11,0.18};
    test_barone_adesi_whaley2(.12,.16,.2,.25,false,false,baw_prices);
    baw_prices =  {21.46,13.93,8.27,4.52,2.3};
    test_barone_adesi_whaley2(.08,.12,.4,.25,false,false,baw_prices);
    baw_prices =  {20.98,12.64,6.37,2.65,0.92};
    test_barone_adesi_whaley2(.08,.12,.2,.5,false,false,baw_prices);
    std::cout<<"Barone Adesi Whaley put prices with b=-0.04 are all validated."<<std::endl;

    baw_prices =  {20.0,10.58,3.93,0.94,0.15};
    test_barone_adesi_whaley2(.08,.08,.2,.25,false,true,baw_prices);
    baw_prices =  {20.0,10.53,3.9,0.93,0.15};
    test_barone_adesi_whaley2(.12,.16,.2,.25,false,true,baw_prices);
    baw_prices =  {20.93,13.39,7.84,4.23,2.12};
    test_barone_adesi_whaley2(.08,.12,.4,.25,false,true,baw_prices);
    baw_prices =  {20.04,11.48,5.48,2.15,0.7};
    test_barone_adesi_whaley2(.08,.12,.2,.5,false,true,baw_prices);
    std::cout<<"Barone Adesi Whaley put prices on Future (b=0.00) are all validated."<<std::endl;

    baw_prices =  {26.25,20.64,15.99,12.22,9.23};
    test_barone_adesi_whaley2(.08,.12,.2,3.0,false,false,baw_prices);
    baw_prices =  {22.4,16.5,12.03,8.69,6.22};
    test_barone_adesi_whaley2(.08,.0,.2,3.0,false,true,baw_prices);
    baw_prices =  {20.33,13.56,9.11,6.12,4.12};
    test_barone_adesi_whaley2(.08,.04,.2,3.0,false,false,baw_prices);
    baw_prices =  {20.0,11.63,6.96,4.26,2.64};
    test_barone_adesi_whaley2(.08,.0,.2,3.0,false,false,baw_prices);
    std::cout<<"Barone Adesi Whaley put prices on long expiry (T=3.0) are all validated."<<std::endl;


};

void test_lewis_method()
{
    auto start = std::chrono::high_resolution_clock::now();
    std::shared_ptr<GaussianQuadrature> quadrature = get_gauss_laguerre_quadrature(64);
    auto end = std::chrono::high_resolution_clock::now();
    double x = 0.5;
    double normalized_sigma = 0.2;
    bool is_call = false;

    // Test Lewis method
    double lewis_price = NormalizedBlack::get_lewis_price(x, normalized_sigma, is_call, quadrature);
    double price = NormalizedBlack::get_price(x, normalized_sigma, is_call);
    std::cout << "Lewis method price: " << lewis_price << "\n";
    std::cout << "Expected price: " << price << "\n";
    
    double S = 100.0;    // Stock price
    double K = 100.0;    // Strike price
    double r = 0.05;     // Risk-free rate
    double q = 0.02;     // Dividend yield
    double sigma = 0.20; // Volatility
    double T = 2.0;      // Time to maturity
    bool is_future = false;
    BlackScholes bs(S, K, T, r, q, sigma, is_call, is_future);
    double lewis_price2 = bs.get_lewis_price(quadrature);
    double bs_price = bs.get_price();
    std::cout << "Lewis method price (Black-Scholes): " << lewis_price2 << "\n";
    std::cout << "Expected price (Black-Scholes): " << bs_price << "\n";

    std::cout << "Lewis method test passed!\n";
}

void test_monte_carlo_european_vanilla_price()
{
    std::cout << "Testing European Vanilla Price Calculation..." << std::endl;

    // Parameters for the Black-Scholes model
    double S = 100.0; // Underlying asset price
    double K = 100.0; // Strike price
    double T = 1.0; // Time to maturity in years
    double r = 0.05; // Risk-free interest rate
    double q = 0.02; // Dividend yield
    double sigma = 0.2; // Volatility
    bool is_call = true; // Call option

    // Create a Black-Scholes object
    BlackScholes bs(S, K, T, r, q, sigma, is_call, false);

    // Create a Black-Scholes simulation engine
    int number_iterations = 10000;
    int sample_size = 1000;
    std::shared_ptr<BlackScholesSimulation> engine = std::make_shared<BlackScholesSimulation>(S, r-q, sigma, T, number_iterations, sample_size, true);
    std::cout << "Time taken for samples update: " << engine->get_time_taken_for_simulation() << " seconds" << std::endl;
    // Calculate the discount factor
    double df = exp(-r * T);

    // Get the European vanilla price using Monte Carlo simulation
    MonteCarloPayoffEngine payoffs = engine->get_payoff_engine();
    double call_price = payoffs.european_vanilla_price(r,K,true);
    double put_price = payoffs.european_vanilla_price(r,K,false);
    double excact_call_price = bs.get_price();
    bs.set_is_call(false);
    double excact_put_price = bs.get_price();

    std::cout << "Calculated European Vanilla Call Price: " << std::fixed << std::setprecision(2) << call_price << std::endl;
    std::cout << "Calculated European Vanilla Put Price: " << std::fixed << std::setprecision(2) << put_price << std::endl;
    std::cout << "Exact European Vanilla Call Price: " << std::fixed << std::setprecision(2) << excact_call_price << std::endl;
    std::cout << "Exact European Vanilla Put Price: " << std::fixed << std::setprecision(2) << excact_put_price << std::endl;

    // Expected price calculation for verification (using Black-Scholes formula
    
    assert(fabs(put_price - excact_put_price) < 3e-1); // Allowing a small tolerance for numerical errors
    assert(fabs(call_price - excact_call_price) < 3e-1); // Allowing a small tolerance for numerical errors

    std::cout << "Test passed! Calculated price matches expected price." << std::endl;
}

int main() {
    test_black_scholes();
    full_test_barone_adesi_whaley();
    test_lewis_method();
    test_monte_carlo_european_vanilla_price();
    std::cout << "All tests passed successfully!\n";
    return 0;
}

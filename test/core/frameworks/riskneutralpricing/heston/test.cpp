#include <iostream>
#include <iomanip>
#include <memory>
#include <cassert>
#include <map>
#include <chrono>
#include "../../../../../src/core/frameworks/riskneutralpricing/heston/heston.h" 

// Helper function for floating-point comparisons
bool is_close(double a, double b, double tol = 1e-3) {
    return std::abs(a - b) <= tol;
}

void test_heston_lewis_price_function()
{
    // Test Heston model
    double S = 100.0;    // Stock price
    double K = 100.0;    // Strike price
    double r = 0.05;     // Risk-free rate
    double q = 0.01;     // Dividend yield
    double T = 1.5;      // Time to maturity
    double F = S * std::exp((r - q) * T); // Forward price

    Heston heston(2.0,0.05,0.3,0.45,0.05);
    heston.set_quadrature(64); // Set the number of quadrature points

    double heston_price = heston.get_lewis_price(F, K, T,r,true);
    std::cout << "Heston model price: " << heston_price << "\n";
    std::cout << "Expected: " << 13.2561 << "\n";
    assert(is_close(heston_price, 13.2561, 1e-4));

}

void test_heston_simulation()
{
    double S = 100.0; 
    double K = 100.0;    // Strike price
    double r = 0.05;     // Risk-free rate
    double q = 0.01;     // Stock price
    double mu = r-q;    // Drift
    double sigma = 0.2;  // Volatility
    double T = 1.0;      // Time to maturity
    double kappa = 2.0;  // Speed of mean reversion
    double theta = 0.05; // Long-term variance
    double eta = 0.3;    // Volatility of volatility
    double rho = 0;  // Correlation between asset and volatility processes
    double v0 = 0.05;    // Initial variance
    double F = S * std::exp((r - q) * T); // Forward price
    double df = std::exp(-r * T); // Discount factor
    std::shared_ptr<HestonSimulation> engine = std::make_shared<HestonSimulation>(S, mu, T, kappa, theta, eta, rho, v0, 1000, 5000);
    Heston model = Heston(kappa, theta, eta, rho, v0);
    MonteCarloPayoffEngine payoff = engine->get_payoff_engine();
    double expected_call_price = model.get_lewis_price(F, K, T, r, true);
    double expected_put_price = model.get_lewis_price(F, K, T, r, false);
    double call_price = payoff.european_vanilla_price(r, K, true);
    double put_price = payoff.european_vanilla_price(r, K, false);
    std::cout << "Heston Simulation Call Price: " << call_price << "\n";
    std::cout << "Expected Call Price: " << expected_call_price << "\n";
    std::cout << "Heston Simulation Put Price: " << put_price << "\n";
    std::cout << "Expected Put Price: " << expected_put_price << "\n";

    std::cout << "Simulation completed successfully.\n";
}

int main() {
    test_heston_lewis_price_function(); 
    test_heston_simulation();
    std::cout << "All tests passed successfully!\n";
    return 0;
}

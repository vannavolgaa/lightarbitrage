#include "montecarlo.h"

double MonteCarloPayoffEngine::european_vanilla_price(double r, double K, bool is_call) const
{
    double df = std::exp(-r * T_); // Discount factor
    Eigen::VectorXd terminal_prices = get_terminal_prices(); 
    double sum = 0.0;
    int M = terminal_prices.size(); // Number of samples
    int c = is_call ? 1 : -1; // Call or Put option
    for (int i = 0; i < M; ++i)
    {
        sum += std::max(c * (terminal_prices(i) - K), 0.0);
    }
    return df * sum / double(M);
}
#include "heston.h"

std::complex<double> Heston::get_characteristic_function(std::complex<double> u, double T)
{
    std::complex<double> i(0.0, 1.0);
    double e = eta_, k = kappa_, ro = rho_;
    std::complex<double> beta = k - ro*e*i*u;
    std::complex<double> d = sqrt(beta*beta + e*e*(i*u + u*u));
    std::complex<double> g = (beta - d)/(beta + d); 
    std::complex<double> D = (beta - d)*(1.0-exp(-d*T))/(e*e*(1.0-g*exp(-d*T)));
    std::complex<double> C = k*(T*(beta-d) - 2.0*log((1.0-g*exp(-d*T))/(1.0-g)))/(e*e); 
    return exp(C*theta_+D*v0_); 
}

double Heston::get_lewis_normalized_price(const double x, const double T, const bool is_call)
{
    std::function<std::complex<double>(std::complex<double>)> cf = [this, T](std::complex<double> u) {
        return get_characteristic_function(u, T);
    };
    return get_lewis_normalized_black_price_gauss_laguerre(x, is_call, cf, quadrature_);
}

double Heston::get_lewis_price(const double F, const double K, const double T, const double r, const bool is_call)
{
    double x = log(F/K);
    return exp(-r*T)*get_lewis_normalized_price(x, T, is_call)*sqrt(F*K);
};

void HestonSimulation::update_euler_milstein(bool milstein)
{
    auto start = std::chrono::high_resolution_clock::now();
    double dt = T_ / double(engine_.get_sample_size());
    std::vector<Eigen::MatrixXd> Z = engine_.get_samples();
    Eigen::MatrixXd Zv = Z[0];
    Eigen::MatrixXd Zs = Z[1];
    double v0 = model_.get_v0(), k = model_.get_kappa(), th = model_.get_theta(), eta = model_.get_eta();
    for (int i = 0; i < Zv.rows(); ++i) {
        for (int j = 0; j < Zv.cols(); ++j) {
            double milstein_term = milstein ? 0.25 * eta * eta * dt * (Zv(i, j) * Zv(i, j) - 1.0) : 0.0;
            if (j==0) 
            {
                double v = get_variance_from_truncation(v0 + k * (th - v0) * dt + eta * sqrt(dt*v0) * Zv(i, j) + milstein_term);
                simulated_variances(i, j) = v;
                simulated_prices(i, j) = S_ * exp((mu_ - 0.5 * v0) * dt + sqrt(v0*dt) * Zs(i, j));
            }
            else
            {
                double v_1 = simulated_variances(i, j-1);
                double v = get_variance_from_truncation(v_1 + k * (th - v_1) * dt + eta * sqrt(dt*v_1) * Zv(i, j) + milstein_term);
                simulated_variances(i, j) = v;
                simulated_prices(i, j) = simulated_prices(i, j-1) * exp((mu_ - 0.5 * v_1) * dt + sqrt(v_1*dt) * Zs(i, j));
            }
        }
    }
    auto end = std::chrono::high_resolution_clock::now();
    time_taken_for_simulation_ = std::chrono::duration<double>(end - start).count();
    is_updated_ = true;
}

void HestonSimulation::update()
{
    switch (discretization_method_)
    {
        case HestonDiscretizationMethod::EULER:
            update_euler_milstein(false);
            break;
        case HestonDiscretizationMethod::MILSTEIN:
            update_euler_milstein(true);
            break;
        default:
            throw FrameworksError(FrameworksErrorCode::HESTON_UNKNOWN_DISCRETIZATION_METHOD);
    }
}
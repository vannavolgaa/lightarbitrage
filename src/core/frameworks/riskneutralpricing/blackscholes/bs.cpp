#include "bs.h"

double BlackScholes::get_future_price()  
{
    check_update();
    return F_; 
}; 

void BlackScholes::update() 
{
    if (S_ < 0 || K_ <= 0 || T_ <= 0 || sigma_ <= 0) {
        throw FrameworksError(FrameworksErrorCode::BLACK_SCHOLES_MODEL_INVALID_PARAM);
    }

    F_ = S_ * exp(future_flag_*(r_ - q_) * T_);
    d1_ = (log(F_ / K_) + 0.5 * sigma_ * sigma_ * T_) / (sigma_ * sqrt(T_)); 
    d2_ = d1_ - sigma_ * sqrt(T_);
    Nd1_ = stdnorm.cdf(call_put_flag_*d1_);
    Nd2_ = stdnorm.cdf(call_put_flag_*d2_);
    nd1_ = stdnorm.pdf(d1_);
    nd2_ = stdnorm.pdf(d2_);
    
    is_updated_ = true;

};

double BlackScholes::get_price() 
{
    check_update(); 
    return exp(-r_*T_)*NormalizedBlack::get_price(log(F_ / K_), sigma_*sqrt(T_), get_is_call())*sqrt(F_*K_);
};

double BlackScholes::get_dual_delta()
{
    check_update();
    return call_put_flag_*exp(-r_*T_)*Nd2_;
}

double BlackScholes::get_dual_gamma()
{
    check_update();
    return exp(-r_*T_)*nd2_/(K_*sigma_*sqrt(T_));
}

double BlackScholes::get_delta()
{
    check_update();
    double mu = future_flag_*(r_-q_);
    return exp(-r_*T_)*call_put_flag_*exp(mu*T_)*Nd1_;
}

double BlackScholes::get_gamma()
{
    check_update();
    double mu = future_flag_*(r_-q_);
    return exp(-r_*T_)*exp(mu*T_)*exp(mu*T_)*nd1_/(F_*sigma_*sqrt(T_));
}

double BlackScholes::get_vega()
{
    check_update();
    return F_*exp(-r_*T_)*nd1_*sqrt(T_);
}

double BlackScholes::get_theta()
{
    check_update();
    double mu = future_flag_*(r_-q_);
    double term1 = -F_*exp(-r_*T_)*nd1_*sigma_/(2*sqrt(T_));
    double term2 = -call_put_flag_*r_*K_*exp(-r_*T_)*Nd2_;
    double term3 = call_put_flag_*(r_-mu)*F_*exp(-r_*T_)*Nd1_;
    return term1+term2+term3;
}

double BlackScholes::get_lewis_price(std::shared_ptr<GaussianQuadrature> quadrature)
{
    check_update();
    double x = log(F_/K_);
    return exp(-r_*T_)*NormalizedBlack::get_lewis_price(x, sigma_*sqrt(T_), get_is_call(), quadrature)*sqrt(F_*K_);
};

void BaroneAdesiWhaley::update()
{
    if (S_ < 0 || K_ <= 0 || T_ <= 0 || sigma_ <= 0) {
        throw FrameworksError(FrameworksErrorCode::BARONE_ADESI_WHALEY_MODEL_INVALID_PARAM);
    }
    double mu = future_flag_*(r_-q_);
    double N = 2*mu/(sigma_*sigma_);
    double M = 2*r_/(sigma_*sigma_);
    qq_inf_ = .5*(1-N+call_put_flag_*sqrt((N-1)*(N-1)+4*M));
    qq_ = .5*(1-N+call_put_flag_*sqrt((N-1)*(N-1)+4*M/(1-exp(-r_*T_))));
    double Sinf = K_/(1 - 1/qq_inf_);
    double payoff = call_put_flag_*(Sinf - K_); 
    double h = -call_put_flag_*(mu*T_+call_put_flag_*2*sigma_*sqrt(T_))*K_/payoff; 
    Se0_ = Sinf-call_put_flag_*exp(h)*payoff;
    Se_ = get_optimal_exercise_price();
    is_updated_ = true;
}; 

BlackScholes BaroneAdesiWhaley::get_bs(double s) const
{
    return BlackScholes(s, K_, T_, r_, q_, sigma_, get_is_call(), get_is_future());
};

double BaroneAdesiWhaley::target_function(const double s)
{
    BlackScholes bs_ = get_bs(s);
    double bs_price = bs_.get_price();
    double bs_delta = bs_.get_delta();
    return call_put_flag_*(s-K_) - bs_price - call_put_flag_*s*(1-call_put_flag_*bs_delta)/qq_;
};

double BaroneAdesiWhaley::derivative_target_function(const double s)
{
    BlackScholes bs_ = get_bs(s);
    double bs_delta = bs_.get_delta();
    double bs_gamma = bs_.get_gamma();
    return call_put_flag_ - bs_delta*(1-1/qq_) - (call_put_flag_-bs_gamma*s)/qq_;
};

double BaroneAdesiWhaley::get_optimal_exercise_price()
{
    if (future_flag_*(r_-q_) >= r_ and get_is_call()){return S_*S_*S_;}
    NewtonRaphson nr = NewtonRaphson(Se0_, [this](double s){return target_function(s);}, [this](double s){return derivative_target_function(s);});
    nr.set_max_iterations(max_iter_);
    nr.set_tolerance_rate(epsilon_);
    std::shared_ptr<OptimizerResult> optim = nr.optimize();
    if (optim->number_iterations_ == optim->maximum_iterations_){
        if(abs(target_function(Se0_))<abs(optim->function_value_))return Se0_; 
        else return optim->x_[0] ;
    }
    else return optim->x_[0];
};

double BaroneAdesiWhaley::get_exercise_premium()
{
    check_update();
    BlackScholes bs_ = get_bs(Se_);
    double bs_delta = bs_.get_delta();
    return call_put_flag_ * Se_ * std::pow(S_ / Se_, qq_) * (1 - call_put_flag_ * bs_delta) / qq_;
};

double BaroneAdesiWhaley::get_price()
{
    check_update();
    if (call_put_flag_*S_ >= call_put_flag_*Se_) return call_put_flag_*(S_-K_);
    BlackScholes bs_ = get_bs(S_);
    return bs_.get_price() + get_exercise_premium();
    
};

void BlackScholesSimulation::update()
{
    auto start_time = std::chrono::high_resolution_clock::now();
    if (S_ < 0 || sigma_ <= 0 || T_ <= 0) {
        throw FrameworksError(FrameworksErrorCode::BLACK_SCHOLES_MODEL_INVALID_PARAM);
    }
    double dt = T_ / double(engine_.get_sample_size());
    double drift = mu_ - 0.5 * sigma_ * sigma_;
    double diffusion = sigma_ * sqrt(dt);
    Eigen::MatrixXd Z = engine_.get_samples();
    for (int i = 0; i < Z.rows(); ++i) {
        for (int j = 0; j < Z.cols(); ++j) {
            double moneyness = exp(drift * dt + diffusion * Z(i,j));
            if (j==0) simulated_prices(i, j) = S_ * moneyness; 
            else simulated_prices(i, j) = simulated_prices(i, j-1) * moneyness; 
        }
    }
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_time = end_time - start_time;
    time_taken_for_simulation_ = elapsed_time.count();
    is_updated_ = true;
}
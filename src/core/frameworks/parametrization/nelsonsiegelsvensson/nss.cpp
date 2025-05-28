#include "nss.h"


double NelsonSiegel::get_function1(double t) const {
    return (1 - exp(-t / tau_)) / (t / tau_);
}

double NelsonSiegel::get_function2(double t) const {
    return get_function1(t) - exp(-t / tau_);
}

double NelsonSiegel::get_rate(const double& t) const {
    return b0_ + b1_ * get_function1(t) + b2_ * get_function2(t);
}

double NelsonSiegel::get_forward_rate(const double& t) const {
    return b0_ + b1_ * exp(-t / tau_) + b2_ * exp(-t / tau_) *  t / tau_;
}

double NelsonSiegelSvensson::get_function3(double t) const {
    double t1 = t / tau2_;
    return (1 - exp(-t1)) / t1 - exp(-t1);;
}

double NelsonSiegelSvensson::get_rate(const double& t) const {
    
    return NelsonSiegel::get_rate(t) + b3_ * get_function3(t);
}

double NelsonSiegelSvensson::get_forward_rate(const double& t) const {
    double t1 = t / tau2_;
    return NelsonSiegel::get_forward_rate(t) + b3_ * exp(-t1) * t / tau2_; 
}
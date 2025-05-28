#pragma once 
#include <iostream>

class NelsonSiegel 
{
    public: 
        NelsonSiegel(double b0, double b1, double b2, double tau) 
            : b0_(b0), b1_(b1), b2_(b2), tau_(tau) {}
        double get_rate(const double& t) const; 
        double get_forward_rate(const double& t) const;
        double get_b1() const { return b1_; }
        double get_b2() const { return b2_; }
        double get_b0() const { return b0_; }
        double get_tau() const { return tau_; }
        void set_b0(double b0) { b0_ = b0; }
        void set_b1(double b1) { b1_ = b1; }
        void set_b2(double b2) { b2_ = b2; }
        void set_tau(double tau) { tau_ = tau; }
    private:
        double b0_;  
        double b1_;  
        double b2_;  
        double tau_;
        double get_function1(double t) const; 
        double get_function2(double t) const; 
}; 

class NelsonSiegelSvensson : public NelsonSiegel
{
    public: 
        NelsonSiegelSvensson(double b0, double b1, double b2, double b3, double tau1, double tau2) 
            : NelsonSiegel(b0, b1, b2, tau1), b3_(b3), tau2_(tau2) {}
        double get_rate(const double& t) const; 
        double get_forward_rate(const double& t) const;
        double get_b3() const { return b3_; }
        double get_tau2() const { return tau2_; }
        void set_b3(double b3) { b3_ = b3; }
        void set_tau2(double tau2) { tau2_ = tau2; }
    private:
        double b3_; 
        double tau2_;
        double get_function3(double t) const;
};
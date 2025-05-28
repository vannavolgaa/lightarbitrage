#pragma once
#include <iostream>
#include <random>
#include <complex>
#include <cfloat>
#include "../../../../src/core/math/constants.h"
#include "../../../../src/core/math/erf_cody.h"
#include "../../../../src/errors/core/math.h"

using namespace constants;

class ProbabilityDistribution
{
    public: 
        virtual double inv_cdf(double p) const = 0;
        virtual double cdf(double x) const = 0; 
        virtual double mgf(double t) const = 0; 
        virtual std::complex<double> cf(double t) const = 0; 
        virtual double r() = 0; 
        ProbabilityDistribution(){}; 
        virtual ~ProbabilityDistribution()=default;
}; 

class DiscreteProbabilityDistribution : public ProbabilityDistribution
{
    public: 
        virtual double inv_cdf(double p) const = 0;
        virtual double cdf(double x) const = 0; 
        virtual double mgf(double t) const = 0; 
        virtual std::complex<double> cf(double t) const = 0; 
        virtual double r() = 0; 
        virtual double pmf(double x) const = 0; 
        DiscreteProbabilityDistribution(){}; 
        virtual ~DiscreteProbabilityDistribution()=default;
}; 

class ContinuousProbabilityDistribution : public ProbabilityDistribution
{
    public: 
        virtual double inv_cdf(double p) const = 0;
        virtual double cdf(double x) const = 0; 
        virtual double mgf(double t) const = 0; 
        virtual std::complex<double> cf(double t) const = 0; 
        virtual double r() = 0; 
        virtual double pdf(double x) const = 0; 
        ContinuousProbabilityDistribution(){}; 
        virtual ~ContinuousProbabilityDistribution()=default;
}; 

class UniformDistribution final: public ContinuousProbabilityDistribution
{
    public: 
        UniformDistribution(double left_bound, double right_bound): 
        a_(left_bound), b_(right_bound), dist_(left_bound, right_bound), gen_(std::random_device{}())   
        {
            check_bounds(a_, b_);
        };
        UniformDistribution(): a_(0.0), b_(1.0), dist_(0.0, 1.0), gen_(std::random_device{}()){};
        virtual ~UniformDistribution()=default;

        double get_left_bound() const {return a_;}; 
        double get_right_bound() const {return b_;};
        double inv_cdf(double p) const override {return a_ + (b_ - a_) * p;}
        double cdf(double x) const override;
        double mgf(double t) const override;
        std::complex<double> cf(double t) const override;
        double pdf(double x) const override; 
        double r() override;

        void set_bounds(double left_bound, double right_bound){a_ = left_bound; b_=right_bound; check_bounds(a_, b_);}; 

    private: 
        double a_; 
        double b_; 
        std::uniform_real_distribution<double> dist_;
        std::mt19937 gen_;
        void check_bounds(double a, double b) const {if (a >= b) throw MathError(MathErrorCode::UNIFORM_INVALID_BOUNDS);}

};

class NormalDistribution final: public ContinuousProbabilityDistribution
{
    public: 
        NormalDistribution(double mu, double sigma): mu_(mu), sigma_(sigma)
        {
            if (sigma<=0) throw MathError(MathErrorCode::NORMAL_DISTRIBUTION_INVALID_SIGMA);
        };
        NormalDistribution(): mu_(0.0), sigma_(1.0){};
        virtual ~NormalDistribution()=default;

        double get_sigma() const {return sigma_;}; 
        double get_mu() const {return mu_;};
        double cdf(double x) const override;
        double inv_cdf(double p) const override {return mu_ + sigma_ * acklam_standard_inv_cdf(p);}
        double mgf(double t) const override{ return exp(mu_*t + sigma_*sigma_*t*t/2);};
        std::complex<double> cf(double t) const override;
        double pdf(double x) const override; 
        double r() override; 

        void set_sigma(double sigma){sigma_ = sigma;}; 
        void set_mu(double mu){mu_ = mu;};
        
    private: 
        double acklam_standard_inv_cdf(double p) const;
        double mu_; 
        double sigma_; 
}; 






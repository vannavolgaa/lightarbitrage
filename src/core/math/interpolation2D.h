#pragma once
#include <iostream>
#include <map>
#include <algorithm>
#include "../../../src/errors/core/math.h" 

class Interpolation2D
{
    public : 
        void set_values(std::map<double, double> mapped_x_y_)
        {   
            int n = mapped_x_y_.size();
            if (n < 2){throw MathError(MathErrorCode::INTERPOLATION_WRONG_VECTOR_SIZE);}
            int i = 0;
            for(auto const& imap: mapped_x_y_)
            {
                if (i == 0){x_min = imap.first;}
                if (i == n-1){x_max = imap.first;}
                i++;
                x_vector.push_back(imap.first);
                y_vector.push_back(imap.second);
                
            }    
        };
        double evaluate(double x_value) const 
        {
            if (x_value < x_min || x_value > x_max) {throw MathError(MathErrorCode::INTERPOLATION_OUT_OF_RANGE);}
            return implemented_evaluate(x_value);
        };
        double get_x_min() const {return x_min;}; 
        double get_x_max() const {return x_max;}; 
        std::vector<double> get_x_vector() const {return x_vector;};
        std::vector<double> get_y_vector() const {return y_vector;};
        Interpolation2D(std::map<double, double> mapped_x_y_){set_values(mapped_x_y_);};
        virtual ~Interpolation2D() = default;

    protected: 
        virtual double implemented_evaluate(double x_value) const = 0;

    private:
        std::vector<double> x_vector; 
        std::vector<double> y_vector; 
        double x_min ;
        double x_max ;
};

class LinearInterpolation2D: public Interpolation2D
{
    public :
        static double linear_interpolate(double x0, double y0, double x1, double y1, double x) {return y0 + (x - x0) * (y1 - y0) / (x1 - x0);};
        LinearInterpolation2D(std::map<double, double> mapped_x_y_): Interpolation2D(mapped_x_y_){};
        ~LinearInterpolation2D() = default;

    protected:
        double implemented_evaluate(double x_value) const override
        {
            std::vector<double> x_vector_ = get_x_vector();
            std::vector<double> y_vector_ = get_y_vector();
            int n = x_vector_.size();
            for (int i = 1; i < n; ++i) {
                if (x_value >= x_vector_[i-1] && x_value <= x_vector_[i]) {
                    double x0 = x_vector_[i-1]; 
                    double x1 = x_vector_[i];
                    return linear_interpolate(x0, y_vector_[i-1],x1, y_vector_[i],x_value);
                }
            }
            return 0.0;
        }
};

class CubicSpline2D: public Interpolation2D
{
    public :
        CubicSpline2D(std::map<double, double> mapped_x_y_): Interpolation2D(mapped_x_y_){get_parameters();};
        ~CubicSpline2D() = default;

    protected:
        double implemented_evaluate(double x_value) const override
        {
            std::vector<double> x_vector_ = get_x_vector();
            std::vector<double> y_vector_ = get_y_vector();
            int n = x_vector_.size();
            for (int i = 1; i < n; ++i) {
                if (x_value >= x_vector_[i-1] && x_value <= x_vector_[i]) {
                    double dx = x_value - x_vector_[i-1];
                    return a_[i-1] + b_[i-1] * dx + c_[i-1] * dx * dx + d_[i-1] * dx * dx * dx;
                }
            }
            return 0.0;
        }

    private : 
        void get_parameters()
        {
            x_ = get_x_vector();
            int n = x_.size() - 1; 
            a_ = get_y_vector();
            
            std::vector<double> h(n);
            std::vector<double> alpha(n, 0.0);

            for (int i = 0; i < n; ++i) {
                h[i] = x_[i + 1] - x_[i];
            }

            for (int i = 1; i < n; ++i) {
                alpha[i] = (3.0 / h[i]) * (a_[i + 1] - a_[i]) - (3.0 / h[i - 1]) * (a_[i] - a_[i - 1]);
            }

            std::vector<double> l(n + 1, 0.0);
            std::vector<double> mu(n, 0.0);
            std::vector<double> z(n + 1, 0.0);
            
            l[0] = 1.0;  
            mu[0] = 0.0;
            z[0] = 0.0;

            for (int i = 1; i < n; ++i) {
                l[i] = 2.0 * (x_[i + 1] - x_[i - 1]) - h[i - 1] * mu[i - 1];
                mu[i] = h[i] / l[i];
                z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i];
            }

            l[n] = 1.0;
            z[n] = 0.0;
            
            c_.resize(n + 1, 0.0);
            b_.resize(n, 0.0);
            d_.resize(n, 0.0);
            
            for (int j = n - 1; j >= 0; --j) {
                c_[j] = z[j] - mu[j] * c_[j + 1];
                b_[j] = (a_[j + 1] - a_[j]) / h[j] - h[j] * (c_[j + 1] + 2.0 * c_[j]) / 3.0;
                d_[j] = (c_[j + 1] - c_[j]) / (3.0 * h[j]);
            }
        }

        std::vector<double> a_, b_, c_, d_, x_;
};

inline std::shared_ptr<Interpolation2D> create_interpolation_2D_object(std::map<double, double> mapped_x_y_, bool linear = true)
{
    if (linear) return std::make_shared<LinearInterpolation2D>(mapped_x_y_);
    else return std::make_shared<CubicSpline2D>(mapped_x_y_);
}

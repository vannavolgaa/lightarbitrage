#include "lbr.h"

namespace LetsBeRational
{
    const double min_r = -(1 - sqrt(DBL_EPSILON));
    const double max_r = 2 / (DBL_EPSILON * DBL_EPSILON);

    double get_rational_cubic_interpolate(double x, double x0, double x1,double y0, double y1,double dy0, double dy1,double r)
    {
        double h = x1 - x0;
        double s = (x - x0) / h;

        double numerator =
            y1 * s * s * s + 
            (r * y1 - h * dy1) * s * s * (1 - s) +
            (r * y0 + h * dy0) * s * (1 - s) * (1 - s) +
            y0 * (1 - s) * (1 - s) * (1 - s);

        double denominator = 1 + (r - 3) * s * (1 - s);

        return numerator / denominator;
    };

    double get_r_left_side(double x0, double x1,double y0, double y1,double dy0, double dy1, double ddy0)
    {
        double h = x1 - x0;
        double delta = (y1 - y0) / h;
        double denominator = (delta - dy0);
        double numerator =  (.5*h*ddy0 + dy1 - dy0);
        if (abs(denominator) < DBL_MIN) return (numerator > 0) ? max_r : min_r;
        else return numerator / denominator;
    };

    double get_r_right_side(double x0, double x1,double y0, double y1,double dy0, double dy1, double ddy1)
    {
        double h = x1 - x0;
        double delta = (y1 - y0) / h;
        double denominator = (dy1 - delta);
        double numerator =  (.5*h*ddy1 + dy1 - dy0);
        if (abs(denominator) < DBL_MIN) return (numerator > 0) ? max_r : min_r;
        else return numerator / denominator;
    };

    double get_right_asymptotic(double sigma)
    {
        return stdnorm.cdf(-sigma/2);
    };

    double get_right_asymptotic_first_derivative(double x, double sigma)
    {
        return (abs(x)<DBL_MIN) ? -.5 : -.5 * exp(.5*x*x/(sigma*sigma));
    };

    double get_right_asymptotic_second_derivative(double x, double sigma)
    {
        double sigsq = sigma*sigma;
        return (abs(x)<DBL_MIN) ? 0 : exp(sigsq / 8 +  x*x/sigsq)*sqrt(PI/2) * x * x / (sigsq * sigma);
    };

    double get_left_asymptotic_z(double x, double sigma)
    {
        return -abs(x)/(sqrt(3)*sigma);
    }

    double get_left_asymptotic(double x, double sigma)
    {
        double cdfz = stdnorm.cdf(get_left_asymptotic_z(x, sigma));
        return (2*PI*cdfz*cdfz*cdfz*abs(x))/(3*sqrt(3));
    };

    double get_left_asymptotic_first_derivative(double x, double sigma)
    {
        double z = get_left_asymptotic_z(x, sigma);
        double cdfz = stdnorm.cdf(z);
        return 2*PI*z*z*cdfz*cdfz*exp(z*z+sigma*sigma/8.0);
    };

    double get_left_asymptotic_second_derivative(double x, double sigma)
    {
        double z = get_left_asymptotic_z(x, sigma);
        double cdfz = stdnorm.cdf(z);
        double pdfz = stdnorm.pdf(z);
        double sigmasq = sigma*sigma;
        double term1 = (PI/6)*(z*z/(sigmasq*sigma))*cdfz*exp(2*z*z+sigmasq/4.0);
        double term2 = 8*sqrt(3)*sigma*abs(x)+(3*sigmasq*(sigmasq-8)-8*x*x)*cdfz/pdfz;
        return term1*term2;
    }

    std::tuple<double, double, double, double, double, double, double, double, double, double, bool> get_initial_data(double beta, double x, double is_call)
    {
        int c = is_call ? 1 : -1;
        if (c*x>0){
            beta = fabs(std::max(beta-NormalizedBlack::get_intrisic(x, is_call),0.));
            is_call = not is_call;
        }
        if (!is_call) x = -x; is_call = true;
        double b_max = exp(.5*x);
        double sigma_c = sqrt(2*abs(x));
        double b_c = NormalizedBlack::get_price(x, sigma_c, is_call);
        double v_c = NormalizedBlack::get_vega(x, sigma_c);
        double sigma_u = (v_c>DBL_MIN) ? sigma_c + (b_max - b_c)/v_c : sigma_c;
        double sigma_l = (v_c>DBL_MIN) ? sigma_c - b_c/v_c : sigma_c;
        double b_u = NormalizedBlack::get_price(x, sigma_u, is_call);
        double b_l = NormalizedBlack::get_price(x, sigma_l, is_call);
        return std::make_tuple(beta,x,b_max,sigma_l, sigma_c, sigma_u, b_l, b_c, b_u, v_c, is_call); 
    }

    double get_initial_guess_normalized_volatility(double beta, double x, bool is_call)
    {
        auto [beta_,x_,b_max,sigma_l, sigma_c, sigma_u, b_l, b_c, b_u, v_c, is_call_] = get_initial_data(beta, x, is_call);
        if (beta_ < b_l && beta_ >= 0) {
            double fl = get_left_asymptotic(x_, sigma_l);
            double dfl = get_left_asymptotic_first_derivative(x_, sigma_l);
            double ddfl = get_left_asymptotic_second_derivative(x_, sigma_l);
            double r = get_r_right_side(0, b_l, 0, fl, 1, dfl, ddfl);
            double frc = get_rational_cubic_interpolate(beta_, 0, b_l, 0, fl, 1, dfl, r);
            double sq3 = sqrt(3), a = frc/(2*PI*abs(x_));
            return abs(x_/(sq3*stdnorm.inv_cdf(sq3*std::pow(a, 1.0/3.0))));
        }
        else if (beta_ <= b_c && beta_ >= b_l){
            double v_l = NormalizedBlack::get_vega(x_, sigma_l);
            double r = get_r_right_side(b_l, b_c, sigma_l, sigma_c, 1/v_l, 1/v_c, 0);
            return get_rational_cubic_interpolate(beta_, b_l, b_c, sigma_l, sigma_c, 1/v_l, 1/v_c, r);
        }
        else if (beta_ <= b_u && beta_ > b_c){
            double v_u = NormalizedBlack::get_vega(x_, sigma_u);
            double r = get_r_left_side(b_c, b_u, sigma_c, sigma_u, 1/v_c, 1/v_u, 0);
            return get_rational_cubic_interpolate(beta_, b_c, b_u, sigma_c, sigma_u, 1/v_c, 1/v_u, r);
        }
        else {
            double fu = get_right_asymptotic(sigma_u);
            double dfu = get_right_asymptotic_first_derivative(x_, sigma_u);
            double ddfu = get_right_asymptotic_second_derivative(x_, sigma_u);
            double r = get_r_left_side(b_u, b_max, fu, 0, dfu, -.5, ddfu);
            double frc = get_rational_cubic_interpolate(beta_, b_u, b_max, fu, 0, dfu, -.5, r);
            return -2.0*stdnorm.inv_cdf(frc);
        }
    }

    std::function<double(double)> get_target_function(double beta, double x, double b_l, double b_u)
    {

        double b_max = exp(.5*x);
        double bt_u = std::max(b_u, b_max/2);
        if (beta >= 0 && beta < b_l) return [beta, x](double s) {return log(1/NormalizedBlack::get_call_price(x, s)) - log(1/beta);};
        else if (beta >= b_l && beta <= bt_u) return [beta, x](double s) {return NormalizedBlack::get_call_price(x, s) - beta;};
        else return [b_max, beta, x](double s) {return log((b_max-beta)/(b_max-NormalizedBlack::get_call_price(x, s)));};

    }

    std::function<double(double)> get_target_function_first_derivative(double beta, double x, double b_l, double b_u)
    {
        double b_max = exp(.5*x);
        double bt_u = std::max(b_u, b_max/2);
        if (beta >= 0 && beta < b_l) return [x](double s) {return -NormalizedBlack::get_vega(x, s)/NormalizedBlack::get_call_price(x, s);};
        else if (beta >= b_l && beta <= bt_u) return [x](double s) {return NormalizedBlack::get_vega(x, s);};
        else return [b_max, x](double s) {return NormalizedBlack::get_vega(x, s)/(b_max - NormalizedBlack::get_call_price(x, s));};

    }
    
    double get_newton_normalized_volatility(double beta, double x, bool is_call)
    {
        auto [beta_,x_,b_max,sigma_l, sigma_c, sigma_u, b_l, b_c, b_u, v_c, is_call_] = get_initial_data(beta, x, is_call);
        if (beta_<=0) return 0.0;
        NewtonRaphson newton(
            get_initial_guess_normalized_volatility(beta, x, is_call),
            get_target_function(beta_, x_, b_l, b_u),
            get_target_function_first_derivative(beta_, x_, b_l, b_u)
        );
        newton.set_tolerance_rate(1e-20);
        newton.set_max_iterations(100);
        std::shared_ptr<OptimizerResult> result = newton.optimize();
        return result->x_[0];
    }

    double get_newton_black_scholes_implied_volatility(double price, double S, double K, double T, double r, double q, bool is_call, bool is_future)
    {
        if (price<0) return NAN;
        double F = S*exp(get_future_flag(is_future)*(r-q)*T);
        double x = log(F/K);
        double beta = exp(r*T)*price/sqrt(F*K);
        if (abs(NormalizedBlack::get_vega(x, get_initial_guess_normalized_volatility(beta, x, is_call)))<1e-5) return NAN;
        double iv = get_newton_normalized_volatility(beta, x, is_call)/sqrt(T);
        if (iv<=0.0) return NAN;
        else return iv; 
    } 
}

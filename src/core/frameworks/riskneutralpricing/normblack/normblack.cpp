#include "normblack.h"

inline double get_lewis_normalized_black_price_gauss_laguerre(
    double x, bool is_call,
    std::function<std::complex<double>(std::complex<double>)> cf, 
    std::shared_ptr<GaussianQuadrature> quadrature)
{
    int n = quadrature->points;
    const std::vector<double>& laguerre_nodes = quadrature->roots;
    const std::vector<double>& laguerre_weights = quadrature->weights;
    const std::complex<double> i(0.0, 1.0);
    double sum = 0.0;
    for (int j = 0; j < n; ++j) {
        double u = laguerre_nodes[j];
        std::complex<double> integrand = std::exp(i * u * x) * cf(std::complex<double>(u, -0.5));
        double rintegrand = abs(integrand.real()) < DBL_MIN ? 0.0 : integrand.real();
        if (rintegrand == 0.0 or abs(laguerre_weights[j])<1e-20) sum += 0.0;
        else sum += laguerre_weights[j] * std::exp(u) * rintegrand / (u * u + 0.25);
    }
    double call = std::exp(x / 2.0) - sum / M_PI;
    return is_call ? call : call - NormalizedBlack::get_intrisic(x, true);
}

inline double get_lewis_normalized_black_price_gauss_laguerre(
    double x, bool is_call,
    std::function<std::complex<double>(std::complex<double>)> cf, 
    int points)
{
    std::shared_ptr<GaussianQuadrature> quadrature = get_gauss_laguerre_quadrature(points);
    return get_lewis_normalized_black_price_gauss_laguerre(x, is_call, cf, quadrature);
}

namespace NormalizedBlack
{
    // Code source from https://github.com/vollib/lets_be_rational/
    double get_call_price_region_1(double h, double t)
    {
        const double e=(t/h)*(t/h), r=((h+t)*(h-t)), q=(h/r)*(h/r);
        const double asymptotic_expansion_sum = (2.0+q*(-6.0E0-2.0*e+3.0*q*(1.0E1+e*(2.0E1+2.0*e)+5.0*q*(-1.4E1+e*(-7.0E1+e*(-4.2E1-2.0*e))+7.0*q*(1.8E1+e*(1.68E2+e*(2.52E2+e*(7.2E1+2.0*e)))+9.0*q*(-2.2E1+e*(-3.3E2+e*(-9.24E2+e*(-6.6E2+e*(-1.1E2-2.0*e))))+1.1E1*q*(2.6E1+e*(5.72E2+e*(2.574E3+e*(3.432E3+e*(1.43E3+e*(1.56E2+2.0*e)))))+1.3E1*q*(-3.0E1+e*(-9.1E2+e*(-6.006E3+e*(-1.287E4+e*(-1.001E4+e*(-2.73E3+e*(-2.1E2-2.0*e))))))+1.5E1*q*(3.4E1+e*(1.36E3+e*(1.2376E4+e*(3.8896E4+e*(4.862E4+e*(2.4752E4+e*(4.76E3+e*(2.72E2+2.0*e)))))))+1.7E1*q*(-3.8E1+e*(-1.938E3+e*(-2.3256E4+e*(-1.00776E5+e*(-1.84756E5+e*(-1.51164E5+e*(-5.4264E4+e*(-7.752E3+e*(-3.42E2-2.0*e))))))))+1.9E1*q*(4.2E1+e*(2.66E3+e*(4.0698E4+e*(2.3256E5+e*(5.8786E5+e*(7.05432E5+e*(4.0698E5+e*(1.08528E5+e*(1.197E4+e*(4.2E2+2.0*e)))))))))+2.1E1*q*(-4.6E1+e*(-3.542E3+e*(-6.7298E4+e*(-4.90314E5+e*(-1.63438E6+e*(-2.704156E6+e*(-2.288132E6+e*(-9.80628E5+e*(-2.01894E5+e*(-1.771E4+e*(-5.06E2-2.0*e))))))))))+2.3E1*q*(5.0E1+e*(4.6E3+e*(1.0626E5+e*(9.614E5+e*(4.08595E6+e*(8.9148E6+e*(1.04006E7+e*(6.53752E6+e*(2.16315E6+e*(3.542E5+e*(2.53E4+e*(6.0E2+2.0*e)))))))))))+2.5E1*q*(-5.4E1+e*(-5.85E3+e*(-1.6146E5+e*(-1.77606E6+e*(-9.37365E6+e*(-2.607579E7+e*(-4.01166E7+e*(-3.476772E7+e*(-1.687257E7+e*(-4.44015E6+e*(-5.9202E5+e*(-3.51E4+e*(-7.02E2-2.0*e))))))))))))+2.7E1*q*(5.8E1+e*(7.308E3+e*(2.3751E5+e*(3.12156E6+e*(2.003001E7+e*(6.919458E7+e*(1.3572783E8+e*(1.5511752E8+e*(1.0379187E8+e*(4.006002E7+e*(8.58429E6+e*(9.5004E5+e*(4.7502E4+e*(8.12E2+2.0*e)))))))))))))+2.9E1*q*(-6.2E1+e*(-8.99E3+e*(-3.39822E5+e*(-5.25915E6+e*(-4.032015E7+e*(-1.6934463E8+e*(-4.1250615E8+e*(-6.0108039E8+e*(-5.3036505E8+e*(-2.8224105E8+e*(-8.870433E7+e*(-1.577745E7+e*(-1.472562E6+e*(-6.293E4+e*(-9.3E2-2.0*e))))))))))))))+3.1E1*q*(6.6E1+e*(1.0912E4+e*(4.74672E5+e*(8.544096E6+e*(7.71342E7+e*(3.8707344E8+e*(1.14633288E9+e*(2.07431664E9+e*(2.33360622E9+e*(1.6376184E9+e*(7.0963464E8+e*(1.8512208E8+e*(2.7768312E7+e*(2.215136E6+e*(8.184E4+e*(1.056E3+2.0*e)))))))))))))))+3.3E1*(-7.0E1+e*(-1.309E4+e*(-6.49264E5+e*(-1.344904E7+e*(-1.4121492E8+e*(-8.344518E8+e*(-2.9526756E9+e*(-6.49588632E9+e*(-9.0751353E9+e*(-8.1198579E9+e*(-4.6399188E9+e*(-1.6689036E9+e*(-3.67158792E8+e*(-4.707164E7+e*(-3.24632E6+e*(-1.0472E5+e*(-1.19E3-2.0*e)))))))))))))))))*q)))))))))))))))));
        return ONE_OVER_SQRT_TWO_PI*exp((-0.5*(h*h+t*t)))*(t/r)*asymptotic_expansion_sum;
    };

    // Code source from https://github.com/vollib/lets_be_rational/
    double get_call_price_region_2(double h, double t)
    {
        const double a = 1+h*(0.5*SQRT_TWO_PI)*erfcx_cody(-ONE_OVER_SQRT_TWO*h), w=t*t, h2=h*h;
        const double expansion = 2*t*(a+w*((-1+3*a+a*h2)/6+w*((-7+15*a+h2*(-1+10*a+a*h2))/120+w*((-57+105*a+h2*(-18+105*a+h2*(-1+21*a+a*h2)))/5040+w*((-561+945*a+h2*(-285+1260*a+h2*(-33+378*a+h2*(-1+36*a+a*h2))))/362880+w*((-6555+10395*a+h2*(-4680+17325*a+h2*(-840+6930*a+h2*(-52+990*a+h2*(-1+55*a+a*h2)))))/39916800+((-89055+135135*a+h2*(-82845+270270*a+h2*(-20370+135135*a+h2*(-1926+25740*a+h2*(-75+2145*a+h2*(-1+78*a+a*h2))))))*w)/6227020800.0))))));
        return ONE_OVER_SQRT_TWO_PI*exp((-0.5*(h*h+t*t)))*expansion;
    };

    double get_call_price_region_3(double h, double t)
    {
        const NormalDistribution norm(0.0,1.0);
        return norm.cdf(h+t)*exp(h*t) - norm.cdf(h-t)/exp(h*t);
    };

    double get_call_price_region_4(double h, double t)
    {
        return 0.5 * exp(-0.5*(h*h+t*t)) * (erfcx_cody(-ONE_OVER_SQRT_TWO*(h+t)) - erfcx_cody(-ONE_OVER_SQRT_TWO*(h-t)));
    };

    double get_intrisic(double x, bool is_call) 
    {
        if (get_call_put_flag(is_call)*x <= 0) return 0;
        double bm = exp(.5*x);
        return std::max(get_call_put_flag(is_call)*(bm - 1/bm), 0.0);
    };

    double get_call_price(double x, double normalized_sigma)
    {
        if (x>0) return get_intrisic(x, true)+get_call_price(-x,normalized_sigma);
        if (normalized_sigma<=0) return get_intrisic(x, true);
        if ( x < normalized_sigma*H_LARGE  &&  0.5*normalized_sigma*normalized_sigma+x < normalized_sigma*(T_SMALL+H_LARGE)) return get_call_price_region_1(x/normalized_sigma, .5*normalized_sigma);
        if ( 0.5*normalized_sigma < T_SMALL) return get_call_price_region_2(x/normalized_sigma, .5*normalized_sigma);
        if (x+0.5*normalized_sigma*normalized_sigma > normalized_sigma*0.85) return get_call_price_region_3(x/normalized_sigma, .5*normalized_sigma);
        return get_call_price_region_4(x/normalized_sigma, .5*normalized_sigma);
    };

    double get_vega(double x, double normalized_sigma)
    {
        if (normalized_sigma<=0) return 0;
        return ONE_OVER_SQRT_TWO_PI*exp(-0.5*(x*x/(normalized_sigma*normalized_sigma)+0.25*normalized_sigma*normalized_sigma));
    };

    double get_volga(double x, double normalized_sigma)
    {
        if (normalized_sigma<=0) return 0;
        return (x/(normalized_sigma*normalized_sigma*normalized_sigma) - .125)*get_vega(x, normalized_sigma);
    };

    double get_price(double x, double normalized_sigma, bool is_call)
    {
        if (normalized_sigma<=0) return get_intrisic(x, is_call);
        if (!is_call) x = -x;
        return get_call_price(x,normalized_sigma);
    }

    std::complex<double> get_characteristic_function(std::complex<double> u, double x, double normalized_sigma)
    {
        return std::exp(-0.5*u*(u+std::complex<double>(0.0,1.0))*normalized_sigma*normalized_sigma);
    }

    double get_lewis_price(double x, double normalized_sigma, bool is_call, std::shared_ptr<GaussianQuadrature> quadrature)
    {
        std::function<std::complex<double>(std::complex<double>)> cf = [normalized_sigma, x](std::complex<double> u) {
            return get_characteristic_function(u, x, normalized_sigma);
        };
        return get_lewis_normalized_black_price_gauss_laguerre(x, is_call, cf, quadrature);
    }
}

#include <iostream>
#include <iomanip>
#include <memory>
#include <cassert>
#include <map>
#include <chrono>
#include "../../../../../src/core/frameworks/riskneutralpricing/letsberational/lbr.h" 
#include "../../../../../src/core/frameworks/riskneutralpricing/blackscholes/bs.h" 

void test_get_newton_normalized_volatility() 
{
    std::cout << "Testing get_newton_normalized_volatility..." << std::endl;
    for (double s = 0.01; s < 7; s+= .5)
    {
        double x = -4;
        for (double x = -10; x<0; x+=0.5)
        {

            double b = NormalizedBlack::get_price(x,s,true); 
            double s_lbr = LetsBeRational::get_newton_normalized_volatility(b,x,true);
            if (s_lbr>0) assert(abs(s_lbr-s)<= 1e-6);

        }
        
    }
        
    std::cout << "All tests passed for get_newton_normalized_volatility!" << std::endl;
}

void test_get_newton_black_scholes_implied_volatility()
{
    std::cout << "Testing get_newton_black_scholes_implied_volatility..." << std::endl;
    double S = 100; 
    double r = 0.01; 
    double q = 0.02; 
    double bs_price, bs_iv;
    for (double K = 1; K <= 500; K+= 10)
    {
        for (double T = 0.001; T <=30; T+= 0.5)
        {
            for (double s= 0.01; s<=5; s+=0.1)
            {
                BlackScholes bs(S,K,T,r,q,s,true,false);
                bs_price = bs.get_price();
                bs_iv = LetsBeRational::get_newton_black_scholes_implied_volatility(bs_price,S,K,T,r,q,true,false);
                if (!isnan(bs_iv)) assert(abs(bs_iv-s)<= 1e-3);
                
                bs.set_is_call(false);
                bs_price = bs.get_price();
                bs_iv = LetsBeRational::get_newton_black_scholes_implied_volatility(bs_price,S,K,T,r,q,false,false);
                if (!isnan(bs_iv)) assert(abs(bs_iv-s)<= 1e-3);

                bs.set_is_future(true);
                bs_price = bs.get_price();
                bs_iv = LetsBeRational::get_newton_black_scholes_implied_volatility(bs_price,S,K,T,r,q,false,true);
                if (!isnan(bs_iv)) assert(abs(bs_iv-s)<= 1e-3);

                bs.set_is_call(true);
                bs_price = bs.get_price();
                bs_iv = LetsBeRational::get_newton_black_scholes_implied_volatility(bs_price,S,K,T,r,q,true,true);
                if (!isnan(bs_iv)) assert(abs(bs_iv-s)<= 1e-3);

            }
        }
    }
    std::cout << "All tests passed for get_newton_black_scholes_implied_volatility!" << std::endl;
}

int main()
{
    test_get_newton_normalized_volatility();
    test_get_newton_black_scholes_implied_volatility();
    return 0; 
}


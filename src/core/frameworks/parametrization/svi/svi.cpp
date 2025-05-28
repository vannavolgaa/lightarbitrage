#include "svi.h"

bool SSVI::check_butterfly_arbitrage(double atm_total_variance) const
{
    double f = get_f(atm_total_variance);
    double cond1 = atm_total_variance*f*(1+abs(rho_));
    double cond2 = cond1*f;
    if (cond1>0 and cond2>0) return false;
    else return true;
};

bool SSVI::check_calendar_spread_arbitrage(double atm_total_variance) const
{
    double f = get_f(atm_total_variance);
    double df = get_df(atm_total_variance);
    double value = f*(1+sqrt(1-rho_*rho_))/(rho_*rho_);
    if (df>=0 and df<=value) return false;
    else return true;
};

double SSVI::get_total_variance(double k, double atm_total_variance) const
{
    double f = get_f(atm_total_variance);
    double term1 = f*k+rho_;
    double term2 = sqrt(term1*term1 + (1-rho_*rho_));
    return .5*atm_total_variance*(1+rho_*k*f+term2);
};

std::shared_ptr<SVI> SSVI::get_svi(double atm_total_variance, double t) const
{
    double f = get_f(atm_total_variance);
    return std::make_shared<SVI>(
        atm_total_variance/t, 
        .5*rho_*f*sqrt(atm_total_variance),
        .5*(1+rho_)*f*sqrt(atm_total_variance), 
        .5*(1-rho_)*f*sqrt(atm_total_variance), 
        atm_total_variance*(1-rho_*rho_)/t, 
        t);
};

double SSVI::get_risk_neutral_density(double k, double atm_total_variance, double t) const
{
    std::shared_ptr<SVI> svi = get_svi(atm_total_variance, t);
    return svi->get_risk_neutral_density(k); 
};

double SSVI::get_local_volatility(double k, double atm_total_variance, double t) const
{
    std::shared_ptr<SVI> svi = get_svi(atm_total_variance, t);
    return svi->get_local_volatility(k); 
};

double SSVI::get_atm_volatility_skew(double k, double atm_total_variance, double t) const
{
    std::shared_ptr<SVI> svi = get_svi(atm_total_variance, t);
    return svi->get_ut(); 
};

double SSVI::get_undiscounted_normalized_black_price(double k, double atm_total_variance, double t, bool is_call) const
{
    return NormalizedBlack::get_price(k,get_implied_volatility(k, atm_total_variance,t)*sqrt(t),is_call);
};

void SVI::update_parameters()
{
    set_b();
    set_p();
    set_beta();
    set_alpha();
    set_m();
    set_a();
    set_s();
    set_dbdt();
    set_dmdt();
    set_dsdt(); 
    set_dadt();
    if(a_+b_*s_*sqrt(1-p_*p_)<0) throw FrameworksError(FrameworksErrorCode::RAW_SVI_CONDITION);
};

void SVI::set_b()
{
    b_ = sqrt(vt_*T_)*(ct_+pt_)/2;
    if (b_<0) throw FrameworksError(FrameworksErrorCode::RAW_SVI_INVALID_PARAM_B);
};

void SVI::set_p()
{
    if (b_==0) p_ = 0.0;
    else p_ = 1 - pt_*sqrt(vt_*T_)/b_;
    if (abs(p_)>1.0) throw FrameworksError(FrameworksErrorCode::RAW_SVI_INVALID_PARAM_P);
}

void SVI::set_beta()
{
    if (b_==0) beta_ = .5;
    else beta_ = p_ - 2*ut_*sqrt(vt_*T_)/b_;
    if (abs(beta_)>1.0) throw FrameworksError(FrameworksErrorCode::RAW_SVI_INVALID_PARAM_BETA);
}

void SVI::set_alpha()
{
    if (beta_<0) alpha_ = -sqrt(1/(beta_*beta_) - 1);
    if (beta_>0) alpha_ = sqrt(1/(beta_*beta_) - 1);
    if (beta_==0) alpha_ = 0.0;
}

void SVI::set_m()
{
    if (b_==0) m_ = 0.0;
    else{
        double intterm;
        if (alpha_<0) intterm = -sqrt(1+alpha_*alpha_);
        else intterm = sqrt(1+alpha_*alpha_);
        m_ = (T_*(vt_-vmt_)/(b_*(-p_+intterm-alpha_*sqrt(1-p_*p_))));
    }
};

void SVI::set_a()
{
    if(m_==0) a_ = T_*(vmt_+vt_*sqrt(1-p_*p_))/(1-sqrt(1-p_*p_));
    else a_ = T_*vmt_-b_*alpha_*m_*sqrt(1-p_*p_);
};

void SVI::set_s()
{
    if(m_==0){
        if (b_==0) s_ = 1.0;
        else s_ = (vt_*T_ - a_)/b_;
    }
    else s_ = alpha_*m_;
    if (s_<=0) throw FrameworksError(FrameworksErrorCode::RAW_SVI_INVALID_PARAM_S);
};

void SVI::set_dbdt(){dbdt_ = sqrt(vt_)*(ct_+pt_)/(4*sqrt(T_));};

void SVI::set_dmdt()
{
    if (b_==0) dmdt_ = 0.0;
    else dmdt_ = m_ * (b_ - T_*dbdt_)/(T_*b_);
};

void SVI::set_dsdt()
{
    if (m_==0){
        if (b_==0) dsdt_ = 0; 
        else dsdt_ = (b_*(vt_-a_/T_)+dbdt_*(vt_*T_-a_))/(b_*b_);
    }
    else dsdt_ = alpha_*dmdt_;
}

void SVI::set_dadt()
{
    if (m_==0) dadt_ = a_/T_;
    else dadt_ = vmt_-sqrt(1-p_*p_)*(dbdt_*s_+dsdt_*b_);
};

bool SVI::check_butterfly_arbitrage() const
{
    double cond1 = sqrt(vt_*T_)*std::max(ct_,pt_);
    double cond2 = (ct_+pt_)*std::max(ct_,pt_);
    if (cond1<2 and cond2 <=2){return false;}
    else{return true;}
};

std::shared_ptr<SSVI> SVI::get_ssvi() const
{
    double rho_ = 1/(1+pt_/ut_);
    return std::make_shared<SSVI>(rho_, 2*ut_/rho_, .5);
};

double SVI::get_risk_neutral_density(double k) const
{
    double w = get_total_variance(k);
    double dwdk = get_dwdk(k);
    double dw2dk2 = get_dw2dk2(k);
    double sqdwdk = dwdk*dwdk;
    double term = (1 - .5*k*dwdk/w);
    return term * term - .25 * sqdwdk * (0.25 + 1/w) + .5 * dw2dk2;
};

double SVI::get_asymptotic_value(bool left, double epsilon) const
{
    int x;
    if (left) x = -1;
    else x = 1;
    double term = (1 - epsilon*(1 + x*p_));
    return m_ + x * s_ * term / sqrt(1 - term*term);
}

bool SVI::check_numerical_grid_cs(std::shared_ptr<SVI> slice, int grid_size, double epsilon) const
{
    if (grid_size<=1) grid_size = 100;
    epsilon = abs(epsilon); 
    int x; 
    double st = slice->get_T();

    if (st > T_) x = 1;
    else if (st < T_) x = -1;
    else return false; 

    double vleft = get_asymptotic_value(true, epsilon);
    double vright = get_asymptotic_value(false, epsilon);
    double increment = (vright - vleft) / (grid_size-1);
    for (int i = 0; i < grid_size; i++)
    {
        double k = vleft + i * increment;
        double sw = slice->get_total_variance(k);
        double w = get_total_variance(k);
        if (x*sw<=x*w) return true;
    }
    return false;
};

bool SVI::check_calendar_spread_arbitrage(std::shared_ptr<SVI> slice, int grid_size, double epsilon) const
{
    if (grid_size<=0) grid_size = 100; 
    epsilon = abs(epsilon); 
    int x; 
    double st = slice->get_T();
    if (st > T_) x = 1;
    else if (st < T_) x = -1;
    else return false;

    double sw = slice->get_atm_total_variance(); 
    double w = get_atm_total_variance();
    double spt = slice->get_pt();
    double sct = slice->get_ct();
    double svmt = slice->get_vmt();
    if (x*sw<=x*w) return true;
    if (x*st*svmt<=x*T_*vmt_) return true;
    if (x*sqrt(sw)*spt<=x*sqrt(w)*pt_) return true;
    if (x*sqrt(sw)*sct<=x*sqrt(w)*ct_) return true;

    return check_numerical_grid_cs(slice, grid_size, epsilon);
};

double SVI::get_undiscounted_normalized_black_price(double k, double t, bool is_call) const
{
    return NormalizedBlack::get_price(k,get_implied_volatility(k)*sqrt(t),is_call);
}
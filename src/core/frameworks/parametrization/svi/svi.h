#pragma once 
#include <iostream>
#include "../../../../../src/core/frameworks/riskneutralpricing/normblack/normblack.h"
#include "../../../../../src/errors/core/frameworks.h"

class SVI; 

class ReducedSVI; 

class SSVI 
{
    public: 
        void set_rho(double rho) {rho_ = rho;};
        void set_nu(double nu) {nu_ = nu;};
        void set_gamma(double gamma) {gamma_ = gamma;};
        double get_rho() const {return rho_;};
        double get_nu() const {return nu_;};
        double get_gamma() const {return gamma_;};
        std::shared_ptr<SVI> get_svi(double atm_total_variance, double t) const;
        bool check_butterfly_arbitrage(double atm_total_variance) const; 
        bool check_calendar_spread_arbitrage(double atm_total_variance) const; 
        double get_total_variance(double k, double atm_total_variance) const;
        double get_implied_variance(double k, double atm_total_variance, double t) const {return get_total_variance(k,atm_total_variance)/t;};
        double get_implied_volatility(double k, double atm_total_variance, double t) const {return sqrt(get_implied_variance(k,atm_total_variance,t));};
        double get_risk_neutral_density(double k, double atm_total_variance, double t) const;
        double get_local_volatility(double k, double atm_total_variance, double t) const; 
        double get_atm_volatility_skew(double k, double atm_total_variance, double t) const;
        double get_undiscounted_normalized_black_price(double k, double atm_total_variance, double t, bool is_call) const;
        SSVI(const double rho, const double nu, const double gamma): 
            rho_(rho), nu_(nu), gamma_(gamma)
            {
                if (std::abs(rho)>1 or nu<0 or gamma<0 or gamma>1) throw FrameworksError(FrameworksErrorCode::SSVI_INVALID_PARAM);
            };
        ~SSVI(){}; 
    private: 
        double get_f(double atm_total_variance) const {return nu_*pow(atm_total_variance, -gamma_);}
        double get_df(double atm_total_variance) const {return (1-gamma_)*get_f(atm_total_variance);}
        double rho_;
        double nu_;
        double gamma_;
};

class SVI
{
    public: 
        
        double get_vt() const {return vt_;};
        double get_ut() const {return ut_;};
        double get_ct() const {return ct_;};
        double get_pt() const {return pt_;};
        double get_vmt() const {return vmt_;};
        double get_T() const {return T_;};
        void set_vt(double vt) {vt_ = vt; update_parameters();};
        void set_ut(double ut) {ut_ = ut; update_parameters();};
        void set_ct(double ct) {ct_ = ct; update_parameters();};
        void set_pt(double pt) {pt_ = pt; update_parameters();};
        void set_vmt(double vmt) {vmt_ = vmt; update_parameters();};
        void set_T(double T) {T_ = T; update_parameters();};
        std::shared_ptr<SSVI> get_ssvi() const;
        bool check_butterfly_arbitrage() const; 
        bool check_calendar_spread_arbitrage(std::shared_ptr<SVI> slice, int grid_size, double epsilon) const; 
        double get_total_variance(double k) const {return a_ + b_*(p_*(k-m_)+sqrt((k-m_)*(k-m_)+s_*s_));};
        double get_atm_total_variance() const {return vt_*T_;};
        double get_minimum_total_variance() const {return vmt_*T_;};
        double get_implied_variance(double k) const {return get_total_variance(k)/T_;};
        double get_implied_volatility(double k) const {return sqrt(get_implied_variance(k));};
        double get_risk_neutral_density(double k) const; 
        double get_local_variance(double k) const {return get_dwdt(k)/get_risk_neutral_density(k);};
        double get_local_volatility(double k) const {return sqrt(get_local_variance(k));};
        double get_asymptotic_value(bool left, double epsilon) const;
        double get_undiscounted_normalized_black_price(double k, double t, bool is_call) const;
        SVI(double vt, double ut, double ct, double pt, double vmt, double t): 
            vt_(vt), ut_(ut), ct_(ct), pt_(pt), vmt_(vmt), T_(t)
            {
                if (T_<=0 or vt_<=0 or vmt_<=0) throw FrameworksError(FrameworksErrorCode::JUMP_WINGS_SVI_INVALID_PARAM);
                update_parameters();
            };
        ~SVI(){};

    private:
        void update_parameters();
        void set_b(); 
        void set_p(); 
        void set_beta();
        void set_alpha();
        void set_m();
        void set_a();
        void set_s();
        void set_dbdt();
        void set_dmdt();
        void set_dsdt();
        void set_dadt();
        double vt_; 
        double ut_; 
        double ct_;
        double pt_;
        double vmt_; 
        double T_; 
        double b_; 
        double p_; 
        double beta_; 
        double alpha_;
        double m_; 
        double a_; 
        double s_; 
        double dbdt_; 
        double dmdt_; 
        double dsdt_;
        double dadt_;   
        double get_g(double k) const {return p_*(k-m_)+sqrt((k-m_)*(k-m_)+s_*s_);};
        double get_dgdt(double k) const {return -p_*dmdt_ + (dsdt_*s_-dmdt_*(k-m_))/sqrt((k-m_)*(k-m_)+s_*s_);};
        double get_dwdt(double k) const {return dadt_ + b_*get_dgdt(k) + dbdt_*get_g(k);};
        double get_dwdk(double k) const {return b_*(p_+(k-m_)/(sqrt((k-m_)*(k-m_) + s_*s_)));};
        double get_dw2dk2(double k) const {return b_*s_*s_/std::pow((k-m_)*(k-m_) + s_*s_,1.5);};
        bool check_numerical_grid_cs(std::shared_ptr<SVI> slice, int grid_size, double epsilon) const;
};

class ReducedSVI: public SVI
{
    public: 
        double get_eta() const {return eta_;};
        double get_rho() const {return rho_;};
        void set_eta(double eta) {eta_ = eta; set_ut(.5*eta*rho_); set_ct(.5*eta*(1+rho_)); set_pt(.5*eta*(1-rho_));};
        void set_rho(double rho) {rho_ = rho; set_ut(.5*eta_*rho); set_ct(.5*eta_*(1+rho)); set_pt(.5*eta_*(1-rho)); set_vmt(get_vt()*(1-rho*rho));};
        ReducedSVI(double vt, double eta, double rho, double t): 
        SVI(vt, .5*eta*rho, .5*eta*(1+rho), .5*eta*(1-rho), vt*(1-rho*rho), t),eta_(eta), rho_(rho){};
        ~ReducedSVI(){};
    private: 
        double eta_;
        double rho_; 
};
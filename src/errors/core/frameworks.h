#pragma once 
#include <iostream>
#include "../../../src/errors/base.h"

enum class FrameworksErrorCode
{
    INVALID_NORMALIZED_BLACK_PRICE = 0,
    JUMP_WINGS_SVI_INVALID_PARAM = 1, 
    RAW_SVI_INVALID_PARAM_B = 2, 
    RAW_SVI_INVALID_PARAM_P = 3, 
    RAW_SVI_INVALID_PARAM_BETA = 4, 
    RAW_SVI_INVALID_PARAM_S = 5, 
    RAW_SVI_CONDITION = 6, 
    SSVI_INVALID_PARAM = 7, 
    VASICEK_INVALID_PARAM = 8, 
    BLACK_SCHOLES_MODEL_INVALID_PARAM = 9,
    BARONE_ADESI_WHALEY_MODEL_INVALID_PARAM = 10,
    LETS_BE_RATIONAL_INVALID_BS_PRICE = 11,
    HESTON_INVALID_PARAM = 12, 
    HESTON_UNKNOWN_DISCRETIZATION_METHOD = 13
};

class FrameworksError: public CoreArbitrageError<FrameworksErrorCode>
{
    public:
        explicit FrameworksError(const FrameworksErrorCode code): CoreArbitrageError<FrameworksErrorCode>(code){};
        ~FrameworksError() override = default;
    protected: 
        std::string get_module_name() const override final {return "frameworks";};
        std::string get_error_message() const override final
        {
            switch (get_code())
            {
                case FrameworksErrorCode::INVALID_NORMALIZED_BLACK_PRICE:
                    return "Normalized black price is invalid. It can't be negative and ccannot exceed its upper bound b_max = e^(x/2).";
                case FrameworksErrorCode::JUMP_WINGS_SVI_INVALID_PARAM:
                    return "Invalid JW-SVI parameters. Either ATM variance, minimum variance or time to maturity is negative of equal to 0.";
                case FrameworksErrorCode::RAW_SVI_INVALID_PARAM_B:
                    return "Invalid raw SVI b parameter. Must have b>=0.";
                case FrameworksErrorCode::RAW_SVI_INVALID_PARAM_P: 
                    return "Invalid raw SVI p parameter. Must have abs(p)<=1.";
                case FrameworksErrorCode::RAW_SVI_INVALID_PARAM_BETA: 
                    return "Invalid raw SVI beta parameter. Must have abs(beta)<=1.";
                case FrameworksErrorCode::RAW_SVI_INVALID_PARAM_S: 
                    return "Invalide raw SVI s parameter. Must have s>0.";
                case FrameworksErrorCode::RAW_SVI_CONDITION: 
                    return "Invalid SVI condition, must have a+b*s*sqrt(1-p^2)>0";
                case FrameworksErrorCode::SSVI_INVALID_PARAM: 
                    return "Invalid SSVI parameters. Either abs(rho)>1, nu<0, or gamma <0 or gamma>1.";
                case FrameworksErrorCode::VASICEK_INVALID_PARAM: 
                    return "Invalid Vasicek parameters. Either Kappa is non-positive or sigma is negative.";
                case FrameworksErrorCode::BLACK_SCHOLES_MODEL_INVALID_PARAM: 
                    return "Invalid black scholes model parameters. Either Spot, Strike time to maturity or sigma are given non-positive.";
                case FrameworksErrorCode::LETS_BE_RATIONAL_INVALID_BS_PRICE: 
                    return "Black scholes price must be a positive member in order to get an implied volatility mapping.";
                case FrameworksErrorCode::BARONE_ADESI_WHALEY_MODEL_INVALID_PARAM: 
                    return "Invalid Barone Adesi Whaley model parameters. Either Spot, Strike time to maturity or sigma are given non-positive.";
                case FrameworksErrorCode::HESTON_INVALID_PARAM:
                    return "Invalid Heston parameters.";
                case FrameworksErrorCode::HESTON_UNKNOWN_DISCRETIZATION_METHOD:
                    return "Unknown Heston discretization method.";
            }
            return "Unknown frameworks error.";
        };
};
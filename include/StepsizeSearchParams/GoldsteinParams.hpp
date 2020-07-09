#ifndef GOLDSTEIN_PARAMS_HPP
#define GOLDSTEIN_PARAMS_HPP

#include <StepsizeSearchParams/InaccurateSearchParams.hpp>

namespace NLOP {
/// @class NLOP::GoldsteinParams
/// @brief Params for Goldstein method
class GoldsteinParams: public InaccurateSearchParams
{
public:
    /// @brief Constructor
    GoldsteinParams() { setDefaults(); }

    /// @brief Use default stepsize searcher params
    void setDefaults() override
    {
        alpha = 1.5;
        beta = 0.5;
        init_lambda_factor = 0.1;
        rho = 0.1;
    }

    /// @brief Getters and setters for params
    void setRho(const double value) { rho = value; }

    double getRho() const { return rho; }

private:
    double rho;
};
}

#endif

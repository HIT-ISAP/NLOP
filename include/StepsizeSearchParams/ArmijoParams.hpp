#ifndef ARMIJO_PARAMS_HPP
#define ARMIJO_PARAMS_HPP

#include <StepsizeSearchParams/InaccurateSearchParams.hpp>

namespace NLOP {
/// @class NLOP::ArmijoParams
/// @brief Params for Armijo method
class ArmijoParams: public InaccurateSearchParams
{
public:
    /// @brief Constructor
    ArmijoParams() { setDefaults(); }

    /// @brief Use default stepsize searcher params
    void setDefaults() override
    {
        alpha = 1.5;
        beta = 0.5;
        rho = 0.1;
        mu = 10;
    }

    /// @brief Getters and setters for params
    void setRho(const double value) { rho = value; }
    void setMu(const double value) { mu = value; }

    double getRho() const { return rho; }
    double getMu() const { return mu; }

private:
    double rho;     // rho <= 0.1
    double mu;      // mu = 5 ~ 10
};
}

#endif

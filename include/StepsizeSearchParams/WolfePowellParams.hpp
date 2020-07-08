#ifndef WOLFEPOWELL_PARAMS_HPP
#define WOLFEPOWELL_PARAMS_HPP

#include <StepsizeSearchParams/InaccurateSearchParams.hpp>

namespace NLOP {
/// @class NLOP::WolfePowellParams
/// @brief Params for Wolfe-Powell method
class WolfePowellParams: public InaccurateSearchParams
{
public:
    /// @brief Constructor
    WolfePowellParams() { setDefaults(); }

    /// @brief Use default stepsize searcher params
    void setDefaults() override
    {
        alpha = 1.5;
        beta = 0.5;
        rho = 0.1;
        sigma = 0.8;
    }

    /// @brief Getters and setters for params
    void setRho(const double value) { rho = value; }
    void setSigma(const double value) { sigma = value; }

    double getRho() const { return rho; }
    double getSigma() const { return sigma; }

private:
    double rho;        // rho <= 0.1
    double sigma;      // sigma ~ (rho, 1)
};
}

#endif
